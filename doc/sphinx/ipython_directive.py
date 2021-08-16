# -*- coding: utf-8 -*-
"""Sphinx directive to support embedded IPython code.

This directive allows pasting of entire interactive IPython sessions, prompts
and all, and their code will actually get re-executed at doc build time, with
all prompts renumbered sequentially.

To enable this directive, simply list it in your Sphinx ``conf.py`` file
(making sure the directory where you placed it is visible to sphinx, as is
needed for all Sphinx directives).

By default this directive assumes that your prompts are unchanged IPython ones,
but this can be customized.  For example, the following code in your Sphinx
config file will configure this directive for the following input/output
prompts ``Yade [1]:`` and ``-> [1]:``::

 import ipython_directive as id
 id.rgxin =re.compile(r'(?:In |Yade )\[(\d+)\]:\s?(.*)\s*')
 id.rgxout=re.compile(r'(?:Out| ->  )\[(\d+)\]:\s?(.*)\s*')
 id.fmtin ='Yade [%d]:'
 id.fmtout=' ->  [%d]:'

 id.rc_override=dict(
   prompt_in1="Yade [\#]:",
   prompt_in2="     .\D..",
   prompt_out=" ->  [\#]:"
 )
 id.reconfig_shell()

 import ipython_console_highlighting as ich
 ich.IPythonConsoleLexer.input_prompt=
    re.compile("(Yade \[[0-9]+\]: )|(   \.\.\.+:)")
 ich.IPythonConsoleLexer.output_prompt=
    re.compile("(( ->  )|(Out)\[[0-9]+\]: )|(   \.\.\.+:)")
 ich.IPythonConsoleLexer.continue_prompt=re.compile("   \.\.\.+:")


ToDo
----

- Turn the ad-hoc test() function into a real test suite.
- Break up ipython-specific functionality from matplotlib stuff into better
  separated code.
- Make sure %bookmarks used internally are removed on exit.


Authors
-------

- John D Hunter: orignal author.
- Fernando Perez: refactoring, documentation, cleanups.
- VáclavŠmilauer <eudoxos-AT-arcig.cz>: Prompt generatlizations.
"""
from __future__ import print_function

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

# Stdlib
from future import standard_library
standard_library.install_aliases()

from builtins import range
from builtins import object
import io
import imp
import os
import re
import shutil
import sys
import warnings

# To keep compatibility with various python versions
try:
	from hashlib import md5
except ImportError:
	from md5 import md5

# Third-party
import matplotlib
import sphinx
from docutils.parsers.rst import directives

matplotlib.use('Agg')

# Our own
import IPython
from IPython.Shell import MatplotlibShell

#-----------------------------------------------------------------------------
# Globals
#-----------------------------------------------------------------------------

sphinx_version = sphinx.__version__.split(".")
# The split is necessary for sphinx beta versions where the string is
# '6b1'
sphinx_version = tuple([int(re.split('[a-z]', x)[0]) for x in sphinx_version[:2]])

COMMENT, INPUT, OUTPUT = list(range(3))
rc_override = {}
rgxin = re.compile('In \[(\d+)\]:\s?(.*)\s*')
rgxcont = re.compile('   \.+:\s?(.*)\s*')
rgxout = re.compile('Out\[(\d+)\]:\s?(.*)\s*')
fmtin = 'In [%d]:'
fmtout = 'Out[%d]:'
fmtcont = '   .\D.:'


#-----------------------------------------------------------------------------
# Functions and class declarations
#-----------------------------------------------------------------------------
def block_parser(part):
	"""
    part is a string of ipython text, comprised of at most one
    input, one ouput, comments, and blank lines.  The block parser
    parses the text into a list of::

      blocks = [ (TOKEN0, data0), (TOKEN1, data1), ...]

    where TOKEN is one of [COMMENT | INPUT | OUTPUT ] and
    data is, depending on the type of token::

      COMMENT : the comment string

      INPUT: the (DECORATOR, INPUT_LINE, REST) where
         DECORATOR: the input decorator (or None)
         INPUT_LINE: the input as string (possibly multi-line)
         REST : any stdout generated by the input line (not OUTPUT)


      OUTPUT: the output string, possibly multi-line
    """

	block = []
	lines = part.split('\n')
	N = len(lines)
	i = 0
	decorator = None
	while 1:

		if i == N:
			# nothing left to parse -- the last line
			break

		line = lines[i]
		i += 1
		line_stripped = line.strip()
		if line_stripped.startswith('#'):
			block.append((COMMENT, line))
			continue

		if line_stripped.startswith('@'):
			# we're assuming at most one decorator -- may need to
			# rethink
			decorator = line_stripped
			continue

		# does this look like an input line?
		matchin = rgxin.match(line)
		if matchin:
			lineno, inputline = int(matchin.group(1)), matchin.group(2)

			# the ....: continuation string
			#continuation = '   %s:'%''.join(['.']*(len(str(lineno))+2))
			#Nc = len(continuation)
			# input lines can continue on for more than one line, if
			# we have a '\' line continuation char or a function call
			# echo line 'print'.  The input line can only be
			# terminated by the end of the block or an output line, so
			# we parse out the rest of the input line if it is
			# multiline as well as any echo text

			rest = []
			while i < N:

				# look ahead; if the next line is blank, or a comment, or
				# an output line, we're done

				nextline = lines[i]
				matchout = rgxout.match(nextline)
				matchcont = rgxcont.match(nextline)
				#print "nextline=%s, continuation=%s, starts=%s"%(nextline, continuation, nextline.startswith(continuation))
				if matchout or nextline.startswith('#'):
					break
				elif matchcont:  #nextline.startswith(continuation):
					inputline += '\n' + matchcont.group(1)  #nextline[Nc:]
				else:
					rest.append(nextline)
				i += 1

			block.append((INPUT, (decorator, inputline, '\n'.join(rest))))
			continue

		# if it looks like an output line grab all the text to the end
		# of the block
		matchout = rgxout.match(line)
		if matchout:
			lineno, output = int(matchout.group(1)), matchout.group(2)
			if i < N - 1:
				output = '\n'.join([output] + lines[i:])

			block.append((OUTPUT, output))
			break

	return block


class EmbeddedSphinxShell(object):
	"""An embedded IPython instance to run inside Sphinx"""

	def __init__(self):

		self.cout = io.StringIO()

		IPython.Shell.Term.cout = self.cout
		IPython.Shell.Term.cerr = self.cout
		argv = ['-autocall', '0']
		self.user_ns = {}
		self.user_glocal_ns = {}

		self.IP = IPython.ipmaker.make_IPython(
		        argv,
		        self.user_ns,
		        self.user_glocal_ns,
		        embedded=True,
		        #shell_class=IPython.Shell.InteractiveShell,
		        shell_class=MatplotlibShell,
		        rc_override=dict(colors='NoColor', **rc_override)
		)

		self.input = ''
		self.output = ''

		self.is_verbatim = False
		self.is_doctest = False
		self.is_suppress = False

		# on the first call to the savefig decorator, we'll import
		# pyplot as plt so we can make a call to the plt.gcf().savefig
		self._pyplot_imported = False

		# we need bookmark the current dir first so we can save
		# relative to it
		self.process_input_line('bookmark ipy_basedir')
		self.cout.seek(0)
		self.cout.truncate(0)

	def process_input_line(self, line):
		"""process the input, capturing stdout"""
		#print "input='%s'"%self.input
		stdout = sys.stdout
		sys.stdout = self.cout
		#self.IP.resetbuffer()
		self.IP.push(self.IP.prefilter(line, 0))
		#self.IP.runlines(line)
		sys.stdout = stdout

	# Callbacks for each type of token
	def process_input(self, data, input_prompt, lineno):
		"""Process data block for INPUT token."""
		decorator, input, rest = data
		image_file = None
		#print 'INPUT:', data
		is_verbatim = decorator == '@verbatim' or self.is_verbatim
		is_doctest = decorator == '@doctest' or self.is_doctest
		is_suppress = decorator == '@suppress' or self.is_suppress
		is_savefig = decorator is not None and decorator.startswith('@savefig')

		input_lines = input.split('\n')

		#continuation = '   %s:'%''.join(['.']*(len(str(lineno))+2))
		#Nc = len(continuation)

		if is_savefig:
			saveargs = decorator.split(' ')
			filename = saveargs[1]
			outfile = os.path.join('_static/%s' % filename)
			# build out an image directive like
			# .. image:: somefile.png
			#    :width 4in
			#
			# from an input like
			# savefig somefile.png width=4in
			imagerows = ['.. image:: %s' % outfile]

			for kwarg in saveargs[2:]:
				arg, val = kwarg.split('=')
				arg = arg.strip()
				val = val.strip()
				imagerows.append('   :%s: %s' % (arg, val))

			image_file = outfile
			image_directive = '\n'.join(imagerows)

		# TODO: can we get "rest" from ipython
		#self.process_input_line('\n'.join(input_lines))

		ret = []
		is_semicolon = False

		for i, line in enumerate(input_lines):
			if line.endswith(';'):
				is_semicolon = True

			if i == 0:
				# process the first input line
				if is_verbatim:
					self.process_input_line('')
				else:
					# only submit the line in non-verbatim mode
					self.process_input_line(line)
				formatted_line = '%s %s' % (input_prompt, line)
			else:
				# process a continuation line
				if not is_verbatim:
					self.process_input_line(line)

				formatted_line = fmtcont.replace('\D', '.' * len(str(lineno))) + line  #'%s %s'%(continuation, line)

			if not is_suppress:
				ret.append(formatted_line)

		if not is_suppress:
			if len(rest.strip()):
				if is_verbatim:
					# the "rest" is the standard output of the
					# input, which needs to be added in
					# verbatim mode
					ret.append(rest)

		self.cout.seek(0)
		output = self.cout.read()
		if not is_suppress and not is_semicolon:
			ret.append(output)

		self.cout.truncate(0)
		return ret, input_lines, output, is_doctest, image_file
		#print 'OUTPUT', output  # dbg

	def process_output(self, data, output_prompt, input_lines, output, is_doctest, image_file):
		"""Process data block for OUTPUT token."""
		if is_doctest:
			submitted = data.strip()
			found = output
			if found is not None:
				ind = found.find(output_prompt)
				if ind < 0:
					raise RuntimeError('output prompt="%s" does not match out line=%s' % (output_prompt, found))
				found = found[len(output_prompt):].strip()

				if found != submitted:
					raise RuntimeError(
					        'doctest failure for input_lines="%s" with found_output="%s" and submitted output="%s"' %
					        (input_lines, found, submitted)
					)
				#print 'doctest PASSED for input_lines="%s" with found_output="%s" and submitted output="%s"'%(input_lines, found, submitted)

	def process_comment(self, data):
		"""Process data block for COMMENT token."""
		if not self.is_suppress:
			return [data]

	def process_block(self, block):
		"""
        process block from the block_parser and return a list of processed lines
        """

		ret = []
		output = None
		input_lines = None

		m = rgxin.match(str(self.IP.outputcache.prompt1).strip())
		lineno = int(m.group(1))

		input_prompt = fmtin % lineno
		output_prompt = fmtout % lineno
		image_file = None
		image_directive = None
		# XXX - This needs a second refactor.  There's too much state being
		# held globally, which makes for a very awkward interface and large,
		# hard to test functions.  I've already broken this up at least into
		# three separate processors to isolate the logic better, but this only
		# serves to highlight the coupling.  Next we need to clean it up...
		for token, data in block:
			if token == COMMENT:
				out_data = self.process_comment(data)
			elif token == INPUT:
				out_data, input_lines, output, is_doctest, image_file = self.process_input(data, input_prompt, lineno)
			elif token == OUTPUT:
				out_data = self.process_output(data, output_prompt, input_lines, output, is_doctest, image_file)
			if out_data:
				ret.extend(out_data)

		if image_file is not None:
			self.ensure_pyplot()
			command = 'plt.gcf().savefig("%s")' % image_file
			#print 'SAVEFIG', command  # dbg
			self.process_input_line('bookmark ipy_thisdir')
			self.process_input_line('cd -b ipy_basedir')
			self.process_input_line(command)
			self.process_input_line('cd -b ipy_thisdir')
			self.cout.seek(0)
			self.cout.truncate(0)
		return ret, image_directive

	def ensure_pyplot(self):
		if self._pyplot_imported:
			return
		self.process_input_line('import matplotlib.pyplot as plt')


# A global instance used below. XXX: not sure why this can't be created inside
# ipython_directive itself.
shell = EmbeddedSphinxShell()


def reconfig_shell():
	"""Called after setting module-level variables to re-instantiate
    with the set values (since shell is instantiated first at import-time
    when module variables have default values)"""
	global shell
	shell = EmbeddedSphinxShell()


def ipython_directive(
        name,
        arguments,
        options,
        content,
        lineno,
        content_offset,
        block_text,
        state,
        state_machine,
):

	debug = ipython_directive.DEBUG
	shell.is_suppress = 'suppress' in options
	shell.is_doctest = 'doctest' in options
	shell.is_verbatim = 'verbatim' in options

	#print 'ipy', shell.is_suppress, options
	parts = '\n'.join(content).split('\n\n')
	lines = ['.. sourcecode:: ipython', '']

	figures = []
	for part in parts:
		block = block_parser(part)

		if len(block):
			rows, figure = shell.process_block(block)
			for row in rows:
				lines.extend(['    %s' % line for line in row.split('\n')])

			if figure is not None:
				figures.append(figure)

	for figure in figures:
		lines.append('')
		lines.extend(figure.split('\n'))
		lines.append('')

	#print lines
	if len(lines) > 2:
		if debug:
			print('\n'.join(lines))
		else:
			#print 'INSERTING %d lines'%len(lines)
			state_machine.insert_input(lines, state_machine.input_lines.source(0))

	return []


ipython_directive.DEBUG = False


# Enable as a proper Sphinx directive
def setup(app):
	setup.app = app
	options = {
	        'suppress': directives.flag,
	        'doctest': directives.flag,
	        'verbatim': directives.flag,
	}

	app.add_directive('ipython', ipython_directive, True, (0, 2, 0), **options)


# Simple smoke test, needs to be converted to a proper automatic test.
def test():

	examples = [
	        r"""
In [9]: pwd
Out[9]: '/home/jdhunter/py4science/book'

In [10]: cd bookdata/
/home/jdhunter/py4science/book/bookdata

In [2]: from pylab import *

In [2]: ion()

In [3]: im = imread('stinkbug.png')

@savefig mystinkbug.png width=4in
In [4]: imshow(im)
Out[4]: <matplotlib.image.AxesImage object at 0x39ea850>
        
""",
	        r"""

In [1]: x = 'hello world'

# string methods can be
# used to alter the string
@doctest
In [2]: x.upper()
Out[2]: 'HELLO WORLD'

@verbatim
In [3]: x.st<TAB>
x.startswith  x.strip
""",
	        r"""

In [130]: url = 'http://ichart.finance.yahoo.com/table.csv?s=CROX\
   .....: &d=9&e=22&f=2009&g=d&a=1&br=8&c=2006&ignore=.csv'

In [131]: print url.split('&')
['http://ichart.finance.yahoo.com/table.csv?s=CROX', 'd=9', 'e=22', 'f=2009', 'g=d', 'a=1', 'b=8', 'c=2006', 'ignore=.csv']

In [60]: import urllib

""",
	        r"""\

In [133]: import numpy.random

@suppress
In [134]: numpy.random.seed(2358)

@doctest
In [135]: np.random.rand(10,2)
Out[135]:
array([[ 0.64524308,  0.59943846],
       [ 0.47102322,  0.8715456 ],
       [ 0.29370834,  0.74776844],
       [ 0.99539577,  0.1313423 ],
       [ 0.16250302,  0.21103583],
       [ 0.81626524,  0.1312433 ],
       [ 0.67338089,  0.72302393],
       [ 0.7566368 ,  0.07033696],
       [ 0.22591016,  0.77731835],
       [ 0.0072729 ,  0.34273127]])

""",
	        r"""
In [106]: print x
jdh

In [109]: for i in range(10):
   .....:     print i
   .....:
   .....:
0
1
2
3
4
5
6
7
8
9


""",
	        r"""

In [144]: from pylab import *

In [145]: ion()

# use a semicolon to suppress the output
@savefig test_hist.png width=4in
In [151]: hist(np.random.randn(10000), 100);


@savefig test_plot.png width=4in
In [151]: plot(np.random.randn(10000), 'o');
   """,
	        r"""
# use a semicolon to suppress the output
In [151]: plt.clf()

@savefig plot_simple.png width=4in
In [151]: plot([1,2,3])

@savefig hist_simple.png width=4in
In [151]: hist(np.random.randn(10000), 100);

""",
	        r"""
# update the current fig
In [151]: ylabel('number')

In [152]: title('normal distribution')


@savefig hist_with_text.png
In [153]: grid(True)

        """,
	]

	ipython_directive.DEBUG = True
	#options = dict(suppress=True)
	options = dict()
	for example in examples:
		content = example.split('\n')
		ipython_directive(
		        'debug',
		        arguments=None,
		        options=options,
		        content=content,
		        lineno=0,
		        content_offset=None,
		        block_text=None,
		        state=None,
		        state_machine=None,
		)


# Run test suite as a script
if __name__ == '__main__':
	test()
