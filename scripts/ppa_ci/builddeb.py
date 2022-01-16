#!/usr/bin/env python3

import argparse, git, shutil
import tempfile, tarfile, fileinput
import email.utils, datetime, os

parser = argparse.ArgumentParser(description='Build yadedaily packages.')
parser.add_argument("-d", help="target distribution", action="store", dest="dist", default='buster', type=str)
parser.add_argument("-p", help="numerical precision", action="store", dest="prec", default='double', type=str)
args = parser.parse_args()
dt = datetime.datetime.now()
precName=('' if args.prec == 'double' else ("-"+args.prec))

# Define variables
dirpath = "./deb"
#dirpath = "./deb_%s/"%args.dist
#dirpath = tempfile.mkdtemp()
dirpathyade = dirpath + '/yadedaily'+precName+'/'

pipeline_ID = ""
if "CI_PIPELINE_IID" in os.environ:
    pipeline_ID = "-" + os.environ['CI_PIPELINE_IID']

repoups = git.Repo('.')
versiondebian = dt.strftime("%Y%m%d") + pipeline_ID + "~" + repoups.head.commit.hexsha[0:7] + "~" + args.dist + "1"
tarballname = 'yadedaily%s_%s.orig.tar.xz'%(precName,versiondebian)

# Create tempdir
os.mkdir(dirpath)
# Copy buildtree into the tmpdir
shutil.copytree('.', dirpathyade, ignore=shutil.ignore_patterns('.git', 'deb', '.ccache'))

# Create tarball
with tarfile.open('%s/%s'%(dirpath,tarballname), mode='w:xz') as out:
    print('Creating tarball... %s'%tarballname)
    out.add(dirpathyade + '/.', arcname='yadedaily', recursive=True)

# Copy debian-directory into the proper place
shutil.copytree('./scripts/ppa_ci/debian', dirpathyade + '/debian/')

with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("VERSION", versiondebian), end='')
with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("DISTRIBUTION", args.dist), end='')
with fileinput.FileInput(dirpathyade + '/debian/changelog', inplace=True) as file:
    for line in file:
        print(line.replace("DATE", email.utils.formatdate(localtime=True)), end='')
with fileinput.FileInput(dirpathyade + '/debian/rules', inplace=True) as file:
    for line in file:
        print(line.replace("VERSIONYADEREPLACE", versiondebian), end='')

# Prepare high precision package.
handlePrec={'long-double':('18',2), 'float128':('33',1), 'mpfr150':('150',1)}
if(args.prec in handlePrec):
    # Replace all instances of text 'yadedaily' with 'yadedaily'+precName. In all files
    print ('Building package with high precision: '+str(args.prec))
    with fileinput.FileInput(dirpathyade + '/debian/rules', inplace=True) as file:
        for line in file:
            print(line.replace('-DSUFFIX="daily"', '-DSUFFIX="daily'+precName+'" -DENABLE_MPFR=ON -DREAL_DECIMAL_PLACES='+handlePrec[args.prec][0]), end='')
    for fname in [  'control', 'rules', 'changelog', 'yadedaily.menu', 'yadedaily.desktop'
                  , 'yadedaily-doc.doc-base', 'yadedaily.install']:
        with fileinput.FileInput(dirpathyade + '/debian/' + fname, inplace=True) as file:
            for line in file:
                print(line.replace('yadedaily', 'yadedaily'+precName), end='')
    # 'python3-yadedaily.pyinstall' is special, because python modules cannot have '-' in their name. Use '_' instead.
    with fileinput.FileInput(dirpathyade + '/debian/python3-yadedaily.pyinstall', inplace=True) as file:
        for line in file: # workaround the problem that string.replace allows only to specify how many first occurrences to replace.
            print(line.replace('yadedaily', 'yadedaily'+precName.replace('-','_')).replace('_','-', handlePrec[args.prec][1] ), end='')
    # Use faster MPFR instead of boost::cpp_bin_float
    # so modify the debian/control file dependencies: add libmpfr-dev libmpfrc++-dev libmpc-dev
    with fileinput.FileInput(dirpathyade + '/debian/control', inplace=True) as file:
        for line in file:
            print(line.replace(' python3-numpy', ' libmpfr-dev,\n libmpfrc++-dev,\n libmpc-dev,\n python3-numpy'), end='')
    # Then rename files which have 'yadedaily' in their name.
    for fname in [ 'libyadedaily.install' ,'python3-yadedaily.pyinstall' ,'yadedaily-doc.doc-base'
                  ,'yadedaily-doc.docs' ,'yadedaily.desktop' ,'yadedaily.docs' ,'yadedaily.examples'
                  ,'yadedaily.install' ,'yadedaily.menu' ,'yadedaily_16x16.xpm' ,'yadedaily_32x32.xpm']:
        os.rename(dirpathyade + '/debian/' + fname, dirpathyade + '/debian/' + fname.replace('yadedaily', 'yadedaily'+precName))

print (versiondebian)
print (dirpath)
#shutil.rmtree(dirpath)
