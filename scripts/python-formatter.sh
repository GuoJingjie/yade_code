#!/bin/bash
##########################################################################
#  2021        Anton Gladky                                              #
#  2021        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
#
# Let us use this script on python scripts added or changed in new Merge
# Requests. Eventually all the scripts will be nicely formatted.
#
# Warning: we are still looking for optimal flags (or a better formatter).
#
# invoke:
#       scripts/python-formatter.sh ./Directory   optionally-yapf3-executable
#       scripts/python-formatter.sh ./SingleFile  optionally-yapf3-executable

FORMATTER=yapf3

git diff-index --quiet HEAD --
CLEAN=$?
if [ ${CLEAN} -ne 0 ]; then
	echo "There are uncommitted changes in git directory. Do 'git stash' or 'git commit'."
	echo "Reformatting on top of uncommited files will only cause trouble."
	exit 1
fi

if [ "$#" -ne 1 ]; then
	if [ "$#" -eq 2 ]; then
		echo "Assuming that the yapf3 executable is ${2}"
		FORMATTER=${2}
	else
		echo "Please invoke this script with either one file or one directory as argument."
		echo "Optionally as second argument pass the yapf3 executable to use."
		exit 1
	fi
fi

function finish-print-stats {
	echo "Formatting finished. If the changes are a mistake you can call: 'git reset --hard'."
	git diff --shortstat
	exit 0
}

if [ ! -d "${1}" ]; then
	if [ -f "${1}" ]; then
		FILE=${1}
		extension=${FILE##*.}
		echo "File extension: "${extension}
		if [ "${extension}" = "py"  ]; then
			echo "Formatting Python with yapf3"
			${FORMATTER} --style='{based_on_style: google; column_limit: 160; use_tabs: True; dedent_closing_brackets: True; indent_width = 8; continuation_indent_width = 8 }' -i ${1}
		else
			echo "Please use scripts/clang-formatter.sh for formatting C++ code"
		fi
		finish-print-stats
	else
		echo "${1} is neither a file nor a directory"
		exit 1
	fi
else
	echo "Running yapf3 (from packages python3-yapf yapf3). It is slow, so use 8 cores with: -P 8"
	find ${1} -iname *.py -print0 | xargs -0 -P 8 -I'{}' ${FORMATTER} --style='{based_on_style: google; column_limit: 160; use_tabs: True; dedent_closing_brackets: True; indent_width = 8; continuation_indent_width = 8 }' -i {}
	finish-print-stats
fi

