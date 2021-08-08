#!/bin/bash
##########################################################################
#  2019        Anton Gladky                                              #
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################
#
# invoke:
#       scripts/clang-format-all.sh ./Directory   optionally-clang-format-executable
#       scripts/clang-format-all.sh ./SingleFile  optionally-clang-format-executable

FORMATTER=clang-format

git diff-index --quiet HEAD --
CLEAN=$?
if [ ${CLEAN} -ne 0 ]; then
	echo "There are uncommitted changes in git directory. Do 'git stash' or 'git commit'."
	echo "Reformatting on top of uncommited files will only cause trouble."
	exit 1
fi

if [ "$#" -ne 1 ]; then
	if [ "$#" -eq 2 ]; then
		echo "Assuming that the clang-format executable is ${2}"
		FORMATTER=${2}
	else
		echo "Please invoke this script with either one file or one directory as argument."
		echo "Optionally as second argument pass the clang-format executable to use."
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
			echo "Formatting Python with: autopep8 --max-line-length 160 --ignore=E401 -i ${1}"
			autopep8 --max-line-length 160 --ignore=E401 -i ${1}
		else
			echo "Formatting C++ with: ${FORMATTER} -i ${1}"
			${FORMATTER} -i ${1}
		fi
		finish-print-stats
	else
		echo "${1} is neither a file nor a directory"
		exit 1
	fi
else
	echo "Formatting directory: ${FORMATTER} -i ${1}"
	find ${1} -iname *.cpp -o -iname *.hpp -o -iname *.ipp -o -iname *.h | xargs ${FORMATTER} -i
	echo "clang-format finished, now running autopep8 (from package python3-autopep8). It is slower, so uses 8 cores with: -P 8"
	find ${1} -iname *.py -print0 | xargs -0 -P 8 -I'{}' autopep8 --max-line-length 160 --ignore=E401 -i {}
	finish-print-stats
fi

# note: reformatting of the files:
#		-o -name '*.hpp.in' \
#		-o -name '*.cpp.in' \
#		-o -name '*.ipp.in' \
# is temporarily disabled. The @TEMPLATE_FLOW_NAME@ does not play well with clang-format

