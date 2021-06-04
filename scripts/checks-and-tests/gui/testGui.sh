#!/bin/bash

# This check is to be called inside xvfb-run, so that it has a working Xserver in which a simple yade GUI session can be started, the
# testGui.py script is only a slightly modified simple-scene-energy-tracking.py example.
# the screenshots have to be taken inside yade session, while the GUI windows are open.

#YADE_EXECUTABLE=install/bin/yade-ci			<- this is for CI, uncomment when finished testing
#GUI_TESTS_PATH=scripts/checks-and-tests/gui  FIXME: Path should be taken out of the call to the script	<- this is for CI, uncomment when finished testing
#YADE_EXECUTABLE=yade-2020-07-27.git-23b17b5  <- non QM version of yade
YADE_EXECUTABLE=~/Programs/yade_QM/install/bin/yade-2020-08-21.git-9ac74cd
GUI_TESTS_PATH=.

# You can test locally using this script, just change YADE_EXECUTABLE into something that works for you:
#YADE_EXECUTABLE=./examples/yade
#
# Also remember that default path for GuiTest files is scripts/checks-and-tests/gui/* so it to something that works for you e.g. GUI_TESTS_PATH=./
#
# then launch this command:
#   xvfb-run -a -s "-screen 0 1600x1200x16" scripts/checks-and-tests/gui/testGui.sh
#
# or just this command:
#   scripts/checks-and-tests/gui/testGui.sh


testTool () {
	WHICHtool=`which $1`
	echo -e "is $1 present? We found this: ${WHICHtool}"
	ls -la $2
	if [[ ${3} == "ERROR_OK" ]] ; then
		echo "OK: ${1} presence is not obligatory."
	else
		if [[ $2 == "${WHICHtool}" ]] ; then
			echo " OK."
		else
			if [[ $3 == "${WHICHtool}" ]] ; then
				echo " OK (second path)"
			else
				echo "ERROR: $2 is not \"${WHICHtool}\", this script is too stupid for that, aborting."
				exit 1
			fi
		fi
		if [[ ! -f $2 ]] ; then
			echo "ERROR: $2 is missing, aborting."
			exit 1
		fi
	fi
}

testTool "xterm"   "/usr/bin/xterm"
testTool "scrot"   "/usr/bin/scrot"
testTool "xdotool" "/usr/bin/xdotool"
testTool "bash"    "/bin/bash"        "/usr/bin/bash"
testTool "gdb"     "/usr/bin/gdb"     "ERROR_OK"

echo -e "\n\n=== Will now test inside xterm, all usefull output, including gdb crash backtrace, will be on screenshots ===\n\n"

mkdir -p screenshots

# FIXME: this should be deduced automatically from the files matching pattern testGui*.py, see also testGui.py
#        currently these names are written manually inside:
#          *  scripts/checks-and-tests/gui/testGuiEmpty.py
#          *  scripts/checks-and-tests/gui/testGuiSimple.py
#	declare -a TESTS=( "Empty" "Simple" )		<- old version, to be deleted if new approach is accepted
#FIX: line below runs python script which gets the names of files and passes them to bash variable
#Path to tests needs to be passed as an argument for it to find correct place with files

TESTS=($(python3 ${GUI_TESTS_PATH}/helper/readNames.py ${GUI_TESTS_PATH} | tr -d '[],'))

for TestFile in ${TESTS[@]}; do

	LOGFILE="${GUI_TESTS_PATH}/screenshots/testGui_${TestFile}.txt"
	tail -F ${LOGFILE} &
	TAIL_PID=$!

	echo -e "******************************************\n*** Testing file testGui${TestFile}.py ***\n******************************************\nLog in file: ${LOGFILE}\ntail pid:${TAIL_PID}\n"

	/usr/bin/xterm -l -xrm "XTerm*logFile:${LOGFILE}" -geometry 100x48+5+560  -e /bin/bash -c "${YADE_EXECUTABLE} ${GUI_TESTS_PATH}/testGui${TestFile}.py"

	# FIXME: the idea is to have a screenshot from outside of yade. But taking a screenshot after it finished (crashed, or by normal exit)
	#        will just produce an empty screenshot. It has to be done in a different way.
	# scrot -z scrBash01.png

	mv scr*.png screenshots
	sleep 0.25
	echo -e "******************************************\n*** Finished file testGui${TestFile}.py ***\n******************************************\n"
	kill -9 ${TAIL_PID}

# FIXME : this number 14 is hardcoded in scripts/checks-and-tests/gui/helper/testGuiHelper.py as self.maxTestNum=14
	if [[ ! -f screenshots/scr_${TestFile}_14.png ]] ; then
	    echo "File screenshots/scr_${TestFile}_14.png is missing, aborting."
	    exit 1
	fi

done

