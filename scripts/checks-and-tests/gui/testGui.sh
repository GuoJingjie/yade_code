#!/bin/bash
# This check is to be called inside xvfb-run, so that it has a working Xserver in which a simple yade GUI session can be started,
# the simplest of those checks, testGuiEmpty.py or testGuiSimple.py are only a slightly modified simple-scene-energy-tracking.py example.
# The screenshots are taken inside yade session, while the GUI windows are open, by calling scrot from helper/testGuiHelper.py

YADE_EXECUTABLE=install/bin/yade-ci
GUI_TESTS_PATH=scripts/checks-and-tests/gui
CREATE_NEW_SCREENSHOTS=screenshots
REFERENCE_SCREENSHOTS=data/checks-and-tests/gui/screenshots

# You can test this locally using this script, just change YADE_EXECUTABLE into something that works for you:
#
#YADE_EXECUTABLE=./examples/yade
#
# then launch this command (inside xvfb):
#   xvfb-run -a -s "-screen 0 1600x1200x16" scripts/checks-and-tests/gui/testGui.sh
#
# or just this command (to see it happening on your desktop):
#   scripts/checks-and-tests/gui/testGui.sh


# This function checks if necessary tools are installed.
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

echo -e "\n\n=== Will now test inside xterm, all useful output, including gdb crash backtrace ===\n=== will be on screenshots and in the xterm logs: /screenshots/*.txt             ===\n\n"

mkdir -p ${CREATE_NEW_SCREENSHOTS} # screenshots

# Loop over files ${GUI_TESTS_PATH}/testGui*py
# TODO: tell the called scripts what name to use "Empty" "Simple", etc. currently these names are written manually inside
#  * scripts/checks-and-tests/gui/testGuiEmpty.py
#  * scripts/checks-and-tests/gui/testGuiSimple.py
#  * etc.
export TIMEFORMAT='[93m Real time spent: %3lR [0m'
for FileName in ${GUI_TESTS_PATH}/testGui*py; do
	TestFile=($(echo $FileName | sed -e "s@${GUI_TESTS_PATH}/testGui@@g" | sed -e 's/.py$//g' ))

	LOGFILE="${CREATE_NEW_SCREENSHOTS}/testGui_${TestFile}.txt"
	tail -F ${LOGFILE} &
	TAIL_PID=$!

	echo -e "******************************************\n*** Testing file testGui${TestFile}.py ***\n******************************************\nLog in file: ${LOGFILE}\ntail pid:${TAIL_PID}\n"

	time /usr/bin/xterm -l -xrm "XTerm*logFile:${LOGFILE}" -geometry 100x48+5+560  -e /bin/bash -c "${YADE_EXECUTABLE} ${GUI_TESTS_PATH}/testGui${TestFile}.py"

	# The idea here is to have a screenshot from outside of yade. But taking a screenshot after it finished (crashed, or by normal exit)
	# will just produce an empty screenshot. It has to be done by calling scrot -z from inside scripts/checks-and-tests/gui/helper/testGuiHelper.py

	# Running is finished, move all created screenshots to the artifacts directory.
	mv scr_*.png ${CREATE_NEW_SCREENSHOTS}
	sleep 0.25
	echo -e "******************************************\n*** Finished file testGui${TestFile}.py ***\n******************************************\n"
	kill -9 ${TAIL_PID}

	# If running was a success then a file testGui_${TestFile}_OK_or_Skipped.txt, it is created by scripts/checks-and-tests/gui/helper/testGuiHelper.py, at the end of screenshotEngine function.
	echo -e "*** Checking if it was a success        ***"
	if [[ ! -f ${CREATE_NEW_SCREENSHOTS}/testGui_${TestFile}_OK_or_Skipped.txt ]] ; then
		echo "File ${CREATE_NEW_SCREENSHOTS}/testGui_${TestFile}_OK_or_Skipped.txt is missing, aborting."
		exit 1
	else
		ls -la ${CREATE_NEW_SCREENSHOTS}/testGui_${TestFile}_OK_or_Skipped.txt
		echo -e "*** OK ***"
	fi
done

sleep 1
echo -e "******************************************\n*** Checking screenshots now ***\n******************************************\n"
echo "This is CI job: ${CI_JOB_NAME}"
if [[ ${CI_JOB_NAME} == tstHP* || ${CI_JOB_NAME} == make_asan_HP || ${CI_JOB_NAME} == test_opposite || ${CI_JOB_NAME} == test_SSE ]]; then
	# FIXME : high precision tests produce different output in terminal. This is related to https://gitlab.com/yade-dev/trunk/-/issues/203
	#         for now I temporarily use higher xterm tolerance for this test. Remember to remove this once #203 is fixed.
	# NOTE  : The test_opposite prints different messages in the terminal. We could add different reference screenshots. For now let's just use higher xterm tolerance.
	# NOTE  : test_SSE also needs higher tolerance, because of Eigen messages about memory alignment.
	python3 ${GUI_TESTS_PATH}/helper/compareScreenshotsParts.py ${REFERENCE_SCREENSHOTS} ${CREATE_NEW_SCREENSHOTS} 45 || { sleep 1 ; exit 1; }
else
	# Smallest xterm tolerance is 5, enough for different session cookie and some small variation in messages.
	python3 ${GUI_TESTS_PATH}/helper/compareScreenshotsParts.py ${REFERENCE_SCREENSHOTS} ${CREATE_NEW_SCREENSHOTS}  5 || { sleep 1 ; exit 1; }
fi
echo -e "******************************************\n*** Checking screenshots finished ***\n******************************************\n"

