import numpy as np
import os, sys
print('\033[93m Checking screenshots, paths: \n', sys.argv[-1], '\n', sys.argv[-2], ' \033[0m')
try:
	import cv2
except:
	print("'\033[93m'+Please install python3-opencv package for automated screenshot comparison. Screenshots for manual comparison are in directory: ", sys.argv[-1], '\033[0m')
	sys.stdout.flush()
	os._exit(0)  # os._exit(1) use 1 for error and 0 for ok


screenshotNames = []
# Number of pixels which have to be different in order to get attention
thresholdDict = dict({'view': 60000, 'console': 110000, 'controller': 110000, 'inspector': 100000})
# FIXME: Current threshold is really high due to inspector and controller window sizes and contents differing between gui tests (and blinkHighlight on yellow)


def printFlushExit(msg, code):
	print(msg)
	sys.stdout.flush()
	os._exit(code)


def getPngNames(path):
	pngList = []
	files = os.listdir(path)

	for png in files:
		if(png[len(png) - 4:] == ".png"):  # we are looking for png's only thus we check if file is a .png
			pngList.append(png)

	pngList.sort()
	print("Path:  ", path)
	print("Files: ", pngList)
	return pngList


# Hhre we get the screenshotNames list
screenshotNames = getPngNames(path=sys.argv[-1])  # path to screenshot's folder should be given as an argument to command calling this script

i = 0
j = 0
while i < len(screenshotNames):
	print('\033[93m Checking screenshot ', screenshotNames[i], ' \033[0m')
	# We must only check if comparable file exists as screenshot file must exist if it is in the screenshotNames list as generation requirement
	if not os.path.isfile(sys.argv[-2] + "/compareSources/" + screenshotNames[i]):
		printFlushExit('\033[91m' + "ERROR: file, " + screenshotNames[i] + " does not exist in compareSources folder."'\033[0m', 1)
	# screenshotNames list is used for comparable generation aswell as compareSources and screenshots both have the exact same filenames (this prevents wrong files comparison)
	comparable = cv2.imread(sys.argv[-2] + "/compareSources/" + screenshotNames[i], cv2.IMREAD_COLOR)
	screenshot = cv2.imread(sys.argv[-1] + "/" + screenshotNames[i], cv2.IMREAD_COLOR)
	if screenshot.shape == comparable.shape:
		view1 = comparable[0:560, 0:550]
		console1 = comparable[560:1600, 0:550]
		controller1 = comparable[0:1600, 550:1050]
		inspector1 = comparable[0:1600, 1050:]
		view2 = screenshot[0:560, 0:550]
		console2 = screenshot[560:1600, 0:550]
		controller2 = screenshot[0:1600, 550:1050]
		inspector2 = screenshot[0:1600, 1050:]
		# not perfect solution for threshold-part synchronization
		# but better than depending on hardcoded order
		comparableParts = [(view1, 'view'), (console1, 'console'), (controller1, 'controller'), (inspector1, 'inspector')]
		screenshotParts = [(view2, 'view'), (console2, 'console'), (controller2, 'controller'), (inspector2, 'inspector')]

		while j < len(screenshotParts):
			if(screenshotParts[j][0].shape == comparableParts[j][0].shape):
				partDiff = cv2.subtract(screenshotParts[j][0], comparableParts[j][0])
				r, g, b = cv2.split(partDiff)
				# arithmetic average of errors from all 3 colours (could try to replace with weighted average)
				error = (cv2.countNonZero(r) + cv2.countNonZero(g) + cv2.countNonZero(b)) / 3

				if(error > thresholdDict[screenshotParts[j][1]]):
					printFlushExit("There are differences between: " + screenshotNames[j] + " " + screenshotParts[j][1]
					               + " and " + screenshotNames[j] + " " + comparableParts[j][1] + "Pixels diffrent: " + error, 1)
			else:
				printFlushExit("can't compare pictures with different sizes", 1)
			j += 1
	else:
		printFlushExit("can't compare pictures with different sizes", 1)

	i += 1
