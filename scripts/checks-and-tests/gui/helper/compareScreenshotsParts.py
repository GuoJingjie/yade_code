import numpy as np
import os, sys
try:
	import cv2
except:
	print("Sorry, you don't have openCV library which is required for automated screenshot comparison. Please analyze screenshots manually.")
	os._exit(0)	#os._exit(1) use 1 for error and 0 for ok


screenshotNames = []
comparableNames = []
thresholdDict = dict({'view': 60000, 'console': 110000, 'controller': 110000, 'inspector': 100000})	#Number of pixels which have to be different in order to draw differentiation
					#FIXME: Current threshold is really high due to inspector and controller window sizes and contents differing between gui tests (and blinkHighlight on yellow)
def getPngNames(path):
	pngList = []
	files=os.listdir(path)
	
	for png in files:
		if(png[len(png)-4:]==".png"):	#we are looking for png's only thus we check if file is a .png
			pngList.append(png)
	
	pngList.sort()
	return pngList

#Hhre we get the screenshotNames list
screenshotNames = getPngNames(path=sys.argv[-1])	#path to screenshot's folder should be given as an argument to command calling this script

#here we get the comparableNames list
comparableNames = getPngNames(path=sys.argv[-1]+"../compareSources")
i=0
j=0
while i < len(screenshotNames):
	comparable = cv2.imread(sys.argv[-1]+"../compareSources/"+comparableNames[i], cv2.IMREAD_COLOR)
	screenshot = cv2.imread(sys.argv[-1]+screenshotNames[i], cv2.IMREAD_COLOR)
	if screenshot.shape == comparable.shape:
		print(f"Comparing: {screenshotNames[i]} with {comparableNames[i]}")
		view1 = comparable[0:560, 0:550]
		console1 = comparable[560:1600, 0:550]
		controller1 = comparable[0:1600, 550:1050]
		inspector1 = comparable[0:1600, 1050:]
		view2 = screenshot[0:560, 0:550]
		console2 = screenshot[560:1600, 0:550]
		controller2 = screenshot[0:1600, 550:1050]
		inspector2 = screenshot[0:1600, 1050:]
		comparableParts = [(view1,'view'), (console1,'console'), (controller1,'controller'), (inspector1,'inspector')]	#not really perfect solution for threshold-part synchronization
		screenshotParts = [(view2,'view'), (console2,'console'), (controller2,'controller'), (inspector2,'inspector')]	#but better than depending on hardcoded order

		while j < len(screenshotParts):
			if(screenshotParts[j][0].shape == comparableParts[j][0].shape):
				partDiff = cv2.subtract(screenshotParts[j][0], comparableParts[j][0])
				r,g,b = cv2.split(partDiff)
				error = (cv2.countNonZero(r)+cv2.countNonZero(g)+cv2.countNonZero(b))/3	#arithmetic average of errors from all 3 colours (could try to replace with weighted average)

				if(error > thresholdDict[screenshotParts[j][1]]):
					print(f"There are differences between: {screenshotNames[j]} {screenshotParts[j][1]} and {comparableNames[j]} {comparableParts[j][1]}")
					print(f"Pixels diffrentiating: {error}")
					os._exit(1)
			else:
				print("can't compare pictures with different sizes")
				os._exit(1)
			j+=1
	else:
		print("can't compare pictures with different sizes")
		os._exit(1)

	i+=1
