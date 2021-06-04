import numpy as np
import os, sys
try:
	import cv2
except:
	print("Sorry, you don't have openCV library which is required for automated screenshot comparison. Please analyze screenshots manually.")
	os._exit(0)	#os._exit(1) use 1 for error and 0 for ok


screenshots = []
comparables = []
threshold = 150000	#Number of pixels which have to be different in order to draw differentiation
					#FIXME: Current threshold is really high due to inspector and controller window sizes and contents differing between gui tests
def getPngNames(path):
	pngList = []
	files=os.listdir(path)
	
	for png in files:
		if(png[len(png)-4:]==".png"):	#we are looking for png's only thus we check if file is a .png
			pngList.append(png)
	
	pngList.sort()
	return pngList

#Hhre we get the screenshots list
screenshots = getPngNames(path=sys.argv[-1])	#path to screenshot's folder should be given as an argument to command calling this script

#here we get the comparables list
comparables = getPngNames(path=sys.argv[-1]+"../compareSources")

#here we compare screenshots and comparables
#NOTE: script presumes that names of files are sorted, no name match checking is done (maybe it should?)
i = 0
while i < len(screenshots):
	comparable = cv2.imread(sys.argv[-1]+"../compareSources/"+comparables[i], 0)
	screenshot = cv2.imread(sys.argv[-1]+screenshots[i], 0)
	
	if(comparable.shape == screenshot.shape):
		diff = cv2.subtract(screenshot, comparable)
		error = cv2.countNonZero(diff)
		if(error > threshold):
			print("There are differences between: ", screenshots[i], " and ", comparables[i])
			print("Pixels diffrentiating: ", error)
			os._exit(1)
		else:
			print(screenshots[i], " ok")

	else:
		print("can't compare pictures with different sizes")
		os._exit(1)
	i+=1


