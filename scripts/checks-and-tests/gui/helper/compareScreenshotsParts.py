import numpy as np
import os, sys
scrDir = sys.argv[-1]
refDir = sys.argv[-2]
print('\033[93m Checking screenshots\n  reference directory : ', refDir, '\n  new screenshots dir : ', scrDir  , ' \033[0m (can be found in job artifacts)')
try:
	import cv2
except:
	print("'\033[93m'+Please install python3-opencv package for automated screenshot comparison. Screenshots for manual comparison are in directory: ", scrDir, '\033[0m')
	sys.stdout.flush()
	os._exit(0)  # os._exit(1) use 1 for error and 0 for ok

screenshotNames = []
# Number of pixels which have to be different in order to get attention
thresholdDict = dict({'view': 60000, 'term': 110000, 'cont': 110000, 'insp': 100000})
maxEncountered= dict({'view': 0    , 'term': 0     , 'cont': 0     , 'insp': 0     })

def printFlushExit(msg, code):
	print(msg)
	sys.stdout.flush()
	os._exit(code)

def getPngNames(path):
	pngList = []
	files = os.listdir(path)
	for png in files:
		if(png[len(png) - 4:] == ".png"):
			pngList.append(png)
	pngList.sort()
	#print("Path:  ", path)
	#print("Files: ", pngList)
	return pngList

screenshotNames = getPngNames(path=refDir)

print('                      ----- MAX number of different pixels -----')
print('                      '+str(thresholdDict['view']).ljust(11)+str(thresholdDict['term']).ljust(11)+str(thresholdDict['cont']).ljust(11)+str(thresholdDict['insp']).ljust(11))
print('                      ------------------------------------------')
print('Screenshot            View       Terminal   Controller Inspector')
hadError    = False

try:
	for refScr in screenshotNames:
		# Check if file exists
		if not os.path.isfile(scrDir   + "/" + refScr):
			printFlushExit('\033[91m' + "ERROR: file, " + refScr + " does not exist in compareSources folder."'\033[0m', 1)
		# Both screenshots have the same filenames
		scr = cv2.imread(scrDir + "/" + refScr, cv2.IMREAD_COLOR)
		ref = cv2.imread(refDir + "/" + refScr, cv2.IMREAD_COLOR)
		print(refScr.ljust(22), end='')
		if ref.shape == scr.shape:
			#             View                        Terminal                       Controller                       Inspector
			scrParts = [ (scr[0:560, 0:550],'view'), (scr[560:1200, 0:550],'term'), (scr[20:1120, 550:1050],'cont'), (scr[20:1120, 1050:1550],'insp') ]
			refParts = [ (ref[0:560, 0:550],'view'), (ref[560:1200, 0:550],'term'), (ref[20:1120, 550:1050],'cont'), (ref[20:1120, 1050:1550],'insp') ]
			j = 0
			rowHadError = False
			while j < len(refParts):
				if(refParts[j][0].shape == scrParts[j][0].shape):
					partDiff = cv2.subtract(refParts[j][0], scrParts[j][0])
					r, g, b = cv2.split(partDiff)
					# arithmetic average of errors from all 3 colours
					error       = int( (cv2.countNonZero(r) + cv2.countNonZero(g) + cv2.countNonZero(b)) / 3. )
					nowHasError = error > thresholdDict[refParts[j][1]]
					if((refScr.endswith('_01.png') or refScr.endswith('_03.png')) and j==0):
						# 01: Controller - before it is moved into correct position it can appear as having different sizes.
						# 03: View       - it may not appear on time when this screenshot is created
						nowHasError = False
					else:
						maxEncountered[refParts[j][1]] = max(maxEncountered[refParts[j][1]],error)
					rowHadError = rowHadError or nowHasError
					hadError    = hadError    or nowHasError
					print(('\033[91m' if nowHasError else '') + str(error).ljust(11) + ('\033[0m' if nowHasError else ''), end=('' if (j!=3 or rowHadError) else '\n'))
				else:
					printFlushExit("can't compare pictures with different sizes", 1)
				j += 1
			if(rowHadError):
				print("Error")
		else:
			printFlushExit("can't compare pictures with different sizes", 1)
	
	print('                      ------------------------------------------')
	print('maximum error found   '+str(maxEncountered['view']).ljust(11)+str(maxEncountered['term']).ljust(11)+str(maxEncountered['cont']).ljust(11)+str(maxEncountered['insp']).ljust(11))
	
	if(hadError):
		printFlushExit("Errors detected!", 1)
except:
	print("Exception encountered during tests")
	os._exit(1)

