import numpy as np
import os, sys
xtermTolerance  = int(sys.argv[-1])
scrDir          = sys.argv[-2]
refDir          = sys.argv[-3]
try:
	import cv2
except:
	print("\033[93m Please install python3-opencv package for automated screenshot comparison.\033[0m")
	sys.stdout.flush()
	os._exit(0)  # os._exit(1) use 1 for error and 0 for ok
print('\033[93m Checking screenshots\n  reference directory : ', refDir, '\n  new screenshots dir : ', scrDir  , ' \033[0m (can be found in job artifacts)')

screenshotNames = []
# Maximum allowed average pixel difference between screenhots.
thresholdDict = dict({'view': 5    , 'term': xtermTolerance , 'cont': 20    , 'insp': 20    })
maxEncountered= dict({'view': 0    , 'term': 0              , 'cont': 0     , 'insp': 0     })

def getAllowedError(viewName,scrNumber,testName):
	skipCheck = 1000000
	assert(viewName  in ['view','term','cont','insp'])
	assert(scrNumber in range(15))
	# Some exceptions are allowed here

	# Before Controller is moved into correct position it can appear as having different sizes, in the View area.
	if((viewName == 'view') and (scrNumber == 1)):
		return (skipCheck,'Start Controller ')

	# The View may not appear on time when this screenshot is created
	if((viewName == 'view') and (scrNumber == 3)):
		return (skipCheck,'Start View ')

	# The testGuiHopper.py and testGuiVideo.py are more random than other tests. Especially in different precisions.
	# Their goal is to test for crashes and when drawing and erasing, not for reproducibility of screenshots.
	if('Hopper' in testName):
		# In different precisions the clumps can go to different positions in 3D view
		if((viewName == 'view') and (scrNumber >= 4)):
			return (skipCheck,'Hpr View ')
		# and the interactions might be different.
		if((viewName == 'insp') and (scrNumber == 7)):
			return (skipCheck,'Hpr Insp ')
	if('Video' in testName):
		# In different precisions the clumps can go to different positions in 3D view
		if((viewName == 'view') and (scrNumber >= 3) and (scrNumber <= 8)):
			return (skipCheck,'Vid View ')
		# and the interactions might be different.
		if((viewName == 'insp') and (scrNumber == 7)):
			return (skipCheck,'Vid Insp ')

	return (thresholdDict[viewName],'')

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

print('                      ---------- MAX allowed mean error --------')
print('                      '+str(thresholdDict['view']).ljust(11)+str(thresholdDict['term']).ljust(11)+str(thresholdDict['cont']).ljust(11)+str(thresholdDict['insp']).ljust(11))
print('                      ------------------------------------------')
print('Screenshot            View       Terminal   Controller Inspector')
hadError    = False

try:
	for refScr in screenshotNames:
		# Check if file exists
		if not os.path.isfile(scrDir   + "/" + refScr):
			testFile = refScr[4:-7]
			skipFile = "testGui_"+testFile+"_OK_or_Skipped.txt"
			if(os.path.isfile(scrDir + "/" + skipFile)):
				# test is normally skipped
				print(refScr+" : "+skipFile)
				continue
			else:
				printFlushExit("\033[91m ERROR: file, " + refScr + " does not exist in compareSources folder.\033[0m", 1)
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
			rowMsg      = ''
			while j < len(refParts):
				if(refParts[j][0].shape == scrParts[j][0].shape):
					partDiff = cv2.subtract(refParts[j][0], scrParts[j][0])
					r, g, b = cv2.split(partDiff)
					# sum of mean error for each channel, see https://docs.opencv.org/3.4/d2/de8/group__core__array.html#ga191389f8a0e58180bb13a727782cd461
					error       = cv2.mean(r)[0] + cv2.mean(g)[0] + cv2.mean(b)[0]
					tol,msg     = getAllowedError(refParts[j][1], int(refScr[-6:-4]), refScr)
					nowHasError = error > tol
					rowMsg     += msg
					maxEncountered[refParts[j][1]] = max(maxEncountered[refParts[j][1]],error)
					rowHadError = rowHadError or nowHasError
					hadError    = hadError    or nowHasError
					# https://misc.flogisoft.com/bash/tip_colors_and_formatting
					print(
						  ('\033[91m' if nowHasError else ('\033[95m' if (msg != '') else ''))
						+ str(round(error,6)).ljust(11)
						+ ('\033[0m' if (nowHasError or rowMsg != '') else '')
						, end=('' if ((j != 3) or rowHadError or (rowMsg != '')) else '\n')
					)
					if(nowHasError):
						cv2.imwrite(scrDir+'/'+refScr+'_diff_'+refParts[j][1]+'.png',partDiff)
				else:
					printFlushExit("can't compare pictures with different sizes", 1)
				j += 1
			if(rowHadError):
				print("Error "+rowMsg)
			elif(rowMsg != ''):
				print(rowMsg)
		else:
			printFlushExit("can't compare pictures with different sizes", 1)
	
	print('                      ------------------------------------------')
	print('maximum error found   '+str(round(maxEncountered['view'],6)).ljust(11)+str(round(maxEncountered['term'],6)).ljust(11)+str(round(maxEncountered['cont'],6)).ljust(11)+str(round(maxEncountered['insp'],6)).ljust(11))
	
	if(hadError):
		printFlushExit("\033[91m Errors detected! The picture differences are saved in the artifacts.\033[0m", 1)
except:
	print("Exception encountered during tests")
	os._exit(1)

