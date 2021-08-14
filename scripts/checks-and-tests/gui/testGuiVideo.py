# -*- encoding=utf-8 -*-

from __future__ import print_function
import os, sys
from testGuiHelper import TestGUIHelper
TriaxialTest(noFiles=True).load()

#yade.qt.View() # cannot open 3D View. Must allow testGuiHelper do this.

scr = TestGUIHelper("Video",True)
guiIterPeriod = 20

from yade import qt

openViewAt       = 61 # must be synchronized with clicking
snapEveryNthIter = 5  # save image every 5 iterations
snapshotFiles    = O.tmpFilename()

O.engines=O.engines+[
	qt.SnapshotEngine(fileBase=snapshotFiles, label='snapshotter', iterPeriod=snapEveryNthIter, firstIterRun=openViewAt, ignoreErrors=False),
	PyRunner(iterPeriod=snapEveryNthIter, firstIterRun=openViewAt, command='checkFile()'),
# Skip mencoder creating a final video. Let us only test that SnapshotEngine works in the GUI.
# This can be improved later.
#	PyRunner(iterPeriod=500,command='finito()')
	PyRunner(iterPeriod=guiIterPeriod, command='scr.screenshotEngine()')
]
rr=qt.Renderer()
rr.shape,rr.intrPhys=False,True

def checkFile():
	# Let's check if snapshotter save the image
	num = (O.iter - openViewAt)//snapEveryNthIter
	fname = snapshotFiles + str(num).rjust(5,'0') + ".png"
	if(not os.path.isfile(fname)):
		# FIXME: throwing exceptions inside GUI test just causes a timeout in the pipeline CI. So we could improve this.
		msg = "SnapshotEngine is not saving files! Tried to find this file: "+fname
		sys.__stderr__.write(msg)
		sys.__stderr__.flush()
		os._exit(1)
	else:
		pass
		#sys.__stderr__.write("File OK : " + fname)
		#sys.__stderr__.flush()

def finito():
	"""This function will be called after 500 steps. Since SnapshotEngine waits for a new 3d view to open,
	it must run after the script has finished and the command line appears
	(see old bug report on old site: https://bugs.launchpad.net/yade/+bug/622669).

	For that reason, O.run() is at the end of the script and this function will be called
	once we want to exit really.
	"""
	# skip mencoder in the test. Let's just test that SnapshotEngine works in the GUI
	# This can be improved later.
	makeVideo(snapshotter.snapshots,out='/tmp/video.avi')
	print("Video saved in /tmp/video.avi")
	import sys
	os._exit(0)

O.run(guiIterPeriod * scr.getTestNum() + 1)

