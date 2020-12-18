# encoding: utf-8
# 2011 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
from __future__ import print_function
from past.builtins import execfile
import yade,math,os,sys
import time, datetime

class YadeCheckError(Exception):
	"""Use ``raise YadeCheckError("message") when check fails. Or throw any other python exception that seems suitable."""
	pass

scriptsToRun=os.listdir(checksPath)
failedScripts=list()
maxElapsedTime=0

#some scripts are singleCore only because of required 100% reproducibility
singleCore= ['checkVTKRecorder.py' ,'checkPotentialVTKRecorders.py' ,'checkJCFpm.py' ,'checkColliderCorrectness.py' ,'checkColliderConstantness.py', 'checkThermalEngine.py']

#checkSpawn.py fails always for now, needs investigations
skipScripts = ['checkList.py','checkSpawn.py']
onlyOneScript = [] # use this if you want to test only one script, it takes precedence over skipScripts.

def multiCore(): # multi core --check is running.
	return ((opts.threads != None and opts.threads != 1) or (opts.cores != None and opts.cores != '1'))

def mustCheck(sc):
	if(len(onlyOneScript)==1): return sc in onlyOneScript
	if(multiCore()):
		return (sc not in singleCore) and (sc not in skipScripts)
	return sc not in skipScripts

for script in scriptsToRun:
	if (script[len(script)-3:]==".py" and mustCheck(script)):
		print("###################################")
		print("running: ",script)
		try:
			t0 = time.time()
			execfile(checksPath+"/"+script)
			t1 = time.time()
			elapsedTime = t1-t0
			maxElapsedTime = max(elapsedTime,maxElapsedTime)
			print("Status:\033[92m success\033[0m, time spent on this check:"+("\033[92m " if elapsedTime<30 else "\033[91m ")+str(datetime.timedelta(seconds=elapsedTime))+"\033[0m")
			print("___________________________________")
		except Exception as e:
			failedScripts.append(script)
			print('\033[91m',script," failure, caught exception ",e.__class__.__name__,": ",e,'\033[0m')
		O.reset()
	elif (not mustCheck(script)):
		print("###################################")
		if(multiCore()):
			print("Skipping %s, because this script only works with -j1 (single core) or is in skipScripts"%script)
		else:
			print("Skipping %s, because it is in skipScripts"%script)

if(maxElapsedTime > 30):
	print("\033[95mWARNING: some checks took longer than 30 seconds.\033[0m")
print("Most time spend on a single check: "+str(datetime.timedelta(seconds=maxElapsedTime)))

if (len(failedScripts)!=0):
	print('\033[91m', len(failedScripts) , " tests are failed"+'\033[0m')
	for s in failedScripts: print("  "+s)
	sys.exit(1)
else:
	# https://misc.flogisoft.com/bash/tip_colors_and_formatting
	print('\033[92m'+"*** ALL CHECKS PASSED ***"+'\033[0m')
	sys.exit(0)

