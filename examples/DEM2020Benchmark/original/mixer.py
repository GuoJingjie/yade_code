#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 
# 2020 © Robert Caulk <rob.caulk@gmail.com>

# Benchmark of basic performance of open-source DEM simulation systems
# Case 2: Rotating Drum

# Units: SI (m, N, Pa, kg, sec)

#####################   1. INPUT/OUTPUT  #####################


angularVelocity = 2

# 5 sec. This is the timeframe of interest, proposed by the organizers 
# pass additional argument to the script to change it: 'yadedaily -n -x Case2_rotating_drum.py 0.0001'

simulatedTime = 5 if len(sys.argv)<=1 else float(sys.argv[1])

# -------------------------------------------------------------------- #
# MPI initialization
numMPIThreads=1

if numMPIThreads > 1:
	from yade import mpy as mp
	mp.initialize(numMPIThreads)
	
try:
    os.mkdir('inputData')
except:
    pass # will pass if folders already exist

try:
    os.mkdir('outputData')
except:
    pass # will pass if folders already exist

# -------------------------------------------------------------------- #
# Materials
Steel = O.materials.append(FrictMat(young=210e9,poisson=0.2,density=7200,label='Steel'))

# Coeff of restitution (e) / Coeff of friction (f)
e_M1_M2=0.45;	f_M1_M2=0.2
e_M1_M1=0.5;	f_M1_M1=0.3
e_M1_St=0.4;	f_M1_St=0.2
e_M2_M2=0.4;	f_M2_M2=0.4
e_M2_St=0.4;	f_M2_St=0.2


M1=O.materials.append(FrictMat(young=1.0e9,poisson=0.2,density=2500,label='M1'))

M2=O.materials.append(FrictMat(young=0.5e9,poisson=0.2,density=2000,label='M2'))

F_g1g2=atan(f_M1_M2)
F_g1g1=atan(f_M1_M1)
F_g2g2=atan(f_M2_M2) # Friction Angle between granular material (g) and granular material (g).
F_gs=atan(f_M1_St) # Friction Angle between granular material (g) and steel (s).Same for both materials


# -------------------------------------------------------------------- #
# Load particles from .txt file

from yade import ymport
wallFile='Case2_Drum_Walls.txt'
spheres_M1='Case2_Drum_PartCoordinates_M1.txt'
spheres_M2='Case2_Drum_PartCoordinates_M2.txt'

hasInputSpheres = os.path.exists('inputData/'+spheres_M1)
if not hasInputSpheres:
	print("Downloading sphere file",spheres_M1)
	try:
		os.system('wget -nc -O inputData/'+spheres_M1+' http://yade-dem.org/publi/data/DEM8/'+spheres_M1)
		os.system('wget -nc -O inputData/'+spheres_M2+' http://yade-dem.org/publi/data/DEM8/'+spheres_M2)
	except:
		print("** probably no internet connection, grab",spheres_M1,"by yourself **")


hasInputWall = os.path.exists('inputData/'+wallFile)
if not hasInputWall:
	print("Downloading mesh file",wallFile)
	try:
		os.system('wget -nc -O inputData/'+wallFile+' http://yade-dem.org/publi/data/DEM8/'+wallFile)
	except:
		print("** probably no internet connection, grab",wallFile,"by yourself **")


#####################   2. YADE PART  #####################

sp_m1=ymport.text('inputData/'+spheres_M1,material=M1,color=(1,0,0))
sp_m2=ymport.text('inputData/'+spheres_M2,material=M2,color=(0,0,1))
facets = ymport.textFacets('inputData/'+wallFile,color=(0,1,0),material=Steel)
#facets = ymport.stl(fileName+'.stl',color=(0,1,0),material=Steel)
drum_ids = range(len(facets))

## Engines 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],label="collider"),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(
			frictAngle = MatchMaker(matches=((M1,M1,F_g1g1),(Steel,M1,F_gs),(Steel,M2,F_gs),(M1,M2,F_g1g2),(M2,M2,F_g2g2))), 
			en         = MatchMaker(matches=((M1,M1,e_M1_M1),(Steel,M1,e_M1_St),(Steel,M2,e_M2_St),(M1,M2,e_M1_M2),(M2,M2,e_M2_M2)))  
		)],
		[Law2_ScGeom_MindlinPhys_Mindlin()],
	),
	NewtonIntegrator(damping=0,gravity=[0,0,-9.81],label="newton"),
	RotationEngine(dead=1,rotateAroundZero=True,zeroPoint=(0,0,0),rotationAxis=(0,1,0),angularVelocity=angularVelocity,ids=drum_ids,label='rotation'),
	#GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()) #FIXME Remember to reinstate parallelMode=True when we use MPI
	#VTKRecorder(virtPeriod=0.01,fileName='tmp/Case2_drum-',recorders=['spheres','facets']),
]


# -------------------------------------------------------------------- #

NSTEPS = 1000

sp=sp_m1+sp_m2

if numMPIThreads > 1:
	mp.mprint("appending bodies, rank", mp.rank)
	if mp.rank==0:
			O.bodies.append(facets)
			mp.mprint("master has",len(O.bodies), "facets")
	else:
			import numpy as np
			layers = np.array_split(sp,numMPIThreads-1)
			mp.mprint("layers",[len(l) for l in layers])
			layerNo = mp.rank-1 #rank zero is for facets
			nextId = int(len(facets) + np.sum([len(x) for x in layers[:layerNo]]))
			for s in layers[layerNo]:
			    s.subdomain = mp.rank
			    O.bodies.insertAtId(s,nextId)
			    nextId += 1
			mp.mprint("s",s,"lastId",nextId-1)
			
			# tune mpi
	mp.VERBOSE_OUTPUT=False
	mp.DISTRIBUTED_INSERT=True
	mp.REALLOCATE_FREQUENCY=4
	mp.ACCUMULATE_FORCES=False
	mp.MAX_RANK_OUTPUT=4
    

else:
	O.bodies.append(facets)
	O.bodies.append(sp_m1)
	O.bodies.append(sp_m2)

collider.verletDist = 2e-4 # 10% of smallest radius
O.dt=8e-7
#O.dt=0.8*PWaveTimeStep() 
O.dynDt=False

# let the unbalanced force settle before?
while 0:
	O.run(1000, True)
	unb=unbalancedForce()
	print('settling particles unb',unb)
	if unb<0.01: break

rotation.dead=0
O.resetTime()
startTime = time.time()

if numMPIThreads>1:
	mp.mpirun(1,numMPIThreads,False) #this is to eliminate initialization overhead in Cundall number and timings.
	mp.YADE_TIMING=True
	t1=time.time()
	mp.mpirun(NSTEPS,withMerge=False)
	t2=time.time()
	mp.mprint("num. bodies:",len([b for b in O.bodies])," ",len(O.bodies))
	if mp.rank==0:
		  mp.mprint("CPU wall time for ",NSTEPS," iterations:",t2-t1,"; Cundall number =",NSTEPS/(t2-t1))
	#mp.mergeScene()

else: 
	O.run(1,True)
	t1=time.time()
	O.run(NSTEPS,True)
	t2=time.time()
	print("num. bodies:",len([b for b in O.bodies])," ",len(O.bodies))
	print("CPU wall time for ",NSTEPS," iterations:",t2-t1,"; Cundall number =",len(sp_m1+sp_m2)*NSTEPS/(t2-t1))

## Check particle quadrant
def getNumParticlesInQuadrants():
	zone1_M1_count = 0
	zone2_M1_count = 0
	zone1_M2_count = 0
	zone2_M2_count = 0
	for b in O.bodies:
		if isinstance(b.shape, Sphere):
			pos = b.state.pos
			if pos[0] > 0 and pos[2]<0: 
				if b.material.label=='M1': zone1_M1_count+=1 # +x-z in yade
				if b.material.label=='M2': zone1_M2_count+=1
			if pos[0] < 0 and pos[2]>0: 
				if b.material.label=='M1': zone2_M1_count+=1 # -x+z in yade
				if b.material.label=='M2': zone2_M2_count+=1
				
	return zone1_M1_count,zone1_M2_count,zone2_M1_count,zone2_M2_count

from yade import plot
plot.plots={'time1':('zone1_M1_count','zone1_M2_count','zone2_M1_count','zone2_M2_count')}

def addPlotData(save=True): 
	zone1_M1_count,zone1_M2_count,zone2_M1_count,zone2_M2_count = getNumParticlesInQuadrants()
	plot.addData(zone1_M1_count=zone1_M1_count,
						zone1_M2_count=zone1_M2_count,
						zone2_M1_count=zone2_M1_count,
						zone2_M2_count=zone2_M2_count, 
						time1=O.time)
	if save:
		plot.saveDataTxt('outputData/Case2_rotating_drum.txt',vars=('time1','zone1_M1_count','zone1_M2_count','zone2_M1_count','zone2_M2_count'))
		plot.plot(noShow=True).savefig('outputData/Case2_rotating_drum.png')

print('particles settled and ready to rotate')

O.stopAtTime=simulatedTime

addPlotData(False) # I use this to record the initial state for O.time=0.0 FIXME: Since we already run NSTEPS, O.time is not exactly zero.
O.engines=O.engines+[PyRunner(virtPeriod=0.1,command='addPlotData()')] # Here I use virtPeriod=0.1, following the provided .xlsx example file.



#####################   3. GUI and timings  #####################

if opts.nogui==False:
	from yade import qt
	v=qt.View()
	v.eyePosition = Vector3(0,-0.4,0)
	#v.upVector    = Vector3(0,0,1)
	v.viewDir     = Vector3(0,1,0)
#	v.grid=(False,True,False)
#	v.ortho       = True
	rndr=yade.qt.Renderer()
	#rndr.shape=False
	#rndr.bound=True
print("stop at ",O.stopAtTime," current=",O.time)

O.run(-1,wait = opts.nogui)

wallTime = time.time() - startTime
f = open("timings.txt","a")
f.write('Case2_rotating_drum '+str(O.time)+' '+str(wallTime)+'\n')
