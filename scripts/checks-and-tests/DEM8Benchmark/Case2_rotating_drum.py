#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 
# 2020 © Robert Caulk <rob.caulk@gmail.com>

# Benchmark of basic performance of open-source DEM simulation systems
# Case 2: Rotating Drum

# Units: SI (m, N, Pa, kg, sec)

angularVelocity = 2

# -------------------------------------------------------------------- #
# MPI initialization
numMPIThreads=0

if numMPIThreads > 1:
	print('in mpi?')
	from yade import mpy as mp
	mp.initialize(numMPIThreads)

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


#	O.dt=1.5e-6	# FIXME: This is the suggested timestep value by the organizers for M1.

M2=O.materials.append(FrictMat(young=0.5e9,poisson=0.2,density=2000,label='M2'))


#	O.dt=2e-6	# FIXME: This is the suggested timestep value by the organizers for M2.

F_g1g2=atan(f_M1_M2)
F_g1g1=atan(f_M1_M1)
F_g2g2=atan(f_M2_M2) # Friction Angle between granular material (g) and granular material (g).
F_gs=atan(f_M1_St) # Friction Angle between granular material (g) and steel (s).Same for both materials


# -------------------------------------------------------------------- #
# Load particles from .txt file
from yade import ymport
import os

# FIXME: Download coordinate files like this
#os.system('wget http://perso.3sr-grenoble.fr/users/bchareyre/yade/input/'+fileName+'.stl')

sp_m2=ymport.text('inputData/Case2_Drum_PartCoordinates_M2.txt',material=M2,color=(0,0,1))
sp_m1=ymport.text('inputData/Case2_Drum_PartCoordinates_M1.txt',material=M1,color=(1,0,0))

# Load facets making the cylinder
facets = ymport.textFacets('inputData/Case2_Drum_Walls.txt',color=(0,1,0),material=Steel)
drum_ids = range(len(facets))


# -------------------------------------------------------------------- #
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
	RotationEngine(rotateAroundZero=True,zeroPoint=(0,0,0),rotationAxis=(0,1,0),angularVelocity=angularVelocity,ids=drum_ids),
	#GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()) #FIXME Remember to reinstate parallelMode=True when we use MPI
	#VTKRecorder(virtPeriod=0.01,fileName='tmp/Case2_drum-',recorders=['spheres','facets']),
]


# -------------------------------------------------------------------- #

NSTEPS = 1000

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
			mp.mprint("s",s,"startId",nextId)
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


#collider.verletDist = 0.5e-3
O.dt=2e-6
#O.dt=0.8*PWaveTimeStep() #FIXME: Consider larger timestep instead?
O.dynDt=False

O.stopAtTime=5 #5 sec. This is the timeframe of interest, proposed by the organizers in the latest update of the data.

if numMPIThreads>1:
	mp.mpirun(1,numMPIThreads,False) #this is to eliminate initialization overhead in Cundall number and timings.
	mp.YADE_TIMING=True
	t1=time.time()
	mp.mpirun(NSTEPS,withMerge=False)
	t2=time.time()
	mp.mprint("num. bodies:",len([b for b in O.bodies])," ",len(O.bodies))
	if mp.rank==0:
		  mp.mprint("CPU wall time for ",NSTEPS," iterations:",t2-t1,"; Cundall number = TODO")
	#mp.mergeScene()

else: 
	O.run(1,True)
	t1=time.time()
	O.run(NSTEPS,True)
	t2=time.time()
	print("num. bodies:",len([b for b in O.bodies])," ",len(O.bodies))
	print("CPU wall time for ",NSTEPS," iterations:",t2-t1,"; Cundall number = TODO")

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
def addPlotData(): 
	zone1_M1_count,zone1_M2_count,zone2_M1_count,zone2_M2_count = getNumParticlesInQuadrants()
	plot.addData(zone1_M1_count=zone1_M1_count,
						zone1_M2_count=zone1_M2_count,
						zone2_M1_count=zone2_M1_count,
						zone2_M2_count=zone2_M2_count, 
						time1=O.time)
	plot.saveDataTxt('Case2_rotating_drum_data.txt',vars=('time1','zone1_M1_count','zone1_M2_count','zone2_M1_count','zone2_M2_count'))

addPlotData() # I use this to record the initial state for O.time=0.0 FIXME: Since we already run NSTEPS, O.time is not exactly zero.
O.engines=O.engines+[PyRunner(virtPeriod=0.01,command='addPlotData()')] # Here I use virtPeriod=0.1, following the provided .xlsx example file.

def stopAt():
	if O.time>=5: O.stop()

O.engines=O.engines+[PyRunner(iterPeriod=1000,command='stopAt()')]

plot.plots={'time1':('zone1_M1_count','zone1_M2_count','zone2_M1_count','zone2_M2_count')}
plot.plot(noShow=False)

# -------------------------------------------------------------------- #
# GUI
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

O.run()
