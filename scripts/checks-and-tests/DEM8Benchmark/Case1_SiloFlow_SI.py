#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 

# Benchmark of basic performance of open-source DEM simulation systems
# Case 1: Silo flow

# Units: SI (m, N, Pa, kg, sec)

# -------------------------------------------------------------------- #
# Input Data -> Define Material and Orifice size. Uncomment the prefered choice. The rest of the script should be automated for all different scenarios.
granularMaterial='M1'
#granularMaterial='M2'

#fileName='LargeOrifice'; z=55.4222/1000.	# "z" corresponds to the lowest point of the funnel, in the coordinate system proposed in the latest update of the data.
fileName='SmallOrifice'; z=59.3008/1000.

# -------------------------------------------------------------------- #
# MPI initialization
numMPIThreads=1

if numMPIThreads > 1:
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

if granularMaterial=='M1':
	M1=O.materials.append(FrictMat(young=1.0e9,poisson=0.2,density=2500,label='M1'))
	e_gg=e_M1_M1	# Coefficient of restitution (e) between granular material (g) and granular material (g).
	f_gg=f_M1_M1	# Coefficient of friction (f)...

	e_gs=e_M1_St	# Coefficient of restitution (e) between granular material (g) and steel (s).
	f_gs=f_M1_St	# Coefficient of friction (f)...

#	O.dt=1.5e-6	# FIXME: This is the suggested timestep value by the organizers for M1.

elif granularMaterial=='M2':
	M2=O.materials.append(FrictMat(young=0.5e9,poisson=0.2,density=2000,label='M2'))
	e_gg=e_M2_M2
	f_gg=f_M2_M2

	e_gs=e_M2_St
	f_gs=f_M2_St

#	O.dt=2e-6	# FIXME: This is the suggested timestep value by the organizers for M2.

F_gg=atan(f_gg) # Friction Angle between granular material (g) and granular material (g).
F_gs=atan(f_gs) # Friction Angle between granular material (g) and steel (s).


# -------------------------------------------------------------------- #
# Load particles from .txt file
from yade import ymport
sp=ymport.text('inputData/Case1_SiloFlow_PartCoordinates_'+fileName+'.txt',material=granularMaterial)
#sp=sp[0:2000] #Remove this; was used for debbuging purposes

# -------------------------------------------------------------------- #
## Engines 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],label="collider"),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(
			frictAngle = MatchMaker(matches=((1,1,F_gg),(0,1,F_gs))), # 0 being the id of Steel and
			en         = MatchMaker(matches=((1,1,e_gg),(0,1,e_gs)))  # 1 being the id of granularMaterial.
		)],
		[Law2_ScGeom_MindlinPhys_Mindlin()],
	),
	NewtonIntegrator(damping=0,gravity=[0,0,-9.81],label="newton"),
	#GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()), #FIXME Remember to reinstate parallelMode=True when we use MPI
	#VTKRecorder(virtPeriod=0.1,fileName='/tmp/Case1_SiloFlow-',recorders=['spheres','facets']),
]


# -------------------------------------------------------------------- #
# Load facets making the cylinder
facets = ymport.textFacets('inputData/Case1_SiloFlow_Walls_'+fileName+'.txt',color=(0,1,0),material=Steel)
fctIds = range(len(facets))

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
    O.bodies.append(sp)


collider.verletDist = 0.5e-3
O.dt=1e-10
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

# -------------------------------------------------------------------- #
# Erase particles flowing out of the silo
def eraseEscapedParticles():
	for i in indSpheres:
		if O.bodies[i].state.pos[2]<-z: # Delete particles when they pass the orifice #FIXME: Maybe -z-radius instead?
			indSpheres.remove(i)
			O.bodies.erase(i)

# -------------------------------------------------------------------- #
# Record time-dependent number of retained particles
indSpheres=[sp[i].id for i in range(0,len(sp))] #Indices of spherical particles

from yade import plot
def addPlotData(): 
	retained=sum(1 for i in indSpheres if O.bodies[i].state.pos[2]>-z)  
	eraseEscapedParticles()
	plot.addData(retained=retained, time1=O.time)

addPlotData() # I use this to record the initial state for O.time=0.0 FIXME: Since we already run NSTEPS, O.time is not exactly zero.
O.engines=O.engines+[PyRunner(virtPeriod=0.1,command='addPlotData()')] # Here I use virtPeriod=0.1, following the provided .xlsx example file.

plot.plots={'time1':'retained'}
plot.plot(noShow=False)

# -------------------------------------------------------------------- #
# GUI
if opts.nogui==False:
	from yade import qt
	v=qt.View()

	v.eyePosition = Vector3(0,-0.6,0.1)
	v.upVector    = Vector3(0,0,1)
	v.viewDir     = Vector3(0,1,0)
#	v.grid=(False,True,False)
	v.ortho       = True

	rndr=yade.qt.Renderer()
	#rndr.shape=False
	#rndr.bound=True

O.run()
plot.saveDataTxt('Case1_SiloFlow_'+fileName+'_'+granularMaterial+'.txt',vars=('time1','retained'))

