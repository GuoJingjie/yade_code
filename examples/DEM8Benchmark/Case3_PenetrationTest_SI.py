#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 

# Benchmark of basic performance of open-source DEM simulation systems
# Case 3: Penetration test

# Units: SI (m, N, Pa, kg, sec)

# -------------------------------------------------------------------- #
# Input Data -> Define Number of particles. Uncomment the prefered choice. The rest of the script should be automated for all different scenarios.
N=25000
#N=50000
#N=100000

# -------------------------------------------------------------------- #
## MPI initialization --> FIXME: I leave the MPI magic to Bruno :)
#numMPIThreads=1

#if numMPIThreads > 1:
#    from yade import mpy as mp
#    mp.initialize(numMPIThreads)

# -------------------------------------------------------------------- #
# Materials
Steel=O.materials.append(FrictMat(young=210e9,poisson=0.2,density=7200,frictionAngle=atan(0.2),label='Steel'))

# -------------------------------------------------------------------- #
# Assign coeff of restitution (e)
e_M1_M1=0.5;
e_M1_St=0.4;

M1=O.materials.append(FrictMat(young=1.0e9,poisson=0.2,density=2500,frictionAngle=atan(0.3),label='M1'))
e_gg=e_M1_M1	# Coefficient of restitution (e) between granular material (g) and granular material (g)
e_gs=e_M1_St	# Coefficient of restitution (e) between granular material (g) and steel (s)
#e_ss=0.4 	# Used only for debugging purposes in the Ip2 functor; not needed

# -------------------------------------------------------------------- #
# Load the initial sphere pack from .txt file
from yade import ymport
sp=ymport.text('inputData/Case3_PenetrationTest_'+str(N)+'_Particles.txt',material=M1)

sp[-1]=sphere(sp[-1].state.pos,sp[-1].shape.radius,material='Steel') # Redefine this with the correct material, as sp[-1] corresponds to the large sphere made of steel


# FIXME: [...] MPI stuff :)

O.bodies.append(sp)

# -------------------------------------------------------------------- #
# Load facets making the box
from yade import ymport
facets = ymport.textFacets('inputData/Case3_PenetrationTest_Walls.txt',color=(0,1,0),material=Steel)
fctIds = range(len(facets))

# FIXME: [...] MPI stuff :)

O.bodies.append(facets)

# -------------------------------------------------------------------- #
# Timestep
#O.dt=1e-10
O.dt=0.8*PWaveTimeStep() # Alternatively, we could use this increased timestep or consider using the GST.

# -------------------------------------------------------------------- #
## Engines 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(
			en         = MatchMaker(matches=((1,1,e_gg),(0,1,e_gs)))  # 0 being the id of Steel and 1 being the id of M1 # ,(0,0,e_ss) # e_ss is not needed, as the steel ball is not supposed to touch the steel box. Was defined for debugging purposes above
		)],
		[Law2_ScGeom_MindlinPhys_Mindlin()],
	),
	NewtonIntegrator(damping=0,gravity=[0,0,-9.81]),
#	GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()), #FIXME Remember to reinstate parallelMode=True when we use MPI 
	#VTKRecorder(virtPeriod=5e-4,fileName='/tmp/Case3_PenetrationTest-',recorders=['spheres','facets']),
]

# -------------------------------------------------------------------- #
# Make sure O.bodies[N] is the large sphere (bSphere)	# Maybe this is an overkill, since we have the initial packing. Feel free to remove this check and assign velocity directly to O.bodies[N] in all cases.
if O.bodies[N].shape.radius==0.01:
	bSphere=O.bodies[N]
	bSphere.state.vel=[0,0,-20]
else:
	for b in sp:
		if b.shape.radius==0.01:
			bSphere=b
			b.state.vel=[0,0,-20]

# -------------------------------------------------------------------- #
# Record time-dependent position of the large sphere (bSphere)
from yade import plot
def addPlotData(): 
	z=bSphere.state.pos[2]
	plot.addData(z=z, time1=O.time)

addPlotData() # I use this to record the initial state for O.time=0.0
O.engines=O.engines+[PyRunner(virtPeriod=5e-4,command='addPlotData()')] # Here I use virtPeriod=0.0005, following the provided .xlsx example file

plot.plots={'time1':'z'}
plot.plot(noShow=False)


# FIXME: For N=25000, the initial position of the large sphere is z=+0.04m instead of z=+0.06m which is shown in the provided .xlsx example file. I think we need to ask for clarifications on this.

# -------------------------------------------------------------------- #
# GUI
if opts.nogui==False:
	from yade import qt
	v=qt.View()

	v.eyePosition=Vector3(0,-0.35,0)
	v.upVector    = Vector3(0,0,1)
	v.viewDir     = Vector3(0,1,0)
#	v.grid        = (False,True,False)
	v.ortho       = True

	rndr=yade.qt.Renderer()
	#rndr.shape=False
	#rndr.bound=True

if N==25000 or N==50000:
	O.stopAtTime=0.04 # These values were taken from the .xlsx example file
elif N==100000:
	O.stopAtTime=0.05

O.run()
plot.saveDataTxt('Case3_PenetrationTest_'+str(N)+'.txt',vars=('time1','z'))

