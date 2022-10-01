# -*- encoding=utf-8 -*-
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>

# Benchmark of basic performance of open-source DEM simulation systems
# Case 1: Silo flow
# Units: SI (m, N, Pa, kg, sec)

# script execution:
# yade Case1SiloFlow.py   (with default arguments, else:)
# yade Case1SiloFlow.py  [small|large] [M1|M2]
# example: yade Case1SiloFlow.py  small M1 (equivalent to default)
#
# third optional argument is simulated time (default 5s), example:
# yade Case1SiloFlow.py small M1 0.3

numThreads = 1
reportTiming = False

# Configure MPI module if needed
mpi = 'MPI' in yade.config.features
if mpi:
	from yade import mpy as mp
else:
	print("yade is compiled without MPI support, numThreads>1 ignored")
	numThreads = 1

#####################   1. INPUT/OUTPUT  #####################

try:
	os.mkdir('inputData')
except:
	pass  # will pass if folders already exist

try:
	os.mkdir('outputData')
except:
	pass  # will pass if folders already exist

# -------------------------------------------------------------------- #
# Input Data -> Define Material and Orifice size.

size = str(sys.argv[1]) if len(sys.argv) > 1 else 'small'
material = str(sys.argv[2]) if len(sys.argv) > 2 else 'M1'
simulationTime = float(sys.argv[3]) if len(sys.argv) > 3 else 5

if size == 'large':
	fileName = 'Case1_large'
	z = 55.4222 / 1000.  # This is the height of the lowest point of the funnel (at the orifice), measuring from the lowest cylindrical cross section of the silo
else:
	z = 59.3008 / 1000.
	fileName = 'Case1_small'

# -------------------------------------------------------------------- #
# Materials
Steel = O.materials.append(FrictMat(young=210e9, poisson=0.2, density=7200, label='Steel'))

# -------------------------------------------------------------------- #
# Asign values based on the Material and Orifice size
# Coeff of restitution (e) / Coeff of friction (f)
e_M1_M2 = 0.45
f_M1_M2 = 0.2
e_M1_M1 = 0.5
f_M1_M1 = 0.3
e_M1_St = 0.4
f_M1_St = 0.2
e_M2_M2 = 0.4
f_M2_M2 = 0.4
e_M2_St = 0.4
f_M2_St = 0.2

if material == 'M1':
	M1 = O.materials.append(FrictMat(young=1.0e9, poisson=0.2, density=2500, label='M1'))
	e_gg = e_M1_M1  # Coefficient of restitution (e) between granular material (g) and granular material (g)
	f_gg = f_M1_M1  # Coefficient of friction (f)...
	e_gs = e_M1_St  # Coefficient of restitution (e) between granular material (g) and steel (s)
	f_gs = f_M1_St  # Coefficient of friction (f)...
elif material == 'M2':
	M2 = O.materials.append(FrictMat(young=0.5e9, poisson=0.2, density=2000, label='M2'))
	e_gg = e_M2_M2
	f_gg = f_M2_M2
	e_gs = e_M2_St
	f_gs = f_M2_St

F_gg = atan(f_gg)  # Friction Angle between granular material (g) and granular material (g)
F_gs = atan(f_gs)  # Friction Angle between granular material (g) and steel (s)

# urls where data should be obtained
urls = {}  # case name will be size+material
# FIXME: particle coordinates for M1 and M2 are the same right now, the original data needs to be retrieved
urls["smallM1"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_PartCoordinates_SmallOrifice.txt"
urls["smallM2"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_PartCoordinates_SmallOrifice.txt"
urls["largeM1"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_PartCoordinates_LargeOrifice.txt"
urls["largeM2"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_PartCoordinates_LargeOrifice.txt"
urls["small"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_Walls_SmallOrifice.txt"
urls["large"] = "https://yade-dem.org/publi/data/DEM8/Case1_SiloFlow_Walls_LargeOrifice.txt"

# This condition is not abolutely necessary but it would be inelegant to
# download *.stl and generate densePack N times when we need it done only on master (centralized scene method)
if not mpi or mp.rank == 0:
	from yade import ymport
	wallFile = 'inputData/Case1_SiloFlow_Walls_' + fileName + '.txt'
	sphereFile = 'inputData/Case1_SiloFlow_PartCoordinates_' + fileName + '.txt'

	hasInputSpheres = os.path.exists(sphereFile)
	if not hasInputSpheres:
		print("Downloading sphere file", sphereFile)
		try:
			print('wget -nc -O ' + sphereFile + '_temp ' + urls[size + material])
			os.system('wget -nc -O ' + sphereFile + '_temp ' + urls[size + material])
		except:
			print("** probably no internet connection, grab", sphereFile, "by yourself **")

		with open(sphereFile + '_temp') as x, open(sphereFile, 'w') as y:
			for line in x:
				columns = line.split()
				if len(columns) < 3:
					continue  # trailing empty lines
				y.write('\t'.join(columns) + ' 0.002\n')
		y.close()

	hasInputWall = os.path.exists(wallFile)
	if not hasInputWall:
		print("Downloading mesh file", wallFile)
		try:
			os.system('wget -nc -O ' + wallFile + ' ' + urls[size])
		except:
			print("** probably no internet connection, grab", wallFile, "by yourself **")

#####################   2. YADE PART  #####################

# -------------------------------------------------------------------- #
## Engines
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()], label="collider"),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(hertzian=True)],
                [
                        Ip2_FrictMat_FrictMat_MindlinPhys(
                                frictAngle=MatchMaker(matches=((1, 1, F_gg), (0, 1, F_gs))),  # 0 being the id of Steel and
                                en=MatchMaker(matches=((1, 1, e_gg), (0, 1, e_gs)))  # 1 being the id of material
                        )
                ],
                [Law2_ScGeom_MindlinPhys_Mindlin(preventGranularRatcheting=False)],
        ),
        NewtonIntegrator(damping=0, gravity=[0, 0, -9.810], label="newton"),
        #GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()), #FIXME Remember to reinstate parallelMode=True when we use MPI
        #VTKRecorder(virtPeriod=0.04,fileName='/tmp/Silo-',recorders=['spheres','facets']),
]

if not mpi or mp.rank == 0:
	sp = ymport.text(sphereFile, material=material)
	facets = ymport.textFacets(wallFile, color=(0, 1, 0), material=Steel)
	#facets = ymport.stl(fileName+'.stl',color=(0,1,0),material=Steel)
	fctIds = range(len(facets))

	O.bodies.append(facets)
	O.bodies.append(sp)

	# -------------------------------------------------------------------- #
	# Count the number of spherical particles to verify sample size. We can comment this out later on.
	numSpheres = 0
	for b in O.bodies:
		if isinstance(b.shape, Sphere):
			numSpheres = numSpheres + 1
	print('The total number of spheres is: ', numSpheres)

	collider.verletDist = 0.1 * O.bodies[-1].shape.radius
	O.dynDt = False
	O.dt = 1.5e-6 if material == 'M1' else 2e-6

# -------------------------------------------------------------------- #
# Erase particles flowing out of the silo


def eraseEscapedParticles():
	global numErased
	count = 0
	ts = time.time()
	ers = []
	for b in O.bodies:
		if isinstance(
		        b.shape, Sphere
		) and b.state.pos[2] < -z:  # I do not delete the particles right after they pass the orifice, to disturb the simulation as little as possible
			ers.append(b.id)
			count += 1
	if mpi:
		mp.bodyErase(ers)
	else:
		for b in ers:
			O.bodies.erase(b.id)
	numErased += count


#-------------------------------------------------------------------- #
#Record time-dependent number of retained particles and vtk export

from yade import plot
plot.plots = {'time': (('retained', 'bo--'), None, ('Cu', "kx--"))}

numErased = 0


def addPlotData(Cu):
	plot.addData(time=O.time, retained=numSpheres - numErased, Cu=Cu)


addPlotData(0)

from yade import export
vtk = export.VTKExporter("spheresFinal")

# -------------------------------------------------------------------- #
# Run iterations

if mpi:  # import and tune MPI module
	mp.DOMAIN_DECOMPOSITION = True
	mp.ACCUMULATE_FORCES = False
	mp.MERGE_W_INTERACTIONS = False
	mp.REALLOCATE_MINIMAL = False
	mp.REALLOCATE_FREQUENCY = 20
	mp.USE_CPP_REALLOC = True
	mp.MINIMAL_INTERSECTIONS = True
	mp.YADE_TIMING = reportTiming
	mp.USE_CPP_INTERS = True
	#mp.VERBOSE_OUTPUT=True
else:
	O.timingEnabled = reportTiming

startTime = time.time()
substeps = 500
while len(O.bodies) - numErased > 0 and O.time < simulationTime:
	#for k in range(4):
	t1 = time.time()
	if mpi:
		mp.mpirun(substeps, numThreads, withMerge=True)  # if numThreads=1 this will fall-back to normal O.run() and mp.rank=0
	else:
		O.run(substeps, True)
	eraseEscapedParticles()
	if mpi and mp.rank > 0:
		continue  # mpi workers do not record
	t2 = time.time()
	addPlotData((numSpheres - numErased) * substeps / (t2 - t1))
	#vtk.exportSpheres(what=dict(particleVelocity='b.state.vel',domain='b.subdomain'))
	plot.plot(noShow=True).savefig('outputData/' + fileName + '_' + material + '_np' + str(numThreads) + '.png')
	plot.saveDataTxt('outputData/' + fileName + '_' + material + '.txt')
	print("iter=", O.iter, ", last substep erased", numErased, "in", t2 - t1, "s")

#####################   3. GUI and timings  #####################

if opts.nogui == False:
	from yade import qt
	v = qt.View()

	v.eyePosition = Vector3(0, -.6, .1)
	v.upVector = Vector3(0, 0, 1)
	v.viewDir = Vector3(0, 1, 0)
	#	v.grid=(False,True,False)

	rndr = yade.qt.Renderer()
	#rndr.shape=False
	#rndr.bound=True

wallTime = time.time() - startTime
f = open("timings.txt", "a")
f.write(fileName + '_' + material + ' ' + str(O.time) + ' ' + str(wallTime) + '\n')

## To play interactively with mpi execution:
## mp.mpirun(100,numThreads,True) #'True' so we see merged scene after the run
## eraseEscapedParticles()
## mp.mpirun(100,numThreads,True)
## etc.
