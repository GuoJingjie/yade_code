#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 

# Benchmark of basic performance of open-source DEM simulation systems
# Case 3: Penetration test

# for efficient production run with the following command (it will terminate automatically when done):
#   'yadedaily -n -x Case3_PenetrationTest_SI.py 25000 M1'
# for keeping GUI run like this:
#   'yadedaily Case3_PenetrationTest_SI.py 25000 M1'
# adapt yade executable if not 'daily' version and pick one in 25000,50000,10000
# the input scripts are downloaded from TUHH as part of the script execution if not available in current path
# provided input should be in ./inputData relative to where yade is executed

# Units: SI (m, N, Pa, kg, sec)
# -------------------------------------------------------------------- #
# Input Data -> Define Number of particles. Uncomment the prefered choice. The rest of the script should be automated for all different scenarios.
# additional arg, if any is used to define N, else 25k is default
if len(sys.argv)>1: N=int(sys.argv[1])
else: N=25000
if N not in [25000,50000,100000]: print("input error, call 'yadedaily -n -x scriptName.py' or 'yadedaily -n -x scriptName.py 25000' (or 50000, etc.)")

material=str(sys.argv[2]) if len(sys.argv)>2 else 'M1'

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
Steel=O.materials.append(FrictMat(young=210e9,poisson=0.2,density=7200,frictionAngle=atan(0.2),label='Steel'))

# -------------------------------------------------------------------- #
# Assign coeff of restitution (e)
e_M1_M1=0.5;
e_M1_St=0.4;

M1=O.materials.append(FrictMat(young=1.0e9,poisson=0.2,density=2500,frictionAngle=atan(0.3),label='M1'))
e_gg=e_M1_M1	# Coefficient of restitution (e) between granular material (g) and granular material (g)
e_gs=e_M1_St	# Coefficient of restitution (e) between granular material (g) and steel (s)

# -------------------------------------------------------------------- #
# Load the initial sphere pack from .txt file

inputFileName = 'inputData/'+str(int(N/1000))+'KParticles.txt'
altName = 'inputData/'+str(int(N/1000))+'KParticlesSwapped.txt'

urls = {}
urls["25000"]="https://cloud.tuhh.de/index.php/s/9BHg3p39FzC7pYL/download"
urls["50000"]="https://cloud.tuhh.de/index.php/s/28ETraDs3fbY3xo/download"
urls["100000"]="https://cloud.tuhh.de/index.php/s/Pw2Tyc7Aw9g2kCB/download"
#urls["walls"]="https://cloud.tuhh.de/index.php/s/49HZwje9zj4R78m/download"  # <======= NOT USED / Need a fix - we are using data from yade-dem  herefafter 
urls["walls"]="http://yade-dem.org/publi/data/DEM8/Case3_PenetrationTest_Walls.txt"

# download from organizer at TUHH
if not os.path.exists(inputFileName):    
    os.system("wget -O "+inputFileName+" "+urls[str(N)])

# convert data with radius in last column
if not os.path.exists(altName):
    skipFirst=True # skip column titles
    with open(inputFileName) as x, open(altName, 'w') as y:
        for line in x:
            if skipFirst:
                skipFirst = False
                continue
            columns=line.split()
            # BEGIN HACK: fix provided input files
            if len(columns)<4: continue
            if len(columns)>4: # must be big sphere line (hopefully)
                    y.write('0 0 0.06 0.01\n')
                    continue
            # END HACK
            # nomal case, radius goes to last column
            columns = columns[1:]+columns[:1]
            y.write('\t'.join(columns)+'\n')
    y.close()

from yade import ymport
sp=ymport.text(altName,material= (M1 if material=="M1" else M2))

sp[-1]=sphere(sp[-1].state.pos,sp[-1].shape.radius,material='Steel') # Redefine this with the correct material, as sp[-1] corresponds to the large sphere made of steel
O.bodies.append(sp)

# -------------------------------------------------------------------- #
# Load facets making the box
wallFileName = 'inputData/Case3Walls'+str(int(N/1000))+'.txt' # a bit pointless to download again with different N but it's to avoid collisions between parallel jobs

# download from organizer at TUHH
if not os.path.exists(wallFileName):    
    os.system("wget -O "+wallFileName+" "+urls["walls"])

from yade import ymport
facets = ymport.textFacets(wallFileName,color=(0,1,0),material=Steel)
fctIds = range(len(facets))
O.bodies.append(facets)

# -------------------------------------------------------------------- #
# Timestep
O.dt=5e-7 #as required by guidelines
O.dynDt=False # Alternatively, we could use this increased timestep or consider using the GST.

# -------------------------------------------------------------------- #
## Engines 
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],verletDist=-0.1),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(
			en         = MatchMaker(matches=((1,1,e_gg),(0,1,e_gs)))  # 0 being the id of Steel and 1 being the id of M1 # ,(0,0,e_ss) # e_ss is not needed, as the steel ball is not supposed to touch the steel box. Was defined for debugging purposes above
		)],
		[Law2_ScGeom_MindlinPhys_Mindlin()],
	),
	NewtonIntegrator(damping=0,gravity=[0,0,-9.81]),
	#VTKRecorder(iterPeriod=200,fileName='/tmp/Case3_PenetrationTest-',recorders=['spheres']),
]

# -------------------------------------------------------------------- #
# Initial condition (assumption: the big sphere is the last one)
bSphere=O.bodies[N]
bSphere.state.vel=[0,0,-5]

# -------------------------------------------------------------------- #
# Record time-dependent position of the large sphere (bSphere)
from yade import plot
plot.plots={'time1':'z'}
def addPlotData(): 
	z=bSphere.state.pos[2]
	plot.addData(z=z, time1=O.time)
	plot.plot(noShow=True).savefig('outputData/Case3_PenetrationTest_'+str(N)+material+'.png')
	plot.saveDataTxt('outputData/Case3_PenetrationTest_'+str(N)+material+'.txt',vars=('time1','z'))

addPlotData() #record the initial state for O.time=0.0
O.engines=O.engines+[PyRunner(virtPeriod=5e-4,command='addPlotData()')] # Here I use virtPeriod=0.0005, following the provided .xlsx example file

# -------------------------------------------------------------------- #
# GUI
if opts.nogui==False:
	Gl1_Sphere.stripes=True
	plot.plot(noShow=False)
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

O.stopAtTime=0.1 # These values were taken from the .xlsx example file
O.run(-1,wait = opts.nogui)

