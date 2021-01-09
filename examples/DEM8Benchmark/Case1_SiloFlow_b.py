#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 

# Benchmark of basic performance of open-source DEM simulation systems
# Case 1: Silo flow


# Configure MPI module if needed
numThreads=4
from yade import mpy as mp
if numThreads>1:
    mp.DOMAIN_DECOMPOSITION =True
    mp.ACCUMULATE_FORCES=False
    mp.MERGE_W_INTERACTIONS=False
    mp.REALLOCATE_MINIMAL = False # if true, intersections are minimized before reallocations, hence minimizing the number of reallocated bodies
    mp.REALLOCATE_FREQUENCY=20
    mp.USE_CPP_REALLOC=True
    mp.MINIMAL_INTERSECTIONS = False
    mp.MASTER_UPDATE_STATES = False
    mp.ERASE_REMOTE_MASTER = True
    #mp.mpirun(1,numThreads,True) #this is to eliminate initialization overhead in Cundall number and timings
    #mp.MINIMAL_INTERSECTIONS = True
    mp.YADE_TIMING=False
else:     
    O.timingEnabled=True


# The stl file is in mm -> We can keep the units below or scale the stl file to meters (ymport.stl does not have a scale parameter, but I can edit the .stl file manually)

# -------------------------------------------------------------------- #
# Units: N/mm/MPa/ton
# 	force:   1 N
# 	length:  1 mm = 1e-3 m
# 	stress:  1 MPa = 1 N/mm^2 = 1e6 Pa = 1e-3 GPa
#	mass:    1 ton = 1e3 kg
#	density: 1 ton/mm^3 = 1e-12 kg/m^3

# -------------------------------------------------------------------- #
# Input Data -> Define Material and Orifice size. Uncomment the prefered choice

#FIXME: most of what follows don't have to be guarded by rank==0
fileName='SiloLargeOrifice'
#fileName='SiloSmallOrifice'
if fileName=='SiloLargeOrifice':
    z=70	# This is the height of the lowest point of the funnel (at the orifice), measuring from the lowest cylindrical cross section of the silo
elif fileName=='SiloSmallOrifice':
    z=80

particleRadius = 3 # 2mm in the benchmark, bigger for less particles

granularMaterial='M1'
#granularMaterial='M2'

if mp.rank==0:

    # -------------------------------------------------------------------- #
    # Materials
    Steel = O.materials.append(FrictMat(young=210e-3,poisson=0.2,density=7200e-12,label='Steel'))

    # -------------------------------------------------------------------- #
    # Asign values based on the Material and Orifice size
    # Coeff of restitution (e) / Coeff of friction (f)
    e_M1_M2=0.45;	f_M1_M2=0.2
    e_M1_M1=0.5;	f_M1_M1=0.3
    e_M1_St=0.4;	f_M1_St=0.2
    e_M2_M2=0.4;	f_M2_M2=0.4
    e_M2_St=0.4;	f_M2_St=0.2

    if granularMaterial=='M1':
        M1=O.materials.append(FrictMat(young=1.0e-3,poisson=0.2,density=2500e-12,label='M1'))
        e_gg=e_M1_M1	# Coefficient of restitution (e) between granular material (g) and granular material (g)
        f_gg=f_M1_M1	# Coefficient of friction (f)...

        e_gs=e_M1_St	# Coefficient of restitution (e) between granular material (g) and steel (s)
        f_gs=f_M1_St	# Coefficient of friction (f)...
    elif granularMaterial=='M2':
        M2=O.materials.append(FrictMat(young=0.5e-3,poisson=0.2,density=2000e-12,label='M2'))
        e_gg=e_M2_M2
        f_gg=f_M2_M2

        e_gs=e_M2_St
        f_gs=f_M2_St

    F_gg=atan(f_gg) # Friction Angle between granular material (g) and granular material (g)
    F_gs=atan(f_gs) # Friction Angle between granular material (g) and steel (s)


    # -------------------------------------------------------------------- #
    ## Engines 
    O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],label="collider"),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_MindlinPhys(
                frictAngle = MatchMaker(matches=((1,1,F_gg),(0,1,F_gs))), # 0 being the id of Steel and
                en         = MatchMaker(matches=((1,1,e_gg),(0,1,e_gs)))  # 1 being the id of granularMaterial
            )],
            [Law2_ScGeom_MindlinPhys_Mindlin()],
        ),
        NewtonIntegrator(damping=0,gravity=[0,0,-9.810],label="newton"),
        #GlobalStiffnessTimeStepper(active=1,timestepSafetyCoefficient=0.8, timeStepUpdateInterval=100, parallelMode=False, label = "ts",defaultDt=PWaveTimeStep()), #FIXME Remember to reinstate parallelMode=True when we use MPI
        #VTKRecorder(virtPeriod=0.04,fileName='/tmp/Silo-',recorders=['spheres','facets']),
    ]



    # -------------------------------------------------------------------- #
    # Generate initial packing. Choose among regularOrtho, regularHexa or randomDensePack (which I think is best).

    ## Using regularOrtho
    #sp=pack.regularOrtho(pack.inCylinder(Vector3(0,0,0),Vector3(0,0,305),radius),radius=2,gap=r*1/10.,material=granularMaterial)

    ## Using regularHexa
    #sp=pack.regularHexa(pack.inCylinder(Vector3(0,0,0),Vector3(0,0,215),radius),radius=2,gap=r*1/10.,material=granularMaterial)  

    # Using randomDensePack
    sp=pack.randomDensePack(pack.inCylinder((0,0,0),(0,0,254),100),radius=particleRadius,spheresInCell=500,returnSpherePack=False,material=granularMaterial,seed=1) # FIXME: We can maybe increase spheresInCell to 5k or 10k. I used 500 to have a fast generation of the sample at this stage. If you prefer, we can load the packing from a db using memoizeDb, to achieve a faster sample generation.

    # -------------------------------------------------------------------- #
    # Sort packing in ascending Z coordinates and delete excess particles to achieve sample size of 122k
    zValues=[]
    for s in sp:
        zValues.append(s.state.pos[2])

    from operator import itemgetter
    indices, zValues_sorted = zip(*sorted(enumerate(zValues), key=itemgetter(1)))
    list(zValues)
    list(indices)

    sp_new=[]
    for i in range(0,len(sp)):
        sp_new.append(sp[indices[i]])

    Nspheres=122000
    sp_new=sp_new[0:Nspheres]

    from yade import ymport
    if not os.path.exists(fileName+'.stl'):
        print("Downloading mesh file")
        try:
            os.system('wget http://perso.3sr-grenoble.fr/users/bchareyre/yade/input/'+fileName+'.stl')
        except:
            print("** probably no internet connection, grab the *.stl files by yourself **")
    facets = ymport.stl(fileName+'.stl',color=(0,1,0),material=Steel)
    fctIds = range(len(facets))
    
    O.bodies.append(facets)
    O.bodies.append(sp_new)
        
    # -------------------------------------------------------------------- #
    # Count the number of spherical particles to verify sample size. We can comment this out later on.
    numSpheres=0
    for b in O.bodies:
        if isinstance(b.shape,Sphere):
            numSpheres=numSpheres+1
    print('The total number of spheres is: ',numSpheres)
        
    collider.verletDist = 0.2*particleRadius
    O.dt=0.6 * PWaveTimeStep()
    O.dynDt=False


#-------------------------------------------------------------------- #
#Record time-dependent number of retained particles and vtk export

from yade import plot
plot.plots={'time':('retained')}

numErased = 0
def addPlotData(): 
	plot.addData(retained=numSpheres-numErased, time=O.time)

# Can't be used very easily with MPI split/merge in between, so I'm calling in a loop below (Bruno)
#O.engines=O.engines+[PyRunner(iterPeriod=100,command='addPlotData()')] # FIXME: Remember to revisit the iterPeriod

from yade import export
vtk = export.VTKExporter("spheres")

# -------------------------------------------------------------------- #
# Erase particles flowing out of the silo

def eraseEscapedParticles():
    global numErased
    count=0
    ts=time.time()
    ers=[]
    for b in O.bodies:
        if isinstance(b.shape,Sphere) and b.state.pos[2]<-z-20: # I do not delete the particles right after they pass the orifice, to disturb the simulation as little as possible
            ers.append(b.id)
            #O.bodies.erase(b.id)
            count+=1
    test = mp.bodyErase(ers)
    numErased+=count
    mp.mprint("erased",test,"/",count," in",time.time()-ts,"s")


# -------------------------------------------------------------------- #
# Run iterations 

t1=time.time()

#while numSpheres-numErased>0:
for k in range(2000):
    #if mp.rank>0: mp.mprint("has",len(O.subD.ids),"particles")
    if mp.rank==0: mp.mprint("running 100+",O.iter)
    mp.mpirun(50,numThreads,True)
    if mp.rank==0:
        addPlotData()
        vtk.exportSpheres(what=dict(particleVelocity='b.state.vel',domain='b.subdomain'))
    #if mp.rank>0: mp.mprint("has",len(O.subD.ids),"particles")
    eraseEscapedParticles()

t2=time.time()
mp.mprint("CPU time:",t2-t1)
print("erased =", numErased)

if mp.rank==0:
    print("Cundall number:", 0.5*numSpheres*O.iter/(t2-t1)) # computed with N*0.5 since we go from N to zero
    plot.plot(noShow=True).savefig(fileName+'_'+granularMaterial+'.png')
    plot.saveDataTxt(fileName+'_'+granularMaterial+'.txt')

## -------------------------------------------------------------------- #
## GUI
#if opts.nogui==False:
	#from yade import qt
	#v=qt.View()

	#v.eyePosition = Vector3(0,-600,100)
	#v.upVector    = Vector3(0,0,1)
	#v.viewDir     = Vector3(0,1,0)
##	v.grid=(False,True,False)

	#rndr=yade.qt.Renderer()
	##rndr.shape=False
	##rndr.bound=True

## To play interactively with mpi execution:
## mp.mpirun(100,numThreads,True) #'True' so we see merged scene after the run
## eraseEscapedParticles()
## mp.mpirun(100,numThreads,True)
## etc.
