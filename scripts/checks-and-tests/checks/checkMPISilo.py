#!/usr/bin/python
# 2020 © Vasileios Angelidakis <v.angelidakis2@ncl.ac.uk>
# 2020 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> 

# This script is directly derived from DEM8 benchmark case no.1, with less bodies
# and less timesteps 

# Configure MPI module if needed
numThreads=4
from yade import mpy as mp
if numThreads>1:
    mp.DOMAIN_DECOMPOSITION =True
    mp.ACCUMULATE_FORCES=False
    mp.MERGE_W_INTERACTIONS=False
    mp.REALLOCATE_MINIMAL = False
    mp.REALLOCATE_FREQUENCY=20
    mp.MINIMAL_INTERSECTIONS = True
    mp.MASTER_UPDATE_STATES = False
    mp.ERASE_REMOTE_MASTER = True
    mp.USE_CPP_REALLOC=True
    mp.VERBOSE_OUTPUT=False
    mp.SEND_BYTEARRAYS=True
    mp.YADE_TIMING=False
else:     
    O.timingEnabled=True

#FIXME: most of what follows don't have to be guarded by rank==0
fileName='SiloLargeOrifice'
#fileName='SiloSmallOrifice'
if fileName=='SiloLargeOrifice':
    z=70	# This is the height of the lowest point of the funnel (at the orifice), measuring from the lowest cylindrical cross section of the silo
elif fileName=='SiloSmallOrifice':
    z=80

particleRadius = 8 # 2mm in the benchmark, bigger for less particles

granularMaterial='M1'
#granularMaterial='M2'

if mp.rank==0:
    Steel = O.materials.append(FrictMat(young=0.1e-2,poisson=0.2,density=7200e-12,label='Steel'))
    e_M1_M2=0.45;	f_M1_M2=0.2
    e_M1_M1=0.5;	f_M1_M1=0.3
    e_M1_St=0.4;	f_M1_St=0.2
    e_M2_M2=0.4;	f_M2_M2=0.4
    e_M2_St=0.4;	f_M2_St=0.2

    if granularMaterial=='M1':
        M1=O.materials.append(FrictMat(young=0.1e-2,poisson=0.2,density=2500e-12,label='M1'))
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
    ]


    sp=pack.randomDensePack(pack.inCylinder((0,0,0),(0,0,254),100),radius=particleRadius,spheresInCell=500,returnSpherePack=False,material=granularMaterial,seed=1)
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
         
    collider.verletDist = 0.2*particleRadius
    O.dt=0.6 * PWaveTimeStep()
    O.dynDt=False


numErased = 0

def eraseEscapedParticles():
    global numErased
    count=0
    ts=time.time()
    ers=[]
    for b in O.bodies:
        if isinstance(b.shape,Sphere) and b.state.pos[2]<-z-20:
            ers.append(b.id)
            count+=1
    mp.sendCommand("all","bodyErase("+str(ers)+")",True)
    numErased+=count

mp.mpirun(1000,numThreads,True)
if mp.rank==0:
    eraseEscapedParticles()
    if numErased==32:
        mp.mprint("Parallel MPI silo -N4 succeeds")
    else:
        raise YadeCheckError("Parallel MPI silo -N4 fails")
