import numpy as np, sys

# Test assertion function
def assert_passed(init_hsize, deltaE, deltaS, interactive_run=hasattr(sys, "ps1")):
    error_str = "Error in flipCell when trying to flip cell with hSize={:}, energy changed by {:.2e} and stress changed by {:.2e}".format(init_hsize, deltaE, deltaS)
    if not interactive_run: raise YadeCheckError(error_str)
    else: print(error_str)

# Test parameters
ucube = Matrix3(np.eye(3)) # unit cube, to be recovered after each flipCell
nSpheres = 4   # number of spheres along a cube side, final number will the cube of this
fail_flag = False

O.periodic = True

dx = 1/nSpheres
for i in range(nSpheres):
    for j in range(nSpheres):
        for k in range(nSpheres):
            O.bodies.append(sphere([i*dx,j*dx,k*dx],0.58*dx,fixed=True))

O.engines= [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb()],verletDist=-0.1,label='collider'),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack(label="law"),],
    ),
    NewtonIntegrator()
]

O.dt=1
O.step()
O.saveTmp()

def testFlip(velGrad,time,reset):
    if reset: O.reload() # start from unit cube
    O.cell.velGrad=velGrad
    O.run(time,True)
    O.cell.velGrad=Matrix3.Zero
    O.step()
    E = law.elasticEnergy()
    stress = getStress()
    init_hsize = O.cell.hSize
    try:
        flip=flipCell() # flip and run one step, see if the forces are unchanged through elastic energy and stress
    except:
        raise YadeCheckError("Internal error while flipping")
    O.step()
    energyChange = (E - law.elasticEnergy())/E
    stressChange = (stress - getStress()).norm() / stress.norm()
    if energyChange>1e-10 or stressChange>1e-10:
        assert_passed(init_hsize, deltaE, deltaS)
        fail_flag = True
    return flip


# Test multidirectional deformations
def runTests():
    #for test in range(int(n_random_tests)):
    hSize = np.array(ucube)
    velGrad = Matrix3.Zero
    strainRate = 0.1
    outers = [(0,1),(0,2),(1,2),(1,0),(2,0),(2,1)]
    for ij1 in outers:
        for ori1 in [-1., 1.]: # get a random orientation
            for ij2 in outers:
                if ij1 == ij2: continue
                for ori2 in [0, -0.5*ori1]:
                    velGrad = Matrix3.Zero
                    velGrad[ij1]=ori1
                    velGrad[ij2]=ori2
                    testFlip(strainRate*velGrad,8,True)
                    testFlip(strainRate*velGrad,8,False)
runTests()
if not fail_flag: print("flipCell keeps the mechanical state unchanged")
 
