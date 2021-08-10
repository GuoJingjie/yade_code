# -*- encoding=utf-8 -*-
################################################################################
# 1) Generate a packing of particles inside an axis-aligned bounding box using
# makeCloud(). 2) Use triaxialStressController() with internal confinement to
# achieve an average isotropic stress on the aabb. 3) Finally, confine the
# aggregate by applying an isotropic Cauchy stress tensor on the model's
# power diagram described by the alpha shape surface reconstruction algorithm
################################################################################
import yade.export
from yade import pack
#######################################
###       FUNCTION DEFINITIONS      ###
#######################################


def averageVelGrad(xq, Aq, vol0):  # Compute Velocity Gradient of the Aggregate
    velGrad = Matrix3.Zero
    for i, el in enumerate(Aq):
        velGrad = velGrad + Vector3.outer(xq[i], Aq[i]) / vol0
    return velGrad


def triStressMicro(stressPt, triThresh, creepThresh):  # Apply Uniform Stress Tensor
    while 1:
        O.run(1000, True)
        meanStress = TW.calcAlphaStress(alpha, alphaShrinked, True)
        triDiff = stressPt.diagonal() - meanStress.diagonal()
        caps = TW.getAlphaCaps(alpha, alphaShrinked, True)
        velArr = [O.bodies[cap[0]].state.vel for cap in caps]
        areaArr = [cap[1] for cap in caps]
        lbar = averageVelGrad(velArr, areaArr, TW.alphaCapsVol)
        meanBVel = 0.5 * lbar.spectralDecomposition()[1].maxAbsCoeff()
        integrator.damping = max(0.2, min(0.99, (meanBVel / boundVel)**4))
        print('Boundary Velocity:', meanBVel, 'Stress Diagonal', meanStress.diagonal())
        if meanBVel < abs(creepThresh * boundVel) \
                and abs(triDiff[0] / stressPt[0][0]) < triThresh \
                and abs(triDiff[1] / stressPt[1][1]) < triThresh \
                and abs(triDiff[2] / stressPt[2][2]) < triThresh:
            break
        print('damping = ', integrator.damping)


def isoLoad(stabThresh, triThresh):  # Use TriaxialStressController for Initial Confinement
    while 1:
        O.run(1000, True)
        unb = unbalancedForce()
        triDiff = abs((triax.goal1 - triax.meanStress) / (triax.goal1))
        print('unb:', unb, ' triDiff: ', triDiff, ' meanS: ', triax.meanStress,
              ' porosity: ', triax.porosity)
        if unb < stabThresh and triDiff < triThresh:
            break


############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
confinement_p = -5.0e6
damp = 0.2  # damping coefficient
fricDegree = 35.0

young = 160.0e9  # contact stiffness
poisson = 0.17

sphereMat = FrictMat(young=young,
                     poisson=poisson,
                     frictionAngle=radians(fricDegree),
                     density=2720,
                     label='spheres')

wallMat = FrictMat(young=young,
                   poisson=poisson,
                   frictionAngle=0,
                   density=0,
                   label='walls')
#######################################
###            MAKE THE PACK       ###
#######################################
num_spheres = 1000   # number of spheres
mn, mx = Vector3(0, 0, 0), Vector3.Ones  # corners of the packing

sp = pack.SpherePack()
sp.makeCloud(mn, mx, num=num_spheres, rMean=0.04, rRelFuzz=0.5,
             distributeMass=True, seed=1)
#######################################
###          DEFINE THE ENGINES     ###
#######################################
engineList = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(),
                           Bo1_Box_Aabb()],
                          label="sorter"),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(),
         Ig2_Facet_Sphere_ScGeom(),
         Ig2_Box_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_MindlinPhys(label='mindlinContact')],
        [Law2_ScGeom_MindlinPhys_Mindlin(label="mindlinLaw2")],),
    GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=100,
                               timestepSafetyCoefficient=0.2, label='stepper'),
    NewtonIntegrator(damping=damp, label='integrator'),
    TriaxialStressController(label='triax'),
    PyRunner(iterPeriod=0,
             command="TW.applyAlphaForces(stressM,alpha,alphaShrinked,True)",
             label='loader')]
##########################################
### INTERNAL COMPACTION TO CONFINEMENT ###
##########################################
O.materials.append(sphereMat)
O.materials.append(wallMat)
walls = aabbWalls([mn, mx], thickness=0, material='walls')
wallIds = O.bodies.append(walls)
sp.toSimulation(material='spheres')
O.engines = engineList
triax.goal1 = triax.goal2 = triax.goal3 = confinement_p
isoLoad(0.01, 0.01)
#######################################
###   APPLY UNIFORM STRESS TENSOR   ###
#######################################
sp.fromSimulation()
O.resetThisScene()
O.materials.append(sphereMat)
sp.toSimulation(material='spheres')

engineList.remove(triax)
O.engines = engineList

stepper.timestepSafetyCoefficient = 0.6
stressM = Matrix3.Identity * (confinement_p)
minRad = min(b.shape.radius for b in O.bodies)
meanDia = 2. * sum([b.shape.radius for b in O.bodies]) / len(O.bodies)
alpha = (22. * minRad) ** 2
alphaShrinked = (21. * minRad) ** 2
loader.iterPeriod = 2
TW = TesselationWrapper()
integrator.damping = 0.5
iNum = 7.9e-5       # Inertial Number Target
boundVel = 0.45 * iNum * sqrt(5.e6 / O.bodies[0].mat.density) / meanDia
print('Quasistatic Strain Rate', boundVel)
triStressMicro(stressM, 0.01, 0.1)
#######################################
###        SAVE CONFINED GZIP       ###
#######################################
O.save('confinedState.yade.gz')
print("###      Confining Pressure Achieved      ###")
#######################################
### VTK SPHERES AND POWER DIAGRAM   ###
#######################################
vtkExporter = yade.export.VTKExporter('vtkConfined')
vtkExporter.exportSpheres()

caps = TW.getAlphaCaps(alpha, alphaShrinked, True)
segments = TW.getAlphaGraph(alpha, alphaShrinked, True)
f = open('LVdiagram.vtk', 'w')
f.write("# vtk DataFile Version 2.0\nGenerated from YADE\n")
f.write("ASCII\nDATASET POLYDATA\n")
f.write("POINTS %d float \n" % (len(segments)))
[f.write("%f %f %f\n" % (point[0], point[1], point[2])) for point in segments]
f.write("LINES %d %d \n" % (len(segments) / 2, 3 * len(segments) / 2))
[f.write("2 %d %d\n" % (2 * k, 2 * k + 1)) for k in range(int(len(segments) / 2))]
f.close()
