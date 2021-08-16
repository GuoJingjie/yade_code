# -*- encoding=utf-8 -*-
################################################################################
# Pack particles inside of a torus imported as an STL file. Then, confine the
# aggregate by applying an isotropic Cauchy stress tensor on the model's
# power diagram described by the alpha shape surface reconstruction algorithm 
################################################################################
import gts
from yade import ymport,pack
#######################################
###       FUNCTION DEFINITIONS      ###
#######################################
def predStl(stlF): 
    facets = ymport.stl(stlF)
    s = gts.Surface()
    for facet in facets: # creates fts.Face for each facet. The vertices and edges are duplicated
        vs = [facet.state.pos + facet.state.ori*v for v in facet.shape.vertices]
        vs = [gts.Vertex(v[0],v[1],v[2]) for v in vs]
        es = [gts.Edge(vs[i],vs[j]) for i,j in ((0,1),(1,2),(2,0))]
        f = gts.Face(es[0],es[1],es[2])
        s.add(f)
    s.cleanup(1e-4) # removes duplicated vertices and edges
    assert s.is_closed()
    return inGtsSurface(s)

def averageVelGrad(xq,Aq,vol0): # Compute Velocity Gradient of the Aggregate
  velGrad = Matrix3.Zero
  for i, el in enumerate(Aq):
    velGrad = velGrad + Vector3.outer(xq[i],Aq[i]) / vol0
  return velGrad
  
def triStressMicro(stressPt,triThresh,creepThresh): # Apply Uniform Stress Tensor
  while 1:	
    O.run(1000, True) 
    meanStress = TW.calcAlphaStress(alpha,alphaShrinked,True)
    triDiff = stressPt.diagonal() - meanStress.diagonal()
    caps=TW.getAlphaCaps(alpha,alphaShrinked,True)
    velArr = [O.bodies[cap[0]].state.vel for cap in caps]
    areaArr = [cap[1] for cap in caps]
    lbar = averageVelGrad(velArr,areaArr,TW.alphaCapsVol)
    meanBVel = 0.5 * lbar.spectralDecomposition()[1].maxAbsCoeff()
    integrator.damping=max(0.2,min(0.99,(meanBVel/boundVel)**4))
    print('Boundary Velocity:',meanBVel,'Stress Diagonal',meanStress.diagonal())
    if meanBVel < abs(creepThresh * boundVel) and abs(triDiff[0]/stressPt[0][0]) < triThresh and abs(triDiff[1]/stressPt[1][1]) < triThresh and abs(triDiff[2]/stressPt[2][2]) < triThresh:
      break
    print('damping = ',integrator.damping)

############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
confinement_p = -5.0e6
damp=0.2 # damping coefficient
fricDegree=35.0

young=160.0e9 # contact stiffness
poisson=0.17

sphereMat=FrictMat( young=young,
                    poisson=poisson,
                    frictionAngle=radians(fricDegree),
                    density=2720,
                    label='spheres')
#######################################
###            MAKE THE PACK       ###
#######################################
os.system("wget -nc https://yade-dem.org/publi/data/AlphaS/torout.stl") 
predOut = predStl('torout.stl')
os.system("wget -nc https://yade-dem.org/publi/data/AlphaS/torin.stl") 
predIn = predStl('torin.stl')

# sphs = regularHexa(predDiff,r,0)
sp=pack.SpherePack()
sp=pack.randomDensePack(predOut-predIn,radius=0.04,useOBB=True,rRelFuzz=0.5,returnSpherePack=True)
#######################################
###          DEFINE THE ENGINES     ###
#######################################
engineList = [
  ForceResetter(),
  InsertionSortCollider([ Bo1_Sphere_Aabb() ],
                          label="sorter"),
  InteractionLoop(
    [ Ig2_Sphere_Sphere_ScGeom(),
      Ig2_Facet_Sphere_ScGeom(),
      Ig2_Box_Sphere_ScGeom()],
    [Ip2_FrictMat_FrictMat_MindlinPhys( label='mindlinContact')],
    [Law2_ScGeom_MindlinPhys_Mindlin( label="mindlinLaw2")],),
  GlobalStiffnessTimeStepper( active=1,timeStepUpdateInterval=100,
                              timestepSafetyCoefficient=0.2, label='stepper'),
  NewtonIntegrator(damping=damp, label='integrator'),
  PyRunner(iterPeriod=0,
            command="TW.applyAlphaForces(stressM,alpha,alphaShrinked,True)",
            label='loader')]
#######################################
###   APPLY UNIFORM STRESS TENSOR   ###
#######################################
O.materials.append(sphereMat)
sp.toSimulation(material='spheres')
O.engines=engineList

stepper.timestepSafetyCoefficient=0.6
stressM = Matrix3.Identity * (confinement_p)
minRad = min(b.shape.radius for b in O.bodies)
meanDia = 2. * sum([b.shape.radius for b in O.bodies]) / len(O.bodies)
alpha = (15. * minRad) ** 2
alphaShrinked = (13.5 * minRad) ** 2
loader.iterPeriod=1
TW=TesselationWrapper()
integrator.damping=0.5
iNum = 7.9e-5       # Inertial Number Target
boundVel = 0.45 * iNum * sqrt(5.e6 /O.bodies[0].mat.density) / meanDia
print('Quasistatic Strain Rate', boundVel)
triStressMicro(stressM,0.01,0.1)
#######################################
###        SAVE CONFINED GZIP       ###
#######################################
O.save('confinedState.yade.gz')
print("###      Confining Pressure Achieved      ###")
#######################################
### VTK SPHERES AND POWER DIAGRAM   ###
#######################################
import yade.export
vtkExporter = yade.export.VTKExporter('vtkConfined')
vtkExporter.exportSpheres()

caps=TW.getAlphaCaps(alpha,alphaShrinked,True)
segments=TW.getAlphaGraph(alpha,alphaShrinked,True)
f=open('LVdiagram.vtk','w')
f.write("# vtk DataFile Version 2.0\nGenerated from YADE\n")
f.write("ASCII\nDATASET POLYDATA\n")
f.write("POINTS %d float \n" %(len(segments)))
[f.write("%f %f %f\n" %(point[0],point[1],point[2])) for point in segments]
f.write("LINES %d %d \n" %(len(segments)/2, 3*len(segments)/2))
[f.write("2 %d %d\n" %(2*k,2*k+1)) for k in range(int(len(segments)/2))]
f.close()
