from yade import pack, ymport, plot, utils, export, timing
import numpy as np

"""
Check test for ThermalEngine

"""
if ('THERMAL' in features):
  young=5e6

  mn,mx=Vector3(0,0,0),Vector3(0.05,0.05,0.05)

  O.materials.append(FrictMat(young=young*100,poisson=0.5,frictionAngle=0,density=2600e20,label='walls'))
  O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(30),density=2600e10,label='spheres'))

  walls=aabbWalls([mn,mx],thickness=0,material='walls')
  wallIds=O.bodies.append(walls)

  sp=pack.SpherePack()
  O.bodies.append(pack.regularHexa(pack.inAlignedBox(mn,mx),radius=0.005,gap=0))

  triax=TriaxialStressController(
  maxMultiplier=1.+2e4/young,
  finalMaxMultiplier=1.+2e3/young,
  thickness = 0,
  stressMask = 7,
  internalCompaction=True,
  )

  O.engines=[
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=1,label='is2aabb'),Bo1_Box_Aabb()]),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=1,label='ss2sc'),Ig2_Box_Sphere_ScGeom()],
    [Ip2_FrictMat_FrictMat_FrictPhys()],
    [Law2_ScGeom_FrictPhys_CundallStrack()],label="iloop"
  ),
  GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.5),
  triax,
  FlowEngine(dead=1,label="flow",multithread=False),
  ThermalEngine(dead=1,label='thermal'),
  NewtonIntegrator(damping=0.5)
  ]

  O.step()
  ss2sc.interactionDetectionFactor=-1
  is2aabb.aabbEnlargeFactor=-1

  tri_pressure = 10000
  triax.goal1=triax.goal2=triax.goal3=-tri_pressure
  triax.stressMask=7
  while 1:
    O.run(100, True)
    unb=unbalancedForce()
    #print('unbalanced force:',unb,' mean stress: ',triax.meanStress)
    if unb<0.1 and abs(-tri_pressure-triax.meanStress)/tri_pressure<0.001:
      break

  triax.internalCompaction=False

  flow.debug=False
  # add flow
  flow.permeabilityMap = False
  flow.pZero = 10
  flow.meshUpdateInterval=2
  flow.fluidBulkModulus=2.2e9
  flow.useSolver=4
  flow.permeabilityFactor=-1e-5
  flow.viscosity= 0.001
  flow.decoupleForces =  False
  flow.bndCondIsPressure=[0,0,0,0,1,1]
  flow.bndCondValue=[0,0,0,0,10,10]

  ## Thermal Stuff
  flow.bndCondIsTemperature=[0,0,0,0,0,0] 
  flow.thermalEngine=True
  flow.thermalBndCondValue=[0,0,0,0,0,0] 
  flow.tZero=25

  flow.dead=0
  thermal.dead=1

  thermal.conduction=True
  thermal.fluidConduction=True
  thermal.debug=0
  thermal.thermoMech=True
  thermal.solidThermoMech = True
  thermal.fluidThermoMech = True
  thermal.advection=True
  thermal.useKernMethod=False
  thermal.bndCondIsTemperature=[0,0,0,0,0,1]
  thermal.thermalBndCondValue=[0,0,0,0,0,45]
  thermal.fluidK = 0.650
  thermal.fluidBeta = 2e-5 # 0.0002
  thermal.particleT0 = 25
  thermal.particleK = 2.0
  thermal.particleCp = 710
  thermal.particleAlpha = 3.0e-5
  thermal.particleDensity = 2700
  thermal.tsSafetyFactor = 0 #0.01
  thermal.uniformReynolds =10
  thermal.minimumThermalCondDist=0

  timing.reset()
  O.dt=1e-2
  O.dynDt=False
  thermal.dead=0
  flow.emulateAction()

  def bodyByPos(x,y,z):
  cBody = O.bodies[1]
  cDist = Vector3(100,100,100)
  for b in O.bodies:
    if isinstance(b.shape, Sphere):
    dist = b.state.pos - Vector3(x,y,z)
    if np.linalg.norm(dist) < np.linalg.norm(cDist):
      cDist = dist
      cBody = b
  return cBody

  bodyOfInterest = bodyByPos(0.025,0.025,0.025)

  O.run(100,1)
  ftemp = flow.getPoreTemperature((0.025,0.025,0.025))
  btemp = O.bodies[bodyOfInterest.id].state.temp

  print('ftemp',ftemp,'btemp',btemp)
  tolerance = 0.01
  ftarget = 25.556741323314558
  btarget = 43.463950268408055
  if (ftemp - ftarget)/ftarget > tolerance:
    raise YadeCheckError('ThermalEngine checktest: fluid temp incorrect')
  if (ftemp - btarget)/btarget > tolerance:
    raise YadeCheckError('ThermalEngine checktest: body temp incorrect')
else:
  print("This checkThermalEngine.py cannot be executed because ENABLE_THERMAL is disabled")

