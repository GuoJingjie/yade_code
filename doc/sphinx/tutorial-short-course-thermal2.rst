.. _tutorial-thermal1:

Part 2 of our Thermal Hands-on session will focus on the full THM coupling.

Day 3 - Thermal Hands-on part 2
===============================

Let's build a triaxially loaded cubic specimen:

.. code-block:: python

    from yade import pack, ymport, plot, utils, export, timing
    import numpy as np

    young=5e6

    mn,mx=Vector3(0,0,0),Vector3(0.05,0.05,0.05)

    O.materials.append(FrictMat(young=young*100,poisson=0.5,frictionAngle=0,density=2600e10,label='walls'))
    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(30),density=2600e10,label='spheres'))

    walls=aabbWalls([mn,mx],thickness=0,material='walls')
    wallIds=O.bodies.append(walls)

    sp=pack.SpherePack()
    sp.makeCloud(mn,mx,rMean=0.0015,rRelFuzz=0.333,num=100,seed=1) 
    sp.toSimulation(color=(0.752, 0.752, 0.752),material='spheres')


Here we see that we are appending a sphere cloud to the simulation (we will compact them after setting the `O.engines` list).

Next, we need to construct our engines list as usual:

.. code-block:: python

    O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=1,label='is2aabb'),Bo1_Box_Aabb()]),
        InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=1,label='ss2sc'),Ig2_Box_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()],label="iloop"
        ),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.5),
        TriaxialStressController(label='triax'),
        FlowEngine(dead=1,label="flow"),
        ThermalEngine(dead=1,label='thermal'),
        VTKRecorder(iterPeriod=500,fileName='./spheres-',recorders=['spheres','thermal','intr'],dead=1,label='VTKrec'),
        NewtonIntegrator(damping=0.5)
    ]

Now we have the full `O.engines` list set, which contains a `TriaxialStressController()` for our stress control, a `FlowEngine()` for the fluid fluxes and heat advection,
and a `ThermalEngine()` for our thermal coupling.

Compacting the specimen
-----------------------

Let's setup the `TriaxialStressController()` for our compaction:

.. code-block:: python

    triax.maxMultiplier=1.+2e4/young
    triax.finalMaxMultiplier=1.+2e3/young
    triax.thickness = 0
    triax.stressMask = 7
    triax.internalCompaction=True
    tri_pressure = 1000
    triax.goal1=triax.goal2=triax.goal3=-tri_pressure
    triax.stressMask=7

    while 1:
        O.run(1000, True)
        unb=unbalancedForce()
        print('unbalanced force:',unb,' mean stress: ',triax.meanStress)
        if unb<0.1 and abs(-tri_pressure-triax.meanStress)/tri_pressure<0.001:
            break

    triax.internalCompaction=False

Here we see that we are running a loop where we run 1000 iterations of `internalCompaction` (particles grow in radius to achieve stress), then 
testing the `unbalancedForce()` and ultimately stopping if our stopping criteria is achieved. We can't forget that our `FlowEngine()` and 
`ThermalEngine()` are both set to `dead=1` in the `O.engines` list, so they will not activate during this compaction stage. 


Setting up the `FlowEngine()`
-----------------------------

.. code-block:: python

    # initial pressure condition 
    flow.pZero = 10
    flow.meshUpdateInterval = 2
    # we will activate compressibility in the fluid
    flow.fluidBulkModulus = 2.2e9
    flow.useSolver = 4
    # enforcing a darcy permeability in the specimen 
    flow.permeabilityFactor = -1e-5
    flow.viscosity = 0.001
    # setting the boundary conditions
    flow.bndCondIsPressure = [0,0,0,0,1,1]
    flow.bndCondValue = [0,0,0,0,10,10]

    ## Thermal Stuff
    flow.bndCondIsTemperature  [0,0,0,0,0,0] 
    # activate the thermal engine
    flow.thermalEngine = True
    flow.thermalBndCondValue = [0,0,0,0,0,0] 
    # initial temperature conditions
    flow.tZero = 25

    flow.dead=0


Setting up the `ThermalEngine()`
--------------------------------

.. code-block:: python

    thermal.dead = 0
    thermal.conduction = True
    thermal.fluidConduction = True
    thermal.thermoMech = True
    thermal.solidThermoMech = True
    thermal.fluidThermoMech = True
    thermal.advection = True
    thermal.useKernMethod = False
    thermal.bndCondIsTemperature = [0,0,0,0,0,1]
    thermal.thermalBndCondValue = [0,0,0,0,0,45]
    thermal.fluidK = 0.650
    thermal.fluidBeta = 2e-5
    thermal.particleT0 = 25
    thermal.particleK = 2.0
    thermal.particleCp = 710
    thermal.particleAlpha = 3.0e-5
    thermal.particleDensity = 2700
    thermal.tsSafetyFactor = 0
    thermal.uniformReynolds = 10

We won't describe each parameter here, those descriptions can be found in the `Class Reference <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.ThermalEngine>`_. 
However, it is clear we are activating conduction, advection, the thermo-fluid mechanical coupling, the solid-fluid mechanical coupling, and fluid conduction. 
Each component can be deactivated in case the user does not need the full THM coupling. We also see a similar assignment of boundary conditions as we saw in the 
previous hands-on sessions. Some additional parameters shown here include the fluid thermal conductivity (`thermal.fluidK`), the coefficient of thermal expansion for the fluid 
(`thermal.fluidBeta`). 

Running the coupled simulation
------------------------------

The simulation is set and ready to run, first we will let `FlowEngine()` detect and assign the boundary conditions by running `flow.emulateAction()`:

.. code-block:: python

    O.dt=1e-4
    O.dynDt=False
    thermal.dead=0
    flow.emulateAction()


Now it is up to you to finish the script

#. collect the temperature at some interesting points in the specimen
#. plot the temperature 
#. export the VTK files for viewing in paraview


Example script
---------------

Please find a full script located in the `examples folder <https://gitlab.com/yade-dev/trunk/-/blob/master/examples/ThermalEngine/thermoHydroMechanical_coupling.py>`_