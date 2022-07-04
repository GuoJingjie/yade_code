.. _tutorial-fluids:

Today we will build a script which will simulate fluid flow through a spherical packing
using Yade's FlowEngine. 

Day 2 - Fluids Hands-on part 1
=======================

First let's import the libraries and set some parameters for future tweaking:

.. code-block:: python

    from builtins import range
    from yade import pack

    num_spheres = 1000  # number of spheres
    young = 1e6
    compFricDegree = 3  # initial contact friction during the confining phase
    finalFricDegree = 30  # contact friction during the deviatoric loading
    mn, mx = Vector3(0, 0, 0), Vector3(1, 1, 1)  # corners of the initial packing


Next, we already know how to add materials and geometry:

.. code-block:: python
    
    # append sphere and wall materials
    O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=radians(compFricDegree), density=2600, label='spheres'))
    O.materials.append(FrictMat(young=young, poisson=0.5, frictionAngle=0, density=0, label='walls'))

    # create and append 4 walls of a cube sized to our mn, mx parameters
    walls = aabbWalls([mn, mx], thickness=0, material='walls')
    wallIds = O.bodies.append(walls)

    # use makeCloud to generate a cloud of spheres inside our mn, mx bounds
    sp = pack.SpherePack()
    sp.makeCloud(mn, mx, -1, 0.3333, num_spheres, False, 0.95, seed=1)  #"seed" make the "random" generation always the same
    sp.toSimulation(material='spheres')


These commands should all look familiar after passing the previous two tutorials. In brief, we are appending the `FrictMat`
material type, then we assign that material to a set of walls which we then append to the scene with `O.bodies.append(walls)`. 
Following the walls, we create and append the spheres. 

Notice how we add the walls **first** and then we add the spheres. `FlowEngine` expects by default to see the walls in the first 6 bodies
(ids 0 through 5). If we need to place the walls in a different location, we can do so but we would need to set 
additional parameters in the `FlowEngine` engine. For now, we append the walls first.

Triaxial Stress Control
------------------------

Next, we will create our `TriaxialStressController` (`full parameter list with descriptions found here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.TriaxialStressController>`_) and set some standard parameters to it:

.. code-block:: python

    triax = TriaxialStressController(
            internalCompaction=True,  
            stressMask=7,
            goal1=-10000,
            goal2=-10000,
            goal3=-10000,
            maxMultiplier=1. + 2e4 / young,  # spheres growing factor (fast growth)
            finalMaxMultiplier=1. + 2e3 / young,  # spheres growing factor (slow growth)
    )

Most of these parameters are geared towards how the stress is applied and achieved inside our specimen. 
`internalCompaction`  tells `TriaxialStressController()` that we want the stress to be achieved by holding
the walls fixed and growing the particles until the desired stress is achieved. `stressMask` is an integer between 
0 and 7 which indicates the loading conditions (stress or strain, or which axis, `more details found here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.TriaxialStressController.stressMask>`_).
`stressMask = 7` tells `TriaxialStressController()` that we want all axes loaded to a constant stress condition. 
`goalX` indicates the value along each of the three axes. So here we are asking for all 3 axes to achieve a constant 
compressive stress of -10000. `maxMultiplier` and `finalMaxMultiplier` control how quickly the particles can grow,
`more details found here. <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.TriaxialStressController.maxMultiplier>`_

Engine list
-----------

Next, we will set up our engine list, as usual:

.. code-block:: python

    O.engines = [
            ForceResetter(),
            InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
            InteractionLoop(
                    [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], 
                    [Ip2_FrictMat_FrictMat_FrictPhys()], 
                    [Law2_ScGeom_FrictPhys_CundallStrack()],
                    label="iloop"
            ),
            FlowEngine(dead=1, label="flow"),  #introduced as a dead engine for the moment, see 2nd section
            GlobalStiffnessTimeStepper(active=1, timeStepUpdateInterval=100, timestepSafetyCoefficient=0.8),
            triax,
            NewtonIntegrator(damping=0.2, label="newton")
    ]


This should look familiar based on the previous two tutorials we completed. In summary, we need to ensure that 
Yade knows to expect collisions between our spheres and our walls (boxes), so we add the `Ig2_Sphere_Sphere_ScGeom()` 
and the `Ig2_Box_sphere_ScGeom()`. Here we will stick to the classic Cundall Strack contact law. Next we add the `FlowEngine`
which is set to `dead=1` so that we can run some non-flow time steps before initiating our flow simulation (see below).
Here we introduce a new engine called the `GlobalStiffnessTimeStepper` which will automatically control the timestep
during the simulation (`see more details here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.GlobalStiffnessTimeStepper>`_).
We then see the placement of our predefined `triax` followed by the familiar `NewtonIntegrator`. Our engine list 
now contains all the engines necessary to run a fluid-coupling simulation in Yade. 


Finding an equilibrated state
-----------------------------

But before running and fluid simulation, we need our spheres to be in a balanced and packed state. In order to 
achieve this, we can run some steps and check the `unbalancedForce()` while the particles grow (remember, we 
set `internalCompaction=True`):

.. code-block:: python

    while 1:
        O.run(1000, True)
        unb = unbalancedForce()
        if unb < 0.001 and abs(-10000 - triax.meanStress) / 10000 < 0.001:
            break


This while loop will start by telling Yade to run 1000 iterations through our `O.engines` list. Next it will check
the total `unbalancedForce()` between all the particles. Finaly, it will ensure that the meanStress is close to our
desired stress. If the unbalanced force and mean stress are not adequate, it will repeat the proces again until the 
break criteria is satisfied. 

When this loop is completed, we know we have achieved a packed state, and we can check this visually by activating
the viewer:

.. code-block:: python

    yade.qt.View()


It is common to keep the friction low to expedite the unbalanced force phase. But once the packing
is acheived, we can simply increase the friction to match our physical properties:

.. code-block:: python

    setContactFriction(radians(finalFricDegree))


Setting up the `FlowEngine`
^^^^^^^^^^^^^^^^^^^^^^^^^

we are almost ready to run a fluid coupled test, but first we want to set up the `FlowEngine` parameters:

.. code-block:: python

    flow.dead = 0
    # boundaries
    flow.bndCondIsPressure = [0, 0, 1, 1, 0, 0]
    flow.bndCondValue = [0, 0, 1, 0, 0, 0]
    flow.boundaryUseMaxMin = [0, 0, 0, 0, 0, 0]
    # permeability control
    flow.permeabilityFactor = 1
    flow.viscosity = 10
    # remeshing criteria


All these parameters, and more, can be found with `full descriptions here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.FlowEngine>`_. 
`flow.dead = 0` tells Yade that we now want to activate the `FlowEngine`. Next we set the boundary conditions using 
`bndCondIsPressure` and `bndCondValue`. These tell `FlowEngine` whch boundaries should have a dirichlet boundary condition
and what that pressure value should be at those boundaries. `boundaryUseMaxMin` tells `FlowEngine` if the boundary should
be set automatically using max min coordinates of the bodies, or if it should use the locations of the appended walls. We
appended walls and thus set all 6 components of this array to False (0). 

Next we are setting the permeability parameters. `permeabilityFactor=1` tells `FlowEngine` that the permeability 
between pores should be set according to the Poisseuille equation. More details associated with this parameter can be 
found `in the class reference <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.FlowEngine.permeabilityFactor>`_. 
Similar to `permeabilityFactor`, `viscosity` sets the viscosity used within the Poisseuille equation as well as the
viscous forces. 

Remeshing parameters
--------------------

Understanding the remeshing methods in `FlowEngine` is integral to using 
the `FlowEngine` properly. During our presentations, you saw how `FlowEngine` uses a Delauay triangulation with
a Voronoi dual to triangulate the pores. However, as the particles are moving, the mesh also needs to be 
re-computed since all the geometrical information associated with each of the pores will change (which changes
permeability and force integrals). This remeshing process is expensive, so we need to find a way to remesh
frequently enough that we capture the deformation, but not too frequently that the computer spends all of its time
remeshing instead of running the simulation. We control the frequency of remeshing using the following parameters:

.. code-block:: python

    flow.defTolerance = 0.3
    flow.meshUpdateInterval = 200

Where the `defTolerance` is a value which detects the maximum volumetric deformation within the system and 
triggers a remesh if the deformation is in excess of this value. Meanwhile, the `meshUpdateInterval` forces a 
remesh every XXX iterations (here we are asking for a new mesh every 200 iterations). Details about these 
parameters can be `found here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.FlowEngine.defTolerance>`_.

There are a few final settings that any `FlowEngine` user should be made aware of:

.. code-block:: python

    # solver 
    flow.useSolver = 3
    # manually setting the timestep
    O.dt = 0.1e-3
    O.dynDt = False

Here we see a `useSolver` parameter which tells `FlowEngine` which of the various solvers we want to employ 
for our simulation. Both 3 and 4 are direct solvers employing a Cholesky decomposition. The difference is that 
4 is more parallelized and ready for GPU acceleration. We also set the time step here manually with `O.dt` and 
`O.dynDt = False`. This is because there is currently no automated way to set a stable timestep for `FlowEngine`
This means the user should use trial and error to find a stable timestep since it depends strongly on the 
dynamics/geometry of the simulation.

Getting the starting permeability
---------------------------------

.. code-block:: python

    O.run(1, 1)
    Qin = flow.getBoundaryFlux(2)
    Qout = flow.getBoundaryFlux(3)
    permeability = abs(Qin) / 1.e-4  #size is one, we compute K=V/∇H
    print("Qin=", Qin, " Qout=", Qout, " permeability=", permeability)


We employ an easy `FlowEngine` method called `getBoundaryFlux()` for obtaining fluxes into and out of the specimen 
for the second and third walls in our model. We can compute the permeability here (remembering that pressure = density * gravity * head).


Starting the oedometer
----------------------

The next part will require your help, we know we need new boundary conditions for the oedometer, 
so complete the `bndCondIsPressure` and `bndCondValue` entries below. 

.. code-block:: python

    flow.bndCondIsPressure = [_, _, _, _, _, _]
    flow.bndCondValue = [_, _, _, _, _, _]
    flow.updateTriangulation = True  #force remeshing to reflect new BC immediately
    newton.damping = 0

Before we start, we need to make sure we can collect data for plotting. 

.. code-block:: python

    def history():
        plot.addData(
                e22=triax.strain[_],
                t=O.time,
                p=flow.getPorePressure((_, _, _)),
                s22=triax.stress(_)[_]
        )

We can add any data collection we wish inside this function. For example, here we will collect
the triaxial strain using the `strain` function in our `TriaxialStressController`. We are also using a `FlowEngine` 
function called `getPorePressure` which lets use obtain the pore pressure at any user defined cooradinate. As we've mentioned
before, you can find a variety of additional functions in the `Yade class reference <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.FlowEngine>`_

Complete the `history()` function above before proceeding to the next code block.

We need Yade to call our `history()` function once per loop. We can do that by creating a `PyRunner`:

.. code-block:: python

    O.engines = O.engines + [PyRunner(iterPeriod=200, command='history()', label='recorder')]

Here we are appending the `PyRunner` to our existing `O.engines` list. We are telling the `PyRunner` that we 
want it to call the command `history()` once every 200 ierations. 

Plotting live data
------------------

Yade has a module for plotting live data, the details of the `plot module can be found here <https://yade-dem.org/doc/yade.plot.html>`_tutorial-fluids

Here is an example of how we can plot the data live:

.. code-block:: python

    from yade import plot
    plot.plots = {'t': (('e22', 'b--'), None, ('s22', 'g--'), ('p', 'g-'))}
    plot.plot()

The plot module is letting us plot `t` vs `e22` using a blue line (`b--`) for the principle y-axis. Meanwhile,
it is plotting `s22` and `p` using a green lines on the secondary y-axis. 

We are now all set to run the fluid coupling simulation. 

Example script
---------------

Please find a full script located in the `examples folder <https://gitlab.com/yade-dev/trunk/-/blob/master/examples/FluidCouplingPFV/oedometer.py>`_