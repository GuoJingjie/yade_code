.. _tutorial-thermal1:

Today we will learn how to build a script that simulates heat conduction through a 
spherical packing and compares the numerical values to Fourier's analytical solution.

Day 3 - Thermal Hands-on part 1
===============================

We know where to start, let's import the necessary libraries and set our variables:

.. code-block:: python

    from yade import pack
    from yade import timing
    import numpy as np

    num_spheres=1000
    young=1e6
    rad=0.003

    mn,mx=Vector3(0,0,0),Vector3(1.0,0.008,0.008) 

These are all recognizable variables from previous hands-on sessions. Next, we append our materials and walls as we've done in the past:

.. code-block:: python

    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=radians(3),density=2600,label='spheres'))
    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=0,density=0,label='walls'))
    walls=aabbWalls([mn,mx],thickness=0,material='walls')
    wallIds=O.bodies.append(walls)

    O.bodies.append(pack.regularOrtho(pack.inAlignedBox(mn,mx),radius=rad,gap=-1e-8,material='spheres'))

Here we see that we are appending a new type of sphere packing called `regularOrtho`. As the name suggests, this creates a regular 
orthogonal packing which will be useful for ensuring that randomness doesn't affect our comparison to the analytic conduction solution 
later. 

Next, we need to construct our engines list as usual:

.. code-block:: python

    O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=intRadius),Bo1_Box_Aabb()]),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=intRadius),Ig2_Box_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_FrictPhys()],
            [Law2_ScGeom_FrictPhys_CundallStrack()],label="iloop"
        ),
        FlowEngine(label="flow"),
        ThermalEngine(label='thermal'),
        GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8),
        VTKRecorder(iterPeriod=500,fileName='VTK'+timeStr+identifier+'/spheres-',recorders=['spheres','thermal','intr'],dead=1,label='VTKrec'),
        NewtonIntegrator(damping=0.2)
    ]

Most of this should look familiar based on our previous hands-on sessions. But we see two important components including `FlowEngine` and `ThermalEngine`. 
These two engines rely intimately on one another for simulating THM processes, and thus `ThermalEngine` *cannot* be used without `FlowEngine`. 
We instantiate both of these engines without setting any parameters so that we can do so in detail in following steps.

We are only interested in validating the thermal conduction scheme in Yade, so we need to turn many default functionalities off 
starting with body dynamics:

.. code-block:: python

    for b in O.bodies:
	if isinstance(b.shape, Sphere): 
		b.dynamic = False

`b.dynamic` is a body parameter which tells Yade to consider it for force calculations or not. Setting this value to false ensures that these spheres
will not move during the entirety of our simulation. 

Next, we set our thermal parameters:

.. code-block:: python

    thermal.conduction = True
    thermal.thermoMech = False
    thermal.advection = False
    thermal.fluidThermoMech = False
    thermal.solidThermoMech = False
    thermal.fluidConduction = False

    thermal.bndCondIsTemperature = [1,1,0,0,0,0]
    thermal.thermalBndCondValue = [0,0,0,0,0,0]
    thermal.particleDensity = 2600 # kg/m^3
    thermal.particleT0 = 400 # K
    thermal.particleCp = 710 #J(kg K) 
    thermal.particleK = 2. #W/(mK)
    thermal.particleAlpha = 11.6e-3
    thermal.useKernMethod = False

The full set of available `ThermalEngine` parameters and all their specific details can be found `here <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.ThermalEngine>`_ inside our Class Reference.
We see that we need to ensure many of the functionalities are set to `False` for the basic conduction example here. Next, we set our boundary conditions in the same way we 
learned how to set boundary conditions during the previous `FlowEngine` hands-on session. Meanwhile, the initial temperature of the particles is set with `particleT0`. Finally, we set the basic thermal conduction parameters such as the particle density (`particleDensity`), 
thermal conductivity (`particleK`), heat capacity (`particleCp`), and diffusivity (`particleAlpha`). 

Now we need to employ the `FlowEngine` for one step so that it can identify the boundaries for our `ThermalEngine`. We do not require the `FlowEngine` beyond this step because 
we are not simulating any fluid fluxes in the present conduction example:

.. code-block:: python
    
    O.dt=1.
    O.dynDt=False

    flow.updateTriangulation=True
    flow.dead=0
    flow.emulateAction()
    flow.dead=1

Here we see that we are forcing `FlowEngine` to update the triangulation in a fake timestep with `flow.emulateAction`. Once this is done, we reset the `FlowEngine` to `dead=1`
So that we do not waste computational effort calculating pressure fields. 

Gathering field data
-----------------------------------

Since we are comparing our numerical conduction to an analytical scheme, we need a way to obtain field data from 
arbitrary coordinates. Here is an example of one way to do so:

.. code-block:: python

    def bodyByPos(x,y,z):
        cBody = O.bodies[1]
        cDist = Vector3(100,100,100)
        for b in O.bodies:
            if isinstance(b.shape, Sphere):
                dist = b.state.pos - Vector3(x,y,z)
                if np.linalg.norm(dist) < np.linalg.norm(cDist):
                    cDist = dist
                    cBody = b
        print('found closest body ', cBody.id, ' at ', cBody.state.pos)
        return cBody

Where we simply feed it arbitrary coordinates and it will return the closest body with which we can extract physical quantities such as 
temperature, velocity, etc.

Let's use this function to grab 10 bodies along the x-axis for us to track during the simulation:

.. code-block:: python

    axis = np.linspace(mn[0], mx[0], num=11)
    axisBodies = [None] * len(axis)
    axisTrue = np.zeros(len(axis))
    for i,x in enumerate(axis):
        axisBodies[i] = bodyByPos(x, mx[1]/2, mx[2]/2)
        axisTrue[i] = axisBodies[i].state.pos[0]


Additionally, we need a way to compute the analytical solution. Here is the solution to the heat equation
for a uniform initial temperature condition and boundary conditions at 0 K:

.. code-block:: python

    k = 2  
    Cp = 710
    rho = 2600
    alpha = 6.*k/(np.pi*Cp*rho)

    def analyticalHeatSolution(x,t,u0,L,alpha):
        ns = np.linspace(1,1000,1000)
        solution = 0
        for i,n in enumerate(ns):
            integral = (-2./L)*u0*L*(np.cos(n*np.pi)-1.) / (n*np.pi)
            solution += integral * np.sin(n*np.pi*x/L)*np.exp((-alpha*(n*np.pi/L)**2)*t)
        return solution

Where `x` is the x coordainte along the x-axis, `t` is the time of measurement, `u0` is the initial temperature of the rod, `L` is the length of the rod, and `k` is the 
thermal diffusivity of the rod. `alpha` is an effective thermal diffusivity which scales the discrete elements to cubical continuum elements. 

Finally, we need to collect and plot the data during the simulation. The temperature can be obtained via the `body state <https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.ThermalState>`_.
And you have the bodies of interest set in `axisBodies`. Using the information from previous hands-on sessions, fill out the following template to collect data for 

.. code-block:: python

    def history():
        plot.addData(
            t = O.time,
            i = O.iter,
            temp1 = ________,
            temp2 = ________,
            temp3 = ________,
            tempAnalytic1 = analyticalHeatSolution(________),
            tempAnalytic2 = analyticalHeatSolution(________),
            tempAnalytic3 = analyticalHeatSolution(________)
    )
        plot.saveDataTxt('conductionAnalyticalComparison.txt',vars=('t','i','temp1','temp2','temp3','tempAnalytic1','tempAnalytic2','tempAnalytic3'))

    O.engines=O.engines+[PyRunner(iterPeriod=500,command='history()',label='recorder')]

Use the lessons we learned from previous hands-on sessions to:

#. plot the comparison between the numerical temperature and the analytical temperature.
#. ensure that our `VTKRecorder` is also collecting and printing files for paraview.
#. start the simulation.

Example script
---------------

Please find a full script located in the `examples folder <https://gitlab.com/yade-dev/trunk/-/blob/master/examples/ThermalEngine/conductionVerification.py>`_