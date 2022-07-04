.. _short-course:

Day 1 - Yade Hands-on part 1
====================

Let's create a bouncing sphere
------------------------------

First we need to define a material for our sphere:

.. code-block:: python

    # Start by defining a material
    O.materials.append(FrictMat(young=1.0e9, poisson=0.2, density=2500, label='frictmat'))

``O`` is our scene, it contains all the information that we need to run a DEM simulation. We can edit various components of the
scene such as ``materials`` here. We use the python function ``.append()`` to add this material to the existing list of materials
inside Python. 

The ``FrictMat`` is a type of material available in Yade. Yade boasts a wide variety of `materials <https://www.yade-dem.org/doc//yade.wrapper.html#yade.wrapper.Material>`_
such as ``CohFrictMat``. These materials all have different constitutive laws associated with them. For now we stick with the simplest one, ``FrictMat``.

Next, we need to create two spheres by appending two ``bodies`` to our scene, ``O``:

.. code-block:: python

    O.bodies.append(
        [
        # fixed: particle's position in space will not change (support)
        sphere(center=(0, 0, 0), radius=.5, fixed=True),
        # this particles is free, subject to dynamics
        sphere((0, 0, 2), .5)
        ]
    )

We see that we appended a ``sphere`` to the scene by designating its center and radius. Yade has a variety of `shapes <https://www.yade-dem.org/doc//yade.wrapper.html#shape>`_ 
that can be appended as bodies such as ``Facet``, ``Box``, and others. 

Now it is time to define how these spheres should move. The scene ``O`` has an "engines" list in ``O.engines`` which is a list of actions that are taken for each
iteration in our simulation. 

.. code-block:: python

    O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom()],  # collision geometry
                [Ip2_FrictMat_FrictMat_FrictPhys()],  # collision "physics"
                [Law2_ScGeom_FrictPhys_CundallStrack()]  # contact law -- apply forces
        ),
        # Apply gravity force to particles. damping: numerical dissipation of energy.
        NewtonIntegrator(gravity=(0, 0, -9.81), damping=0.1)
    ]

Some of this may look foreign to you, but there is a logic to it. The ``ForceResetter()`` resets all forces stored in the scene, ``InsertionSortCollider`` is simply creating a 
sorted list of possible body interactions. ``InteractionLoop()`` is where we assign the interaction geometry (``Ig2_Sphere_Sphere_ScGeom()``) which conveniently 
matches our appended body shape (``Sphere``), the collision physics ``Ip2_FrictMat_FrictMat_FrictPhys()`` conveniently matches our material assignment (``FrictMat``). Finally,
the constitutive law is defined with ``Law2_ScGeom_FrictPhys_CundallStrack()``, so we are using the classical CundallStrack contact law here. The timeintegration occurs in the 
last component ``NewtonIntegrator`` where we can add a gravitational force and damping.

Although we have already added our material, body, and engines to the scene, we should still take care to define the time step in ``O.dt``:

.. code-block:: python 

    O.dt = .5 * PWaveTimeStep()

Where the ``PWaveTimeStep()`` automatically estimates the critical time step associated with the stiffness of the packing. We factor that down by 1/2 to be safe, 
since it is but an estimate.

Starting the Script
------

Now that we have the full script written for our bouncing ball, it is time to start it by executing:

.. code-block:: bash

    yade bouncing_sphere.py

in our terminal. 