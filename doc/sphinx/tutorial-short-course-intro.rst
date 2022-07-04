.. _Intro:

Introduction to Bash and Python
============

Terminal
-------------

The terminal is a shell designed to let us interact with the computer and its filing system with a basic set of commands::

	user@machine:~\$                         # user operating at machine, in the directory ~ (= user's home directory)
	user@machine:~\$ ls .                    # list contents of the current directory
	user@machine:~\$ ls foo                  # list contents of directory foo, relative to the dcurrent directory ~ (= ls ~/foo = ls /home/user/foo)
	user@machine:~\$ ls /tmp                 # list contents of /tmp
	user@machine:~\$ cd foo                  # change directory to foo
	user@machine:~/foo\$ ls ~                # list home directory (= ls /home/user)
	user@machine:~/foo\$ cd bar              # change to bar (= cd ~/foo/bar)
	user@machine:~/foo/bar\$ cd ../../foo2   # go to the parent directory twice, then to foo2 (cd ~/foo/bar/../../foo2 = cd ~/foo2 = cd /home/user/foo2)
	user@machine:~/foo2\$ cd                 # go to the home directory (= ls ~ = ls /home/user)
	user@machine:~\$

Keys
----

Useful keys on the command-line are:

============= =========================
<tab>         show possible completions of what is being typed (use abundantly!)
^C (=Ctrl+C)  delete current line
^D            exit the shell
↑↓            move up and down in the command history
^C            interrupt currently running program
^\\           kill currently running program
Shift-PgUp    scroll the screen up (show past output)
Shift-PgDown  scroll the screen down (show future output; works only on quantum computers)
============= =========================

Starting yade
-------------

If yade is installed on the machine, it can be (roughly speaking) run as any other program; without any arguments, it runs in the "dialog mode", where a command-line is presented:

::

	user@machine:~\$ yade
	Welcome to Yade 2022.01a
	TCP python prompt on localhost:9002, auth cookie `adcusk'
	XMLRPC info provider on http://localhost:21002
	[[ ^L clears screen, ^U kills line. F12 controller, F11 3d view, F10 both, F9 generator, F8 plot. ]]
	Yade [1]:                                            #### hit ^D to exit
	Do you really want to exit ([y]/n)?
	Yade: normal exit.

The command-line is in fact ``python``, enriched with some yade-specific features. (Pure python interpreter can be run with ``python`` or ``ipython`` commands).

Instead of typing commands one-by-one on the command line, they can be be written in a file (with the .py extension) and given as argument to Yade::

	user@machine:~\$ yade simulation.py

For a complete help, see ``man yade``


.. rubric:: Exercises

#. Open the terminal, navigate to your home directory
#. Create a new empty file and save it in ``~/first.py``
#. Change directory to ``/tmp``; delete the file ``~/first.py``
#. Run program ``xeyes``
#. Look at the help of Yade.
#. Look at the *manual page* of Yade
#. Run Yade, exit and run it again.

Yade basics
------------

Yade objects are constructed in the following manner (this process is also called "instantiation", since we create concrete instances of abstract classes: one individual sphere is an instance of the abstract :yref:`Sphere`, like Socrates is an instance of "man"):

.. ipython::

	@suppress
	Yade [1]: from yade import *

	Yade [1]: Sphere           # try also Sphere?
	
	Yade [1]: s=Sphere()       # create a Sphere, without specifying any attributes

	Yade [1]: s.radius         # 'nan' is a special value meaning "not a number" (i.e. not defined)

	Yade [1]: s.radius=2       # set radius of an existing object

	Yade [1]: s.radius            

	Yade [1]: ss=Sphere(radius=3)   # create Sphere, giving radius directly

	Yade [1]: s.radius, ss.radius     # also try typing s.<tab> to see defined attributes

Particles
--------

Particles are the "data" component of simulation; they are the objects that will undergo some processes, though do not define those processes yet.

Singles
------

There is a number of pre-defined functions to create particles of certain type; in order to create a sphere, one has to (see the source of :yref:`yade.utils.sphere` for instance):

#. Create :yref:`Body`
#. Set :yref:`Body.shape` to be an instance of :yref:`Sphere` with some given radius
#. Set :yref:`Body.material` (last-defined material is used, otherwise a default material is created)
#. Set position and orientation in :yref:`Body.state`, compute mass and moment of inertia based on :yref:`Material` and :yref:`Shape`

In order to avoid such tasks, shorthand functions are defined in the :yref:`yade.utils` module; to mention a few of them, they are :yref:`yade.utils.sphere`, :yref:`yade.utils.facet`, :yref:`yade.utils.wall`.

.. ipython::

	@suppress
	Yade [1]: from yade import utils
	
	Yade [1]: s=utils.sphere((0,0,0),radius=1)    # create sphere particle centered at (0,0,0) with radius=1

	Yade [1]: s.shape                       # s.shape describes the geometry of the particle

	Yade [1]: s.shape.radius                # we already know the Sphere class

	Yade [1]: s.state.mass, s.state.inertia # inertia is computed from density and geometry

	Yade [1]: s.state.pos                   # position is the one we prescribed

	Yade [1]: s2=utils.sphere((-2,0,0),radius=1,fixed=True)     # explanation below

In the last example, the particle was fixed in space by the ``fixed=True`` parameter to :yref:`yade.utils.sphere`; such a particle will not move, creating a primitive boundary condition.

A particle object is not yet part of the simulation; in order to do so, a special function :yref:`O.bodies.append<BodyContainer::append>` (also see :yref:`Omega::bodies` and :yref:`Scene`) is called:

.. ipython::

	Yade [1]: O.bodies.append(s)            # adds particle s to the simulation; returns id of the particle(s) added


Packs
-----

There are functions to generate a specific arrangement of particles in the :yref:`yade.pack` module; for instance, cloud (random loose packing) of spheres can be generated with the :yref:`yade._packSpheres.SpherePack` class:

.. ipython::

	Yade [1]: from yade import pack

	Yade [1]: sp=pack.SpherePack()                   # create an empty cloud; SpherePack contains only geometrical information

	Yade [1]: sp.makeCloud((1,1,1),(2,2,2),rMean=.2) # put spheres with defined radius inside box given by corners (1,1,1) and (2,2,2)

	Yade [1]: for c,r in sp: print(c,r)               # print center and radius of all particles (SpherePack is a sequence which can be iterated over)
	   ...:

	Yade [1]: sp.toSimulation()                      # create particles and add them to the simulation

Boundaries
----------

:yref:`yade.utils.facet` (triangle :yref:`Facet`) and :yref:`yade.utils.wall` (infinite axes-aligned plane :yref:`Wall`) geometries are typically used to define boundaries. For instance, a "floor" for the simulation can be created like this:
 
.. ipython::

	Yade [1]: O.bodies.append(utils.wall(-1,axis=2))

There are other conveinence functions (like :yref:`yade.utils.facetBox` for creating closed or open rectangular box, or family of :yref:`yade.ymport` functions)

Look inside
-----------

The simulation can be inspected in several ways. All data can be accessed from python directly:

.. ipython::
	:okexcept:

	Yade [1]: len(O.bodies)

	Yade [1]: O.bodies[10].shape.radius   # radius of body #10 (will give error if not sphere, since only spheres have radius defined)

	Yade [1]: O.bodies[12].state.pos      # position of body #12

Besides that, Yade says this at startup (the line preceding the command-line)::

	[[ ^L clears screen, ^U kills line. F12 controller, F11 3d view, F10 both, F9 generator, F8 plot. ]]

:guilabel:`Controller`
	Pressing ``F12`` brings up a window for controlling the simulation. Although typically no human intervention is done in large simulations (which run "headless", without any graphical interaction), it can be handy in small examples. There are basic information on the simulation (will be used later).
:guilabel:`3d view`
	The 3d view can be opened with F11 (or by clicking on button in the *Controller* -- see below). There is a number of keyboard shortcuts to manipulate it (press ``h`` to get basic help), and it can be moved, rotated and zoomed using mouse.  Display-related settings can be set in the "Display" tab of the controller (such as whether particles are drawn).
:guilabel:`Inspector`
	*Inspector* is opened by clicking on the appropriate button in the *Controller*. It shows (and updates) internal data of the current simulation. In particular, one can have a look at engines, particles (*Bodies*) and interactions (*Interactions*). Clicking at each of the attribute names links to the appropriate section in the documentation.


Engines
-------

Engines define processes undertaken by particles. As we know from the theoretical introduction, the sequence of engines is called *simulation loop*. Let us define a simple interaction loop:

.. ipython::
	
	@suppress
	Yade [1]: O.reset()

	Yade [1]: O.engines=[                   # newlines and indentations are not important until the brace is closed
	   ...:    ForceResetter(),
	   ...:    InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
	   ...:    InteractionLoop(           # dtto for the parenthesis here
	   ...:        [Ig2_Sphere_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom()],
	   ...:        [Ip2_FrictMat_FrictMat_FrictPhys()],
	   ...:        [Law2_ScGeom_FrictPhys_CundallStrack()]
	   ...:    ),
	   ...:    NewtonIntegrator(damping=.2,label='newtonCustomLabel')      # define a label newtonCustomLabel under which we can access this engine easily
	   ...: ]
	   ...:

	Yade [1]: O.engines

	Yade [1]: O.engines[-1]==newtonCustomLabel    # is it the same object?

	Yade [1]: newtonCustomLabel.damping

Instead of typing everything into the command-line, one can describe simulation in a file (*script*) and then run yade with that file as an argument. We will therefore no longer show the command-line unless necessary; instead, only the script part will be shown. Like this::

	O.engines=[                   # newlines and indentations are not important until the brace is closed
		ForceResetter(),
		InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
		InteractionLoop(           # dtto for the parenthesis here
			 [Ig2_Sphere_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom()],
			 [Ip2_FrictMat_FrictMat_FrictPhys()],
			 [Law2_ScGeom_FrictPhys_CundallStrack()]
		),
		GravityEngine(gravity=(0,0,-9.81)),                    # 9.81 is the gravity acceleration, and we say that
		NewtonIntegrator(damping=.2,label='newtonCustomLabel') # define a label under which we can access this engine easily
	]

Besides engines being run, it is likewise important to define how often they will run. Some engines can run only sometimes (we will see this later), while most of them will run always; the time between two successive runs of engines is *timestep* ($\Dt$). There is a mathematical limit on the timestep value, called *critical timestep*, which is computed from properties of particles. Since there is a function for that, we can just set timestep using :yref:`yade.utils.PWaveTimeStep`::

	O.dt=utils.PWaveTimeStep()

Each time when the simulation loop finishes, time ``O.time`` is advanced by the timestep ``O.dt``:

.. ipython::

	Yade [1]: O.dt=0.01

	Yade [1]: O.time

	Yade [1]: O.step()

	Yade [1]: O.time

For experimenting with a single simulations, it is handy to save it to memory; this can be achieved, once everything is defined, with::

	O.saveTmp()



