# -*- encoding=utf-8 -*-
# 2021 J.Kozicki, A.Gladky, K.Thoeni "Implementation of high-precision computation capabilities into the dynamic simulation framework YADE" (to be submitted)
# The script is tested with yade_2021.01a
from yade import plot, qt, math

# initialize all floating point variables with Real(arg) to avoid precision loss
from yade.math import toHP1 as Real
# after release yade_2021.01a math.toHP1 has an alias Real, same with radiansHP1

### set parameters ###
L = Real('0.1')      # length [m]
n = 4                # number of nodes for the length [-]
r = L / 100            # radius [m]
g = Real('9.81')     # gravity
inclination = math.radiansHP1(20)  # Initial inclination of rods [degrees]
color = [1, 0.5, 0]  # Define a color for bodies

O.dt = Real('1e-05')  # time step
damp = Real('1e-1')  # damping. It is interesting to examine  damp = 0

O.engines = [  # define engines, main functions for simulation
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb()],
                          label='ISCollider', avoidSelfInteractionMask=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom6D()],
        [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
            setCohesionNow=True, setCohesionOnNewContacts=False)],
        [Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
    ),
    NewtonIntegrator(gravity=(0, -g, 0), damping=damp, label='newton'),
]

# define material:
O.materials.append(CohFrictMat(young=1e5, poisson=0, density=1e1,
                               frictionAngle=math.radiansHP1(0), normalCohesion=1e7,
                               shearCohesion=1e7, momentRotationLaw=False, label='mat'))

# create spheres
nodeIds = []
for i in range(0, n):
    nodeIds.append(O.bodies.append(sphere([i * L / n * math.cos(inclination),
                                           i * L / n * math.sin(inclination), 0], r, wire=False, fixed=False,
                                          material='mat', color=color)))

# create rods
for i, j in zip(nodeIds[:-1], nodeIds[1:]):
    inter = createInteraction(i, j)
    inter.phys.unp = -(O.bodies[j].state.pos - O.bodies[i].state.pos).norm() + \
        O.bodies[i].shape.radius + O.bodies[j].shape.radius

O.bodies[0].dynamic = False  # set a fixed upper node
qt.View()                    # create a GUI view
Gl1_Sphere.stripes = True    # mark spheres with stripes
rr = qt.Renderer()           # get instance of the renderer
rr.intrAllWire = True        # draw wires
rr.intrPhys = True           # draw the normal forces between the spheres.
