#this example demonstrates the difference between including gravity in the numerical damping versus not
#
# two particles are created, one with numerical damping applied, and one without
# The particles start at the same vertical position with an initial upward velocity in free ball, and then eventually fall
# down under gravity to impact a bottom wall.
#
# If the below variable is set to true (default behavior), then gravity is damped, so the two particles have different behaviors
# both under free fall and when in contact with the bottom wall.
#
# If the below variable is set to false, then gravity is not damped, so the two particles have identical motion under free fall
# but different behavior when contacting the bottom wall, as the contact forces are damped for one but not the other.

dampGrav = True

###############################################################################3

import matplotlib.pyplot as pyplot
from yade import qt, plot
qt.View()  #open the controlling and visualization interfaces

box_x = 0.05
box_y = 0.05
box_z = 0.05

particle_dia = 0.005

mn = Vector3(0, box_y, 0)
mx = Vector3(box_x, 2 * box_y, box_z)  #corners of the initial packing
thick = 2 * particle_dia  # the thickness of the walls

bigmx = (mx[0], 3 * mx[1], mx[2])

O.materials.append(FrictMat(young=2e8, poisson=0.3, frictionAngle=0, density=8000, label='spheres'))

ball1 = sphere((mx[0] / 2, 2 * mx[1], mx[2] / 2), particle_dia / 2, material='spheres')
ball1Id = O.bodies.append(ball1)

ball2 = sphere((mx[0] / 2, 2 * mx[1], mx[2] / 2 + 2 * particle_dia), particle_dia / 2, material='spheres')
ball2Id = O.bodies.append(ball2)
O.bodies[ball2Id].state.isDamped = False

O.bodies[ball1Id].state.vel[1] = 1
O.bodies[ball2Id].state.vel[1] = 1

walls = utils.aabbWalls([mn, bigmx], thickness=thick, oversizeFactor=1.5, material='spheres')
wallIds = O.bodies.append(walls)

#turn on gravity and let it settle
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()]),
        InteractionLoop([Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()]),
        NewtonIntegrator(gravity=(0, -9.81, 0), damping=0.2, dampGravity=dampGrav),
        PyRunner(command='addPlotData()', iterPeriod=2),
]
O.dt = .2 * PWaveTimeStep()


def addPlotData():
	plot.addData(time=O.time, DampedParticlePosition=O.bodies[ball1Id].state.pos[1], UndampedParticlePosition=O.bodies[ball2Id].state.pos[1])


plot.plots = {'time': ('DampedParticlePosition', 'UndampedParticlePosition')}
plot.plot()

#O.run(-1, True)
