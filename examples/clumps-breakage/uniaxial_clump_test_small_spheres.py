from yade import qt
from yade import pack
import numpy as np
from yade import plot
from bf import stressTensor, checkFailure, replaceSphere, evalClump

########################
# FUNCTIONS
########################


def clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P):
    for i in range(len(O.bodies)):
        b = O.bodies[i]
        if not b == None:
            if isinstance(b.shape, Clump):
                outer_predicate = pack.inCylinder(
                    w1.state.pos, w2.state.pos, 1)
                evalClump(b.id, radius_ratio, tension_strength, compressive_strength,
                          outer_predicate=outer_predicate, discretization=discretization)
                # refresh time step
                O.dt = time_step_sf*PWaveTimeStep()


########################
# SIMULATION
########################
# contstants
time_step_sf = 0.8
strain_rate = 0.01
old_spheres_color = (0.7, 0.2, 0.7)
new_spheres_color = (0.2, 0.7, 0.2)
tension_strength = 1.0e6
compressive_strength = 5.0e6
radius_ratio = 10
wei_V0 = 1e-6
wei_P = 0.9
young = 10e9
discretization = 10
# MATERIAL
mat = FrictMat(label='grain')
wall_mat = FrictMat(label='wall')
mat.young = young
wall_mat.young = 10*young
O.materials.append([mat, wall_mat])

# clump to be broken
# three spheres
sphere_list_with_small_sphere = [
    [[0, 0, 0], 10e-3], [[0, 10e-3, 0], 10e-3], [[0, 0, 9e-3], 1.1e-3]]

O.bodies.appendClumped([sphere(s[0], s[1], material=mat)
                        for s in sphere_list_with_small_sphere], discretization=discretization)
# walls
z_min, z_max = aabbExtrema()[0][2], aabbExtrema()[1][2]

w1 = utils.wall(z_min, axis=2, sense=1, material=wall_mat)
w2 = utils.wall(z_max, axis=2, sense=-1, material=wall_mat)
w1_id, w2_id = O.bodies.append([w1, w2])

w2.state.vel = (0, 0, -strain_rate*abs(z_max-z_min))

# engines
O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()],
    ),
    NewtonIntegrator(damping=0.2),
    PyRunner(command="clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P)", iterPeriod=1000)
]

O.dt = time_step_sf*PWaveTimeStep()


qt.Controller()
v = qt.View()
v.ortho = True
v.viewDir = (-1, 0, 0)
v.upVector = (0, 0, 1)
v.showEntireScene()

O.run()
