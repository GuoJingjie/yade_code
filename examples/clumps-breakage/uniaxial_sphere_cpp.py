from yade import qt
from bf import stressTensor, checkFailure, replaceSphere, evalClump
from yade import pack
import numpy as np
from yade import plot
import sys
sys.path.append('./')

########################
# FUNCTIONS
########################
# functions for only for example simulation


def stateUpdate(sphere_id, radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P):
    if O.bodies[sphere_id] != None:
        #mass = O.bodies[sphere_id].state.mass
        sigma = stressTensor(O.bodies[sphere_id], stress_correction=False)
        effort = checkFailure(
            O.bodies[sphere_id], tension_strength, compressive_strength, wei_V0, wei_P, wei_m=3)
        plot.addData(t=O.time, sigma_z=sigma[2, 2])
        if effort >= 1:
            # here I removed the option "contact_with_flat_bodies" but I used "outer_predicate" instead.
            # the outer predicate is a cylinder with a large radius (top and bottom position on the walls)
            outer_predicate = pack.inCylinder(w1.state.pos, w2.state.pos, 1)
            replaceSphere(sphere_id, radius_ratio=radius_ratio, relative_gap=relative_gap,
                          grow_radius=grow_radius, outer_predicate=outer_predicate, max_scale=max_scale)
            # refresh time step
            O.dt = time_step_sf*PWaveTimeStep()
            O.pause()
        for b in O.bodies:  # another loop for coloring
            # if not b == None:
            if isinstance(b.shape, Sphere) and b.iterBorn > 1:
                b.shape.color = new_spheres_color


########################
# SIMULATION
########################
# contstants
time_step_sf = 0.8
strain_rate = 0.01

old_spheres_color = (0.7, 0.2, 0.7)
new_spheres_color = (0.2, 0.7, 0.2)
tension_strength = 5.0e6
compressive_strength = 50.0e6
radius_ratio = 5
wei_V0 = 1e-6
wei_P = 0.5
young = 10e9

# subparticles packing
relative_gap = 0
grow_radius = 1
max_scale = 5.

# MATERIAL
mat = FrictMat(label='grain')
wall_mat = FrictMat(label='wall')
mat.young = young
wall_mat.young = 10*young
O.materials.append([mat, wall_mat])


# sphere to check splitting algorithm
O.bodies.append(sphere((0, 0, 0), 1e-2, color=old_spheres_color, material=mat))

# walls
z_min, z_max = aabbExtrema()[0][2], aabbExtrema()[1][2]

w1 = utils.wall(z_min, axis=2, sense=1, material=wall_mat)
w2 = utils.wall(z_max, axis=2, sense=-1, material=wall_mat)
w1_id, w2_id = O.bodies.append([w1, w2])

w2.state.vel = (0, -1, -strain_rate*abs(z_max-z_min))


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
    PyRunner(command="stateUpdate(0,radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P)", iterPeriod=1),
    ClumpBreakage(iterPeriod=1, stressCorrection=True, label='CB')
]

O.dt = 0.8*PWaveTimeStep()


qt.Controller()
v = qt.View()
v.ortho = True
v.viewDir = (0, 1, 0)
v.showEntireScene()
# move camera back a little bit
v.eyePosition = Vector3(0., 1.3*v.eyePosition[1], 0.)
O.step()
O.step()
# v.saveSnapshot('uniaxial_sphere_before_compression.jpg')

#O.run(wait = True)

# v.saveSnapshot('uniaxial_sphere_compression-radius_ratio_{:.2f}-grow_radius_{:.2f}-relative_gap_{:.2f}.jpg'.format(radius_ratio,grow_radius,relative_gap))
