from yade import pack
from bf import evalClump

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
                clump_broken = evalClump(b.id, radius_ratio, tension_strength, compressive_strength, outer_predicate=outer_predicate,
                                         max_scale=max_scale, initial_packing_scale=initial_packing_scale, grow_radius=grow_radius, relative_gap=relative_gap)
                if clump_broken and pause_after_first_break:
                    O.pause()
                    O.dt = time_step_sf*PWaveTimeStep()
    for b in O.bodies:  # another loop just for coloring
        if isinstance(b.shape, Sphere) and b.iterBorn > 1 and b.isStandalone:
            b.shape.color = new_spheres_color
        elif isinstance(b.shape, Sphere) and b.iterBorn > 1 and b.isClumpMember:
            b.shape.color = new_clumps_color


########################
# SIMULATION
########################
# contstants
pause_after_first_break = False
time_step_sf = 0.8
strain_rate = 0.01
old_clumps_color = (0.1, 0.3, 0.5)
new_clumps_color = (0.5, 0.3, 0.1)
old_spheres_color = (0.7, 0.2, 0.7)
new_spheres_color = (0.2, 0.7, 0.2)
tension_strength = 5.0e6
compressive_strength = 5.0e6
radius_ratio = 5
wei_V0 = 1e-6
wei_P = 0.5
young = 10e9
density = 1000

# params related to subparticles packing
relative_gap = 0
grow_radius = 1
max_scale = 4
initial_packing_scale = 1.5

# MATERIAL
mat = FrictMat(label='grain')
# the second material just to test whether the particle's material will becorrectly applied to the sub-particle material (I deliberately give the same label to both materials)
mat2 = FrictMat(label='grain')
wall_mat = FrictMat(label='wall')
mat.young = young
mat2.young = 2*young
mat.density = density
wall_mat.young = 10*young
O.materials.append([mat, mat2, wall_mat])

# clumps to be broken
# two spheres
sphere_list_two = [[[0, 4e-2, 3.7e-3], 6e-3], [[0, 4e-2, -2e-3], 8e-3]]
# four spheres
sphere_list_four = [[[-5e-3, -1e-3, 0], 9.9e-3], [[5e-3, -
                                                   1e-3, 0], 8e-3], [[0, 1e-2, 0], 1e-2], [[0, 2e-2, 0], 7e-3]]

O.bodies.appendClumped([sphere(s[0], s[1], color=old_clumps_color, material=mat)
                        for s in sphere_list_two], discretization=20)
O.bodies.appendClumped([sphere(s[0], s[1], color=old_clumps_color, material=mat2)
                        for s in sphere_list_four], discretization=20)

# additional "neighbour" sphere to check splitting algorithm
O.bodies.append(sphere((2e-2, 1e-2, 0), 0.7e-2,
                       color=old_spheres_color, material=mat))

# walls
z_min, z_max = aabbExtrema()[0][2], aabbExtrema()[1][2]

if (yade.math.usesHP() or yade.libVersions.getArchitecture() == 'ppc64el'):
    print("This is HP, we do not support it! Or architecture is pp64el. Skip")
else:
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

    errors = 0

    if (len(O.bodies) != 11):
        print("Number of bodies before the test is not 11 but " + len(O.bodies))
        errors += 1

    if pause_after_first_break:
        O.run(wait=True)
        v.saveSnapshot('example_clump_compression-radius_ratio_{:.2f}-grow_radius_{:.2f}-relative_gap_{:.2f}.jpg'.format(
            radius_ratio, grow_radius, relative_gap))

    else:
        O.run(200000, True)

    if (len(O.bodies) != 93):
        print("Number of bodies before the test is not 93 but " + len(O.bodies))
        errors += 1

    if (errors):
        raise YadeCheckError("Test is failed")
