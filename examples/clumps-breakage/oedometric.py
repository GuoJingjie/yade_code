from yade import pack, plot, export
import numpy as np
import pandas as pd
import os
from bf import stressTensor, checkFailure, replaceSphere, evalClump

########################
# FUNCTIONS
########################

# functions for only for example simulation


def clump_psd(density, skip_Ids=[], percentage_result=True):
    """
    Finds psd of clumps and standalone spheres in simulation. The simplest version (the only so far) takes the mass of the clump, density of the material, and computes equivalent radius.
    skip_Ids - Ids of the clumps that should be omitted in simulation.
    density - one density is assumed for all the clumps since the clump has no density and only members does. Density needs to be passed by user.
    Function returns diameters and cumulated mass 'as is' (not exactly matching sieves.
    percentage_result - if true - the result is returned in percent instead of actual mass.
    """
    df = pd.DataFrame(columns=['r', 'm'])
    for b in O.bodies:
        if (isinstance(b.shape, Clump) or (isinstance(b.shape, Sphere) and b.isStandalone)) and not(b.id in skip_Ids):
            m = b.state.mass
            eq_radius = (3*m/(4*np.pi*density))**(1/3)
            df = df.append(pd.DataFrame(
                {'r': [eq_radius], 'm': [m]}), ignore_index=True)
    df_sorted = df.sort_values(by=['r'])
    df_grouped = df_sorted.groupby('r').sum().reset_index()
    df_grouped['D'] = 2*df_grouped['r']
    df_grouped['cum_mass'] = df_grouped['m'].cumsum()
    if percentage_result:
        df_grouped['cum_mass'] = 100 * \
            df_grouped['cum_mass']/max(df_grouped['cum_mass'])

    return df_grouped

# FUNCTIONS FOR PREDICATE PREPARING


def gts_cylinder_from_facet_cylinder(center, radius, height, segmentsNumber=10, result='pred'):
    """
    I used source code ( https://yade-dem.org/doc/_modules/yade/geom.html#facetCylinderConeGenerator )
    to retrieve the points created for cylinder generator.

    In the part 2, I create two polylines and use yade.pack.sweptPolylines2gtsSurface to create GTS surface for predicate exactly equal to the facet cylinder.
    Based on the option 'result' it can return 'pred', 'pts' or 'surface'.
    """
    #######################
    radius
    angleRange = (0, 2*math.pi)
    anglesInRad = numpy.linspace(
        angleRange[0], angleRange[1], segmentsNumber+1, endpoint=True)

    PTop = []
    PBottom = []

    for i in anglesInRad:
        XTop = radius*math.cos(i)
        YTop = radius*math.sin(i)
        PTop.append(Vector3(XTop, YTop, center[2]+height/2))

        XBottom = radius*math.cos(i)
        YBottom = radius*math.sin(i)
        PBottom.append(Vector3(XBottom, YBottom, center[2]-height/2))

    # part 2
    pts = [PBottom, PTop]
    surface = pack.sweptPolylines2gtsSurface(
        pts=pts, capStart=True, capEnd=True)
    surface.cleanup(1e-6)
    assert surface.is_closed()
    if result == 'pred':
        pred = pack.inGtsSurface(surface)
        return pred
    elif result == 'pts':
        return pts
    elif result == 'surface':
        return surface


# end of functions for predicate preparing

def initial_clump_coloring():
    """
    Random, uniform coloring of clumps.
    """
    for b in O.bodies:
        if isinstance(b.shape, Clump):
            color = randomColor()
            members_ids = b.shape.members.keys()
            for member_id in members_ids:
                O.bodies[member_id].shape.color = color


def clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P):
    # prepare cylinder predicate
    global predicate_fault_iters
    global outer_predicate
    global count_broken_particles
    try:
        outer_predicate = gts_cylinder_from_facet_cylinder(center=(
            0.0, 0.0, w2.state.pos[2]/2), radius=cyl_radius, height=w2.state.pos[2], segmentsNumber=segmentsNumber)
    except:
        predicate_fault_iters += [O.iter]
    for i in range(len(O.bodies)):  # We can't iterate directly on the O.bodies, because I would get "segmentation fault core dumped error" (the loop changes O. bodies list)
        b = O.bodies[i]
        if not b == None:
            if isinstance(b.shape, Clump):
                clump_broken = evalClump(b.id, radius_ratio, tension_strength, compressive_strength, outer_predicate=outer_predicate,
                                         max_scale=max_scale, grow_radius=grow_radius,  relative_gap=relative_gap, discretization=discretization)
                # refresh time step
                if clump_broken:
                    O.dt = time_step_sf*PWaveTimeStep()
                    count_broken_particles += 1
            elif isinstance(b.shape, Sphere) and b.isStandalone:
                evalSphere(b.id, radius_ratio, tension_strength, compressive_strength, wei_V0,
                           wei_P, outer_predicate, grow_radius=grow_radius,  relative_gap=relative_gap)
    for b in O.bodies:  # another loop just for coloring
        # if not b == None:
        if isinstance(b.shape, Sphere) and b.iterBorn > 1 and b.isStandalone:
            b.shape.color = new_spheres_color
        elif isinstance(b.shape, Sphere) and b.iterBorn > 1 and b.isClumpMember:
            b.shape.color = new_clumps_color


def evalSphere(sphere_id, radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P, outer_predicate, grow_radius, relative_gap):
    global count_broken_particles
    if O.bodies[sphere_id].shape.radius > min_radius_to_break:
        sigma = stressTensor(O.bodies[sphere_id], stress_correction=True)
        effort = checkFailure(
            O.bodies[sphere_id], tension_strength, compressive_strength, wei_V0, wei_P, wei_m=3)
        if effort >= 1:
            replaceSphere(sphere_id, radius_ratio=radius_ratio, outer_predicate=outer_predicate,
                          max_scale=max_scale, grow_radius=grow_radius, relative_gap=relative_gap)
            count_broken_particles += 1
            # refresh time step
            O.dt = time_step_sf*PWaveTimeStep()


########################
# SIMULATION
########################
# contstants
simulation_prefix = 'only_clumps_'
global predicate_fault_iters
global outer_predicate
global count_broken_particles
predicate_fault_iters = []
count_broken_particles = 0
# simulation
time_step_sf = 0.8
strain_rate = 0.01
new_clumps_color = (0.5, 0.3, 0.1)
new_spheres_color = (0.2, 0.7, 0.2)
# material properties
tension_strength = 40.0e6
compressive_strength = 200.0e6
radius_ratio = 2.0
wei_V0 = 5.236e-7
wei_P = 0.6
young = 10e9
density = 2600

# subparticles packing
relative_gap = 0
grow_radius = 1.0
max_scale = 5.
# geometry
zMin = 0
cyl_radius = 0.1
segmentsNumber = 20
packing_height = 1.0
# do not evaluate standalone spheres smaller than this radius
min_radius_to_break = 0.009
grain_mean_radius = 0.015
radius_RelFuzz = 0.2
clump_templates_list = [[[1, 0.9, 1], [[0, 0, 0], [0, 0.5, 0.5], [0, 0, 1]]],
                        [[1, 1, 1, 1, 1, 1], [[0, 0, -0.8], [0, 0, 0.8],
                                              [-0.8, 0, 0], [0.8, 0, 0], [0, -0.8, 0], [0, 0.8, 0]]],
                        [[1, 1, 1, 1, 1, 1], [[0, 0, 0], [0, 0, 0.7], [0.7, 0.35, 0], [1, 0, 0], [1, 1, 0], [1.7, 0.35, 0]]]]
clump_fractions = [0.2, 0.4, 0.4]
discretization = 20
# MATERIAL
mat = FrictMat(label='grain')
wall_mat = FrictMat(label='wall')
mat.young = young
mat.density = density
wall_mat.young = 10*young
O.materials.append([mat, wall_mat])

# cylinder from facets
# cylinder
cylinder = O.bodies.append(geom.facetCylinder(center=(0.0, 0.0, 0.0+packing_height/2), radius=cyl_radius,
                                              height=packing_height, wallMask=4, color=(1, 0, 0), material=wall_mat, segmentsNumber=segmentsNumber))


# PACKING
sp = pack.SpherePack()
sp.makeCloud((-cyl_radius, -cyl_radius, 0), (cyl_radius, cyl_radius,
                                             packing_height), num=300, rMean=grain_mean_radius, rRelFuzz=radius_RelFuzz)

filtered_sp = pack.filterSpherePack(predicate=gts_cylinder_from_facet_cylinder(center=(0.0, 0.0, 0.0+packing_height/2),
                                                                               radius=cyl_radius*0.99, height=packing_height, segmentsNumber=segmentsNumber), spherePack=sp, returnSpherePack=True)
filtered_sp.toSimulation(material=O.materials['grain'])
# replace spheres by clumps
clump_templates = []
for ct in clump_templates_list:
    clump_templates.append(clumpTemplate(ct[0], ct[1]))
O.bodies.replaceByClumps(clump_templates, clump_fractions,
                         discretization=discretization)
initial_clump_coloring()

# WALLS
z_min, z_max = aabbExtrema()[0][2], aabbExtrema()[1][2]

w1 = utils.wall(z_min, axis=2, sense=1, material=wall_mat)
w2 = utils.wall(z_max, axis=2, sense=-1, material=wall_mat)
w1_id, w2_id = O.bodies.append([w1, w2])

# ENGINES
O.engines = [
    ForceResetter(),
    InsertionSortCollider(
        [Bo1_Sphere_Aabb(), Bo1_Wall_Aabb(), Bo1_Facet_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Wall_Sphere_ScGeom(),
         Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()],
    ),
    NewtonIntegrator(damping=0.2, gravity=(0, 0, -9.81)),
]

O.dt = time_step_sf*PWaveTimeStep()
O.trackEnergy = True


O.run(int(2.5e5), wait=True)

# second part of the simulation
# prepare directories for storingresults
directories = ['output/oedometric/vtk/',
               'output/oedometric/psd/', 'output/oedometric/plots/']
for directory in directories:
    if not(os.path.exists(directory)):
        os.makedirs(directory)
# shift w2 closer to the grains
_, new_z_max = aabbExtrema()[0][2], aabbExtrema()[1][2]
shiftBodies([w2.id], (0, 0, new_z_max-z_max))
setRefSe3()

# predicate preparation (outer boundary for the new particles)
outer_predicate = gts_cylinder_from_facet_cylinder(center=(
    0.0, 0.0, w2.state.pos[2]/2), radius=cyl_radius, height=w2.state.pos[2], segmentsNumber=segmentsNumber)

w2.state.vel = (0, 0, -strain_rate*abs(z_max-z_min))


def count_particles():
    count = 0
    for b in O.bodies:
        if isinstance(b.shape, Clump):
            count += 1
        elif isinstance(b.shape, Sphere) and b.isStandalone:
            count += 1
    return count


disp_increment = 0.005  # PSD will be stored for each disp incremet
global last_disp
last_disp = 0.
global PSD
PSD = clump_psd(density=density)
PSD['disp'] = np.repeat(last_disp, len(PSD))

A = np.pi*cyl_radius**2  # cylinder crossection


def dataMine():
    global last_disp
    global PSD
    global count_broken_particles
    wall_force = O.forces.f(w2.id)[2]
    wall_disp = w2.state.displ()[2]
    N = count_particles()
    particles_mass = 0
    for i in range(len(O.bodies)):
        b = O.bodies[i]
        if b != None:
            if isinstance(b.shape, Clump) or (isinstance(b.shape, Sphere) and b.isStandalone):
                particles_mass += b.state.mass
    particles_volume = particles_mass/density
    sample_volume = A*(w2.state.pos[2] - w1.state.pos[2])
    porosity = 1.-particles_volume/sample_volume

    plot.addData(i=O.iter, t=O.time, F=wall_force, d=wall_disp, N=N, broken_particles=count_broken_particles,
                 porosity=porosity, particles_mass=particles_mass, **O.energy)
    # storing PSD information
    if wall_disp > last_disp+disp_increment:
        last_disp += disp_increment
        PSD_temp = clump_psd(density=density)
        PSD_temp['disp'] = np.repeat(wall_disp, len(PSD_temp))
        PSD = PSD.append(PSD_temp, ignore_index=True)
        PSD.to_csv('output/oedometric/psd/'+simulation_prefix+'.csv')
        # save plot data temporarily and after termination condition
        plot.saveDataTxt('output/oedometric/plots/' +
                         simulation_prefix+'_plot_data.txt.bz2')
    # termination condition
    if wall_disp > 0.10*(new_z_max-z_min):
        O.pause()
        plot.saveDataTxt('output/oedometric/plots/' +
                         simulation_prefix+'_plot_data.txt.bz2')


exporter = export.VTKExporter('output/oedometric/vtk/'+simulation_prefix)


def exp_vtk():
    exporter.exportSpheres(what=dict(color='b.shape.color', cId='b.clumpId'))
    exporter.exportFacets()
    inters = []  # prepare interactions for export
    for ii in O.interactions:
        b1 = O.bodies[ii.id1]
        b2 = O.bodies[ii.id2]
        # export only interactions between spheres
        if isinstance(b1.shape, Sphere) and isinstance(b2.shape, Sphere):
            # ommit interactions between the same clump members
            if b1.clumpId == b2.clumpId and not(b1.clumpId == -1):
                continue
            else:
                inters += [(ii.id1, ii.id2)]
    exporter.exportInteractions(ids=inters, what=dict(
        forceN='i.phys.normalForce.norm()'))


O.engines = O.engines + [PyRunner(command="clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P)", iterPeriod=10000),
                         PyRunner(command="dataMine()", virtPeriod=0.01),
                         PyRunner(command="exp_vtk()", virtPeriod=0.02)]


O.run()
