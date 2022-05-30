from bf import stressTensor, checkFailure, replaceSphere, evalClump
from yade import pack, plot, export, qt
import numpy as np
import os


# FUNCTIONS
def clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P):
    outer_predicate = pack.inCylinder(
        (0, 0, -1*sphere_radius), (0, 0, 1*sphere_radius), 100*sphere_radius)
    for i in range(len(O.bodies)):
        b = O.bodies[i]
        if b != None:
            if isinstance(b.shape, Clump) and b.id != saw_clump:
                evalClump(b.id, radius_ratio, tension_strength, compressive_strength,
                          wei_V0=wei_V0, wei_P=wei_P, wei_m=wei_m, outer_predicate=outer_predicate)

    # another loop for coloring and fix clumps
    for b in O.bodies:
        if not(b is None) and b.isClump and b.id != saw_clump:
            b.state.blockedDOFs = 'xyzXYZ'
            for m in b.shape.members.keys():
                O.bodies[m].shape.color = (0, 0, 1)

        if b.isStandalone:
            if remove_spheres:
                O.bodies.erase(b.id)
            else:
                b.shape.color = (1, 1, 0)
    # refresh time step
    O.dt = time_step_sf*PWaveTimeStep()


def addPlotData():
    torque = O.forces.t(saw_clump)[2]
    vertical_force = O.forces.f(saw_clump)[1]
    displ = -O.bodies[saw_clump].state.displ()[1]
    plot.addData(
        t=O.time,
        i=O.iter,
        displacement=displ,
        torque=torque,
        vertical_force=vertical_force,
    )
    if O.time > 10:
        O.pause()
        plot.saveDataTxt(main_dir+'force_plot_rot_speed_{:.2f}_weibull_V0_{:.9f}_compressive_{:.0f}.txt'.format(
            rotations_per_sec, wei_V0, compressive_strength))


########################
# SIMULATION
########################
# params
# default parameters or from table
readParamsFromTable(noTableOk=True,
                    time_step_sf=0.8,
                    tension_strength=1.0e6,
                    compressive_strength=2.5e6,
                    sphere_radius=1e-2,
                    radius_ratio=3,
                    wei_V0=5.236e-7,  # 5.236e-7 is for radius 5 mm
                    wei_P=0.63,
                    wei_m=3,
                    young=10e6,
                    main_dir='abbrasive_results/',
                    relative_gap=0,
                    grow_radius=1,
                    max_scale=5.0,
                    remove_spheres=True,
                    rotations_per_sec=10,
                    )
from yade.params.table import *

for directory in [main_dir, main_dir+'vtk/']:
    if not(os.path.exists(directory)):
        os.makedirs(directory)

# MATERIAL
mat = FrictMat(label='grain')
saw_mat = FrictMat(label='saw')
mat.young = young
saw_mat.young = young
O.materials.append([mat, saw_mat])

# BODIES
# abrasive wall
# prepare spheres
wall_clump_spheres = []
for x in np.arange(0, 20*sphere_radius, 2*sphere_radius):
    for y in np.arange(0, 10*sphere_radius, 0.5*sphere_radius):
        for z in [0]:
            wall_clump_spheres += [sphere((x, y, z),
                                          sphere_radius, material=mat, color=(0, 0, 1))]

# clump spheres
c_id, _ = O.bodies.appendClumped(wall_clump_spheres, discretization=20)
O.bodies[c_id].state.blockedDOFs = 'xyzXYZ'

# saw clump
center = (9*sphere_radius, 20.5*sphere_radius, 0)
saw_radius = 9*sphere_radius
# I multiplied it by some number in order to control distances between saw members.
number_of_spheres = int(np.pi*saw_radius/sphere_radius*0.75)
#ang_increment = 2*np.pi/number_of_spheres

saw_clump, _ = O.bodies.appendClumped([sphere((center[0]+saw_radius*np.cos(phi), center[1]+saw_radius*np.sin(
    phi), 0), sphere_radius, material=saw_mat, color=(1, 0, 0)) for phi in np.linspace(0, 2*np.pi, number_of_spheres)], discretization=20)

O.bodies[saw_clump].state.blockedDOFs = 'xyzXYZ'
O.bodies[saw_clump].state.angVel = (0, 0, 2*pi*saw_radius*rotations_per_sec)
O.bodies[saw_clump].state.vel = (0, -0.01, 0)


# engines
O.engines = [
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb()]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()],
    ),
    NewtonIntegrator(damping=0.2, label='newton'),
    PyRunner(command="clumpUpdate(radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P)", iterPeriod=1000),
]

O.dt = time_step_sf*PWaveTimeStep()

if not(runningInBatch()):
    remove_spheres = False
    newton.gravity = (0, -10, 0)
    exporter = export.VTKExporter(main_dir+'vtk/rot_speed_{:.2f}_weibull_V0_{:.9f}_compressive_{:.0f}'.format(
        rotations_per_sec, wei_V0, compressive_strength))
    O.engines += [PyRunner(command="exporter.exportSpheres(what = dict(color='b.shape.color[0]+b.shape.color[1]', cId='b.clumpId', born = 'b.timeBorn'))", virtPeriod=0.025),
                  PyRunner(command="addPlotData()", virtPeriod=0.001)]

    plot.plots = {'displacement': ('vertical_force')}
    plot.plot()

    qt.Controller()
    v = qt.View()
    v.ortho = True
    v.viewDir = (0, 0, 1)
    v.showEntireScene()
    # move camera back a little bit
    v.eyePosition = Vector3(
        v.eyePosition[0], v.eyePosition[1], 1.3*v.eyePosition[2])

setRefSe3()
O.run()
waitIfBatch()
