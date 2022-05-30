from bf import stressTensor, checkFailure, replaceSphere, evalClump
from yade import pack, plot
import numpy as np
import os
sys.path.append('.')
########################
# FUNCTIONS
########################


def stateUpdate(sphere_id, radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P, wei_m):
    if O.bodies[sphere_id] != None:
        effort = checkFailure(
            O.bodies[sphere_id], tension_strength, compressive_strength, wei_V0, wei_P, wei_m)
        if effort >= 1:
            outer_predicate = pack.inCylinder(w1.state.pos, w2.state.pos, 1)
            replaceSphere(sphere_id, radius_ratio=radius_ratio, relative_gap=relative_gap,
                          grow_radius=grow_radius, outer_predicate=outer_predicate, max_scale=max_scale)
            # refresh time step
            O.dt = time_step_sf*PWaveTimeStep()
        for b in O.bodies:  # another loop just for coloring
            # if not b == None:
            if isinstance(b.shape, Sphere) and b.iterBorn > 1:
                b.shape.color = new_spheres_color


def addPlotData():
    if O.bodies[sphere_id] == None:  # pause after particle breakage
        O.pause()
        if runningInBatch():  # if running in batch - save to file
            try:
                # create a new file, fails if file exists
                results = open(main_dir+'strength-uniaxial.txt', 'x')
                results.mode = 'w'
                results.write(
                    'batch_description,sphere_radius,wei_P,wei_m,max_stress,max_force \n')
                results.close()

            except:
                pass
            results = open(main_dir+'strength-uniaxial.txt',
                           'a')  # open file in append mode
            results_line = O.tags['description']+','+str(sphere_radius)+','+str(wei_P)+','+str(
                wei_m)+','+str(plot.data['sigma_z'][-1])+','+str(plot.data['f'][-1])+'\n'
            results.write(results_line)
            results.close()
    else:
        f1, f2 = [O.forces.f(i)[2] for i in [w1_id, w2_id]]
        f2 *= -1
        # average force
        f = .5*(f1+f2)
        # displacement of the moving wall
        dspl = 2*w2.state.displ()[2]
        # stress - according to the stress tensor
        sigma = stressTensor(O.bodies[sphere_id], stress_correction=True)
        sigma_z = sigma[2, 2]
        # quasi stress - simplistic evaluation: force/sphere crossection
        quasi_stress = f/(np.pi*sphere_radius**2)
        # store values
        yade.plot.addData(
            t=O.time,
            i=O.iter,
            displacement=dspl,
            sigma_z=sigma_z,
            quasi_stress=quasi_stress,
            f=f,
        )


########################
# SIMULATION
########################
# params
# default parameters or from table
readParamsFromTable(noTableOk=True,  # unknownOk=True,
                    time_step_sf=0.8,
                    strain_rate=0.001,

                    old_spheres_color=(0.7, 0.2, 0.7),
                    new_spheres_color=(0.2, 0.7, 0.2),
                    tension_strength=1.0e6,
                    compressive_strength=5.0e6,
                    sphere_radius=1e-2,
                    radius_ratio=3,
                    wei_V0=5.236e-7,  # 4.189e-9 is volume of sphere of radius 1 mm, 5.236e-7 is for radius 5 mm
                    wei_P=0.63,
                    wei_m=3,
                    main_dir='uniaxial_sphere_results/',
                    #simName = 'aaa',
                    # subparticles packing
                    relative_gap=0,
                    grow_radius=1,
                    max_scale=5.,
                    )
from yade.params.table import *
if runningInBatch():
    if not(os.path.exists(main_dir)):
        os.makedirs(main_dir)

# MATERIAL
mat = FrictMat(label='grain')
wall_mat = FrictMat(label='wall')
wall_mat.young = 10*mat.young
O.materials.append([mat, wall_mat])


# sphere to check splitting algorithm
sphere_id = O.bodies.append(
    sphere((0, 0, 0), sphere_radius, color=old_spheres_color, material=mat))

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
    PyRunner(command="stateUpdate(sphere_id,radius_ratio, tension_strength, compressive_strength, wei_V0, wei_P, wei_m )", iterPeriod=100),
    PyRunner(command="addPlotData()", iterPeriod=100)
]

O.dt = time_step_sf*PWaveTimeStep()


# plot.plots={'displacement':('sigma_z')}
# plot.plot()
O.run()
waitIfBatch()
