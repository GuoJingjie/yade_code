# -*- encoding=utf-8 -*-
# 2021 J.Kozicki, A.Gladky, K.Thoeni "Implementation of high-precision computation capabilities into the dynamic simulation framework YADE" (to be submitted)
# The script is tested with yade commit 7f645d64f (April 2021)
from yade import plot, math

# initialize all floating point variables with Real(arg) to avoid precision loss
from yade.math import Real
import sys
import csv
import lzma

if(math == sys.modules['math']):
    raise RuntimeError(
        "To avoid precision problems it is better to avoid the python math (double only) module and use yade.math (faster) or mpmath (slower) instead.")

### set parameters ###
L = Real('0.1')      # length [m]
n = 4                # number of nodes for the length [-]
r = L/100            # radius [m]
g = Real('9.81')     # gravity
inclination = math.radiansHP1(20)  # Initial inclination of rods [degrees]
color = [1, 0.5, 0]  # Define a color for bodies

O.dt = Real('1e-05') # time step
damp = Real('1e-1')  # damping. It is interesting to examine  damp = 0

totalSecDuration = 100
totalIterDuration = int(totalSecDuration/O.dt)
iterPyRunnerPeriod = int(1e-4/O.dt)

O.engines = [  # define engines, main functions for simulation
    ForceResetter(),
    PyRunner(command='addPlotData()', initRun=True,
             iterPeriod=iterPyRunnerPeriod, nDo=1),  # run this at the beginning to get Eref=Epot(t=0)
    InsertionSortCollider([Bo1_Sphere_Aabb()],
                          label='ISCollider', avoidSelfInteractionMask=True),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom6D()],
        [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
            setCohesionNow=True, setCohesionOnNewContacts=False)],
        [Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
    ),
    NewtonIntegrator(gravity=(0, -g, 0), damping=damp, label='newton'),
    PyRunner(command='addPlotData()', initRun=True,
             iterPeriod=iterPyRunnerPeriod, firstIterRun=0),
    PyRunner(command='getPrevVel()', initRun=True,
             iterPeriod=iterPyRunnerPeriod, firstIterRun=iterPyRunnerPeriod-1),
    PyRunner(command='saveResults()', initRun=True,
             iterPeriod=iterPyRunnerPeriod, firstIterRun=totalIterDuration-1)
]

# define material:
O.materials.append(CohFrictMat(young=1e5, poisson=0, density=1e1,
                               frictionAngle=math.radiansHP1(0), normalCohesion=1e7,
                               shearCohesion=1e7, momentRotationLaw=False, label='mat'))

# create spheres
nodeIds = []
for i in range(0, n):
    nodeIds.append(O.bodies.append(sphere([i*L/n*math.cos(inclination),
        i*L/n*math.sin(inclination), 0], r, wire=False, fixed=False,
        material='mat', color=color)))

# create rods
for i, j in zip(nodeIds[:-1], nodeIds[1:]):
    inter = createInteraction(i, j)
    inter.phys.unp = -(O.bodies[j].state.pos-O.bodies[i].state.pos).norm() + O.bodies[i].shape.radius+O.bodies[j].shape.radius

# Set a fixed node
O.bodies[0].dynamic = False

# Calculate angle between three points: a, b, c
def calculateAngle(a, b, c):
    ba = a - b
    bc = c - b
    cosineAngle = ba.dot(bc) / (ba.norm() * bc.norm())

    if (cosineAngle > 1):
        angle = math.acos(1)
    elif (cosineAngle < -1):
        angle = math.acos(-1)
    else:
        angle = math.acos(cosineAngle)
    return (math.degrees(angle))

# Calculate potential energy
def calculatePotentialEnergy():
    potEnergy = Real(0)
    for i in range(0, n):
        h1 = O.bodies[i].state.pos[1]  # mass is lumped into nodes
        h2 = L/n*i
        potEnergy += O.bodies[i].state.mass * g * (h1 + h2)
    return potEnergy

# Calculate elastic energy
def calculateElasticEnergy():  # normal force only
    elEnergy = Real(0)
    for i in O.interactions:
        elEnergy += i.phys.normalForce.squaredNorm()/(2 * i.phys.kn)
    return elEnergy


dataToSave = []

iter_prev = 0
pos_prev = {'1': [0, 0],
            '2': [0, 0],
            '3': [0, 0]}
vel_prev = {'1': [0, 0],
            '2': [0, 0],
            '3': [0, 0]}
ang1_prev = 0
ang2_prev = 0

kin_prev = 0
pot_prev = 0
elast_prev = 0

# Get velocities for the previous step


def getPrevVel():
    global iter_prev
    global pos_prev
    global vel_prev
    global ang1_prev
    global ang2_prev
    global kin_prev
    global pot_prev
    global elast_prev
    iter_prev = O.iter
    pos_prev = {'1': [O.bodies[1].state.pos[0], O.bodies[1].state.pos[1]],
                '2': [O.bodies[2].state.pos[0], O.bodies[2].state.pos[1]],
                '3': [O.bodies[3].state.pos[0], O.bodies[3].state.pos[1]]}
    vel_prev = {'1': [O.bodies[1].state.vel[0], O.bodies[1].state.vel[1]],
                '2': [O.bodies[2].state.vel[0], O.bodies[2].state.vel[1]],
                '3': [O.bodies[3].state.vel[0], O.bodies[3].state.vel[1]]}
    ang1_prev = calculateAngle(O.bodies[0].state.pos,
                               O.bodies[1].state.pos, O.bodies[2].state.pos)
    ang2_prev = calculateAngle(O.bodies[1].state.pos,
                               O.bodies[2].state.pos, O.bodies[3].state.pos)
    kin_prev = kineticEnergy() # not used in the paper, only used for cross-checking
    pot_prev = calculatePotentialEnergy()
    elast_prev = calculateElasticEnergy()

# Function to accumulate and add data to plot


def addPlotData():
    ang1 = calculateAngle(O.bodies[0].state.pos,
                          O.bodies[1].state.pos, O.bodies[2].state.pos)
    ang2 = calculateAngle(O.bodies[1].state.pos,
                          O.bodies[2].state.pos, O.bodies[3].state.pos)

    curTime = O.iter*O.dt
    kin = kineticEnergy()
    pot = calculatePotentialEnergy()
    elast = calculateElasticEnergy()

    dataCurrent  = {'time': curTime, 'iter_prev': iter_prev, 'iter': O.iter,
                    'pos1_prev_x': pos_prev['1'][0], 'pos1_prev_y': pos_prev['1'][1],
                    'pos2_prev_x': pos_prev['2'][0], 'pos2_prev_y': pos_prev['2'][1],
                    'pos3_prev_x': pos_prev['3'][0], 'pos3_prev_y': pos_prev['3'][1],
                    'pos1_x': O.bodies[1].state.pos[0], 'pos1_y': O.bodies[1].state.pos[1],
                    'pos2_x': O.bodies[2].state.pos[0], 'pos2_y': O.bodies[2].state.pos[1],
                    'pos3_x': O.bodies[3].state.pos[0], 'pos3_y': O.bodies[3].state.pos[1],
                    'vel1_prev_x': vel_prev['1'][0],    'vel1_prev_y': vel_prev['1'][1],
                    'vel2_prev_x': vel_prev['2'][0],    'vel2_prev_y': vel_prev['2'][1],
                    'vel3_prev_x': vel_prev['3'][0],    'vel3_prev_y': vel_prev['3'][1],
                    'vel1_x': O.bodies[1].state.vel[0], 'vel1_y': O.bodies[1].state.vel[1],
                    'vel2_x': O.bodies[2].state.vel[0], 'vel2_y': O.bodies[2].state.vel[1],
                    'vel3_x': O.bodies[3].state.vel[0], 'vel3_y': O.bodies[3].state.vel[1],
                    'ang1_prev': ang1_prev,
                    'ang2_prev': ang2_prev,
                    'ang1': ang1,
                    'ang2': ang2,
                    'kin_prev': kin_prev,
                    'pot_prev': pot_prev,
                    'elast_prev': elast_prev,
                    'kin': kin,
                    'pot': pot,
                    'elast': elast
                    }
    dataToSave.append(dataCurrent)
    plot.addData(
        ang1=ang1, ang2=ang2, time1=curTime, time2=curTime, kin=kin, pot=pot, elast=elast
    )

def infoString():
    return '_digits_'+str(yade.math.getDigits10(1))+'_dt_'+str(O.dt)

def saveResults():
    fname1 = 'data_angles'+infoString()
    fname2 = 'data_energy'+infoString()+'.csv.xz'
    print('\n'+'-'*70)
    print(int(totalIterDuration), "iterations finished. Saving results to:\n",fname1,"\n",fname2)
    print('-'*70)
    plot.saveGnuplot(fname1)
    with lzma.open(fname2, mode='wt') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(dataToSave[0].keys()))
        writer.writeheader()
        for i in dataToSave:
            writer.writerow(i)
    print('\n'+'-'*70)
    print('data saved')
    print('-'*70)

plot.plots = {'time1': ('ang1', 'ang2'), 'time2': ('kin', 'pot', 'elast')}
plot.plot()

O.saveTmp()

if(opts.nogui==False):
    from yade import qt
    qt.View()                    # create a GUI view
    Gl1_Sphere.stripes = True    # mark spheres with stripes
    rr = qt.Renderer()           # get instance of the renderer
    rr.intrAllWire = True        # draw wires
    rr.intrPhys = True           # draw the normal forces between the spheres.


print('-'*70)
print("Will now run",int(totalIterDuration),"iterations, which is",totalSecDuration,"seconds.")
print("Using time step = ",O.dt)
print('-'*70)

O.run(totalIterDuration)

