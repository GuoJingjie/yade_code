# -*- coding: utf-8

# 3 Spheres have an initial velociries: 0, +5, -5
# Their positions and velocities are checked during free fall
# checkGravity.py                     Checks the correctness of NewtonIntegrator and GravityEngine
# checkGravityRungeKuttaCashKarp54.py Checks the correctness of RungeKuttaCashKarp54Integrator and GravityEngine

# Omega
from __future__ import print_function
from yade.math import toHP1, radiansHP1
o = Omega()

# PhysicalParameters
# In present test these parameters are in fact not used, because the spheres are only falling down.
# However for good measure let us use the full precision Real here also.
Density = toHP1(2400)
frictionAngle = radiansHP1(35)
sphereRadius = toHP1('0.05')
tc = toHP1('0.001')
en = toHP1('0.3')
et = toHP1('0.3')

O.dt = 0.02 * tc

sphereMat = O.materials.append(ViscElMat(density=Density, frictionAngle=frictionAngle, tc=tc, en=en, et=et))


v_down = toHP1(-5)
v_up = toHP1(5)

# Important  !!!! when using yade high-precision be careful because any *quickly* *typed* ad-hoc floating point number can "destroy" the calculations.
# Important  !!!! For example in yade-mpfr150 declaring      g = -9.81 in fact produces number which is not useful to
#            ↓↓↓↓ do a simulation with 150 decimal places :  g = -9.8100000000000004973799150320701301097869873046875
#            ↓↓↓↓ make sure to use, yade.math.toHP1(………), so that the numbers on python side use native yade rpecision.
g = yade.math.toHP1("-9.81")

tolerance = 1e-10
print("checkGravity.py : yade precision is ", yade.math.getDigits10(1), " decimal places. Will use error tolerance of: ", tolerance)

print("Note: for high precision calculations use yade.math.toHP1(...) see this:")
print("g = yade.math.toHP1('-9.81') ## produces: ", (yade.math.toHP1('-9.81')).__repr__())
print("g = -9.81                    ## produces: ", (-9.81 * yade.math.toHP1('1')).__repr__())
print("See https://yade-dem.org/doc/HighPrecisionReal.html#string-conversions for more info.")


def calcPos(v, t):
    return v * t + g * t * t / 2


id_0 = o.bodies.append(sphere((0.0, calcPos(0, O.dt / 2), 0), 0.2, material=sphereMat))  # The body has no initial vertical Velocity
id_down = o.bodies.append(sphere((1.0, calcPos(v_down, O.dt / 2), 0), 0.2, material=sphereMat))  # The body has an initial vertical Velocity -5
id_up = o.bodies.append(sphere((2.0, calcPos(v_up, O.dt / 2), 0), 0.2, material=sphereMat))  # The body has an initial vertical Velocity +5

# Instead of setting initial positions (above) to be declared at O.dt/2 one can also set initial velocities to be at -O.dt/2
# see https://gitlab.com/yade-dev/trunk/-/merge_requests/555#note_462560944
# Or use a different integrator, see checkGravityRungeKuttaCashKarp54.py
O.bodies[id_down].state.vel[1] = v_down
O.bodies[id_up].state.vel[1] = v_up

# Usual engines
o.engines = [
    ForceResetter(),
    InsertionSortCollider([
        Bo1_Sphere_Aabb(),
        Bo1_Facet_Aabb(),
    ]),
    InteractionLoop(
        [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
        [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
        [Law2_ScGeom_ViscElPhys_Basic()],
    ),
    NewtonIntegrator(damping=0, gravity=[0, g, 0]),
    PyRunner(command='checkPos()', iterPeriod=10000),
]


def checkPos():
    #print("Iter=%i, time=%.15f, dt=%.15f" % (O.iter, O.time, O.dt))
    if (abs((O.bodies[id_0].state.pos[1] - getCurrentPos(0)) / O.bodies[id_0].state.pos[1]) > tolerance):
        warningMessagePos(0, O.bodies[id_0].state.pos[1], getCurrentPos(0))

    if (abs((O.bodies[id_down].state.pos[1] - getCurrentPos(v_down)) / O.bodies[id_down].state.pos[1]) > tolerance):
        warningMessagePos(v_down, O.bodies[id_down].state.pos[1], getCurrentPos(v_down))

    if (abs((O.bodies[id_up].state.pos[1] - getCurrentPos(v_up)) / O.bodies[id_up].state.pos[1]) > tolerance):
        warningMessagePos(v_up, O.bodies[id_up].state.pos[1], getCurrentPos(v_up))

    if (abs((O.bodies[id_0].state.vel[1] - getCurrentVel(0)) / O.bodies[id_0].state.vel[1]) > tolerance):
        warningMessageVel(0, O.bodies[id_0].state.vel[1], getCurrentVel(0))

    if (abs((O.bodies[id_down].state.vel[1] - getCurrentVel(v_down)) / O.bodies[id_down].state.vel[1]) > tolerance):
        warningMessageVel(v_down, O.bodies[id_down].state.vel[1], getCurrentVel(v_down))

    if (abs((O.bodies[id_up].state.vel[1] - getCurrentVel(v_up)) / O.bodies[id_up].state.vel[1]) > tolerance):
        warningMessageVel(v_up, O.bodies[id_up].state.vel[1], getCurrentVel(v_up))


def getCurrentPos(inVel=0):
    t = O.time + O.dt * 1.5
    return inVel * t + g * t * t / 2


def getCurrentVel(inVel=0):
    t = O.time + O.dt
    return inVel + g * t


def warningMessagePos(inVel, y_pos, y_pos_need):
    raise YadeCheckError("The body with the initial velocity %.3f, has an y-position %.19f, but it should be %.19f. Iter=%i, time=%.15f, dt=%.15f" %
                         (inVel, y_pos, y_pos_need, O.iter, O.time, O.dt))


def warningMessageVel(inVel, y_vel, y_vel_need):
    raise YadeCheckError("The body with the initial velocity %.3f, has an y-velocity %.19f, but it should be %.19f. Iter=%i, time=%.15f, dt=%.15f" %
                         (inVel, y_vel, y_vel_need, O.iter, O.time, O.dt))


O.run(1000000, True)
