#This example demonstrates the use of HarmonicForceEngine
#A simple harmonic oscillator (spring-mass-damper) excited by an external force is demonstrated
#The mass of the oscillator is a single spherical particle, stiffness comes from contact with a wall,
#damping is provided by numerical damping in the Newton integrator and the external force is applied with HarmonicForceEngine
#with the values used below, the interaction stiffness for each wall is 100,000 N/m, (so total 200,000 N/m
#considering top & bottom walls), the ball mass is 0.00431 kg. This gives a natural frequency of 1084 Hz.
#The x, y, and z directions are excited at frequencies of 900 Hz (below resonance), 1084 (near resonance) and 1300 Hz
#(above resonance).  As expected the largest response is when the excitation frequency is near resonance.

import matplotlib.pyplot as pyplot
from yade import qt, plot
qt.View() #open the controlling and visualization interfaces

#box dimensions intentionally smaller than the particle diameter to create an interference fit (clearance would result in
#a nonlinear stiffness behavior)
box_x = 0.004 
box_y = 0.004
box_z = 0.004

particle_dia = 0.010

young2 = 2e7 
rho    = 8230 

mn = Vector3(-box_x, -box_y, -box_z)
mx = Vector3( box_x,  box_y,  box_z) #corners of the initial packing
thick = 1.1*particle_dia # the thickness of the walls
               

O.materials.append(FrictMat(young=young2, frictionAngle=0, density=rho ,   label='mat'))

        
ball = sphere( (0, 0 , 0 ) , particle_dia/2, material='mat' )
ballIds = O.bodies.append(ball)

ballWeight = O.bodies[ballIds].state.mass * 9.81

walls=utils.aabbWalls([mn,mx],thickness=thick,oversizeFactor=1.1,material='mat')        
wallIds=O.bodies.append(walls)


swing = HarmonicForceEngine( ids = [ballIds],  A = (ballWeight, ballWeight, ballWeight), f = (900, 1084 , 1300) )


O.engines=[
	ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()] ),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()] ),
        swing,
	NewtonIntegrator(gravity=(0,0 ,0),damping=0.05),
        PyRunner(command='addPlotData()', iterPeriod=1),
]

O.dt=.5*PWaveTimeStep()


def addPlotData():
        plot.addData(time=O.time ,posx = O.bodies[ballIds].state.pos[0] , posy = O.bodies[ballIds].state.pos[1], posz = O.bodies[ballIds].state.pos[2]   )
                

plot.plots={'time':('posx', 'posy', 'posz')}
plot.plot()


#O.run(-1, True) 





        
