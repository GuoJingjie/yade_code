#########################################################################################################################################################################
# Author: Raphael Maurin, raphael.maurin@imft.fr
# 12/07/2021
#
# Script to validate the spring-dashpot law (Law2_ScGeom_ViscElPhys_Basic) and the lubrication force calculation from HydroForceEngine.
#
# Consider a simple case of two spheres colliding with each other without gravity, each with a given opposite velocity of magnitude v/2
# Measure the ratio between the relative particles velocity after collision and the relative particles velocity before collision
# Without lubrication or at Stokes > 10^3, this represents the dry restitution coefficient and one should recover the imposed value
#
# All the quantities are made dimensionless
#
#########################################################################################################################################################################
from yade import plot
import numpy as np
effectiveRestitutionApproach = True
lubricationForceApproach = False

#
## Parameters of the configuration
#

## Input parameters
# Particles
normalRestit = 0.9  # Imposed dry restitution coefficient
kappa = 1e-6  # Stiffness dimensionless number, ratio between the particle stiffness and the granular pressure. Should be greater to 1e-4 to be in the rigid grain limit.
partMass = 1.  # Dimensionless mass
diameter = 1.  # Dimensionless diameters of the spheres
relativeVel = 1.  # Dimensionless initial relative velocity
delta = 1e-3 * diameter  # Dimensionless initial relative position
roughnessRatio = 1e-2  # Roughness particle scale, dimensionless

#Fluid
densityRatio = 100.  # Particle to fluid density ratio
Stokes = 100.  # Stokes number associated with the collision
effectiveRestit = normalRestit * (1. + 1. / Stokes * log(2. * roughnessRatio))

## Calculated parameters
# Particles
partDensity = partMass / (pi / 6. * pow(diameter, 3.))  # Particle density from its mass
O.dt = kappa * diameter / relativeVel  # Time step in order to have a maximum overlap of kappa diameter
tContact = O.dt * 300.  # Contact time necessary to have a good resolution
stiffnessContact = partMass / 2. / tContact * (
        log(effectiveRestit)**2 + pi**2
)  #See Schwager & Poschel 2007 on spring-dashpot contact law, estimated according contact stiffness

# Fluid
viscoDyn = partDensity * relativeVel * diameter / Stokes  # Fluid viscosity calculated to match the stokes number required
fluidDensity = partDensity / densityRatio  # Fluid density calculated to match the density ratio

gravityVector = Vector3(0., 0., 0.)
#
## Creation of the framework
#

#Material
if effectiveRestitutionApproach:
	O.materials.append(
	        ViscElMat(
	                en=normalRestit,
	                et=1.,
	                frictionAngle=0.4,
	                young=stiffnessContact / diameter,
	                poisson=0.5,
	                density=partDensity,
	                lubrication=True,
	                viscoDyn=viscoDyn,
	                roughnessScale=roughnessRatio * diameter,
	                label='Mat1'
	        )
	)
else:
	O.materials.append(
	        ViscElMat(en=normalRestit, et=1., frictionAngle=0.4, young=stiffnessContact / diameter, poisson=0.5, density=partDensity, label='Mat1')
	)

#Create a free sphere at the origin with velocity in the direction of the other particle
O.bodies.append(sphere((0, 0, 0), diameter / 2.0, material='Mat1'))
O.bodies[0].state.vel[2] = relativeVel / 2.

#Create a free sphere to send on the other sphere/the wall with a velocity v
O.bodies.append(sphere((0, 0, diameter + delta), diameter / 2.0, material='Mat1'))
O.bodies[-1].state.vel[2] = -relativeVel / 2.

#
## Simulation loop
#

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb(), Bo1_Facet_Aabb(), Bo1_Box_Aabb()], label='contactDetection'),
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom(),
                 Ig2_Wall_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()], [Ip2_ViscElMat_ViscElMat_ViscElPhys()], [Law2_ScGeom_ViscElPhys_Basic()]
        ),
        HydroForceEngine(
                densFluid=fluidDensity,
                viscoDyn=viscoDyn,
                zRef=0.,
                gravity=gravityVector,
                deltaZ=diameter + delta,
                nCell=1,
                vCell=1,
                ids=[0, 1],
                label='hydroEngine',
                dead=True
        ),
        NewtonIntegrator(damping=0.0, gravity=gravityVector),
        PyRunner(command='Plot()', iterPeriod=100, label='plotFunction')
]

if lubricationForceApproach == True:
	hydroEngine.dead = False
	hydroEngine.lubrication = True
	hydroEngine.roughnessPartScale = roughnessRatio * diameter
	hydroEngine.initialization()

O.saveTmp()
O.run()


def Plot():
	plot.addData(ratioVel=abs(O.bodies[1].state.vel[2] - O.bodies[0].state.vel[2]) / relativeVel, dimensionlessTime=O.time / tContact)
	if O.time > delta / relativeVel + tContact * 500:
		print("Effective restitution coefficient measured:", abs(O.bodies[1].state.vel[2] - O.bodies[0].state.vel[2]) / relativeVel)
		print("Effective restitution coefficient expected:", effectiveRestit)
		O.pause()


#Plot the dimensionless sediment transport rate as a function of time during the simulation
plot.plots = {'dimensionlessTime': ('ratioVel')}
plot.plot()
