#########################################################################################################################################################################
# Author: Raphael Maurin, raphael.maurin@imft.fr
# 27/08/2021
#
# Script to validate the lubrication component of the spring-dashpot law (Law2_ScGeom_ViscElPhys_Basic)
#
# Reproduce N times the simpleBallRebound script where two particles collides, varying the imposed effective restitution coefficient as a function of the stokes number
# Allows to reproduce curves of effective restitution coefficient as a function of the stokes number.
#
# All the quantities are made dimensionless
#
# THERE IS STILL A PROBLEM IN THE RESULT: AN UNEXPECTED CHANGE OF BEHAVIOR AT ST = 2e-2
#
#########################################################################################################################################################################
from yade import plot
import numpy as np

#
## Parameters of the configuration
#

stokesRange = np.logspace(
        0, 3, 200
)  #Stokes range explored for the collisions created. The length of stokes range correspond to the number of collisions that will be modeled.

## Input parameters
# Particles
normalRestit = 0.9  # Imposed dry restitution coefficient
kappa = 1e-6  # Stiffness dimensionless number, ratio between the particle stiffness and the granular pressure. Should be greater to 1e-4 to be in the rigid grain limit.
partMass = 1.  # Dimensionless mass
diameter = 1.  # Dimensionless diameters of the spheres
relativeVel = 1.  # Dimensionless initial relative velocity
delta = 1e-4 * diameter  # Dimensionless initial relative position
roughnessRatio = 1e-2  # Roughness particle scale, dimensionless

#Fluid
densityRatio = 1.  # Particle to fluid density ratio
Stokes = 100.  # Stokes number associated with the collision

## Calculated parameters
# Particles
partDensity = partMass / (pi / 6. * pow(diameter, 3.))  # Particle density from its mass
O.dt = kappa * diameter / relativeVel  # Time step in order to have a maximum overlap of kappa diameter
tContact = O.dt * 300.  # Contact time necessary to have a good resolution

# Fluid
fluidDensity = partDensity / densityRatio  # Fluid density calculated to match the density ratio

gravityVector = Vector3(0., 0., 0.)
#
## Creation of the framework
#

#Material
n = 0
effectiveRestit = np.zeros(len(stokesRange))
beadsCoupleID = np.zeros([len(stokesRange), 2], dtype=int)
for Stokes in stokesRange:
	nameMat = 'Mat' + str(n)
	viscoDyn = partDensity * relativeVel * diameter / Stokes
	effectiveRestit[n] = max(1e-3, normalRestit * (1. + 1. / Stokes * log(2. * roughnessRatio)))
	stiffnessContact = partMass / 2. / tContact * (
	        log(effectiveRestit[n])**2 + pi**2
	)  #See Schwager & Poschel 2007 on spring-dashpot contact law, estimated according contact stiffness
	O.materials.append(
	        ViscElMat(
	                en=normalRestit,
	                et=1.,
	                frictionAngle=0.5,
	                young=stiffnessContact / diameter,
	                poisson=0.5,
	                density=partDensity,
	                lubrication=True,
	                viscoDyn=viscoDyn,
	                roughnessScale=roughnessRatio * diameter,
	                label=nameMat
	        )
	)
	#Create a free sphere at the origin with velocity in the direction of the other particle
	beadsCoupleID[n, 0] = O.bodies.append(sphere((2 * n * diameter, 0, 0), diameter / 2.0, material=nameMat))
	beadsCoupleID[n, 1] = O.bodies.append(sphere((2 * n * diameter, 0, diameter + delta), diameter / 2.0, material=nameMat))
	for b in O.bodies:
		if b.id == beadsCoupleID[n, 0]:
			b.state.vel[2] = relativeVel / 2.
		if b.id == beadsCoupleID[n, 1]:
			b.state.vel[2] = -relativeVel / 2.
	n += 1

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
        NewtonIntegrator(damping=0.0, gravity=gravityVector),
        PyRunner(command='Meas()', iterPeriod=100)
]
O.saveTmp()
O.run()

effectiveRestitMeas = np.zeros(len(stokesRange))
stokesRangeMeas = np.zeros(len(stokesRange))


def Meas():
	if O.time > delta / relativeVel + tContact * 500:
		n = 0
		k = 0
		for b1 in O.bodies:
			for b2 in O.bodies:
				n = 0
				for idCouple in beadsCoupleID:
					if b1.id == idCouple[0] and b2.id == idCouple[1]:
						effectiveRestitMeas[k] = abs(b1.state.vel[2] - b2.state.vel[2]) / relativeVel
						stokesRangeMeas[k] = stokesRange[n]
						k += 1
					n += 1
		print('finished, execute Plot() function now to see the results')
		O.pause()


import matplotlib.pyplot as plt


def Plot():
	plt.figure()
	plt.semilogx(stokesRangeMeas, effectiveRestitMeas / normalRestit, 'ok', label='measured')
	plt.semilogx(stokesRange, effectiveRestit / normalRestit, '-r', label='imposed')
	plt.xlabel('Stokes number')
	plt.ylabel('Effective restitution coefficient/dry restitution coefficient')
	plt.show()
