# -*- encoding=utf-8 -*-
from __future__ import print_function
#########################################################################################################################################################################
# Authors: 	Remi Chassagne, remi.chassagne@univ-grenoble-alpes.fr
#		Raphael Maurin, raphael.maurin@imft.fr
# 03/09/2021
#
# Same as sedimentTransportExampleCoupling with two particle size, to study size segregation.
# Impose a given number of layers of smaller particles at the top of the bed and let it evolve under bedload transport at imposed fluid shear stress
# One recovers the logarithmic descent found by Chassagne et al (2020), JFM 895 A30
#
# Data can be saved (hdf5 format), to be further plotted with the file postProcess.py
#
############################################################################################################################################################################

#Import libraries
from builtins import range

#Import libraries
from yade import pack, plot
import math
import random as rand
import numpy as np
import h5py

##
## Main parameters of the simulation
##

#Particles
diameterPart1 = 6e-3  #Diameter of the particles initially at the bottom, in m
diameterPart2 = 4e-3  #Diameter of the particles initially at the top, in m
densPart = 2500  #density of the particles, in kg/m3
phiPartMax = 0.61  #Value of the dense packing solid volume fraction, dimensionless
restitCoef = 0.5  #Restitution coefficient of the particles, dimensionless
partFrictAngle = atan(0.4)  #friction angle of the particles, in radian

#Fluid
densFluid = 1000.  #Density of the fluid, in kg/m^3
kinematicViscoFluid = 1e-6  #kinematic viscosity of the fluid, in m^2/s
waterDepth = 9.  #Water depth, in diamete PHF:3.6 vs 3.1 pour comparer au cas mono qui ajoute 0.5Dl
dtFluid = 1e-5  #Time step for the fluid resolution, in s
fluidResolPeriod = 1e-2  #Time between two fluid resolution, in s

#Configuration: inclined channel
slope = 0.05  #Inclination angle of the channel slope in radian
lengthCell = 10  #Streamwise length of the periodic cell, in diameter
widthCell = 10  #Spanwise length of the periodic cell, in diameter
Nlayer = 10.  #Total nb of layer of particle, in 1st class particle diameter
Nlayer2 = 1.  #Nb of 2nd class particles, in 2nd class particle diameter
Nlayer1 = Nlayer - Nlayer2 * diameterPart2 / diameterPart1
print("Number of layer of large particles : {}".format(Nlayer1))
fluidHeight = (0.5 + Nlayer1) * diameterPart1 + Nlayer2 * diameterPart2 + waterDepth * diameterPart1  #Height of the flow from the bottom of the sample, in m

saveData = 1  #If put to 1, at each execution of function measure() save the sediment transport rate, fluid velocity, solid volume fraction and velocity profiles for post-processing
dtSave = 0.5  #Time step for saving results
endTime = 300  #Time simulated (in seconds)
#Save data details
if saveData == 1:  #If saveData option is activated, requires a folder data
	scriptPath = os.path.abspath(os.path.dirname(sys.argv[-1]))  #Path where the script is stored

##
## Secondary parameters of the simulation
##

expoDrag = 3.1  # Richardson Zaki exponent for the hindrance function of the drag force applied to the particles

#Discretization of the sample in ndimz wall-normal (z) steps of size dz, between the bottom of the channel and the position of the water free-surface. Should be equal to the length of the imposed fluid profile. Mesh used for HydroForceEngine.
dz = diameterPart1 / 30  # Fluid discretization step in the wall-normal direction
ndimz = int(1 + fluidHeight / dz)  # Number of cells in the height
dz = fluidHeight / (1.0 * (ndimz - 1))  #Correct the value of dz due to int() in ndimz
print("Nb cells = {}, dz={}".format(ndimz, dz))

# Initialization of the main vectors
vxFluid = np.zeros(ndimz)  # Vertical fluid velocity profile: u^f = u_x^f(z) e_x, with x the streamwise direction and z the wall-normal
phiPart = np.zeros(ndimz)  # Vertical particle volume fraction profile
vxPart = np.zeros(ndimz)  # Vertical average particle velocity profile
phiPart1 = np.zeros(ndimz)  # Vertical volume fraction profile of 1st particle class
phiPart2 = np.zeros(ndimz)  # Vertical volume fraction profile of 2nd particle class
vxPart1 = np.zeros(ndimz)  # Vertical velocity profile of 1st particle class
vxPart2 = np.zeros(ndimz)  # Vertical velocity profile of 2nd particle class

#Geometrical configuration, define useful quantities
height = 5 * fluidHeight  #heigth of the periodic cell, in m (bigger than the fluid height to take into particles jumping above the latter)
length = lengthCell * diameterPart1  #length of the stream, in m
width = widthCell * diameterPart1  #width of the stream, in m
groundPosition = height / 4.0  #Definition of the position of the ground, in m
gravityVector = Vector3(9.81 * sin(slope), 0.0, -9.81 * cos(slope))  #Gravity vector to consider a channel inclined with slope angle 'slope'

#Particles contact law/material parameters
maxPressure = (densPart - densFluid) * phiPartMax * (Nlayer1 * diameterPart1 +
                                                     Nlayer2 * diameterPart2) * abs(gravityVector[2])  #Estimated max particle pressure from the static load
normalStiffness = maxPressure * diameterPart1 * 1e4  #Evaluate the minimal normal stiffness to be in the rigid particle limit (cf Roux and Combe 2002)
youngMod1 = normalStiffness / diameterPart1  #Young modulus of the 1st particle class from the stiffness wanted.
youngMod2 = normalStiffness / diameterPart2  #Young modulus of the 2nd particle class from the stiffness wanted.
poissonRatio = 0.5  #poisson's ratio of the particles. Classical values, does not have much influence
O.materials.append(ViscElMat(en=restitCoef, et=0., young=youngMod1, poisson=poissonRatio, density=densPart, frictionAngle=partFrictAngle, label='Mat1'))
O.materials.append(ViscElMat(en=restitCoef, et=0., young=youngMod2, poisson=poissonRatio, density=densPart, frictionAngle=partFrictAngle, label='Mat2'))

########################
## FRAMEWORK CREATION ##
########################

#Definition of the semi-periodic cell
O.periodic = True
O.cell.setBox(length, width, height)

# Reference walls: build two planes at the ground and free-surface to have a reference for the eyes in the 3D view
lowPlane = box(center=(length / 2.0, width / 2.0, groundPosition), extents=(200, 200, 0), fixed=True, wire=False, color=(0., 1., 0.), material='Mat1')
WaterSurface = box(
        center=(length / 2.0, width / 2.0, groundPosition + fluidHeight),
        extents=(2000, width / 2.0, 0),
        fixed=True,
        wire=False,
        color=(0, 0, 1),
        material='Mat1',
        mask=0
)
O.bodies.append([lowPlane, WaterSurface])  #add to simulation

# Regular arrangement of spheres sticked at the bottom with random height
L = range(0, int(length / (diameterPart1)))  #The length is divided in particle diameter
W = range(0, int(width / (diameterPart1)))  #The width is divided in particle diameter

for x in L:  #loop creating a set of sphere sticked at the bottom with a (uniform) random altitude comprised between 0.5 (diameter/12) and 5.5mm (11diameter/12) with steps of 0.5mm. The repartition along z is made around groundPosition.
	for y in W:
		n = rand.randrange(
		        0, 12, 1
		) / 12.0 * diameterPart1  #Define a number between 0 and 11/12 diameter with steps of 1/12 diameter (0.5mm in the experiment)
		O.bodies.append(
		        sphere(
		                (x * diameterPart1, y * diameterPart1, groundPosition - 11 * diameterPart1 / 12.0 / 2.0 + n),
		                diameterPart1 / 2.,
		                color=(0, 0, 0),
		                fixed=True,
		                material='Mat1'
		        )
		)

#Create a loose cloud of particle inside the cell
partCloud1 = pack.SpherePack()
partCloud2 = pack.SpherePack()
partVolume1 = pi / 6. * pow(diameterPart1, 3)  #Volume of a particle 1
partVolume2 = pi / 6. * pow(diameterPart2, 3)  #Volume of a particle 2
partNumber1 = int(Nlayer1 * phiPartMax * diameterPart1 * length * width / partVolume1)  #Number of beads to obtain Nlayer1 layers of 1st class particles
partNumber2 = int(Nlayer2 * phiPartMax * diameterPart2 * length * width / partVolume2)  #Number of beads to obtain Nlayer2 layers of 2nd class particles
#Define the deposition height considering that the packing realised by make cloud is 0.2
depositionHeight1 = max(Nlayer1 * phiPartMax / 0.2 * diameterPart1, diameterPart1)  #Consider that the packing realised by make cloud is 0.2
depositionHeight2 = max(Nlayer2 * phiPartMax / 0.2 * diameterPart2, diameterPart2)  #Consider that the packing realised by make cloud is 0.2
partCloud1.makeCloud(
        minCorner=(0, 0., groundPosition + diameterPart1 + 1e-4),
        maxCorner=(length, width, groundPosition + depositionHeight1),
        rRelFuzz=0.,
        rMean=diameterPart1 / 2.0,
        num=partNumber1
)
partCloud2.makeCloud(
        minCorner=(0, 0., groundPosition + depositionHeight1),
        maxCorner=(length, width, groundPosition + depositionHeight1 + depositionHeight2),
        rRelFuzz=0.,
        rMean=diameterPart2 / 2.0,
        num=partNumber2
)
partCloud1.toSimulation(material='Mat1')  #Send this packing to simulation with material Mat1
partCloud2.toSimulation(material='Mat2', color=(0, 0, 1))  #Send this packing to simulation with material Mat2
#Evaluate the deposition time considering the free-fall time of the highest particle to the ground
depoTime = sqrt(fluidHeight * 2 / abs(gravityVector[2]))

# Collect the ids of the spheres which are dynamic to add a fluid force through HydroForceEngines
idApplyForce = []
for b in O.bodies:
	if isinstance(b.shape, Sphere) and b.dynamic:
		idApplyForce += [b.id]

#########################
#### SIMULATION LOOP#####
#########################

O.engines = [
        # Reset the forces
        ForceResetter(),
        # Detect the potential contacts
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb(), Bo1_Facet_Aabb(), Bo1_Box_Aabb()], label='contactDetection', allowBiggerThanPeriod=True),
        # Calculate the different interactions
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_ViscElMat_ViscElMat_ViscElPhys()], [Law2_ScGeom_ViscElPhys_Basic()],
                label='interactionLoop'
        ),
        #Apply an hydrodynamic force to the particles
        HydroForceEngine(
                densFluid=densFluid,
                viscoDyn=kinematicViscoFluid * densFluid,
                zRef=groundPosition,
                gravity=gravityVector,
                deltaZ=dz,
                expoRZ=expoDrag,
                nCell=ndimz,
                vCell=length * width * dz,
                ids=idApplyForce,
                label='hydroEngine',
                phiMax=phiPartMax,
                dead=True
        ),
        #Solve the fluid volume-averaged 1D momentum balance, RANS 1D
        PyRunner(command='fluidModel()', virtPeriod=fluidResolPeriod, label='fluidRes', dead=True),
        #Measurement, output files
        PyRunner(command='measure()', virtPeriod=dtSave, label='measurement', dead=True),
        # Check if the packing is stabilized, if yes activate the hydro force on the grains and the slope.
        PyRunner(command='gravityDeposition(depoTime)', virtPeriod=0.01, label='gravDepo'),
        #GlobalStiffnessTimeStepper, determine the time step
        GlobalStiffnessTimeStepper(defaultDt=1e-4, viscEl=True, timestepSafetyCoefficient=0.7, label='GSTS'),
        # Integrate the equation and calculate the new position/velocities...
        NewtonIntegrator(damping=0.2, gravity=gravityVector, label='newtonIntegr')
]
#save the initial configuration to be able to recharge the simulation starting configuration easily
O.saveTmp()
#run
O.run()

####################################################################################################################################
####################################################  FUNCTION DEFINITION  #########################################################
####################################################################################################################################

#Fluid coupling parametrization
hydroEngine.initialization()  #Initialize the variables
hydroEngine.enableMultiClassAverage = True  #Activate the averaging over different class/size of particles


######                                                                     ######
### LET THE TIME FOR THE GRAVITY DEPOSITION AND ACTIVATE THE FLUID AT THE END ###
######                                                                     ######
def gravityDeposition(lim):
	if O.time < lim:
		return
	else:
		print('\n Gravity deposition finished, apply fluid forces !\n')
		newtonIntegr.damping = 0.0  # Set the artificial numerical damping to zero
		gravDepo.dead = True  # Remove the present engine for the following
		hydroEngine.dead = False  # Activate the HydroForceEngine
		hydroEngine.vxFluid = vxFluid  # Send the fluid velocity vector used to apply the drag fluid force on particles in HydroForceEngine (see c++ code)
		measurement.dead = False  # Activate the measure() PyRunner
		fluidRes.dead = False  # Activate the 1D fluid resolution
		hydroEngine.averageProfile()  #Evaluate the solid volume fraction, velocity and drag, necessary for the fluid resolution.
		hydroEngine.fluidResolution(1., dtFluid)  #Initialize the fluid resolution, run the fluid resolution for 1s

	return


#######                       ########
###         FLUID RESOLUTION       ###
#######                       ########
def fluidModel():
	global vxFluid, taufsi
	#Evaluate the average vx,vy,vz,phi,drag profiles and store it in hydroEngine, to prepare the fluid resolution
	hydroEngine.averageProfile()
	#Fluid resolution
	hydroEngine.fluidResolution(fluidResolPeriod, dtFluid)  #Solve the fluid momentum balance for a time of fluidResolPeriod s with a time step dtFluid
	#update the fluid velocity for later save
	vxFluid = np.array(hydroEngine.vxFluid)


#######               ########
###         OUTPUT         ###
#######               ########
#Initialization
qsMean = 0  #Mean dimensionless sediment transport rate
fileNumber = 0
zAxis = np.zeros(ndimz)  #z scale, in diameter
for i in range(0, ndimz):  #z scale used for the possible plot at the end
	zAxis[i] = i * dz / diameterPart1


# Averaging/Save
def measure():
	global qsMean, vxPart, vxPart1, vxPart2, vyPart, vyPart1, vyPart2, vzPart, vzPart1, vzPart2, phiPart, phiPart1, phiPart2, turbulentViscosity, turbStress, fileNumber, stress1, stress2
	#Evaluate the average depth profile of streamwise, spanwise and wall-normal particle velocity, particle volume fraction (and drag force for coupling with RANS fluid resolution), and store it in hydroEngine variables vxPart, phiPart, vyPart, vzPart, averageDrag.
	hydroEngine.averageProfile()
	#Extract the calculated vector. They can be saved and plotted afterwards.
	vxPart = np.array(hydroEngine.vPart)[:, 0]
	vyPart = np.array(hydroEngine.vPart)[:, 1]
	vzPart = np.array(hydroEngine.vPart)[:, 2]
	vxPart1 = np.array(hydroEngine.multiVxPart)[0]
	vxPart2 = np.array(hydroEngine.multiVxPart)[1]
	vyPart1 = np.array(hydroEngine.multiVyPart)[0]
	vyPart2 = np.array(hydroEngine.multiVyPart)[1]
	vzPart1 = np.array(hydroEngine.multiVzPart)[0]
	vzPart2 = np.array(hydroEngine.multiVzPart)[1]
	phiPart = np.array(hydroEngine.phiPart)
	phiPart1 = np.array(hydroEngine.multiPhiPart[0])
	phiPart2 = np.array(hydroEngine.multiPhiPart[1])
	turbulentViscosity = np.array(hydroEngine.turbulentViscosity)
	turbStress = 1. / densFluid * np.array(hydroEngine.ReynoldStresses)

	#Function to get stress profile per class, not implemented for the moment. TO DO for Remi.
	#stress1 = getStressProfile_perClass(diameterPart1/2, dz*width*length, ndimz, dz, groundPosition, vxPart, vyPart, vzPart)[0]
	#stress2 = getStressProfile_perClass(diameterPart2/2, dz*width*length, ndimz, dz, groundPosition, vxPart, vyPart, vzPart)[0]
	#for i in range(len(stress1)):
	#   stress1[i] = np.array(stress1[i]).reshape(1,9)[0]
	#   stress2[i] = np.array(stress2[i]).reshape(1,9)[0]

	#Evaluate the dimensionless sediment transport rate for information
	qsMean = sum(phiPart * vxPart) * dz / sqrt((densPart / densFluid - 1) * abs(gravityVector[2]) * pow(diameterPart1, 3))

	#Condition to stop the simulation after endTime seconds
	if O.time >= endTime:
		print('\n End of the simulation, simulated {0}s as required !\n '.format(endTime))
		O.pause()

	#Evaluate the Shields number from the maximum of the Reynolds stresses evaluated in the fluid resolution
	dummy = np.where(phiPart[:ndimz] < 0.8 * phiPartMax)
	phi1Sum = sum(phiPart1[dummy])
	phi2Sum = sum(phiPart2[dummy])
	diameterMix = (diameterPart1 * phi1Sum + diameterPart2 * phi2Sum) / (phi1Sum + phi2Sum)
	shieldsNumber = max(hydroEngine.ReynoldStresses) / ((densPart - densFluid) * diameterMix * abs(gravityVector[2]))
	print('Shields number', shieldsNumber)

	if saveData == 1:  #Save data for postprocessing
		nameFile = scriptPath + '/data.hdf5'  # Name of the file that will be saved
		globalParam = [
		        'qsMean', 'phiPart', 'phiPart1', 'phiPart2', 'vxPart', 'vxPart1', 'vxPart2', 'vyPart', 'vyPart1', 'vyPart2', 'vzPart', 'vzPart1',
		        'vzPart2', 'vxFluid', 'zAxis', 'turbulentViscosity', 'turbStress'
		]  #, 'stress1', 'stress2']   # Variables to save
		Save(nameFile, globalParam)  #Save
		fileNumber += 1  #Increment the file number


##########################################
##########################################


#Function to save global variables in a python file which can be re-executed for post-processing
def Save(filePathName, globalVariables):
	h5File = h5py.File(filePathName, 'a')
	grp = h5File.create_group(str(globals()['fileNumber']))
	for varName in globalVariables:
		dataSet = grp.create_dataset(name=varName, data=globals()[varName])
	h5File.close()
