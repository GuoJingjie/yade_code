# -*- encoding=utf-8 -*-
from __future__ import print_function
#########################################################################################################################################################################
# Author: Raphael Maurin, raphael.maurin@imft.fr
# 25/08/2021
#
# Validation of the 1D fluid resolution from HydroForceEngine, perform 3 different cases: 
#	- 1D laminar fluid resolution without DEM coupling: reproduce a poiseuille flow and compares to analytical solution (inspired from fluidValidation/poiseuille.py)
#	- 1D turbulent fluid resolution without DEM coupling: reproduce a logarithmic profile and compares to analytical solution (inspired from fluidValidation/logProfile.py)
#	- fluid-DEM coupling: perform a coupled simulation and compare to reference data obtained with current Yade version (inspired from twoWayCoupling/sedimentTransportExample1DCoupling.py)
#
# The script send back a message confirming if the test have been passed or not, as a function of the comparison. 
#
############################################################################################################################################################################
import numpy as np
from math import *
#Function to evaluate the relative RMS error btwn two vectors
def errorRMS(ref,a):
	if len(ref)!= len(a): 
		print("errorRMS: you can only compare two vectors of the same size")
		exit()
	err = 0.
	N = 0.
	for i in range(len(ref)):
		err += sqrt(pow(ref[i]-a[i],2.))
		N+=1.
	return err/N/ref.mean()


##################
# Validation options
fullValidation = False	#Fluid + DEM coupling
fluidValidation = False	#Fluid only validation (short)
DEMcouplingValidation = True	#DEM coupling validation only (long, for comparison with another run/another yade version)



if fullValidation ==True or fluidValidation==True:
	##########################################################
	# POISEUILLE FLOW VALIDATION
	# see fluidValidations/poiseuille.py for details about the script

	[fluidHeight,dpdx,densFluid,kinematicViscoFluid,tfin,dtFluid,ndimz] = [1.,-0.10,1e3,1e-6,1e7,1e4,101]
	dz =  fluidHeight/float(ndimz-1)
	flResP = HydroForceEngine(densFluid=densFluid,viscoDyn=kinematicViscoFluid*densFluid,deltaZ=dz,nCell=ndimz,dpdx=dpdx,gravity = Vector3(0,0,-9.81),phiMax=0.61,iturbu=0,iusl=0,uTop=0.)
	flResP.initialization()
	flResP.fluidResolution(tfin,dtFluid)
	vxFluidP = np.array(flResP.vxFluid)

	# Re-create the fluid resolution mesh & evaluate the analytical solution
	[zScale,zScale[0],zScale[1]]=[np.zeros(ndimz+1),0.,0.5*dz]
	for i in range(2,ndimz):
		zScale[i]=zScale[i-1]+dz
	zScale[ndimz] = zScale[ndimz-1] + 0.5*dz
	analyticalVxFluidP=-dpdx*1/2./(kinematicViscoFluid*densFluid)*(fluidHeight*zScale-zScale*zScale) # Analytical solution of the fluid velocity
	errorP = errorRMS(analyticalVxFluidP,vxFluidP)
	print('\nPoiseuille: error wrt analytical solution',errorP)
	if errorP<0.0015: 
		print('Poiseuille profile reproduced with less than 0.15% error!\nVALIDATED\n')
	else:
		print('!!!!!!!!!!!!!!!!!!\n NOT VALIDATED: the error is larger than the one expected with classical version of YADE, be careful with the fluid resolution formulation, there might be a problem!\n!!!!!!!!!!!!!!!!!!!!!!')


	##########################################################
	# POISEUILLE FLOW VALIDATION
	# see fluidValidations/logProfile.py for details about the script

	## Analytical Solution
	[alpha,fluidHeight,tfin,dtFluid,kappa] = [0.5e-1,1e-1,1e3,2e-1,0.41]
	gravityVector = Vector3(9.81*sin(alpha),0,-9.81*cos(alpha))
	ustar = sqrt(abs(gravityVector[0])*fluidHeight)
	dz = .1*kinematicViscoFluid/ustar
	[ndimz, ztrans] = [int(round(fluidHeight/dz)),11.3]
	[utrans, z_analytic, analyticalVxFluidLP] = [ztrans*ustar, np.linspace(0,1e0,ndimz)*fluidHeight, np.zeros(ndimz)]
	zplus=z_analytic*ustar/kinematicViscoFluid
	for i in range(ndimz):
		if (zplus[i]<=ztrans): 
			analyticalVxFluidLP[i]=ustar*zplus[i]	#Linear fluid profile in the viscous sub layer
		else:
			analyticalVxFluidLP[i]=ustar/kappa*np.log(zplus[i]/ztrans)+utrans;	#Logarythmic fluid profile above
	#Fluid resolution
	ndimz-=1
	flResL = HydroForceEngine(densFluid=densFluid,viscoDyn=kinematicViscoFluid*densFluid,deltaZ=dz,nCell=ndimz,gravity=gravityVector,phiMax=0.61,iturbu=1,ilm=0,iusl=1,viscousSubLayer=1,kappa = kappa)
	flResL.initialization()
	flResL.fluidResolution(tfin,dtFluid)
	vxFluidLP = np.array(flResL.vxFluid)


	errorLP = errorRMS(analyticalVxFluidLP,vxFluidLP)
	print('\nLog profile: error wrt analytical solution',errorLP)
	if errorLP<0.035: 
		print('Log profile reproduced with less than 3.5% error! \nVALIDATED\n')
	else:
		print('!!!!!!!!!!!!!!!!!!\n NOT VALIDATED: the error is larger than the one expected with classical version of YADE, be careful with the fluid resolution formulation, there might be a problem!\n!!!!!!!!!!!!!!!!!!!!!!')



if fullValidation or DEMcouplingValidation==True:
	##########################################################
	# Fluid-DEM coupling validation
	# see twoWayCoupling/sedimentTransportExample_1DRANSCoupling.py for details about the script
	#Import libraries
	from yade import pack, plot
	import random as rand
	rand.seed(0)

	exactReproductionOnSameComputer = True	# To activate this option: run a simulation with the present script on the version of Yade you want and save the data. 
	#Then, run another simulation with a different version of Yade 
	dataBaseRef = False
	if os.path.exists('dataRef')==False: 
		print('Create a reference data base for exact comparison')
		dataBaseRef = True
		

	oldVersion = True	#Version of yade before 09/2021, requires some small adaptations
	if oldVersion==False: 
		from builtins import range

	##
	## Main parameters of the simulation
	##

	#Particles
	diameterPart = 6e-3	#Diameter of the particles, in m
	densPart = 2500		#density of the particles, in kg/m3
	phiPartMax = 0.61	#Value of the dense packing solid volume fraction, dimensionless
	restitCoef = 0.7	#Restitution coefficient of the particles, dimensionless
	partFrictAngle = atan(0.5)	#friction angle of the particles, in radian

	#Fluid
	densFluidPY = 1000.	#Density of the fluid, in kg/m^3
	kinematicViscoFluid = 1e-6	#kinematic viscosity of the fluid, in m^2/s
	waterDepth = 17.5	#Water depth, in diameter
	dtFluid = 1e-5 	#Time step for the fluid resolution, in s
	fluidResolPeriod = 1e-2	#Time between two fluid resolution, in s

	#Configuration: inclined channel
	slope = 0.05	#Inclination angle of the channel slope in radian
	lengthCell = 10	#Streamwise length of the periodic cell, in diameter
	widthCell = 10	#Spanwise length of the periodic cell, in diameter
	Nlayer = 10.	#nb of layer of particle, in diameter
	fluidHeight = (Nlayer+waterDepth)*diameterPart	#Height of the flow from the bottom of the sample, in m

	saveData = 0	#If put to 1, at each execution of function measure() save the sediment transport rate, fluid velocity, solid volume fraction and velocity profiles for post-processing
	endTime = 100	#Time simulated (in seconds)

	if exactReproductionOnSameComputer:
		endTime = 1
		saveData = 1

	##
	## Secondary parameters of the simulation
	##

	expoDrag_PY = 3.1	# Richardson Zaki exponent for the hindrance function of the drag force applied to the particles

	#Discretization of the sample in ndimz wall-normal (z) steps of size dz, between the bottom of the channel and the position of the water free-surface. Should be equal to the length of the imposed fluid profile. Mesh used for HydroForceEngine.
	ndimz = 300	#Number of cells in the height
	dz =  fluidHeight/(1.0*(ndimz-1))	# Fluid discretization step in the wall-normal direction	

	# Initialization of the main vectors
	vxFluidPY = np.zeros(ndimz+1)	# Vertical fluid velocity profile: u^f = u_x^f(z) e_x, with x the streamwise direction and z the wall-normal. Fluid velocity defined in between the mesh nodes and at the node at the two boundaries, i.e. at ndimz-1 + 2 location = ndimz+1.
	phiPartPY = np.zeros(ndimz-1)	# Vertical particle volume fraction profile, evaluated in between the cells, i.e. at ndimz-1 locations
	vxPartPY = np.zeros(ndimz-1)	# Vertical average particle velocity profile, evaluated in between the cells, i.e. at ndimz-1 locations

	#Geometrical configuration, define useful quantities
	height = 5*fluidHeight	#heigth of the periodic cell, in m (bigger than the fluid height to take into particles jumping above the latter)
	length = lengthCell*diameterPart #length of the stream, in m
	width = widthCell*diameterPart  #width of the stream, in m
	groundPosition = height/4.0 #Definition of the position of the ground, in m
	gravityVector = Vector3(9.81*sin(slope),0.0,-9.81*cos(slope)) #Gravity vector to consider a channel inclined with slope angle 'slope'

	#Particles contact law/material parameters
	maxPressure = (densPart-densFluidPY)*phiPartMax*Nlayer*diameterPart*abs(gravityVector[2]) #Estimated max particle pressure from the static load
	normalStiffness = maxPressure*diameterPart*1e4 #Evaluate the minimal normal stiffness to be in the rigid particle limit (cf Roux and Combe 2002)
	youngMod = normalStiffness/diameterPart	#Young modulus of the particles from the stiffness wanted.
	poissonRatio = 0.5	#poisson's ratio of the particles. Classical values, does not have much influence
	O.materials.append(ViscElMat(en=restitCoef, et=0., young=youngMod, poisson=poissonRatio, density=densPart, frictionAngle=partFrictAngle, label='Mat'))  


	########################
	## FRAMEWORK CREATION ##
	########################

	#Definition of the semi-periodic cell
	O.periodic = True 
	O.cell.setBox(length,width,height)

	# Reference walls: build two planes at the ground and free-surface to have a reference for the eyes in the 3D view
	lowPlane = box(center= (length/2.0, width/2.0,groundPosition),extents=(200,200,0),fixed=True,wire=False,color = (0.,1.,0.),material = 'Mat')
	WaterSurface = box(center= (length/2.0, width/2.0,groundPosition+fluidHeight),extents=(2000,width/2.0,0),fixed=True,wire=False,color = (0,0,1),material = 'Mat',mask = 0)
	O.bodies.append([lowPlane,WaterSurface]) #add to simulation


	# Regular arrangement of spheres sticked at the bottom with random height
	L = list(range(0,int(length/(diameterPart)))) #The length is divided in particle diameter
	W = list(range(0,int(width/(diameterPart)))) #The width is divided in particle diameter

	for x in L: #loop creating a set of sphere sticked at the bottom with a random altitude comprised between -0.5 and 0.5 diameter around groundPosition.
		for y in W:
			n =  rand.random()     #Define a number between 0 and 1
			O.bodies.append(sphere((x*diameterPart, y*diameterPart,groundPosition + (-0.5+n)*diameterPart),diameterPart/2.,color=(0,0,0),fixed = True,material = 'Mat'))

	#Create a loose cloud of particle inside the cell
	partCloud = pack.SpherePack()
	partVolume = pi/6.*pow(diameterPart,3) #Volume of a particle
	partNumber = int(Nlayer*phiPartMax*diameterPart*length*width/partVolume) #Volume of beads to obtain Nlayer layers of particles
	partCloud.makeCloud(minCorner=(0,0.,groundPosition+diameterPart),maxCorner=(length,width,groundPosition+fluidHeight),rRelFuzz=0., rMean=diameterPart/2.0, num = partNumber,seed = 0)
	partCloud.toSimulation(material='Mat') #Send this packing to simulation with material Mat
	#Evaluate the deposition time considering the free-fall time of the highest particle to the ground
	depoTime = sqrt(fluidHeight*2/abs(gravityVector[2]))

	# Collect the ids of the spheres which are dynamic to add a fluid force through HydroForceEngines
	idApplyForce = []
	for b in O.bodies: 
		if isinstance(b.shape,Sphere) and b.dynamic:
			idApplyForce+=[b.id]



	#########################
	#### SIMULATION LOOP#####
	#########################

	O.engines = [
		# Reset the forces
		ForceResetter(),
		# Detect the potential contacts
		InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Wall_Aabb(),Bo1_Facet_Aabb(),Bo1_Box_Aabb()],label='contactDetection',allowBiggerThanPeriod = True),
		# Calculate the different interactions
		InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()],
		[Ip2_ViscElMat_ViscElMat_ViscElPhys()],
		[Law2_ScGeom_ViscElPhys_Basic()]
		,label = 'interactionLoop'),				
		#Apply an hydrodynamic force to the particles
		HydroForceEngine(densFluid = densFluidPY,viscoDyn = kinematicViscoFluid*densFluidPY,zRef = groundPosition,gravity = gravityVector,deltaZ = dz,expoRZ = expoDrag_PY,lift = False,nCell = ndimz,vCell = length*width*dz,radiusPart=diameterPart/2.,phiMax = 0.61,ids = idApplyForce, label = 'hydroEngine', dead = True),
		#Solve the fluid volume-averaged 1D momentum balance, RANS 1D
		PyRunner(command = 'fluidModel()', virtPeriod = fluidResolPeriod, label = 'fluidRes', dead = True),
		#Measurement, output files
		PyRunner(command = 'measure()', virtPeriod = 0.1, label = 'measurement', dead = True),
		# Check if the packing is stabilized, if yes activate the hydro force on the grains and the slope.
		PyRunner(command='gravityDeposition(depoTime)',virtPeriod = 0.01,label = 'gravDepo'),
		#GlobalStiffnessTimeStepper, determine the time step
		GlobalStiffnessTimeStepper(defaultDt = 1e-4, viscEl = True,timestepSafetyCoefficient = 0.7,  label = 'GSTS'),
		# Integrate the equation and calculate the new position/velocities...
		NewtonIntegrator(damping=0.2, gravity=gravityVector, label='newtonIntegr')
		]
	#save the initial configuration to be able to recharge the simulation starting configuration easily
	O.saveTmp()
	#Initialize HydroForceEngine variables to zero (fluid velocity, fluctuations,...)
	if oldVersion==True: 
		hydroEngine.vxFluid = vxFluidPY
		hydroEngine.vxPart  = np.zeros(ndimz)
		hydroEngine.phiPart  = np.zeros(ndimz)
		hydroEngine.vFluctX = np.zeros(len(O.bodies))
		hydroEngine.vFluctY = np.zeros(len(O.bodies))
		hydroEngine.vFluctZ = np.zeros(len(O.bodies))
	else: 
		hydroEngine.initialization()

	#run
	O.run()

	####################################################################################################################################
	####################################################  FUNCTION DEFINITION  #########################################################
	####################################################################################################################################



	######								           ######
	### LET THE TIME FOR THE GRAVITY DEPOSITION AND ACTIVATE THE FLUID AT THE END ###
	######								           ######
	def gravityDeposition(lim):
		if O.time<lim : return
		else :
			print('\n Gravity deposition finished, apply fluid forces !\n')
			newtonIntegr.damping = 0.0	# Set the artificial numerical damping to zero
			gravDepo.dead = True	# Remove the present engine for the following
			measurement.dead = False	# Activate the measure() PyRunner
			fluidRes.dead = False		# Activate the 1D fluid resolution

			hydroEngine.dead = False	# Activate the HydroForceEngine
			hydroEngine.ReynoldStresses = np.ones(ndimz)*1e-4 # Send the simplified fluid Reynolds stresses Rxz/\rho^f used to account for the fluid velocity fluctuations in HydroForceEngine (see c++ code)
			hydroEngine.averageProfile()	#Evaluate the solid volume fraction, velocity and drag, necessary for the fluid resolution. 
			hydroEngine.fluidResolution(1.,dtFluid)	#Initialize the fluid resolution, run the fluid resolution for 1s
		return
	###############
	#########################################


	#######			      ########
	###	    FLUID RESOLUTION	   ###
	#######			      ########
	def fluidModel():
		global vxFluidPY
		#Evaluate the average vx,vy,vz,phi,drag profiles and store it in hydroEngine, to prepare the fluid resolution
		hydroEngine.averageProfile()
		#Fluid resolution
		hydroEngine.fluidResolution(fluidResolPeriod,dtFluid)	#Solve the fluid momentum balance for a time of fluidResolPeriod s with a time step dtFluid
		#update the fluid velocity for later save
		vxFluidPY = np.array(hydroEngine.vxFluid)
		

	#######		      ########
	###	    OUTPUT	   ###
	#######		      ########
	#Initialization
	qsMean = 0		#Mean dimensionless sediment transport rate
	zAxis = np.zeros(ndimz)	#z scale, in diameter
	for i in range(0,ndimz):#z scale used for the possible plot at the end
		zAxis[i] = i*dz/diameterPart

	# Averaging/Save
	def measure():
		global qsMean,vxPartPY,phiPartPY
		#Evaluate the average depth profile of streamwise, spanwise and wall-normal particle velocity, particle volume fraction (and drag force for coupling with RANS fluid resolution), and store it in hydroEngine variables vxPart, phiPart, vyPart, vzPart, averageDrag.
		hydroEngine.averageProfile()
		#Extract the calculated vector. They can be saved and plotted afterwards. 
		if oldVersion: 
			vxPartPY = np.array(hydroEngine.vxPart)
		else:
			vxPartPY = np.array(hydroEngine.vPart)[:,0]
		phiPartPY = np.array(hydroEngine.phiPart)

		#Evaluate the dimensionless sediment transport rate for information
		qsMean = sum(phiPartPY*vxPartPY)*dz/sqrt((densPart/densFluidPY - 1)*abs(gravityVector[2])*pow(diameterPart,3))

		#Evaluate the Shields number from the maximum of the Reynolds stresses evaluated in the fluid resolution
		shieldsNumber = max(hydroEngine.ReynoldStresses)/((densPart-densFluidPY)*diameterPart*abs(gravityVector[2]))	
		print('Shields number', shieldsNumber)

		if saveData==1:	#Save data for postprocessing
			global fileNumber
			nameFile = scriptPath + dataPath + str(fileNumber)+'.py'	# Name of the file that will be saved
			globalParam =  ['qsMean','phiPartPY','vxPartPY','vxFluidPY','zAxis']	# Variables to save
			Save(nameFile, globalParam)	#Save
			fileNumber+=1	#Increment the file number

		#Condition to stop the simulation after endTime seconds
		if O.time>=endTime:
			print('\n End of the simulation, simulated {0}s as required !\n '.format(endTime))
			O.pause()
	################
	##########################################


	#Save data details
	fileNumber = 0	# Counter for the file saved
	if saveData==1:	#If saveData option is activated, requires a folder data
		scriptPath = os.path.abspath(os.path.dirname(sys.argv[-1])) #Path where the script is stored
		dataPath = '/data/'
		if dataBaseRef:
			dataPath = '/dataRef/'
		if os.path.exists(scriptPath +dataPath)==False:
			os.mkdir(scriptPath +dataPath)
		else:
			print('\n!! Save data: overwrite the files contains in the folder data/ !!\n')
	#Function to save global variables in a python file which can be re-executed for post-processing
	def Save(filePathName, globalVariables):
		f = open(filePathName,'w')
		f.write('from numpy import *\n')
		for i in globalVariables:
			f.write(i + ' = '+repr(globals()[i]) + '\n')
		f.close()
