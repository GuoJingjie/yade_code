#########################################################################################################################################################################
# Author : Remi Monthiller, remi.monthiller@etu.enseeiht.fr
# Adapted from the code of Raphael Maurin, raphael.maurin@imft.fr
# 30/10/2018
#
# Incline plane simulations
#
#########################################################################################################################################################################

import numpy as np


class sim:  # Simulation

	@staticmethod
	def simulation():
		"""Creates a simulation.
		
		"""
		sim.engineCreation()
		frameworkCreation()
		sim.init()
		# Save tmp so that the user can reload the simulation.
		O.saveTmp()

	@staticmethod
	def init():
		""" Initialise the simulation.

		"""
		### Initialisation of fluid
		if pF.enable:
			partsIds = []
			for i in range(len(O.bodies)):
				b = O.bodies[i]
				if not b.isClump and isinstance(b.shape, Sphere) and b.dynamic:
					partsIds.append(i)
			hydroEngine.ids = partsIds
			hydroEngine.ReynoldStresses = np.ones(pN.n_z) * 0.0
			hydroEngine.turbulentFluctuation()
			hydroEngine.initialization()

	####---------------####
	#### Engines Def
	####---------------####

	@staticmethod
	def cylinderContactDetectionCreation(engines):
		""" Add to the engines the engine necessary to cylinder contact detection.
		
		"""
		### Detect the potential contacts
		engines.append(InsertionSortCollider([Bo1_GridConnection_Aabb(), Bo1_Sphere_Aabb()], label='contactDetection', allowBiggerThanPeriod=True))

	@staticmethod
	def sphereContactDetectionCreation(engines):
		""" Add to the engines the engine necessary to sphere contact detection.
		
		"""
		### Detect the potential contacts
		engines.append(InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Box_Aabb()], label='contactDetection', allowBiggerThanPeriod=True))

	@staticmethod
	def cylinderInteractionsCreation(engines):
		""" Add to the engines the engine necessary to the cylinder interactions.
		
		"""
		### Calculate different interactions
		engines.append(
		        InteractionLoop(
		                [
		                        # Sphere interactions
		                        Ig2_Sphere_Sphere_ScGeom(),
		                        # Sphere-Cylinder interactions
		                        Ig2_Sphere_GridConnection_ScGridCoGeom(),
		                        # Cylinder interactions
		                        Ig2_GridNode_GridNode_GridNodeGeom6D(),
		                        # Cylinder-Cylinder interaction :
		                        Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
		                ],
		                [
		                        # Internal Cylinder Physics
		                        Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True, setCohesionOnNewContacts=False),
		                        # Physics for External Interactions (cylinder-cylinder)
		                        Ip2_FrictMat_FrictMat_FrictPhys()
		                ],
		                [
		                        # Contact law for sphere-sphere interactions
		                        Law2_ScGeom_FrictPhys_CundallStrack(),
		                        # Contact law for cylinder-sphere interactions
		                        Law2_ScGridCoGeom_FrictPhys_CundallStrack(),
		                        # Contact law for "internal" cylider forces :
		                        Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		                        # Contact law for cylinder-cylinder interaction :
		                        Law2_GridCoGridCoGeom_FrictPhys_CundallStrack()
		                ],
		                label='interactionLoop'
		        )
		)

	@staticmethod
	def sphereInteractionsCreation(engines):
		""" Add to the engines the engine necessary to the sphere interactions.
		
		"""
		engines.append(
		        InteractionLoop(
		                [Ig2_Sphere_Sphere_ScGeom(), Ig2_Box_Sphere_ScGeom()], [Ip2_ViscElMat_ViscElMat_ViscElPhys()], [Law2_ScGeom_ViscElPhys_Basic()],
		                label='interactionLoop'
		        )
		)

	####---------------####
	#### Hydro engine
	####---------------####

	@staticmethod
	def hydroEngineCreation(engines):
		""" Add to the engines the engine necessary to the application of hydrodynamic forces.
		
		"""
		# Creating Hydro Engine
		engines.append(
		        HydroForceEngine(
		                densFluid=pF.rho,
		                viscoDyn=pF.nu * pF.rho,
		                zRef=pM.z_ground,
		                gravity=pM.g,
		                deltaZ=pF.dz,
		                expoRZ=pF.expoDrag,
		                lift=False,
		                nCell=pN.n_z,
		                vCell=pM.l * pM.w * pF.dz,
		                phiPart=pP.phi,
		                vxFluid=pF.vx,
		                vPart=pP.v,
		                ids=[],
		                phiMax=pF.turb_phi_max,
		                ilm=pF.turbulence_model_type,
		                nbAverageT=pF.nb_average_over_time,
		                phiBed=pF.phi_bed,
		                enableMultiClassAverage=pF.enable_poly_average,
		                fluidWallFriction=pF.enable_wall_friction,
		                dead=True,
		                label='hydroEngine'
		        )
		)
		# Fluid resolution
		if pF.solve:
			engines.append(PyRunner(command='pyRuns.solveFluid()', virtPeriod=pF.t, label='fluidSolve'))
		# Turbulent fluctuations
		if pF.enable_fluctuations:
			engines.append(PyRunner(command='pyRuns.computeTurbulentFluctuations()', virtPeriod=pF.t_fluct, label='turbFluct'))
		# Display fluid velocity profile
		if pF.display_enable:
			engines.append(PyRunner(command='pyRuns.updateFluidDisplay()', virtPeriod=pF.t, label='fluidDisplay'))

	####---------------####
	#### Engines Creation
	####---------------####

	@staticmethod
	def engineCreation():
		"""Creates the engine.
		
		"""
		engines = []
		### Reset the forces
		engines.append(ForceResetter())
		### Detect the potential contacts
		if pP.kind == "cylinder":
			sim.cylinderContactDetectionCreation(engines)
		else:
			sim.sphereContactDetectionCreation(engines)
		### Calculate different interactions
		if pP.kind == "cylinder":
			sim.cylinderinteractionsCreation(engines)
		else:
			sim.sphereInteractionsCreation(engines)
		#### Apply hydrodynamic forces
		if pF.enable:
			sim.hydroEngineCreation(engines)
		### GlobalStiffnessTimeStepper, determine the time step for a stable integration
		engines.append(GlobalStiffnessTimeStepper(defaultDt=1e-4, viscEl=True, timestepSafetyCoefficient=0.7, label='GSTS'))
		### Integrate the equation and calculate the new position/velocities...
		engines.append(NewtonIntegrator(damping=0.0, gravity=pM.g, label='newtonIntegr'))
		### PyRunner Calls
		# End of simulation
		engines.append(PyRunner(command='pyRuns.exitWhenFinished()', virtPeriod=0.1, label='exit'))
		# Save data
		if pSave.yadeSavePeriod:
			engines.append(PyRunner(command='pyRuns.save()', virtPeriod=pSave.yadeSavePeriod, label='save'))
		# Shaker
		if pM.shake_enable:
			engines.append(PyRunner(command='pyRuns.shaker(pM.shake_axis)', virtPeriod=pM.shake_dt, label='shaker'))
		### Recorder
		if pSave.vtkRecorderVirtPeriod > 0:
			engines.append(VTKRecorder(virtPeriod=pSave.vtkRecorderVirtPeriod, recorders=['spheres', 'velocity', 'colors'], fileName='./vtk/sim-'))
		### Adding engines to Omega
		O.engines = engines

	@staticmethod
	def resetEngine():
		"""Creates the engine.
		
		"""
		#### Apply hydrodynamic forces
		if pF.enable:
			if 'hydroEngine' in globals():
				hydroEngine.dead = False
				# Creating Hydro Engine
				hydroEngine.densFluid = pF.rho
				hydroEngine.viscoDyn = pF.nu * pF.rho
				hydroEngine.zRef = pM.z_ground
				hydroEngine.gravity = pM.g
				hydroEngine.deltaZ = pF.dz
				hydroEngine.expoRZ = pF.expoDrag
				hydroEngine.lift = False
				hydroEngine.nCell = pN.n_z
				hydroEngine.vCell = pM.l * pM.w * pF.dz
				hydroEngine.phiPart = pP.phi
				hydroEngine.vxFluid = pF.vx
				hydroEngine.vPart = pP.v
				hydroEngine.ids = []
				hydroEngine.phiMax = pF.turb_phi_max
				hydroEngine.ilm = pF.turbulence_model_type
				ydroEngine.nbAverageT = pF.nb_average_over_time
				hydroEngine.phiBed = pF.phi_bed
				hydroEngine.enableMultiClassAverage = pF.enable_poly_average
				hydroEngine.fluidWallFriction = pF.enable_wall_friction
				hydroEngine.dead = True
				# Fluid resolution
				if pF.solve:
					if 'fluidSolve' in globals():
						fluidSolve.dead = False
						fluidSolve.virtPeriod = pF.t
					else:
						print("fluidSolve engine does not already exist and can't be added")
				else:
					if 'fluidSolve' in globals():
						fluidSolve.dead = True
				# Turbulent fluctuations
				if pF.enable_fluctuations:
					if 'turbFluct' in globals():
						turbFluct.dead = False
						turbFluct.virtPeriod = pF.t_fluct
					else:
						print("turbFluct engine does not already exist and can't be added")
				else:
					if 'turbFluct' in globals():
						turbFluct.dead = True
				# Display fluid velocity profile
				if pF.display_enable:
					if 'fluidDisplay' in globals():
						fluidDisplay.dead = False
						fluidDisplay.virtPeriod = pF.t
					else:
						print("fluidDisplay engine does not already exist and can't be added")
				else:
					if 'fluidDisplay' in globals():
						fluidDisplay.dead = True
			else:
				print("hydroEngine engine does not already exist and can't be added")
		else:
			if 'hydroEngine' in globals():
				hydroEngine.dead = True
		### Newton Integrator
		newtonIntegr.gravity = pM.g
		### PyRunner Calls
		# End of simulation
		exit.dead = False
		# Save data
		if pSave.yadeSavePeriod:
			if 'save' in globals():
				save.dead = False
				save.virtPeriod = pSave.yadeSavePeriod
			else:
				print("save engine does not already exist and can't be added")
		else:
			if 'save' in globals():
				save.dead = True
		# Shaker
		if pM.shake_enable:
			if 'shaker' in globals():
				shaker.dead = False
				shaker.virtPeriod = pM.shake_dt
			else:
				print("shaker engine does not already exist and can't be added")
		else:
			if 'shaker' in globals():
				shaker.dead = True


### Reading pyRunners
execfile("../common/simulationPyRunners.py")
