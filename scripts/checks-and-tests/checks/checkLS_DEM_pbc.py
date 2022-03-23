if ('LS_DEM' in features):
	# Handy function for comparing obtained quantities with expected ones
	######################################################################
	def equalNbr(x, xRef, tol=0.02):
		'''Whether x and xRef are equal numbers, up to a relative tolerance tol'''
		if abs(x - xRef) / abs(xRef) > tol:
			return False
		else:
			return True

	# Comparting the compression and decompression of a box of spheres to test if
	# periodic boundary conditions behave as expected.
	################################################################################

	# Settings and necessities
	import numpy as np
	from numpy import random
	ydim = 1.0
	xdim = 1.0
	zdim = 1.0
	nPartPerDim = 4.0
	spac = xdim / nPartPerDim / 2.0
	r = spac * 1.1

	#
	# Create a simulation with normal spheres as reference
	#
	O.periodic = True
	# Create a box with unequal dimension to create unequal stresses in each direction
	O.cell.hSize = Matrix3(xdim, 0.0, 0.0, 0.0, 0.95 * ydim, 0.0, 0.0, 0.0, 0.90 * zdim)

	# Alternate the step sizes such that the strain is differen and the particles correctly fit into the box
	for i in np.arange(spac + spac / 100, xdim - spac + 2 * spac / 100, 2 * spac):
		for j in np.arange(spac + spac / 100, ydim - spac + 2 * spac / 100, 1.9 * spac):
			for k in np.arange(spac + spac / 100, zdim - spac + 2 * spac / 100, 1.8 * spac):
				O.bodies.append(sphere(center=(i, j, k), radius=r))

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Sphere_Aabb()], label='collider'),
	        InteractionLoop(
	                [Ig2_Sphere_Sphere_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()], label='interactionLoop'
	        ),
	        PeriTriaxController(label='triaxPbc'),
	        NewtonIntegrator(damping=.2, label='newton')
	]

	timestep = 0.25 * PWaveTimeStep()
	O.dt = timestep
	O.run(1, wait=True)

	# Record reference values
	refStressXcomp = triaxPbc.stress[0]
	refStressYcomp = triaxPbc.stress[1]
	refStressZcomp = triaxPbc.stress[2]

	#
	# Create a scene for the simulation with level set spheres
	#
	lsSim = O.addScene()
	O.switchToScene(lsSim)

	O.periodic = True
	# Create a box with unequal dimension to create unequal stresses in each direction
	O.cell.hSize = Matrix3(xdim, 0.0, 0.0, 0.0, 0.95 * ydim, 0.0, 0.0, 0.0, 0.90 * zdim)

	# Alternate the step sizes such that the strain is differen and the particles correctly fit into the box
	for i in np.arange(spac + spac / 100, xdim - spac + 2 * spac / 100, 2 * spac):
		for j in np.arange(spac + spac / 100, ydim - spac + 2 * spac / 100, 1.9 * spac):
			for k in np.arange(spac + spac / 100, zdim - spac + 2 * spac / 100, 1.8 * spac):
				O.bodies.append(
				        levelSetBody(
				                "sphere",
				                center=(i, j, k),
				                radius=r,
				                nSurfNodes=10000,  #Need a high number of nodes to approach analytical sphere result
				                nodesPath=2,
				                dynamic=True,
				                spacing=0.01
				        )
				)

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_LevelSet_Aabb()], verletDist=0, label='colliderLS'),
	        InteractionLoop(
	                [Ig2_LevelSet_LevelSet_ScGeom()], [Ip2_FrictMat_FrictMat_FrictPhys()], [Law2_ScGeom_FrictPhys_CundallStrack()],
	                label='interactionLoopLS'
	        ),
	        PeriTriaxController(label='triaxLS'),
	        NewtonIntegrator(damping=.2, label='newtonLS')
	]

	O.dt = timestep  # The same time step as would be used for normal spheres
	O.run(1, wait=True)

	# Test the stresses and strains against the reference values, allow for an error of 2%
	if not equalNbr(triaxLS.stress[0], refStressXcomp, 0.02):
		raise YadeCheckError("Failed, periodic compact in LS-DEM gives a x stress of", triaxLS.stress[0], "vs", refStressXcomp, "expected")
	if not equalNbr(triaxLS.stress[1], refStressYcomp, 0.02):
		raise YadeCheckError("Failed, periodic compact in LS-DEM gives a y stress of", triaxLS.stress[1], "vs", refStressYcomp, "expected")
	if not equalNbr(triaxLS.stress[2], refStressZcomp, 0.02):
		raise YadeCheckError("Failed, periodic compact in LS-DEM gives a z stress of", triaxLS.stress[2], "vs", refStressZcomp, "expected")

	print('Periodic boundary conditions for LS-DEM as correct as expected')
else:
	print("Skip checkLS_DEM_pbc, LS-DEM feature not available")
