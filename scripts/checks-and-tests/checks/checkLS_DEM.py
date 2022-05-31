if ('LS_DEM' in features):
	# Handy functions for comparing obtained quantities with expected ones
	######################################################################
	def equalVectors(n, nRef, tol=0.02):
		'''Whether n and nRef are equal or opposite vectors, up to a relative tolerance tol (that applies to the norm)'''
		if (n - nRef).norm() > tol * nRef.norm() and (n + nRef).norm() > tol * nRef.norm():
			return False
		else:
			return True

	def equalNbr(x, xRef, tol=0.02):
		'''Whether x and xRef are equal numbers, up to a relative tolerance tol'''
		if abs(x - xRef) / abs(xRef) > tol:
			return False
		else:
			return True

	# Starting with some volume and surface comparisons
	###################################################
	# sphere case:
	lsSph = levelSetBody('sphere', radius=1, spacing=0.05, nodesPath = 1)
	if not equalNbr(lsSph.shape.volume(), 4. / 3. * pi, 4.e-3):
		raise YadeCheckError("Failed because of an incorrect sphere volume in LS-DEM:", lsSph.shape.volume(), "vs", 4. / 3. * pi, "expected")
	lNorm = [nod.norm() for nod in lsSph.shape.surfNodes]
	if not equalNbr(max(lNorm), 1, 1.e-3) or not equalNbr(min(lNorm), 1, 1.e-3):
		raise YadeCheckError("Failed because of incorrect boundary nodes on a sphere in LS-DEM")
	if not equalNbr(lsSph.shape.getSurface(), 4*pi, 0.025): # tolerance could be 1.1e-3 with 2502 surface nodes
        	raise YadeCheckError("Failed because of incorrect getSurface() for a sphere in LS-DEM")

	# superellipsoid case:
	rx, ry, rz, epsE, epsN = 0.5, 1.2, 1., 0.1, 0.5
	lsSe = levelSetBody('superellipsoid', extents=(rx, ry, rz), epsilons=(epsE, epsN), spacing=0.05)

	# see A.H. Barr, in Graphics Gems III, D. Kirk (1995) for the following expressions of a superellipsoid volume. With log therein = natural logarithm
	def funG(x):  # see also http://people.math.sfu.ca/~cbm/aands/abramowitz_and_stegun.pdf p. 75
		if x <= 0:
			raise YadeCheckError("Failed because Gamma function does not apply to", x, "< 0")
		gam0, gam1, gam2, gam3 = 1. / 12, 1. / 30, 53. / 210, 195. / 371
		gam4, gam5 = 22999. / 22737, 29944523. / 19733142
		gam6 = 109535241009. / 48264275462
		contFrac = gam0 / (x + gam1 / (x + gam2 / (x + gam3 / (x + gam4 / (x + gam5 / (x + gam6 / x))))))
		return 0.5 * log(2 * pi) - x + (x - 0.5) * log(x) + contFrac

	def funGamma(x):  # exp(funG) is already the gamma function. Using nevertheless the factorial property for a more precise evaluation as below
		return exp(funG(x + 5)) / (x * (x + 1) * (x + 2) * (x + 3) * (x + 4))

	def beta(x, y):
		return funGamma(x) * funGamma(y) / funGamma(x + y)

	volExp = 2. / 3. * rx * ry**rz * epsE * epsN * beta(epsE / 2, epsE / 2) * beta(epsN, epsN / 2)
	if not equalNbr(lsSe.shape.volume(), volExp, 0.05):
		raise YadeCheckError("Failed because of an incorrect superellipsoid volume in LS-DEM:", lsSe.shape.volume(), "vs", volExp, "expected")
	print('LS-DEM volume and surface descriptions as correct as expected')

	# Now looking at the relative movements of 2 spheres and 2 LevelSet-shaped twins
	################################################################################

	rad = 1  # the smallest sphere
	rRatio = 1.8  # rBig / rSmall
	centrSmall, centrBig = (0, 0, 0), (0, 0, rad * (1 + rRatio))
	prec = 80  # grid fineness
	nSurfNodes = 2502
	# the 2 true spheres, along z-axis
	O.bodies.append(sphere(centrSmall, rad, dynamic=False))
	O.bodies.append(sphere(centrBig, rRatio * rad, dynamic=False))
	movSph = O.bodies[1]

	# their LevelSet counterparts:
	O.bodies.append(levelSetBody('sphere', centrSmall, rad, spacing=2 * rad / prec, nSurfNodes=nSurfNodes, nodesPath=1, dynamic=False))
	O.bodies.append(levelSetBody('sphere', centrBig, rRatio * rad, spacing=2 * rRatio * rad / prec, nSurfNodes=nSurfNodes, nodesPath=1, dynamic=False))
	movLS = O.bodies[3]

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_LevelSet_Aabb()]),
	        InteractionLoop(
	                [Ig2_Sphere_Sphere_ScGeom(avoidGranularRatcheting=False),
	                 Ig2_LevelSet_LevelSet_ScGeom()],
	                [Ip2_FrictMat_FrictMat_FrictPhys(kn=MatchMaker(algo='val', val=1.e7), ks=MatchMaker(algo='val', val=1.e7))],
	                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)]
	        ),
	        NewtonIntegrator()
	]
	O.dt = 5.e-4

	# 1. Normal relative displacement
	movSph.state.vel = movLS.state.vel = (0, 0, -1.)
	O.run(100, True)

	if not O.interactions[2, 3]:
		raise YadeCheckError("Failed, we do not have LS interaction after first stage of relative normal displacement")
	lsCont = O.interactions[2, 3]
	sphCont = O.interactions[0, 1]

	if not equalNbr(lsCont.geom.penetrationDepth, sphCont.geom.penetrationDepth, 1.e-12):  # 2.2e-13 is a feasible goal on that ideal case
		raise YadeCheckError(
		        "Failed, normal overlap is too wrong in LS-DEM after first stage:", lsCont.geom.penetrationDepth, "vs", sphCont.geom.penetrationDepth,
		        "in DEM"
		)

	if not equalVectors(lsCont.geom.normal, sphCont.geom.normal):
		raise YadeCheckError("Failed, the two normals are too different after first stage:", sphCont.geom.normal, "vs", lsCont.geom.normal)
	movSph.state.vel = movLS.state.vel = Vector3.Zero

	#2. Circular relative displacement (~ pure shear)
	arc = rad * (1 + rRatio)
	dAlpha = 1.2 * pi / 4.
	nIt = 1500
	normV = arc * dAlpha / (nIt * O.dt)
	lsShearDisp = sphShearDisp = Vector3.Zero
	for i in range(nIt):
		movSph.state.vel = movLS.state.vel = normV * Vector3(cos(dAlpha / nIt * i), 0, -sin(dAlpha / nIt * i))
		O.step()
		lsShearDisp += lsCont.geom.shearInc
		sphShearDisp += sphCont.geom.shearInc

	if not equalNbr(
	        lsCont.geom.penetrationDepth, sphCont.geom.penetrationDepth, 0.03
	):  # 0.0284 error is expected here, would be eg 0.007 with 6402 nodes and grid precision 80
		raise YadeCheckError(
		        "Failed, normal overlaps are too different after 2nd stage:", lsCont.geom.penetrationDepth, "vs", sphCont.geom.penetrationDepth
		)
	if not equalVectors(lsShearDisp, sphShearDisp):
		raise YadeCheckError("Failed, the two shear displacements are too different after second stage:", sphShearDisp, "vs", lsShearDisp)
	if not equalVectors(
	        lsCont.geom.normal, sphCont.geom.normal, 0.03
	):  # allowing here 3 % of error. 6402 nodes and grid precision 80 would allow to go under 2 %
		raise YadeCheckError("Failed, the two normals are too different after second stage:", sphCont.geom.normal, "vs", lsCont.geom.normal)
	print('LS-DEM contact description as correct as expected')

	# Consideration of Fast Marching Method to finish
	#################################################
	grid = RegularGrid(-1.1, 1.1, 23)  # a cubic grid from -1.1 to 1.1 with 23 gp ie a 0.1 step
	fmm = FastMarchingMethod(
	        phiIni=distIniSE(radii=[1, 1, 1], epsilons=[1, 1], grid=grid), grid=grid
	)  # checking fast marching method when applied to the distance to the unit sphere
	phiField = fmm.phi()
	error = 0
	for i in range(23):
		for j in range(23):
			for k in range(23):
				phi_ijk = phiField[i][j][k]
				error += abs(grid.gridPoint(i, j, k).norm() - 1 - phi_ijk)
	error /= grid.nGP.prod()  # average (and dimensionless with respect to unit radius)
	errorExpected = 0.009393853398395624  # e.g. on Ubuntu 20.04.3 and jduriez axp17* while would be 0.009419100794945902 on Debian Bullseye, see https://gitlab.com/yade-dev/trunk/-/jobs/1832563583
	if not equalNbr(error, errorExpected, 5.e-3):
		raise YadeCheckError("Failed, Fast Marching Method gives an error of", error, "vs", errorExpected, "expected")
	print('Fast Marching Method as correct as expected')
else:
	print("Skip checkLSdem, LS-DEM feature not available")
