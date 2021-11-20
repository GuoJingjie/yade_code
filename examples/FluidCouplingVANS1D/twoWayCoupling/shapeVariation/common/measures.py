#########################################################################################################################################################################
# Author : Remi Monthiller, remi.monthiller@etu.enseeiht.fr
# Adapted from the code of Raphael Maurin, raphael.maurin@imft.fr
# 30/10/2018
#
# Incline plane simulations
#
#########################################################################################################################################################################

import numpy as np

# --------------------------------- #
# Shields
# --------------------------------- #


def getShields():
	"""Returns the current shields number.
	
	"""
	if pF.enable:
		return max(hydroEngine.ReynoldStresses) / ((pP.rho - pF.rho) * pP.dvs * abs(pM.g[2]))
	else:
		return 0


# --------------------------------- #
# Orientation
# --------------------------------- #


def cart2sph(x, y, z):
	XsqPlusYsq = x**2 + y**2
	r = sqrt(XsqPlusYsq + z**2)
	theta = atan2(sqrt(XsqPlusYsq), z)
	phi = atan2(y, x)
	return r, theta, phi


def getOrientationHist(binsNb=3, kind="half"):
	rs = []
	thetas = []
	phis = []
	for body in O.bodies:
		if body.dynamic == True and body.isClump:
			u = (O.bodies[body.id - 1].state.pos - O.bodies[body.id - 3].state.pos).normalized()
			if kind != "full":
				if u.dot(Vector3(1.0, 0.0, 0.0)) < 0.0:
					u = -u
				if kind == "quart":
					if u[1] > 0.0:
						u[1] = -u[1]
			r, theta, phi = cart2sph(u[0], u[1], u[2])
			rs.append(r)
			thetas.append(theta)
			phis.append(phi)
	rs, thetas, phis = np.histogram2d(thetas, phis, bins=binsNb, normed=False)
	thetas = (thetas[:-1] + thetas[1:]) / 2.0
	phis = (phis[:-1] + phis[1:]) / 2.0
	X = []
	Y = []
	Z = []
	C = []
	n_sum = rs.sum()
	n_max = rs.max()
	for j in range(rs.shape[1]):
		for i in range(rs.shape[0]):
			n = rs[i, j] / n_sum
			n2 = rs[i, j] / n_max
			X.append(sin(thetas[i]) * cos(phis[j]) * n2)
			Y.append(sin(thetas[i]) * sin(phis[j]) * n2)
			Z.append(cos(thetas[i]) * n2)
			C.append(n)
	return [X, Y, Z, C]


def getMeanOrientation(kind="half"):
	rs = []
	thetas = []
	phis = []
	for body in O.bodies:
		if body.dynamic == True and body.isClump:
			u = (O.bodies[body.id - 1].state.pos - O.bodies[body.id - 3].state.pos).normalized()
			if kind != "full":
				if u.dot(Vector3(1.0, 0.0, 0.0)) < 0.0:
					u = -u
				if kind == "quart":
					if u[1] > 0.0:
						u[1] = -u[1]
			r, theta, phi = cart2sph(u[0], u[1], u[2])
			rs.append(r)
			thetas.append(theta)
			phis.append(phi)
	theta_mean = np.array(thetas).mean()
	theta_var = sqrt(np.array(thetas).var())
	phi_mean = np.array(phis).mean()
	phi_var = sqrt(np.array(phis).var())
	return [theta_mean, theta_var, phi_mean, phi_var]


def getVectorMeanOrientation(kind="half"):
	[theta_mean, theta_var, phi_mean, phi_var] = getMeanOrientation(kind)
	X = []
	Y = []
	Z = []
	C = []
	tmp = Vector3()
	tmp[0] = sin(theta_mean) * cos(phi_mean)
	tmp[1] = sin(theta_mean) * sin(phi_mean)
	tmp[2] = cos(theta_mean)
	tmp = tmp.normalized()
	X.append(tmp[0])
	Y.append(tmp[1])
	Z.append(tmp[2])
	C.append(1.0)
	tmp = Vector3()
	tmp[0] = sin(theta_mean + theta_var / 2.0) * cos(phi_mean)
	tmp[1] = sin(theta_mean + theta_var / 2.0) * sin(phi_mean)
	tmp[2] = cos(theta_mean + theta_var / 2.0)
	tmp = tmp.normalized()
	X.append(tmp[0])
	Y.append(tmp[1])
	Z.append(tmp[2])
	C.append(0.5)
	tmp[0] = sin(theta_mean - theta_var / 2.0) * cos(phi_mean)
	tmp[1] = sin(theta_mean - theta_var / 2.0) * sin(phi_mean)
	tmp[2] = cos(theta_mean - theta_var / 2.0)
	tmp = tmp.normalized()
	X.append(tmp[0])
	Y.append(tmp[1])
	Z.append(tmp[2])
	C.append(0.5)
	tmp[0] = sin(theta_mean) * cos(phi_mean + phi_var / 2.0)
	tmp[1] = sin(theta_mean) * sin(phi_mean + phi_var / 2.0)
	tmp[2] = cos(theta_mean)
	tmp = tmp.normalized()
	X.append(tmp[0])
	Y.append(tmp[1])
	Z.append(tmp[2])
	C.append(0.5)
	tmp[0] = sin(theta_mean) * cos(phi_mean - phi_var / 2.0)
	tmp[1] = sin(theta_mean) * sin(phi_mean - phi_var / 2.0)
	tmp[2] = cos(theta_mean)
	tmp = tmp.normalized()
	X.append(tmp[0])
	Y.append(tmp[1])
	Z.append(tmp[2])
	C.append(0.5)
	return [X, Y, Z, C]


def getOrientationProfiles(z_ground, dz, n_z, kind="half"):
	zs = [(i + 0.5) * dz for i in range(n_z)]
	thetas = [[] for i in range(n_z)]
	phis = [[] for i in range(n_z)]
	mean_theta = [0 for i in range(n_z)]
	var_theta = [0 for i in range(n_z)]
	mean_phi = [0 for i in range(n_z)]
	var_phi = [0 for i in range(n_z)]

	for body in O.bodies:
		if body.dynamic == True and body.isClump:
			i_pos = int((body.state.pos[2] - z_ground) / dz)
			if i_pos > pN.n_z and i_pos > 0:
				u = (O.bodies[body.id - 1].state.pos - O.bodies[body.id - 3].state.pos).normalized()
				if kind != "full":
					if u.dot(Vector3(1.0, 0.0, 0.0)) < 0.0:
						u = -u
					if kind == "quart":
						if u[1] > 0.0:
							u[1] = -u[1]
				r, theta, phi = cart2sph(u[0], u[1], u[2])
				thetas[i_pos].append(theta)
				phis[i_pos].append(phi)

	for i in range(n_z):
		mean_theta[i] = np.array(thetas[i]).mean()
		var_theta[i] = sqrt(np.array(thetas[i]).var())
		mean_phi[i] = np.array(phis[i]).mean()
		var_phi[i] = sqrt(np.array(phis[i]).var())

	return [zs, mean_theta, var_theta, mean_phi, var_phi]


# --------------------------------- #
# Get profiles utils.
# --------------------------------- #


def getDryProfile(prop):
	"""Should be removed at the end.
	
	"""
	partsIds = []
	for i in range(len(O.bodies)):
		b = O.bodies[i]
		if not b.isClump:
			partsIds.append(i)
	hydroEngineTmp = HydroForceEngine(
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
	        ids=partsIds,
	        label='hydroEngine',
	        dead=True
	)
	hydroEngineTmp.ReynoldStresses = np.ones(pN.n_z) * 0.0
	hydroEngineTmp.turbulentFluctuation()
	hydroEngineTmp.newInit()
	hydroEngineTmp.newAverageProfile()
	return eval("hydroEngineTmp." + prop)


def getVxPartProfile():
	"""Get the average velocity profile.
	
	"""
	if (pF.enable):
		vxPart = []
		if pF.method == "new":
			for i in range(0, len(hydroEngine.vPart)):
				vxPart.append(hydroEngine.vPart[i][0])
		elif pF.method == "old":
			for i in range(0, len(hydroEngine.vxPart)):
				vxPart.append(hydroEngine.vxPart[i])
		return vxPart
	else:
		partsIds = []
		for i in range(len(O.bodies)):
			b = O.bodies[i]
			if not b.isClump:
				partsIds.append(i)
		hydroEngineTmp = HydroForceEngine(
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
		        ids=partsIds,
		        label='hydroEngine',
		        dead=True
		)
		hydroEngineTmp.ReynoldStresses = np.ones(pN.n_z) * 0.0
		hydroEngineTmp.turbulentFluctuation()
		hydroEngineTmp.newInit()
		hydroEngineTmp.newAverageProfile()
		vxPart = []
		for v in hydroEngineTmp.vPart:
			vxPart.append(v[0])
		return vxPart


def getInertialProfile():
	"""Returns the current mean velocity value of the dynamic objects as a Vector3
	
	"""
	if (pF.enable):
		# Computation of \dot{\gamma}
		## Starting with vxPart
		vxPart = []
		for i in range(0, len(hydroEngine.vPart)):
			vxPart.append(hydroEngine.vPart[i][0])
		gamma_dot = []
		for i in range(0, len(vxPart)):
			gamma_dot.append((vxPart[i + 1][0] - vxPart[i][0]) / 2.0)
		# Computation of the Granular Pressure :
		# P_p = (\rho_p - \rho_f) * g * d_{vs} * \int{\phi dz}_z^\inf
		pressure = [0] * len(hydroEngine.phiPart)
		integ = 0
		for i in range(len(hydroEngine.phiPart) - 1, -1, -1):
			integ += hydroEngine.phiPart[i] * pF.dz
			pressure[i] = (pP.rho - pF.rho) * pM.g[2] * pP.dvs * integ
		# Finishing by computing the inertial number profile
		# I = \frac{\dot{\gamma} d}{\sqrt{\frac{P}{\rho_p}}}
		inertial = []
		for i in range(0, len(pressure)):
			inertial.append(gamma_dot[i] * dvs / sqrt(pressure / pP.rho))
		return [[i * pF.dz for i in range(pN.n_z)], inertial]
	else:
		partsIds = []
		for i in range(len(O.bodies)):
			b = O.bodies[i]
			if not b.isClump:
				partsIds.append(i)
		hydroEngineTmp = HydroForceEngine(
		        densFluid=pF.rho,
		        viscoDyn=pF.nu * pF.rho,
		        zRef=pM.z_ground,
		        gravity=pM.g,
		        deltaZ=pF.dz,
		        expoRZ=pF.expoDrag,
		        lift=False,
		        nCell=pM.n_z,
		        vCell=pM.l * pM.w * pF.dz,
		        radiusPart=pP.r,
		        phiPart=pP.phi,
		        vxFluid=pF.vx,
		        vPart=pP.v,
		        ids=partsIds,
		        label='hydroEngine',
		        dead=True
		)
		hydroEngineTmp.ReynoldStresses = np.ones(pM.n_z) * 0.0
		hydroEngineTmp.turbulentFluctuation()
		hydroEngineTmp.averageProfile()
		vxPart = []
		for v in hydroEngine.vPart:
			vxPart.append(v[0])
		return [[i * pF.dz for i in range(pM.n_z)], hydroEngineTmp.phiPart, vxPart, hydroEngineTmp.vxFluid[0:-1]]


# --------------------------------- #
# Get Old profile : Should be useless.
# --------------------------------- #


def getOldProfiles():
	dz = h / float(n_z)
	zs = [dz * i for i in range(n_z)]
	phi = [0 for i in range(n_z)]
	vpx = [0 for i in range(n_z)]
	for body in O.bodies:
		if body.dynamic == True and not body.isClump:
			z = body.state.pos[2]
			r = body.shape.radius
			z_min = z - r
			z_max = z + r
			n_min = int(z_min * n_z / h)
			n_max = int(z_max * n_z / h)
			vel = body.state.vel
			for i in range(n_min, n_max + 1):
				z_inf = max(zs[i], z_min) - z
				z_sup = min(zs[i + 1], z_max) - z
				vol = math.pi * pow(r, 2) * ((z_sup - z_inf) + (pow(z_inf, 3) - pow(z_sup, 3)) / (3 * pow(r, 2)))
				phi[i] += vol
				vpx[i] += vol * vel[0]
	for i in range(n_z):
		if (phi[i] > 0):
			vpx[i] /= phi[i]
			phi[i] /= dz * l * w
	return [zs, phi, vpx]
