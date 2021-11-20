class pyRuns:

	@staticmethod
	def exitWhenFinished():
		if (O.time > pN.t_max):
			print("INFO: Simulation finished at time : " + str(O.time))
			O.pause()

	@staticmethod
	def save():
		if pF.enable:
			hydroEngine.averageProfile()
		O.save("data/" + str(O.time) + ".yade")

	count = 0

	@staticmethod
	def solveFluid():
		if O.time < pF.solve_begin_time:
			return
		elif pyRuns.count < 1 and O.time < pSave.yadeSavePeriod:
			print("INFO: Starting to apply fluid.\n")
			pyRuns.count = 1
			for e in O.engines:
				if e.label == "hydroEngine":
					e.dead = False
					break
			hydroEngine.averageProfile()
			hydroEngine.fluidResolution(1., pF.dt)
		# Computes average vx, vy, vz, phi, drag profiles
		hydroEngine.averageProfile()
		hydroEngine.fluidResolution(pF.t, pF.dt)

	@staticmethod
	def computeTurbulentFluctuations():
		if O.time < pF.solve_begin_time:
			return
		# Computes average vx, vy, vz, phi, drag profiles
		### Evaluate nBed, the position of the bed which is assumed to be located around the first maximum of concentration when considering decreasing z.
		phi = hydroEngine.phiPart
		nBed = pN.n_z - 2
		while nBed > 0 and not (phi[nBed] < phi[nBed - 1] and phi[nBed] > 0.5):
			# If there is a peak and its value is superior to 0.5,
			# we consider it to be the position of the bed
			nBed -= 1
		waterDepth = (pN.n_z - 1 - nBed) * pF.dz
		### Evaluate the bed elevation for the following
		bedElevation = pF.h - waterDepth
		### (Re)Define the bed elevation over which fluid turbulent fluctuations will be applied.
		hydroEngine.bedElevation = bedElevation
		### Impose a unique constant lifetime for the turbulent fluctuation, flucTimeScale
		# Fluid velocity scale in the water depth
		vMeanAboveBed = sum(hydroEngine.vxFluid[nBed:]) / (pN.n_z - nBed)
		# TODO : Very stange that it can be 0
		if vMeanAboveBed > 0:
			flucTimeScale = waterDepth / vMeanAboveBed  # time scale of the fluctuation w_d/v, eddy turn over time
			# New evaluation of the random fluid velocity fluctuation for each particle.
			hydroEngine.turbulentFluctuation()
			# Actualize when will be calculated the next fluctuations.
			turbFluct.virtPeriod = flucTimeScale

	shake_old_time = 0.0
	stop = False

	@staticmethod
	def shaker(axis=2):
		if pM.shake_wait_f > 0.0 and O.time > pyRuns.shake_old_time + 1.0 / pM.shake_wait_f:
			pyRuns.shake_old_time = O.time
			pyRuns.stop = False
		shaker_vel = 2.0 * pi * pM.shake_f * 0.5 * pM.shake_a * sin((O.time - pyRuns.shake_old_time) * 2.0 * pi * pM.shake_f)
		if (O.time - pyRuns.shake_old_time) > (1 / pM.shake_f):
			for b in O.bodies:
				if not b.dynamic:
					b.state.vel[axis] = 0.0
			if O.time > pM.shake_time:
				for e in O.engines:
					if e.label == "shaker":
						e.dead = True
						break
			if pM.shake_wait_f != 0.0:
				pyRuns.stop = True
		elif not pyRuns.stop:
			for b in O.bodies:
				if not b.dynamic:
					b.state.vel[axis] = shaker_vel

	fluidDisplayIds = []

	@staticmethod
	def updateFluidDisplay():
		if len(yade.qt.views()) > 0:
			# Managing Renderer
			renderer = yade.qt._GLViewer.Renderer()
			renderer.ghosts = False
			# Drawing profile
			dz = (pF.z - pM.z_ground) / pF.display_n
			nn = pN.n_z / pF.display_n
			vx = []
			vmax = 0.001
			for i in pyRuns.fluidDisplayIds:
				O.bodies.erase(i)
			mult = 0.0
			for i in range(pF.display_n):
				vx.append(0.0)
				for j in range(nn):
					vx[i] += hydroEngine.vxFluid[i * nn + j]
				vx[i] /= nn
				if vx[i] > vmax:
					vmax = vx[i]
			if pF.display_mult == 0.0:
				mult = (pM.l - pP.S) / vmax
			else:
				mult = pF.display_mult
			for i in range(pF.display_n):
				z = pM.z_ground + i * dz
				b = vx[i] / vmax
				v = b / 2
				r = 0.0
				pyRuns.fluidDisplayIds.append(
				        frCrea.createBox(
				                center=(max(pM.l / 100.0, vx[i] * mult / 2.0), pM.w / 2.0, z),
				                extents=(max(pM.l / 100.0, vx[i] * mult / 2.0), pM.w * 0.6, dz / 2.0),
				                color=(r, v, b),
				                wire=True,
				                mask=0
				        )
				)
		else:
			for i in pyRuns.fluidDisplayIds:
				O.bodies.erase(i)
			#fluidDisplay.dead = True
