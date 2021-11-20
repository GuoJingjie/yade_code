import numpy as np
import math

m_A = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
m_phi = [0.61, 0.61, 0.62, 0.61, 0.60, 0.58, 0.56, 0.545, 0.54, 0.535, 0.54]
phi_max = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(m_A_1, m_phi, 3))

#m_A_2 = [1.0, 1.5, 2.0, 2.5, 3.0]
#m_mu = [tan(0.475), tan(0.545), tan(0.615), tan(0.615), tan(0.615)]
#mu_A = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(m_A_2, m_mu, 3))


### Numerical Parameters
class pN:
	### Creates new engines after loading. Useful to start from an previously computed time step.
	### For example, this parameter must be true if you want to run a dry simulation from the result of a wet simulation.
	enable_new_engines = False
	### Time of the simulation
	t_max = 400.0
	### Number of cells of the mesh
	n_z = 800
	### Verbose
	verbose = True


### Particle Parameters
class pP:
	### Shape and length parameters of the particles
	# Kind of simulated grains : "clump" or "cylinder". Cylinders do not work yet.
	kind = "clump"
	# Equivalent diameter of the particles.
	dvs = 6.0e-3
	# Shape factors
	A = 1.0
	B = 1.0
	## Number of spheres according to the parameters A and B
	## The total number of spheres is given by na*nb
	na = 3
	nb = 1
	# Characteristic lengths, the valueis are computed automatically in the following
	S = 1.0
	L = A * S
	I = B * S
	# Volume and frontale surface that will be computed differently depending on the kind
	vol = 0
	surf = 0
	### Clump parameters
	if kind == "clump":
		### Creating a row :
		# List of the diameters of the spheres in the clump
		dsa = []
		if L > S:
			# Computing diameter slope
			ca = ((2.0 * int(na / 2) + 1) * S - L) / (int(na / 2) * (int(na / 2) + 1))
			# Creating diameter row
			for i in range(int(na / 2)):
				dsa.append(S - ca * (int(na / 2) - i))
			dsa.append(S)
			for i in range(int(na / 2)):
				dsa.append(S - ca * (i + 1))
		else:
			dsa = [S]
		### Creating columns
		# List of the diameters of the spheres in the clump
		ds = []
		if I > S:
			# Computing diameter slope
			cb = ((2.0 * int(nb / 2) + 1) * S - I) / (int(nb / 2) * (int(nb / 2) + 1))
			# Creating diameters columns
			for i in range(int(nb / 2)):
				d = S - cb * (int(nb / 2) - i)
				ds.append([d / S * da for da in dsa])
			ds.append(dsa)
			for i in range(int(nb / 2)):
				d = S - cb * (i + 1)
				ds.append([d / S * da for da in dsa])
		else:
			ds = [dsa]
		# Computation of the volume and the surface of the particles
		vol = 0.0
		surf = 0.0
		for row in ds:
			for d in row:
				vol += math.pi * pow(d, 3) / 6.0
				surf += math.pi * pow(d, 2) / 4.0
		c = dvs / (3.0 / 2.0 * vol / surf)
		# Computation of the volume and the surface of the particles
		for i in range(len(ds)):
			for j in range(len(ds[i])):
				ds[i][j] *= c
		# Computation of the volume and the surface of the particles
		vol = 0.0
		surf = 0.0
		for row in ds:
			for d in row:
				vol += math.pi * pow(d, 3) / 6.0
				surf += math.pi * pow(d, 2) / 4.0
		S = max([max(row) for row in ds])
		L = A * S
		I = B * S
	# Verification
	dvs_calc = 3.0 / 2.0 * vol / surf
	### Cylinder parameters
	# TODO
	### Density of particles
	rho = 2.5e3
	### Coefficient of restitution
	c_r = 0.5
	### Maximum volume fraction
	phi_max = phi_max(A)
	### Friction angle
	mu = math.atan(0.5)
	### Initial particle velocity and volume fraction that are given to the HydroEngine
	v = [Vector3(0, 0, 0)] * (pN.n_z - 1)
	phi = [0] * (pN.n_z - 1)


### Macroscopic Parameters
class pM:
	### Sediment height
	hs = 12.0 * pP.dvs
	### Framework parameters
	# Slope
	alpha = 0.05
	# Cell dimensions
	l = 10 * pP.dvs * pP.A
	w = 10 * pP.dvs * pP.A
	h = 2000.0 * pP.dvs
	# Ground position
	z_ground = h / 2.0
	### Number of Particles
	n = pP.phi_max * l * w * hs / pP.vol
	# Number of particles "layers"
	n_l = n / (pP.phi_max * l * w * pP.S / pP.vol)
	### Gravity parameters
	g_scale = 9.81
	g = Vector3(g_scale * math.sin(alpha), 0, -g_scale * math.cos(alpha))
	### Ground Rugosity
	d_rug = pP.S
	### Shake Parameters
	shake_enable = True
	# Frequency of the sinusoide
	shake_f = 20.0
	# Time discretisation of the shake
	shake_dt = 0.05 / shake_f
	# Shake amplitude
	shake_a = pP.dvs / 2.0
	# Frequency of oscillatons :
	# 1/shake_wait_f time will be waited between each oscillations
	shake_wait_f = 19.0
	# Time until the shake stops
	shake_time = 0.6


### Param Save
class pSave:
	# Data will be saved all "yadeSavePeriod" simulation (virtual) time. Disable saving by setting it to 0.
	yadeSavePeriod = 2.0
	# Data will be saved as vtk (for Paraview for example) all "vtkRecorderIterPeriod" times. Disable saving by setting it to 0.
	vtkRecorderVirtPeriod = 0


### Param Fluid
class pF:
	enable = True
	solve = True
	# Can be "old" or "new"
	# "old" corresponds to the code before my internship (doesn't work for clumps)
	# "new" corresponds to the code after my internship (works for clumps)
	method = "new"
	# Will store data separately for each size of particles
	enable_poly_average = True
	# The fluid will start to be solved after solve_begin_time
	solve_begin_time = 0.8
	## Physics
	# Fluid density
	rho = 1e3
	# Fluid cinematic viscosity
	nu = 1e-6
	# Guess of the shields number of the simulation.
	init_shields = 0.5
	# Height of the fluid set automatically from the value of init_shields.
	h = 0.0
	if pM.alpha != 0 and enable and rho != 0:
		h = init_shields * (pP.rho / rho - 1) * pP.dvs / math.sin(pM.alpha)
	# Turbulent model type:
	# 1=Prandtl, 2=Li&Sawamoto, 3=Monthiller(turbulent structures pentration based), 4=Rousseau(pore-length based), 5=Ni&Capart(diameter based)
	# The parameter enable_poly_average must be set to True to use the turbulent model 3
	turbulence_model_type = 3
	# Model 2 parameter
	turb_phi_max = pP.phi_max
	# Model 5 parameters
	phi_bed = 0.08
	nb_average_over_time = 0
	## Numeric
	# Fluid resolution time discretisation
	dt = 1e-5
	# Fluid resolution period
	t = 1e-2
	## Fluid mesh
	# Richardson Zaki exponent for the hindrance function of the drag force applied to the particles ???
	expoDrag = 3.1
	# Other computed parameters
	z = pM.hs + h + pM.z_ground
	dz = (z - pM.z_ground) / float(pN.n_z - 1)
	# Attributes of the fluid
	vx = [0] * (pN.n_z + 1)
	# Display parameters
	display_enable = False
	display_n = 100
	display_mult = 0
	# Mostly useless parameters
	enable_wall_friction = False
	enable_fluctuations = False
	t_fluct = 1e-1


if pN.verbose:
	print("\n")

	print("INFO: A : " + str(pP.A))
	print("INFO: B : " + str(pP.B))
	print("INFO: phi_max : " + str(pP.phi_max))

	print("INFO: Particle volume : " + str(pP.vol))
	print("INFO: Particle surface : " + str(pP.surf))
	print("INFO: dvs:" + str(pP.dvs) + " dvs_calc:" + str(pP.dvs_calc))
	print("INFO: Diameter of the spheres: " + str(pP.ds))

	print("INFO: Number of particles : " + str(pM.n))
	print("INFO: Estimated number of particle layers : " + str(pM.n_l))

	if pF.enable:
		print("INFO: Estimated fluid height : " + str(pF.h))

	print("\n")
