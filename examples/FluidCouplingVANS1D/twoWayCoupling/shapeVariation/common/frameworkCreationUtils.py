#########################################################################################################################################################################
# Author : Remi Monthiller, remi.monthiller@gmail.com
#
# Incline plane simulations
#
#########################################################################################################################################################################

import math
import random

from yade import pack
from yade.gridpfacet import *
"""Call frameworkCreation to create the framework in your simulation. You will need to define the parameters in the frameworkCreation.py file.

"""


class frCrea:  # Framework Creation

	####---------------####
	#### Immaterial
	####---------------####

	@staticmethod
	def createPeriodicCell(l=pM.l, w=pM.w, h=pM.h):
		""" Defines the domain as periodic

		"""
		# Defining the domain periodic
		O.periodic = True
		# Defining cell with good dimensions
		O.cell.setBox(l, w, h)

	@staticmethod
	def defineMaterials():
		""" Defines materials.

		"""
		## Estimated max particle pressure from the static load
		if pF.enable:
			p_max = -(pP.rho - pF.rho) * pP.phi_max * pM.n_l * pP.S * pM.g[2]
		else:
			p_max = -pP.rho * pP.phi_max * pM.n_l * pP.S * pM.g[2]

		## Evaluate the minimal normal stiffness to be in the rigid particle limit (cf Roux and Combe 2002)
		N_s = p_max * pP.S * 1.0e4
		## Young modulus of the particles from the stiffness wanted.
		E = N_s / pP.S
		## Poisson's ratio of the particles. Classical values, does not have much influence.
		nu = 0.5
		## Finaly difining material
		O.materials.append(ViscElMat(en=pP.c_r, et=0., young=E, poisson=nu, density=pP.rho, frictionAngle=pP.mu, label='mat'))

	####---------------####
	#### Macro
	####---------------####

	@staticmethod
	def createGround():
		""" Create the ground using a very big plane.

		"""
		cell_center = O.cell.size / 2.0
		ground = box(
		        center=(cell_center[0], cell_center[1], pM.z_ground),
		        extents=(1.0e5 * pM.l, 1.0e5 * pM.w, 0.0),
		        fixed=True,
		        wire=True,
		        color=(0., 0.5, 0.),
		        material='mat'
		)  # Flat bottom plane
		groundDisplay = box(
		        center=(cell_center[0], cell_center[1], pM.z_ground),
		        extents=(pM.l / 2.0, pM.w / 2.0, 0.0),
		        fixed=True,
		        wire=False,
		        color=(0., 0.5, 0.),
		        material='mat',
		        mask=0
		)  # Display of flat bottom plane
		O.bodies.append(ground)
		O.bodies.append(groundDisplay)
		return (len(O.bodies) - 2, len(O.bodies) - 1)

	@staticmethod
	def createRugosity():
		"""Creates a rugous ground based on spheres.

		"""
		cell_center = O.cell.size / 2.0

		randRangeZ = pM.d_rug / 2.0
		randRangeXY = 0.0
		d_eff = pM.d_rug

		# Create particles
		z = pM.z_ground
		x = -pM.l / 2.0 + cell_center[0] + d_eff / 2.0
		while x < pM.l / 2.0 + cell_center[0]:
			y = -pM.w / 2.0 + cell_center[1] + d_eff / 2.0
			while y < pM.w / 2.0 + cell_center[1]:
				O.bodies.append(
				        sphere(
				                center=(
				                        x + random.uniform(-1, 1) * randRangeXY, y + random.uniform(-1, 1) * randRangeXY,
				                        z + random.uniform(-1, 1) * randRangeZ
				                ),
				                radius=pM.d_rug / 2.0,
				                color=(0.0, 0.0, 0.0),
				                fixed=True,
				                material='mat',
				                wire=False
				        )
				)
				y += d_eff
			x += d_eff

	@staticmethod
	def createWalls():
		"""Creates walls.

		"""
		cell_center = O.cell.size / 2.0
		wall_1 = box(
		        center=(cell_center[0], cell_center[1] - pM.w / 2.0, pM.z_ground),
		        extents=(1.0e5 * pM.l, 0.0, 1.0e5 * pM.h),
		        fixed=True,
		        wire=True,
		        color=(0., 0.5, 0.),
		        material='mat'
		)
		display_1 = box(
		        center=(cell_center[0], cell_center[1] - pM.w / 2.0, pM.z_ground),
		        extents=(pM.l / 2.0, 0.0, pM.h / 2.0),
		        fixed=True,
		        wire=True,
		        color=(0., 0.5, 0.),
		        material='mat',
		        mask=0
		)
		wall_2 = box(
		        center=(cell_center[0], cell_center[1] + pM.w / 2.0, pM.z_ground),
		        extents=(1.0e5 * pM.l, 0.0, 1.0e5 * pM.h),
		        fixed=True,
		        wire=True,
		        color=(0., 0.5, 0.),
		        material='mat'
		)
		display_2 = box(
		        center=(cell_center[0], cell_center[1] + pM.w / 2.0, pM.z_ground),
		        extents=(pM.l / 2.0, 0.0, pM.h / 2.0),
		        fixed=True,
		        wire=True,
		        color=(0., 0.5, 0.),
		        material='mat',
		        mask=0
		)
		wall_1_id = O.bodies.append(wall_1)
		display_1_id = O.bodies.append(display_1)
		wall_2_id = O.bodies.append(wall_2)
		display_2_id = O.bodies.append(display_2)
		return (wall_1_id, display_1_id, wall_2_id, display_2_id)

	####---------------####
	#### Micro
	####---------------####

	@staticmethod
	def createBox(**kwargs):
		""" Creates a simple box

		Parameters :
		- center : center of the box : Vector3
		- extents : extents of the box : Vector3
		- orientation : orientation of the box : Quaternion
		- wire : display as wire : bool
		- color : color of the box : Vector3
		- mask : collision mask : int
		"""
		# Parameters
		center = kwargs.get('center', Vector3(0.0, 0.0, 0.0))
		extents = kwargs.get('extents', Vector3(0.0, 0.0, 0.0))
		orientation = kwargs.get('orientation', Quaternion(Vector3(0.0, 1.0, 0.0), 0.0))
		wire = kwargs.get('wire', False)
		color = kwargs.get('color', Vector3(1.0, 1.0, 1.0))
		mask = kwargs.get('mask', 1)

		# Creating ground at (0, 0, 0) with right length, width and orientation
		new_box = box(
		        center=center, extents=extents, orientation=orientation, fixed=True, wire=wire, color=color, material='mat', mask=mask
		)  # Flat bottom plane
		O.bodies.append(new_box)
		return len(O.bodies) - 1

	@staticmethod
	def createSphere(center, radius, color, wire, mask):
		""" Creates a simple sphere

		Parameters :
		- center : center of the sphere : Vector3
		- radius : radius of the sphere : float
		- color : color of the sphere : Vector3
		- wire : display as wire : bool
		- mask : collision mask : int
		"""
		# Creating ground at (0, 0, 0) with right length, width and orientation
		new_sphere = sphere(center=center, radius=radius, fixed=True, wire=wire, color=color, material='mat', mask=mask)  # Flat bottom plane
		O.bodies.append(new_sphere)
		return len(O.bodies) - 1

	####---------------####
	#### Display
	####---------------####

	@staticmethod
	def createDisplayFluidHeight():
		display = box(
		        center=(pM.l / 2.0, pM.w / 2.0, pF.z),
		        extents=(pM.l / 2.0, pM.w / 2.0, 0.0),
		        fixed=True,
		        wire=False,
		        color=(0., 0., 0.5),
		        material='mat',
		        mask=0
		)
		O.bodies.append(display)
		return len(O.bodies) - 1

	@staticmethod
	def createDisplayG():
		display = box(
		        center=(pM.l / 2.0, pM.w / 2.0, pF.z),
		        extents=(pM.w / 50.0, pM.w / 50.0, pM.h / 10.0),
		        orientation=Quaternion((0, -1, 0), pM.alpha),
		        fixed=True,
		        wire=False,
		        color=(0.5, 0., 0.),
		        material='mat',
		        mask=0
		)
		O.bodies.append(display)
		return len(O.bodies) - 1

	####---------------####
	#### Clump
	####---------------####

	@staticmethod
	def addClump(**kwargs):
		""" Add a 'fiber' clump.
			
			Parameters:
			- center : center of the clump : Vector3
			- ds : diameter of each spheres of the clump : float list
			- iter_vects : direction vectors : Vector3 list
		"""
		# Parameters
		center = kwargs.get('center', Vector3(0.0, 0.0, 0.0))
		ds = kwargs.get('ds', [1.0e-6])
		iter_vects = kwargs.get('iter_vects', [Vector3(1.0, 0, 0), Vector3(0, 1.0, 0)])
		color = kwargs.get('color', [Vector3(0.5, 0.0, 0.0)])

		# Sphere list
		ss = []
		d_max_tot = sum([max(row) for row in ds])
		center = center - d_max_tot / 2.0 * iter_vects[1]
		for row in ds:
			d_tot = sum(row)
			d_max = max(row)
			center += d_max / 2.0 * iter_vects[1]
			pos = center - d_tot / 2.0 * iter_vects[0]
			for d in row:
				pos += d / 2.0 * iter_vects[0]
				ss.append(sphere(center=pos, radius=d / 2.0, color=color, fixed=False, material='mat'))
				pos += d / 2.0 * iter_vects[0]
			center += d_max / 2.0 * iter_vects[1]
		# Adding clump to simulation
		result = O.bodies.appendClumped(ss)
		(id_clump, ids_clumped) = result

		# Random orientation
		O.bodies[id_clump].state.ori = Quaternion((0, 0, 1), random.uniform(0, 2 * math.pi)) * Quaternion(
		        (0, 1, 0), random.uniform(0, 2 * math.pi)
		) * Quaternion((1, 0, 0), random.uniform(0, 2 * math.pi))

		return result

	@staticmethod
	def createParticles(colors=[Vector3(0.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0)]):
		"""Creates a simple clump cloud.
		
		"""
		cell_center = O.cell.size / 2.0

		d_eff = sqrt(pow(pP.L, 2.0) + pow(pP.I, 2.0))
		n_i = 0

		iter_vects = [Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0)]

		# Create particles
		z = pM.z_ground + d_eff + d_eff / 2.0
		while n_i < pM.n:
			x = -pM.l / 2.0 + cell_center[0] + d_eff / 2.0
			while x < pM.l / 2.0 + cell_center[0] - d_eff / 2.0 and n_i < pM.n:
				y = -pM.w / 2.0 + cell_center[1] + d_eff / 2.0
				if y > pM.w / 2.0 + cell_center[1] - d_eff / 2.0:
					y = cell_center[1]
				y_init = y
				while (y == y_init) or (y < pM.w / 2.0 + cell_center[1] - d_eff / 2.0 and n_i < pM.n):
					col = random.uniform(0.0, 1.0)
					col = col * colors[0] + (1.0 - col) * colors[1]
					if pP.kind == "sphere":
						O.bodies.append(sphere(center=Vector3(x, y, z), radius=pP.dvs / 2.0, color=col, fixed=False, material='mat'))
					elif pP.kind == "clump":
						(id_clump, ids_clumped) = frCrea.addClump(center=Vector3(x, y, z), ds=pP.ds, iter_vects=iter_vects, color=col)
					n_i += 1
					#d_eff = frCrea.computeOutD(id_clump, ids_clumped)
					y += d_eff
				x += d_eff
			z += d_eff

	####---------------####
	#### Cylinder
	####---------------####

	@staticmethod
	def addCylinder(**kwargs):
		""" Creates a simple cylinder.

		Parameters :
		- center : center of the cylinder : Vector3
		- l : length of the cylinder : float
		- s : diameter of the cylinder : float
		"""
		# Parameters
		center = kwargs.get('center', Vector3(0.0, 0.0, 0.0))
		l = kwargs.get('l', [1.0e-6])
		s = kwargs.get('s', l / 4.0)

		# Sphere list
		col = (random.uniform(0.0, 0.5), random.uniform(0.0, 0.5), random.uniform(0.0, 0.5))
		nodesIds = []
		cylIds = []
		cylinder(
		        center - Vector3((l - s) / 2.0, 0.0, 0.0),
		        center + Vector3((l - s) / 2.0, 0.0, 0.0),
		        radius=s / 2.0,
		        nodesIds=nodesIds,
		        cylIds=cylIds,
		        fixed=False,
		        color=col,
		        intMaterial='cMat',
		        extMaterial='fMat'
		)

		# Random orientation
		# c_body = O.bodies[id_clump].state.ori = Quaternion((0, 0, 1), random.uniform(0, 2 * math.pi)) * Quaternion((0, 1, 0), random.uniform(0, 2 * math.pi)) * Quaternion((1, 0, 0), random.uniform(0, 2 * math.pi))
		return (cylIds[0], nodesIds)

	@staticmethod
	def createCylinderCloud():
		"""Creates a simple cylinder cloud.
		
		"""

		d_eff = pS.L
		r_eff = pS.L
		n_i = 0

		# Create particles
		z = pM.z_ground + d_eff + d_eff / 2.0
		while n_i < pM.n:
			x = -pM.l / 2.0
			while x < pM.l / 2.0 - d_eff / 2.0 and n_i < pM.n:
				y = -pM.w / 2.0
				while y < pM.w / 2.0 - d_eff / 2.0 and n_i < pM.n:
					(id_cyl, ids_nodes) = frCrea.addCylinder(
					        center=Vector3(x, y, z),
					        l=pS.L,
					        s=pS.S,
					)
					n_i += 1
					y += d_eff
				x += d_eff
			z += d_eff

	####---------------####
	#### Misc
	####---------------####

	@staticmethod
	def addTetra(**kwargs):
		# Parameters
		center = kwargs.get('center', Vector3(0.0, 0.0, 0.0))
		d_tetra = kwargs.get('d', 1.0e-6)
		l_tetra = kwargs.get('l', d_tetra)

		x = center[0]
		y = center[1]
		z = center[2]
		r_tetra = d_tetra / 2.0

		s_1 = sphere(center=(
		        x + l_tetra,
		        y,
		        z - l_tetra * sin(math.pi / 6.0),
		), radius=r_tetra, color=(1.0, 0.0, 0.0), fixed=False, material='mat')
		s_2 = sphere(
		        center=(x - l_tetra * sin(math.pi / 6.0), y - l_tetra * cos(math.pi / 6.0), z - l_tetra * sin(math.pi / 6.0)),
		        radius=r_tetra,
		        color=(1.0, 0.0, 0.0),
		        fixed=False,
		        material='mat'
		)
		s_3 = sphere(
		        center=(x - l_tetra * sin(math.pi / 6.0), y + l_tetra * cos(math.pi / 6.0), z - l_tetra * sin(math.pi / 6.0)),
		        radius=r_tetra,
		        color=(1.0, 0.0, 0.0),
		        fixed=False,
		        material='mat'
		)
		s_4 = sphere(center=(
		        x,
		        y,
		        z + l_tetra,
		), radius=r_tetra, color=(1.0, 0.0, 0.0), fixed=False, material='mat')
		return O.bodies.appendClumped([s_1, s_2, s_3, s_4])
