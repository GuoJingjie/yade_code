"""
Import geometry from various formats ('import' is python keyword, hence the name 'ymport').
"""

from builtins import range
from builtins import object
from yade.wrapper import *
from yade import utils

from yade.minieigenHP import *


def textExt(fileName, format='x_y_z_r', shift=Vector3.Zero, scale=1.0, attrs=[], **kw):
	"""Load sphere coordinates from file in a format selected by the ``format`` argument, returns a list of corresponding bodies; that may be inserted to the simulation with O.bodies.append().

	:param str filename: file name
	:param str format: selected input format. Supported ``'x_y_z_r'``(default), ``'x_y_z_r_matId'``, ``'x_y_z_r_attrs'``
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param list attrs: attrs read from file if export.textExt(format='x_y_z_r_attrs') were used ('passed by reference' style)
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""
	infile = open(fileName, "r")
	lines = infile.readlines()
	infile.close()
	ret = []
	for line in lines:
		data = line.split()
		if (data[0] == "#format"):
			format = data[1]
			continue
		elif (data[0][0] == "#"):
			continue

		if (format == 'x_y_z_r'):
			pos = Vector3(float(data[0]), float(data[1]), float(data[2]))
			ret.append(utils.sphere(shift + scale * pos, scale * float(data[3]), **kw))
		elif (format == 'x_y_z_r_matId'):
			pos = Vector3(float(data[0]), float(data[1]), float(data[2]))
			ret.append(utils.sphere(shift + scale * pos, scale * float(data[3]), material=int(data[4]), **kw))

		elif (format == 'id_x_y_z_r_matId'):
			pos = Vector3(float(data[1]), float(data[2]), float(data[3]))
			ret.append(utils.sphere(shift + scale * pos, scale * float(data[4]), material=int(data[5]), **kw))

		elif (format == 'x_y_z_r_attrs'):
			pos = Vector3(float(data[0]), float(data[1]), float(data[2]))
			s = utils.sphere(shift + scale * pos, scale * float(data[3]), **kw)
			ret.append(s)
			attrs.append(data[4:])

		else:
			raise RuntimeError("Please, specify a correct format output!")
	return ret


def textFacets(fileName, format='x1_y1_z1_x2_y2_z2_x3_y3_z3', shift=Vector3.Zero, scale=1.0, attrs=[], **kw):
	"""Load facet coordinates from file in a format selected by the ``format`` argument, returns a list of corresponding bodies; that may be inserted to the simulation with O.bodies.append().
	
	:param str filename: file name
	:param str format: selected input format. Supported ``'x1_y1_z1_x2_y2_z2_x3_y3_z3'``(default), ``'x1_y1_z1_x2_y2_z2_x3_y3_z3_matId'``, ``'id_x1_y1_z1_x2_y2_z2_x3_y3_z3_matId'`` or ``'x1_y1_z1_x2_y2_z2_x3_y3_z3_attrs'``
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param list attrs: attrs read from file ('passed by reference' style)
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.facet`
	:returns: list of facets.

	Lines starting with # are skipped
	"""
	infile = open(fileName, "r")
	lines = infile.readlines()
	infile.close()
	ret = []
	for line in lines:
		data = line.split()
		if (data[0] == "#format"):
			format = data[1]
			continue
		elif (data[0][0] == "#"):
			continue

		if (format == 'x1_y1_z1_x2_y2_z2_x3_y3_z3'):
			V1 = Vector3(float(data[0]), float(data[1]), float(data[2]))
			V2 = Vector3(float(data[3]), float(data[4]), float(data[5]))
			V3 = Vector3(float(data[6]), float(data[7]), float(data[8]))
			ret.append(utils.facet((shift + V1 * scale, shift + V2 * scale, shift + V3 * scale), **kw))

		elif (format == 'x1_y1_z1_x2_y2_z2_x3_y3_z3_matId'):
			V1 = Vector3(float(data[0]), float(data[1]), float(data[2]))
			V2 = Vector3(float(data[3]), float(data[4]), float(data[5]))
			V3 = Vector3(float(data[6]), float(data[7]), float(data[8]))
			ret.append(utils.facet((shift + V1 * scale, shift + V2 * scale, shift + V3 * scale), material=int(data[9]), **kw))

		elif (format == 'id_x1_y1_z1_x2_y2_z2_x3_y3_z3_matId'):
			V1 = Vector3(float(data[1]), float(data[2]), float(data[3]))
			V2 = Vector3(float(data[4]), float(data[5]), float(data[6]))
			V3 = Vector3(float(data[7]), float(data[8]), float(data[9]))
			ret.append(utils.facet((shift + V1 * scale, shift + V2 * scale, shift + V3 * scale), material=int(data[10]), **kw))

		elif (format == 'x1_y1_z1_x2_y2_z2_x3_y3_z3_attrs'):
			V1 = Vector3(float(data[0]), float(data[1]), float(data[2]))
			V2 = Vector3(float(data[3]), float(data[4]), float(data[5]))
			V3 = Vector3(float(data[6]), float(data[7]), float(data[8]))
			ret.append(utils.facet((shift + V1 * scale, shift + V2 * scale, shift + V3 * scale), **kw))
			attrs.append(data[4:])

		else:
			raise RuntimeError("Please, specify a correct format output!")
	return ret


def textClumps(fileName, shift=Vector3.Zero, discretization=0, orientation=Quaternion((0, 1, 0), 0.0), scale=1.0, **kw):
	"""Load clumps-members from file in a format selected by the ``format`` argument, insert them to the simulation.

	:param str filename: file name
	:param str format: selected input format. Supported ``'x_y_z_r'``(default), ``'x_y_z_r_clumpId'``
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""
	infile = open(fileName, "r")
	lines = infile.readlines()
	infile.close()
	ret = []

	curClump = []
	newClumpId = -1

	for line in lines:
		data = line.split()
		if (data[0][0] == "#"):
			continue
		pos = orientation * Vector3(float(data[0]), float(data[1]), float(data[2]))

		if (newClumpId < 0 or newClumpId == int(data[4])):
			idD = curClump.append(utils.sphere(shift + scale * pos, scale * float(data[3]), **kw))
			newClumpId = int(data[4])
		else:
			newClumpId = int(data[4])
			ret.append(O.bodies.appendClumped(curClump, discretization=discretization))
			curClump = []
			idD = curClump.append(utils.sphere(shift + scale * pos, scale * float(data[3]), **kw))

	if (len(curClump) != 0):
		ret.append(O.bodies.appendClumped(curClump, discretization=discretization))

	# Set the mask to a clump the same as the first member of it
	for i in range(len(ret)):
		O.bodies[ret[i][0]].mask = O.bodies[ret[i][1][0]].mask
	return ret


def text(fileName, shift=Vector3.Zero, scale=1.0, **kw):
	"""Load sphere coordinates from file, returns a list of corresponding bodies; that may be inserted to the simulation with O.bodies.append().

	:param string filename: file which has 4 colums [x, y, z, radius].
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param \*\*kw: (unused keyword arguments)	is passed to :yref:`yade.utils.sphere`
	:returns: list of spheres.

	Lines starting with # are skipped
	"""

	return textExt(fileName=fileName, format='x_y_z_r', shift=shift, scale=scale, **kw)


def stl(file, dynamic=None, fixed=True, wire=True, color=None, highlight=False, noBound=False, material=-1, scale=1.0, shift=Vector3.Zero):
	""" Import a .stl geometry in the form of a set of :yref:`Facet`-shaped bodies.
	
	:param string file: the .stl file serving as geometry input
	:param bool dynamic: controls :yref:`Body.dynamic`
	:param bool fixed: controls :yref:`Body.dynamic` (with fixed = True imposing :yref:`Body.dynamic` = False) if *dynamic* attribute is not given
	:param bool wire: rendering option, passed to :yref:`Facet.wire`
	:param color: rendering option, passed to :yref:`Facet.color`
	:param bool highlight: rendering option, passed to :yref:`Facet.highlight`
	:param bool noBound: sets :yref:`Body.bounded` to False if True, preventing collision detection (and vice-versa)
	:param material: defines :yref:`material<Body.material>` properties, see :ref:`DefiningMaterials` for usage
	:param float scale: scaling factor to e.g. dilate the geometry if > 1
	:param Vector3 shift: for translating the geometry
	:returns: a corresponding list of :yref:`Facet`-shaped bodies"""
	imp = STLImporter()
	facets = imp.ymport(file)
	for b in facets:
		b.shape.setVertices(b.shape.vertices[0] * scale, b.shape.vertices[1] * scale, b.shape.vertices[2] * scale)
		b.shape.color = color if color else utils.randomColor()
		b.shape.wire = wire
		b.shape.highlight = highlight
		pos = b.state.pos * scale + shift
		utils._commonBodySetup(b, 0, Vector3(0, 0, 0), material=material, pos=pos, noBound=noBound, dynamic=dynamic, fixed=fixed)
		b.aspherical = False
	return facets


def gts(meshfile, shift=Vector3.Zero, scale=1.0, **kw):
	""" Read given meshfile in gts format.

	:Parameters:
		`meshfile`: string
			name of the input file.
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.facet`
	:Returns: list of facets.
	"""
	import gts, yade.pack
	surf = gts.read(open(meshfile))
	surf.scale(scale, scale, scale)
	surf.translate(shift[0], shift[1], shift[2])
	yade.pack.gtsSurface2Facets(surf, **kw)


def gmsh(meshfile="file.mesh", shift=Vector3.Zero, scale=1.0, orientation=Quaternion((0, 1, 0), 0.0), **kw):
	""" Imports geometry from .mesh file and creates facets.

	:Parameters:
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`orientation`: quaternion
			orientation of the imported mesh
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.facet`
	:Returns: list of facets forming the specimen.
	
	mesh files can easily be created with `GMSH <http://www.geuz.org/gmsh/>`_.
	Example added to :ysrc:`examples/packs/packs.py`
	
	Additional examples of mesh-files can be downloaded from 
	http://www-roc.inria.fr/gamma/download/download.php
	"""
	infile = open(meshfile, "r")
	lines = infile.readlines()
	infile.close()

	nodelistVector3 = []
	elementlistVector3 = []  # for deformable elements
	findVerticesString = 0

	while (lines[findVerticesString].split()[0] != 'Vertices'):  #Find the string with the number of Vertices
		findVerticesString += 1
	findVerticesString += 1
	numNodes = int(lines[findVerticesString].split()[0])

	for i in range(numNodes):
		nodelistVector3.append(Vector3(0.0, 0.0, 0.0))
	id = 0

	for line in lines[findVerticesString + 1:numNodes + findVerticesString + 1]:
		data = line.split()
		nodelistVector3[id] = orientation * Vector3(float(data[0]) * scale, float(data[1]) * scale, float(data[2]) * scale) + shift
		id += 1

	findTriangleString = findVerticesString + numNodes
	while (lines[findTriangleString].split()[0] != 'Triangles'):  #Find the string with the number of Triangles
		findTriangleString += 1
	findTriangleString += 1
	numTriangles = int(lines[findTriangleString].split()[0])

	triList = []
	for i in range(numTriangles):
		triList.append([0, 0, 0, 0])

	tid = 0
	for line in lines[findTriangleString + 1:findTriangleString + numTriangles + 1]:
		data = line.split()
		id1 = int(data[0]) - 1
		id2 = int(data[1]) - 1
		id3 = int(data[2]) - 1
		triList[tid][0] = tid
		triList[tid][1] = id1
		triList[tid][2] = id2
		triList[tid][3] = id3
		tid += 1
		ret = []
	for i in triList:
		a = nodelistVector3[i[1]]
		b = nodelistVector3[i[2]]
		c = nodelistVector3[i[3]]
		ret.append(utils.facet((nodelistVector3[i[1]], nodelistVector3[i[2]], nodelistVector3[i[3]]), **kw))
	return ret


def gengeoFile(fileName="file.geo", shift=Vector3.Zero, scale=1.0, orientation=Quaternion((0, 1, 0), 0.0), **kw):
	""" Imports geometry from LSMGenGeo .geo file and creates spheres. 
	Since 2012 the package is available in Debian/Ubuntu and known as python-demgengeo
	http://packages.qa.debian.org/p/python-demgengeo.html
	
	:Parameters:
		`filename`: string
			file which has 4 colums [x, y, z, radius].
		`shift`: Vector3
			Vector3(X,Y,Z) parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`orientation`: quaternion
			orientation of the imported geometry
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	:Returns: list of spheres.
	
	LSMGenGeo library allows one to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.
	
	Example: :ysrc:`examples/packs/packs.py`, usage of LSMGenGeo library in :ysrc:`examples/test/genCylLSM.py`.
	
	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	from yade.utils import sphere

	infile = open(fileName, "r")
	lines = infile.readlines()
	infile.close()

	numSpheres = int(lines[6].split()[0])
	ret = []
	for line in lines[7:numSpheres + 7]:
		data = line.split()
		pos = orientation * Vector3(float(data[0]), float(data[1]), float(data[2]))
		ret.append(utils.sphere(shift + scale * pos, scale * float(data[3]), **kw))
	return ret


def gengeo(mntable, shift=Vector3.Zero, scale=1.0, **kw):
	""" Imports geometry from LSMGenGeo library and creates spheres.
	Since 2012 the package is available in Debian/Ubuntu and known as python-demgengeo
	http://packages.qa.debian.org/p/python-demgengeo.html

	:Parameters:
		`mntable`: mntable
			object, which creates by LSMGenGeo library, see example
		`shift`: [float,float,float]
			[X,Y,Z] parameter moves the specimen.
		`scale`: float
			factor scales the given data.
		`**kw`: (unused keyword arguments)
				is passed to :yref:`yade.utils.sphere`
	
	LSMGenGeo library allows one to create pack of spheres
	with given [Rmin:Rmax] with null stress inside the specimen.
	Can be useful for Mining Rock simulation.
	
	Example: :ysrc:`examples/packs/packs.py`, usage of LSMGenGeo library in :ysrc:`examples/test/genCylLSM.py`.
	
	* https://answers.launchpad.net/esys-particle/+faq/877
	* http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
	* https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/"""
	try:
		from GenGeo import MNTable3D, Sphere
	except ImportError:
		from gengeo import MNTable3D, Sphere
	ret = []
	sphereList = mntable.getSphereListFromGroup(0)
	for i in range(0, len(sphereList)):
		r = sphereList[i].Radius()
		c = sphereList[i].Centre()
		ret.append(
		        utils.sphere(
		                [shift[0] + scale * float(c.X()), shift[1] + scale * float(c.Y()), shift[2] + scale * float(c.Z())], scale * float(r), **kw
		        )
		)
	return ret


def unv(fileName, shift=(0, 0, 0), scale=1.0, returnConnectivityTable=False, **kw):
	""" Import geometry from unv file, return list of created facets.

		:param string fileName: name of unv file
		:param (float,float,float)|Vector3 shift: (X,Y,Z) parameter moves the specimen.
		:param float scale: factor scales the given data.
		:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.facet`
		:param bool returnConnectivityTable: if True, apart from facets returns also nodes (list of (x,y,z) nodes coordinates) and elements (list of (id1,id2,id3) element nodes ids). If False (default), returns only facets
	
	unv files are mainly used for FEM analyses (are used by `OOFEM <http://www.oofem.org/>`_ and `Abaqus <http://www.simulia.com/products/abaqus_fea.html>`_), but triangular elements can be imported as facets.
	These files cen be created e.g. with open-source free software `Salome <http://salome-platform.org>`_.
	
	Example: :ysrc:`examples/test/unv-read/unvRead.py`."""

	class UNVReader(object):
		# class used in ymport.unv function
		# reads and evaluate given unv file and extracts all triangles
		# can be extended to read tetrahedrons as well
		def __init__(self, fileName, shift=(0, 0, 0), scale=1.0, returnConnectivityTable=False, **kw):
			self.shift = shift
			self.scale = scale
			self.unvFile = open(fileName, 'r')
			self.flag = 0
			self.line = self.unvFile.readline()
			self.lineSplit = self.line.split()
			self.nodes = []
			self.elements = []
			self.read(**kw)

		def readLine(self):
			self.line = self.unvFile.readline()
			self.lineSplit = self.line.split()

		def read(self, **kw):
			while self.line:
				self.evalLine()
				self.line = self.unvFile.readline()
			self.unvFile.close()
			self.createFacets(**kw)

		def evalLine(self):
			self.lineSplit = self.line.split()
			if len(self.lineSplit) <= 1:  # eval special unv format
				if self.lineSplit[0] == '-1':
					pass
				elif self.lineSplit[0] == '2411':
					self.flag = 1
					# nodes
				elif self.lineSplit[0] == '2412':
					self.flag = 2
					# edges (lines)
				else:
					self.flag = 4
					# volume elements or other, not interesting for us (at least yet)
			elif self.flag == 1:
				self.evalNodes()
			elif self.flag == 2:
				self.evalEdge()
			elif self.flag == 3:
				self.evalFacet()
			#elif self.flag == 4: self.evalGroup()
		def evalNodes(self):
			self.readLine()
			self.nodes.append(
			        (
			                self.shift[0] + self.scale * float(self.lineSplit[0]), self.shift[1] + self.scale * float(self.lineSplit[1]),
			                self.shift[2] + self.scale * float(self.lineSplit[2])
			        )
			)

		def evalEdge(self):
			if self.lineSplit[1] == '41':
				self.flag = 3
				self.evalFacet()
			else:
				self.readLine()
				self.readLine()

		def evalFacet(self):
			if self.lineSplit[1] == '41':  # triangle
				self.readLine()
				self.elements.append((int(self.lineSplit[0]) - 1, int(self.lineSplit[1]) - 1, int(self.lineSplit[2]) - 1))
			else:  # is not triangle
				self.readLine()
				self.flag = 4
				# can be added function to handle tetrahedrons
		def createFacets(self, **kw):
			self.facets = [utils.facet(tuple(self.nodes[i] for i in e), **kw) for e in self.elements]

	#
	unvReader = UNVReader(fileName, shift, scale, returnConnectivityTable, **kw)
	if returnConnectivityTable:
		return unvReader.facets, unvReader.nodes, unvReader.elements
	return unvReader.facets


def iges(fileName, shift=(0, 0, 0), scale=1.0, returnConnectivityTable=False, **kw):
	""" Import triangular mesh from .igs file, return list of created facets.

		:param string fileName: name of iges file
		:param (float,float,float)|Vector3 shift: (X,Y,Z) parameter moves the specimen.
		:param float scale: factor scales the given data.
		:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.facet`
		:param bool returnConnectivityTable: if True, apart from facets returns also nodes (list of (x,y,z) nodes coordinates) and elements (list of (id1,id2,id3) element nodes ids). If False (default), returns only facets
	"""
	nodes, elems = [], []
	f = open(fileName)
	for line in f:
		if line.startswith('134,'):  # read nodes coordinates
			ls = line.split(',')
			v = Vector3(float(ls[1]) * scale + shift[0], float(ls[2]) * scale + shift[1], float(ls[3]) * scale + shift[2])
			nodes.append(v)
		if line.startswith('136,'):  # read elements
			ls = line.split(',')
			i1, i2, i3 = int(ls[3]) / 2, int(ls[4]) / 2, int(ls[5]) / 2  # the numbering of nodes is 1,3,5,7,..., hence this int(ls[*])/2
			elems.append((i1, i2, i3))
	facets = [utils.facet((nodes[e[0]], nodes[e[1]], nodes[e[2]]), **kw) for e in elems]
	if returnConnectivityTable:
		return facets, nodes, elems
	return facets


def ele(nodeFileName, eleFileName, shift=(0, 0, 0), scale=1.0, **kw):
	""" Import tetrahedral mesh from .ele file, return list of created tetrahedrons.

		:param string nodeFileName: name of .node file
		:param string eleFileName: name of .ele file
		:param (float,float,float)|Vector3 shift: (X,Y,Z) parameter moves the specimen.
		:param float scale: factor scales the given data.
		:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.polyhedron`
	"""
	f = open(nodeFileName)
	line = f.readline()
	while line.startswith('#'):
		line = f.readline()
	ls = line.split()
	nVertices = int(ls[0])
	if int(ls[1]) != 3:
		raise RuntimeError("wrong .node file, number of dimensions should be 3")
	vertices = [None for i in range(nVertices)]
	shift = Vector3(shift)
	for i in range(nVertices):
		line = f.readline()
		while line.startswith('#'):
			line = f.readline()
		ls = line.split()
		if not ls:
			continue
		v = shift + scale * Vector3(tuple(float(ls[j]) for j in (1, 2, 3)))
		vertices[int(ls[0]) - 1] = v
	f.close()
	#
	f = open(eleFileName)
	line = f.readline()
	while line.startswith('#'):
		line = f.readline()
	ls = line.split()
	if int(ls[1]) != 4:
		raise RuntimeError("wrong .ele file, unsupported tetrahedra's number of nodes")
	nTetras = int(ls[0])
	tetras = [None for i in range(nTetras)]
	for i in range(nTetras):
		ls = f.readline().split()
		tetras[int(ls[0]) - 1] = utils.polyhedron([vertices[int(ls[j]) - 1] for j in (1, 2, 3, 4)], **kw)
	f.close()
	return tetras


def textPolyhedra(fileName, material, shift=Vector3.Zero, scale=1.0, orientation=Quaternion((0, 1, 0), 0.0), **kw):
	"""Load polyhedra from a text file.
	
	:param str filename: file name. Expected file format is the one output by export.textPolyhedra.
	:param [float,float,float] shift: [X,Y,Z] parameter moves the specimen.
	:param float scale: factor scales the given data.
	:param quaternion orientation:  orientation of the imported polyhedra
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.polyhedra_utils.polyhedra`
	:returns: list of polyhedras.

	Lines starting with # are skipped
	"""
	from yade import polyhedra_utils
	infile = open(fileName, "r")
	lines = infile.readlines()
	infile.close()
	ret = []
	i = -1
	while (i < (len(lines) - 1)):
		i += 1
		line = lines[i]
		data = line.split()
		if (data[0][0] == "#"):
			continue

		if (len(data) != 3):
			raise RuntimeError("Check polyhedra input file! Number of parameters in the first line is not 3!")
		else:
			vertLoad = []
			#			ids = int(data[0])
			verts = int(data[1])
			surfs = int(data[2])
			i += 1
			for d in range(verts):
				dataV = lines[i].split()
				pos = orientation * Vector3(float(dataV[0]) * scale, float(dataV[1]) * scale, float(dataV[2]) * scale) + shift
				vertLoad.append(pos)
				i += 1
			polR = polyhedra_utils.polyhedra(material=material, v=vertLoad, **kw)
			if (polR != -1):
				ret.append(polR)
			i = i + surfs - 1
	return ret


# blockMeshDict ==================================================

from dataclasses import dataclass, field
from typing import NamedTuple

BOUNDARY_ERROR = 0
BOUNDARY_PATCH = 1
BOUNDARY_WALL = 2
BOUNDARY_EMPTY = 3


class BoundaryFace3(NamedTuple):
	v0: int
	v1: int
	v2: int

	def __str__(self) -> str:
		return f"({self.v0}, {self.v1}, {self.v2})"


class BoundaryFace4(NamedTuple):
	v0: int
	v1: int
	v2: int
	v3: int

	def __str__(self) -> str:
		return f"({self.v0}, {self.v1}, {self.v2}, {self.v3})"


@dataclass
class Boundary:
	name: str = ""
	typ: int = BOUNDARY_ERROR
	faces4: [BoundaryFace4] = field(default_factory=list)

	def __str__(self) -> str:
		s = f"{self.name} ({self.typ}):\n"
		s += f"\tfaces4({len(self.faces4)})\n"
		for f in self.faces4:
			s += f"\t  {f}\n"
		return s


def tryParseFloat(x: str, n: int) -> float:
	val: float
	try:
		val = float(x)
	except ValueError:
		assert False, f"{n}: Expected 'float', got: {x}"

	return val


def tryParseInt(x: str, n: int) -> int:
	val: int
	try:
		val = int(x)
	except ValueError:
		assert False, f"{n}: Expected 'int', got: {x}"

	return val


def blockMeshDict(path, patchasWall=True, emptyasWall=True, **kw):
	"""Load openfoam's blockMeshDict file's "boundary" section into facets.

	:param str path: file name. Typical value is: "system/blockMeshDic".
	:param bool patchasWall: load "patch"-es as walls.
	:param bool emptyasWall: load "empty"-es as walls.
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.facet`
	:returns: list of facets.
    	"""
	convertToMeters = 1.0

	verticesBlock = False
	vertices = []
	boundariesBlock = False
	boundaries = []
	facesBlock = False
	facets = []

	lines = []
	with open(path) as fp:
		lines = fp.readlines()

	lineNumber = 0
	blockDepth = 0

	currentBoundary = Boundary()

	for line in lines:
		line = line.strip()
		lineNumber += 1

		if len(line) == 0:
			continue

		if line.startswith("convertToMeters") and blockDepth == 0:
			convertToMetersStr = line.split(" ")[-1]
			assert convertToMetersStr[-1] == ';', f"{lineNumber}: Expected ';' at the end of convertToMeters"
			convertToMeters = tryParseFloat(convertToMetersStr[:-1], lineNumber)
		elif line.startswith("vertices") and blockDepth == 0:
			verticesBlock = True
		elif line.startswith("boundary") and blockDepth == 0:
			boundariesBlock = True
		elif line == "(" or line == "{":
			blockDepth += 1
			continue
		elif line.endswith(");"):
			blockDepth -= 1

			if blockDepth == 0:
				verticesBlock = False
				boundariesBlock = False
			elif facesBlock == 2:
				facesBlock = False

			continue
		elif line.endswith("}"):
			blockDepth -= 1

			if blockDepth == 1:
				assert boundariesBlock, f"{lineNumber}: Parser error, only boundaries are supported."
				assert currentBoundary.name != "", f"{lineNumber}: Empty name for boundary"
				assert currentBoundary.typ != BOUNDARY_ERROR, f"{lineNumber}: Invalid type for boundary, supported: patch, wall, empty"
				assert len(currentBoundary.faces4) > 0, f"{lineNumber}: Boundary must contain at least one face"

				boundaries.append(currentBoundary)

				currentBoundary = Boundary()

			continue
		else:
			if verticesBlock:
				if line.find("(") != -1 and line.find(")") != -1:
					vStr = line[1:-1]
					vStrs = vStr.split(" ")
					assert len(vStrs) == 3, f"{lineNumber}: Vertice format expected: (x y z), got: {vStr}"
					x = tryParseFloat(vStrs[0], lineNumber)
					y = tryParseFloat(vStrs[1], lineNumber)
					z = tryParseFloat(vStrs[2], lineNumber)

					vertices.append((x, y, z))
				else:
					assert False, f"{lineNumber}: Vertice format expected: (x y z), got: {line}"
			elif boundariesBlock:
				if blockDepth == 1:
					assert len(line.split(" ")) == 1, f"{lineNumber}: Expected surface name, got: {line}"

					currentBoundary.name = line
				elif blockDepth == 2:
					if line.startswith("type"):
						typStrs = line.split(" ")
						assert len(typStrs) == 2, f"{lineNumber}: Expected boundary type definition, got: {line}"
						assert typStrs[1][-1] == ';', f"{lineNumber}: Expected ';' at the end of boundary type"

						typStr = typStrs[1][:-1]
						if typStr == "patch":
							currentBoundary.typ = BOUNDARY_PATCH
						elif typStr == "wall":
							currentBoundary.typ = BOUNDARY_WALL
						elif typStr == "empty":
							currentBoundary.typ = BOUNDARY_EMPTY
						else:
							assert False, f"{lineNumber}: Unknown boundary type: {typStr}"

					elif line.startswith("faces"):
						facesBlock = True
					else:
						assert False, f"{lineNumber}: Expected 'type' or 'faces', got: {line}"
				elif blockDepth == 3:
					assert facesBlock, f"{lineNumber}: Parser error, only faces supported at level 3"

					if line.find("(") != -1 and line.find(")") != -1:
						fStr = line[1:-1]
						fStrs = fStr.split(" ")
						assert len(fStrs) == 4, f"{lineNumber}: Face format expected: (v0 v1 v2 v3), got: {fStr}"

						v0 = tryParseInt(fStrs[0], lineNumber)
						v1 = tryParseInt(fStrs[1], lineNumber)
						v2 = tryParseInt(fStrs[2], lineNumber)
						v3 = tryParseInt(fStrs[3], lineNumber)

						assert v0 >= 0, f"{lineNumber}: Face index must be greater or equal to 0, not {v0}"
						assert v1 >= 0, f"{lineNumber}: Face index must be greater or equal to 0, not {v1}"
						assert v2 >= 0, f"{lineNumber}: Face index must be greater or equal to 0, not {v2}"
						assert v3 >= 0, f"{lineNumber}: Face index must be greater or equal to 0, not {v3}"

						face = BoundaryFace4(v0, v1, v2, v3)
						currentBoundary.faces4.append(face)
					else:
						assert False, f"{lineNumber}: Face format expected: (v0 v1 v2 v3), got: {line}"
				else:
					assert blockDepth == 0, f"{lineNumber}: Parser error, only levels 0, 1, 2 and 3 are supported"

	for b in boundaries:
		if b.typ == BOUNDARY_PATCH and patchasWall == False:
			continue

		if b.typ == BOUNDARY_EMPTY and emptyasWall == False:
			continue

		for f4 in b.faces4:
			f0 = BoundaryFace3(f4.v0, f4.v1, f4.v2)
			f1 = BoundaryFace3(f4.v2, f4.v3, f4.v0)

			f0v0 = tuple([x * convertToMeters for x in vertices[f0.v0]])
			f0v1 = tuple([x * convertToMeters for x in vertices[f0.v1]])
			f0v2 = tuple([x * convertToMeters for x in vertices[f0.v2]])

			f1v0 = tuple([x * convertToMeters for x in vertices[f1.v0]])
			f1v1 = tuple([x * convertToMeters for x in vertices[f1.v1]])
			f1v2 = tuple([x * convertToMeters for x in vertices[f1.v2]])

			facets.append(utils.facet((f0v0, f0v1, f0v2), **kw))
			facets.append(utils.facet((f1v0, f1v1, f1v2), **kw))

	return facets
