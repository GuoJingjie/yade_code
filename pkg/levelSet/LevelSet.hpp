/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <lib/computational-geometry/MarchingCube.hpp>
#include <core/Shape.hpp>
#include <pkg/levelSet/RegularGrid.hpp>

namespace yade {

class LevelSet : public Shape {
private:
	Real     minRad, maxRad; // for sphericity
	Vector3r center;
	Real     volume, lengthChar;
	Vector3r inertia; // the eigenvalues of the inertia matrix: its diagonal expression in localAxes basis (in a Vector3r form here)
	bool     initDone;
	bool     initDoneMarchingCubes;
	int      nVoxInside;
	void     init();          // compute nVoxInside, center, volume, and inertia and calls initSurfNodes
	void     initSurfNodes(); // fills surfNodes
	bool rayTraceInCell(const Vector3r&, const Vector3r&, const Vector3r&, const Vector3i&); // handles the ray tracing from a given point in a given cell
	void rayTrace(const Vector3r&); // recursively calls rayTraceInCell, walking accross the whole grid along a ray starting from center
	struct mcData {                 // Structure for holding marching cubes triangulation of level set particle
		vector<Vector3r> triangles;
		vector<Vector3r> normals;
		int              nbTriangles;
	};
	mcData marchingCubesData; // Actual marching cubes data holder
public:
	Real             distance(const Vector3r&) const; // gives through interpolation the distance from a point to the surface
	Vector3r         normal(const Vector3r&) const;   // gives the outwards normal at some point
	Real             getVolume();                     // these 3 get*() may call init() if not already done, they can not be const-declared
	Vector3r         getCenter();
	Vector3r         getInertia();
	Real             getSurface() const;           // this one can be const-declared
	void             computeMarchingCubes();       // Compute the marching cube triangulation for the LS shape
	vector<Vector3r> getMarchingCubeTriangles();   // Retrieve marching cube triangles
	vector<Vector3r> getMarchingCubeNormals();     // Retrieve marching cube normals
	int              getMarchingCubeNbTriangles(); // Retrieve marching cube number of triangles
	virtual ~LevelSet() {};
	// clang-format off
  YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(LevelSet,Shape,"A level set description of particle shape based on a :yref:`discrete distance field<LevelSet.distField>` and :yref:`surface nodes<LevelSet.surfNodes>` [Duriez2021a]_ [Duriez2021b]_.",
		((vector< vector< vector<Real> > >,distField,,Attr::readonly,"The signed (< 0 when inside) distance-to-surface function as a discrete scalar field on :yref:`lsGrid<LevelSet.lsGrid>`, with distField[i][j][k] corresponding to lsGrid.gridPoint(i,j,k). From Python, slice this multi-dimensional list with care: while distField[i][:][:] corresponds to values on a x-cst plane, distField[:][:][k] is not at z-constant (use [[distField[i][j][k] for j in ..] for i in ..] instead)"))
		((vector<Vector3r>,corners,,Attr::readonly,"The 8 corners of an axis-aligned bounding box, in local axes. It is computed once for all by :yref:`Bo1_LevelSet_Aabb` and used by the same Functor to get :yref:`Body.bound`."))
		((vector<Vector3r>,surfNodes,,Attr::readonly,"Surface discretization nodes (the list of) used for exact contact treatment in :yref:`Ig2_LevelSet_LevelSet_ScGeom`, previously coined boundNodes in [Duriez2021b]_. Expressed in local frame. Getting them back after a save/load cycle requires to launch one iteration or to first ask for shape.center.")) // NB: just nodes as a name would "conflict" with many PFacet variables
		((int,nSurfNodes,100,,"The number of boundary nodes in :yref:`surfNodes<LevelSet.surfNodes>`, previously coined nNodes in [Duriez2021b]_. Usually set through utils levelSetBody() function (has to be set at instantiation in all cases). Please use a perfect square + 2 if not :yref:`twoD<LevelSet.twoD>` and if :yref:`nodesPath<LevelSet.nodesPath>` = 1."))
		((int,nodesPath,2,,"Defines how the space of spherical coordinates $(\\theta \\in [0;\\pi] ,\\varphi\\in [0;2 \\pi])$ is discretized when ray tracing the boundary nodes: 1 gives a rectangular partition of that space, plus two nodes at $\\theta = 0 [\\pi]$; 2 locates the nodes along a spiral path [Duriez2021a]_")) // Elias'polyhedralBall code; and Rakhmanov1994
		((Real,nodesTol,50,,"Tolerance coefficient for accepting (if $|\\phi| / L <$ nodesTol $\\times$ numeric precision with $\\phi$ the return value of :yref:`distance<LevelSet.distance>` and $L$ a body-characteristic length taken as $\\sqrt[3]{V}$ with $V$ the :yref:`volume<LevelSet.volume>`, or $\\sqrt{V/g}$ with $g$ the grid :yref:`spacing<RegularGrid.spacing>` if :yref:`twoD<LevelSet.twoD>`) boundary nodes proposed by the ray tracing algorithm.")) // phi will be a varphi, because of a let\phi\varphi in doc/sphinx/conf.py ?
		((Real,sphericity,-1,Attr::readonly,"Shape sphericity computed from boundary nodes and assuming both largest inscribed sphere and smallest circumscribed sphere have the origin (of local axes) as center."))
		((shared_ptr<RegularGrid>,lsGrid,new RegularGrid,Attr::readonly,"The :yref:`regular grid<RegularGrid>` carrying :yref:`distField<LevelSet.distField>`, in local axes."))
		((bool,twoD,false,Attr::readonly,"True for z-invariant shapes. Serves to restrict the definition of :yref:`surfNodes<LevelSet.surfNodes>` in the (x,y) plane."))
		,
		minRad = std::numeric_limits<Real>::infinity();
		maxRad = 0;
		inertia = Vector3r::Zero();
		initDone = false; // after hesitation, it is finally chosen to save the least of data, but to recall init() after a save/load
		initDoneMarchingCubes = false;
		lengthChar = -1;
		volume = -1;
		nVoxInside = -1;
		center = Vector3r(std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::infinity());
		createIndex(); // necessary for such a Shape-derived class, see https://yade-dem.org/doc/prog.html#indexing-dispatch-types
 		,
		.def("volume",&LevelSet::getVolume,"The volume defined by the negative domain of the :yref:`level set function<LevelSet.distField>`, in a voxellised fashion. A voxel is said to be inside according to the level set value at its minimum grid point.")
		.def("center",&LevelSet::getCenter,"The center of mass of the :yref:`volume<LevelSet.volume>` (considering obviously an uniform density for this volume), in local axes (for verification purposes, by comparison with the origin).")
 		.def("inertia",&LevelSet::getInertia,"The eigenvalues of the geometric inertia matrix (the one considering the infinitesimal volume as the integrand, instead of infinitesimal mass) as a Vector3r.")
// 		.def("nodesInCell",&LevelSet::getNodesInCellCube,(boost::python::args("i", "j", "k")),"Which boundary nodes belong to a given grid cube (given by its i,j,k indices)")
		.def("distance",&LevelSet::distance,(boost::python::arg("pt")),"Distance to surface value at pt, pt being expressed in local frame.")
		.def("normal",&LevelSet::normal,(boost::python::arg("pt")),"Normal vector to the surface, at some pt. Local frame applies to both output normal and input pt.")
		.def("rayTrace",&LevelSet::rayTrace,(boost::python::arg("ray")),"Performs one ray tracing, possibly modifying :yref:`surfNodes<LevelSet.surfNodes>`. Provided for debugging purposes")
		.def("getSurface",&LevelSet::getSurface,"Returns particle surface as computed from numeric integration over the :yref:`surface nodes<LevelSet.surfNodes>`. Requires :yref:`nodesPath<LevelSet.nodesPath>` = 1.")
		.def("computeMarchingCubes",&LevelSet::computeMarchingCubes,"Compute or recompute the triangulation of the particle surface after using the Marching Cubes algorithm on :yref:`distField<LevelSet.distField>`.")
		.def("marchingCubesVertices",&LevelSet::getMarchingCubeTriangles,"Returns the vertices for a surface triangulation obtained after executing the Marching Cubes algorithm on :yref:`distField<LevelSet.distField>`.")
		.def("marchingCubesNormals",&LevelSet::getMarchingCubeNormals,"Returns the normals for a surface triangulation obtained after executing the Marching Cubes algorithm on :yref:`distField<LevelSet.distField>`.")
		.def("marchingCubesNbTriangles",&LevelSet::getMarchingCubeNbTriangles,"Returns the number of triangles forming the surface triangulation as per the Marching Cubes algorithm (executed on :yref:`distField<LevelSet.distField>`).")
	)
	// clang-format on
	REGISTER_CLASS_INDEX(LevelSet, Shape); // necessary for such a Shape-derived class, see https://yade-dem.org/doc/prog.html#indexing-dispatch-types
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(LevelSet);
} // namespace yade
#endif // YADE_LS_DEM
