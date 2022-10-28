/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once

#include <lib/base/Logging.hpp>
#include <lib/base/Math.hpp>
#include <core/Clump.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/levelSet/LevelSet.hpp>

namespace yade { // Cannot have #include directive inside.
class ShopLS {
public:
	DECLARE_LOGGER;
	//****** Helpful functions for RegularGrid creation *******
	//NB: I could overload in C++ the two following functions but then I would not know how to expose them to Python. Here, with two different C++ names, and one Python name (and two .def in _utils.cpp), it works...
	//See e.g. https://www.boost.org/doc/libs/1_66_0/libs/python/doc/html/tutorial/tutorial/functions.html#tutorial.functions.overloading https://stackoverflow.com/questions/26200998/how-to-wrap-functions-overloaded-by-type for some reading
	static int      nGPr(Real, Real, Real);                              // to have a good number of gridpoints
	static Vector3i nGPv(const Vector3r&, const Vector3r&, const Real&); // "overloaded" (not in C++, will be in Python) version of the above
	//****** Cartesian <-> spherical coordinates functions *******
	static Vector3r cart2spher(const Vector3r&);
	static Vector3r spher2cart(const Vector3r&);
	//****** Functions for LS bodies creation *******
	static shared_ptr<LevelSet> lsSimpleShape(
	        int,
	        const AlignedBox3r&,
	        const Real&,
		const Real&,
	        const Vector2r&,
	        shared_ptr<
	                Clump>); // will be passed to Python in eg levelSetBody() and "returning shared_ptr<…> objects is the preferred way of passing objects from c++ to python" according to https://yade-dem.org/doc/prog.html#reference-counting
	static Real
	distApproxSE(const Vector3r& pt, const Vector3r& extents, const Vector2r& epsilons); // the approximated distance function to a superellipsoid
	static vector<vector<vector<Real>>>
	        distIniClump(shared_ptr<Clump>, shared_ptr<RegularGrid>); // appropriate FastMarchingMethod.phiIni for the distance to a clump
	static vector<vector<vector<Real>>>
	                distIniSE(const Vector3r&, const Vector2r&, shared_ptr<RegularGrid>); // appropriate FastMarchingMethod.phiIni for the distance to a se
	static Real     distToCircle(const Vector3r&, const Real&);                           // the distance to a circle's contour, in a (x,y) plane
	static Real     distToClumpSpheres(Vector3r, shared_ptr<Clump>);                      // minimum distance to all Clump members
	static Real     distToInterval(Real, Real, Real);                                     // the distance to an interval in 1D space
	static Real     distToRecParallelepiped(Vector3r pt, Vector3r extents);               // the distance to a box
	static Real     distToSph(const Vector3r&, const Real&, Vector3r center = Vector3r::Zero());          // the distance to a sphere surface
	static bool     insideClump(const Vector3r&, shared_ptr<Clump>);                                      // could fit into core/Clump ?
	static Real     insideOutsideSE(const Vector3r& pt, const Vector3r& radii, const Vector2r& epsilons); // inside-outside function to a superellipsoid
	static Real     fioRose(Vector3r);
	static Vector3r grad_fioRose(Vector3r);
	static Real     distApproxRose(Vector3r);
	static vector<vector<vector<Real>>>
	                                    phiIni(int,
	                                           Vector3r,
	                                           Vector2r,
	                                           shared_ptr<Clump>,
	                                           shared_ptr<RegularGrid>); // for constructing FastMarchingMethod.phiIni in various cases, see implementation. Not intended to be Python-exposed
	static vector<vector<vector<Real>>> phiIniCppPy(shared_ptr<RegularGrid>); // handy function for constructing a Py-related FastMarchingMethod.phiIni
	//****** Miscellaneous functions ********
	static Vector3r rigidMapping(const Vector3r&, const Vector3r&, const Vector3r&, const Quaternionr&); // a rigid body mapping
	static Real     biInterpolate(
                std::array<Real, 2>,
                std::array<Real, 2>,
                std::array<Real, 2>,
                std::array<std::array<Real, 2>, 2>); // bi-interpolation in a plane, used in LevelSet.distance()
	static shared_ptr<ScGeom>
	                          geomPtr(Vector3r                       ctctPt,
	                                  Real                           un,
	                                  Real                           rad1,
	                                  Real                           rad2,
	                                  const State&                   rbp1,
	                                  const State&                   rbp2,
	                                  const shared_ptr<Interaction>& c,
	                                  const Vector3r&                currentNormal,
	                                  const Vector3r&                shift2);
	static shared_ptr<ScGeom> geomPtrForLaterRemoval(const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c);
};
} // namespace yade
#endif //YADE_LS_DEM
