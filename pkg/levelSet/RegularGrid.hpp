/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once
#include <lib/serialization/Serializable.hpp>

namespace yade {

class RegularGrid : public Serializable {
public:
	void pyHandleCustomCtorArgs(
	        boost::python::tuple& t,
	        boost::python::dict&) override; // way to go for constructor variants: adding a .def(boost::python::init<..>) does not seem to work
	//	adopting const member functions for the following, so that the compiler knows they should not change anything of *this:
	Vector3i closestCorner(const Vector3r&) const; // the i,j,k grid indices of the closest "smaller" grid point
	Vector3r getDims() const;
	Vector3r gridPoint(int, int, int) const;
	//	Vector3r gridPoint(Vector3i) const; // would make sense for the C++ world, but would require extra work for the python exposure to work (as an overloaded C++ function)
	Vector3r max() const;

	// clang-format off
  YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(RegularGrid,Serializable,"A rectilinear (aka uniform or regular) grid, for :yref:`LevelSet` shapes or other purposes. A cubic regular grid extending from a :yref:`min<RegularGrid.min>` = (m,m,m) to a max at (M,M,M), with :yref:`nGP<RegularGrid.nGP>` = (n,n,n) ie a :yref:`spacing<RegularGrid.spacing>` = (M-m)/(n-1), can be conveniently obtained from RegularGrid(m,M,n)",
		((Vector3r,min,Vector3r(NaN,NaN,NaN),Attr::readonly,"The minimum corner of the grid."))
		((Vector3i,nGP,Vector3i::Zero(),Attr::readonly,"The number of grid points along the three axes as a Vector3i."))
		((Real,spacing,-1,Attr::readonly,"The (uniform and isotropic) grid spacing between two axis-consecutive grid points."))
		,
		,
		.def("gridPoint",&RegularGrid::gridPoint,(boost::python::args("i", "j", "k")),"Returns the Vector3 position of any grid point, given its indices *i* (along the X-axis), *j* (Y-axis), *k* (Z-axis).")
		.def("max",&RegularGrid::max,"Returns the maximum corner of the grid.")
		.def("dims",&RegularGrid::getDims,"Returns the grid dimensions along the three axes, as a Vector3.")
		.def("closestCorner",&RegularGrid::closestCorner,(boost::python::arg("pt")),"Returns the Vector3i indices of the closest gridpoint which is smaller (for all components) than *pt*.")
				    );
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(RegularGrid);
};     // namespace yade
#endif // YADE_LS_DEM
