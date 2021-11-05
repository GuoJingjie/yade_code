/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#pragma once

#include <lib/base/Logging.hpp>
#include <lib/serialization/Serializable.hpp>
#include <pkg/levelSet/RegularGrid.hpp>

namespace yade {

class FastMarchingMethod : public Serializable { // let s derive from Serializable for Python access
private:
	enum gpState { knownState, trialState, farState };
	// Definitions:
	// - knownState: that gp has a distance value we're sure about / can not modify anymore anyway
	// - trialState: for a gp in the narrow band, with at least one "known" neighbour. It carries some finite (no longer infinite) distance value we re unsure of
	// - farState: just the initial state for all gridpoints
	vector<vector<vector<gpState>>> gpStates;
	vector<Vector3i>                knownTmp, // a temporary version of known (below), which we will empty at some point for the purpose of the algorithm
	        trials;                           // the narrow band = all gps (their indices) with a trialState
	vector<vector<vector<Real>>> phiField;
	void                         printNeighbValues(int, int, int) const;
	//	bool compFnIn(Vector3i, Vector3i), compFnOut(Vector3i, Vector3i);
	Real                          eikDiscr(Real, Real, Real) const;
	Real                          eikDiscr(Real, Real, Real, Real) const;
	Real                          phiAtNgbr(int idx, int i, int j, int k) const;
	Real                          phiWhenKnown(int i, int j, int k, bool exterior) const;
	Real                          phiFromEik(Real, Real, Real, bool) const;
	Real                          phiFromEik(Real, Real, Real, Real, bool) const;
	void                          iniFront(bool);
	void                          iniStates();
	void                          updateFastMarchingMethod(int i, int j, int k, bool);
	void                          loopTrials(bool);
	std::pair<vector<Real>, Real> surroundings(int, int, int, bool) const;
	void                          trializeFromKnown(int xInd, int yInd, int zInd, bool); // indices here refer to some known gridpoint
	void                          trialize(int xInd, int yInd, int zInd, bool);          // indices here refer to the actual trial-to-be gridpoint
	void                          confirm(int xInd, int yInd, int zInd, Real phiVal, bool, bool checkState);
	//	see https://stackoverflow.com/questions/1902311/problem-sorting-using-member-function-as-comparator for the use of struct below
	struct downwindSort { // will give a sorted items with the 1st element being the closest
		downwindSort(const FastMarchingMethod& info)
		        : m_info(info)
		{
		} // only if you really need the object state
		const FastMarchingMethod& m_info;
		bool                      operator()(const Vector3i& gp1, const Vector3i& gp2)
		{
			Real phi1(m_info.phiField[gp1[0]][gp1[1]][gp1[2]]), phi2(m_info.phiField[gp2[0]][gp2[1]][gp2[2]]);
			return math::abs(phi1) < math::abs(phi2); // true if gp1 is the closest one to the surface, in both (inside and outside) cases
		}
	};
	struct upwindSort { // will sort such that 1st element is the farthest from surface => closest is the last one
		upwindSort(const FastMarchingMethod& info)
		        : m_info(info)
		{
		} // only if you really need the object state
		const FastMarchingMethod& m_info;
		bool                      operator()(const Vector3i& gp1, const Vector3i& gp2)
		{
			Real phi1(m_info.phiField[gp1[0]][gp1[1]][gp1[2]]), phi2(m_info.phiField[gp2[0]][gp2[1]][gp2[2]]);
			return math::abs(phi1) > math::abs(phi2);
		}
	};

public:
	DECLARE_LOGGER;                     // see https://yade-dem.org/doc/prog.html#debug-macros
	vector<vector<vector<Real>>> phi(); // kind-of the "main"
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(FastMarchingMethod,Serializable,"Executes a Fast Marching Method (FMM) to solve $||\\vec \\nabla \\phi|| = c$ for a discrete field $\\phi$ defined on :yref:`grid<FastMarchingMethod.grid>`, with :yref:`phiIni<FastMarchingMethod.phiIni>` serving as boundary condition. Typically, $c=1$ (see :yref:`speed<FastMarchingMethod.speed>`) and $\\phi$ is a distance field. Note that the minimum search inherent to the FMM is not yet optimal in terms of execution speed and faster implementations of the FMM may be found elsewhere. See [Duriez2021b]_ for more details, where the class was coined DistFMM.",
		((vector<Vector3i>,known,,Attr::readonly,"Gridpoints (indices) with distance known for good: they have been at some point the shortest gp to the surface while executing the FMM."))
		((vector<vector<vector<Real>>>,phiIni,,,"Initial discrete field defined on the :yref:`grid<FastMarchingMethod.grid>` that will serve as a boundary condition for the FMM. Field values have to be - inf (resp. inf) for points being far inside (resp. outside) and correct (finite) on each side of the interface. Built-in functions *distIniSE* (for superellipsoids), *phiIniCppPy* (for a Python user function, through a mixed C++-Py internal implementation) or *phiIniPy* (for a Python user function through a pure Py internal implementation) may be used for such a purpose."))
		((shared_ptr<RegularGrid>,grid,,,"The underlying :yref:`regular grid<RegularGrid>`."))
		((Real,speed,1,,"Keep to 1 for a true distance, 2 for the flake-like rose verification of [Duriez2021b]_."))
		,
		,.def("phi",&FastMarchingMethod::phi,"Executes the FMM and returns its solution as a list of list of list, with the [i][j][k] element corresponding to grid.gridPoint(i,j,k)")
		)
	// clang-format on
};
REGISTER_SERIALIZABLE(FastMarchingMethod);
} // namespace yade
#endif
