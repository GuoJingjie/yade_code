/*************************************************************************
*  Copyright (C) 2010 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "FlowBoundingSphere.hpp"

#ifdef FLOW_ENGINE

namespace yade { // Cannot have #include directive inside.


namespace CGT {

	template <class _Tesselation> class PeriodicFlowLinSolv : public FlowBoundingSphereLinSolv<_Tesselation, PeriodicFlow<_Tesselation>> {
	public:
		typedef _Tesselation         Tesselation;
		typedef Network<Tesselation> _N;
		DECLARE_TESSELATION_TYPES(Network<Tesselation>)
		typedef FlowBoundingSphereLinSolv<_Tesselation, PeriodicFlow<_Tesselation>> BaseFlowSolver;
		typedef typename BaseFlowSolver::ETriplet                                   ETriplet;

		///painfull, but we need that for templates inheritance...
		using _N::boundaries;
		using _N::boundingCells;
		using _N::boundsIds;
		using _N::cornerMax;
		using _N::cornerMin;
		using _N::currentTes;
		using _N::debugOut;
		using _N::facetNFictious;
		using _N::facetVertices;
		using _N::Height;
		using _N::idOffset;
		using _N::nOfSpheres;
		using _N::num_particles;
		using _N::Rmoy;
		using _N::sectionArea;
		using _N::sSolidTot;
		using _N::T;
		using _N::vPoral;
		using _N::vPoralPorosity;
		using _N::VSolidTot;
		using _N::vtkInfiniteCells;
		using _N::vtkInfiniteVertices;
		using _N::vTotal;
		using _N::Vtotalissimo;
		using _N::vTotalPorosity;
		using _N::xMax;
		using _N::xMaxId;
		using _N::xMin;
		using _N::xMinId;
		using _N::yMax;
		using _N::yMaxId;
		using _N::yMin;
		using _N::yMinId;
		using _N::zMax;
		using _N::zMaxId;
		using _N::zMin;
		using _N::zMinId;
		//same for functions
		using _N::addBoundingPlanes;
		using _N::boundary;
		using _N::defineFictiousCells;
		;
		using _N::surfaceSolidThroatInPore;
		using _N::tesselation;

		using BaseFlowSolver::checkSphereFacetOverlap;
		using BaseFlowSolver::clampKValues;
		using BaseFlowSolver::distanceCorrection;
		using BaseFlowSolver::factorizeOnly;
		using BaseFlowSolver::fluidBulkModulus;
		using BaseFlowSolver::kFactor;
		using BaseFlowSolver::KOptFactor;
		using BaseFlowSolver::maxKdivKmean;
		using BaseFlowSolver::meanKStat;
		using BaseFlowSolver::minKdivKmean;
		using BaseFlowSolver::minPermLength;
		using BaseFlowSolver::noCache;
		using BaseFlowSolver::permeabilityMap;
		using BaseFlowSolver::rAverage;
		using BaseFlowSolver::relax;
		using BaseFlowSolver::resetRHS;
		using BaseFlowSolver::tolerance;
		using BaseFlowSolver::viscosity;
		/// More members from LinSolv variant
		using BaseFlowSolver::A;
		using BaseFlowSolver::areCellsOrdered;
		using BaseFlowSolver::bodv;
		using BaseFlowSolver::errorCode;
		using BaseFlowSolver::fullAcolumns;
		using BaseFlowSolver::fullAvalues;
		using BaseFlowSolver::gsB;
		using BaseFlowSolver::gsdV;
		using BaseFlowSolver::gsP;
		using BaseFlowSolver::isFullLinearSystemGSSet;
		using BaseFlowSolver::isLinearSystemSet;
		using BaseFlowSolver::ncols;
		using BaseFlowSolver::orderedCells;
		using BaseFlowSolver::T_b;
		using BaseFlowSolver::T_bv;
		using BaseFlowSolver::T_cells;
		using BaseFlowSolver::T_index;
		using BaseFlowSolver::T_nnz;
		using BaseFlowSolver::T_x;
		using BaseFlowSolver::tripletList;
		using BaseFlowSolver::useSolver;
		using BaseFlowSolver::xodv;

		vector<int> indices; //redirection vector containing the rank of cell so that T_cells[indices[cell->info().index]]=cell

		virtual ~PeriodicFlowLinSolv();
		PeriodicFlowLinSolv();

		///Linear system solve
		int setLinearSystem(Real dt = 0) override;
		int setLinearSystemFullGS(Real dt = 0) override;
	};

} //namespace CGTF

}; // namespace yade

#include "PeriodicFlowLinSolv.ipp"
#endif //FLOW_ENGINE
