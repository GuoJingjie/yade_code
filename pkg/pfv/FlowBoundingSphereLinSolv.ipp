/*************************************************************************
*  Copyright (C) 2010 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef FLOW_ENGINE

// #define XVIEW
#include "FlowBoundingSphereLinSolv.hpp" //include after #define XVIEW
#if CGAL_VERSION_NR < CGAL_VERSION_NUMBER(4, 11, 0)
#include "CGAL/constructions/constructions_on_weighted_points_cartesian_3.h"
#endif
#include "vector"
#include <CGAL/Width_3.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <new>
#include <sys/stat.h>
#include <sys/types.h>
#include <utility>

#ifdef XVIEW
//#include "Vue3D.h" //FIXME implicit dependencies will look for this class (out of tree) even ifndef XVIEW
#endif

#ifdef YADE_OPENMP
#include <omp.h>
// #define GS_OPEN_MP //It should never be defined if Yade is not using openmp
#endif

// #define PARDISO //comment this if pardiso lib is not available

#ifdef LINSOLV
extern "C" {
void openblas_set_num_threads(int num_threads);
}
#endif

namespace yade { // Cannot have #include directive inside.
namespace CGT {

#ifdef PARDISO
#ifdef AIX
#define F77_FUNC(func) func
#else
#define F77_FUNC(func) func##_
#endif
	/* PARDISO prototype. */
	// Note about using high precisision Real, for start move these prototypes to file lib/compatibility/PardisoCompatibility.hpp,
	// and do similar wrapper/converter like lib/compatibility/LapackCompatibility.cpp is done
	// then try to replace it with native Eigen high-precision Real roooutines.
	extern "C" int F77_FUNC(pardisoinit)(void*, int*, int*, int*, double*, int*);

	extern "C" int F77_FUNC(pardiso)(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);
#endif

#ifdef XVIEW
	Vue3D Vue1;
#endif
	template <class _Tesselation, class FlowType> FlowBoundingSphereLinSolv<_Tesselation, FlowType>::~FlowBoundingSphereLinSolv()
	{
#ifdef TAUCS_LIB
		if (Fccs) taucs_ccs_free(Fccs);
#endif
#ifdef SUITESPARSE_VERSION_4
		if (useSolver == 4) {
			if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);
			CHOLMOD(free_sparse)(&Achol, &com);
			CHOLMOD(free_factor)(&L, &com);
			CHOLMOD(finish)(&com);
			if (getCHOLMODPerfTimings) {
				gettimeofday(&end, NULL);
				cout << "CHOLMOD Time to finalize multithreaded com "
				     << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) << endl;
			}
		}
#endif
	}
	template <class _Tesselation, class FlowType>
	FlowBoundingSphereLinSolv<_Tesselation, FlowType>::FlowBoundingSphereLinSolv()
	        : FlowType()
	{
		useSolver = 0;
		isLinearSystemSet = 0;
		isFullLinearSystemGSSet = 0;
		areCellsOrdered = 0; //true when orderedCells is filled, turn it false after retriangulation
		updatedRHS = false;
		ZERO = 0;
#ifdef TAUCS_LIB
		T_A = &SystemMatrix;
		F = NULL;    //The taucs factor
		Fccs = NULL; //The taucs factor in CCS format
#endif
		pardisoInitialized = false;
		pTimeInt = 0;
		pTime1N = 0;
		pTime2N = 0;
		pTime1 = 0;
		pTime2 = 0;
#ifdef LINSOLV
		factorizedEigenSolver = false;
		numFactorizeThreads = 1;
		numSolveThreads = 1;
#endif
#ifdef SUITESPARSE_VERSION_4
		CHOLMOD(start)(&com);
		//CHOLMOD(wildcard)();
		factorExists = false;
		com.nmethods = 1;                       // nOrderingMethods; //1;
		com.method[0].ordering = CHOLMOD_METIS; // orderingMethod; //CHOLMOD_METIS;
#ifdef PFV_GPU
		if (multithread) com.maxGpuMemFraction = 0.4; //using (less than) half of the available memory for each solver
#endif
		com.supernodal = CHOLMOD_AUTO; //CHOLMOD_SUPERNODAL;
#endif
		reuseOrdering = false;
	}


	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::swapFwd(Real* v, int i2)
	{
		// FIXME: why not write std::swap(v[i2] , v[i2+1]); ? // it can be optimized even more than these three lines.
		Real temp = v[i2];
		v[i2] = v[i2 + 1];
		v[i2 + 1] = temp;
	}
	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::swapFwd(int* v, int i2)
	{
		Real temp = v[i2];
		v[i2] = v[i2 + 1];
		v[i2 + 1] = temp;
	}

//spatial sort traits to use with a pair of CGAL::sphere pointers and integer.
//template<class _Triangulation>
// usage : spatial_sort(pointsPtrs.begin(),pointsPtrs.end(), RTraits_for_spatial_sort()/*, CGT::RTriangulation::Weighted_point*/);
#include <utility>
	template <class Triangulation> struct CellTraits_for_spatial_sort : public Triangulation::Geom_traits {
		typedef typename Triangulation::Geom_traits                 Gt;
		typedef const typename Triangulation::Finite_cells_iterator Point_3;
		struct Less_x_3 {
			bool operator()(const Point_3& p, const Point_3& q) const { return typename Gt::Less_x_3()(p->info(), q->info()); }
		};
		struct Less_y_3 {
			bool operator()(const Point_3& p, const Point_3& q) const { return typename Gt::Less_y_3()(p->info(), q->info()); }
		};
		struct Less_z_3 {
			bool operator()(const Point_3& p, const Point_3& q) const { return typename Gt::Less_z_3()(p->info(), q->info()); }
		};
		Less_x_3 less_x_3_object() const { return Less_x_3(); }
		Less_y_3 less_y_3_object() const { return Less_y_3(); }
		Less_z_3 less_z_3_object() const { return Less_z_3(); }
	};

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::resetLinearSystem()
	{
		FlowType::resetLinearSystem();
		isLinearSystemSet = false;
		isFullLinearSystemGSSet = false;
		areCellsOrdered = false;
#ifdef TAUCS_LIB
		if (F) taucs_supernodal_factor_free(F);
		F = NULL;
		if (Fccs) taucs_ccs_free(Fccs);
		Fccs = NULL;
#endif
#ifdef LINSOLV
		factorizedEigenSolver = false;
#endif
#ifdef PARDISO
		if (pardisoInitialized) {
			phase = -1;
			F77_FUNC(pardiso)
			(pt, &maxfct, &mnum, &mtype, &phase, &ncols, &ddum, NULL, NULL, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
			pardisoInitialized = false;
		}
#endif
	}

	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::setLinearSystem(Real dt)
	{
#ifdef SUITESPARSE_VERSION_4
		if (!multithread && factorExists && useSolver == 4) {
			if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);
			CHOLMOD(free_sparse)(&Achol, &com);
			//cholmod_l_free_triplet(&cholT, &com);
			if (!reuseOrdering) {
				CHOLMOD(free_factor)(&L, &com);
				CHOLMOD(finish)(&com);
				if (getCHOLMODPerfTimings) {
					gettimeofday(&end, NULL);
					cout << "CHOLMOD Time to finalize singlethreaded com "
					     << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) << endl;
				}
				CHOLMOD(start)(&com);
			}
			com.nmethods = 1;                       // nOrderingMethods; //1;
			com.method[0].ordering = CHOLMOD_METIS; // orderingMethod; //CHOLMOD_METIS;
			factorExists = false;
		}
#endif

		if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);

		RTriangulation& Tri = T[currentTes].Triangulation();
		int             n_cells = Tri.number_of_finite_cells();
		vector<int>     clen;
		vector<int>     is;
		vector<int>     js;
		vector<Real>    vs;
		if (!areCellsOrdered) {
			T_nnz = 0;
			ncols = 0;
			///Ordered cells
			orderedCells.clear();
			const FiniteCellsIterator cellEnd = Tri.finite_cells_end();
			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
				orderedCells.push_back(cell);
				cell->info().index = 0;
				if (!cell->info().Pcondition && !cell->info().blocked) ++ncols;
			}
			orderedCells.shrink_to_fit();
			//		//Segfault on 14.10, and useless overall since SuiteSparse has preconditionners (including metis)
			// 		spatial_sort(orderedCells.begin(),orderedCells.end(), CellTraits_for_spatial_sort<RTriangulation>());
			T_cells.clear();
			T_index = 0;
			isLinearSystemSet = false;
			areCellsOrdered = true;
		}
		if (!isLinearSystemSet) {
#ifdef TAUCS_LIB
			if (Fccs) taucs_ccs_free(Fccs); //delete the old factor
#endif
			int n = 3 * (ncols + 1); //number of non-zero in triangular matrix
			is.resize(n);
			js.resize(n);
			vs.resize(n);
			T_x.resize(ncols);
			T_x.shrink_to_fit();
			T_b.resize(ncols);
			T_b.shrink_to_fit();
			T_bv.resize(ncols);
			T_bv.shrink_to_fit();
			bodv.resize(ncols);
			xodv.resize(ncols);
			//gsB.resize(ncols+1);
			T_cells.resize(ncols + 1);
			T_cells.shrink_to_fit();
			T_nnz = 0;
		}
		for (int kk = 0; kk < ncols; kk++)
			T_b[kk] = 0;
		///Ordered cells
		int        index = 0, nIndex = 0;
		CellHandle neighbourCell;
		for (int i2 = 0; i2 < n_cells; i2++) {
			FiniteCellsIterator& cell = orderedCells[i2];
			///Non-ordered cells
			// 	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cell_end; cell++) {
			if (!cell->info().Pcondition && !cell->info().blocked) {
				index = cell->info().index;
				if (index == 0) {
					T_cells[++T_index] = cell;
					cell->info().index = index = T_index;
				}
				if (!isLinearSystemSet) {
					//Add diagonal term
					is[T_nnz] = index;
					js[T_nnz] = index;
					vs[T_nnz] = 0;
					for (int j = 0; j < 4; j++)
						if (!cell->neighbor(j)->info().blocked) vs[T_nnz] += (cell->info().kNorm())[j];
					// 				vs[T_nnz] = (cell->info().kNorm())[0]+ (cell->info().kNorm())[1]+ (cell->info().kNorm())[2]+ (cell->info().kNorm())[3];
					if (fluidBulkModulus > 0) {
						if (cell->info().isCavity && phiZero > 0)
							vs[T_nnz] += (1.f * equivalentCompressibility / (dt * cell->info().invVoidVolume()));
						else
							vs[T_nnz] += (1.f / (dt * fluidBulkModulus * cell->info().invVoidVolume()));
					}
					++T_nnz;
				}
				for (int j = 0; j < 4; j++) {
					neighbourCell = cell->neighbor(j);
					nIndex = neighbourCell->info().index;
					if (Tri.is_infinite(neighbourCell)) continue;
					if (!isLinearSystemSet && !(neighbourCell->info().Pcondition || neighbourCell->info().blocked)) {
						if (nIndex == 0) {
							T_cells[++T_index] = neighbourCell;
							neighbourCell->info().index = nIndex = T_index;
						} else if (index > nIndex) {
							is[T_nnz] = index;
							js[T_nnz] = nIndex;
							vs[T_nnz] = -(cell->info().kNorm())[j];
							T_nnz++;
						}
					} else if (neighbourCell->info().Pcondition && !neighbourCell->info().blocked) {
						//ADD TO b, FIXME : use updated volume change
						T_b[index - 1] += cell->info().kNorm()[j] * neighbourCell->info().p();
					}
				}
			}
		}

		if (getCHOLMODPerfTimings) {
			gettimeofday(&end, NULL);
			cout << "CHOLMOD Time to build linear equations " << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec))
			     << endl;
		}
		updatedRHS = true;
		if (!isLinearSystemSet) {
			if (useSolver == 1 || useSolver == 2) {
#ifdef TAUCS_LIB
				clen.resize(ncols + 1);
				T_jn.resize(ncols + 1);
				T_A->colptr = &T_jn[0];
				T_ia.resize(T_nnz);
				T_A->rowind = &T_ia[0];
				T_A->flags = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
				T_an.resize(T_nnz);
				T_A->values.d = &T_an[0];
				T_A->n = ncols;
				T_A->m = ncols;
				int i, j, k;
				for (j = 0; j < ncols; j++)
					clen[j] = 0;
				for (k = 0; k < T_nnz; k++) {
					i = is[k] - 1; /* make it 1-based */
					j = js[k] - 1; /* make it 1-based */
					(clen[j])++;
				}
				/* now compute column pointers */
				k = 0;
				for (j = 0; j < ncols; j++) {
					int tmp;
					tmp = clen[j];
					clen[j] = (T_A->colptr[j]) = k;
					k += tmp;
				}
				clen[ncols] = (T_A->colptr[ncols]) = k;

				/* now read matrix into data structure */
				for (k = 0; k < T_nnz; k++) {
					i = is[k] - 1; /* make it 1-based */
					j = js[k] - 1; /* make it 1-based */
					assert(i < ncols);
					assert(j < ncols);
					(T_A->taucs_values)[clen[j]] = vs[k];
					(T_A->rowind)[clen[j]] = i;
					clen[j]++;
					// 			cerr<<"i="<< i <<" j="<< j<<" v="<<vs[k]<<" clen[j]="<<clen[j]-1<<endl;
				}
#endif //TAUCS_LIB
#ifdef LINSOLV
			} else if (useSolver == 3) {
				tripletList.clear();
				tripletList.resize(T_nnz);
				tripletList.shrink_to_fit();
				for (int k = 0; k < T_nnz; k++)
					tripletList[k] = ETriplet(is[k] - 1, js[k] - 1, vs[k]);
				A.resize(ncols, ncols);
				A.data().squeeze();
				A.setFromTriplets(tripletList.begin(), tripletList.end());
#endif
#ifdef SUITESPARSE_VERSION_4
			} else if (useSolver == 4) {
				if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);
				cholmod_triplet* trip = CHOLMOD(allocate_triplet)(ncols, ncols, T_nnz, 1, CHOLMOD_REAL, &com);
				// set all the values for the cholmod triplet matrix
				for (int k = 0; k < T_nnz; k++) {
					add_T_entry(trip, is[k] - 1, js[k] - 1, vs[k]);
				}
				Achol = CHOLMOD(triplet_to_sparse)(trip, trip->nnz, &com);
				CHOLMOD(free_triplet)(&trip, &com);
				//trip=0;
				if (getCHOLMODPerfTimings) {
					CHOLMOD(print_sparse)(Achol, "Achol", &com);
					gettimeofday(&end, NULL);
					cout << "CHOLMOD Time to allocate matrix "
					     << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) << endl;
				}
#endif
			}
			isLinearSystemSet = true;
		}
		return ncols;
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::copyGsToCells()
	{
		for (int ii = 1; ii <= ncols; ii++)
			T_cells[ii]->info().p() = gsP[ii];
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::copyCellsToGs(Real dt)
	{
		for (int ii = 1; ii <= ncols; ii++) {
			gsP[ii] = T_cells[ii]->info().p();
			gsdV[ii] = T_cells[ii]->info().dv();
			if (fluidBulkModulus > 0) { gsdV[ii] -= T_cells[ii]->info().p() / (fluidBulkModulus * dt * T_cells[ii]->info().invVoidVolume()); }
		}
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::copyLinToCells()
	{
		for (int ii = 1; ii <= ncols; ii++) {
			T_cells[ii]->info().p() = T_x[ii - 1];
		}
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::copyCellsToLin(Real dt)
	{
		// trunk
		/*	for (int ii=1; ii<=ncols; ii++) {*/
		/*		T_bv[ii-1]=T_b[ii-1]-T_cells[ii]->info().dv();*/
		/*		if (fluidBulkModulus>0) T_bv[ii-1] += T_cells[ii]->info().p()/(fluidBulkModulus*dt*T_cells[ii]->info().invVoidVolume());*/
		/*	}*/


// update T_b for cells abutting Pconditions since our Pcondition value might want to change sub remesh
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (int ii = 1; ii <= ncols; ii++) {
			if (controlCavityPressure) {
				T_b[ii - 1] = 0;
				for (int j = 0; j < 4; j++) {
					CellHandle neighborCell = T_cells[ii]->neighbor(j);
					if (neighborCell->info().Pcondition && !neighborCell->info().blocked) {
						T_b[ii - 1] += T_cells[ii]->info().kNorm()[j] * (neighborCell->info().p());
					}
				}
			}
			T_bv[ii - 1] = T_b[ii - 1] - T_cells[ii]->info().dv();
			if (fluidBulkModulus > 0) {
				if (phiZero > 0 && T_cells[ii]->info().isCavity) { // consider air compressibility in cavity
					T_bv[ii - 1] += T_cells[ii]->info().p() * equivalentCompressibility / (dt * T_cells[ii]->info().invVoidVolume());
					if (controlCavityVolumeChange) T_bv[ii - 1] += cavityDV;
				} else { // use normal bulkmodulus
					T_bv[ii - 1] += T_cells[ii]->info().p() / (fluidBulkModulus * dt * T_cells[ii]->info().invVoidVolume());
				}
			}
		}
	}

	/// For Gauss Seidel, we need the full matrix
	template <class _Tesselation, class FlowType>
	// int FlowBoundingSphereLinSolv<_Tesselation,FlowType>::SetLinearSystemFullGS()
	int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::setLinearSystemFullGS(Real dt)
	{
		//WARNING : boundary conditions (Pcondition, p values) must have been set for a correct definition
		RTriangulation& Tri = T[currentTes].Triangulation();
		int             n_cells = Tri.number_of_finite_cells();
		if (!areCellsOrdered) {
			T_cells.clear();
			T_index = 0;
			T_nnz = 0;
			ncols = 0;
			const FiniteCellsIterator cellEnd = Tri.finite_cells_end();
			orderedCells.clear();

			for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
				orderedCells.push_back(cell);
				if (!cell->info().Pcondition && !cell->info().blocked) ++ncols;
			}
			//FIXME: does it really help? test by commenting this "sorting" line
			spatial_sort(orderedCells.begin(), orderedCells.end(), CellTraits_for_spatial_sort<RTriangulation>());

			// 		Real pZero=0;
			// 		if (yMinId>=0 and yMaxId>yMinId) pZero = abs((boundary(yMinId).value-boundary(yMaxId).value)/2);
			gsP.resize(ncols + 1);
			// 		_gsP.resize(ncols+1);
			gsB.resize(ncols + 1);
			T_b.resize(ncols + 1);
			gsdV.resize(ncols + 1);
			fullAcolumns.resize(ncols + 1);
			fullAvalues.resize(ncols + 1);
			T_cells.resize(ncols + 1);
			for (int k = 0; k <= ncols; k++) {
				fullAcolumns[k].resize(4);
				fullAvalues[k].resize(5);
				gsdV[k] = 0;
				// 			gsP[k]=pZero;
			}
			// 		_gsP[0]= &ZERO;
			gsP[0] = 0;
			areCellsOrdered = true;
			isFullLinearSystemGSSet = false;
		}
		for (int k = 0; k <= ncols; k++)
			gsB[k] = 0;

		///we build the full matrix + RHS here, else only the RHS in the other loop
		if (!isFullLinearSystemGSSet)
			///Ordered cells
			for (int i2 = 0; i2 < n_cells; i2++) {
				FiniteCellsIterator& cell = orderedCells[i2];
				///Non-ordered cells
				// 	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cell_end; cell++) {
				if (!cell->info().Pcondition && !cell->info().blocked) {
					if (cell->info().index == 0) {
						T_cells[++T_index] = cell;
						cell->info().index = T_index;
					}
					gsP[cell->info().index] = cell->info().pression;
					//Add diagonal term
					Real num
					        = (cell->info().kNorm())[0] + (cell->info().kNorm())[1] + (cell->info().kNorm())[2] + (cell->info().kNorm())[3];
					if (fluidBulkModulus > 0) num += (1.f / (dt * fluidBulkModulus * cell->info().invVoidVolume()));
					fullAvalues[cell->info().index][4] = 1.f / num;
					++T_nnz;

					for (int j = 0; j < 4; j++) {
						CellHandle neighbourCell = cell->neighbor(j);
						if (Tri.is_infinite(neighbourCell)) {
							fullAvalues[cell->info().index][j] = 0;
							fullAcolumns[cell->info().index][j] = &gsP[0];
							continue;
						}
						if (!neighbourCell->info().Pcondition) {
							if (neighbourCell->info().index == 0) {
								T_cells[++T_index] = neighbourCell;
								neighbourCell->info().index = T_index;
							}
							++T_nnz;
							fullAvalues[cell->info().index][j] = (cell->info().kNorm())[j];
							fullAcolumns[cell->info().index][j] = &gsP[neighbourCell->info().index];
						} else {
							fullAvalues[cell->info().index][j] = 0;
							fullAcolumns[cell->info().index][j] = &gsP[0];
							gsB[cell->info().index] += cell->info().kNorm()[j] * neighbourCell->info().p();
						}
					}
				}
			}
		///define only the new RHS, accouting for new imposed pressures
		else
			///Ordered cells
			for (int i2 = 0; i2 < n_cells; i2++) {
				FiniteCellsIterator& cell = orderedCells[i2];
				///Non-ordered cells
				// 	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cell_end; cell++) {
				if (!cell->info().Pcondition && !cell->info().blocked)
					for (int j = 0; j < 4; j++) {
						CellHandle neighbourCell = cell->neighbor(j);
						if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().Pcondition)
							gsB[cell->info().index] += cell->info().kNorm()[j] * neighbourCell->info().p();
					}
			}
		isFullLinearSystemGSSet = true;
		return ncols;
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::vectorizedGaussSeidel(Real dt)
	{
		using math::max;
		using math::min;

		// 	cout<<"VectorizedGaussSeidel"<<endl;
		if (!isFullLinearSystemGSSet || (isFullLinearSystemGSSet && reApplyBoundaryConditions())) setLinearSystemFullGS(dt);
		copyCellsToGs(dt);

		int  j = 0;
		Real dp_max, p_max, sum_p, p_moy, dp_moy, sum_dp;

#ifdef GS_OPEN_MP
		const int num_threads = 1;
		omp_set_num_threads(num_threads);
		vector<Real> t_sum_p, t_dp_max, t_sum_dp, t_p_max;
		t_sum_dp.resize(num_threads);
		t_dp_max.resize(num_threads);
		t_p_max.resize(num_threads);
		t_sum_p.resize(num_threads);
#endif
		int j2 = -1;
		dp_max = 0;
		p_max = 0;
		p_moy = 0;
		dp_moy = 0;
		sum_p = 0;
		sum_dp = 0;
		do {
			if (++j2 >= 10) j2 = 0; //compute max/mean only each 10 iterations
			if (j2 == 0) {
				dp_max = 0;
				p_max = 0;
				p_moy = 0;
				dp_moy = 0;
				sum_p = 0;
				sum_dp = 0;
#ifdef GS_OPEN_MP
				for (int ii = 0; ii < num_threads; ii++)
					t_p_max[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_dp_max[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_sum_p[ii] = 0;
				for (int ii = 0; ii < num_threads; ii++)
					t_sum_dp[ii] = 0;
#endif
			}
#ifdef GS_OPEN_MP
#pragma omp parallel for schedule(dynamic, 2000)
#endif

			for (int ii = 1; ii <= ncols; ii++) {
				double** Acols = &(fullAcolumns[ii][0]);
				double*  Avals = &(fullAvalues[ii][0]);
				double   dp = (((gsB[ii] - gsdV[ii] + Avals[0] * (*Acols[0]) + Avals[1] * (*Acols[1]) + Avals[2] * (*Acols[2])
                                               + Avals[3] * (*Acols[3]))
                                              * Avals[4])
                                             - gsP[ii])
				        * relax;

				gsP[ii] = dp + gsP[ii];
				if (j2 == 0) {
#ifdef GS_OPEN_MP
					const int tn = omp_get_thread_num();
					t_sum_dp[tn] += math::abs(dp);
					t_dp_max[tn] = max(t_dp_max[tn], math::abs(dp));
					t_p_max[tn] = max(t_p_max[tn], gsP[ii]);
					t_sum_p[tn] += math::abs(gsP[ii]);
#else
					dp_max = max(dp_max, math::abs(dp));
					p_max = max(p_max, math::abs(gsP[ii]));
					sum_p += math::abs(gsP[ii]);
					sum_dp += math::abs(dp);
#endif
				}
			}
#ifdef GS_OPEN_MP
			if (j2 == 0) {
				for (int jj = 0; jj < num_threads; jj++)
					p_max = max(p_max, t_p_max[jj]);
				for (int jj = 0; jj < num_threads; jj++)
					dp_max = max(dp_max, t_dp_max[jj]);
				for (int jj = 0; jj < num_threads; jj++)
					sum_p += t_sum_p[jj];
				for (int jj = 0; jj < num_threads; jj++)
					sum_dp += t_sum_dp[jj];
			}
#endif
			if (j2 == 0) {
				p_moy = sum_p / ncols;
				dp_moy = sum_dp / ncols;
				if (debugOut) cerr << "GS : j=" << j << " p_moy=" << p_moy << " dp_moy=" << dp_moy << endl;
			}
#ifdef GS_OPEN_MP
#pragma omp master
#endif
			j++;
		} while ((dp_max / p_max) > tolerance && j < 20000 /*&& ( dp_max > tolerance )*/ /* &&*/ /*( j<50 )*/);
		copyGsToCells();
		if (j >= 20000) cerr << "GS did not converge in 20k iterations (maybe because the reference pressure is 0?)" << endl;
		if (debugOut) cerr << "GS iterations : " << j - 1 << endl;
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::sortV(int k1, int k2, int* is, Real* ds)
	{
		for (int k = k1; k < k2; k++) {
			int kk = k;
			while (kk >= k1 && is[kk] > is[kk + 1]) {
				swapFwd(is, kk);
				swapFwd(ds, kk);
				--kk;
			}
		}
	}

	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::eigenSolve(Real dt)
	{
#ifdef LINSOLV
		if (!isLinearSystemSet || (isLinearSystemSet && reApplyBoundaryConditions()) || !updatedRHS) ncols = setLinearSystem(dt);
		copyCellsToLin(dt);
		//FIXME: we introduce new Eigen vectors, then we have to copy from/to c-arrays, can be optimized later
		VectorXr eb(ncols);
		VectorXr ex(ncols);
		for (int k = 0; k < ncols; k++)
			eb[k] = T_bv[k];
		if (!factorizedEigenSolver) {
			eSolver.setMode(Eigen::CholmodSupernodalLLt);
			openblas_set_num_threads(numFactorizeThreads);
			eSolver.compute(A);
			//Check result
			if (eSolver.cholmod().status > 0) {
				cerr << "something went wrong in Cholesky factorization, use LDLt as fallback this time" << eSolver.cholmod().status << endl;
				eSolver.setMode(Eigen::CholmodLDLt);
				eSolver.compute(A);
			}
			factorizedEigenSolver = true;
		}
		// backgroundAction only wants to factorize, no need to solve and copy to cells.
		if (!factorizeOnly) {
			openblas_set_num_threads(numSolveThreads);
			ex = eSolver.solve(eb);
			for (int k = 0; k < ncols; k++)
				T_x[k] = ex[k];
			copyLinToCells();
		}
#else
		cerr << "Flow engine not compiled with eigen, nothing computed if useSolver=3" << endl;
#endif
		return 0;
	}

	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::cholmodSolve(Real dt)
	{
#ifdef SUITESPARSE_VERSION_4
		if (!isLinearSystemSet || (isLinearSystemSet && reApplyBoundaryConditions()) || !updatedRHS) ncols = setLinearSystem(dt);
		copyCellsToLin(dt);
		cholmod_dense* B = CHOLMOD(zeros)(ncols, 1, Achol->xtype, &com); //cholmod_l_zeros(ncols, 1, Achol->xtype, &com);
		Real*          B_x = (Real*)B->x;
		for (int k = 0; k < ncols; k++)
			B_x[k] = T_bv[k];
		if (!factorizedEigenSolver) {
			openblas_set_num_threads(numFactorizeThreads);
			if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);
			if (!reuseOrdering) {
				L = CHOLMOD(analyze)(Achol, &com); //cholmod_l_analyze(Achol, &com);
			} else {
				N = CHOLMOD(copy_factor)(L, &com);
			}
			if (getCHOLMODPerfTimings) {
				gettimeofday(&end, NULL);
				cout << "Reusing reordering? " << reuseOrdering << ". CHOLMOD Time to Analyze "
				     << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) << endl;
			}
			if (getCHOLMODPerfTimings) gettimeofday(&start, NULL);
			if (!reuseOrdering) {
				CHOLMOD(factorize)(Achol, L, &com); //cholmod_l_factorize(Achol, L, &com);
			} else {
				CHOLMOD(factorize)(Achol, N, &com);
			}

			if (getCHOLMODPerfTimings) {
				gettimeofday(&end, NULL);
				cout << "CHOLMOD Time to factorize " << ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec))
				     << endl;
			}
			factorExists = true;
			factorizedEigenSolver = true;
		}

		if (!factorizeOnly) {
			openblas_set_num_threads(numSolveThreads);
			cholmod_dense* ex = CHOLMOD(solve)(CHOLMOD_A, L, B, &com); // cholmod_l_solve(CHOLMOD_A, L, B, &com);
			Real*          e_x = (Real*)ex->x;
			for (int k = 0; k < ncols; k++) {
				T_x[k] = e_x[k];
			}
			copyLinToCells();
			CHOLMOD(free_dense)(&ex, &com); //cholmod_l_free_dense(&ex, &com);
		}
		CHOLMOD(free_dense)(&B, &com); //cholmod_l_free_dense(&B, &com);
#else
		cerr << "Flow engine not compiled with CHOLMOD, nothing computed if useSolver=4" << endl;
#endif
		return 0;
	}

	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::initializeInternalEnergy()
	{
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i2 = 0; i2 < sizeCells; i2++) {
			CellHandle& cell = Tes.cellHandles[i2];
			if (!cell->info().isFictious && !cell->info().blocked && !cell->info().isCavity) {
				Real volume = thermalPorosity > 0 ? thermalPorosity / cell->info().invVoidVolume() : 1. / cell->info().invVoidVolume();
				cell->info().internalEnergy = fluidCp * fluidRho * cell->info().temp() * volume;
			}
			if (cell->info().isCavity)
				cell->info().internalEnergy = fluidCp * fluidRho * cell->info().temp()
				        * (cell->info().volume()); // ignore particles used for fluid discr. in cavity (i.e. use volume())
		}
	}


	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::augmentConductivityMatrix(Real dt)
	{
		Real            energyFlux;
		Real            upwindTemp;
		Real            facetFlowRate;
		RTriangulation& Tri = T[currentTes].Triangulation();
		// cycle through facets instead of cells to avoid duplicate math
		/*   	Tesselation& Tes = T[currentTes];*/
		/*	const long sizeFacets = Tes.facetHandles.size();*/
		/*	#pragma omp parallel for*/
		/*    	for (long i=0; i<sizeFacets; i++){*/
		/*		Facet& facet = Tes.facetHandles[i];*/
		/*		const CellHandle& cell = facet->first;*/
		/*		const CellHandle& neighborCell = facet->first->neighbor(facet->second);*/
		/*		facetFlowRate = cell->info().kNorm()[facet->second] * (cell->info().p() - cell->neighbor(facet->second)->info().p());*/
		/*		if (facetFlowRate>0){*/
		/*			upwindTemp = cell->info().temp();*/
		/*		} else { */
		/*			upwindTemp = neighborCell->info().temp();*/
		/*		}*/
		/*		energyFlux = fluidCp*fluidRho*dt*upwindTemp*facetFlowRate;*/
		/*		if (!cell->info().Tcondition && !cell->info().isFictious && !cell->info().blocked) cell->info().internalEnergy -= energyFlux;*/
		/*		if (!neighborCell->info().Tcondition && !neighborCell->info().isFictious && !neighborCell->info().blocked) neighborCell->info().internalEnergy += energyFlux;*/
		/*	}*/


		for (FiniteFacetsIterator f_it = Tri.finite_facets_begin(); f_it != Tri.finite_facets_end(); f_it++) {
			const CellHandle& cell = f_it->first;
			const CellHandle& neighborCell = f_it->first->neighbor(f_it->second);
			if (cell->info().blocked || neighborCell->info().blocked || (cell->info().Pcondition && neighborCell->info().Pcondition)) continue;
			facetFlowRate = cell->info().kNorm()[f_it->second] * (cell->info().p() - cell->neighbor(f_it->second)->info().p());
			if (facetFlowRate > 0) {
				upwindTemp = cell->info().temp();
			} else {
				upwindTemp = neighborCell->info().temp();
			}
			energyFlux = fluidCp * fluidRho * dt * upwindTemp * facetFlowRate;
			if (!cell->info().Tcondition && !cell->info().isFictious) cell->info().internalEnergy -= energyFlux;
			if (!neighborCell->info().Tcondition && !neighborCell->info().isFictious) neighborCell->info().internalEnergy += energyFlux;
		}
	}


	template <class _Tesselation, class FlowType> void FlowBoundingSphereLinSolv<_Tesselation, FlowType>::setNewCellTemps(bool addToDeltaTemp)
	{
		Tesselation& Tes = T[currentTes];
		const long   sizeCells = Tes.cellHandles.size();
		Real         cavityInternalEnergy = 0;
		Real         cavityVolume = 0;
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
		for (long i2 = 0; i2 < sizeCells; i2++) {
			CellHandle& cell = Tes.cellHandles[i2];
			if (!cell->info().isFictious && !cell->info().blocked) {
				Real oldTemp = cell->info().temp();
				//cell->info().temp()=cell->info().internalEnergy/(cell->info().volume()*fluidCp*fluidRho);
				if (!cell->info().isCavity) {
					Real volume = thermalPorosity > 0 ? thermalPorosity / cell->info().invVoidVolume() : 1. / cell->info().invVoidVolume();
					cell->info().temp() = cell->info().internalEnergy
					        / (volume * fluidCp
					           * fluidRho); //FIXME: invVoidVolume depends on volumeSolidPore() which uses CGAL points only updated each remesh. We might need our own volumeSolidPore().
				} else {
					cell->info().temp() = cell->info().internalEnergy / ((cell->info().volume()) * fluidCp * fluidRho);
				}
				if (!addToDeltaTemp) cell->info().dtemp() = cell->info().temp() - oldTemp;
				else
					cell->info().dtemp() += cell->info().temp()
					        - oldTemp; // fluid conduction is a midstep proces that uses a midsteptemp, in this case we want to add to the existing deltatemp so that fluid expansion is computed based on full temp step
			}
			if (controlCavityPressure && cell->info().isCavity && !cell->info().blocked) {
				cavityInternalEnergy += cell->info().internalEnergy;
				cavityVolume += 1. / cell->info().invVoidVolume();
			}
		}
		Real cavityTemp;
		if (controlCavityPressure) {
			cavityTemp = cavityInternalEnergy / (cavityVolume * fluidCp * fluidRho); //use cavityFluidDensity?
#ifdef YADE_OPENMP
#pragma omp parallel for
#endif
			for (long i2 = 0; i2 < sizeCells; i2++) {
				CellHandle& cell = Tes.cellHandles[i2];
				if (!cell->info().isCavity) continue;
				Real oldTemp = cell->info().temp();
				cell->info().temp() = cavityTemp;
				cell->info().dtemp() = cell->info().temp() - oldTemp;
			}
		}
	}

	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::taucsSolve(Real /*dt*/)
	{
#ifdef TAUCS_LIB
		if (debugOut) cerr << endl << "TAUCS solve" << endl;
		Real t = taucs_ctime();  //timer
		Real t2 = taucs_ctime(); //global timer
		if (!isLinearSystemSet || (isLinearSystemSet && reApplyBoundaryConditions())) {
			ncols = setLinearSystem(dt);
			if (debugOut) cerr << "Assembling the matrix : " << taucs_ctime() - t << endl;
			t = taucs_ctime();
		}

		copyCellsToLin(dt);
		if (debugOut) cerr << "Updating dv's (Yade->LinSolver) : " << taucs_ctime() - t << endl;
		t = taucs_ctime();
		//taucs_logfile("stdout");//! VERY USEFULL!!!!!!!!!!! (disable to exclude output time from taucs_ctime() measurments)

		taucs_double* x = &T_x[0]; // the unknown vector to solve Ax=b
		taucs_double* bod = &bodv[0];
		taucs_double* xod = &xodv[0];

		if (Fccs == NULL) {
			if (debugOut) cerr << "_entering taucs_" << endl;
			// 1) Reordering
			taucs_ccs_order(T_A, &perm, &invperm, (char*)"metis");
			if (debugOut) cerr << "_entering taucs_" << endl;
			taucs_ccs_matrix* Aod;
			if (debugOut) cerr << "_entering taucs_" << endl;
			Aod = taucs_ccs_permute_symmetrically(T_A, perm, invperm);
			if (debugOut) cerr << "Reordering : " << taucs_ctime() - t << endl;
			t = taucs_ctime();

			// 2) Factoring
			F = taucs_ccs_factor_llt_mf(Aod);
			if (F == NULL) cerr << "factorization failed" << endl;
			taucs_dccs_free(Aod);
			Aod = NULL;
			//convert F to ccs for faster solve
			Fccs = taucs_supernodal_factor_to_ccs(F);
			//... then delete F
			taucs_supernodal_factor_free(F);
			F = NULL;
			if (debugOut) cerr << "Factoring : " << taucs_ctime() - t << endl;
			t = taucs_ctime();
		}
		taucs_vec_permute(ncols, TAUCS_DOUBLE, &T_bv[0], bod, perm);
		// 3) Back substitution and reodering the solution back
		taucs_ccs_solve_llt(Fccs, xod, bod); //the ccs format (faster)
		                                     // 	taucs_supernodal_solve_llt(F, xod, bod);//the blackbox format (slower)
		if (debugOut) cerr << "Solving : " << taucs_ctime() - t << endl;
		t = taucs_ctime();
		t = taucs_ctime();
		taucs_vec_ipermute(ncols, TAUCS_DOUBLE, xod, x, perm);
		//     	cerr << "Deordering : " <<  taucs_ctime()-t << endl; t = taucs_ctime();
		// 4) Copy back to the triangulation
		copyLinToCells();
		// 	cerr << "Updating P (LinSolver->Yade) : " <<  taucs_ctime()-t << endl;
		if (debugOut) cerr << "Total TAUCS time ................ : " << taucs_ctime() - t2 << endl;
#else
		cerr << "Flow engine not compiled with taucs, nothing computed if useSolver=1" << endl;
#endif
		return 0;
	}
	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::pardisoSolve(Real /*dt*/)
	{
		cerr << endl << "PardisoSolve solve" << endl;
#ifndef PARDISO
		return 0;
#else
		Real iniT = taucs_ctime();

		if (debugOut) cerr << "_entering pardiso_" << endl;
		/* Matrix data. */
		Real t = taucs_ctime(); //timer
		bool wasLSystemSet = isLinearSystemSet;
		if (!isLinearSystemSet || (isLinearSystemSet && reApplyBoundaryConditions())) {
			ncols = setLinearSystem(dt);
			if (debugOut) cerr << "Assembling the matrix : " << taucs_ctime() - t << endl;
			t = taucs_ctime();
		}

		if (debugOut) cerr << taucs_ctime() - t << "s : set system" << endl;
		t = taucs_ctime();
		ia = T_A->colptr;
		ja = T_A->rowind;
		a = T_A->values.d;
		if (debugOut) cerr << taucs_ctime() - t << "s : set system" << endl;
		if (!wasLSystemSet)
			for (int k = 0; k < ncols; k++)
				sortV(ia[k], ia[k + 1] - 1, ja, a);
		if (debugOut) cout << taucs_ctime() - t << "s for ordering CCS format" << endl;
		t = taucs_ctime();

		nnz = ia[ncols];
		//    int mtype = -2;        /* Real symmetric matrix */
		mtype = 2; /* Real symmetric positive def. matrix */
		/* RHS and solution vectors. */
		copyCellsToLin(dt);
		b = &T_bv[0];
		// 	P_x.resize(n);
		x = &T_x[0]; // the unknown vector to solve Ax=b
		nrhs = 1;    /* Number of right hand sides. */

		error = 0;
		solver = 0; /* use sparse direct solver */

		/* Numbers of processors, value of OMP_NUM_THREADS */

		if (!pardisoInitialized) {
			var = getenv("OMP_NUM_THREADS");
			if (var != NULL) sscanf(var, "%d", &num_procs);
			else {
				num_procs = 1;
				cerr << "Set environment OMP_NUM_THREADS to something. Pardiso needs it defined!" << endl;
			}
			if (debugOut) cerr << taucs_ctime() - t << "pardisoinit" << endl;
			F77_FUNC(pardisoinit)(pt, &mtype, &solver, iparm, dparm, &error);
			if (debugOut) cerr << taucs_ctime() - t << "pardisoinit'ed" << endl;
			pardisoInitialized = true;
			if (error != 0) {
				if (error == -10) printf("No license file found \n");
				if (error == -11) printf("License is expired \n");
				if (error == -12) printf("Wrong username or hostname \n");
				return 1;
			}
			iparm[2] = num_procs;
			maxfct = 1; /* Maximum number of numerical factorizations.  */
			mnum = 1;   /* Which factorization to use. */
			msglvl = 0; /* Print statistical information  */
			error = 0;  /* Initialize error flag */

			/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
			/*     notation.                                                        */
			if (debugOut) cout << taucs_ctime() - t << "tuning" << endl;
			t = taucs_ctime();
			for (i = 0; i < ncols + 1; i++) {
				ia[i] += 1;
			}
			for (i = 0; i < nnz; i++) {
				ja[i] += 1;
			}
			if (debugOut) cout << taucs_ctime() - t << "s : Convert matrix from 0-based" << endl;
			t = taucs_ctime();
			/* ..  Reordering and Symbolic Factorization.  This step also allocates */
			/*     all memory that is necessary for the factorization.              */
			phase = 11;
			F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &ncols, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
			if (error != 0) {
				printf("\nERROR during symbolic factorization: %d", error);
				exit(1);
			}
			if (debugOut) cout << taucs_ctime() - t << "s : Reordering and Symbolic Factorization" << endl;
			t = taucs_ctime();

			/* ..  Numerical factorization.                                         */
			phase = 22;
			F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &ncols, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);

			if (error != 0) {
				printf("\nERROR during numerical factorization: %d", error);
				exit(2);
			}
			if (debugOut) cerr << taucs_ctime() - t << "s : Numerical factorization. " << endl;
			t = taucs_ctime();
		}
		/* ..  Back substitution and iterative refinement.                      */
		phase = 33;
		iparm[7] = 0; /* Max numbers of iterative refinement steps. */
		t = taucs_ctime();
		F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &ncols, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error, dparm);
		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			exit(3);
		}
		if (debugOut) cerr << taucs_ctime() - t << "s : Back substitution and iterative refinement." << endl;
		t = taucs_ctime();
		copyLinToCells();
		if (debugOut) cerr << taucs_ctime() - t << "s : Copy back." << endl;

		if (wasLSystemSet) {
			pTime1N++;
			pTime1 += (taucs_ctime() - iniT);
		} else {
			pTime2N++;
			pTime2 += (taucs_ctime() - iniT);
		}
		if (pTimeInt > 99) {
			cout << "Pardiso.....  " << pTime1 / (double)pTime1N << " s/iter for " << pTime1N << "/" << (pTime2N + pTime1N) << " std iter., "
			     << pTime2 / pTime2N << " s/iter for " << pTime2N << "/" << (pTime2N + pTime1N) << " retriangulation iter." << endl;
			pTime1 = 0;
			pTime2 = 0;
			pTime1N = 0;
			pTime2N = 0;
			pTimeInt = 0;
		}
		pTimeInt++;
		/* ..  Convert matrix back to 0-based C-notation.                       */
		// 	for (i = 0; i < n+1; i++) {
		// 		ia[i] -= 1;
		// 	}
		// 	for (i = 0; i < nnz; i++) {
		// 		ja[i] -= 1;
		// 	}
		return 0;
#endif
	}


	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::pardisoSolveTest()
	{
#ifndef PARDISO
		return 0;
#else
		/* Matrix data. */
		Real t = taucs_ctime(); //timer
		bool wasLSystemSet = isLinearSystemSet;
		int n = setLinearSystem();
		// 	ncols=n;//for VectorizesGS
		cout << taucs_ctime() - t << "s : set system" << endl;
		t = taucs_ctime();
		int* ia = T_A->colptr;
		int* ja = T_A->rowind;
		double* a = T_A->values.d;

		if (!wasLSystemSet)
			for (int k = 0; k < n; k++)
				sortV(ia[k], ia[k + 1] - 1, ja, a);
		cout << taucs_ctime() - t << "s for ordering CCS format" << endl;
		t = taucs_ctime();

		int nnz = ia[n];
		//    int mtype = -2;        /* Real symmetric matrix */
		int mtype = 2; /* Real symmetric positive def. matrix */
		/* RHS and solution vectors. */
		double* b = &T_b[0];
		P_x.resize(n);
		double* x = &P_x[0]; // the unknown vector to solve Ax=b
		int nrhs = 1;        /* Number of right hand sides. */

		/* Internal solver memory pointer pt,                  */
		/* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
		/* or void *pt[64] should be OK on both architectures  */
		void* pt[64];
		/* Pardiso control parameters. */
		int iparm[64];
		double dparm[64];
		int maxfct, mnum, phase, error, msglvl, solver;
		/* Number of processors. */
		int num_procs;
		/* Auxiliary variables. */
		char* var;
		int i;

		double ddum; /* Double dummy */
		int idum;    /* Integer dummy. */
		/* ..  Setup Pardiso control parameters.                     */
		error = 0;
		solver = 0; /* use sparse direct solver */

		/* Numbers of processors, value of OMP_NUM_THREADS */
		var = getenv("OMP_NUM_THREADS");
		if (var != NULL) sscanf(var, "%d", &num_procs);
		else {
			num_procs = 1;
			printf("Set environment OMP_NUM_THREADS to something. Pardiso needs it defined. \n");
			// 	exit(0);
		}

		F77_FUNC(pardisoinit)(pt, &mtype, &solver, iparm, dparm, &error);
		if (error != 0) {
			if (error == -10) printf("No license file found \n");
			if (error == -11) printf("License is expired \n");
			if (error == -12) printf("Wrong username or hostname \n");
			return 1;
		}
		iparm[2] = num_procs;
		maxfct = 1; /* Maximum number of numerical factorizations.  */
		mnum = 1;   /* Which factorization to use. */
		msglvl = 0; /* Print statistical information  */
		error = 0;  /* Initialize error flag */

		/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
		/*     notation.                                                        */
		cout << taucs_ctime() - t << "tuning" << endl;
		t = taucs_ctime();
		for (i = 0; i < n + 1; i++) {
			ia[i] += 1;
		}
		for (i = 0; i < nnz; i++) {
			ja[i] += 1;
		}
		cout << taucs_ctime() - t << "s : Convert matrix from 0-based" << endl;
		t = taucs_ctime();
		/* ..  Reordering and Symbolic Factorization.  This step also allocates */
		/*     all memory that is necessary for the factorization.              */
		phase = 11;
		F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
		if (error != 0) {
			printf("\nERROR during symbolic factorization: %d", error);
			exit(1);
		}
		cout << taucs_ctime() - t << "s : Reordering and Symbolic Factorization" << endl;
		t = taucs_ctime();

		/* ..  Numerical factorization.                                         */
		phase = 22;
		F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);

		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			exit(2);
		}
		cerr << taucs_ctime() - t << "s : Numerical factorization. " << endl;
		t = taucs_ctime();

		/* ..  Back substitution and iterative refinement.                      */
		phase = 33;
		iparm[7] = 0; /* Max numbers of iterative refinement steps. */
		t = taucs_ctime();
		F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error, dparm);
		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			exit(3);
		}
		cerr << taucs_ctime() - t << "s : Back substitution and iterative refinement." << endl;
		t = taucs_ctime();

		/* ..  Convert matrix back to 0-based C-notation.                       */
		for (i = 0; i < n + 1; i++) {
			ia[i] -= 1;
		}
		for (i = 0; i < nnz; i++) {
			ja[i] -= 1;
		}
		/* ..  Termination and release of memory.                               */
		phase = -1; /* Release internal memory. */
		            // 	F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase,
		            // 	                  &n, &ddum, ia, ja, &idum, &nrhs,
		            // 	                  iparm, &msglvl, &ddum, &ddum, &error,  dparm);
		//Release only what has to be
		F77_FUNC(pardiso)(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, NULL, NULL, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
		return 0;
#endif
	}
	template <class _Tesselation, class FlowType> int FlowBoundingSphereLinSolv<_Tesselation, FlowType>::taucsSolveTest()
	{
#ifdef TAUCS_LIB
		cout << endl << "TAUCS solve test" << endl;

		double t = taucs_ctime(); //timer
		ncols = setLinearSystem();

		//taucs_logfile("stdout");//! VERY USEFULL!!!!!!!!!!! (disable to exclude output time from taucs_ctime() measurments)

		// allocate TAUCS solution vector
		T_x.resize(ncols);
		double* x = &*T_x.begin(); // the unknown vector to solve Ax=b
		cout << "Assembling the matrix2 : " << taucs_ctime() - t << endl;
		t = taucs_ctime();
		// solve the linear system
		void* F = NULL;

		//Allocate reoredered x and b
		vector<double> bodv(ncols);
		taucs_double*  bod = &*bodv.begin();
		vector<double> xodv(ncols);
		taucs_double*  xod = &*xodv.begin();


		int*              perm;
		int*              invperm;
		taucs_ccs_matrix* Aod;

		t = taucs_ctime();
		double t2 = taucs_ctime();
		// 1) Reordering
		taucs_ccs_order(T_A, &perm, &invperm, "metis");
		Aod = taucs_ccs_permute_symmetrically(T_A, perm, invperm);
		taucs_vec_permute(ncols, TAUCS_DOUBLE, &T_b[0], bod, perm);
		cout << "Reordering : " << taucs_ctime() - t << endl;
		t = taucs_ctime();
		// 2) Factoring

		F = taucs_ccs_factor_llt_mf(Aod);
		//F = taucs_ccs_factor_llt_mf(T_A);
		cout << "Factoring : " << taucs_ctime() - t << endl;
		t = taucs_ctime();
		// 3) Back substitution and reodering the solution back
		double t4 = taucs_ctime();
		// for (int k=0;k<10;k++){

		taucs_supernodal_solve_llt(F, xod, bod);
		//     cout << "B3) Solving : " <<  taucs_ctime()-t << endl;
		t = taucs_ctime();
		taucs_vec_ipermute(ncols, TAUCS_DOUBLE, xod, x, perm);
		//     cout << "B4) Deordering : " <<  taucs_ctime()-t << endl;
		t = taucs_ctime();
		// }
		double T4 = taucs_ctime() - t4;
		cout << "Solving : " << T4 << endl;
		cout << "Low level reordered total time : " << taucs_ctime() - t2 << endl;
		t2 = taucs_ctime();
		taucs_supernodal_factor_free(F);
		taucs_ccs_free(Aod);

		/// Using TAUCS inverse factoring
		// 	t = taucs_ctime(); t2 = taucs_ctime();
		// 	F=taucs_ccs_factor_xxt(Aod);
		// 	cout << "C1_) inverse factoring (reordered) : " <<  taucs_ctime()-t << endl; t = taucs_ctime();
		// 	taucs_ccs_solve_xxt (F,xod,bod);///REALLY TOO SLOW!!
		// 	cout << "C2) solve : " <<  taucs_ctime()-t << endl; t = taucs_ctime();
		// 	taucs_vec_ipermute(ncols, TAUCS_DOUBLE, xod, x, perm);
		// 	cout << "C3) reordering : " <<  taucs_ctime()-t << endl; t = taucs_ctime();
		// 	cout << "Inverse factoring total time : " <<  taucs_ctime()-t2 << endl;
		// 	taucs_linsolve(NULL, &F, 0, NULL, NULL, NULL, NULL);
		// 	taucs_ccs_free(Aod);

		// 	ofstream file("result.dat");
		//   	if (!file.is_open()) cout << "problem opening file";
		// 	for  (unsigned i = 0; i < M ; i++) {
		// 		for (unsigned j = 0; j < L ; j++) file << i << " "<< j << " "<< x[index(i,j)] << endl;
		//         	file << endl << endl;
		//         }


		//     const FiniteCellsIterator cell_end = T[currentTes].Triangulation().finite_cells_end();
		//     if (debugOut) for (FiniteCellsIterator cell = T[currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++)
		//         {
		//
		//             if (cell->info().index>0) {
		//                 cerr <<"ix "<< cell->info().index<<" ";
		//                 cerr <<cell->info().p()<<" "<<gsP[cell->info().index]<<" "<<P_x[max((unsigned int)0,cell->info().index-1)] <<" "<<x[max((unsigned int)0,cell->info().index-1)]<<" err="<<100*(cell->info().p()-x[max((unsigned int)0,cell->info().index-1)])/x[max((unsigned int)0,cell->info().index-1)]<<"%"<<endl;
		//             }
		//             else if (cell->info().index==0) {
		//                 cerr <<"ix "<< cell->info().index<<" ";
		//                 cerr <<cell->info().p()<<" "<<"NaN"<<endl;
		//             }
		//         }

#endif
		return 0;
	}


} //namespace CGT
} // namespace yade

#endif //FLOW_ENGINE
