/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>     *
*  Copyright (C) 2013 by T. Sweijen (T.sweijen@uu.nl)                    *
*  Copyright (C) 2012 by Chao Yuan <chao.yuan@3sr-grenoble.fr>           *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef TWOPHASEFLOW
#include "TwoPhaseFlowEngine.hpp"
#include <boost/range/algorithm_ext/erase.hpp>

namespace yade { // Cannot have #include directive inside.

using math::max;
using math::min; // using inside .cpp file is ok.

YADE_PLUGIN((TwoPhaseFlowEngineT)(TwoPhaseFlowEngine)(PhaseCluster));

CREATE_LOGGER(TwoPhaseFlowEngine);
CREATE_LOGGER(PhaseCluster);

PhaseCluster::~PhaseCluster()
{
#ifdef LINSOLV
	resetSolver();
#endif
}

void PhaseCluster::solvePressure()
{
	if (pores.size() == 0) {
		LOG_WARN("nothing to solve for cluster " << label);
		return;
	}
	vector<int>  clen;
	vector<int>  is;
	vector<int>  js;
	int          ncols = 0; //number of unknown
	vector<Real> vs;
	vector<Real> RHS;
	vector<Real> RHSvol;

	vector<CellHandle> pCells; //the pores in which pressure will be solved
#ifdef LINSOLV
	for (vector<CellHandle>::iterator cellIt = pores.begin(); cellIt != pores.end(); cellIt++) {
		CellHandle cell = *cellIt;
		if ((!cell->info().Pcondition) && !cell->info().blocked) {
			cell->info().index = ncols++;
			pCells.push_back(cell);
		} else
			cell->info().index = -1;
	}
	is.reserve(ncols * 3);
	js.reserve(ncols * 3);
	vs.reserve(ncols * 3);
	RHS.resize(ncols, 0);
	RHSvol.resize(ncols, 0);
	unsigned T_nnz = 0;
	for (auto cellIt = pCells.begin(); cellIt != pCells.end(); cellIt++) {
		CellHandle cell = *cellIt;
		Real       diag = 0;
		for (int j = 0; j < 4; j++)
			if ((not tes->Triangulation().is_infinite(cell->neighbor(j))) and !cell->neighbor(j)->info().blocked) {
				diag += (cell->info().kNorm())[j];
				if ((not cell->neighbor(j)->info().Pcondition) and not cell->neighbor(j)->info().isNWRes) {
					if (not factorized) {
						if (cell->info().label == cell->neighbor(j)->info().label) {
							// off-diag coeff, only if neighbor cell is part of same cluster and in upper triangular part of the matrix
							if (cell->info().index < cell->neighbor(j)->info().index) {
								T_nnz++;
								is.push_back(cell->info().index);
								js.push_back(cell->neighbor(j)->info().index);
								vs.push_back(-(cell->info().kNorm())[j]);
							}
						} else {
							LOG_WARN(
							        "adjacent pores from different W-clusters:" << cell->info().id << " "
							                                                    << cell->neighbor(j)->info().id);
						}
					}
				} else { //imposed pressure can be in the W-phase or the NW-phase but capillary pressure will be added in another loop
					// for the moment add neighbor pressure regadless of the phase
					RHS[cell->info().index] += (cell->info().kNorm())[j] * cell->neighbor(j)->info().p();
				}
			} else {
				if (tes->Triangulation().is_infinite(cell->neighbor(j))) LOG_WARN("infinite neighbour");
			}
		// define the diag coeff
		if (not factorized) {
			T_nnz++;
			is.push_back(cell->info().index);
			js.push_back(cell->info().index);
			vs.push_back(diag);
		}

		// source term from volume change, to be updated later
		RHSvol[cell->info().index] -= cell->info().dv();
	}
	for (vector<Interface>::iterator it = interfaces.begin(); it != interfaces.end(); it++) {
		if (not tes) LOG_WARN("no tes!!");
		const CellHandle& innerCell = tes->cellHandles[it->first.first];
		if (innerCell->info().Pcondition) continue;
		RHS[innerCell->info().index] -= innerCell->info().kNorm()[it->outerIndex] * it->capillaryP;
	}
	//comC.useGPU=useGPU; //useGPU;
	//FIXME: is it safe to share "comC" among parallel cluster resolution?
	if (not factorized) {
		cholmod_triplet* T = cholmod_l_allocate_triplet(ncols, ncols, T_nnz, 1, CHOLMOD_REAL, &(comC));
		for (unsigned k = 0; k < T_nnz; k++) {
			((long*)T->i)[k] = is[k];
			((long*)T->j)[k] = js[k];
			((Real*)T->x)[k] = vs[k];
		}
		T->nnz = T_nnz;
		// convert triplet list into a cholmod sparse matrix, then factorize it
		cholmod_sparse* AcholC = cholmod_l_triplet_to_sparse(T, T->nnz, &(comC));
		LC = cholmod_l_analyze(AcholC, &(comC));
		cholmod_l_factorize(AcholC, LC, &(comC));
		// clean
		cholmod_l_free_triplet(&T, &(comC));
		cholmod_l_free_sparse(&AcholC, &(comC));
		factorized = true;
	}
	cholmod_dense* B = cholmod_l_zeros(ncols, 1, LC->xtype, &(comC));
	Real*          B_x = (Real*)B->x;
	for (int k = 0; k < ncols; k++)
		B_x[k] = RHS[k] + RHSvol[k];

	ex = cholmod_l_solve(CHOLMOD_A, LC, B, &(comC));
	Real* e_x = (Real*)ex->x;

	for (auto cellIt = pCells.begin(); cellIt != pCells.end(); cellIt++) {
		const CellHandle& cell = *cellIt;
		cell->info().p() = e_x[cell->info().index];
	}
	//clean
	cholmod_l_free_dense(&B, &(comC));

#endif
}

void TwoPhaseFlowEngine::initialization()
{
	scene = Omega::instance().getScene().get(); //here define the pointer to Yade's scene -necessary if the engine is used outside O.engines
	setPositionsBuffer(true);                   //copy sphere positions in a buffer...
	if (!keepTriangulation) {
		buildTriangulation(0.0, *solver);
	} //create a triangulation and initialize pressure in the elements (connecting with W-reservoir), everything will be contained in "solver"
	  // 		initializeCellIndex();//initialize cell index
	  // 		if(isInvadeBoundary) {computePoreThroatRadius();}
	// // 		else {computePoreThroatRadiusTrickyMethod1();}//save pore throat radius before drainage. Thomas, here you can also revert this to computePoreThroatCircleRadius().

	// Determine the entry-pressure
	if (entryPressureMethod == 1 && isInvadeBoundary) {
		computePoreThroatRadiusMethod1();
	} //MS-P method
	else if (entryPressureMethod == 1 && isInvadeBoundary == false) {
		computePoreThroatRadiusTrickyMethod1();
	} //MS-P method
	else if (entryPressureMethod == 2) {
		computePoreThroatRadiusMethod2();
	} //Inscribed circle}
	else if (entryPressureMethod == 3) {
		computePoreThroatRadiusMethod3();
	} //Area equivalent circle}
	else if (entryPressureMethod > 3) {
		cout << endl << "ERROR - Method for determining the entry pressure does not exist";
	}

	computePoreBodyRadius(); //save pore body radius before imbibition
	computePoreBodyVolume(); //save capillary volume of all cells, for fast calculating saturation. Also save the porosity of each cell.
	computeSolidLine();      //save cell->info().solidLine[j][y]
	initializeReservoirs();  //initial pressure, reservoir flags and local pore saturation
	if (isCellLabelActivated) updateCellLabel();
	solver->noCache = true;
}

void TwoPhaseFlowEngine::computePoreBodyVolume()
{
	initializeVolumes(*solver);
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		cell->info().poreBodyVolume = math::abs(cell->info().volume()) - math::abs(solver->volumeSolidPore(cell));
		cell->info().porosity = cell->info().poreBodyVolume / math::abs(cell->info().volume());
	}
}

void TwoPhaseFlowEngine::computePoreThroatRadiusMethod2()
{
	//Calculate the porethroat radii of the inscribed sphere in each pore-body.
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (unsigned int i = 0; i < 4; i++) {
			cell->info().poreThroatRadius[i] = math::abs(solver->computeEffectiveRadius(cell, i));
		}
	}
}

void TwoPhaseFlowEngine::computePoreThroatRadiusMethod3()
{
	//Calculate the porethroat radii of the surface equal circle of a throat
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (unsigned int i = 0; i < 4; i++) {
			cell->info().poreThroatRadius[i] = solver->computeEquivalentRadius(cell, i);
		}
	}
}

void TwoPhaseFlowEngine::computePoreBodyRadius()
{
	// This routine finds the radius of the inscribed sphere within each pore-body
	// Following Mackay et al., 1972.
	Real         d01 = 0.0, d02 = 0.0, d03 = 0.0, d12 = 0.0, d13 = 0.0, d23 = 0.0, Rin = 0.0, r0 = 0.0, r1 = 0.0, r2 = 0.0, r3 = 0.0;
	bool         check = false;
	unsigned int i = 0;
	Real         dR = 0.0, tempR = 0.0;
	bool         initialSign = false; //False = negative, true is positive
	bool         first2 = true;

	MatrixXr M(6, 6);

	FOREACH(CellHandle & cell, solver->T[solver->currentTes].cellHandles)
	{
		//Distance between multiple particles, can be done more efficient
		d01 = d02 = d03 = d12 = d13 = d23 = r0 = r1 = r2 = r3 = 0.0;

		d01 = pow((cell->vertex(0)->point().x() - cell->vertex(1)->point().x()), 2)
		        + pow((cell->vertex(0)->point().y() - cell->vertex(1)->point().y()), 2)
		        + pow((cell->vertex(0)->point().z() - cell->vertex(1)->point().z()), 2);

		d02 = pow((cell->vertex(0)->point().x() - cell->vertex(2)->point().x()), 2)
		        + pow((cell->vertex(0)->point().y() - cell->vertex(2)->point().y()), 2)
		        + pow((cell->vertex(0)->point().z() - cell->vertex(2)->point().z()), 2);

		d03 = pow((cell->vertex(0)->point().x() - cell->vertex(3)->point().x()), 2)
		        + pow((cell->vertex(0)->point().y() - cell->vertex(3)->point().y()), 2)
		        + pow((cell->vertex(0)->point().z() - cell->vertex(3)->point().z()), 2);

		d12 = pow((cell->vertex(1)->point().x() - cell->vertex(2)->point().x()), 2)
		        + pow((cell->vertex(1)->point().y() - cell->vertex(2)->point().y()), 2)
		        + pow((cell->vertex(1)->point().z() - cell->vertex(2)->point().z()), 2);

		d13 = pow((cell->vertex(1)->point().x() - cell->vertex(3)->point().x()), 2)
		        + pow((cell->vertex(1)->point().y() - cell->vertex(3)->point().y()), 2)
		        + pow((cell->vertex(1)->point().z() - cell->vertex(3)->point().z()), 2);

		d23 = pow((cell->vertex(2)->point().x() - cell->vertex(3)->point().x()), 2)
		        + pow((cell->vertex(2)->point().y() - cell->vertex(3)->point().y()), 2)
		        + pow((cell->vertex(2)->point().z() - cell->vertex(3)->point().z()), 2);

		//Radii of the particles
		r0 = sqrt(cell->vertex(0)->point().weight());
		r1 = sqrt(cell->vertex(1)->point().weight());
		r2 = sqrt(cell->vertex(2)->point().weight());
		r3 = sqrt(cell->vertex(3)->point().weight());

		//Fill coefficient matrix
		M(0, 0) = 0.0;
		M(1, 0) = d01;
		M(2, 0) = d02;
		M(3, 0) = d03;
		M(4, 0) = pow((r0 + Rin), 2);
		M(5, 0) = 1.0;
		M(0, 1) = d01;
		M(1, 1) = 0.0;
		M(2, 1) = d12;
		M(3, 1) = d13;
		M(4, 1) = pow((r1 + Rin), 2);
		M(5, 1) = 1.0;
		M(0, 2) = d02;
		M(1, 2) = d12;
		M(2, 2) = 0.0;
		M(3, 2) = d23;
		M(4, 2) = pow((r2 + Rin), 2);
		M(5, 2) = 1.0;
		M(0, 3) = d03;
		M(1, 3) = d13;
		M(2, 3) = d23;
		M(3, 3) = 0.0;
		M(4, 3) = pow((r3 + Rin), 2);
		M(5, 3) = 1.0;
		M(0, 4) = pow((r0 + Rin), 2);
		M(1, 4) = pow((r1 + Rin), 2);
		M(2, 4) = pow((r2 + Rin), 2);
		M(3, 4) = pow((r3 + Rin), 2);
		M(4, 4) = 0.0;
		M(5, 4) = 1.0;
		M(0, 5) = 1.0;
		M(1, 5) = 1.0;
		M(2, 5) = 1.0;
		M(3, 5) = 1.0;
		M(4, 5) = 1.0;
		M(5, 5) = 0.0;

		i = 0;
		check = false;
		dR = Rin = 0.0 + (min(r0, min(r1, min(r2, r3))) / 50.0); //Estimate an initial dR
		first2 = true;
		//Iterate untill check = true, such that an accurate answer as been found
		while (check == false) {
			i = i + 1;
			tempR = Rin;
			Rin = Rin + dR;

			M(4, 0) = pow((r0 + Rin), 2);
			M(4, 1) = pow((r1 + Rin), 2);
			M(4, 2) = pow((r2 + Rin), 2);
			M(4, 3) = pow((r3 + Rin), 2);
			M(0, 4) = pow((r0 + Rin), 2);
			M(1, 4) = pow((r1 + Rin), 2);
			M(2, 4) = pow((r2 + Rin), 2);
			M(3, 4) = pow((r3 + Rin), 2);

			if (first2) {
				first2 = false;
				if (M.determinant() < 0.0) { initialSign = false; } //Initial D is negative
				if (M.determinant() > 0.0) { initialSign = true; }  // Initial D is positive
			}

			if (math::abs(M.determinant()) < 1E-100) { check = true; } //TODO:M.determinant should be converted to dimensionless.

			if ((initialSign == true) && (check == false)) {
				if (M.determinant() < 0.0) {
					Rin = Rin - dR;
					dR = dR / 2.0;
				}
			}

			if ((initialSign == false) && (check == false)) {
				if (M.determinant() > 0.0) {
					Rin = Rin - dR;
					dR = dR / 2.0;
				}
			}

			if (solver->debugOut) { cout << endl << i << " " << Rin << " " << dR << " " << M.determinant(); }
			if (i > 4000) {
				cout << endl << "error, finding solution takes too long cell:" << cell->info().id;
				check = true;
			}
			if (math::abs(tempR - Rin) / Rin < 0.001) { check = true; }
		}
		cell->info().poreBodyRadius = Rin;
	}
}

void TwoPhaseFlowEngine::computePoreThroatRadiusMethod1()
{
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	CellHandle          neighbourCell;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (int j = 0; j < 4; j++) {
			neighbourCell = cell->neighbor(j);
			if (!tri.is_infinite(neighbourCell)) {
				cell->info().poreThroatRadius[j] = computeEffPoreThroatRadius(cell, j);
				neighbourCell->info().poreThroatRadius[tri.mirror_index(cell, j)] = cell->info().poreThroatRadius[j];
			}
		}
	}
}
Real TwoPhaseFlowEngine::computeEffPoreThroatRadius(CellHandle cell, int j)
{
	Real       rInscribe = math::abs(solver->computeEffectiveRadius(cell, j));
	CellHandle cellh = CellHandle(cell);
	int        facetNFictious = solver->detectFacetFictiousVertices(cellh, j);
	Real       r;
	if (facetNFictious == 0) {
		r = computeEffPoreThroatRadiusFine(cell, j);
	} else
		r = rInscribe;
	return r;
}
Real TwoPhaseFlowEngine::computeEffPoreThroatRadiusFine(CellHandle cell, int j)
{
	RTriangulation& tri = solver->T[solver->currentTes].Triangulation();
	if (tri.is_infinite(cell->neighbor(j))) return 0;

	Vector3r pos[3]; //solid pos
	Real     r[3];   //solid radius

	for (int i = 0; i < 3; i++) {
		pos[i] = makeVector3r(cell->vertex(facetVertices[j][i])->point().point());
		r[i] = sqrt(cell->vertex(facetVertices[j][i])->point().weight());
	}
	return computeMSPRcByPosRadius(pos[0], r[0], pos[1], r[1], pos[2], r[2]);
}
Real TwoPhaseFlowEngine::computeMSPRcByPosRadius(
        const Vector3r& posA, const Real& rA, const Vector3r& posB, const Real& rB, const Vector3r& posC, const Real& rC)
{
	Real e[3]; //edges of triangulation
	Real g[3]; //gap radius between solid

	e[0] = (posB - posC).norm();
	e[1] = (posC - posA).norm();
	e[2] = (posB - posA).norm();
	g[0] = ((e[0] - rB - rC) > 0) ? 0.5 * (e[0] - rB - rC) : 0;
	g[1] = ((e[1] - rC - rA) > 0) ? 0.5 * (e[1] - rC - rA) : 0;
	g[2] = ((e[2] - rA - rB) > 0) ? 0.5 * (e[2] - rA - rB) : 0;

	Real rmin = (math::max(g[0], math::max(g[1], g[2])) == 0) ? 1.0e-11 : math::max(g[0], math::max(g[1], g[2]));
	Real rmax = computeEffRcByPosRadius(posA, rA, posB, rB, posC, rC);
	if (rmin > rmax) { cerr << "WARNING! rmin>rmax. rmin=" << rmin << " ,rmax=" << rmax << endl; }

	Real deltaForceRMin = computeDeltaForce(posA, rA, posB, rB, posC, rC, rmin);
	Real deltaForceRMax = computeDeltaForce(posA, rA, posB, rB, posC, rC, rmax);
	Real effPoreRadius;

	if (deltaForceRMin > deltaForceRMax) {
		effPoreRadius = rmax;
	} else if (deltaForceRMax < 0) {
		effPoreRadius = rmax;
	} else if (deltaForceRMin > 0) {
		effPoreRadius = rmin;
	} else {
		effPoreRadius = bisection(posA, rA, posB, rB, posC, rC, rmin, rmax);
	}
	return effPoreRadius;
}
Real TwoPhaseFlowEngine::bisection(
        const Vector3r& posA, const Real& rA, const Vector3r& posB, const Real& rB, const Vector3r& posC, const Real& rC, Real a, Real b)
{
	Real m = 0.5 * (a + b);
	if (math::abs(b - a) > computeEffRcByPosRadius(posA, rA, posB, rB, posC, rC) * 1.0e-6) {
		if (computeDeltaForce(posA, rA, posB, rB, posC, rC, m) * computeDeltaForce(posA, rA, posB, rB, posC, rC, a) < 0) {
			b = m;
			return bisection(posA, rA, posB, rB, posC, rC, a, b);
		} else {
			a = m;
			return bisection(posA, rA, posB, rB, posC, rC, a, b);
		}
	} else {
		return m;
	}
}
Real TwoPhaseFlowEngine::computeDeltaForce(
        const Vector3r& posA, const Real& rA, const Vector3r& posB, const Real& rB, const Vector3r& posC, const Real& rC, Real r)
{
	Real rRc[3];    //r[i] + r (r: capillary radius)
	Real e[3];      //edges of triangulation
	Real rad[4][3]; //angle in radian

	rRc[0] = rA + r;
	rRc[1] = rB + r;
	rRc[2] = rC + r;

	e[0] = (posB - posC).norm();
	e[1] = (posC - posA).norm();
	e[2] = (posB - posA).norm();

	rad[3][0] = acos(((posB - posA).dot(posC - posA)) / (e[2] * e[1]));
	rad[3][1] = acos(((posC - posB).dot(posA - posB)) / (e[0] * e[2]));
	rad[3][2] = acos(((posA - posC).dot(posB - posC)) / (e[1] * e[0]));

	rad[0][0] = computeTriRadian(e[0], rRc[1], rRc[2]);
	rad[0][1] = computeTriRadian(rRc[2], e[0], rRc[1]);
	rad[0][2] = computeTriRadian(rRc[1], rRc[2], e[0]);

	rad[1][0] = computeTriRadian(rRc[2], e[1], rRc[0]);
	rad[1][1] = computeTriRadian(e[1], rRc[0], rRc[2]);
	rad[1][2] = computeTriRadian(rRc[0], rRc[2], e[1]);

	rad[2][0] = computeTriRadian(rRc[1], e[2], rRc[0]);
	rad[2][1] = computeTriRadian(rRc[0], rRc[1], e[2]);
	rad[2][2] = computeTriRadian(e[2], rRc[0], rRc[1]);

	Real lNW = (rad[0][0] + rad[1][1] + rad[2][2]) * r;
	Real lNS = (rad[3][0] - rad[1][0] - rad[2][0]) * rA + (rad[3][1] - rad[2][1] - rad[0][1]) * rB + (rad[3][2] - rad[1][2] - rad[0][2]) * rC;
	Real lInterface = lNW + lNS;

	Real sW0 = 0.5 * rRc[1] * rRc[2] * sin(rad[0][0]) - 0.5 * rad[0][0] * pow(r, 2) - 0.5 * rad[0][1] * pow(rB, 2) - 0.5 * rad[0][2] * pow(rC, 2);
	Real sW1 = 0.5 * rRc[2] * rRc[0] * sin(rad[1][1]) - 0.5 * rad[1][1] * pow(r, 2) - 0.5 * rad[1][2] * pow(rC, 2) - 0.5 * rad[1][0] * pow(rA, 2);
	Real sW2 = 0.5 * rRc[0] * rRc[1] * sin(rad[2][2]) - 0.5 * rad[2][2] * pow(r, 2) - 0.5 * rad[2][0] * pow(rA, 2) - 0.5 * rad[2][1] * pow(rB, 2);
	Real sW = sW0 + sW1 + sW2;

	CVector facetSurface = 0.5 * CGAL::cross_product(makeCgVect(posA - posC), makeCgVect(posB - posC));
	Real    sVoid = sqrt(facetSurface.squared_length()) - (0.5 * rad[3][0] * pow(rA, 2) + 0.5 * rad[3][1] * pow(rB, 2) + 0.5 * rad[3][2] * pow(rC, 2));
	Real    sInterface = sVoid - sW;

	Real deltaF = lInterface - sInterface / r; //deltaF=surfaceTension*(perimeterPore - areaPore/rCap)
	return deltaF;
}
//calculate radian with law of cosines. (solve $\alpha$)
Real TwoPhaseFlowEngine::computeTriRadian(Real a, Real b, Real c)
{
	Real cosAlpha = (pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * b * c);
	if (cosAlpha > 1.0) { cosAlpha = 1.0; }
	if (cosAlpha < -1.0) { cosAlpha = -1.0; }
	Real alpha = acos(cosAlpha);
	return alpha;
}

void TwoPhaseFlowEngine::savePhaseVtk(const char* folder, bool withBoundaries)
{
	vector<int>
	            allIds; //an ordered list of cell ids (from begin() to end(), for vtk table lookup), some ids will appear multiple times since boundary cells are splitted into multiple tetrahedra
	vector<int> fictiousN;
	bool        initNoCache = solver->noCache;
	solver->noCache = false;

	static unsigned int number = 0;
	char                filename[250];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	sprintf(filename, "%s/out_%d.vtk", folder, number++);
	basicVTKwritter vtkfile(0, 0);
	solver->saveMesh(vtkfile, withBoundaries, allIds, fictiousN, filename);
	solver->noCache = initNoCache;

	vtkfile.begin_data("Pressure", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().p());
	vtkfile.end_data();

	vtkfile.begin_data("fictious", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(fictiousN[kk]);
	vtkfile.end_data();

	vtkfile.begin_data("id", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(allIds[kk]);
	vtkfile.end_data();

#define SAVE_CELL_INFO(INFO)                                                                                                                                   \
	vtkfile.begin_data(#INFO, CELL_DATA, SCALARS, FLOAT);                                                                                                  \
	for (unsigned kk = 0; kk < allIds.size(); kk++)                                                                                                        \
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().INFO);                                                                \
	vtkfile.end_data();
	SAVE_CELL_INFO(saturation)
	SAVE_CELL_INFO(hasInterface)
	SAVE_CELL_INFO(Pcondition)
	SAVE_CELL_INFO(flux)
	SAVE_CELL_INFO(mergedID)
	SAVE_CELL_INFO(accumulativeDV)
	SAVE_CELL_INFO(porosity)
	SAVE_CELL_INFO(label)
}

void TwoPhaseFlowEngine::computePoreThroatRadiusTrickyMethod1()
{
	computePoreThroatRadiusMethod1();
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	CellHandle          neighbourCell;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (int j = 0; j < 4; j++) {
			neighbourCell = cell->neighbor(j);
			if (cell->info().isFictious && neighbourCell->info().isFictious) {
				cell->info().poreThroatRadius[j] = -1.0;
				neighbourCell->info().poreThroatRadius[tri.mirror_index(cell, j)] = cell->info().poreThroatRadius[j];
			}
		}
	}
}

void TwoPhaseFlowEngine::computeSolidLine()
{
	RTriangulation&     Tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (int j = 0; j < 4; j++) {
			solver->lineSolidPore(cell, j);
		}
	}
	if (solver->debugOut) { cout << "----computeSolidLine-----." << endl; }
}

void TwoPhaseFlowEngine::initializeReservoirs()
{
	boundaryConditions(*solver);
	solver->pressureChanged = true;
	solver->reApplyBoundaryConditions();
	///keep boundingCells[2] as W-reservoir.
	for (FlowSolver::VCellIterator it = solver->boundingCells[2].begin(); it != solver->boundingCells[2].end(); it++) {
		(*it)->info().isWRes = true;
		(*it)->info().isNWRes = false;
		(*it)->info().saturation = 1.0;
	}
	///keep boundingCells[3] as NW-reservoir.
	for (FlowSolver::VCellIterator it = solver->boundingCells[3].begin(); it != solver->boundingCells[3].end(); it++) {
		(*it)->info().isNWRes = true;
		(*it)->info().isWRes = false;
		(*it)->info().saturation = 0.0;
	}

	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	///if we start from drainage
	if (drainageFirst) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().Pcondition) continue;
			cell->info().p() = bndCondValue[2];
			cell->info().isWRes = true;
			cell->info().isNWRes = false;
			cell->info().saturation = 1.0;
		}
	}
	///if we start from imbibition
	if (!drainageFirst) {
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().Pcondition) continue;
			cell->info().p() = bndCondValue[3];
			cell->info().isWRes = false;
			cell->info().isNWRes = true;
			cell->info().saturation = 0.0;
		}
	}
	if (solver->debugOut) { cout << "----initializeReservoirs----" << endl; }
}


void TwoPhaseFlowEngine::savePoreNetwork(const char* folder)
{
	//Open relevant files
	std::ofstream filePoreBodyRadius;
	std::cout << "Opening File: "
	          << "PoreBodyRadius" << std::endl;
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	filePoreBodyRadius.open(string(folder) + "/PoreBodyRadius.txt", std::ios::trunc);
	if (!filePoreBodyRadius.is_open()) {
		std::cerr << "Error opening file ["
		          << "PoreBodyRadius" << ']' << std::endl;
		return;
	}

	std::ofstream filePoreBoundary;
	std::cout << "Opening File: "
	          << "PoreBoundary" << std::endl;
	filePoreBoundary.open(string(folder) + "/PoreBoundaryIndex.txt", std::ios::trunc);
	if (!filePoreBoundary.is_open()) {
		std::cerr << "Error opening file ["
		          << "PoreBoundary" << ']' << std::endl;
		return;
	}

	std::ofstream filePoreBodyVolume;
	std::cout << "Opening File: "
	          << "PoreBodyVolume" << std::endl;
	filePoreBodyVolume.open(string(folder) + "/PoreBodyVolume.txt", std::ios::trunc);
	if (!filePoreBodyVolume.is_open()) {
		std::cerr << "Error opening file ["
		          << "PoreBodyVolume" << ']' << std::endl;
		return;
	}
	std::ofstream fileLocation;
	std::cout << "Opening File: "
	          << "Location" << std::endl;
	fileLocation.open(string(folder) + "/PoreBodyLocation.txt", std::ios::trunc);
	if (!fileLocation.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileLocation" << ']' << std::endl;
		return;
	}

	std::ofstream fileNeighbor;
	std::cout << "Opening File: "
	          << "fileNeighbor" << std::endl;
	fileNeighbor.open(string(folder) + "/PoreBodyNeighbor.txt", std::ios::trunc);
	if (!fileNeighbor.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileNeighbor" << ']' << std::endl;
		return;
	}

	std::ofstream fileThroatRadius;
	std::cout << "Opening File: "
	          << "fileThroatRadius" << std::endl;
	fileThroatRadius.open(string(folder) + "/throatRadius.txt", std::ios::trunc);
	if (!fileThroatRadius.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileThroatRadius" << ']' << std::endl;
		return;
	}
	std::ofstream fileThroats;
	std::cout << "Opening File: "
	          << "fileThroats" << std::endl;
	fileThroats.open(string(folder) + "/throatConnectivityPoreBodies.txt", std::ios::trunc);
	if (!fileThroats.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileThroats" << ']' << std::endl;
		return;
	}
	std::ofstream fileThroatFluidArea;
	std::cout << "Opening File: "
	          << "fileThroatFluidArea" << std::endl;
	fileThroatFluidArea.open(string(folder) + "/fileThroatFluidArea.txt", std::ios::trunc);
	if (!fileThroatFluidArea.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileThroatFluidArea" << ']' << std::endl;
		return;
	}
	std::ofstream fileHydraulicRadius;
	std::cout << "Opening File: "
	          << "fileHydraulicRadius" << std::endl;
	fileHydraulicRadius.open(string(folder) + "/fileHydraulicRadius.txt", std::ios::trunc);
	if (!fileHydraulicRadius.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileHydraulicRadius" << ']' << std::endl;
		return;
	}
	std::ofstream fileConductivity;
	std::cout << "Opening File: "
	          << "fileConductivity" << std::endl;
	fileConductivity.open(string(folder) + "/fileConductivity.txt", std::ios::trunc);
	if (!fileConductivity.is_open()) {
		std::cerr << "Error opening file ["
		          << "fileConductivity" << ']' << std::endl;
		return;
	}

	//Extract pore network based on triangulation
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isGhost == false && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
			filePoreBodyRadius << cell->info().poreBodyRadius << '\n';
			filePoreBodyVolume << cell->info().poreBodyRadius << '\n';
			CVector center(0, 0, 0);
			Real    count = 0.0;
			for (int k = 0; k < 4; k++) {
				if (cell->vertex(k)->info().id() > 5) {
					center = center + (cell->vertex(k)->point().point() - CGAL::ORIGIN);
					count = count + 1.0;
				}
			}
			if (count != 0.0) { center = center * (1. / count); }
			fileLocation << center << '\n';
			for (unsigned int i = 0; i < 4; i++) {
				if (cell->neighbor(i)->info().isGhost == false
				    && cell->neighbor(i)->info().id < solver->T[solver->currentTes].cellHandles.size()
				    && (cell->info().id < cell->neighbor(i)->info().id)) {
					fileNeighbor << cell->neighbor(i)->info().id << '\n';
					fileThroatRadius << cell->info().poreThroatRadius[i] << '\n';
					fileThroats << cell->info().id << " " << cell->neighbor(i)->info().id << '\n';
					const CVector& Surfk = cell->info().facetSurfaces[i];
					Real           area = sqrt(Surfk.squared_length());
					fileThroatFluidArea << cell->info().facetFluidSurfacesRatio[i] * area << '\n';
					fileHydraulicRadius << 2.0 * solver->computeHydraulicRadius(cell, i) << '\n';
					fileConductivity << cell->info().kNorm()[i] << '\n';
				}
			}

			if (cell->info().isFictious == 1 && cell->info().isGhost == false
			    && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
				//add boundary condition
				if (cell->info().isFictious == 1
				    && (cell->vertex(0)->info().id() == 3 || cell->vertex(1)->info().id() == 3 || cell->vertex(2)->info().id() == 3
				        || cell->vertex(3)->info().id() == 3)) {
					filePoreBoundary << "3" << '\n';
				} else if (
				        cell->info().isFictious == 1
				        && (cell->vertex(0)->info().id() == 2 || cell->vertex(1)->info().id() == 2 || cell->vertex(2)->info().id() == 2
				            || cell->vertex(3)->info().id() == 2)) {
					filePoreBoundary << "2" << '\n';
				} else {
					filePoreBoundary << "2" << '\n';
				}
			}
			if (cell->info().isFictious == 0 && cell->info().isGhost == false
			    && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
				filePoreBoundary << "0" << '\n';
			}
		}
	}
	fileThroatFluidArea.close();
	fileHydraulicRadius.close();
	fileConductivity.close();
	filePoreBodyRadius.close();
	filePoreBoundary.close();
	filePoreBodyVolume.close();
	fileLocation.close();
	fileNeighbor.close();
	fileThroatRadius.close();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//..............................................................Library............................................................//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Real TwoPhaseFlowEngine::getKappa(int numberFacets) const
{
	if (numberFacets == 0) {
		return 0;
		cout << endl << "Pore with zero throats? Check your data";
	} else {
		Real kappa = 0.0;
		if (numberFacets == 4) {
			kappa = 3.8716;
		} //Tetrahedra
		else if (numberFacets == 6) {
			kappa = 8.7067;
		} //Octahedra
		else if (numberFacets == 8) {
			kappa = 6.7419;
		} //Cube
		else if (numberFacets == 10) {
			kappa = 5.150;
		} //Octahedron + hexahedron
		else if (numberFacets == 12) {
			kappa = 24.105;
		} //Icosahedra
		else if (numberFacets == 20) {
			kappa = 22.866;
		} // Dodecahedron
		else {
			kappa = 1.2591 * float(numberFacets) - 1.1041;
		} //Other pore shapes
		return kappa;
	}
}

Real TwoPhaseFlowEngine::getChi(int numberFacets) const
{
	if (numberFacets == 0) {
		return 0;
		cout << endl << "Pore with zero throats? Check your data";
	} else {
		Real chi = 0.0;
		if (numberFacets == 4) {
			chi = 0.416;
		} //Tetrahedra
		else if (numberFacets == 6) {
			chi = 0.525;
		} //Octahedra
		else if (numberFacets == 8) {
			chi = 0.500;
		} //Cube
		else if (numberFacets == 10) {
			chi = 0.4396;
		} //Octahedron + hexahedron
		else if (numberFacets == 12) {
			chi = 0.583;
		} //Icosahedra
		else if (numberFacets == 20) {
			chi = 0.565;
		} // Dodecahedron
		else {
			chi = 0.0893 * math::log(numberFacets) + 0.326;
		} //Other pore shapes
		return chi;
	}
}


Real TwoPhaseFlowEngine::getLambda(int numberFacets) const
{
	if (numberFacets == 0) {
		return 0;
		cout << endl << "Pore with zero throats? Check your data";
	} else {
		Real lambda = 0.0;
		if (numberFacets == 4) {
			lambda = 2.0396;
		} //Tetrahedra
		else if (numberFacets == 6) {
			lambda = 1.2849;
		} //Octahedra
		else if (numberFacets == 8) {
			lambda = 1;
		} //Cube
		else if (numberFacets == 10) {
			lambda = 0.77102;
		} //Octahedron + hexahedron
		else if (numberFacets == 12) {
			lambda = 0.771025;
		} //Icosahedra
		else if (numberFacets == 20) {
			lambda = 0.50722;
		} // Dodecahedron
		else {
			lambda = 7.12 * math::pow(numberFacets, -0.89);
		} //Other pore shapes
		return lambda;
	}
}

Real TwoPhaseFlowEngine::getN(int numberFacets) const
{
	if (numberFacets == 0) {
		return 0;
		cout << endl << "Pore with zero throats? Check your data";
	} else {
		Real n = 0.0;
		if (numberFacets == 4) {
			n = 6.0;
		} //Tetrahedra
		else if (numberFacets == 6) {
			n = 12.0;
		} //Octahedra
		else if (numberFacets == 8) {
			n = 8.0;
		} //Cube
		else if (numberFacets == 10) {
			n = 12.0; /*cout << endl << "number of edges requested for octa + hexahedron!";*/
		}                 //Octahedron + hexahedron NOTE this should not be requested for calculations!
		else if (numberFacets == 12) {
			n = 30.0;
		} //Icosahedra
		else if (numberFacets == 20) {
			n = 30.0;
		} // Dodecahedron
		else {
			n = 1.63 * Real(numberFacets);
		} //Other pore shapes
		return n;
	}
}

Real TwoPhaseFlowEngine::getDihedralAngle(int numberFacets) const
{ //given in radians which is reported as tetha in manuscript Sweijen et al.,
	if (numberFacets == 0) {
		return 0;
		cout << endl << "Pore with zero throats? Check your data";
	} else {
		Real DihedralAngle = 0.0;
		if (numberFacets == 4) {
			DihedralAngle = 1.0 * math::atan(2.0 * math::sqrt(2.0));
		} //Tetrahedra
		else if (numberFacets == 6) {
			DihedralAngle = 1.0 * math::acos(-1.0 / 3.0);
		} //Octahedra
		else if (numberFacets == 8) {
			DihedralAngle = 0.5 * 3.1415926535;
		} //Cube
		else if (numberFacets == 10) {
			DihedralAngle = (1. / 4.) * 3.1415926535; /*cout << endl << "dihedral angle requested for octa + hexahedron!";*/
		}                                                 //Octahedron + hexahedron NOTE this should not be requested for calculations!
		else if (numberFacets == 12) {
			DihedralAngle = math::acos((-1.0 / 3.0) * math::sqrt(5.0));
		} //Icosahedra
		else if (numberFacets == 20) {
			DihedralAngle = math::acos((-1.0 / 5.0) * math::sqrt(5.0));
		} // Dodecahedron
		else {
			DihedralAngle = (1. / 4.) * 3.1415926535;
		} //Other pore shapes
		return DihedralAngle;
	}
}

Real TwoPhaseFlowEngine::getConstantC3(CellHandle cell) const
{
	Real c1 = 54.92 * math::pow(Real(cell->info().numberFacets), -1.14);
	if (cell->info().numberFacets == 4) { c1 = 8.291; }
	if (cell->info().numberFacets == 6) { c1 = 2.524; }
	if (cell->info().numberFacets == 8) { c1 = 2.524; }
	if (cell->info().numberFacets == 10) { c1 = 6.532; }
	if (cell->info().numberFacets == 12) { c1 = 6.087; }
	if (cell->info().numberFacets == 20) { c1 = 0.394; }


	Real c3 = c1 * math::pow(2.0 * surfaceTension, 3) / cell->info().mergedVolume;
	return c3;
}

Real TwoPhaseFlowEngine::getConstantC4(CellHandle cell) const
{
	Real c2 = 4.85 * math::pow(Real(cell->info().numberFacets), -1.19);
	if (cell->info().numberFacets == 4) { c2 = 1.409; }
	if (cell->info().numberFacets == 6) { c2 = 0.353; }
	if (cell->info().numberFacets == 8) { c2 = 0.644; }
	if (cell->info().numberFacets == 10) { c2 = 0.462; }
	if (cell->info().numberFacets == 12) { c2 = 0.0989; }
	if (cell->info().numberFacets == 20) { c2 = 0.245; }
	Real c4 = c2 * math::pow(2.0 * surfaceTension, 3) / math::pow(Real(cell->info().mergedVolume), 2. / 3.);
	return c4;
}

Real TwoPhaseFlowEngine::dsdp(CellHandle cell, Real pw)
{
	if (pw == 0) { std::cout << endl << "Error! water pressure is zero, while computing capillary pressure ... cellId= " << cell->info().id; }
	Real exp = math::exp(-1 * getKappa(cell->info().numberFacets) * cell->info().saturation);
	Real dsdp2 = (1.0 / cell->info().thresholdPressure) * math::pow((1.0 - exp), 2.0) / (getKappa(cell->info().numberFacets) * exp);
	//       if(math::abs(dsdp2) > 1e10){ std::cerr << "Huge dsdp! : "<< dsdp2 << " " << exp << " "<< cell->info().thresholdPressure << " " << getKappa(cell->info().numberFacets);}


	//     Real dsdp2 = (3.0 * getConstantC3(cell) - 2.0 * getConstantC4(cell) * pw) / math::pow(pw,4);
	if (dsdp2 != dsdp2) {
		std::cerr << endl
		          << "Error! sat in dsdp is nan: " << cell->info().saturation << " kappa:" << getKappa(cell->info().numberFacets) << " exp: " << exp
		          << " mergedVolume=" << cell->info().mergedVolume << " pthreshold=" << cell->info().thresholdPressure;
	}
	if (dsdp2 < 0.0) {
		std::cerr << endl << "Error! dsdp is negative!" << dsdp2;
		dsdp2 = 0.0;
	}
	//     if(dsdp2 >1e6){std::cerr<<endl<< "Error! dsdp is huge!" << dsdp2; dsdp2 = 1e6;}
	return dsdp2;
}

Real TwoPhaseFlowEngine::poreSaturationFromPcS(CellHandle cell, Real pw)
{
	//Using equation: Pc = 2*surfaceTension / (Chi * PoreBodyVolume^(1/3) * (1-exp(-kappa * S)))
	Real s = truncationPrecision;
	if (-1 * pw > cell->info().thresholdPressure) {
		s = math::log(1.0 + cell->info().thresholdPressure / pw) / (-1.0 * getKappa(cell->info().numberFacets));
	}
	if (-1 * pw == cell->info().thresholdPressure) { s = cell->info().thresholdSaturation; }
	if (-1 * pw < cell->info().thresholdPressure) {
		if (!remesh && !firstDynTPF) {
			std::cerr << endl
			          << "Error! Requesting saturation while capillary pressure is below threshold value? " << pw << " "
			          << cell->info().thresholdPressure;
		}
		s = cell->info().thresholdSaturation;
	}

	if (s > 1.0 || s < 0.0) {
		std::cout << "Error, saturation from Pc(S) curve is not correct: " << s << " " << cell->info().poreId
		          << " log:" << math::log(1.0 + cell->info().thresholdPressure / pw) << " " << (-1.0 * getKappa(cell->info().numberFacets))
		          << " pw=" << pw << " " << cell->info().thresholdPressure;
		s = 1.0;
	}
	if (s != s) {
		std::cerr << endl
		          << "Error! sat in PcS is nan: " << s << "  " << pw << " " << getConstantC4(cell) << " " << getConstantC3(cell)
		          << " mergedVolume=" << cell->info().mergedVolume << " pthreshold=" << cell->info().thresholdPressure;
	}
	return s;
}


Real TwoPhaseFlowEngine::porePressureFromPcS(CellHandle cell, Real /*saturation*/)
{
	Real pw = -1.0 * cell->info().thresholdPressure / (1.0 - math::exp(-1 * getKappa(cell->info().numberFacets) * cell->info().saturation));
	if (math::exp(-1 * getKappa(cell->info().numberFacets) * cell->info().saturation) == 1.0) {
		std::cerr << endl << "Error! pw = -inf!" << cell->info().saturation;
	}
	if (pw > 0) {
		std::cout << "Pw is above 0! - error: " << pw << " id=" << cell->info().poreId << " pthr=" << cell->info().thresholdPressure
		          << " sat:" << cell->info().saturation << " kappa: " << getKappa(cell->info().numberFacets) << " "
		          << (1.0 - math::exp(-1 * getKappa(cell->info().numberFacets) * cell->info().saturation));
		pw = -1 * cell->info().thresholdPressure;
	}
	if (pw != pw) { std::cout << "Non existing capillary pressure!"; }
	//   if(pw < 100 * waterBoundaryPressure){std::cout << "huge PC!" << pw << " saturation=" << cell->info().saturation << " hasIFace=" << cell->info().hasInterface << " NWRES=" << cell->info().isNWRes; pw = 100 * waterBoundaryPressure;}
	return pw;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//....................................Merging of tetrahedra to find PUA............... ............................................//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TwoPhaseFlowEngine::actionMergingAlgorithm()
{
	//(1) merge tetrahedra together
	mergeCells();
	//(2) count the facets and update the mergedVolume
	countFacets();
	computeMergedVolumes();
	//(3) resolve unsolved pore throats that are too big
	adjustUnresolvedPoreThroatsAfterMerging();
	//(4) export statistics on the merging algorithm
	getMergedCellStats();
	//(5) check for volume conservation
	checkVolumeConservationAfterMergingAlgorithm();
}


void TwoPhaseFlowEngine::mergeCells()
{
	//This function finds the tetrahedra that belong together (i.e. pore throat is too big compared to pore body)
	//Start with worst case scenario, Rij / Ri > 200 (defined as criterion)towards maximumRatioPoreThroatoverPoreBody.
	//If Rij/Ri is too large, than tetrahedra are merged. Then the merged volume and nrFacets is updated.
	//When checkign a subsequent criterion, the updated values of merged volume and nr of facets is used.
	//The remaining unsolved Rij/Ri ratios are fixed later in the program.
	//A limitation of max 20 merged tetrahedra is used to prevent huge pores.

	int                 number = 1, ID = 0;
	Real                dC = 0.0, criterion = 200.0;
	bool                check = false;
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	maxIDMergedCells = 0;

	//Initialize Merged Volumes
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		cell->info().mergedVolume = cell->info().poreBodyVolume;
		cell->info().mergednr = 1;
		cell->info().mergedID = 0;
		cell->info().numberFacets = 4;
	}
	for (unsigned int i = 0; i < 110; i++) {
		if (i < 10) { dC = (200.0 - 50.0) / 9.0; }
		if (i >= 10) { dC = (50.0 - maximumRatioPoreThroatoverPoreBody) / 100.0; }
		if (i == 0) { dC = 0.0; }
		//Decrease the criteria for throat over body radius, so first merge the worst case scenarios.
		criterion = criterion - dC;
		if (debugTPF) { cout << endl << "criterion=" << criterion; }
		for (unsigned int j = 0; j < 5; j++) {
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().isGhost == false && cell->info().mergedID < solver->T[solver->currentTes].cellHandles.size()
				    && cell->info().isFictious == false && cell->info().mergednr < 20) {
					for (unsigned int ngb = 0; ngb < 4; ngb++) {
						if (cell->neighbor(ngb)->info().mergednr < 20) {
							if (cell->neighbor(ngb)->info().isGhost == false
							    && cell->neighbor(ngb)->info().mergedID < solver->T[solver->currentTes].cellHandles.size()
							    && cell->neighbor(ngb)->info().isFictious == false
							    && ((cell->info().mergedID == cell->neighbor(ngb)->info().mergedID && cell->info().mergedID != 0)
							        == false)) {
								if ((cell->info().poreThroatRadius[ngb]
								     / (getChi(cell->info().numberFacets) * math::pow(cell->info().mergedVolume, (1. / 3.))))
								    > criterion) {
									if (cell->info().mergedID == 0 && cell->neighbor(ngb)->info().mergedID == 0) {
										cell->info().mergedID = number;
										cell->neighbor(ngb)->info().mergedID = number;
										number = number + 1;
										countFacets();
										computeMergedVolumes();

									} else if (cell->info().mergedID == 0 && cell->neighbor(ngb)->info().mergedID != 0) {
										cell->info().mergedID = cell->neighbor(ngb)->info().mergedID;
										countFacets();
										computeMergedVolumes();
									} else if (cell->info().mergedID != 0 && cell->neighbor(ngb)->info().mergedID == 0) {
										cell->neighbor(ngb)->info().mergedID = cell->info().mergedID;
										countFacets();
										computeMergedVolumes();
									}
								}
							}
						}
					}
				}
			}
			countFacets();
			computeMergedVolumes();
		}
	}
	maxIDMergedCells = number;

	// RENUMBER THE MERGED CELLS
	for (unsigned int k = 1; k < maxIDMergedCells; k++) {
		check = false;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().mergedID == k) { check = true; }
		}
		if (check) {
			ID = ID + 1;
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().mergedID == k) { cell->info().mergedID = ID; }
			}
		}
	}
	maxIDMergedCells = ID + 1;
	if (debugTPF) { cout << endl << "EFFICIENT - RENUMBER MERGEDCELLS -- FROM: " << number << " TO: " << ID + 1; }
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		cell->info().mergedVolume = cell->info().poreBodyVolume;
	}
}


void TwoPhaseFlowEngine::computeMergedVolumes()
{
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	Real                volume = 0.0, summ = 0.0;
	for (unsigned int mergeID = 1; mergeID < maxIDMergedCells; mergeID++) {
		volume = 0.0;
		summ = 0.0;
		for (FiniteCellsIterator Mergecell = tri.finite_cells_begin(); Mergecell != cellEnd; Mergecell++) {
			if (Mergecell->info().mergedID == mergeID && Mergecell->info().isFictious == false && Mergecell->info().isGhost == false
			    && Mergecell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
				volume = volume + Mergecell->info().poreBodyVolume;
				summ = summ + 1.0;
			}
		}
		if (summ > 1.0) {
			for (FiniteCellsIterator Mergecell = tri.finite_cells_begin(); Mergecell != cellEnd; Mergecell++) {
				if (Mergecell->info().mergedID == mergeID && Mergecell->info().isFictious == false && Mergecell->info().isGhost == false
				    && Mergecell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
					Mergecell->info().poreBodyRadius = getChi(Mergecell->info().numberFacets) * math::pow(volume, (1. / 3.));
					Mergecell->info().mergedVolume = volume;
					Mergecell->info().mergednr = int(math::round(summ));
				}
			}
		}
		if (summ <= 1.0) {
			for (FiniteCellsIterator Mergecell = tri.finite_cells_begin(); Mergecell != cellEnd; Mergecell++) {
				if (Mergecell->info().mergedID == mergeID && Mergecell->info().isFictious == false && Mergecell->info().isGhost == false
				    && Mergecell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
					cout << endl
					     << "isMerged set to -1: " << Mergecell->info().id << " " << Mergecell->info().poreBodyRadius << " "
					     << Mergecell->info().poreThroatRadius[0] << " " << Mergecell->info().poreThroatRadius[1] << " "
					     << Mergecell->info().poreThroatRadius[2] << " " << Mergecell->info().poreThroatRadius[3];
					Mergecell->info().mergednr = 1;
					Mergecell->info().mergedID = 0;
				}
			}
		}
	}
}


void TwoPhaseFlowEngine::countFacets()
{
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	int                 summngb = 0;
	for (unsigned int k = 1; k < maxIDMergedCells; k++) {
		summngb = 0;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().mergedID == k && cell->info().isGhost == false && (cell->info().isFictious == false)
			    && (cell->info().id < solver->T[solver->currentTes].cellHandles.size())) {
				for (unsigned int i = 0; i < 4; i++) {
					if (cell->neighbor(i)->info().mergedID != cell->info().mergedID && cell->neighbor(i)->info().isGhost == false
					    && (cell->neighbor(i)->info().isFictious == false)
					    && (cell->neighbor(i)->info().id < solver->T[solver->currentTes].cellHandles.size())) {
						summngb = summngb + 1;
					}
				}
			}
		}
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().mergedID == k) {
				if (summngb < 4) { summngb = 4; } //Less than 4 throats is not supported -> boundary problems.
				cell->info().numberFacets = summngb;
			}
		}
	}
}


void TwoPhaseFlowEngine::getMergedCellStats() const
{
	std::array<Real, 26> countFacets = { 0 };
	std::array<Real, 30> countMergedNR = { 0 };
	int                  count = 0, countTot = 0;
	RTriangulation&      tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator  cellEnd = tri.finite_cells_end();
	//In this function the amount of tetrahedra per pore is counted and reported in the terminal.

	std::string NameDistributionInFacets = modelRunName;
	std::string NamedistibutionInMergedPoreUnits = modelRunName;
	NameDistributionInFacets.append("-distributionInFacets.txt");
	NamedistibutionInMergedPoreUnits.append("-distributionInMergedPoreUnits.txt");

	std::ofstream distributionInFacets;
	distributionInFacets.open(NameDistributionInFacets, std::ios::trunc);
	if (!distributionInFacets.is_open()) {
		std::cerr << "Error opening file ["
		          << "PoreBodyRadius" << ']' << std::endl;
		return;
	}

	std::ofstream distibutionInMergedPoreUnits;
	distibutionInMergedPoreUnits.open(NamedistibutionInMergedPoreUnits, std::ios::trunc);
	if (!distibutionInMergedPoreUnits.is_open()) {
		std::cerr << "Error opening file ["
		          << "PoreBoundary" << ']' << std::endl;
		return;
	}
	distributionInFacets << "The distribution in the number of pore throats per pore unit - table shows in the first column the number of pore throats and "
	                        "in the second column the total count"
	                     << '\n';
	distibutionInMergedPoreUnits << "The distribution in the number of tetrahedra per merged pore unit - table shows in the first column the number of "
	                                "merged tetrahedra and in the second column the total count"
	                             << '\n';

	countTot = 0;
	count = 0;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isFictious == false && cell->info().isGhost == false && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
			if (cell->info().numberFacets == 4) { count = count + 1; }
			countTot = countTot + 1;
		}
	}

	if (debugTPF) {
		cout << endl
		     << "Number of merged cells is:" << count << "of the total number" << countTot << " which is: " << (float(count) * 100.0 / float(countTot));
	}


	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isFictious == false && cell->info().isGhost == false && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
			if (cell->info().numberFacets < 30) {
				countFacets[cell->info().numberFacets - 4] = countFacets[cell->info().numberFacets - 4] + (1.0 / float(cell->info().mergednr));
			}
		}
	}

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isFictious == false && cell->info().isGhost == false && cell->info().id < solver->T[solver->currentTes].cellHandles.size()) {
			if (cell->info().mergednr < 30) {
				countMergedNR[cell->info().mergednr - 1] = countMergedNR[cell->info().mergednr - 1] + (1.0 / float(cell->info().mergednr));
			}
		}
	}

	for (unsigned int i = 0; i < countFacets.size(); i++) {
		if (debugTPF) { cout << endl << "nrFacets: " << (i + 4) << "-count:" << countFacets[i]; }
		distributionInFacets << (i + 4) << " " << countFacets[i] << '\n';
	}
	for (unsigned int i = 0; i < countMergedNR.size(); i++) {
		if (debugTPF) { cout << endl << "nrMergedUnits: " << i + 1 << "-count:" << countMergedNR[i]; }
		distibutionInMergedPoreUnits << (i + 1) << " " << countMergedNR[i] << '\n';
	}
	distributionInFacets.close();
	distibutionInMergedPoreUnits.close();
}


void TwoPhaseFlowEngine::adjustUnresolvedPoreThroatsAfterMerging()
{
	//Adjust the remaining pore throats, such that all throats are smaller than the pore units.
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	int                 count = 0, countTot = 0;
	for (unsigned int p = 0; p < 5; p++) {
		countTot = 0;
		count = 0;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().isGhost == false && cell->info().isFictious == false) {
				for (unsigned int i = 0; i < 4; i++) {
					if ((cell->info().mergedID != cell->neighbor(i)->info().mergedID
					     || (cell->info().mergedID == 0 && cell->neighbor(i)->info().mergedID == 0))
					    && cell->neighbor(i)->info().isGhost
					            == false /*&& cell->neighbor(i)->info().mergedID < solver->T[solver->currentTes].cellHandles.size()*/) {
						countTot = countTot + 1;
						if (cell->info().poreThroatRadius[i] >= maximumRatioPoreThroatoverPoreBody
						            * (getChi(cell->info().numberFacets)
						               * math::pow(
						                       cell->info().mergedVolume,
						                       (1.
						                        / 3.)))) { // if throat is larger than maximumRatioPoreThroatoverPoreBody time the pore body volume, then adjust pore throat radii
							count = count + 1;
							cell->info().poreThroatRadius[i] = math::min(
							        (maximumRatioPoreThroatoverPoreBody * getChi(cell->info().numberFacets)
							         * math::pow(cell->info().mergedVolume, (1. / 3.))),
							        cell->neighbor(i)->info().poreThroatRadius[i]);
						}
					}
				}
			}
		}
		if (debugTPF) {
			cout << endl
			     << "Total nr Throats = " << countTot << "total throats that are too large: " << count
			     << "that is : " << (float(count) * 100.0 / float(countTot)) << "%";
		}
		if ((float(count) / float(countTot)) > 0.1) {
			cout << endl
			     << "Error! Too many pore throats have been adjusted, more than 10%. Simulation is stopped" << count
			     << " tot:" << countTot; /*stopSimulation = true;*/
		}
	}
}


void TwoPhaseFlowEngine::checkVolumeConservationAfterMergingAlgorithm()
{
	//Check volume of the merging of pores, especially required for truncated pore shapes.
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	Real                volumeSingleCells = 0.0, volumeTotal = 0.0, volumeMergedCells = 0.0;
	bool                check = false;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isFictious == 0) {
			volumeTotal = volumeTotal + cell->info().poreBodyVolume;
			if (cell->info().mergedID == 0) { volumeSingleCells = volumeSingleCells + cell->info().poreBodyVolume; }
		}
	}

	for (unsigned int k = 1; k < maxIDMergedCells; k++) {
		check = false;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if ((cell->info().mergedID == k) && (check == false)) {
				volumeMergedCells = volumeMergedCells + cell->info().mergedVolume;
				check = true;
			}
		}
	}
	//if volume is not conserved give error message
	if (math::abs((volumeTotal - volumeMergedCells - volumeSingleCells) / volumeTotal) > 1e-6) {
		std::cerr << endl
		          << "Error! Volume of pores is not conserved between merged pores and total pores: "
		          << "Total pore volume = " << volumeTotal << "Volume of merged cells = " << volumeMergedCells
		          << "Volume of single cells =" << volumeSingleCells;
		stopSimulation = true;
	}
}

void TwoPhaseFlowEngine::calculateResidualSaturation()
{
	//This function computes the entry pressures of the pore throats, as well as the saturation at which an event occurs. This saturation is used as target saturation for determining dt.
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		//Calculate all the pore body radii based on their volumes
		cell->info().poreBodyRadius = getChi(cell->info().numberFacets) * math::pow(cell->info().mergedVolume, (1. / 3.));
		if (cell->info().poreBodyRadius != 0) { cell->info().thresholdPressure = 2.0 * surfaceTension / cell->info().poreBodyRadius; }
		cell->info().thresholdSaturation = 1.0 - (4. / 3.) * 3.14159265359 * math::pow(getChi(cell->info().numberFacets), 3);

		//First check all the macro pores
		if (cell->info().mergedID > 0) {
			for (unsigned int ngb = 0; ngb < 4; ngb++) {
				//Option (1): Throat radii smaller than pore body radii
				if ((cell->info().poreThroatRadius[ngb] < cell->info().poreBodyRadius)
				    && (cell->info().mergedID != cell->neighbor(ngb)->info().mergedID)) {
					cell->info().entryPressure[ngb] = entryMethodCorrection * surfaceTension / cell->info().poreThroatRadius[ngb];
					cell->info().entrySaturation[ngb] = poreSaturationFromPcS(cell, -1.0 * cell->info().entryPressure[ngb]);

					if (cell->info().entrySaturation[ngb] < 0.0 || cell->info().entrySaturation[ngb] > 1.0) {
						cout << endl
						     << "Error With the entrySaturation of a pore throat! " << cell->info().entrySaturation[ngb] << " "
						     << cell->info().poreThroatRadius[ngb] << " and " << cell->info().poreBodyRadius << " "
						     << getKappa(cell->info().numberFacets) << " "
						     << math::log(
						                1.0
						                - (cell->info().poreThroatRadius[ngb] / (cell->info().poreBodyRadius * entryMethodCorrection)))
						     << " CellID=" << cell->info().id << " MergedID =" << cell->info().mergedID
						     << " Facets=" << cell->info().numberFacets;
						cout << endl << "Simulation is terminated because of an error in entry saturation";
						stopSimulation = true;
					}
				}

				//Option (2): Throat radii bigger than pore body radii (this is an error).
				if ((cell->info().poreThroatRadius[ngb] >= cell->info().poreBodyRadius) && (cell->neighbor(ngb)->info().isFictious == 0)
				    && (cell->info().mergedID != cell->neighbor(ngb)->info().mergedID)) {
					cout << endl
					     << "Error, throat radius is larger than the pore body radius for a merged pores: " << cell->info().id
					     << " MergedID" << cell->info().mergedID << " ThroatRadius: " << cell->info().poreThroatRadius[ngb]
					     << " BodyRadius: " << cell->info().poreBodyRadius << " nr facets = " << cell->info().numberFacets;
					cout << endl << "Simulation is terminated because of an pore throat is larger than pore body!";
					stopSimulation = true;
					cell->info().entrySaturation[ngb] = 1.0;
					cell->info().entryPressure[ngb] = 0.0;
				}

				//Option (3): Two tetrahedra from the same pore body, thus an artificial pore throat that is deactivated here.
				if (cell->info().mergedID == cell->neighbor(ngb)->info().mergedID) {
					cell->info().entrySaturation[ngb] = 1.0;
					cell->info().entryPressure[ngb] = 0.0;
				}

				//Option (4): Neighboring Tetrahedron is a boundary cell.
				if (cell->neighbor(ngb)->info().isFictious == true) {
					cell->info().entrySaturation[ngb] = 1.0;
					cell->info().entryPressure[ngb] = 0.0;
				}
			}
		}
		//check all the non-merged pores
		if (cell->info().mergedID == 0) {
			for (unsigned vert = 0; vert < 4; vert++) {
				//Deactive all boundary cells - which are infact individual tetrahedra.
				if (cell->neighbor(vert)->info().isFictious == 1) {
					cell->info().entrySaturation[vert] = 1.0;
					cell->info().entryPressure[vert] = 0.0;
				}
				if (cell->neighbor(vert)->info().isFictious == 0) {
					//calculate the different entry pressures and so on.
					cell->info().entrySaturation[vert]
					        = math::log(
					                  1.0
					                  - (2.0 * cell->info().poreThroatRadius[vert] / (cell->info().poreBodyRadius * entryMethodCorrection)))
					        / (-1.0 * getKappa(cell->info().numberFacets));
					cell->info().entryPressure[vert] = entryMethodCorrection * surfaceTension / cell->info().poreThroatRadius[vert];
					if ((cell->info().entrySaturation[vert] > 1.0 && !cell->info().isFictious)
					    || ((cell->info().entrySaturation[vert] < 0.0))) {
						cout << endl
						     << "entry saturation error!" << cell->info().entrySaturation[vert] << " " << cell->info().id << " "
						     << cell->info().poreBodyRadius << " " << cell->info().poreThroatRadius[vert];
						cell->info().entrySaturation[vert] = 1.0;
						cout << endl << "Simulation is terminated because of an error in entry saturation!";
						stopSimulation = true;
					}
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//....................................Triangulation while maintaining saturation field ............................................//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void TwoPhaseFlowEngine::reTriangulate()
{
	//Governing function to apply triangulation while maintaining saturation distribution.
	if (debugTPF) { std::cerr << endl << "Apply retriangulation"; }
	initializationTriangulation();
	readTriangulation();
	keepTriangulation = false;
	initialization();
	assignWaterVolumesTriangulation();
	actionMergingAlgorithm();
	equalizeSaturationOverMergedCells();
}

void TwoPhaseFlowEngine::initializationTriangulation()
{
	//Resize all relevant functions

	//per sphere
	leftOverVolumePerSphere.resize(scene->bodies->size(), 0);
	untreatedAreaPerSphere.resize(scene->bodies->size(), 0);
	leftOverDVPerSphere.resize(scene->bodies->size(), 0);
	//per tetrahedra
	finishedUpdating.resize(solver->T[solver->currentTes].cellHandles.size(), 0);
	waterVolume.resize(solver->T[solver->currentTes].cellHandles.size(), 0);
	deltaVoidVolume.resize(solver->T[solver->currentTes].cellHandles.size(), 0);
	tetrahedra.resize(solver->T[solver->currentTes].cellHandles.size());
	solidFractionSpPerTet.resize(solver->T[solver->currentTes].cellHandles.size());
	for (unsigned int i = 0; i < solver->T[solver->currentTes].cellHandles.size(); i++) {
		tetrahedra[i].resize(4, 0);
		solidFractionSpPerTet[i].resize(4, 0);
	}
}

void TwoPhaseFlowEngine::readTriangulation()
{
	//Read all relevant information from old assembly of tetrahedra
	for (unsigned int i = 0; i < scene->bodies->size(); i++) {
		untreatedAreaPerSphere[i] = 0.0;
		leftOverVolumePerSphere[i] = 0.0;
		leftOverDVPerSphere[i] = 0.0;
	}

	for (unsigned int i = 0; i < solver->T[solver->currentTes].cellHandles.size(); i++) {
		tetrahedra[i][0] = tetrahedra[i][1] = tetrahedra[i][2] = tetrahedra[i][3] = 1e6;
		solidFractionSpPerTet[i][0] = solidFractionSpPerTet[i][1] = solidFractionSpPerTet[i][2] = solidFractionSpPerTet[i][3] = 0.0;
		waterVolume[i] = 0.0;
		deltaVoidVolume[i] = 0.0;
		finishedUpdating[i] = 0;
	}

	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		waterVolume[cell->info().id] = cell->info().saturation * cell->info().poreBodyVolume;
		deltaVoidVolume[cell->info().id] = cell->info().dv();
		if (cell->info().isFictious) { finishedUpdating[cell->info().id] = -1; }
		if (!cell->info().isFictious) {
			std::pair<int, Real> pairs[4];
			for (unsigned int i = 0; i < 4; i++) {
				pairs[i] = std::make_pair(cell->vertex(i)->info().id(), math::abs(solver->fractionalSolidArea(cell, i)));
			}

			sort(std::begin(pairs), std::end(pairs));
			for (unsigned int j = 0; j < 4; j++) {
				tetrahedra[cell->info().id][j] = pairs[j].first;
				solidFractionSpPerTet[cell->info().id][j] = pairs[j].second;
			}
		}
	}
}

void TwoPhaseFlowEngine::assignWaterVolumesTriangulation()
{
	//Assign saturation to new assembly of tetrahedra
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	unsigned int        saveID = 1e6;
	static unsigned int index = waterVolume.size();

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious) {
			saveID = 1e6;
			unsigned int vert[4]
			        = { cell->vertex(0)->info().id(), cell->vertex(1)->info().id(), cell->vertex(2)->info().id(), cell->vertex(3)->info().id() };
			std::sort(std::begin(vert), std::end(vert));
			for (unsigned int i = 0; i < index; i++) {
				if (tetrahedra[i][0] == vert[0] && tetrahedra[i][1] == vert[1] && tetrahedra[i][2] == vert[2] && tetrahedra[i][3] == vert[3]) {
					saveID = i;
					break;
				}
			}
			if (saveID != 1e6) {
				cell->info().saturation = waterVolume[saveID] / cell->info().poreBodyVolume;
				cell->info().dv() = deltaVoidVolume[saveID];
				if (cell->info().saturation < 0.0) {
					std::cout << endl
					          << "Negative Sat in subFunction1 :" << cell->info().saturation << " " << waterVolume[saveID] << " "
					          << cell->info().poreBodyVolume;
				}
				finishedUpdating[saveID] = 1;
			}
			if (saveID == 1e6) {
				cell->info().saturation = -1;
				for (unsigned int i = 0; i < 4; i++) {
					untreatedAreaPerSphere[cell->vertex(i)->info().id()] += math::abs(solver->fractionalSolidArea(cell, i));
				}
			}
		}
	}
	for (unsigned int i = 0; i < index; i++) {
		if (finishedUpdating[i] == 0) {
			Real totalArea = solidFractionSpPerTet[i][0] + solidFractionSpPerTet[i][1] + solidFractionSpPerTet[i][2] + solidFractionSpPerTet[i][3];
			for (unsigned int j = 0; j < 4; j++) {
				leftOverVolumePerSphere[tetrahedra[i][j]] += (solidFractionSpPerTet[i][j] / totalArea) * waterVolume[i];
				leftOverDVPerSphere[tetrahedra[i][j]] += (solidFractionSpPerTet[i][j] / totalArea) * deltaVoidVolume[i];
			}
		}
	}

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().saturation == -1) {
			Real vol = 0.0, dv = 0.0;
			for (unsigned int j = 0; j < 4; j++) {
				vol += leftOverVolumePerSphere[cell->vertex(j)->info().id()]
				        * (math::abs(solver->fractionalSolidArea(cell, j)) / untreatedAreaPerSphere[cell->vertex(j)->info().id()]);
				dv += leftOverDVPerSphere[cell->vertex(j)->info().id()]
				        * (math::abs(solver->fractionalSolidArea(cell, j)) / untreatedAreaPerSphere[cell->vertex(j)->info().id()]);
			}
			cell->info().saturation = vol / cell->info().poreBodyVolume;
			cell->info().dv() = dv;
			if (cell->info().saturation < 0.0) {
				std::cout << endl
				          << "Error! Negative Sat in sphere allocation: " << cell->info().saturation << " " << vol << " "
				          << cell->info().poreBodyVolume;
			}
		}
	}
}


void TwoPhaseFlowEngine::equalizeSaturationOverMergedCells()
{
	Real                waterVolume2 = 0.0, volume = 0.0, leftOverVolume = 0.0, workVolume = 0.0;
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();

	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().saturation < 0.0) { std::cout << endl << "Error! Before Starting! Negative sat! ... " << cell->info().saturation; }
	}


	for (unsigned int k = 1; k < maxIDMergedCells; k++) {
		waterVolume2 = 0.0;
		leftOverVolume = 0.0;
		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().mergedID == k) {
				waterVolume2 += cell->info().saturation * cell->info().poreBodyVolume;
				volume = cell->info().mergedVolume;
			}
		}

		if (waterVolume2 > volume) {
			leftOverVolume = waterVolume2 - volume;
			waterVolume2 = volume;
		}
		if (leftOverVolume > 0.0) {
			for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
				if (cell->info().mergedID == k && leftOverVolume > 0.0) {
					for (unsigned int j = 0; j < 4; j++) {
						if (cell->info().mergedID != cell->neighbor(j)->info().mergedID && cell->neighbor(j)->info().saturation < 1.0
						    && !cell->neighbor(j)->info().isFictious && leftOverVolume > 0.0) {
							workVolume = (1.0 - cell->neighbor(j)->info().saturation) * cell->neighbor(j)->info().poreBodyVolume;
							std::cout << endl << workVolume << " " << leftOverVolume << " " << cell->neighbor(j)->info().saturation;
							if (workVolume <= leftOverVolume) {
								leftOverVolume -= workVolume;
								cell->neighbor(j)->info().saturation = 1.0;
								std::cout << "inOne";
							}
							std::cout << endl << workVolume << " " << leftOverVolume << " " << cell->neighbor(j)->info().saturation;
							if (workVolume > leftOverVolume && leftOverVolume > 0.0) {
								cell->neighbor(j)->info().saturation
								        = (leftOverVolume
								           + cell->neighbor(j)->info().saturation * cell->neighbor(j)->info().poreBodyVolume)
								        / cell->neighbor(j)->info().poreBodyVolume;
								leftOverVolume = 0.0;
								std::cout << "inOne";
							}
							std::cout << endl << workVolume << " " << leftOverVolume << " " << cell->neighbor(j)->info().saturation;
						}
					}
				}
			}

			if (leftOverVolume > 0.0) { std::cout << endl << "Error! Left over water volume: " << leftOverVolume; }
		}


		for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
			if (cell->info().mergedID == k) { cell->info().saturation = waterVolume2 / cell->info().mergedVolume; }
		}
	}
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (!cell->info().isFictious && cell->info().saturation > 1.0 && cell->info().mergedID == 0) {
			leftOverVolume = (cell->info().saturation - 1.0) * cell->info().poreBodyVolume;
			cell->info().saturation = 1.0;
			for (unsigned int j = 0; j < 4; j++) {
				if (cell->neighbor(j)->info().saturation < 1.0 && !cell->neighbor(j)->info().isFictious && leftOverVolume > 0.0) {
					workVolume = (1.0 - cell->neighbor(j)->info().saturation) * cell->neighbor(j)->info().poreBodyVolume;
					if (workVolume <= leftOverVolume) {
						leftOverVolume -= workVolume;
						cell->neighbor(j)->info().saturation = 1.0;
						std::cout << "inOne-sec";
					}
					std::cout << endl << " sec-" << workVolume << " " << leftOverVolume << " " << cell->neighbor(j)->info().saturation;
					if (workVolume > leftOverVolume && leftOverVolume > 0.0) {
						cell->neighbor(j)->info().saturation
						        = (leftOverVolume + cell->neighbor(j)->info().saturation * cell->neighbor(j)->info().poreBodyVolume)
						        / cell->neighbor(j)->info().poreBodyVolume;
						leftOverVolume = 0.0;
						std::cout << "inOne-sec";
					}
				}
			}
			if (leftOverVolume > 0.0) { std::cout << "Mass left during remeshing" << leftOverVolume; }
		}
	}

	bool redo = false;
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().saturation > 1.0) { redo = true; }
	}
	if (redo) {
		std::cout << "redo calculation";
		equalizeSaturationOverMergedCells();
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//....................................Set pore network from PUA....................... ............................................//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TwoPhaseFlowEngine::setInitialConditions()
{
	if (debugTPF) { std::cerr << endl << "Set initial condition"; }

	//four possible initial configurations are allowed: primary drainage, primary imbibition, secondary drainage, secondary imbibition
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		//make backup of saturated hydraulic conductivity
		for (unsigned int ngb = 0; ngb < 4; ngb++) {
			cell->info().kNorm2[ngb] = cell->info().kNorm()[ngb];
		}

		cell->info().isFictiousId = -1;
		cell->info().isWRes = false;
		cell->info().isNWRes = false;


		if (cell->info().isFictious) {
			//boundary cells are not used here
			cell->info().p() = 0.0;
			cell->info().saturation = 1.0;
			cell->info().hasInterface = false;
		}
		if (!cell->info().isFictious) {
			//Primary drainage
			if (drainageFirst && primaryTPF) {
				cell->info().p() = -1 * initialPC;
				cell->info().saturation = 1.0;
				cell->info().hasInterface = false;
			}
			//Secondary drainage (using saturation field as input parameter)
			if (drainageFirst && !primaryTPF) {
				cell->info().p() = -1 * initialPC;
				if (cell->info().saturation <= cell->info().thresholdSaturation) {
					cell->info().p() = porePressureFromPcS(cell, cell->info().saturation);
					cell->info().hasInterface = true;
				}
				if (cell->info().saturation > cell->info().thresholdSaturation) {
					cell->info().p() = -1 * initialPC;
					cell->info().saturation = 1.0;
					cell->info().hasInterface = false;
					std::cerr << "Warning: local saturation changed for compatibility of local Pc(S)";
				}
			}
			//Primary imbibition
			if (!drainageFirst && primaryTPF) {
				cell->info().p() = -1 * initialPC;
				cell->info().saturation = poreSaturationFromPcS(cell, -1 * initialPC);
				cell->info().hasInterface
				        = true; //FIXME: hasInterface should be false, but until an imbibition criteria is implementend into solvePressure, this should remain true for testing purposes
			}
			//Secondary imbibition
			if (!drainageFirst && !primaryTPF) {
				cell->info().p() = -1 * initialPC;
				if (cell->info().saturation <= cell->info().thresholdSaturation) {
					cell->info().p() = porePressureFromPcS(cell, cell->info().saturation);
					cell->info().hasInterface = true;
				}
				if (cell->info().saturation > cell->info().thresholdSaturation) {
					cell->info().p() = -1 * initialPC;
					cell->info().saturation = 1.0;
					cell->info().hasInterface = false;
					std::cerr << "Warning: local saturation changed for compatibility of local Pc(S)";
				}
			}
		}
	}
}

void TwoPhaseFlowEngine::transferConditions()
{
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		for (unsigned int ngb = 0; ngb < 4; ngb++) {
			cell->info().kNorm2[ngb] = cell->info().kNorm()[ngb];
		}

		if (cell->info().saturation == 1.0) { cell->info().hasInterface = false; }
		if (cell->info().saturation < 1.0) {
			cell->info().hasInterface = true;
			cell->info().p() = porePressureFromPcS(cell, cell->info().saturation);
		}
	}
}

void TwoPhaseFlowEngine::setBoundaryConditions()
{
	RTriangulation&     tri = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().isFictious) {
			for (unsigned int j = 0; j < 4; j++) {
				for (unsigned int i = 0; i < 6; i++) {
					if (cell->vertex(j)->info().id() == i) {
						cell->info().isFictiousId = i;
						if (bndCondIsPressure[cell->info().isFictiousId] && bndCondIsWaterReservoir[cell->info().isFictiousId]) {
							cell->info().p() = bndCondValue[cell->info().isFictiousId];
							cell->info().isWRes = true;
							waterBoundaryPressure = bndCondValue[cell->info().isFictiousId];
						}
						if (bndCondIsPressure[cell->info().isFictiousId] && !bndCondIsWaterReservoir[cell->info().isFictiousId]) {
							cell->info().p() = bndCondValue[cell->info().isFictiousId];
							cell->info().isNWRes = true;
							airBoundaryPressure = bndCondValue[cell->info().isFictiousId];
						}
					}
				}
			}
		}
	}
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		//set initial interface in simulations

		//Drainage
		if (drainageFirst && cell->info().isNWRes) {
			for (unsigned int i = 0; i < 4; i++) {
				if (!cell->neighbor(i)->info().isFictious) {
					if (!deformation) {
						cell->neighbor(i)->info().hasInterface = true;
						cell->neighbor(i)->info().saturation = poreSaturationFromPcS(cell->neighbor(i), -1 * initialPC);
						cell->neighbor(i)->info().p() = -1 * initialPC;
						cell->neighbor(i)->info().airBC = true;
						if (cell->neighbor(i)->info().saturation != cell->neighbor(i)->info().saturation
						    || cell->neighbor(i)->info().saturation > 1.0 || cell->neighbor(i)->info().saturation < 0.0) {
							std::cout << "Error with initial BC saturation: " << cell->neighbor(i)->info().saturation;
						}
					}
					if (deformation) {
						cell->neighbor(i)->info().hasInterface = true;
						cell->neighbor(i)->info().p() = -1 * initialPC;
						cell->neighbor(i)->info().saturation = poreSaturationFromPcS(cell->neighbor(i), -1 * initialPC);
					}

					//Update properties of other cells in pore unit
					if (cell->neighbor(i)->info().mergedID > 0) {
						for (FiniteCellsIterator Mcell = tri.finite_cells_begin(); Mcell != cellEnd; Mcell++) {
							if (Mcell->info().mergedID == cell->neighbor(i)->info().mergedID) {
								Mcell->info().hasInterface = cell->neighbor(i)->info().hasInterface;
								Mcell->info().saturation = cell->neighbor(i)->info().saturation;
								Mcell->info().p() = cell->neighbor(i)->info().p();
								Mcell->info().isNWRes = cell->neighbor(i)->info().isNWRes;
							}
						}
					}
				}
			}
		}
		//Imbibition                                  FIXME(thomas): Needs to be tested for both rigid packings and deforming packings
		if (!drainageFirst) {
			for (unsigned int i = 0; i < 4; i++) {
				if (cell->info().isNWRes) { cell->neighbor(i)->info().airBC = true; }
				if (cell->info().isWRes) {
					cell->neighbor(i)->info().hasInterface = true;
					cell->neighbor(i)->info().saturation = poreSaturationFromPcS(cell->neighbor(i), -1 * initialPC);
					cell->neighbor(i)->info().p() = -1 * initialPC;
					if (cell->neighbor(i)->info().saturation != cell->neighbor(i)->info().saturation
					    || cell->neighbor(i)->info().saturation > 1.0 || cell->neighbor(i)->info().saturation < 0.0) {
						std::cout << "Error with initial BC saturation: " << cell->neighbor(i)->info().saturation;
					}
					if (cell->neighbor(i)->info().mergedID > 0) {
						for (FiniteCellsIterator Mcell = tri.finite_cells_begin(); Mcell != cellEnd; Mcell++) {
							if (Mcell->info().mergedID == cell->neighbor(i)->info().mergedID) {
								Mcell->info().hasInterface = cell->neighbor(i)->info().hasInterface;
								Mcell->info().saturation = cell->neighbor(i)->info().saturation;
								Mcell->info().p() = cell->neighbor(i)->info().p();
								Mcell->info().isNWRes = cell->neighbor(i)->info().isNWRes;
							}
						}
					}
				}
			}
		}
	}


	if (waterBoundaryPressure == 0.0) { waterBoundaryPressure = -1 * initialPC; }
}

} // namespace yade

#endif //TwoPhaseFLOW
