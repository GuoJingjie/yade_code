/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              *
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
// Experimental engine under development
#ifdef FLOW_ENGINE
#ifdef PARTIALSAT
#include "PartialSatClayEngine.hpp"
#include <lib/high-precision/Constants.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/CohesiveFrictionalContactLaw.hpp>
#include <pkg/dem/HertzMindlin.hpp>
#include <boost/range/algorithm_ext/erase.hpp>

#ifdef YADE_VTK

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wsuggest-override"
//#include<vtkSmartPointer.h>
#include <lib/compatibility/VTKCompatibility.hpp>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
//#include<vtkDoubleArray.h>
//#include<vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#pragma GCC diagnostic pop

#endif

namespace yade { // Cannot have #include directive inside.

// clang-format off
void PartialSatClayEngine::computeFracturePerm(RTriangulation::Facet_circulator& facet, Real aperture, RTriangulation::Finite_edges_iterator& ed_it,const Real openingPressure,bool gasPermFlag, FlowSolver& flow)
{
	const RTriangulation::Facet& currentFacet = *facet; /// seems verbose but facet->first was declaring a junk cell and crashing program (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1666339)
	const RTriangulation& Tri   = flow.T[flow.currentTes].Triangulation(); //solver->T[solver->currentTes].Triangulation();
	const CellHandle&     cell1 = currentFacet.first;
	const CellHandle&     cell2 = currentFacet.first->neighbor(facet->second);
	Real fracturePerm;
	if (Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr << "Infinite cell found in trickPermeability, should be handled somehow, maybe" << endl;
	if (cell1->info().initiallyCracked) return;
	if (!gasPermFlag){
		fracturePerm = apertureFactor * pow(aperture, 3.) / (12. * viscosity);
	} else {
		fracturePerm = apertureFactor * pow(aperture, 3.) / (12. * airViscosity);
	}
	const Real entryPressure = waterSurfaceTension / (aperture * apertureFactor);
	bool entered = false;

	if (cell1->info().Pcondition or cell2->info().Pcondition) return;

	const Real localSuction = ( ( pAir-cell1->info().p() ) + ( pAir - cell2->info().p() ) ) / 2;
	const bool cellOpened = cell1->info().opened[currentFacet.second];


	if (useOpeningPressure and (localSuction < openingPressure) and !cellOpened){
		return;
	}


	numCracks += 1;

	if (localSuction < entryPressure){
		entered = true; // should both cells need to pass condition?
		// cell1->info().entered=true;
		// cell2->info().entered=true;
		cell1->info().entry[currentFacet.second] = 1;
		cell2->info().entry[Tri.mirror_index(cell1, currentFacet.second)] = 1;

	}

	//cout << "Opening pressure " << openingPressure << " local suction " << localSuction << " aperture " << aperture << " entered " << (cell1->info().entered or cell1->info().entered) << endl;
	// bool visited = (cell1->info().visited[currentFacet.second] || cell2->info().visited[Tri.mirror_index(cell1, currentFacet.second)]);
	// if (visited) return;
	if ((changeCrackSaturation and !entered and !gasPermFlag) or (entered and gasPermFlag)) {
		cell1->info().crack = 1;
		cell2->info().crack = 1;
		cell1->info().visited[currentFacet.second] = 1;
		cell2->info().visited[Tri.mirror_index(cell1, currentFacet.second)] = 1;
		cell1->info().opened[currentFacet.second] = 1;
		cell2->info().opened[Tri.mirror_index(cell1, currentFacet.second)] = 1;
		//cell1->info().blocked=1;
		//cell2->info().blocked=1;
		//cell1->info().saturation = SrM;
		//cell2->info().saturation = SrM; //set low saturation to keep some minimum cohesion
		cell1->info().kNorm()[currentFacet.second]                         *= permAreaFactor;
		cell2->info().kNorm()[Tri.mirror_index(cell1, currentFacet.second)] *= permAreaFactor;
		sumOfApertures += aperture;
		numCracks += 1;


		// cell1->info().kNorm()[currentFacet.second]                          = manualCrackPerm > 0 ? manualCrackPerm : solver->averageK * 0.01;
		// cell2->info().kNorm()[Tri.mirror_index(cell1, currentFacet.second)] = manualCrackPerm > 0 ? manualCrackPerm : solver->averageK * 0.01;
		//cell1->info().p() =
		//cout << "tricked perm on cell " << cell1->info().id << endl;
		// for (int i=0;i<4;i++){ // block this cell from all neighbors, cracked or not cracked.
		//         cell1->info().kNorm()[i] = 0;
		//         cell1->neighbor(i)->info().kNorm()[Tri.mirror_index(cell1,currentFacet.second)]=0;
		// }
	} else { // only using parallel plate approximation if the crack is saturated
		cell1->info().visited[currentFacet.second] = 1;
		cell2->info().visited[Tri.mirror_index(cell1, currentFacet.second)] = 1;
		// cell1->info().kNorm()[currentFacet.second]                         *= permAreaFactor;
		// cell2->info().kNorm()[Tri.mirror_index(cell1, currentFacet.second)] *= permAreaFactor;


		if (!onlyFractureExposedCracks or (onlyFractureExposedCracks and cell1->info().isExposed)) {
			cell1->info().entry[currentFacet.second] = 1;

			cell1->info().crack = 1;
			cell1->info().kNorm()[currentFacet.second] += fracturePerm;
		} //
		if (!onlyFractureExposedCracks or (onlyFractureExposedCracks and cell2->info().isExposed)) {
			cell2->info().entry[Tri.mirror_index(cell1, currentFacet.second)] = 1;
			cell2->info().crack = 1;
			cell2->info().kNorm()[Tri.mirror_index(cell1, currentFacet.second)] +=fracturePerm;
		}

		sumOfApertures += aperture;
		numCracks += 1;
	}

	//cout << "crack set to true in pore"<<endl;
	// cell2->info().blocked = cell1->info().blocked = cell2->info().Pcondition = cell1->info().Pcondition = false; /// those ones will be included in the flow problem
	Point&  CellCentre1             = cell1->info();
	Point&  CellCentre2             = cell2->info();
	CVector networkFractureLength   = CellCentre1 - CellCentre2;                    /// Trying to get fracture's surface
	Real    networkFractureDistance = sqrt(networkFractureLength.squared_length()); /// Trying to get fracture's surface
	Real    networkFractureArea     = pow(networkFractureDistance, 2);              /// Trying to get fracture's surface
	totalFractureArea += networkFractureArea;                                       /// Trying to get fracture's surface
	// 	cout <<" ------------------ The total surface area up to here is --------------------" << totalFractureArea << endl;
	// 	printFractureTotalArea = totalFractureArea; /// Trying to get fracture's surface
	if (calcCrackArea and !cell1->info().isFictious and !cell1->info().isAlpha) {
		CVector edge  = ed_it->first->vertex(ed_it->second)->point().point() - ed_it->first->vertex(ed_it->third)->point().point();
		CVector unitV = edge * (1. / sqrt(edge.squared_length()));
		Point   p3    = ed_it->first->vertex(ed_it->third)->point().point()
		        + unitV * (cell1->info() - ed_it->first->vertex(ed_it->third)->point().point()) * unitV;
		const CVector crackVector = cross_product(CellCentre1 - p3, CellCentre2 - p3);
		const Real crackVectorLength = std::abs(sqrt(crackVector.squared_length()));
		Real halfCrackArea = crackAreaFactor * 0.5 * crackVectorLength; //
		cell1->info().crackArea += halfCrackArea;
		cell2->info().crackArea += halfCrackArea;
		crack_fabric_vector += makeVector3r(crackVector/crackVectorLength)*halfCrackArea;
		crack_fabric_area += halfCrackArea;
		crackArea += halfCrackArea;
		crackVolume += halfCrackArea * aperture;
		if (entered) waterVolume += halfCrackArea*aperture;
	}
}

void PartialSatClayEngine::circulateFacets(RTriangulation::Finite_edges_iterator& edge, Real aperture, const Real openingPressure,bool gasPermFlag, FlowSolver& flow)
{
	const RTriangulation&            Tri    = flow.T[flow.currentTes].Triangulation(); //solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0 = facet1++;
	computeFracturePerm(facet0, aperture, edge, openingPressure,gasPermFlag,flow);
	while (facet1 != facet0) {
		computeFracturePerm(facet1, aperture, edge, openingPressure,gasPermFlag,flow);
		facet1++;
	}
	/// Needs the fracture surface for this edge?
	// Real edgeArea = solver->T[solver->currentTes].computeVFacetArea(edge); cout<<"edge area="<<edgeArea<<endl;
}

void PartialSatClayEngine::trickPermeability(FlowSolver& flow, bool gasPermFlag)
{
	leakOffRate               = 0;
	const RTriangulation& Tri = flow.T[flow.currentTes].Triangulation();
	//	if (!first) interpolateCrack(solver->T[solver->currentTes],flow->T[flow->currentTes]);

	const shared_ptr<InteractionContainer> interactions                         = scene->interactions;
	const shared_ptr<BodyContainer>& bodies					    = scene->bodies;
	//int                                    numberOfCrackedOrJointedInteractions = 0;
	sumOfApertures            					            = 0.;
	averageAperture                                                             = 0;
	maxAperture                                                                 = 0;
	crackArea                                                                   = 0;
	crackVolume                                                                 = 0;
	waterVolume								    = 0;
	numCracks								    = 0;
	crack_fabric_area							    = 0;
	crack_fabric_vector 							    = Vector3r(0,0,0);
	//Real totalFractureArea=0; /// Trying to get fracture's surface
	// 	const shared_ptr<IGeom>& ig;
	// 	const ScGeom* geom; // = static_cast<ScGeom*>(ig.get());
	FiniteEdgesIterator edge = Tri.finite_edges_begin();
	for (; edge != Tri.finite_edges_end(); ++edge) {
		const VertexInfo&              vi1         = (edge->first)->vertex(edge->second)->info();
		const VertexInfo&              vi2         = (edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction = interactions->find(vi1.id(), vi2.id());

		if (interaction && interaction->isReal()) {
			if (edge->first->info().isFictious or edge->first->info().isAlpha)
				continue; /// avoid trick permeability for fictitious

			if (displacementBasedCracks) {
				//const shared_ptr<Clump> clump=YADE_PTR_CAST<Clump>(clumpBody->shape);
				auto mindlingeom = YADE_PTR_CAST<ScGeom>(interaction->geom);
				//ScGeom* mindlingeom = YADE_CAST<ScGeom*>(interaction->geom.get());
				auto mindlinphys = YADE_PTR_CAST<MindlinPhys>(interaction->phys);
				//MindlinPhys* mindlinphys = YADE_CAST<MindlinPhys*>(interaction->phys.get());
				Real crackAperture;
				if (useForceForCracks) {
					const Real forceN = mindlinphys->normalForce.norm();
					if (forceN != 0) continue;
					const shared_ptr<Body>& b1  = (*bodies)[vi1.id()];
					const shared_ptr<Body>& b2  = (*bodies)[vi2.id()];
					const auto state1 = YADE_PTR_CAST<PartialSatState>(b1->state);
					const auto state2 = YADE_PTR_CAST<PartialSatState>(b2->state);
					const Real separation = (state1->pos - state2->pos).norm();
					crackAperture = -(separation - (mindlingeom->radius1 + mindlingeom->radius2)); // setting it negative because the algorithm expects a negative value to indicate opening
					//cout << "force of 0 found" << "crackaperture" << crackAperture << endl;
				} else {
					crackAperture = mindlingeom->penetrationDepth - mindlinphys->initD; // if negative, it has opened up
				}
				const Real openingPressure = waterSurfaceTension / (-crackAperture/apertureFactor);
					// if inf!
				if (-crackAperture<=0) continue;
				//if ( crackAperture < 0 ) std::cout << "-crackAp" << -mindlingeom->penetrationDepth << std::endl;
				// shared_ptr< ScGeom > mindlingeom = std::dynamic_pointer_cast< ScGeom >(std::make_shared(interaction->geom.get()));
				if (-crackAperture < residualAperture and !mindlinphys->isBroken and !useOpeningPressure) continue;
				if (!useOpeningPressure) mindlinphys->isBroken = true; //now even if the displacement reduces back below residAp, we keep tricking this edge in the future
				circulateFacets(edge, -crackAperture, openingPressure, gasPermFlag, flow);

			} else {
				cout << "cohfrict phys partial sat integration not enabled in this version" << endl;
				return;
				// CohFrictPhys* cohfrictphys = YADE_CAST<CohFrictPhys*>(interaction->phys.get());
				// //shared_ptr< CohFrictPhys > cohfrictphys = std::dynamic_pointer_cast< CohFrictPhys >(std::make_shared(interaction->phys.get()));
				//
				// if (!cohfrictphys->isBroken) continue;
				// Real aperture = (cohfrictphys->crackAperture <= residualAperture)? residualAperture : cohfrictphys->crackAperture;
				// if (aperture > maxAperture) maxAperture = aperture;
				// SumOfApertures += aperture;
				// circulateFacets(edge,aperture);
			}
		}
	}
	averageAperture = sumOfApertures / numCracks; /// DEBUG
	// 	cout << " Average aperture in joint ( -D ) = " << AverageAperture << endl; /// DEBUG
}

void PartialSatClayEngine::determineFracturePaths(FlowSolver& flow)
{
	RTriangulation&     tri     = flow.T[flow.currentTes].Triangulation(); //solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		if (cell->info().Pcondition)
			continue;
		cell->info().isExposed = false;
	}
	totalCracks = 0;
	// add logic for handling alpha cells
	if (alphaBound >= 0) {
		for (FlowSolver::VCellIterator it = solver->alphaBoundingCells.begin(); it != solver->alphaBoundingCells.end(); it++) {
			if ((*it) == NULL)
				continue;
			// exposureRecursion(*it); FIXME: add the correct bndPressure argument for alpha shape
		}
	} else {
		for (int i = 0; i < 6; i++) {
			for (FlowSolver::VCellIterator it = solver->boundingCells[i].begin(); it != solver->boundingCells[i].end(); it++) {
				if ((*it) == NULL)
					continue;
				Real bndPressure = solver->boundary(wallIds[i]).value;
				exposureRecursion(*it, bndPressure);
			}
		}
	}
}

void PartialSatClayEngine::exposureRecursion(CellHandle cell, Real bndPressure)
{
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell)) continue;
		if (nCell->info().Pcondition) continue;
		//         if ( (nCell->info().isFictious) && (!isInvadeBoundary) ) continue;
		if (!nCell->info().crack) continue;
		if (nCell->info().isExposed) continue; // another recursion already found it

		if (cell->info().crackNum == 0) nCell->info().crackNum = ++totalCracks; // enable visualization of discretely connected cracks
		else nCell->info().crackNum = cell->info().crackNum;

		nCell->info().isExposed = true;
		//imposePressureFromId(nCell->info().id,bndPressure); // make this a boundary condition now
		nCell->info().Pcondition = true;
		nCell->info().p()        = bndPressure;

		exposureRecursion(nCell, bndPressure);
	}
}


// FIXME it is copying the entire vector ↓ better to take    const vector<Body::id_t>& ids2
Body::id_t PartialSatClayEngine::clump(vector<Body::id_t> ids2, unsigned int discretization)
{
	// create and add clump itself
	//Scene*            scene(Omega::instance().getScene().get());
	shared_ptr<Body>  clumpBody = shared_ptr<Body>(new Body());
	shared_ptr<Clump> clump     = shared_ptr<Clump>(new Clump());
	clumpBody->shape            = clump;
	clumpBody->setBounded(false);
	scene->bodies->insert(clumpBody);
	// add clump members to the clump
	FOREACH(Body::id_t id, ids2)
	{
		if (Body::byId(id, scene)->isClumpMember()) {                                                 //Check, whether the body is clumpMember
			Clump::del(Body::byId(Body::byId(id, scene)->clumpId, scene), Body::byId(id, scene)); //If so, remove it from there
		}
	};

	FOREACH(Body::id_t id, ids2) Clump::add(clumpBody, Body::byId(id, scene));
	Clump::updateProperties(clumpBody, discretization);
	return clumpBody->getId();
}


bool PartialSatClayEngine::findInscribedRadiusAndLocation(CellHandle& cell, std::vector<Real>& coordAndRad)
{
	//cout << "using least sq to find inscribed radius " << endl;
	const Real      prec = 1e-5;
	Eigen::MatrixXd A(4, 3);
	Eigen::Vector4d b2;
	Eigen::Vector3d x;
	Eigen::Vector4d r;
	//std::vector<Real> r(4);
	//std::vector<Real> coordAndRad(4);
	CVector baryCenter(0, 0, 0); // use cell barycenter as initial guess
	for (int k = 0; k < 4; k++) {
		baryCenter = baryCenter + 0.25 * (cell->vertex(k)->point().point() - CGAL::ORIGIN);
		if (cell->vertex(k)->info().isFictious)
			return 0;
	}
	Real xo, yo, zo;
	int  count = 0;
	Real rMean;
	xo            = baryCenter[0];
	yo            = baryCenter[1];
	zo            = baryCenter[2];
	bool finished = false;
	while (finished == false) {
		count += 1;
		if (count > 1000) {
			cerr << "too many iterations during sphere inscription, quitting" << endl;
			return 0;
		}
		// build A matrix (and part of b2)
		for (int i = 0; i < 4; i++) {
			Real xi, yi, zi;
			xi               = cell->vertex(i)->point().x();
			yi               = cell->vertex(i)->point().y();
			zi               = cell->vertex(i)->point().z();
			A(i, 0)          = xo - cell->vertex(i)->point().x();
			A(i, 1)          = yo - cell->vertex(i)->point().y();
			A(i, 2)          = zo - cell->vertex(i)->point().z();
			const Real sqrdD = pow(xo - xi, 2) + pow(yo - yi, 2) + pow(zo - zi, 2);
			r(i)             = sqrt(sqrdD) - sqrt(cell->vertex(i)->point().weight());
		}
		rMean = r.sum() / 4.;

		// build b2
		for (int i = 0; i < 4; i++) {
			Real xi, yi, zi;
			xi   = cell->vertex(i)->point().x();
			yi   = cell->vertex(i)->point().y();
			zi   = cell->vertex(i)->point().z();
			b2(i) = (pow(rMean + sqrt(cell->vertex(i)->point().weight()), 2.) - (pow(xo - xi, 2.) + pow(yo - yi, 2.) + pow(zo - zi, 2.))) / 2.;
		}

		// use least squares (normal equation) to minimize residuals
		x = (A.transpose() * A).ldlt().solve(A.transpose() * b2);
		// if the values are greater than precision, update the guess and repeat
		if (abs(x(0)) > prec || abs(x(1)) > prec || abs(x(2)) > prec) {
			xo += x(0);
			yo += x(1);
			zo += x(2);
		} else {
			coordAndRad[0] = xo + x(0);
			coordAndRad[1] = yo + x(1);
			coordAndRad[2] = zo + x(2);
			coordAndRad[3] = rMean;
			if (rMean > sqrt(cell->vertex(0)->point().weight()))
				return 0; // inscribed sphere might be excessively large if it is in a flat boundary cell
			finished = true;
		}

	} // end while finished == false

	return 1;
}

void PartialSatClayEngine::insertMicroPores(const Real fracMicroPore)
{
	cout << "Inserting micro pores in " << fracMicroPore << " perc. of existing pores " << endl;
	Eigen::MatrixXd M(6, 6);
	//if (!solver->T[solver->currentTes]){cerr << "No triangulation, not inserting micropores" << endl; return}
	Tesselation& Tes = solver->T[solver->currentTes];
	//cout << "Tes set" << endl;
	const long        size = Tes.cellHandles.size();
	std::vector<bool> visited(size);
	std::vector<int>  poreIndices(int(ceil(fracMicroPore * size)));
	bool              numFound;
// randomly select the pore indices that we will turn into micro pores
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (unsigned int i = 0; i < poreIndices.size(); i++) {
		numFound = false;
		while (!numFound) {
			const long num = rand() % size; // + 1?
			if (!visited[num] && !Tes.cellHandles[num]->info().isFictious) {
				visited[num]   = true;
				poreIndices[i] = num;
				numFound       = true;
			}
		}
	}
	cout << "find inscribed sphere radius" << endl;

	// find inscribed sphere radius in selected pores and add body
	// FIXME How do we deal with inscribed spheres that might be overlapping after inscription?
	std::vector<Real> coordsAndRad(4);
	//#pragma omp parallel for
	for (unsigned int i = 0; i < poreIndices.size(); i++) {
		const int   idx  = poreIndices[i];
		CellHandle& cell = Tes.cellHandles[idx];
		for (int j = 0; j < 4; j++)
			if (cell->neighbor(j)->info().isFictious)
				continue; // avoid inscribing spheres in very flat exterior cells (FIXME: we can avoid this by using a proper alpha shape)
		//if (cell->info().Pcondition) continue;
		bool inscribed = findInscribedRadiusAndLocation(cell, coordsAndRad);
		if (!inscribed)
			continue; // sphere couldn't be inscribed, continue loop
		bool contained = checkSphereContainedInTet(cell, coordsAndRad);
		if (!contained)
			continue;
		//cout << "converting to Vector3r" << endl;
		Vector3r coords;
		coords[0]         = Real(coordsAndRad[0]);
		coords[1]         = Real(coordsAndRad[1]);
		coords[2]         = Real(coordsAndRad[2]);
		const Real radius = coordsAndRad[3];
		//cout << "adding body" << endl;
		shared_ptr<Body> body;
		createSphere(body, coords, radius);
		scene->bodies->insert(body);
	}
}

//bool PartialSatClayEngine::checkSphereContainedInTet(CellHandle& cell,std::vector<Real>& coordsAndRad)
//{
//	Eigen::Vector3d inscSphere(coordsAndRad[0],coordsAndRad[1],coordsAndRad[2]);
//	Eigen::Vector3d cellLoc(cell->info()[0],cell->info()[1],cell->info()[2]);
//	Real radius = coordsAndRad[3];
//	 for ( int i=0; i<4; i++ ) {
//		Eigen::Vector3d neighborCellLoc(cell->neighbor(i)->info()[0],cell->neighbor(i)->info()[1],cell->neighbor(i)->info()[2]);
//		Eigen::Vector3d vertexLoc(cell->vertex(facetVertices[i][0])->point().x(),cell->vertex(facetVertices[i][0])->point().y(),cell->vertex(facetVertices[i][0])->point().z());

//		Eigen::Vector3d Surfk = cellLoc-neighborCellLoc;
//		Eigen::Vector3d SurfkNormed = Surfk.normalized();
//		Eigen::Vector3d branch = vertexLoc - inscSphere;
//		Real distToFacet = branch.dot(SurfkNormed);
//		if (distToFacet<0){
//			cerr<< "sphere center outside tet, skipping insertion"<<endl;
//			return 0;
//		} else if (distToFacet<radius) {
//			cerr << "inscribed sphere too large for tetrahedral, decreasing size from "<< radius <<" to "<<distToFacet<<endl;
//			coordsAndRad[3] = distToFacet;
//			radius = distToFacet;
//		}
//	}
//	return 1;
//}

bool PartialSatClayEngine::checkSphereContainedInTet(CellHandle& cell, std::vector<Real>& coordsAndRad)
{
	Eigen::Vector3d inscSphere(coordsAndRad[0], coordsAndRad[1], coordsAndRad[2]);
	//const Real origRadius = coordsAndRad[3];
	//Eigen::Vector3d neighborCellLoc(cell->neighbor(i)->info()[0],cell->neighbor(i)->info()[1],cell->neighbor(i)->info()[2]);
	//	Eigen::MatrixXd A(3,4);
	//	Eigen::Vector4d x;
	//	Eigen::Vector3d bvec(0,0,0);
	//	Real a,b,c,d;
	Real radius = coordsAndRad[3];
	for (int i = 0; i < 4; i++) {
		// using same logic as above but more explicit
		Eigen::Vector3d nhat(cell->info().facetSurfaces[i][0], cell->info().facetSurfaces[i][1], cell->info().facetSurfaces[i][2]);
		nhat = nhat / sqrt(cell->info().facetSurfaces[i].squared_length());
		Eigen::Vector3d xi(
		        cell->vertex(facetVertices[i][0])->point().x(),
		        cell->vertex(facetVertices[i][0])->point().y(),
		        cell->vertex(facetVertices[i][0])->point().z());
		Real distToFacet  = nhat.dot(inscSphere - xi);
		Real exampleScale = sqrt(cell->vertex(facetVertices[i][0])->point().weight());
		Real scale        = exampleScale * minMicroRadFrac;
		// even more explicit, creating plane out of 3 verticies!
		//		for (int j=0;j<3;j++){
		//			A(j,0) = cell->vertex(facetVertices[i][j])->point().x();
		//			A(j,1) = cell->vertex(facetVertices[i][j])->point().y();
		//			A(j,2) = cell->vertex(facetVertices[i][j])->point().z();
		//			A(j,3) = 1;
		//		}
		if (!(distToFacet >= scale)) {
			cout << "minimum radius size doesn't fit,in tet skipping" << endl;
			return 0;
		}
		//		x = A.colPivHouseholderQr().solve(bvec);
		//		a=x(0);b=x(1);c=x(2);d=x(3);
		//		Real sqrtSum = sqrt(a*a+b*b+c*c);
		//		Real distToFacet = (a*coordsAndRad[0]+b*coordsAndRad[1]+c*coordsAndRad[2]+d)/sqrtSum;

		if (distToFacet < 0) {
			cerr << "sphere center outside tet, skipping insertion" << endl;
			return 0;
		} else if (distToFacet < radius) {
			cerr << "inscribed sphere too large for tetrahedral, decreasing size from " << radius << " to " << distToFacet << endl;
			coordsAndRad[3] = distToFacet; //*(1.-minMicroRadFrac);
			radius          = distToFacet; //*(1.-minMicroRadFrac);
		}                                      //else {
		//cerr << "inscribed sphere too small, skipping insertion, btw rad*minMicro= " << exampleScale*minMicroRadFrac << " while dist to facet = " << distToFacet << " and the logic " << (distToFacet>=scale) << endl;
		//	return 0;
		//}
	}
	return 1;
}

void PartialSatClayEngine::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius)
{
	body                     = shared_ptr<Body>(new Body);
	body->groupMask          = 2;
	PartialSatState*   state = dynamic_cast<PartialSatState*>(body->state.get());
	shared_ptr<Aabb>   aabb(new Aabb);
	shared_ptr<Sphere> iSphere(new Sphere);
	state->blockedDOFs = State::DOF_NONE;
	const Real volume  = 4. / 3. * M_PI * pow(radius, 3.);
	state->mass        = volume * microStructureRho;
	//body->state->inertia	= Vector3r(2.0/5.0*body->state->mass*radius*radius,
	//			2.0/5.0*body->state->mass*radius*radius,
	//			2.0/5.0*body->state->mass*radius*radius);
	state->pos            = position;
	state->volumeOriginal = volume;
	state->radiiOriginal  = radius;
	shared_ptr<CohFrictMat> mat(new CohFrictMat);
	mat->young          = microStructureE;
	mat->poisson        = microStructureNu;
	mat->frictionAngle  = microStructurePhi * Mathr::PI / 180.0; //compactionFrictionDeg * Mathr::PI/180.0;
	mat->normalCohesion = mat->shearCohesion = microStructureAdh;
	aabb->color                              = Vector3r(0, 1, 0);
	iSphere->radius                          = radius;
	//iSphere->color	= Vector3r(0.4,0.1,0.1);
	//iSphere->color           = Vector3r(math::unitRandom(),math::unitRandom(),math::unitRandom());
	//iSphere->color.normalize();
	body->shape    = iSphere;
	body->bound    = aabb;
	body->material = mat;
}

void PartialSatClayEngine::printPorosityToFile(string file)
{
	RTriangulation&     tri     = solver->T[solver->currentTes].Triangulation();
	FiniteCellsIterator cellEnd = tri.finite_cells_end();
	std::ofstream            myfile;
	myfile.open(file + boost::lexical_cast<string>(scene->iter) + ".txt");
	for (FiniteCellsIterator cell = tri.finite_cells_begin(); cell != cellEnd; cell++) {
		myfile << cell->info().id << " " << cell->info().porosity << " " << cell->info().crack << "\n";
	}
	myfile.close();
}

void PartialSatClayEngine::simulateConfinement() // TODO: needs to be updated for alpha boundary usage
{
	RTriangulation&                  Tri    = solver->T[solver->currentTes].Triangulation();
	const shared_ptr<BodyContainer>& bodies = scene->bodies;
	for (int bound = 0; bound < 6; bound++) {
		int& id = *solver->boundsIds[bound];
		//solver->conductionBoundingCells[bound].clear();
		if (id < 0)
			continue;

		VectorCell tmpCells;
		tmpCells.resize(10000);
		VCellIterator cells_it  = tmpCells.begin();
		VCellIterator cells_end = Tri.incident_cells(solver->T[solver->currentTes].vertexHandles[id], cells_it);

		for (VCellIterator it = tmpCells.begin(); it != cells_end; it++) {
			CellHandle& cell = *it;
			for (int v = 0; v < 4; v++) {
				if (!cell->vertex(v)->info().isFictious) {
					const long int          id2 = cell->vertex(v)->info().id();
					const shared_ptr<Body>& b2  = (*bodies)[id2];
					if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2)
						continue;
					//auto* state = b2->state.get();
					b2->setDynamic(false);
				}
			}
		}
	}
	forceConfinement = false;
}

void PartialSatClayEngine::computeVertexSphericalArea() // TODO: update for alpha boundary
{
	Tesselation& Tes = solver->T[solver->currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		//	#else
		if (cell->info().blocked) //(cell->info().isFictious) ||
			continue;

		VertexHandle W[4];
		for (int k = 0; k < 4; k++)
			W[k] = cell->vertex(k);
		if (cell->vertex(0)->info().isFictious)
			cell->info().sphericalVertexSurface[0] = 0;
		else
			cell->info().sphericalVertexSurface[0]
			        = solver->fastSphericalTriangleArea(W[0]->point(), W[1]->point().point(), W[2]->point().point(), W[3]->point().point());
		if (cell->vertex(1)->info().isFictious)
			cell->info().sphericalVertexSurface[1] = 0;
		else
			cell->info().sphericalVertexSurface[1]
			        = solver->fastSphericalTriangleArea(W[1]->point(), W[0]->point().point(), W[2]->point().point(), W[3]->point().point());
		if (cell->vertex(2)->info().isFictious)
			cell->info().sphericalVertexSurface[2] = 0;
		else
			cell->info().sphericalVertexSurface[2]
			        = solver->fastSphericalTriangleArea(W[2]->point(), W[1]->point().point(), W[0]->point().point(), W[3]->point().point());
		if (cell->vertex(3)->info().isFictious)
			cell->info().sphericalVertexSurface[3] = 0;
		else
			cell->info().sphericalVertexSurface[3]
			        = solver->fastSphericalTriangleArea(W[3]->point(), W[1]->point().point(), W[2]->point().point(), W[0]->point().point());
	}
	solver->sphericalVertexAreaCalculated = true;
}



/// UNUSED EXTRAS ////


void PartialSatClayEngine::interpolateCrack(Tesselation& Tes, Tesselation& NewTes)
{
	RTriangulation& Tri = Tes.Triangulation();
//RTriangulation& newTri = NewTes.Triangulation();
//FiniteCellsIterator cellEnd = newTri.finite_cells_end();
#ifdef YADE_OPENMP
	const long size = NewTes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads > 0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& newCell = NewTes.cellHandles[i];
#else
	FOREACH(CellHandle & newCell, NewTes.cellHandles)
	{
#endif
		if (newCell->info().isGhost or newCell->info().isAlpha)
			continue;
		CVector center(0, 0, 0);
		if (newCell->info().fictious() == 0)
			for (int k = 0; k < 4; k++)
				center = center + 0.25 * (Tes.vertex(newCell->vertex(k)->info().id())->point().point() - CGAL::ORIGIN);
		else {
			Real boundPos = 0;
			int  coord    = 0;
			for (int k = 0; k < 4; k++) {
				if (!newCell->vertex(k)->info().isFictious)
					center = center
					        + (1. / (4. - newCell->info().fictious()))
					                * (Tes.vertex(newCell->vertex(k)->info().id())->point().point() - CGAL::ORIGIN);
			}
			for (int k = 0; k < 4; k++) {
				if (newCell->vertex(k)->info().isFictious) {
					coord    = solver->boundary(newCell->vertex(k)->info().id()).coordinate;
					boundPos = solver->boundary(newCell->vertex(k)->info().id()).p[coord];
					center   = CVector(
                                                coord == 0 ? boundPos : center[0], coord == 1 ? boundPos : center[1], coord == 2 ? boundPos : center[2]);
				}
			}
		}
		CellHandle oldCell    = Tri.locate(CGT::Sphere(center[0], center[1], center[2]));
		newCell->info().crack = oldCell->info().crack;
		//		For later commit newCell->info().fractureTip = oldCell->info().fractureTip;
		//		For later commit newCell->info().cellHalfWidth = oldCell->info().cellHalfWidth;

		/// compute leakoff rate by summing the flow through facets abutting non-cracked neighbors
		//		if (oldCell->info().crack && !oldCell->info().fictious()){
		//			Real facetFlowRate=0;
		//			facetFlowRate -= oldCell->info().dv();
		//			for (int k=0; k<4;k++) {
		//				if (!oldCell->neighbor(k)->info().crack){
		//					facetFlowRate = oldCell->info().kNorm()[k]*(oldCell->info().shiftedP()-oldCell->neighbor(k)->info().shiftedP());
		//					leakOffRate += facetFlowRate;
		//				}
		//			}
		//		}
	}
}

// void crackCellAbovePoroThreshold(CellHandle& cell)
// {
//         cell->info().crack = 1;
//         for (int j=0; j<4; j++){
//                 CellHandle& ncell = cell->neighbor(i);
//                 Real fracturePerm = apertureFactor*pow(residualAperture,3.)/(12.*viscosity);
//                  //nCell->info().crack=1;
//                 cell->info().kNorm()[i] += fracturePerm; //
//                 nCell->info().kNorm()[Tri.mirror_index(cell,i)] += fracturePerm;
//         }
// }

// void PartialSatClayEngine::trickPermOnCrackedCells(FlowSolver& flow)
// {
//         Tesselation& Tes = flow.T[flow.currentTes];
//         //	#ifdef YADE_OPENMP
//         const long size = Tes.cellHandles.size();
//         Real fracturePerm = apertureFactor*pow(residualAperture,3.)/(12.*viscosity);
//         //#pragma omp parallel for
//         //cout << "blocking low poro regions" << endl;
//         for (long i=0; i<size; i++){
//                 CellHandle& cell = Tes.cellHandles[i];
//                 if ( cell->info().initiallyCracked ){
//                         for (int j=0; j<4; j++){
//                                 const CellHandle& nCell = cell->neighbor(j);
//                                 if (!changeCrackSaturation or (changeCrackSaturation and cell->info().saturation>=SsM) or nCell->info().isFictious) {
//                                         cell->info().kNorm()[j] = fracturePerm;
//                                         nCell->info().kNorm()[Tes.Triangulation().mirror_index(cell,j)] = fracturePerm;
//                                 } else { // block cracked cell if it isnt saturated
//                                         cell->info().crack=1;
//                                         nCell->info().crack=1;
//                                         cell->info().blocked=1;
//                                         nCell->info().blocked=1;
//                                         cell->info().saturation=0;
//                                         nCell->info().saturation=0;
//                                 }
//                         }
//                 }
//         }
// }


void PartialSatClayEngine::crackCellsAbovePoroThreshold(FlowSolver& flow)
{
        Tesselation& Tes = flow.T[flow.currentTes];
	const RTriangulation& Tri   = flow.T[flow.currentTes].Triangulation();
        //	#ifdef YADE_OPENMP
        const long size = Tes.cellHandles.size();
        //#pragma omp parallel for
        //cout << "blocking low poro regions" << endl;
        for (long i=0; i<size; i++){
                CellHandle& cell = Tes.cellHandles[i];
                if ( ( first and cell->info().porosity > crackCellPoroThreshold ) or ( cell->info().initiallyCracked ) ){
                        cell->info().crack = true; cell->info().initiallyCracked = true;
                        //Real fracturePerm = apertureFactor*pow(residualAperture,3.)/(12.*viscosity);

                        for (int j=0; j<4; j++){
                                const CellHandle& nCell = cell->neighbor(j);
                                if (changeCrackSaturation) { // or (changeCrackSaturation and cell->info().saturation>=SsM) or nCell->info().isFictious) {
					cell->info().kNorm()[j]                          = manualCrackPerm > 0 ? manualCrackPerm : solver->averageK * 0.01;
					nCell->info().kNorm()[Tri.mirror_index(cell, j)] = manualCrackPerm > 0 ? manualCrackPerm : solver->averageK * 0.01;
                                // } else { // block cracked cell if it isnt saturated
                                //         cell->info().crack=1;
                                //         nCell->info().crack=1;
                                //         cell->info().blocked=1;
                                //         nCell->info().blocked=1;
                                //         cell->info().saturation=0;
                                //         nCell->info().saturation=0;
                                // }
				}
                        }
                }
        }
}


/******************** Ip2_PartialSatMat_PartialSatMat_MindlinPhys *******/
CREATE_LOGGER(Ip2_PartialSatMat_PartialSatMat_MindlinPhys);

void Ip2_PartialSatMat_PartialSatMat_MindlinPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return; // no updates of an already existing contact necessary
	shared_ptr<MindlinPhys> contactPhysics(new MindlinPhys());
	interaction->phys = contactPhysics;
	const auto mat1   = YADE_CAST<FrictMat*>(b1.get());
	const auto mat2   = YADE_CAST<FrictMat*>(b2.get());

	/* from interaction physics */
	const Real Ea = mat1->young;
	const Real Eb = mat2->young;
	const Real Va = mat1->poisson;
	const Real Vb = mat2->poisson;
	const Real fa = mat1->frictionAngle;
	const Real fb = mat2->frictionAngle;


	/* from interaction geometry */
	const auto scg = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	const Real Da  = scg->refR1 > 0 ? scg->refR1 : scg->refR2;
	const Real Db  = scg->refR2;
	//Vector3r normal=scg->normal;        //The variable set but not used


	/* calculate stiffness coefficients */
	const Real Ga            = Ea / (2 * (1 + Va));
	const Real Gb            = Eb / (2 * (1 + Vb));
	const Real G             = (Ga + Gb) / 2;                                                           // average of shear modulus
	const Real V             = (Va + Vb) / 2;                                                           // average of poisson's ratio
	const Real E             = Ea * Eb / ((1. - math::pow(Va, 2)) * Eb + (1. - math::pow(Vb, 2)) * Ea); // Young modulus
	const Real R             = Da * Db / (Da + Db);                                                     // equivalent radius
	const Real Rmean         = (Da + Db) / 2.;                                                          // mean radius
	const Real Kno           = 4. / 3. * E * sqrt(R);                                                   // coefficient for normal stiffness
	const Real Kso           = 2 * sqrt(4 * R) * G / (2 - V);                                           // coefficient for shear stiffness
	const Real frictionAngle = (!frictAngle) ? math::min(fa, fb) : (*frictAngle)(mat1->id, mat2->id, mat1->frictionAngle, mat2->frictionAngle);

	const Real Adhesion = 4. * Mathr::PI * R * gamma; // calculate adhesion force as predicted by DMT theory

	/* pass values calculated from above to MindlinPhys */
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	//contactPhysics->prevNormal = scg->normal; // used to compute relative rotation
	contactPhysics->kno           = Kno; // this is just a coeff
	contactPhysics->kso           = Kso; // this is just a coeff
	contactPhysics->adhesionForce = Adhesion;

	contactPhysics->kr        = krot;
	contactPhysics->ktw       = ktwist;
	contactPhysics->maxBendPl = eta * Rmean; // does this make sense? why do we take Rmean?

	/* compute viscous coefficients */
	if (en && betan) throw std::invalid_argument("Ip2_PartialSatMat_PartialSatMat_MindlinPhys: only one of en, betan can be specified.");
	if (es && betas) throw std::invalid_argument("Ip2_PartialSatMat_PartialSatMat_MindlinPhys: only one of es, betas can be specified.");

	// en or es specified, just compute alpha, otherwise alpha remains 0
	if (en || es) {
		const Real logE       = log((*en)(mat1->id, mat2->id));
		contactPhysics->alpha = -sqrt(5 / 6.) * 2 * logE / sqrt(pow(logE, 2) + pow(Mathr::PI, 2))
		        * sqrt(2 * E * sqrt(R)); // (see Tsuji, 1992), also [Antypov2011] eq. 17
	}

	// betan specified, use that value directly; otherwise give zero
	else {
		contactPhysics->betan = betan ? (*betan)(mat1->id, mat2->id) : 0;
		contactPhysics->betas = betas ? (*betas)(mat1->id, mat2->id) : contactPhysics->betan;
	}
}

} //namespace yade

// clang-format on
#endif //PartialSat
#endif //FLOW_ENGINE
