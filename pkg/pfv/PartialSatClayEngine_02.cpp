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
//////// Post processing tools //////

// void PartialSatClayEngine::saveFractureNetworkVTK(string fileName) // Superceded by savePermeabilityNetworkVTK
// {
// 	vtkSmartPointer<vtkPoints>    intrCellPos = vtkSmartPointer<vtkPoints>::New();
// 	vtkSmartPointer<vtkCellArray> fracCells   = vtkSmartPointer<vtkCellArray>::New();
//
// 	boost::unordered_map<int, int> cIdVector;
// 	int                            curId = 0;
// 	Tesselation&                   Tes   = solver->tesselation(); //flow.T[flow.currentTes];
// 	const long                     size  = Tes.cellHandles.size();
// 	//#pragma omp parallel for
// 	for (long i = 0; i < size; i++) {
// 		CellHandle& cell = Tes.cellHandles[i];
// 		if (solver->T[solver->currentTes].Triangulation().is_infinite(cell))
// 			continue;
// 		if (cell->info().Pcondition)
// 			continue;
// 		if (cell->info().isFictious)
// 			continue;
// 		Point& p2 = cell->info();
// 		intrCellPos->InsertNextPoint(p2[0], p2[1], p2[2]);
// 		cIdVector.insert(std::pair<int, int>(cell->info().id, curId));
// 		curId++;
// 	}
//
// 	//Tesselation& Tes = solver->tesselation(); //flow.T[flow.currentTes];
// 	const long sizeFacets = Tes.facetCells.size();
// 	//#pragma omp parallel for
// 	for (long i = 0; i < sizeFacets; i++) {
// 		std::pair<CellHandle, int> facetPair = Tes.facetCells[i];
// 		const CellHandle&          cell      = facetPair.first;
// 		const CellHandle&          nCell     = cell->neighbor(facetPair.second);
// 		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell) or solver->T[solver->currentTes].Triangulation().is_infinite(cell))
// 			continue;
// 		if (nCell->info().Pcondition or cell->info().Pcondition) continue;
// 		if (nCell->info().isFictious or cell->info().isFictious) continue;
// 		if (cell->info().crack and nCell->info().crack) {
// 			const auto iterId1    = cIdVector.find(cell->info().id);
// 			const auto iterId2    = cIdVector.find(nCell->info().id);
// 			const auto setId1Line = iterId1->second;
// 			const auto setId2Line = iterId2->second;
//
// 			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
// 			line->GetPointIds()->SetId(0, setId1Line);
// 			line->GetPointIds()->SetId(1, setId2Line);
// 			fracCells->InsertNextCell(line);
// 		}
// 	}
// 	vtkSmartPointer<vtkPolyData> intrPd = vtkSmartPointer<vtkPolyData>::New();
// 	intrPd->SetPoints(intrCellPos);
// 	intrPd->SetLines(fracCells);
//
// 	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
// 	//if(compress) writer->SetCompressor(compressor);
// 	//if(ascii) writer->SetDataModeToAscii();
// 	string fn = fileName + "fracNet." + boost::lexical_cast<string>(scene->iter) + ".vtp";
// 	writer->SetFileName(fn.c_str());
// 	writer->SetInputData(intrPd);
// 	writer->Write();
// }

#ifdef YADE_VTK
void PartialSatClayEngine::savePermeabilityNetworkVTK(string fileName)
{
	vtkSmartPointer<vtkPoints>              intrCellPos = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArrayFromReal> permArray   = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	vtkSmartPointer<vtkDoubleArrayFromReal> fracArray   = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	vtkSmartPointer<vtkDoubleArrayFromReal> enteredArray= vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	permArray->SetNumberOfComponents(1);
	permArray->SetName("permeability");
	fracArray->SetName("fractured");
	enteredArray->SetName("entered");
	vtkSmartPointer<vtkCellArray> permCells = vtkSmartPointer<vtkCellArray>::New();

	boost::unordered_map<int, int> cIdVector;
	int                            curId = 0;
	Tesselation&                   Tes   = solver->tesselation(); //flow.T[flow.currentTes];
	const long                     size  = Tes.cellHandles.size();
	//#pragma omp parallel for
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (solver->T[solver->currentTes].Triangulation().is_infinite(cell)) continue;
		//if (cell->info().Pcondition) continue;
		// if (cell->info().isFictious) continue;
		Point& p2 = cell->info();
		intrCellPos->InsertNextPoint(p2[0], p2[1], p2[2]);
		cIdVector.insert(std::pair<int, int>(cell->info().id, curId));
		curId++;
	}

	//Tesselation& Tes = solver->tesselation(); //flow.T[flow.currentTes];
	const long sizeFacets = Tes.facetCells.size();
	//#pragma omp parallel for
	for (long i = 0; i < sizeFacets; i++) {
		std::pair<CellHandle, int> facetPair = Tes.facetCells[i];
		const CellHandle&          cell      = facetPair.first;
		const CellHandle&          nCell     = cell->neighbor(facetPair.second);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell) or solver->T[solver->currentTes].Triangulation().is_infinite(cell))
			continue;
		// if (nCell->info().Pcondition or cell->info().Pcondition)
		// 	continue;
		// if (nCell->info().isFictious or cell->info().isFictious)
		// 	continue;

		const auto iterId1    = cIdVector.find(cell->info().id);
		const auto iterId2    = cIdVector.find(nCell->info().id);
		const auto setId1Line = iterId1->second;
		const auto setId2Line = iterId2->second;

		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		line->GetPointIds()->SetId(0, setId1Line);
		line->GetPointIds()->SetId(1, setId2Line);
		permCells->InsertNextCell(line);
		permArray->InsertNextValue(cell->info().kNorm()[facetPair.second]);
		if (cell->info().opened[facetPair.second] and nCell->info().opened[facetPair.second]) fracArray->InsertNextValue(1);
		else fracArray->InsertNextValue(0);
		if (cell->info().entry[facetPair.second] and nCell->info().entry[facetPair.second]) enteredArray->InsertNextValue(1);
		else enteredArray->InsertNextValue(0);
	}
	vtkSmartPointer<vtkPolyData> intrPd = vtkSmartPointer<vtkPolyData>::New();
	intrPd->SetPoints(intrCellPos);
	intrPd->SetLines(permCells);
	intrPd->GetCellData()->AddArray(permArray);
	intrPd->GetCellData()->AddArray(fracArray);
	intrPd->GetCellData()->AddArray(enteredArray);

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	//if(compress) writer->SetCompressor(compressor);
	//if(ascii) writer->SetDataModeToAscii();
	string fn = fileName + "permNet." + boost::lexical_cast<string>(scene->iter) + ".vtp";
	writer->SetFileName(fn.c_str());
	writer->SetInputData(intrPd);
	writer->Write();
}


void PartialSatClayEngine::saveUnsatVtk(const char* folder, bool withBoundaries)
{
	vector<int>
	            allIds; //an ordered list of cell ids (from begin() to end(), for vtk table lookup), some ids will appear multiple times since boundary cells are splitted into multiple tetrahedra
	vector<int> fictiousN;
	bool        initNoCache = solver->noCache;
	solver->noCache         = false;

	static unsigned int number = 0;
	char                filename[250];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	sprintf(filename, "%s/out_%d.vtk", folder, number++);
	basicVTKwritter vtkfile(0, 0);
	solver->saveMesh(vtkfile, withBoundaries, allIds, fictiousN, filename);
	solver->noCache = initNoCache;

	vtkfile.begin_data("Porosity", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().porosity);
	vtkfile.end_data();

	vtkfile.begin_data("Saturation", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().saturation);
	vtkfile.end_data();

	vtkfile.begin_data("Alpha", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().isAlpha);
	vtkfile.end_data();

	vtkfile.begin_data("Pressure", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().p());
	vtkfile.end_data();

	vtkfile.begin_data("fictious", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(fictiousN[kk]);
	vtkfile.end_data();

	vtkfile.begin_data("blocked", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().blocked);
	vtkfile.end_data();

	vtkfile.begin_data("id", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(allIds[kk]);
	vtkfile.end_data();

	vtkfile.begin_data("fracturedCells", CELL_DATA, SCALARS, INT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().crack);
	vtkfile.end_data();

	vtkfile.begin_data("porosityChange", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(
		        solver->tesselation().cellHandles[allIds[kk]]->info().porosity - solver->tesselation().cellHandles[allIds[kk]]->info().initialPorosity);
	vtkfile.end_data();

	vtkfile.begin_data("saturationChange", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(
		        solver->tesselation().cellHandles[allIds[kk]]->info().saturation
		        - solver->tesselation().cellHandles[allIds[kk]]->info().initialSaturation);
	vtkfile.end_data();


	vtkfile.begin_data("Permeability", CELL_DATA, SCALARS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++) {
		std::vector<Real> perm    = solver->tesselation().cellHandles[allIds[kk]]->info().kNorm();
		Real              permSum = 0;
		for (unsigned int i = 0; i < perm.size(); i++)
			permSum += perm[i];
		vtkfile.write_data(permSum / perm.size());
	}
	vtkfile.end_data();

	solver->averageRelativeCellVelocity();
	vtkfile.begin_data("Velocity", CELL_DATA, VECTORS, FLOAT);
	for (unsigned kk = 0; kk < allIds.size(); kk++)
		vtkfile.write_data(
		        solver->tesselation().cellHandles[allIds[kk]]->info().averageVelocity()[0],
		        solver->tesselation().cellHandles[allIds[kk]]->info().averageVelocity()[1],
		        solver->tesselation().cellHandles[allIds[kk]]->info().averageVelocity()[2]);
	vtkfile.end_data();

#define SAVE_CELL_INFO(INFO)                                                                                                                                   \
	vtkfile.begin_data(#INFO, CELL_DATA, SCALARS, FLOAT);                                                                                                  \
	for (unsigned kk = 0; kk < allIds.size(); kk++)                                                                                                        \
		vtkfile.write_data(solver->tesselation().cellHandles[allIds[kk]]->info().INFO);                                                                \
	vtkfile.end_data();
	//SAVE_CELL_INFO(saturation)
	SAVE_CELL_INFO(Po)
	SAVE_CELL_INFO(lambdao)
	SAVE_CELL_INFO(Pcondition)
	SAVE_CELL_INFO(isExposed)
	SAVE_CELL_INFO(crackNum)
	//	SAVE_CELL_INFO(porosity)
}
#endif

void PartialSatClayEngine::initializeVolumes(FlowSolver& flow)
{
	totalSpecimenVolume = 0;
	//typedef typename Solver::FiniteVerticesIterator FiniteVerticesIterator;

	Solver::FiniteVerticesIterator vertices_end = flow.tesselation().Triangulation().finite_vertices_end();
	CGT::CVector           Zero(0, 0, 0);
	for (Solver::FiniteVerticesIterator V_it = flow.tesselation().Triangulation().finite_vertices_begin(); V_it != vertices_end; V_it++)
		V_it->info().forces = Zero;
	//cout << "in initialize volumes, after force reset" <<endl;
	FOREACH(CellHandle & cell, flow.tesselation().cellHandles)
	{
		//if (cell->info().isAlpha) continue; // not concerned about volumes for alpha cells?
		switch (cell->info().fictious()) {
			case (0): cell->info().volume() = volumeCell(cell); break;
			case (1): cell->info().volume() = volumeCellSingleFictious(cell); break;
			case (2): cell->info().volume() = volumeCellDoubleFictious(cell); break;
			case (3): cell->info().volume() = volumeCellTripleFictious(cell); break;
			default: break;
		}
		//cout << "made it past the volume calc" <<endl;
		if (flow.fluidBulkModulus > 0 || thermalEngine || iniVoidVolumes) {
			cell->info().invVoidVolume() = 1 / (std::abs(cell->info().volume()) - volumeCorrection * flow.volumeSolidPore(cell));
		} else if (partialSatEngine) {
			if (cell->info().volume() <= 0 and cell->info().isAlpha and debug)
				cerr << "cell volume zero, bound to be issues" << endl;
			cell->info().invVoidVolume() = 1 / std::abs(cell->info().volume());
			//cell->info().invVoidVolume() = 1. / std::max(minCellVol, math::abs(cell->info().volume())); // - flow.volumeSolidPore(cell)));
			//if (cell->info().invVoidVolume() == 1. / minCellVol) {
			//	cell->info().blocked = 1;
			//	cout << "using minCellVolume, blocking cell" << endl;
			//}
		}
		//cout << "made it past the invvoidvolume assignment" <<endl;
		if (!cell->info().isAlpha and !cell->info().isFictious)
			totalSpecimenVolume += cell->info().volume();
	}
	if (debug) cout << "Volumes initialised." << endl;
}

void PartialSatClayEngine::updateVolumes(FlowSolver& flow)
{
	if (debug) cout << "Updating volumes.............." << endl;
	Real invDeltaT      = 1 / (partialSatDT == 0 ? scene->dt : solverDT);
	epsVolMax           = 0;
	Real totVol         = 0;
	Real totDVol        = 0;
	totalSpecimenVolume = 0;
#ifdef YADE_OPENMP
	const long size = flow.tesselation().cellHandles.size();
#pragma omp parallel for num_threads(ompThreads > 0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = flow.tesselation().cellHandles[i];
#else
	FOREACH(CellHandle & cell, flow.tesselation().cellHandles)
	{
#endif
		Real newVol, dVol;
		//if (cell->info().isAlpha) continue; // dont care about volumes of alpha cells?
		switch (cell->info().fictious()) {
			case (3): newVol = volumeCellTripleFictious(cell); break;
			case (2): newVol = volumeCellDoubleFictious(cell); break;
			case (1): newVol = volumeCellSingleFictious(cell); break;
			case (0): newVol = volumeCell(cell); break;
			default: newVol = 0; break;
		}
		dVol = cell->info().volumeSign * (newVol - cell->info().volume());
		if (!thermalEngine) cell->info().dv() = dVol * invDeltaT;
		else cell->info().dv() += dVol * invDeltaT; // thermalEngine resets dv() to zero and starts adding to it before this.
		cell->info().volume() = newVol;
		if (!cell->info().isFictious and !cell->info().isAlpha)
			totalSpecimenVolume += newVol;
		if (defTolerance > 0) { //if the criterion is not used, then we skip these updates a save a LOT of time when Nthreads > 1
#pragma omp atomic
			totVol += cell->info().volumeSign * newVol;
#pragma omp atomic
			totDVol += dVol;
		}
	}
	if (defTolerance > 0) epsVolMax = totDVol / totVol;
	//FIXME: move this loop to FlowBoundingSphere
	for (unsigned int n = 0; n < flow.imposedF.size(); n++) {
		flow.IFCells[n]->info().dv() += flow.imposedF[n].second;
		flow.IFCells[n]->info().Pcondition = false;
	}
	if (debug) cout << "Updated volumes, total =" << totVol << ", dVol=" << totDVol << endl;
}


/////// Discrete Fracture Network Functionality ////////

void PartialSatClayEngine::blockCellsAbovePoroThreshold(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	//#pragma omp parallel for
	//cout << "blocking low poro regions" << endl;
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().porosity > crackCellPoroThreshold) {
			cell->info().blocked = true;
			for (int j = 0; j < 4; j++) {
				const CellHandle& nCell = cell->neighbor(j);
				nCell->info().blocked   = true;
			}
		}
	}
}

// void PartialSatClayEngine::blockIsolatedCells(FlowSolver& flow)
// {
//         Tesselation& Tes = flow.T[flow.currentTes];
//         //	#ifdef YADE_OPENMP
//         const long size = Tes.cellHandles.size();
//         //#pragma omp parallel for
//         //cout << "blocking low poro regions" << endl;
//         for (long i=0; i<size; i++){
//                 CellHandle& cell = Tes.cellHandles[i];
//                 if (cell->info().blocked) continue;
//                 for (int j=0; j<4; j++){
//                         const CellHandle& nCell = cell->neighbor(j);
//                         if (!nCell->info().blocked) break;
//                         nCell->info().blocked=true; //cell is surrounded by blocked cells, and therefore needs to be blocked itself.
//                 }
//         }
//
// }

void PartialSatClayEngine::blockIsolatedCells(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	//#pragma omp parallel for
	//cout << "blocking low poro regions" << endl;
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().blocked)
			continue;
		int numBlocked = 0;
		for (int j = 0; j < 4; j++) {
			const CellHandle& nCell = cell->neighbor(j);
			if (nCell->info().blocked)
				numBlocked++;
		}
		if (numBlocked == 4)
			cell->info().blocked = true;
		cell->info().Pcondition = false;
	}
}

void PartialSatClayEngine::removeForcesOnBrokenBonds()
{
	const RTriangulation&                  Tri          = solver->T[solver->currentTes].Triangulation();
	const shared_ptr<InteractionContainer> interactions = scene->interactions;
	FiniteEdgesIterator                    edge         = Tri.finite_edges_begin();
	for (; edge != Tri.finite_edges_end(); ++edge) {
		const VertexInfo&              vi1         = (edge->first)->vertex(edge->second)->info();
		const VertexInfo&              vi2         = (edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction = interactions->find(vi1.id(), vi2.id());

		if (interaction && interaction->isReal()) {
			if (edge->first->info().isFictious) continue;
			auto mindlinphys = YADE_PTR_CAST<MindlinPhys>(interaction->phys);
			if (!mindlinphys->isBroken) continue;
			circulateFacetstoRemoveForces(edge);
		}
	}
}

void PartialSatClayEngine::circulateFacetstoRemoveForces(RTriangulation::Finite_edges_iterator& edge)
{
	const RTriangulation&            Tri    = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0 = facet1++;
	removeForceOnVertices(facet0, edge);
	while (facet1 != facet0) {
		removeForceOnVertices(facet1, edge);
		facet1++;
	}
	/// Needs the fracture surface for this edge?
	// Real edgeArea = solver->T[solver->currentTes].computeVFacetArea(edge); cout<<"edge area="<<edgeArea<<endl;
}

void PartialSatClayEngine::removeForceOnVertices(RTriangulation::Facet_circulator& facet, RTriangulation::Finite_edges_iterator& ed_it)
{
	const RTriangulation::Facet& currentFacet
	        = *facet; /// seems verbose but facet->first was declaring a junk cell and crashing program (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1666339)
	//const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	const CellHandle& cell1 = currentFacet.first;
	const CellHandle& cell2 = currentFacet.first->neighbor(facet->second);
	VertexInfo&       vi1   = (ed_it->first)->vertex(ed_it->second)->info();
	VertexInfo&       vi2   = (ed_it->first)->vertex(ed_it->third)->info();

	// compute area
	Point&  CellCentre1 = cell1->info(); /// Trying to get fracture's surface
	Point&  CellCentre2 = cell2->info(); /// Trying to get fracture's surface
	CVector edge        = ed_it->first->vertex(ed_it->second)->point().point() - ed_it->first->vertex(ed_it->third)->point().point();
	CVector unitV       = edge * (1. / sqrt(edge.squared_length()));
	Point p3 = ed_it->first->vertex(ed_it->third)->point().point() + unitV * (cell1->info() - ed_it->first->vertex(ed_it->third)->point().point()) * unitV;
	Real  halfCrackArea = crackAreaFactor * 0.5 * sqrt(std::abs(cross_product(CellCentre1 - p3, CellCentre2 - p3).squared_length()));

	// modify forces to remove since it is broken
	CVector capillaryForce = edge * halfCrackArea * ((cell1->info().p() + cell2->info().p()) / 2.) * ((cell1->info().sat() + cell2->info().sat()) / 2.);
	//cout << "total force on body"<<vi1.forces[0]<<" "<<vi1.forces[1]<<" "<<vi1.forces[2]<<endl;
	//cout << "capillary force computed" << capillaryForce[0] << " "<<capillaryForce[1]<<" "<<capillaryForce[2]<<endl;
	vi1.forces = vi1.forces + capillaryForce;
	vi2.forces = vi2.forces - capillaryForce;
	//cell1->vertex(facetVertices[j][y])->info().forces = cell1->vertex(facetVertices[j][y])->info().forces -facetNormal*pAir*crossSections[j][y];
}

} //namespace yade

// clang-format on
#endif //PartialSat
#endif //FLOW_ENGINE
