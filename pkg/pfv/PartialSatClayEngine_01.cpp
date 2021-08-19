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

namespace yade { // Cannot have #include directive inside.
CREATE_LOGGER(PartialSatClayEngine);
YADE_PLUGIN((PartialSatClayEngineT)(PartialSatClayEngine)(PartialSatMat)(PartialSatState)(Ip2_PartialSatMat_PartialSatMat_MindlinPhys));

PartialSatClayEngine::~PartialSatClayEngine() { }

// clang-format off
void PartialSatClayEngine::action()
{
	if (debug) cout << "Entering partialSatEngineAction "<<endl;
	if (partialSatDT != 0) timeStepControl();

	if (!isActivated) return;
	timingDeltas->start();
	if (desiredPorosity != 0) {
		Real actualPorosity = Shop::getPorosityAlt();
		volumeCorrection    = desiredPorosity / actualPorosity;
	}
	if (debug) cout << "about to set positions" << endl;
	setPositionsBuffer(true);
	if (!first and alphaBound >= 0)
		addAlphaToPositionsBuffer(true);
	timingDeltas->checkpoint("Position buffer");
	if (first) {
		if (debug) cout << "about to build triangulation" << endl;
		if (multithread)
			setPositionsBuffer(false);
		buildTriangulation(pZero, *solver);
		//if (debug) cout <<"about to add alphatopositionsbuffer" << endl;
		if (alphaBound >= 0) addAlphaToPositionsBuffer(true);
		if (debug) cout << "about to initializevolumes" << endl;
		initializeVolumes(*solver);
		backgroundSolver    = solver;
		backgroundCompleted = true;
		if (partialSatEngine) {
			cout << "setting initial porosity" << endl;
			if (imageryFilePath.compare("none") == 0)
				setInitialPorosity(*solver);
			else {
				cout << "using imagery to set porosity" << endl;
				setPorosityWithImageryGrid(imageryFilePath, *solver);
			}
			if ( crackCellPoroThreshold>0 ) crackCellsAbovePoroThreshold(*solver); // set initial crack network using porosity threshold from imagery
			// if (blockCellPoroThreshold > 0) blockCellsAbovePoroThreshold(*solver);
			if (mineralPoro > 0) {
				cout << "about to block low poros" << endl;
				blockLowPoroRegions(*solver);
			}
			cout << "initializing saturations" << endl;
			initializeSaturations(*solver);
			solver->computePermeability(); //computing permeability again since they depend on the saturations. From here on, we only compute perm once per triangulation
			if (particleSwelling && volumes) setOriginalParticleValues();
			if (useKeq) computeEquivalentBulkModuli(*solver);
			if (debug) cout << "Particle swelling model active, original particle volumes set" << endl;
		}
	}
#ifdef YADE_OPENMP
	solver->ompThreads = ompThreads > 0 ? ompThreads : omp_get_max_threads();
#endif
	timingDeltas->checkpoint("Triangulating");
	updateVolumes(*solver);
	if (!first && partialSatEngine) {
		if (resetOriginalParticleValues) {
			setOriginalParticleValues();
			resetOriginalParticleValues = false;
		}
		if (forceConfinement) simulateConfinement();
		//initializeSaturations(*solver);
		if (!resetVolumeSolids) updatePorosity(*solver);
		else resetPoresVolumeSolids(*solver);

		// if (resetPorosity) { // incase we want to reach a certain state and then set porosity again...chicken and egg problem
		//         resetVolumeSolids=true;
		//         if (imageryFilePath.compare("none")==0) setInitialPorosity(*solver);
		//         else {
		//                 cout << "using imagery to set porosity" << endl;
		//                 setPorosityWithImageryGrid(imageryFilePath, *solver);
		//         }
		// }
		if (useKeq) computeEquivalentBulkModuli(*solver);
	}
	timingDeltas->checkpoint("Update_Volumes");

	epsVolCumulative += epsVolMax;
	if (partialSatDT == 0) retriangulationLastIter++;
	if (!updateTriangulation) updateTriangulation = // If not already set true by another function of by the user, check conditions
		        (defTolerance > 0 && epsVolCumulative > defTolerance) || (meshUpdateInterval > 0 && retriangulationLastIter > meshUpdateInterval);

	// remesh everytime a bond break occurs (for DFNFlow-JCFPM coupling)
	if (breakControlledRemesh) remeshForFreshlyBrokenBonds();

	///compute flow and and forces here
	if (pressureForce) {
#ifdef LINSOLV
		permUpdateIters++;
		if (fixTriUpdatePermInt > 0 && permUpdateIters >= fixTriUpdatePermInt) updateLinearSystem(*solver);
#endif
		if (partialSatEngine) setCellsDSDP(*solver);
		timeDimension = viscosity / maxDSDPj ;
		if (homogeneousSuctionValue==0) solver->gaussSeidel(partialSatDT == 0 ? scene->dt : solverDT);
		else setHomogeneousSuction(*solver);
		timingDeltas->checkpoint("Factorize + Solve");
		if (partialSatEngine) {
			//initializeSaturations(*solver);
			if (!freezeSaturation) updateSaturation(*solver);
			if (debug) cout << "finished initializing saturations" << endl;
		}
		if (!decoupleForces) solver->computeFacetForcesWithCache();
		if (debug) cout << "finished computing facet forces" << endl;
	}

	if (particleSwelling and (retriangulationLastIter == 1 or partialSatDT != 0)) {
		if (suction) {
			if (fracBasedPointSuctionCalc) computeVertexSphericalArea();
			collectParticleSuction(*solver);
		}
		if (swelling) swellParticles();
	}
	if (debug) cout << "finished collecting suction and swelling " << endl;
	//if (freeSwelling && crackModelActive) determineFracturePaths();
	if (brokenBondsRemoveCapillaryforces) removeForcesOnBrokenBonds();
	timingDeltas->checkpoint("compute_Forces");
	///Application of vicscous forces
	scene->forces.sync();
	timingDeltas->checkpoint("forces.sync()");
	//if (!decoupleForces) computeViscousForces ( *solver );
	timingDeltas->checkpoint("viscous forces");
	if (!decoupleForces) {
		if (partialSatDT != 0) addPermanentForces(*solver);
		else applyForces(*solver);
	}
	timingDeltas->checkpoint("Applying Forces");
	if (debug) cout << "finished computing forces and applying" << endl;
///End compute flow and forces
#ifdef LINSOLV
	int sleeping = 0;
	if (multithread && !first) {
		while (updateTriangulation && !backgroundCompleted) { /*cout<<"sleeping..."<<sleeping++<<endl;*/
			sleeping++;
			boost::this_thread::sleep(boost::posix_time::microseconds(1000));
		}
		if (debug && sleeping)
			cerr << "sleeping..." << sleeping << endl;
		if (updateTriangulation || ((meshUpdateInterval > 0 && ellapsedIter > (0.5 * meshUpdateInterval)) && backgroundCompleted)) {
			if (debug) cerr << "switch flow solver" << endl;
			if (useSolver == 0) LOG_ERROR("background calculations not available for Gauss-Seidel");
			if (!fluxChanged) {
				if (fluidBulkModulus > 0 || doInterpolate || partialSatEngine)
					solver->interpolate(solver->T[solver->currentTes], backgroundSolver->T[backgroundSolver->currentTes]);


				//Copy imposed pressures/flow from the old solver
				backgroundSolver->imposedP = vector<pair<CGT::Point, Real>>(solver->imposedP);
				backgroundSolver->imposedF = vector<pair<CGT::Point, Real>>(solver->imposedF);
				solver                     = backgroundSolver;
			} else {
				fluxChanged = false;
			}

			backgroundSolver = shared_ptr<FlowSolver>(new FlowSolver);
			if (metisForced) {
				backgroundSolver->eSolver.cholmod().nmethods           = 1;
				backgroundSolver->eSolver.cholmod().method[0].ordering = CHOLMOD_METIS;
			}
			backgroundSolver->imposedP = vector<pair<CGT::Point, Real>>(solver->imposedP);
			backgroundSolver->imposedF = vector<pair<CGT::Point, Real>>(solver->imposedF);
			if (debug)
				cerr << "switched" << endl;
			setPositionsBuffer(false); //set "parallel" buffer for background calculation
			backgroundCompleted     = false;
			retriangulationLastIter = ellapsedIter;
			updateTriangulation     = false;
			epsVolCumulative        = 0;
			ellapsedIter            = 0;
			boost::thread workerThread(&PartialSatClayEngine::backgroundAction, this);
			workerThread.detach();
			if (debug) cerr << "backgrounded" << endl;
			initializeVolumes(*solver);
			//computeViscousForces(*solver);
			if (debug) cerr << "volumes initialized" << endl;
		} else {
			if (debug && !backgroundCompleted) cerr << "still computing solver in the background, ellapsedIter=" << ellapsedIter << endl;
			ellapsedIter++;
		}
	} else
#endif
	{
		if (updateTriangulation && !first) {
			if (debug) cout << "building tri" <<endl;
			buildTriangulation(pZero, *solver);
			if (debug) cout << "adding alpha to posbuf" <<endl;
			if (alphaBound >= 0) addAlphaToPositionsBuffer(true);
			if (debug) cout << "initializing volumes" <<endl;
			initializeVolumes(*solver);
			if (debug) cout << "out of init volumes" <<endl;
			//resetVolumeSolids=true; //sets new vSolid and (invVoidVolume factored by porosity)
			//computeViscousForces(*solver);
			updateTriangulation     = false;
			epsVolCumulative        = 0;
			retriangulationLastIter = 0;
			ReTrg++;
		}
	}
	first = false;
	timingDeltas->checkpoint("triangulate + init volumes");
	if (getGasPerm) getGasPermeability();
}


/////// Partial Sat Tools /////////
void PartialSatClayEngine::getGasPermeability()
{
	gas_solver->getGasPerm = true;

	// if (gasPermFirst)
	// {
		buildTriangulation(pZero,*gas_solver,true);
		initializeVolumes(*gas_solver);
		//gasPermFirst = false;
	//}
	//setCellsDSDP(*gas_solver);
	gas_solver->gaussSeidel(scene->dt);
	//gas_solver->gaussSeidel(scene->dt);
	//cout << "going into relative cell vel" << endl;
	gas_solver->averageRelativeCellVelocity();
	//cout << "came out of relative cell vel" << endl;
	//getGasPerm = false
}


void PartialSatClayEngine::setHomogeneousSuction(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size  = Tes.cellHandles.size();
	crackedCellTotal = 0;
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell     = Tes.cellHandles[i];
		cell->info().p() = homogeneousSuctionValue;
	}
}

void PartialSatClayEngine::timeStepControl()
{
	if (((elapsedIters > int(partialSatDT / scene->dt)) and partialSatDT != 0) or first) {
		isActivated = true;
		retriangulationLastIter += elapsedIters;
		elapsedIters = 0;
		if (first) {
			collectedDT = scene->dt;
			solverDT    = scene->dt;

		} else {
			solverDT    = collectedDT;
			collectedDT = 0;
		}
		if (debug)
			cout << "using flowtime step =" << solverDT << endl;
	} else {
		if (partialSatDT != 0) {
			elapsedIters++;
			collectedDT += scene->dt;
			isActivated = true;
		}
		isActivated = emulatingAction ? true : false;
		solverDT    = scene->dt;
	}
}

void PartialSatClayEngine::addPermanentForces(FlowSolver& flow)
{
	//typedef typename Solver::FiniteVerticesIterator FiniteVerticesIterator;

	Solver::FiniteVerticesIterator vertices_end = flow.tesselation().Triangulation().finite_vertices_end();

	for (Solver::FiniteVerticesIterator V_it = flow.tesselation().Triangulation().finite_vertices_begin(); V_it != vertices_end; V_it++) {
		scene->forces.setPermForce(V_it->info().id(), makeVector3r(V_it->info().forces));
	}
}

void PartialSatClayEngine::blockLowPoroRegions(FlowSolver& flow)
{
	int          numClumps = 0;
	Tesselation& Tes       = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	//#pragma omp parallel for
	//cout << "blocking low poro regions" << endl;
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().porosity <= mineralPoro) {
			//cout << "found a cell with mineral poro" << endl;
			std::vector<Body::id_t> clumpIds;
			cell->info().blocked = true;
			cell->info().clumped = true;
			//cout << "adding incident particle ids to clump list" << endl;
			addIncidentParticleIdsToClumpList(cell, clumpIds);
			blockMineralCellRecursion(cell, clumpIds); //now cycle on the neighbors
			//cout << "creating clump" << endl;
			if (clumpIds.size() > 0) {
				this->clump(clumpIds, 0);
				numClumps++;
			}
		}
	}
	cout << "clumps created " << numClumps << endl;
}

void PartialSatClayEngine::blockMineralCellRecursion(CellHandle cell, std::vector<Body::id_t>& clumpIds)
{
	for (int facet = 0; facet < 4; facet++) {
		CellHandle nCell = cell->neighbor(facet);
		if (solver->T[solver->currentTes].Triangulation().is_infinite(nCell))
			continue;
		if (nCell->info().Pcondition)
			continue;
		if (nCell->info().clumped)
			continue;
		if (nCell->info().blocked)
			continue;
		if (nCell->info().porosity > mineralPoro)
			continue;

		nCell->info().blocked = true;
		nCell->info().clumped = true;
		//cout << "adding incident particle ids to clump list" << endl;
		addIncidentParticleIdsToClumpList(nCell, clumpIds);
		blockMineralCellRecursion(nCell, clumpIds);
	}
}

void PartialSatClayEngine::addIncidentParticleIdsToClumpList(CellHandle nCell, std::vector<Body::id_t>& clumpIds)
{
	for (int v = 0; v < 4; v++) {
		//if (cell->vertex(v)->info().isFictious) continue;
		const Body::id_t id = Body::id_t(nCell->vertex(v)->info().id());
		if (std::find(clumpIds.begin(), clumpIds.end(), id) != clumpIds.end())
			continue;
		else
			clumpIds.push_back(id);
	}
}

void PartialSatClayEngine::computeEquivalentBulkModuli(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell                   = Tes.cellHandles[i];
		Real        waterFrac              = (cell->info().sat() * cell->info().porosity) / Kw;
		Real        airFrac                = cell->info().porosity * (1. - cell->info().sat()) / Ka;
		Real        solidFrac              = (1. - cell->info().porosity) / Ks;
		Real        Keq                    = 1 / (waterFrac + airFrac + solidFrac);
		cell->info().equivalentBulkModulus = Keq;
	}
}

void PartialSatClayEngine::resetPoresVolumeSolids(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size  = Tes.cellHandles.size();
	crackedCellTotal = 0;
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell     = Tes.cellHandles[i];
		cell->info().vSolids = cell->info().volume() * (1. - cell->info().porosity);
	}
	resetVolumeSolids = false;
}

void PartialSatClayEngine::updatePorosity(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size  = Tes.cellHandles.size();
	crackedCellTotal = 0;
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().Pcondition) continue; //if (cell->info().isAlpha) continue;
		if (!freezePorosity) {
			// if ((!onlyFractureExposedCracks and cell->info().crack)){ // or cell->info().isExposed) {
			// 	crackedCellTotal++; //, cell->info().porosity=fracPorosity;
			// 	                    //crackCellAbovePoroThreshold(cell);
			// }                           //maxPoroClamp;
			// else {
				Real poro = 1. - cell->info().vSolids / cell->info().volume();
				//cout << "old poro" << cell->info().porosity << "new poro" << poro << endl;
				if (poro < minPoroClamp)
					poro = minPoroClamp;
				if (poro > maxPoroClamp)
					poro = maxPoroClamp;
				if (!freezeSaturation and directlyModifySatFromPoro) { // updatesaturation with respect to volume change
					Real dt2        = partialSatDT == 0 ? scene->dt : solverDT;
					Real volWater_o = (cell->info().volume() - cell->info().dv() * dt2) * cell->info().porosity * cell->info().saturation;
					cell->info().saturation = volWater_o / (poro * cell->info().volume()); // update the saturation with respect to new porosity and volume
				}
				cell->info().porosity = poro;
			// }
		} // if we dont want to modify porosity during genesis, but keep the cell curve params updated between triangulations
		// update the parameters for the unique pcs curve in this cell (keep this for interpolation purposes):
		cell->info().Po = Po * exp(a * (meanInitialPorosity - cell->info().porosity)); // use unique cell initial porosity or overall average porosity (mu)?
		if (cell->info().Po > maxPo) cell->info().Po = maxPo;
		cell->info().lambdao = lmbda * exp(b * (meanInitialPorosity - cell->info().porosity));
		if (cell->info().lambdao < minLambdao) cell->info().lambdao = minLambdao;
	}
}


void PartialSatClayEngine::setPorosityWithImageryGrid(string imageryFilePath2, FlowSolver& flow)
{
	std::ifstream file;
	file.open(imageryFilePath2);
	if (!file) {
		cerr << "Unable to open imagery grid reverting to weibull porosity distribution" << endl;
		setInitialPorosity(flow);
		return;
	}
	std::vector<Vector3r> gridCoords;
	std::vector<Real>     porosities;
	int                   l = 0;
	Real                  x, y, z, porosity2;
	while (file >> x >> y >> z >> porosity2) {
		gridCoords.push_back(Vector3r(x, y, z));
		//gridCoords[l][1] = y;
		//gridCoords[l][2] = z;
		porosities.push_back(porosity2);
		l++;
	}
	cout << "finished creating coords vec and porosities" << endl;
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell      = Tes.cellHandles[i];
		Real        finalDist = 1e10;
		Real        finalPoro = meanInitialPorosity;
		CVector     bc        = flow.cellBarycenter(cell);
		// Real xi = cell->info()[0];
		// Real yi = cell->info()[1];
		// Real zi = cell->info()[2];
		//Vector3r cellCoords(xi,yi,zi);
		Vector3r cellCoords(bc[0], bc[1], bc[2]);
		if (!resetVolumeSolids) {
			for (unsigned int k = 0; k < gridCoords.size(); k++) {
				Vector3r vec = cellCoords - gridCoords[k];
				if (vec.norm() < finalDist) {
					finalDist = vec.norm();
					finalPoro = porosities[k];
				}
			}

			// finalPoro = meanInitialPorosity;
			if (finalPoro <= minPoroClamp)
				finalPoro = minPoroClamp;
			if (finalPoro >= maxPoroClamp)
				finalPoro = maxPoroClamp;

			cell->info().porosity = cell->info().initialPorosity = finalPoro;
			if (cell->info().Pcondition) {
				//cout << "setting boundary cell porosity to mean initial" << endl;
				cell->info().porosity = cell->info().initialPorosity = meanInitialPorosity;
			}
			if (finalPoro > maxPorosity)
				maxPorosity = finalPoro;
		}

		cell->info().vSolids = cell->info().volume() * (1. - cell->info().porosity);
		if (!resetVolumeSolids) {
			cell->info().Po = Po* exp(a * (meanInitialPorosity - cell->info().porosity)); // use unique cell initial porosity or overall average porosity (mu)?
			cell->info().lambdao = lmbda * exp(b * (meanInitialPorosity - cell->info().porosity));
		}
	}
	if (resetVolumeSolids)
		resetVolumeSolids = false;
}

} //namespace yade

// clang-format on
#endif //PartialSat
#endif //FLOW_ENGINE
