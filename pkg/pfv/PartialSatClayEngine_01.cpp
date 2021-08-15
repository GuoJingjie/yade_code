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


void PartialSatClayEngine::setInitialPorosity(FlowSolver& flow)
{ // assume that the porosity is weibull distributed
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		//if (cell->info().isAlpha) continue;

		//Real voidSpace = cell->info().volume()*meanInitialPorosity; // assume voidSpace is equivalent to mean porosity
		//		Real oneSixth = (1./6.); // speed up calcs
		//unsigned long long int numPoreCounts = ceil(voidSpace/(oneSixth*pow(meanPoreSizeDiameter,3)*M_PI)); // determine number of pores we expect to encounter for this cell
		//cout << "numPoreCounts "<< numPoreCounts << " voidSpace " << voidSpace << " volPore " << (oneSixth*pow(meanPoreSizeDiameter,3)*M_PI) <<  endl;

		// generate a pore size based on exponential distribution for each expected pore (using MIP values to fit proper parameters)
		//		Real cellPoreVol = 0; int numPoreCounts = 1e6;
		//		#pragma omp parallel
		//		for (int j=0;j<numPoreCounts;j++)
		//		{
		//			Real poreDiam = exponentialDeviate(alphaExpRate,betaExpRate);
		//			//cout << "pore diam generated " << poreDiam <<endl;
		//			cellPoreVol += oneSixth * pow(poreDiam*1e-6,3) * M_PI;
		//		}
		//
		//		// determine new porosity value
		//		Real poro = cellPoreVol / cell->info().volume();
		if (!resetVolumeSolids) {
			Real poro;
			if (!constantPorosity) {
				poro = meanInitialPorosity * weibullDeviate(lambdaWeibullShape, kappaWeibullScale);
				if (poro < minPoroClamp)
					poro = minPoroClamp;
				if (poro > maxPoroClamp)
					poro = maxPoroClamp;
			} else {
				poro = meanInitialPorosity;
			}
			cell->info().porosity = cell->info().initialPorosity = poro;
			if (poro > maxPorosity)
				maxPorosity = poro;
		}

		cell->info().vSolids = cell->info().volume() * (1. - cell->info().porosity);
		// cell->info().invVoidVolume() = 1./(cell->info().volume()*cell->info().porosity); // do this at each triangulation?
		// set parameters for unique PcS curve on this cell:
		if (!resetVolumeSolids) {
			cell->info().Po = Po * exp(a * (meanInitialPorosity - cell->info().porosity)); // use unique cell initial porosity or overall average porosity (mu)?
			cell->info().lambdao = lmbda * exp(b * (meanInitialPorosity - cell->info().porosity));
		}
	}
	if (resetVolumeSolids)
		resetVolumeSolids = false;
}

Real PartialSatClayEngine::weibullDeviate(Real lambda, Real k)
{
	std::random_device              rd;
	std::mt19937                    e2(rd());
	std::weibull_distribution<Real> weibullDistribution(lambda, k);
	Real                            correction = weibullDistribution(e2);
	return correction;
}

Real PartialSatClayEngine::exponentialDeviate(Real a2, Real b2)
{
	std::random_device                   dev;
	std::mt19937                         rng(dev());
	std::uniform_real_distribution<Real> dist(0, 1.);
	Real                                 x = dist(rng);
	if (x == 1.)
		return 9.999999999999e-1; // return value to avoid undefined behavior
	Real deviate = -(1. / b2) * (log(1. - x) - log(a2));
	return exp(deviate); // linearized value, so we convert back using exp(y)
}

Real PartialSatClayEngine::laplaceDeviate(Real mu, Real b2)
{
	std::random_device                   dev;
	std::mt19937                         rng(dev());
	std::uniform_real_distribution<Real> dist(-0.5, 0.5);
	Real                                 x   = dist(rng);
	Real                                 sgn = x > 0 ? 1. : -1.;
	return mu - b2 * sgn * log(1. - 2. * fabs(x)); // inverse of laplace CDF
	                                              //	if (x<mu) return  1./2. * exp((x-mu)/b);
	                                              //	else return 1. - 1./2.*exp(-(x-mu)/b);
}

void PartialSatClayEngine::setCellsDSDP(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	maxDSDPj = 0;
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		Real        deriv;
		if (cell->info().isAlpha) continue;
		deriv = dsdp(cell);
		if ( deriv > maxDSDPj ) maxDSDPj = deriv;
		if (freezeSaturation) deriv = 0;
		if (!math::isnan(deriv)) cell->info().dsdp = deriv;
		else cell->info().dsdp = 0;
		cell->info().oldPressure = cell->info().p();
		// use this value to determine critical timestep
		//cout << "dsdp " << deriv << endl;
	}
}

Real PartialSatClayEngine::dsdp(CellHandle& cell)
{
	//	Real pc1 = pAir - cell->info().p();
	//	// derivative estimate
	//	Real saturation1 = pow((1. + pow(pc1/cell->info().Po,1./(1.-cell->info().lambdao))),-cell->info().lambdao);
	//	Real pc2 = pc1+100.; // small increment
	//	Real saturation2 = pow((1. + pow(pc2/cell->info().Po,1./(1.-cell->info().lambdao))),-cell->info().lambdao);
	//	Real dsdp = (saturation1 - saturation2) / (pc1-pc2);
	//	return dsdp;
	// analytical derivative of van genuchten
	const Real pc = pAir - cell->info().p(); // suction
	if (pc <= 0) return 0;
	const Real term1 = pow(pow(pc / cell->info().Po, 1. / (1. - cell->info().lambdao)) + 1., (-cell->info().lambdao - 1.));
	const Real term2 = cell->info().lambdao * pow(pc / cell->info().Po, 1. / (1. - cell->info().lambdao) - 1.);
	const Real term3 = cell->info().Po * (1. - cell->info().lambdao);
	return term1 * term2 / term3; // techncially this derivative should be negative, but we use a van genuchten fit for suction, not water pressure. Within the numerical formulation, we want the change of saturation with respect to water pressure (not suction). Which essentially reverses the sign of the derivative.

	// alternate form of analytical derivative from VG 1908 https://www.nrc.gov/docs/ML0330/ML033070005.pdf
	//	Real term1 = -cell->info().lambdao/(pc*(1.-cell->info().lambdao));
	//	Real term2 = pow(vanGenuchten(cell,pc),1./cell->info().lambdao);
	//	Real term3 = pow(pow(1.-vanGenuchten(cell,pc),1./cell->info().lambdao),cell->info().lambdao);
	//	return term1*term2*term3;
}

Real PartialSatClayEngine::vanGenuchten(CellHandle& cell, Real pc)
{
	return pow((1. + pow(pc / cell->info().Po, 1. / (1. - cell->info().lambdao))), -cell->info().lambdao);
}

void PartialSatClayEngine::initializeSaturations(FlowSolver& flow)
{
	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		//if (cell->info().blocked) continue;
		setSaturationFromPcS(cell);
	}
}

void PartialSatClayEngine::updateBoundarySaturation(FlowSolver& flow)
{
	if (alphaBound >= 0) {
		for (FlowSolver::VCellIterator it = flow.alphaBoundingCells.begin(); it != flow.alphaBoundingCells.end(); it++) {
			if ((*it) == NULL) continue;
			setSaturationFromPcS(*it);
		}
	} else {
		for (int i = 0; i < 6; i++) {
			for (FlowSolver::VCellIterator it = flow.boundingCells[i].begin(); it != flow.boundingCells[i].end(); it++) {
				if ((*it) == NULL) continue;
				setSaturationFromPcS(*it);
			}
		}
	}
}

Real PartialSatClayEngine::getTotalVolume()
{
	Tesselation& Tes = solver->T[solver->currentTes];
	//	#ifdef YADE_OPENMP
	totalSpecimenVolume = 0;
	const long size     = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().isAlpha or cell->info().isFictious)
			continue;
		totalSpecimenVolume += cell->info().volume();
	}
	return totalSpecimenVolume;
}

void PartialSatClayEngine::setSaturationFromPcS(CellHandle& cell)
{
	// using van genuchten
	Real pc = pAir - cell->info().p(); // suction in macrostructure
	Real saturation;
	if (pc >= 0)
		saturation = vanGenuchten(cell, pc); //pow((1. + pow(pc/cell->info().Po,1./(1.-cell->info().lambdao))),-cell->info().lambdao);
	else
		saturation = 1.;
	if (saturation < SrM)
		saturation = SrM;
	if (saturation > SsM)
		saturation = SsM;
	cell->info().saturation        = saturation;
	cell->info().initialSaturation = saturation;
}


void PartialSatClayEngine::updateSaturation(FlowSolver& flow)
{
	// 	Tesselation& Tes = flow.T[flow.currentTes];
	// //	#ifdef YADE_OPENMP
	// 	const long size = Tes.cellHandles.size();
	// //	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	//     	for (long i=0; i<size; i++){
	// 		CellHandle& cell = Tes.cellHandles[i];
	// 		if (cell->info().Pcondition or cell->info().isAlpha) continue;
	// 		Real sum=0;
	// 		for (int j=0; j<4; j++){
	// 			CellHandle neighborCell = cell->neighbor(j);
	// 			//if (neighborCell->info().isAlpha) continue;
	// 			sum += cell->info().kNorm()[j] * (cell->info().p() - neighborCell->info().p());
	// 		}
	//         	cell->info().saturation = cell->info().saturation - scene->dt * cell->info().invVoidVolume() * sum;
	//                 if (cell->info().saturation<1e-6) {
	//                         cell->info().saturation = 1e-6;
	//                         cerr << "cell saturation dropped below threshold" << endl;
	//                 }
	//                 if (cell->info().saturation>1) {
	//                         cell->info().saturation = 1;
	//                         cerr << "cell saturation exceeded 1" << endl;
	//                 }
	// 	}


	Tesselation& Tes = flow.T[flow.currentTes];
	//	#ifdef YADE_OPENMP
	const long size = Tes.cellHandles.size();
	//	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];
		if (cell->info().Pcondition or cell->info().isAlpha or cell->info().blocked)
			continue;
		cell->info().saturation = cell->info().saturation + cell->info().dsdp * (cell->info().p() - cell->info().oldPressure);
		if (cell->info().saturation < SrM)
			cell->info().saturation = SrM;
		if (cell->info().saturation > SsM)
			cell->info().saturation = SsM;

		// if (cell->info().saturation < 1e-5) {
		// 	cell->info().saturation = 1e-5;
		// 	//cerr << "cell saturation dropped below threshold" << endl;
		// } // keep an ultra low minium value to avoid numerical issues?
	}


	//         Tesselation& Tes = flow.T[flow.currentTes];
	// 	const long sizeFacets = Tes.facetCells.size();
	// //	#pragma omp parallel for  //FIXME: does not like running in parallel for some reason
	//     	for (long i=0; i<sizeFacets; i++){
	// 		std::pair<CellHandle,int> facetPair = Tes.facetCells[i];
	// 		const CellHandle& cell = facetPair.first;
	// 		const CellHandle& neighborCell = cell->neighbor(facetPair.second);
	//                 const Real satFlux = cell->info().invVoidVolume()*cell->info().kNorm()[facetPair.second] * (cell->info().p() - neighborCell->info().p());
	//                 if (!cell->info().Pcondition) cell->info().saturation -= satFlux*scene->dt;
	//                 if (!cell->info().Pcondition) neighborCell->info().saturation += satFlux*scene->dt;
	//         }
}

void PartialSatClayEngine::resetParticleSuctions()
{
	YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b2, scene->bodies)
	{

		if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2)
			continue;
		if (!b2->isStandalone())
			continue;

		PartialSatState* state = dynamic_cast<PartialSatState*>(b2->state.get());
		state->suction         = 0;
	}
	YADE_PARALLEL_FOREACH_BODY_END();
}

void PartialSatClayEngine::collectParticleSuction(FlowSolver& flow)
{
	resetParticleSuctions();
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	Tesselation&               Tes    = flow.T[flow.currentTes];
	const long                 size   = Tes.cellHandles.size();
#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
	for (long i = 0; i < size; i++) {
		CellHandle& cell = Tes.cellHandles[i];

		if (cell->info().isGhost || cell->info().blocked || cell->info().isAlpha) //|| cell->info().Pcondition || cell->info().isFictious || cell->info().isAlpha
			continue; // Do we need special cases for fictious cells? May need to consider boundary contriubtion to node saturation...
		for (int v = 0; v < 4; v++) {
			//if (cell->vertex(v)->info().isFictious) continue;
			if (cell->vertex(v)->info().isAlpha) continue; // avoid adding alpha vertex to suction?
			const long int          id = cell->vertex(v)->info().id();
			const shared_ptr<Body>& b2  = (*bodies)[id];
			if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2)
				continue;
			PartialSatState* state  = dynamic_cast<PartialSatState*>(b2->state.get());
			Sphere*          sphere = dynamic_cast<Sphere*>(b2->shape.get());
			//if (cell->info().isExposed) state->suctionSum+= pAir; // use different pressure for exposed cracks?
			if (!fracBasedPointSuctionCalc) {
				state->suctionSum += pAir - cell->info().p();
				state->incidentCells += 1;
			} else {
				Real areaFrac  = cell->info().sphericalVertexSurface[v];
				Real totalArea = 4 * M_PI * pow(sphere->radius, 2);
				state->suction += (areaFrac / totalArea) * (pAir - cell->info().p());
			}
		}
	}
}

void PartialSatClayEngine::setOriginalParticleValues()
{
	YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b2, scene->bodies)
	{

		if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2)
			continue;
		Sphere*          sphere = dynamic_cast<Sphere*>(b2->shape.get());
		const Real       volume = 4. / 3. * M_PI * pow(sphere->radius, 3.);
		PartialSatState* state  = dynamic_cast<PartialSatState*>(b2->state.get());
		state->volumeOriginal   = volume;
		state->radiiOriginal    = sphere->radius;
	}
	YADE_PARALLEL_FOREACH_BODY_END();
}


void PartialSatClayEngine::swellParticles()
{
	const Real                       suction0 = pAir - pZero;
	totalVolChange                            = 0;

	YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b2, scene->bodies)
	{
		if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2) continue;
		if (!b2->isStandalone()) continue;
		Sphere* sphere = dynamic_cast<Sphere*>(b2->shape.get());
		//const Real volume = 4./3. * M_PI*pow(sphere->radius,3);
		PartialSatState* state = dynamic_cast<PartialSatState*>(b2->state.get());

		if (!fracBasedPointSuctionCalc) {
			state->lastIncidentCells = state->incidentCells;
			state->suction           = state->suctionSum / state->incidentCells;
			state->incidentCells     = 0; // reset to 0 for next time step
			state->suctionSum        = 0; //
		}

		const Real volStrain = betam / alpham * (exp(-alpham * state->suction) - exp(-alpham * suction0));
		//		const Real rOrig = pow(state->volumeOriginal * 3. / (4.*M_PI),1./3.);
		//
		const Real vNew = state->volumeOriginal * (volStrain + 1.);
		const Real rNew = pow(3. * vNew / (4. * M_PI), 1. / 3.);
		if (rNew <= minParticleSwellFactor * state->radiiOriginal) continue;  // dont decrease size too much
		totalVolChange += (pow(rNew, 3.) - pow(sphere->radius, 3.)) * 4. / 3. * M_PI;
		state->radiiChange = rNew - state->radiiOriginal;
		sphere->radius     = rNew;
		//		cout << "volStrain "<<volStrain<<" avgSuction "<<avgSuction<<" suction0 " <<suction0<<" rDel "<<rDel<<" rNew "<< rNew << " rOrig "<< rOrig << endl;
	}
	YADE_PARALLEL_FOREACH_BODY_END();
}

void PartialSatClayEngine::artificialParticleSwell(const Real volStrain){

	YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b2, scene->bodies)
	{
		if (b2->shape->getClassIndex() != Sphere::getClassIndexStatic() || !b2) continue;
		if (!b2->isStandalone()) continue;
		Sphere* sphere = dynamic_cast<Sphere*>(b2->shape.get());
		//const Real volume = 4./3. * M_PI*pow(sphere->radius,3);
		PartialSatState* state = dynamic_cast<PartialSatState*>(b2->state.get());

		// if (!fracBasedPointSuctionCalc) {
		// 	state->lastIncidentCells = state->incidentCells;
		// 	state->suction           = state->suctionSum / state->incidentCells;
		// 	state->incidentCells     = 0; // reset to 0 for next time step
		// 	state->suctionSum        = 0; //
		// }

		// const Real volStrain = betam / alpham * (exp(-alpham * state->suction) - exp(-alpham * suction0));
		//		const Real rOrig = pow(state->volumeOriginal * 3. / (4.*M_PI),1./3.);
		//
		const Real vNew = state->volumeOriginal * (volStrain + 1.);
		const Real rNew = pow(3. * vNew / (4. * M_PI), 1. / 3.);
		if (rNew <= minParticleSwellFactor * state->radiiOriginal) continue;  // dont decrease size too much
		// totalVolChange += (pow(rNew, 3.) - pow(sphere->radius, 3.)) * 4. / 3. * M_PI;
		state->radiiChange = rNew - state->radiiOriginal;
		sphere->radius     = rNew;
		//		cout << "volStrain "<<volStrain<<" avgSuction "<<avgSuction<<" suction0 " <<suction0<<" rDel "<<rDel<<" rNew "<< rNew << " rOrig "<< rOrig << endl;
	}
	YADE_PARALLEL_FOREACH_BODY_END();
}

Real PartialSatClayEngine::getEnteredRatio() const
{
	Tesselation&                   Tes   = solver->tesselation(); //flow.T[flow.currentTes];
	Real enteredThroats = 0;
	Real numCrackedThroats = 0;
	//Tesselation& Tes = solver->tesselation(); //flow.T[flow.currentTes];
	const long sizeFacets = Tes.facetCells.size();
	//#pragma omp parallel for
	for (long i = 0; i < sizeFacets; i++) {
		std::pair<CellHandle, int> facetPair = Tes.facetCells[i];
		const CellHandle&          cell      = facetPair.first;
		const CellHandle&          nCell     = cell->neighbor(facetPair.second);
		if (cell->info().entry[facetPair.second] and nCell->info().entry[facetPair.second]) enteredThroats+=1;
		numCrackedThroats+=1;
	}
	if (numCrackedThroats==0) return 0;
	else return enteredThroats/numCrackedThroats;

}

void PartialSatClayEngine::initSolver(FlowSolver& flow)
{
	flow.Vtotalissimo = 0;
	flow.VSolidTot    = 0;
	flow.vPoral       = 0;
	flow.sSolidTot    = 0;
	flow.slipBoundary = slipBoundary;
	flow.kFactor      = permeabilityFactor;
	flow.debugOut     = debug;
	flow.useSolver    = useSolver;
#ifdef LINSOLV
	flow.numSolveThreads     = numSolveThreads;
	flow.numFactorizeThreads = numFactorizeThreads;
#endif
	flow.factorizeOnly         = false;
	flow.meanKStat             = meanKStat;
	flow.viscosity             = viscosity;
	flow.tolerance             = tolerance;
	flow.relax                 = relax;
	flow.clampKValues          = clampKValues;
	flow.maxKdivKmean          = maxKdivKmean;
	flow.minKdivKmean          = minKdivKmean;
	flow.meanKStat             = meanKStat;
	flow.permeabilityMap       = permeabilityMap;
	flow.fluidBulkModulus      = fluidBulkModulus;
	flow.fluidRho              = fluidRho;
	flow.fluidCp               = fluidCp;
	flow.thermalEngine         = thermalEngine;
	flow.multithread           = multithread;
	flow.getCHOLMODPerfTimings = getCHOLMODPerfTimings;
	flow.tesselation().maxId   = -1;
	flow.blockedCells.clear();
	flow.sphericalVertexAreaCalculated = false;
	flow.xMin = 1000.0, flow.xMax = -10000.0, flow.yMin = 1000.0, flow.yMax = -10000.0, flow.zMin = 1000.0, flow.zMax = -10000.0;
	flow.partialSatEngine   = partialSatEngine;
	flow.pAir               = pAir;
	flow.freeSwelling       = freeSwelling;
	flow.matricSuctionRatio = matricSuctionRatio;
	flow.nUnsatPerm         = nUnsatPerm;
	flow.SrM                = SrM;
	flow.SsM                = SsM;
	flow.tesselation().vertexHandles.clear();
	flow.tesselation().vertexHandles.resize(scene->bodies->size() + 6, NULL);
	flow.tesselation().vertexHandles.shrink_to_fit();
	flow.alphaBound          = alphaBound;
	flow.alphaBoundValue     = alphaBoundValue;
	flow.freezePorosity      = freezePorosity;
	flow.useKeq              = useKeq;
	flow.useKozeny           = useKozeny;
	flow.bIntrinsicPerm      = bIntrinsicPerm;
	flow.meanInitialPorosity = meanInitialPorosity;
	flow.freezeSaturation    = freezeSaturation;
	flow.permClamp           = permClamp;
	flow.manualCrackPerm     = manualCrackPerm;
	flow.getGasPerm 	 = 0;
}

void PartialSatClayEngine::buildTriangulation(Real pZero2, Solver& flow,bool oneTes)
{
	//cout << "retriangulation" << endl;
	if (first or oneTes) flow.currentTes = 0;
	else {
		flow.currentTes = !flow.currentTes;
		if (debug) cout << "--------RETRIANGULATION-----------" << endl;
	}
	if (debug) cout << "about to reset network" << endl;
	flow.resetNetwork();
	initSolver(flow);
	if (alphaBound < 0) addBoundary(flow);
	else flow.alphaBoundingCells.clear();
	if (debug) cout << "about to add triangulate" << endl;
	//   std::cout << "About to triangulate, push button...";
	//   std::cin.get();
	triangulate(flow);
	//    std::cout << "just triangulated, size" << sizeof(flow.tesselation().Triangulation()) << " push button..." << endl;
	//    std::cin.get();
	if (debug) cout << endl << "Tesselating------" << endl << endl;
	flow.tesselation().compute();
	if (alphaBound < 0) flow.defineFictiousCells();
	// For faster loops on cells define this vector
	flow.tesselation().cellHandles.clear();
	flow.tesselation().cellHandles.reserve(flow.tesselation().Triangulation().number_of_finite_cells());
	FiniteCellsIterator cell_end = flow.tesselation().Triangulation().finite_cells_end();
	int k = 0;
	for (FiniteCellsIterator cell = flow.tesselation().Triangulation().finite_cells_begin(); cell != cell_end; cell++) {
		if (cell->info().isAlpha) continue; // do we want alpha cells in the general cell handle list?
		flow.tesselation().cellHandles.push_back(cell);
		cell->info().id = k++;
	} //define unique numbering now, corresponds to position in cellHandles
	flow.tesselation().cellHandles.shrink_to_fit();

	if (thermalEngine || partialSatEngine) {
		flow.tesselation().facetCells.clear();
		flow.tesselation().facetCells.reserve(flow.tesselation().Triangulation().number_of_finite_facets());
		for (FiniteCellsIterator cell = flow.tesselation().Triangulation().finite_cells_begin(); cell != cell_end; cell++) {
			for (int i = 0; i < 4; i++) {
				if (cell->info().id < cell->neighbor(i)->info().id) {
					flow.tesselation().facetCells.push_back(std::pair<CellHandle, int>(cell, i));
				}
			}
		}
		flow.tesselation().facetCells.shrink_to_fit();
	}
	flow.displayStatistics();
	if (!blockHook.empty()) {
		LOG_INFO("Running blockHook: " << blockHook);
		pyRunString(blockHook);
	}
	//flow.computePermeability(); // move to after interpolate since perm now depends on saturation, and saturation is interpolated value
	//	std::cout << "computed perm, about tto initialize pressure, push button...";
	//   std::cin.get();
	if (multithread && (fluidBulkModulus > 0 || partialSatEngine))
		initializeVolumes(flow); // needed for multithreaded compressible flow (old site, fixed bug https://bugs.launchpad.net/yade/+bug/1687355)
		                         //	if (crackModelActive) trickPermeability(&flow);
	porosity = flow.vPoralPorosity / flow.vTotalPorosity;

	if (alphaBound < 0) boundaryConditions(flow);
	if (debug) cout << "about to initialize pressure" << endl;
	flow.initializePressure(pZero2);

	if (thermalEngine) {
		//initializeVolumes(flow);
		thermalBoundaryConditions(flow);
		flow.initializeTemperatures(tZero);
		flow.sphericalVertexAreaCalculated = false;
	}

	if (debug) cout << "about to interpolate" << endl;
	if (!first && !multithread && (useSolver == 0 || fluidBulkModulus > 0 || doInterpolate || thermalEngine || partialSatEngine) && !oneTes) {
		flow.interpolate(flow.T[!flow.currentTes], flow.tesselation());

		//if (mineralPoro>0) blockLowPoroRegions(*solver);
	} else if (partialSatEngine && oneTes) { // only for the steps where we are building an checking the gas permeability of the specimen

		flow.interpolate(solver->T[!solver->currentTes],flow.T[flow.currentTes]);
		flow.imposedP = vector<pair<CGT::Point, Real>>(solver->imposedP);
		//cout << "interpolating to gas matrix" << endl;
		flow.getGasPerm = true;

	}
	if (debug) cout << "made it out of interpolate" << endl;
	flow.computePermeability();
	if (crackCellPoroThreshold>0) crackCellsAbovePoroThreshold(*solver); //trickPermOnCrackedCells(*solver);
	if (crackModelActive) trickPermeability(flow,getGasPerm);
	if (freeSwelling && crackModelActive && computeFracturePaths) determineFracturePaths(flow);
	if (blockIsoCells) blockIsolatedCells(*solver);
	if (partialSatEngine) {
		updateBoundarySaturation(flow);
		flow.sphericalVertexAreaCalculated = false;
	}
	if (waveAction)
		flow.applySinusoidalPressure(flow.tesselation().Triangulation(), sineMagnitude, sineAverage, 30);
	else if (boundaryPressure.size() != 0)
		flow.applyUserDefinedPressure(flow.tesselation().Triangulation(), boundaryXPos, boundaryPressure);
	if (normalLubrication || shearLubrication || viscousShear)
		flow.computeEdgesSurfaces();
}

} //namespace yade

// clang-format on
#endif //PartialSat
#endif //FLOW_ENGINE
