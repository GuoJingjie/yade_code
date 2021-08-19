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

// clang-format off
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
