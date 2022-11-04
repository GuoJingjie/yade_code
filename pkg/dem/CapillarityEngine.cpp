/*************************************************************************
*  Copyright (C) 2013 by Bruno Chareyre    <bruno.chareyre@grenoble-inp.fr>  *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_CGAL

#include <pkg/common/ElastMat.hpp>
#include <pkg/dem/CapillarityEngine.hpp>

#include <lib/base/Math.hpp>
#include <lib/high-precision/Constants.hpp>
#include <lib/high-precision/MathApprox.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/HertzMindlin.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <fstream>
#include <iostream>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((CapillarityEngine));

DelaunayInterpolator::Dt CapillarityEngine::dtVbased;
DelaunayInterpolator::Dt CapillarityEngine::dtPbased;

const unsigned DelaunayInterpolator::comb[] = { 1, 2, 3, 0, 1, 2 };

Real CapillarityEngine::intEnergy()
{
	Real energy = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		ScGeom*                currentGeometry = static_cast<ScGeom*>(I->geom.get());
		CapillaryPhysDelaunay* phys            = static_cast<CapillaryPhysDelaunay*>(I->phys.get());
		if (phys) {
			if (phys->SInterface != 0) {
				energy += liquidTension
				        * (phys->SInterface - 4 * Mathr::PI * (pow(currentGeometry->radius1, 2))
				           - 4 * Mathr::PI * (pow(currentGeometry->radius2, 2)));
			}
		}
	}
	return energy;
}


Real CapillarityEngine::wnInterface()
{
	Real wn = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		ScGeom*                currentGeometry = static_cast<ScGeom*>(I->geom.get());
		CapillaryPhysDelaunay* phys            = static_cast<CapillaryPhysDelaunay*>(I->phys.get());
		if (phys) {
			if (phys->SInterface != 0) {
				wn += (phys->SInterface
				       - 2 * Mathr::PI
				               * (pow(currentGeometry->radius1, 2) * (1 + cos(phys->Delta1))
				                  + pow(currentGeometry->radius2, 2) * (1 + cos(phys->Delta2))));
			}
		}
	}
	return wn;
}

Real CapillarityEngine::swInterface()
{
	Real sw = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		ScGeom*                currentGeometry = static_cast<ScGeom*>(I->geom.get());
		CapillaryPhysDelaunay* phys            = static_cast<CapillaryPhysDelaunay*>(I->phys.get());
		if (phys) {
			sw += (2 * Mathr::PI
			       * (pow(currentGeometry->radius1, 2) * (1 - cos(phys->Delta1)) + pow(currentGeometry->radius2, 2) * (1 - cos(phys->Delta2))));
		}
	}
	return sw;
}


Real CapillarityEngine::waterVolume()
{
	Real volume = 0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
	{
		if (!I->isReal()) continue;
		CapillaryPhysDelaunay* phys = static_cast<CapillaryPhysDelaunay*>(I->phys.get());
		if (phys) { volume += phys->vMeniscus; }
	}
	return volume;
}

void CapillarityEngine::triangulateData()
{
	/// We get data from a file and input them in triangulations
	if (solutions.size() > 0) {
		LOG_WARN("CapillarityEngine asking triangulation for the second time. Ignored.");
		return;
	}
	std::ifstream file(inputFilename.c_str());
	if (!file.is_open()) {
		LOG_ERROR("No data file found for capillary law. Check path and inputFilename.");
		return;
	}

	MeniscusPhysicalData dat;
	double               ending;
	while (file.good()) {
		file >> dat.succion >> dat.force >> dat.distance >> dat.volume >> dat.surface >> dat.arcLength >> dat.delta1 >> dat.delta2 >> dat.R >> ending;
		dat.ending = (bool)ending;
		solutions.push_back(dat);
	}
	file.close();
	// Make lists of points with index, so we can use range insertion, more efficient
	// see http://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation_3SettingInformationWhileInserting
	std::vector<std::pair<K::Point_3, unsigned>> pointsP, pointsV;
	for (unsigned int k = 0; k < solutions.size(); k++) {
		pointsP.push_back(std::make_pair(K::Point_3(solutions[k].R, solutions[k].succion, solutions[k].distance), k));
		pointsV.push_back(std::make_pair(K::Point_3(solutions[k].R, solutions[k].volume, solutions[k].distance), k));
	}
	// and now range insertion
	dtPbased.insert(pointsP.begin(), pointsP.end());
	dtVbased.insert(pointsV.begin(), pointsV.end());
}

int firstIteration = 1;
int x              = 0;

template <
        typename IPhysT> // the template parameter should be an IPhys with capillary attributes, see the selection of correct type below, in CapillarityEngine::action()
void CapillarityEngine::solveBridgesT(Real suction, bool reset)
{
	if (!scene) cerr << "scene not defined!";
	shared_ptr<BodyContainer>& bodies = scene->bodies;
	if (dtPbased.number_of_vertices() < 1) triangulateData();
	if (fusionDetection && !bodiesMenisciiList.initialized) bodiesMenisciiList.prepare(scene);
	timingDeltas->checkpoint("init solver");

	const long size = scene->interactions->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& interaction = (*scene->interactions)[i];
		auto                           phys        = static_cast<IPhysT*>(interaction->phys.get());

		if (interaction->isReal() /*&& phys->computeBridge == true FIXME: not available in Mindlin */) {
			const unsigned int& id1 = interaction->getId1();
			const unsigned int& id2 = interaction->getId2();

			/// interaction geometry search (this test is to compute capillarity only between spheres (probably a better way to do that)
			const int& geometryIndex1 = (*bodies)[id1]->shape->getClassIndex(); // !!!
			const int& geometryIndex2 = (*bodies)[id2]->shape->getClassIndex();
			/// definition of interacting objects (not necessarily in contact)
			ScGeom* currentContactGeometry = static_cast<ScGeom*>(interaction->geom.get());

			if (!(geometryIndex1 == geometryIndex2)) {
				if (currentContactGeometry->penetrationDepth < 0) scene->interactions->requestErase(interaction);
				continue;
			}

			/// Interacting Grains:
			// If you want to define a ratio between YADE sphere size and real sphere size
			Real alpha = 1;
			Real R1    = alpha * std::max(currentContactGeometry->radius2, currentContactGeometry->radius1);
			Real R2    = alpha * std::min(currentContactGeometry->radius2, currentContactGeometry->radius1);

			/// intergranular distance
			Real D = -alpha * ((currentContactGeometry->penetrationDepth) + epsilonMean * (R1 + R2));

			bool wasMeniscus = phys->meniscus;
			if (D <= 0 || createDistantMeniscii) {
				if (fusionDetection && !phys->meniscus) bodiesMenisciiList.insert(interaction);
				phys->meniscus = true;
			}
			Real Dinterpol = D / R1;

			/// Suction (Capillary pressure):
			if (imposePressure || (!imposePressure && totalVolumeConstant)) {
				Real Pinterpol          = 0;
				Pinterpol               = phys->isBroken ? 0 : suction * R1 / liquidTension;
				phys->capillaryPressure = suction;
				/// Capillary solution finder:
				if ((Pinterpol >= 0) && phys->meniscus) {
					// NOTE: here, the parameter phys->m will be missing for MindlinCapillaryPhys, or CapillaryPhys
					// Possible solution for merging the capillary flavors: give the same name to CapillaryPhysDelaunay::m  and CapillaryPhys::currentIndexes
					// and make the interpolator a template parameter of this entire function, one interpolator will use a cell handle (from phys->m) or table indices (from CapillaryPhys)
					MeniscusPhysicalData solution = DelaunayInterpolator::interpolate1(
					        dtPbased, K::Point_3(R2 / R1, Pinterpol, Dinterpol), phys->m, solutions, reset);
					/// capillary adhesion force
					Real Finterpol = solution.force;
					phys->fCap     = Finterpol * R1 * liquidTension * currentContactGeometry->normal;
					/// meniscus volume
					phys->vMeniscus = solution.volume * pow(R1, 3);

					if (phys->vMeniscus > 0) {
						phys->meniscus   = true;
						phys->SInterface = solution.surface * pow(R1, 2);
					} else {
						phys->meniscus   = false;
						phys->SInterface = 4 * Mathr::PI * (pow(R1, 2) + pow(R2, 2));
						if (fusionDetection and wasMeniscus) bodiesMenisciiList.remove(interaction);
						if (currentContactGeometry->penetrationDepth < 0) scene->interactions->requestErase(interaction);
						continue;
					}
					/// wetting angles
					phys->Delta1 = math::max(solution.delta1, solution.delta2);
					phys->Delta2 = math::min(solution.delta1, solution.delta2);
				} else
					scene->interactions->requestErase(interaction);

			} else {
				if (phys->vMeniscus == 0 and phys->capillaryPressure != 0) { //FIXME: test capillaryPressure consistently (!0 or >0 or >=0?!)
					Real Pinterpol = 0;
					Pinterpol      = phys->isBroken ? 0 : phys->capillaryPressure * R1 / liquidTension;
					if ((Pinterpol >= 0) && phys->meniscus) {
						MeniscusPhysicalData solution = DelaunayInterpolator::interpolate1(
						        dtPbased, K::Point_3(R2 / R1, Pinterpol, Dinterpol), phys->m, solutions, reset);
						/// capillary adhesion force
						phys->fCap = solution.force * R1 * liquidTension * currentContactGeometry->normal;
						/// meniscus volume
						Real Vinterpol   = solution.volume * pow(R1, 3);
						Real SInterface  = solution.surface * pow(R1, 2);
						phys->vMeniscus  = Vinterpol;
						phys->SInterface = SInterface;
						if (Vinterpol > 0) phys->meniscus = true;
						else
							phys->meniscus = false;

						if (phys->meniscus == false) phys->SInterface = 4 * Mathr::PI * (pow(R1, 2)) + 4 * Mathr::PI * (pow(R2, 2));
						if (!Vinterpol) {
							if (fusionDetection || phys->isBroken) bodiesMenisciiList.remove(interaction);
							if (D > ((interactionDetectionFactor - 1)
							         * (currentContactGeometry->radius2 + currentContactGeometry->radius1)))
								scene->interactions->requestErase(interaction);
						}
						/// wetting angles
						phys->Delta1 = math::max(solution.delta1, solution.delta2);
						phys->Delta2 = math::min(solution.delta1, solution.delta2);
					}
				} else {
					Real Vinterpol = phys->vMeniscus / pow(R1, 3);
					/// Capillary solution finder:
					if (phys->meniscus) {
						MeniscusPhysicalData solution = DelaunayInterpolator::interpolate1(
						        dtVbased, K::Point_3(R2 / R1, Vinterpol, Dinterpol), phys->m, solutions, reset);
						/// capillary adhesion force
						Real     Finterpol = solution.force;
						Vector3r fCap      = Finterpol * R1 * liquidTension * currentContactGeometry->normal;
						phys->fCap         = fCap;
						/// suction and interfacial area
						Real Pinterpol          = solution.succion * liquidTension / R1;
						Real SInterface         = solution.surface * pow(R1, 2);
						phys->capillaryPressure = Pinterpol;
						phys->SInterface        = SInterface;
						if (Finterpol > 0) phys->meniscus = true;
						else {
							phys->vMeniscus = 0;
							phys->meniscus  = false;
						}
						if (!Vinterpol) {
							if ((fusionDetection) || phys->isBroken) bodiesMenisciiList.remove(interaction);
							if (D > ((interactionDetectionFactor - 1)
							         * (currentContactGeometry->radius2 + currentContactGeometry->radius1)))
								scene->interactions->requestErase(interaction);
						}
						/// wetting angles
						phys->Delta1 = math::max(solution.delta1, solution.delta2);
						phys->Delta2 = math::min(solution.delta1, solution.delta2);
					}
				}
			}
			///interaction is not real //If the interaction is not real, it should not be in the list
		} else {
			if (fusionDetection) bodiesMenisciiList.remove(interaction);
			scene->interactions->requestErase(interaction);
		}
	}
	timingDeltas->checkpoint("compute all");
	if (fusionDetection) checkFusion();
	timingDeltas->checkpoint("check fusion");
}

shared_ptr<CapillaryPhysDelaunay> CapillarityEngine::solveStandalone(
        Real R1, Real R2, Real pressure, Real gap, shared_ptr<CapillaryPhysDelaunay> bridge = shared_ptr<CapillaryPhysDelaunay>() /*, bool pBased = true*/)
{
	if (dtPbased.number_of_vertices() < 1) triangulateData();
	if (not bridge) bridge = shared_ptr<CapillaryPhysDelaunay>(new CapillaryPhysDelaunay());

	Real Dinterpol = gap / R1;
	Real Pinterpol = pressure * R1 / liquidTension;

	MeniscusPhysicalData solution = DelaunayInterpolator::interpolate1(dtPbased, K::Point_3(R2 / R1, Pinterpol, Dinterpol), bridge->m, solutions, false);
	bridge->fCap                  = solution.force * Vector3r(1, 0, 0);
	bridge->vMeniscus             = solution.volume * pow(R1, 3);
	bridge->SInterface            = solution.surface * pow(R1, 2);
	bridge->meniscus              = solution.volume > 0;
	return bridge;
}

void CapillarityEngine::checkPhysType()
{
	// TODO: add more options here if we want more contact models compatible with capillarity
	// NOTE: We are assuming that only one type at a time is used in one simulation here
	for (auto ii = scene->interactions->begin(); ii != scene->interactions->end(); ++ii) {
		if (not(*ii)->phys) continue;
		if (CapillaryPhysDelaunay::getClassIndexStatic() == (*ii)->phys->getClassIndex()) {
			hertzOn          = false;
			hertzInitialized = true;
			break;
		} else if (MindlinCapillaryPhys::getClassIndexStatic() == (*ii)->phys->getClassIndex()) {
			hertzOn          = true;
			hertzInitialized = true;
			break;
		}
	}
	if (not hertzInitialized) LOG_ERROR("The capillary law is not implemented for the interaction physics used here, or no interactions were found.")
}


void CapillarityEngine::action()
{
	timingDeltas->start();
	if (fusionDetection and ompThreads > 1) {
		LOG_WARN("Fusion detection is not thread-safe. Turning ompThreads=1. (It is in fact 'most-likely' unsafe, not tested, comment-out this check "
		         "if you want to see by yourself)")
		ompThreads = 1;
	}
	if (!hertzInitialized) checkPhysType();
	bool switched = (pressureBased != (imposePressure or totalVolumeConstant));
	pressureBased = (imposePressure or totalVolumeConstant);

	void (CapillarityEngine::*solveBridges)(Real capillaryPressure, bool switched);
	solveBridges = hertzOn ? &CapillarityEngine::solveBridgesT<CapillaryMindlinPhysDelaunay> : &CapillarityEngine::solveBridgesT<CapillaryPhysDelaunay>;

	if (imposePressure) {
		(this->*solveBridges)(capillaryPressure, switched);
	} else {
		if (((totalVolumeConstant || (!totalVolumeConstant && firstIteration == 1)) && totalVolumeofWater != -1)
		    || (totalVolumeConstant && totalVolumeofWater == -1)) {
			if (!totalVolumeConstant) x = 1;
			totalVolumeConstant = 1;
			Real p0             = capillaryPressure;
			Real slope;
			Real eps = 0.0000000001;
			(this->*solveBridges)(p0, switched);
			Real V0 = waterVolume();
			if (totalVolumeConstant && totalVolumeofWater == -1 && firstIteration == 1) {
				totalVolumeofWater = V0;
				firstIteration += 1;
			}
			Real p1 = capillaryPressure + 0.1;
			(this->*solveBridges)(p1, switched);
			Real V1 = waterVolume();
			while (abs((totalVolumeofWater - V1) / totalVolumeofWater) > eps) {
				slope = (p1 - p0) / (V1 - V0);
				p0    = p1;
				V0    = V1;
				p1    = p1 - slope * (V1 - totalVolumeofWater);
				if (p1 < 0) {
					cout << "The requested volume of water is quite big, the simulation will continue at constant suction.:" << p0 << endl;
					capillaryPressure = p0;
					imposePressure    = 1;
					break;
				}
				(this->*solveBridges)(p1, switched);
				V1                = waterVolume();
				capillaryPressure = p1;
			}
			if (x == 1) {
				totalVolumeConstant = 0;
				firstIteration += 1;
			}
		} else {
			if ((!totalVolumeConstant && firstIteration == 1) && totalVolumeofWater == -1) {
				totalVolumeConstant = 1;
				(this->*solveBridges)(capillaryPressure, switched);
				firstIteration += 1;
				totalVolumeConstant = 0;
			} else {
				(this->*solveBridges)(capillaryPressure, switched);
			}
		}
		timingDeltas->checkpoint("closing loop");
	}

	const long size = scene->interactions->size();
#ifdef YADE_OPENMP
#pragma omp parallel for schedule(guided) num_threads(ompThreads > 0 ? std::min(ompThreads, omp_get_max_threads()) : omp_get_max_threads())
#endif
	for (long i = 0; i < size; i++) {
		const shared_ptr<Interaction>& I    = (*scene->interactions)[i];
		CapillaryPhysDelaunay*         phys = static_cast<CapillaryPhysDelaunay*>(I->phys.get());
		if (I->isReal() && phys->computeBridge == true) {
			CapillaryPhysDelaunay* cundallContactPhysics = NULL;
			MindlinCapillaryPhys*  mindlinContactPhysics = NULL;
			if (!hertzOn) cundallContactPhysics = static_cast<CapillaryPhysDelaunay*>(I->phys.get()); //use CapillaryPhysDelaunay for linear model
			else
				mindlinContactPhysics = static_cast<MindlinCapillaryPhys*>(I->phys.get()); //use MindlinCapillaryPhys for hertz model

			if ((hertzOn && mindlinContactPhysics->meniscus) || (!hertzOn && cundallContactPhysics->meniscus)) {
				if (fusionDetection) { //version with effect of fusion
					               //BINARY VERSION : if fusionNumber!=0 then no capillary force
					short int& fusionNumber = hertzOn ? mindlinContactPhysics->fusionNumber : cundallContactPhysics->fusionNumber;
					if (binaryFusion) {
						if (fusionNumber != 0) { //cerr << "fusion" << endl;
							hertzOn ? mindlinContactPhysics->fCap : cundallContactPhysics->fCap = Vector3r::Zero();
							continue;
						}
					}
					//LINEAR VERSION : capillary force is divided by (fusionNumber + 1) - NOTE : any decreasing function of fusionNumber can be considered in fact
					else if (fusionNumber != 0)
						hertzOn ? mindlinContactPhysics->fCap : cundallContactPhysics->fCap /= (fusionNumber + 1.);
				}
				scene->forces.addForce(I->getId1(), hertzOn ? mindlinContactPhysics->fCap : cundallContactPhysics->fCap);
				scene->forces.addForce(I->getId2(), -(hertzOn ? mindlinContactPhysics->fCap : cundallContactPhysics->fCap));
			}
		}
	}
	timingDeltas->checkpoint("addForces");
}
void CapillarityEngine::checkFusion()
{
	//Reset fusion numbers
	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for (; ii != iiEnd; ++ii) {
		if ((*ii)->isReal()) {
			if (!hertzOn) static_cast<CapillaryPhysDelaunay*>((*ii)->phys.get())->fusionNumber = 0;
			else
				static_cast<MindlinCapillaryPhys*>((*ii)->phys.get())->fusionNumber = 0;
		}
	}
	std::list<shared_ptr<Interaction>>::iterator firstMeniscus, lastMeniscus, currentMeniscus;
	Real                                         angle1 = -1.0;
	Real                                         angle2 = -1.0;

	for (int i = 0; i < bodiesMenisciiList.size(); ++i) { // i is the index (or id) of the body being tested
		CapillaryPhysDelaunay* cundallInteractionPhysics1 = NULL;
		MindlinCapillaryPhys*  mindlinInteractionPhysics1 = NULL;
		CapillaryPhysDelaunay* cundallInteractionPhysics2 = NULL;
		MindlinCapillaryPhys*  mindlinInteractionPhysics2 = NULL;
		if (!bodiesMenisciiList[i].empty()) {
			lastMeniscus = bodiesMenisciiList[i].end();
			for (firstMeniscus = bodiesMenisciiList[i].begin(); firstMeniscus != lastMeniscus;
			     ++firstMeniscus) { //FOR EACH MENISCUS ON THIS BODY...
				currentMeniscus = firstMeniscus;
				++currentMeniscus;
				if (!hertzOn) {
					cundallInteractionPhysics1 = YADE_CAST<CapillaryPhysDelaunay*>((*firstMeniscus)->phys.get());
					if (i == (*firstMeniscus)->getId1()) angle1 = cundallInteractionPhysics1->Delta1; //get angle of meniscus1 on body i
					else
						angle1 = cundallInteractionPhysics1->Delta2;
				} else {
					mindlinInteractionPhysics1 = YADE_CAST<MindlinCapillaryPhys*>((*firstMeniscus)->phys.get());
					if (i == (*firstMeniscus)->getId1()) angle1 = mindlinInteractionPhysics1->Delta1; //get angle of meniscus1 on body i
					else
						angle1 = mindlinInteractionPhysics1->Delta2;
				}
				for (; currentMeniscus != lastMeniscus; ++currentMeniscus) { //... CHECK FUSION WITH ALL OTHER MENISCII ON THE BODY
					if (!hertzOn) {
						cundallInteractionPhysics2 = YADE_CAST<CapillaryPhysDelaunay*>((*currentMeniscus)->phys.get());
						if (i == (*currentMeniscus)->getId1())
							angle2 = cundallInteractionPhysics2->Delta1; //get angle of meniscus2 on body i
						else
							angle2 = cundallInteractionPhysics2->Delta2;
					} else {
						mindlinInteractionPhysics2 = YADE_CAST<MindlinCapillaryPhys*>((*currentMeniscus)->phys.get());
						if (i == (*currentMeniscus)->getId1())
							angle2 = mindlinInteractionPhysics2->Delta1; //get angle of meniscus2 on body i
						else
							angle2 = mindlinInteractionPhysics2->Delta2;
					}
					if (angle1 == 0 || angle2 == 0) cerr << "THIS SHOULD NOT HAPPEN!!" << endl;

					Vector3r normalFirstMeniscus   = YADE_CAST<ScGeom*>((*firstMeniscus)->geom.get())->normal;
					Vector3r normalCurrentMeniscus = YADE_CAST<ScGeom*>((*currentMeniscus)->geom.get())->normal;

					Real normalDot = 0;
					if ((*firstMeniscus)->getId1() == (*currentMeniscus)->getId1()
					    || (*firstMeniscus)->getId2() == (*currentMeniscus)->getId2())
						normalDot = normalFirstMeniscus.dot(normalCurrentMeniscus);
					else
						normalDot = -(normalFirstMeniscus.dot(normalCurrentMeniscus));

					Real normalAngle = 0;
					if (normalDot >= 0) normalAngle = math::fastInvCos0(normalDot);
					else
						normalAngle = ((Mathr::PI)-math::fastInvCos0(-(normalDot)));

					if ((angle1 + angle2) * Mathr::DEG_TO_RAD > normalAngle) {
						if (!hertzOn) {
							++(cundallInteractionPhysics1->fusionNumber); //count +1 if 2 meniscii are overlaping
							++(cundallInteractionPhysics2->fusionNumber);
						} else {
							++(mindlinInteractionPhysics1->fusionNumber);
							++(mindlinInteractionPhysics2->fusionNumber);
						}
					};
				}
			}
		}
	}
}


BodiesMenisciiList1::BodiesMenisciiList1(Scene* scene)
{
	initialized = false;
	prepare(scene);
}

bool BodiesMenisciiList1::prepare(Scene* scene)
{
	interactionsOnBody.clear();
	shared_ptr<BodyContainer>& bodies = scene->bodies;

	Body::id_t              MaxId = -1;
	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();
	for (; bi != biEnd; ++bi) {
		MaxId = math::max(MaxId, (*bi)->getId());
	}
	interactionsOnBody.resize(MaxId + 1);
	for (unsigned int i = 0; i < interactionsOnBody.size(); ++i) {
		interactionsOnBody[i].clear();
	}

	InteractionContainer::iterator ii    = scene->interactions->begin();
	InteractionContainer::iterator iiEnd = scene->interactions->end();
	for (; ii != iiEnd; ++ii) {
		if ((*ii)->isReal()) {
			if (static_cast<CapillaryPhysDelaunay*>((*ii)->phys.get())->meniscus) insert(*ii);
		}
	}

	return initialized = true;
}

bool BodiesMenisciiList1::insert(const shared_ptr<Interaction>& interaction)
{
	interactionsOnBody[interaction->getId1()].push_back(interaction);
	interactionsOnBody[interaction->getId2()].push_back(interaction);
	return true;
}


bool BodiesMenisciiList1::remove(const shared_ptr<Interaction>& interaction)
{
	interactionsOnBody[interaction->getId1()].remove(interaction);
	interactionsOnBody[interaction->getId2()].remove(interaction);
	return true;
}

std::list<shared_ptr<Interaction>>& BodiesMenisciiList1::operator[](int index) { return interactionsOnBody[index]; }

int BodiesMenisciiList1::size() { return interactionsOnBody.size(); }

void BodiesMenisciiList1::display()
{
	std::list<shared_ptr<Interaction>>::iterator firstMeniscus;
	std::list<shared_ptr<Interaction>>::iterator lastMeniscus;
	for (unsigned int i = 0; i < interactionsOnBody.size(); ++i) {
		if (!interactionsOnBody[i].empty()) {
			lastMeniscus = interactionsOnBody[i].end();
			for (firstMeniscus = interactionsOnBody[i].begin(); firstMeniscus != lastMeniscus; ++firstMeniscus) {
				if (*firstMeniscus) {
					if (firstMeniscus->get()) cerr << "(" << (*firstMeniscus)->getId1() << ", " << (*firstMeniscus)->getId2() << ") ";
					else
						cerr << "(void)";
				}
			}
			cerr << endl;
		} else
			cerr << "empty" << endl;
	}
}

BodiesMenisciiList1::BodiesMenisciiList1() { initialized = false; }

} // namespace yade

#endif // YADE_CGAL
