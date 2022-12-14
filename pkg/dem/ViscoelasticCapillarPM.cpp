#include "ViscoelasticCapillarPM.hpp"
#include <lib/high-precision/Constants.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((ViscElCapMat)(ViscElCapPhys)(Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys)(Law2_ScGeom_ViscElCapPhys_Basic));

/* ViscElCapMat */
ViscElCapMat::~ViscElCapMat() { }

/* ViscElCapPhys */
ViscElCapPhys::~ViscElCapPhys() { }

/* Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys */
void Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	// no updates of an existing contact
	if (interaction->phys) return;

	TIMING_DELTAS_START();
	TIMING_DELTAS_CHECKPOINT("setup");

	shared_ptr<ViscElCapPhys> phys(new ViscElCapPhys());
	Calculate_ViscElMat_ViscElMat_ViscElPhys(b1, b2, interaction, phys);

	ViscElCapMat* mat1 = static_cast<ViscElCapMat*>(b1.get());
	ViscElCapMat* mat2 = static_cast<ViscElCapMat*>(b2.get());

	TIMING_DELTAS_CHECKPOINT("collide_materials");

	if (mat1->Capillar and mat2->Capillar) {
		if (mat1->Vb == mat2->Vb) {
			phys->Vb = mat1->Vb;
		} else {
			throw runtime_error("Vb should be equal for both particles!.");
		}

		if (mat1->gamma == mat2->gamma) {
			phys->gamma = mat1->gamma;
		} else {
			throw runtime_error("Gamma should be equal for both particles!.");
		}

		if (mat1->theta == mat2->theta) {
			phys->theta = (mat1->theta * M_PI / 180.0);
		} else {
			throw runtime_error("Theta should be equal for both particles!.");
		}


		if (mat1->dcap == mat2->dcap) {
			phys->dcap = mat1->dcap;
		} else {
			throw runtime_error("Theta should be equal for both particles!.");
		}

		if (mat1->CapillarType == mat2->CapillarType and mat2->CapillarType != "") {
			if (mat1->CapillarType == "Willett_numeric") {
				phys->CapillarType = Willett_numeric;
			} else if (mat1->CapillarType == "Willett_analytic") {
				phys->CapillarType = Willett_analytic;
			} else if (mat1->CapillarType == "Weigert") {
				phys->CapillarType = Weigert;
			} else if (mat1->CapillarType == "Rabinovich") {
				phys->CapillarType = Rabinovich;
			} else if (mat1->CapillarType == "Lambert") {
				phys->CapillarType = Lambert;
			} else if (mat1->CapillarType == "Soulie") {
				phys->CapillarType = Soulie;
			} else {
				phys->CapillarType = None_Capillar;
			}
		} else {
			throw runtime_error("CapillarType should be equal for both particles!.");
		}
		phys->Capillar = true;
	}
#ifdef YADE_LIQMIGRATION
	phys->LiqMigrEnabled = mat1->LiqMigrEnabled && mat2->LiqMigrEnabled;
#endif
	interaction->phys = phys;
}

/* Law2_ScGeom_ViscElCapPhys_Basic */
bool Law2_ScGeom_ViscElCapPhys_Basic::go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I)
{
	TIMING_DELTAS_START();
	TIMING_DELTAS_CHECKPOINT("setup");

	Vector3r force = Vector3r::Zero();

	const id_t id1 = I->getId1();
	const id_t id2 = I->getId2();

	const ScGeom&        geom   = *static_cast<ScGeom*>(_geom.get());
	Scene*               scene2 = Omega::instance().getScene().get();
	ViscElCapPhys&       phys   = *static_cast<ViscElCapPhys*>(_phys.get());
	const BodyContainer& bodies = *scene2->bodies;

	/*
   * This part for implementation of the capillar model.
   * All main equations are in calculateCapillarForce function. 
   * There is only the determination of critical distance between spheres, 
   * after that the liquid bridge will be broken.
   */

	TIMING_DELTAS_CHECKPOINT("create_liq_bridge");

	if (phys.Capillar and not(phys.liqBridgeCreated) and geom.penetrationDepth >= 0) {
		phys.liqBridgeCreated = true;
		phys.liqBridgeActive  = false;
#ifdef YADE_LIQMIGRATION
		if (phys.LiqMigrEnabled) { scene2->addIntrs.push_back(I); }
#endif
		Sphere* s1 = dynamic_cast<Sphere*>(bodies[id1]->shape.get());
		Sphere* s2 = dynamic_cast<Sphere*>(bodies[id2]->shape.get());
		if (s1 and s2) {
			phys.R = 2 * s1->radius * s2->radius / (s1->radius + s2->radius);
		} else if (s1 and not(s2)) {
			phys.R = s1->radius;
		} else {
			phys.R = s2->radius;
		}
	}

	TIMING_DELTAS_CHECKPOINT("calculate_scrit");

	phys.sCrit = this->critDist(phys.Vb, phys.R, phys.theta);

	TIMING_DELTAS_CHECKPOINT("force_calculation_liquid");
	if (geom.penetrationDepth < 0) {
		if (phys.liqBridgeCreated and -geom.penetrationDepth < phys.sCrit and phys.Capillar) {
			if (not(phys.liqBridgeActive)) {
				phys.liqBridgeActive = true;
				VLiqBridg += phys.Vb;
				NLiqBridg += 1;
			}

			const auto normalCapForceScalar = CapFunctionsPool[phys.CapillarType](geom, phys);
			Real       dampCapForceScalar   = 0.0;

			if (phys.dcap) {
				const State& de1 = *static_cast<State*>(bodies[id1]->state.get());
				const State& de2 = *static_cast<State*>(bodies[id2]->state.get());

				Vector3r& shearForce = phys.shearForce;
				if (I->isFresh(scene2)) shearForce = Vector3r(0, 0, 0);
				shearForce = geom.rotate(shearForce);

				const Vector3r shift2   = scene2->isPeriodic ? scene2->cell->intrShiftPos(I->cellDist) : Vector3r::Zero();
				const Vector3r shiftVel = scene2->isPeriodic ? scene2->cell->intrShiftVel(I->cellDist) : Vector3r::Zero();

				const Vector3r c1x = (geom.contactPoint - de1.pos);
				const Vector3r c2x = (geom.contactPoint - de2.pos - shift2);

				const Vector3r relativeVelocity = (de1.vel + de1.angVel.cross(c1x)) - (de2.vel + de2.angVel.cross(c2x)) + shiftVel;
				const auto     normalVelocity   = geom.normal.dot(relativeVelocity);

				dampCapForceScalar = -phys.dcap * normalVelocity;
			}

			phys.normalForce = -(normalCapForceScalar + dampCapForceScalar) * geom.normal;

			if (I->isActive) {
				addForce(id1, -phys.normalForce, scene2);
				addForce(id2, phys.normalForce, scene2);
			};
			return true;
		} else {
			if (phys.liqBridgeActive) {
				VLiqBridg -= phys.Vb;
				NLiqBridg -= 1;
			}
#ifdef YADE_LIQMIGRATION
			if (phys.Vb > 0.0 and ((phys.Vf1 + phys.Vf2) == 0.0)) {
				phys.Vf1 = phys.Vb / 2.0;
				phys.Vf2 = phys.Vb / 2.0;
			}
			const std::pair<id_t, Real> B1 = { id1, phys.Vf1 };
			const std::pair<id_t, Real> B2 = { id2, phys.Vf2 };
			scene2->delIntrs.push_back(B1);
			scene2->delIntrs.push_back(B2);
#endif
			return false;
		};
	};

	if (phys.liqBridgeActive) {
		phys.liqBridgeActive = false;
		VLiqBridg -= phys.Vb;
		NLiqBridg -= 1;
	}

	if (I->isActive) {
		Vector3r torque1 = Vector3r::Zero();
		Vector3r torque2 = Vector3r::Zero();

		TIMING_DELTAS_CHECKPOINT("force_calculation_penetr");
		if (computeForceTorqueViscEl(_geom, _phys, I, force, torque1, torque2)) {
			addForce(id1, -force, scene2);
			addForce(id2, force, scene2);
			addTorque(id1, torque1, scene2);
			addTorque(id2, torque2, scene2);
			return true;
		} else {
			return false;
		}
	}
	return true;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::critDist(const Real& Vb, const Real& R, const Real& Theta)
{
	const Real Vstar = Vb / (R * R * R);
	const Real Sstar
	        = (1 + 0.5 * Theta) * (pow(Vstar, 1 / 3.0) + 0.1 * pow(Vstar, 2.0 / 3.0)); // [Willett2000], equation (15), use the full-length e.g 2*Sc
	const Real critDist = Sstar * R;
	return critDist;
}

//=========================================================================================
//======================Capillary bridge models============================================
//=========================================================================================

Real Law2_ScGeom_ViscElCapPhys_Basic::Willett_numeric_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/* 
   * Capillar model from [Willett2000]
   */

	const Real R  = phys.R;
	const Real s  = -geom.penetrationDepth;
	const Real Vb = phys.Vb;

	const Real VbS   = Vb / (R * R * R);
	const Real Th1   = phys.theta;
	const Real Th2   = phys.theta * phys.theta;
	const Real Gamma = phys.gamma;

	/*
   * [Willett2000], equations in Attachment
  */
	const Real f1 = (-0.44507 + 0.050832 * Th1 - 1.1466 * Th2) + (-0.1119 - 0.000411 * Th1 - 0.1490 * Th2) * log(VbS)
	        + (-0.012101 - 0.0036456 * Th1 - 0.01255 * Th2) * log(VbS) * log(VbS)
	        + (-0.0005 - 0.0003505 * Th1 - 0.00029076 * Th2) * log(VbS) * log(VbS) * log(VbS);

	const Real f2 = (1.9222 - 0.57473 * Th1 - 1.2918 * Th2) + (-0.0668 - 0.1201 * Th1 - 0.22574 * Th2) * log(VbS)
	        + (-0.0013375 - 0.0068988 * Th1 - 0.01137 * Th2) * log(VbS) * log(VbS);

	const Real f3 = (1.268 - 0.01396 * Th1 - 0.23566 * Th2) + (0.198 + 0.092 * Th1 - 0.06418 * Th2) * log(VbS)
	        + (0.02232 + 0.02238 * Th1 - 0.009853 * Th2) * log(VbS) * log(VbS)
	        + (0.0008585 + 0.001318 * Th1 - 0.00053 * Th2) * log(VbS) * log(VbS) * log(VbS);

	const Real f4 = (-0.010703 + 0.073776 * Th1 - 0.34742 * Th2) + (0.03345 + 0.04543 * Th1 - 0.09056 * Th2) * log(VbS)
	        + (0.0018574 + 0.004456 * Th1 - 0.006257 * Th2) * log(VbS) * log(VbS);

	const Real sPl = (s / 2.0) / sqrt(Vb / R);

	const Real lnFS = f1 - f2 * exp(f3 * log(sPl) + f4 * log(sPl) * log(sPl));
	const Real FS   = exp(lnFS);

	const Real fC = FS * 2.0 * M_PI * R * Gamma;
	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::Willett_analytic_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/* 
   * Capillar model from Willet [Willett2000] (analytical solution), but 
   * used also in the work of Herminghaus [Herminghaus2005]
   */

	const Real R     = phys.R;
	const Real Gamma = phys.gamma;
	const Real s     = -geom.penetrationDepth;
	const Real Vb    = phys.Vb;

	/*
  Real sPl = s/sqrt(Vb/R);                                                            // [Herminghaus2005], equation (sentence between (7) and (8))
  fC = 2.0 * M_PI* R * Gamma * cos(phys.theta)/(1 + 1.05*sPl + 2.5 *sPl * sPl);       // [Herminghaus2005], equation (7)
  */

	const Real sPl    = (s / 2.0) / sqrt(Vb / R); // [Willett2000], equation (sentence after (11)), s - half-separation, so s*2.0
	const Real f_star = cos(phys.theta) / (1 + 2.1 * sPl + 10.0 * pow(sPl, 2.0)); // [Willett2000], equation (12)
	const Real fC     = f_star * (2 * M_PI * R * Gamma);                          // [Willett2000], equation (13), against F

	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::Weigert_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/*
  *  Capillar model from [Weigert1999]
  */
	const Real R  = phys.R;
	const Real a  = -geom.penetrationDepth;
	const Real Ca = (1.0 + 6.0 * a / (R * 2.0));   // [Weigert1999], equation (16)
	const Real Ct = (1.0 + 1.1 * sin(phys.theta)); // [Weigert1999], equation (17)

	/*
  Real Eps = 0.36;                                                                          // Porosity
  Real fi = phys.Vb/(2.0*M_PI/6.0*pow(R*2.0,3.));                                           // [Weigert1999], equation (13)
  Real S = M_PI*(1-Eps)/(pow(Eps, 2.0))*fi;                                                 // [Weigert1999], equation (14)
  Real beta = asin(pow(((S/0.36)*(pow(Eps, 2.0)/(1-Eps))*(1.0/Ca)*(1.0/Ct)), 1.0/4.0));     // [Weigert1999], equation (19)
  */

	const Real beta = asin(pow(phys.Vb / (0.12 * Ca * Ct * pow(2.0 * R, 3.0)), 1.0 / 4.0)); // [Weigert1999], equation (15), against Vb

	const Real r1 = (2.0 * R * (1 - cos(beta)) + a) / (2.0 * cos(beta + phys.theta)); // [Weigert1999], equation (5)
	const Real r2 = R * sin(beta) + r1 * (sin(beta + phys.theta) - 1);                // [Weigert1999], equation (6)
	const Real Pk = phys.gamma * (1 / r1 - 1 / r2);                                   // [Weigert1999], equation (22),
	                                                                                  // see also a sentence over the equation
	                                                                                  // "R1 was taken as positive and R2 was taken as negative"

	//fC = M_PI*2.0*R*phys.gamma/(1+tan(0.5*beta));                                           // [Weigert1999], equation (23), [Fisher]

	const Real fC = M_PI / 4.0 * pow((2.0 * R), 2.0) * pow(sin(beta), 2.0) * Pk
	        + phys.gamma * M_PI * 2.0 * R * sin(beta) * sin(beta + phys.theta); // [Weigert1999], equation (21)

	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::Rabinovich_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/* 
   * Capillar model from Rabinovich [Rabinov2005]
   *
   * This formulation from Rabinovich has been later verified and corrected
   * by Lambert [Lambert2008]. So we can calculate both formulations
   * 
   */

	const Real R     = phys.R;
	const Real Gamma = phys.gamma;
	const Real H     = -geom.penetrationDepth;
	const Real V     = phys.Vb;

	Real fC  = 0.0;
	Real dsp = 0.0;

	if (H != 0.0) {
		dsp              = H / 2.0 * (-1.0 + sqrt(1.0 + 2.0 * V / (M_PI * R * H * H)));       // [Rabinov2005], equation (20)
		fC               = -(2 * M_PI * R * Gamma * cos(phys.theta)) / (1 + (H / (2 * dsp))); // [Lambert2008], equation (65), taken from [Rabinov2005]
		const Real alpha = sqrt(H / R * (-1 + sqrt(1 + 2.0 * V / (M_PI * R * H * H))));       // [Rabinov2005], equation (A3)
		fC -= 2 * M_PI * R * Gamma * sin(alpha) * sin(phys.theta + alpha);                    // [Rabinov2005], equation (19)
	} else {
		fC               = -(2 * M_PI * R * Gamma * cos(phys.theta));
		const Real alpha = 0.0;
		fC -= 2 * M_PI * R * Gamma * sin(alpha) * sin(phys.theta + alpha); // [Rabinov2005], equation (19)
	}

	fC *= -1;
	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::Lambert_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/* 
   * Capillar model from Rabinovich [Rabinov2005]
   *
   * This formulation from Rabinovich has been later verified and corrected
   * by Lambert [Lambert2008]. So we can calculate both formulations
   * 
   */

	const Real R     = phys.R;
	const Real Gamma = phys.gamma;
	const Real H     = -geom.penetrationDepth;
	const Real V     = phys.Vb;

	Real fC  = 0.0;
	Real dsp = 0.0;

	if (H != 0.0) {
		dsp = H / 2.0 * (-1.0 + sqrt(1.0 + 2.0 * V / (M_PI * R * H * H)));       // [Rabinov2005], equation (20)
		fC  = -(2 * M_PI * R * Gamma * cos(phys.theta)) / (1 + (H / (2 * dsp))); // [Lambert2008], equation (65), taken from [Rabinov2005]
	} else {
		fC = -(2 * M_PI * R * Gamma * cos(phys.theta));
	}

	fC *= -1;
	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::Soulie_f(const ScGeom& geom, ViscElCapPhys& phys)
{
	/* 
   * Capillar model from Soulie [Soulie2006]
   *
   * !!! In this implementation the radiis of particles are taken equal
   * to get the symmetric forces.
   * 
   * Please, use this model only for testing purposes.
   * 
   */

	const Real R     = phys.R;
	const Real Gamma = phys.gamma;
	const Real D     = -geom.penetrationDepth;
	const Real V     = phys.Vb;
	const Real Theta = phys.theta;

	const Real a = -1.1 * pow((V / (R * R * R)), -0.53);
	const Real b = (-0.148 * log(V / (R * R * R)) - 0.96) * Theta * Theta - 0.0082 * log(V / (R * R * R)) + 0.48;
	const Real c = 0.0018 * log(V / (R * R * R)) + 0.078;

	const Real fC = Mathr::PI * Gamma * sqrt(R * R) * (c + exp(a * D / R + b));

	return fC;
}

Real Law2_ScGeom_ViscElCapPhys_Basic::None_f(const ScGeom& /*geom*/, ViscElCapPhys& /*phys*/) { return 0; }

#ifdef YADE_LIQMIGRATION
YADE_PLUGIN((LiqControl));
void LiqControl::action()
{
	// This function implements liquid migration model, introduced here [Mani2013]
	mapBodyInt bI;
	mapBodyInt bodyNeedUpdate;

	// Calculate, how much new contacts will be at each body
	for (unsigned int i = 0; i < scene->addIntrs.size(); i++) {
		shared_ptr<Body> b1 = Body::byId(scene->addIntrs[i]->getId1(), scene);
		shared_ptr<Body> b2 = Body::byId(scene->addIntrs[i]->getId2(), scene);

		if (not(mask != 0 && ((b1->groupMask & b2->groupMask & mask) == 0))) {
			addBodyMapInt(bI, scene->addIntrs[i]->getId1());
			addBodyMapInt(bI, scene->addIntrs[i]->getId2());
		}
	}

	// Update volume water at each deleted interaction for each body
	for (unsigned int i = 0; i < scene->delIntrs.size(); i++) {
		shared_ptr<Body> b = Body::byId(scene->delIntrs[i].first, scene);
		b->state->Vf += scene->delIntrs[i].second;
		addBodyMapInt(bodyNeedUpdate, scene->delIntrs[i].first);
		liqVolRup += scene->delIntrs[i].second;
	}
	scene->delIntrs.clear();

	// Update volume bridge at each new added interaction
	mapBodyReal bodyUpdateLiquid;
	for (unsigned int i = 0; i < scene->addIntrs.size(); i++) {
		shared_ptr<Body> b1 = Body::byId(scene->addIntrs[i]->getId1(), scene);
		shared_ptr<Body> b2 = Body::byId(scene->addIntrs[i]->getId2(), scene);

		const id_t id1 = b1->id;
		const id_t id2 = b2->id;

		ViscElCapPhys* Vb = dynamic_cast<ViscElCapPhys*>(scene->addIntrs[i]->phys.get());

		if (mask != 0 && ((b1->groupMask & b2->groupMask & mask) == 0)) {
			Vb->Vb   = 0.0;
			Vb->Vf1  = 0.0;
			Vb->Vf2  = 0.0;
			Vb->Vmax = 0.0;
		} else {
			const Real Vmax = vMax(b1, b2);
			Vb->Vmax        = Vmax;

			Real Vf1 = 0.0;
			Real Vf2 = 0.0;

			if ((b1->state->Vmin) < b1->state->Vf) { Vf1 = (b1->state->Vf - b1->state->Vmin) / bI[id1]; }
			if ((b2->state->Vmin) < b2->state->Vf) { Vf2 = (b2->state->Vf - b2->state->Vmin) / bI[id2]; }

			Real Vrup = Vf1 + Vf2;

			if (Vrup > Vmax) {
				Vf1 *= Vmax / Vrup;
				Vf2 *= Vmax / Vrup;
				Vrup = Vf1 + Vf2;
			}

			liqVolShr += Vrup;
			addBodyMapReal(bodyUpdateLiquid, id1, -Vf1);
			addBodyMapReal(bodyUpdateLiquid, id2, -Vf2);

			Vb->Vb = Vrup;
			if (particleconserve) {
				Vb->Vf1 = Vf1;
				Vb->Vf2 = Vf2;
			} else {
				Vb->Vf1 = Vrup / 2.0;
				Vb->Vf2 = Vrup / 2.0;
			}
		}
	}

	scene->addIntrs.clear();

	// Update water volume in body
	for (mapBodyReal::const_iterator it = bodyUpdateLiquid.begin(); it != bodyUpdateLiquid.end(); ++it) {
		Body::byId(it->first)->state->Vf += it->second;
	}

	// Update contacts around body
	for (mapBodyInt::const_iterator it = bodyNeedUpdate.begin(); it != bodyNeedUpdate.end(); ++it) {
		updateLiquid(Body::byId(it->first));
	}
}

void LiqControl::updateLiquid(shared_ptr<Body> b)
{
	if (b->state->Vf <= b->state->Vmin) {
		return;
	} else {
		// How much liquid can body share
		const Real LiqCanBeShared = b->state->Vf - b->state->Vmin;

		// Check how much liquid can accept contacts
		Real         LiqContactsAccept = 0.0;
		unsigned int contactN          = 0;
		for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
			if (!((*it).second) or !(((*it).second)->isReal())) continue;
			ViscElCapPhys* physT = dynamic_cast<ViscElCapPhys*>(((*it).second)->phys.get());
			if ((physT->Vb < physT->Vmax) and (physT->Vmax > 0)) {
				LiqContactsAccept += physT->Vmax - physT->Vb;
				contactN++;
			}
		}

		if (contactN > 0) {
			//There are some contacts, which can be filled
			Real FillLevel = 0.0;
			if (LiqContactsAccept > LiqCanBeShared) { // Share all available liquid from body to contacts
				const Real LiquidWillBeShared = b->state->Vf - b->state->Vmin;
				b->state->Vf                  = b->state->Vmin;
				FillLevel                     = LiquidWillBeShared / LiqContactsAccept;
			} else { // Not all available liquid from body can be shared
				b->state->Vf -= LiqContactsAccept;
				FillLevel = 1.0;
			}

			for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
				if (!((*it).second) or !(((*it).second)->isReal())) continue;
				ViscElCapPhys* physT = dynamic_cast<ViscElCapPhys*>(((*it).second)->phys.get());
				if ((physT->Vb < physT->Vmax) and (physT->Vmax > 0)) {
					const Real addVolLiq = (physT->Vmax - physT->Vb) * FillLevel;
					liqVolShr += addVolLiq;
					physT->Vb += addVolLiq;
					if (particleconserve) {
						if (((*it).second)->getId1() == (*it).first) {
							physT->Vf2 += addVolLiq;
						} else if (((*it).second)->getId2() == (*it).first) {
							physT->Vf1 += addVolLiq;
						}
					} else {
						physT->Vf1 += addVolLiq / 2.0;
						physT->Vf2 += addVolLiq / 2.0;
					}
				}
			}
			return;
		} else {
			return;
		}
	}
}

void LiqControl::addBodyMapInt(mapBodyInt& m, Body::id_t b)
{
	mapBodyInt::const_iterator got;
	got = m.find(b);
	if (got == m.end()) {
		m.insert(mapBodyInt::value_type(b, 1));
	} else {
		m[b] += 1;
	}
}

void LiqControl::addBodyMapReal(mapBodyReal& m, Body::id_t b, Real addV)
{
	mapBodyReal::const_iterator got;
	got = m.find(b);
	if (got == m.end()) {
		m.insert(mapBodyReal::value_type(b, addV));
	} else {
		m[b] += addV;
	}
}

Real LiqControl::vMax(shared_ptr<Body> const b1, shared_ptr<Body> const b2)
{
	Sphere* s1   = dynamic_cast<Sphere*>(b1->shape.get());
	Sphere* s2   = dynamic_cast<Sphere*>(b2->shape.get());
	Real    minR = 0.0;
	if (s1 and s2) {
		minR = std::min(s1->radius, s2->radius);
	} else if (s1 and not(s2)) {
		minR = s1->radius;
	} else {
		minR = s2->radius;
	}
	return vMaxCoef * minR * minR * minR;
}

Real liqVolIterBody(shared_ptr<Body> b)
{
	Real LiqVol = 0.0;
	if (!b) {
		return LiqVol;
	} else {
		for (Body::MapId2IntrT::iterator it = b->intrs.begin(), end = b->intrs.end(); it != end; ++it) {
			if (!((*it).second) or !(((*it).second)->isReal())) continue;
			ViscElCapPhys* physT = dynamic_cast<ViscElCapPhys*>(((*it).second)->phys.get());
			if (physT and physT->Vb and physT->Vb > 0) {
				if (((*it).second)->id1 == b->id) {
					if (physT->Vf1 > 0 or physT->Vf2 > 0) {
						LiqVol += physT->Vf1;
					} else {
						LiqVol += physT->Vb / 2.0;
					}
				} else {
					if (physT->Vf1 > 0 or physT->Vf2 > 0) {
						LiqVol += physT->Vf2;
					} else {
						LiqVol += physT->Vb / 2.0;
					}
				}
			}
		}
		return LiqVol;
	}
}

Real LiqControl::liqVolBody(id_t id) const
{
	Scene*               scene2 = Omega::instance().getScene().get();
	const BodyContainer& bodies = *scene2->bodies;
	if (bodies[id]) {
		if (bodies[id]->state->Vf > 0) {
			return bodies[id]->state->Vf + liqVolIterBody(bodies[id]);
		} else {
			return liqVolIterBody(bodies[id]);
		}
	} else
		return -1;
}

Real LiqControl::totalLiqVol(int mask2 = 0) const
{
	Scene* scene2      = Omega::instance().getScene().get();
	Real   totalLiqVol = 0.0;
	for (const auto& b : *scene2->bodies) {
		if ((mask2 > 0 && (b->groupMask & mask2) == 0) or (!b)) continue;
		totalLiqVol += liqVolIterBody(b);
		if (b->state->Vf > 0) { totalLiqVol += b->state->Vf; }
	}
	return totalLiqVol;
}

bool LiqControl::addLiqInter(id_t id1, id_t id2, Real liq)
{
	if (id1 == id2 or liq <= 0) return false;

	Scene*                            scene2 = Omega::instance().getScene().get();
	shared_ptr<InteractionContainer>& intrs  = scene2->interactions;
	const shared_ptr<Interaction>&    I      = intrs->find(id1, id2);
	if (I->isReal()) {
		ViscElCapPhys* physT = dynamic_cast<ViscElCapPhys*>(I->phys.get());
		if (physT and physT->Vb <= physT->Vmax and liq <= (physT->Vmax - physT->Vb)) {
			physT->Vb += liq;
			physT->Vf1 += liq / 2.0;
			physT->Vf2 += liq / 2.0;
			return true;
		}
	}
	return false;
}
#endif

} // namespace yade
