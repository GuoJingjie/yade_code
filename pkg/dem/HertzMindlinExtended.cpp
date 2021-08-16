/*
 2016 - Bettina Suhr 
 extension of Hertz-Mindlin model (see HertzMindlin.cpp): 

 normal direction: conical damage model
 tangential direction: stress dependent interparticle friction coefficient
 both models can be switched on/off separately

 references:
 Harkness, Zervos, Le Pen, Aingaran, Powrie: Discrete element simulation of railway ballast: modelling cell pressure effects in triaxial tests, Granular Matter, (2016) 18:65
 Suhr & Six 2017: Parametrisation of a DEM model for railway ballast under different load cases, Granular Matter, (2017) 19:64
 Suhr & Six 2016: On the effect of stress dependent interparticle friction in direct shear tests Powder Technology , (2016) 294:211 - 220
*/

#include "HertzMindlinExtended.hpp"
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>
#include<core/Scene.hpp>
#include <lib/high-precision/Constants.hpp>


namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((FrictMatCDM)(MindlinPhysCDM)(Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM)(Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM)(Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM));

CREATE_LOGGER(Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM);


void Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return; // no updates of an already existing contact necessary

	interaction->phys = shared_ptr<MindlinPhysCDM>(new MindlinPhysCDM());
	const shared_ptr<MindlinPhysCDM>& contactPhysics = YADE_PTR_CAST<MindlinPhysCDM>(interaction->phys);

	const FrictMatCDM* mat1 = YADE_CAST<FrictMatCDM*>(b1.get());
	const FrictMatCDM* mat2 = YADE_CAST<FrictMatCDM*>(b2.get());
	
	/* from interaction physics */
	const Real Ea = mat1->young;
	const Real Eb = mat2->young;
	const Real Va = mat1->poisson;
	const Real Vb = mat2->poisson;
	const Real fa = mat1->frictionAngle;
	const Real fb = mat2->frictionAngle;


	/* from interaction geometry */
	const GenericSpheresContact* scg = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());		

	/* calculate stiffness coefficients */
	const Real Ga = Ea/(2*(1+Va));
	const Real Gb = Eb/(2*(1+Vb));
	//Real V = (Va+Vb)/2; // average of poisson's ratio
	const Real E = Ea*Eb/((1.-std::pow(Va,2))*Eb+(1.-std::pow(Vb,2))*Ea); // equivalent Young's modulus
	
	//CHANGED use equivalent radius from geometry!!
	const Real Da = scg->refR1>0 ? scg->refR1 : scg->refR2; 
	const Real Db = scg->refR2; 
	const Real R = Da*Db/(Da+Db); // equivalent radius

	const Real Kno = 4./3.*E*sqrt(R); // coefficient for normal stiffness
	//Real Kso = 2*sqrt(4*R)*G/(2-V); // coefficient for shear stiffness
	//CHANGED
	const Real Kso = 8*sqrt(R)/( (2-Va)/Ga + (2-Vb)/Gb); // coefficient for shear stiffness
	const Real frictionAngle = (!frictAngle) ? std::min(fa,fb) : (*frictAngle)(mat1->id,mat2->id,mat1->frictionAngle,mat2->frictionAngle);


	/* pass values calculated from above to MindlinPhys */
	//contactPhysics->prevNormal = scg->normal; // used to compute relative rotation
	contactPhysics->E = E; // equiv young modulus
	contactPhysics->G = 1.0/( (2-Va)/Ga + (2-Vb)/Gb);//equiv shear modulus
	contactPhysics->kno = Kno; // this is just a coeff, will be changed
	contactPhysics->kso = Kso; // this is just a coeff, will be changed
	//parameters for conical damage model
	contactPhysics->R=R;//save current contact radius, will increase during yielding
	contactPhysics->radius=R;//HERE Rmin is stored!! will not be changed during contact lifetime
	contactPhysics->sigmaMax=std::min(mat1->sigmaMax,mat2->sigmaMax);
	contactPhysics->alphaFac=(1.0-std::sin(std::min(mat1->alpha,mat2->alpha)))/std::sin(std::min(mat1->alpha,mat2->alpha));
    //parameter for stress dependent interparticle friction coefficient
	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle); //current value of friction coefficient is stored here, will change
	contactPhysics->mu0=std::tan(frictionAngle); // parameter for stress dependent interparticle friction coefficient, will not change during contact lifetime
	contactPhysics->c1=std::min(mat1->c1,mat2->c1); 
	contactPhysics->c2=std::min(mat1->c2,mat2->c2);

	
	if (std::min(mat1->alpha,mat2->alpha) <= 0 or std::min(mat1->alpha,mat2->alpha)>=Mathr::PI/2.0) throw  std::invalid_argument("Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM: alpha must in (0,pi/2) radians ,NOT equal to 0 or pi/2");
	if (contactPhysics->mu0 <= 0) throw  std::invalid_argument("Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM: mu0/friction angle must be > 0");
	if (contactPhysics->sigmaMax <= 0) throw std::invalid_argument("Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM: sigmaMax must be > 0");
	if (contactPhysics->sigmaMax >= E) throw  std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: sigmaMax must be < Young's modulus!");
	if (contactPhysics->c1 <0) throw  std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: c1 must be >=0!");
	if (contactPhysics->c2 <0) throw  std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: c2 must be >=0!");

	//damping currently cot used-----------------------------------------------------
	contactPhysics->alpha=0.0;
	contactPhysics->betan=0.0;
	contactPhysics->betas=0.0;
	//adhesion, bending currently cot used-----------------------------------------------------	
	contactPhysics->adhesionForce = 0;//Skipped at the moment! Adhesion;
	contactPhysics->kr = 0;//Skipped at the moment! krot;
	contactPhysics->ktw = 0;//Skipped at the moment! ktwist;
	contactPhysics->maxBendPl = 0;//Skipped at the moment! 
  
};

CREATE_LOGGER(Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM);

void Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM::go( const shared_ptr<Material>& b1
					, const shared_ptr<Material>& b2
					, const shared_ptr<Interaction>& interaction)
{
	if(interaction->phys) return; // no updates of an already existing contact necessary
	shared_ptr<MindlinPhysCDM> contactPhysics(new MindlinPhysCDM());
	interaction->phys = contactPhysics;
    const FrictMat* matFrictMat;
    const FrictMatCDM* matFrictMatCDM;
    // check which interaction partner is of which material class
    // try to cast b1 to FrictMatCDM, returns null if b1 is FrictMat
    const FrictMatCDM* matFrictMatCDM1= dynamic_cast<FrictMatCDM*>(b1.get());
    if  (matFrictMatCDM1){
        // b1 is FrictMatCDM, b2 is FrictMat
        matFrictMatCDM = YADE_CAST<FrictMatCDM*>(b1.get());
        matFrictMat = YADE_CAST<FrictMat*>(b2.get());
    }else{
        // b1 is FrictMat, b2 is FrictMatCDM
        matFrictMatCDM = YADE_CAST<FrictMatCDM*>(b2.get());
        matFrictMat = YADE_CAST<FrictMat*>(b1.get());
    }
    
	/* from interaction physics */
	const Real Ea = matFrictMat->young;
	const Real Eb = matFrictMatCDM->young;
	const Real Va = matFrictMat->poisson;
	const Real Vb = matFrictMatCDM->poisson;
	const Real fa = matFrictMat->frictionAngle;
	const Real fb = matFrictMatCDM->frictionAngle;


	/* from interaction geometry */
	const GenericSpheresContact* scg = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());		
	if (Va <= 0 or Vb<=0) throw std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: Poisson's ratio must be > 0");

	/* calculate stiffness coefficients */
	const Real Ga = Ea/(2*(1+Va));
	const Real Gb = Eb/(2*(1+Vb));
	//Real V = (Va+Vb)/2; // average of poisson's ratio
	const Real E = Ea*Eb/((1.-std::pow(Va,2))*Eb+(1.-std::pow(Vb,2))*Ea); //equivalent Young's modulus
	
	const Real Da = scg->refR1>0 ? scg->refR1 : scg->refR2; 
	const Real Db = scg->refR2; 
	const Real R = Da*Db/(Da+Db); // equivalent radius

	const Real Kno = 4./3.*E*sqrt(R); // coefficient for normal stiffness
	const Real Kso = 8*sqrt(R)/( (2-Va)/Ga + (2-Vb)/Gb); // coefficient for shear stiffness
	const Real frictionAngle = (!frictAngle) ? std::min(fa,fb) : (*frictAngle)(matFrictMat->id,matFrictMatCDM->id,matFrictMat->frictionAngle,matFrictMatCDM->frictionAngle);


	/* pass values calculated from above to MindlinPhys */
	contactPhysics->E = E; // equiv Young's modulus
	contactPhysics->G = 1.0/( (2-Va)/Ga + (2-Vb)/Gb);//equiv shear modulus
	contactPhysics->kno = Kno; // this is just a coeff, will be changed
	contactPhysics->kso = Kso; // this is just a coeff, will be changed
	//parameters for conical damage model
	contactPhysics->R=R;//save current contact radius, will increase during yielding
	contactPhysics->radius=R;//HERE Rmin is stored!! will not be changed during contact lifetime
	contactPhysics->sigmaMax=matFrictMatCDM->sigmaMax;
	contactPhysics->alphaFac=(1.0-std::sin(matFrictMatCDM->alpha))/std::sin(matFrictMatCDM->alpha);
    //parameter for stress dependent interparticle friction coefficient
	//friction coeff is ALWAYS constant!! 
	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle); 
	contactPhysics->mu0=std::tan(frictionAngle); 
	contactPhysics->c1=0.0;
	contactPhysics->c2=0.0;

	if (matFrictMatCDM->alpha <= 0 or matFrictMatCDM->alpha>=Mathr::PI/2.0) throw std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: alpha must in (0,pi/2) radians ,NOT equal to 0 or pi/2");
	if (contactPhysics->mu0 <= 0) throw std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: mu0/frictionAngle must be > 0");
	if (contactPhysics->sigmaMax <= 0) throw  std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: sigmaMax must be > 0");
	if (contactPhysics->sigmaMax >= E) throw  std::invalid_argument("Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM: sigmaMax must be < Young's modulus!");
	//damping currently cot used-----------------------------------------------------
	contactPhysics->alpha=0.0;
	contactPhysics->betan=0.0;
	contactPhysics->betas=0.0;
	//adhesion, bending currently cot used-----------------------------------------------------
	contactPhysics->adhesionForce = 0;//Skipped at the moment! Adhesion;
	contactPhysics->kr = 0;//Skipped at the moment! krot;
	contactPhysics->ktw = 0;//Skipped at the moment! ktwist;
	contactPhysics->maxBendPl = 0;//Skipped at the moment! 
	
};


/* Function which returns the ratio between the number of sliding contacts to the total number at a given time */
Real Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM::ratioSlidingContacts()
{
	Real ratio(0); int count(0);
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		const MindlinPhysCDM* phys = dynamic_cast<MindlinPhysCDM*>(I->phys.get());
		if (phys->isSliding) {ratio+=1;}
		count++;
	}  
	ratio/=count;
	return ratio;
}

/* Function which returns the ratio between the number of yielding contacts to the total number at a given time */
Real Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM::ratioYieldingContacts()
{
	Real ratio(0); int count(0);
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		const MindlinPhysCDM* phys = dynamic_cast<MindlinPhysCDM*>(I->phys.get());
		if (phys->isYielding) {ratio+=1;}
		count++;
	}  
	ratio/=count;
	return ratio;
}

/******************** Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM *********/
CREATE_LOGGER(Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM);

bool Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
  
	
	//const Body::id_t id1 = contact->getId1(); // get id body 1
 	//const Body::id_t id2 = contact->getId2(); // get id body 2
    const auto id1 = contact->getId1();
    const auto id2 = contact->getId2();
    
	const State* de1 = Body::byId(id1,scene)->state.get();
	const State* de2 = Body::byId(id2,scene)->state.get();	

	const ScGeom* scg = static_cast<ScGeom*>(ig.get());
	MindlinPhysCDM* phys = static_cast<MindlinPhysCDM*>(ip.get());	
 
	/****************/
	/* NORMAL FORCE */
	/****************/
	const Real uDEM = scg->penetrationDepth;// get DEM overlap
	if (uDEM<0){
	  if (neverErase) {phys->shearForce = phys->normalForce = Vector3r::Zero(); phys->kn=phys->ks=0; return true;}
 	  else return false;
	}
  

	// Hertz-Mindlin's formulation + conical damage model 
	//conical damage model: Harkness et al. 2016 (see header for full reference)
	//modification in Suhr&Six 2017 (see header for full reference)
	//DIFFERENT INTERPRETATION OF OVERLAP DEFINITION
	//phys->radius: R_eq as in Hertz law
	//phys->R: current contact radius
	Real uN = uDEM + (phys->radius-phys->R) * phys->alphaFac;//elastic overlap
	if (uN<0) {
		phys->shearForce = phys->normalForce = Vector3r::Zero();
		phys->kn=phys->ks=0;
		return true;
	}
	//reformulated formula
	phys->isYielding=false;
	//check yield condition
	if (2.0 * phys->E / Mathr::PI *std::pow(uN / phys->R, 0.5) > phys->sigmaMax){
	  phys->isYielding=true;
	  Real hfac = std::pow(  Mathr::PI * phys->sigmaMax / 2.0 / phys->E , 2.0);
	  phys->R =  (uDEM + phys->radius*phys->alphaFac) / ( hfac + phys->alphaFac);
	  uN = uDEM + (phys->radius-phys->R) * phys->alphaFac;//adapted elastic overlap
	}
	
	// here we store the value of kn to compute the time step
	phys->kn = 4./3.0 * phys->E * std::pow(phys->R * uN, 0.5); 
	Real Fn = phys->kn *uN; // normal Force (scalar)	
	phys->normalForce = Fn*scg->normal; // normal Force (vector)

	/***************/
	/* SHEAR FORCE */
	/***************/
	
        phys->ks = 8.0 * phys->G *std::pow(phys->R * uN, 0.5);//adapted tangential stiffness
	// 1. Rotate shear force
	Vector3r& shearElastic = scg->rotate(phys->shearElastic);
	// 2. Get shear force (incrementally)
	const Vector3r& shearDisp = scg->shearIncrement();
	shearElastic = shearElastic - phys->ks*(shearDisp);

	
	/********************/
	/* MOHR-COULOMB law */
	/********************/
	phys->isSliding=false;
	phys->shearViscous=Vector3r::Zero(); // reset so that during sliding, the previous values is not there
	// NO ADHESION
    //-------Suhr&Six16---------------------------------------------------------
    //change mu to be pressure dependent!
    Real pm;
    if (Fn>0.0){
        //divide force with contact area (circular), radius: sqrt(phys->R*uN)
        pm=Fn/(phys->R*uN*M_PI); //mean pressure of contact 
    }else{
        pm=0.0;
    }
    Real mu=phys->mu0+phys->c1/(1.0+phys->c2*pm);
    phys->tangensOfFrictionAngle=mu;
		
    Real maxFs = Fn*phys->tangensOfFrictionAngle;

    if (shearElastic.squaredNorm() > maxFs*maxFs){
        phys->isSliding=true;
        Real ratio = maxFs/shearElastic.norm();
        shearElastic *= ratio; phys->shearForce = shearElastic; /*store only elastic shear displacement*/			
    }else {phys->shearForce = shearElastic;} // update the shear force at the elastic value if no damping is present and if we passed MC

	/****************/
	/* APPLY FORCES */
	/****************/
	
	if (!scene->isPeriodic)
		applyForceAtContactPoint(-phys->normalForce - phys->shearForce, scg->contactPoint , id1, de1->se3.position, id2, de2->se3.position);
	else { // in scg we do not wrap particles positions, hence "applyForceAtContactPoint" cannot be used
		Vector3r force = -phys->normalForce - phys->shearForce;
		scene->forces.addForce(id1,force);
		scene->forces.addForce(id2,-force);
		scene->forces.addTorque(id1,(scg->radius1-0.5*scg->penetrationDepth)* scg->normal.cross(force));
		scene->forces.addTorque(id2,(scg->radius2-0.5*scg->penetrationDepth)* scg->normal.cross(force));
	}
	

	return true;
}



} // namespace yade


