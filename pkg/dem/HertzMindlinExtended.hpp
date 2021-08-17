/*
 2016 - Bettina Suhr 
 extension of Hertz-Mindlin model (see HertzMindlin.cpp): 

 normal direction: conical damage model (Harkness et al 2016, Suhr & Six 2017)
 tangential direction: stress dependent interparticle friction coefficient (Suhr & Six 2016) 
 both models can be switched on/off separately

 references:
 Harkness, Zervos, Le Pen, Aingaran, Powrie: Discrete element simulation of railway ballast: modelling cell pressure effects in triaxial tests, Granular Matter, (2016) 18:65
 Suhr & Six 2017: Parametrisation of a DEM model for railway ballast under different load cases, Granular Matter, (2017) 19:64
 Suhr & Six 2016: On the effect of stress dependent interparticle friction in direct shear tests, Powder Technology , (2016) 294:211 - 220
*/

#pragma once

#include <core/Dispatching.hpp>
#include <core/Material.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/dem/HertzMindlin.hpp>
#include <pkg/dem/ScGeom.hpp>

#include <lib/base/openmp-accu.hpp>
#include <boost/tuple/tuple.hpp>

namespace yade { // Cannot have #include directive inside.

class FrictMatCDM : public FrictMat {
public:
	virtual ~FrictMatCDM() {};
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(
	        FrictMatCDM,
	        FrictMat,
	        "Material to be used for extended Hertz-Mindlin contact law. Normal direction:  parameters for Conical Damage Model  (Harkness et al. 2016, "
	        "Suhr & Six 2017). Tangential direction: parameters for stress dependent interparticle friction coefficient (Suhr & Six 2016). Both models can "
	        "be switched on/off separately. ",
	        ((Real, sigmaMax, 1e99, , " >0 [Pa] max compressive strength of material, choose 1e99 to switch off conical damage model"))(
	                (Real, alpha, 1e-6, , "[rad] angle of conical asperities, alpha in (0, pi/2)"))
	        //((Real,mu0,0.15,,"[-] parameter of pressure dependent friction model mu0"))
	        ((Real, c1, 0.0, , "[-] parameter of pressure dependent friction model c1, choose 0 for constant interparticle friction coefficient"))(
	                (Real, c2, 0.0, , "[-] parameter of pressure dependent friction model c2, choose 0 for constant interparticle friction coefficient")),
	        createIndex(););
	REGISTER_CLASS_INDEX(FrictMatCDM, FrictMat);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(FrictMatCDM);


class MindlinPhysCDM : public MindlinPhys {
public:
	virtual ~MindlinPhysCDM() = default; //{};
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(
	        MindlinPhysCDM,
	        MindlinPhys,
	        "Representation of an interaction of an extended  Hertz-Mindlin type. Normal direction:  parameters for Conical Damage Model  (Harkness et al. "
	        "2016, Suhr & Six 2017). Tangential direction: parameters for stress dependent interparticle friction coefficient (Suhr & Six 2016). Both "
	        "models can be switched on/off separately, see FrictMatCDM.",
	        //((Real,tangensOfFrictionAngle,NaN,,"tan of angle of friction"))
	        ((Real, E, 0.0, , " [Pa] equiv. Young's modulus"))((Real, G, 0.0, , " [Pa] equiv. shear modulus"))(
	                (Real, sigmaMax, 0.0, , " [Pa] max compressive strength of material"))(
	                (Real, alphaFac, 0.0, , "factor considering angle of conical asperities"))(
	                (Real, R, 0.0, , "[m] contact radius in conical damage model"))((bool, isYielding, false, , "bool: is contact currently yielding?"))(
	                (Real, mu0, 0.0, , "[-] parameter of pressure dependent friction model mu0"))(
	                (Real, c1, 0.0, , "[-] parameter of pressure dependent friction model c1"))(
	                (Real, c2, 0.0, , "[-] parameter of pressure dependent friction model c2")),
	        createIndex());
	REGISTER_CLASS_INDEX(MindlinPhysCDM, MindlinPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(MindlinPhysCDM);


class Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	FUNCTOR2D(FrictMatCDM, FrictMatCDM);
	YADE_CLASS_BASE_DOC_ATTRS(
	        Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM,
	        IPhysFunctor,
	        "Create a :yref:`MindlinPhysCDM` from two :yref:`FrictMatCDMsExts<FrictMatCDMExt>`. ",
	        ((shared_ptr<MatchMaker>,
	          frictAngle,
	          ,
	          ,
	          "Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used.")));
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM);


class Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM : public IPhysFunctor {
public:
	void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;
	FUNCTOR2D(FrictMat, FrictMatCDM);
	YADE_CLASS_BASE_DOC_ATTRS(
	        Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM,
	        IPhysFunctor,
	        "Create a :yref:`MindlinPhysCDM` from one FrictMat and one FrictMatCDM instance. ",
	        ((shared_ptr<MatchMaker>,
	          frictAngle,
	          ,
	          ,
	          "Instance of :yref:`MatchMaker` determining how to compute interaction's friction angle. If ``None``, minimum value is used.")));
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM);


/******************** Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM *********/
class Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	Real ratioSlidingContacts();
	Real ratioYieldingContacts();

	FUNCTOR2D(ScGeom, MindlinPhysCDM);
	YADE_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(
	        Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM,
	        LawFunctor,
	        "Hertz-Mindlin model extended: Normal direction: conical damage model from Harkness et al. 2016./ Suhr & Six 2017. Tangential direction: "
	        "stress dependent interparticle friction coefficient, Suhr & Six 2016. Both models can be switched on/off separately. In this version there is "
	        "NO damping (neither viscous nor linear), NO adhesion and NO calc_energy, NO includeMoment, NO preventGranularRatcheting. NOT tested for "
	        "periodic simulations.",
	        ((bool,
	          neverErase,
	          false,
	          ,
	          "Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. "
	          ":yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))

	                , /*deprec*/
	        ,         /* init */
	        ,         /* ctor */
	        ,         /* py */
	        .def("ratioSlidingContacts",
	             &Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM::ratioSlidingContacts,
	             "Return the ratio between the number of contacts sliding to the total number at a given time.")
	                .def("ratioYieldingContacts",
	                     &Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM::ratioYieldingContacts,
	                     "Return the ratio between the number of contacts yielding to the total number at a given time."));
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM);


} // namespace yade
