// © 2021-2022 Anton Gladky
// © 2021-2022 Karol Brzeziński <brzezink.ubuntu@gmail.com>


#pragma once

#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <pkg/dem/Ig2_Facet_Sphere_ScGeom.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <unordered_map>


namespace yade {

//*********************************************************************************
/* Clump breakage algorithm */
class ClumpBreakage : public PeriodicEngine {
private:
	struct bodyProp {
		Matrix3r stress = Matrix3r::Zero();
		Real     v      = 0; // volume
		Real     effort = 0; // effort
		bodyProp()      = default;
		bodyProp(Matrix3r _stress, Real _v)
		        : stress(_stress)
		        , v(_v) {};
	};

	// using MapBody2Stress = std::unordered_map<Body::id_t, Matrix3r>;
	// using MapBody2StressTuple = std::unordered_map<Body::id_t, std::tuple<Matrix3r, Real, Real>>; // <stress, volume, effort
	using MapBody2StressTuple = std::unordered_map<Body::id_t, bodyProp>; // <stress, volume, effort

public:
	void action() override;

	MapBody2StressTuple getStressForEachBody();
	void                checkFailure(MapBody2StressTuple&);
	bool                replaceSphere(
	                       MapBody2StressTuple&,
	                       Body::id_t sphere_id              = 0,
	                       Real       subparticles_mass      = 0.,
	                       Real       radius_ratio           = 2,
	                       Real       relative_gap           = 0,
	                       Real       initial_packing_scale  = 1.5,
	                       Real       max_scale              = 3,
	                       bool       search_for_neighbourse = true,
	                       Real       grow_radius            = 1.0,
	                       Real       max_grow_radius        = 2.0,
	                       Vector3r   shift_vector           = Vector3r::Zero());
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(
		ClumpBreakage,PeriodicEngine,"Engine that splits clumps.\n\n.."
		,
		((bool,stress_correction,true,,"boolean controling the applying of a stress correction in a stress calculation"))
		((bool,weibull,true,,"boolean controlling the usage of weibull formulation"))
   		((Real,compressive_strength,-1,,"Stress at which polyhedra of volume 4/3*pi [mm] breaks."))
		((Real,tension_strength,-1,,"Tangential stress at which polyhedra of volume 4/3*pi [mm] breaks."))
		((int,Wei_m,3,,"Weibull Formulation, Weibull modulus, m, (if negative - disabled), [Gladky2017]_"))
		((Real,Wei_V0,0.01,,"Weibull Formulation, V0, m^3, representative volume, [Gladky2017]_."))
		((Real,Wei_P,0.63,,"Weibull Formulation, failure  probability, P, [Gladky2017]_."))
		((Real,young,1e8,,"Young modulus"))
        ,
		/*ctor*/
		,/*py*/
	);
	// clang-format on


	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ClumpBreakage);

} // namespace yade
