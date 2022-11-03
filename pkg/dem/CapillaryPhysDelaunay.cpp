//keep this #ifdef as long as you don't really want to realize a final version publicly, it will save compilation time for everyone else
//when you want it compiled, you can just uncomment the following line

#ifdef YADE_CGAL
#define CAPILLARYPHYS1
#ifdef CAPILLARYPHYS1

#include <pkg/dem/CapillaryPhysDelaunay.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <core/Omega.hpp>
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

// CapillaryPhysDelaunay::~CapillaryPhysDelaunay() { }
// CapillaryMindlinPhysDelaunay::~CapillaryMindlinPhysDelaunay() { }
YADE_PLUGIN((CapillaryPhysDelaunay)(CapillaryMindlinPhysDelaunay)(Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay)(Ip2_FrictMat_FrictMat_CapillaryMindlinPhysDelaunay));

// void Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay::go(
//         const shared_ptr<Material>& b1 //FrictMat
//         ,
//         const shared_ptr<Material>& b2 // FrictMat
//         ,
//         const shared_ptr<Interaction>& interaction)
// {
// 	ScGeom* geom = YADE_CAST<ScGeom*>(interaction->geom.get());
// 	if (geom) {
// 		if (!interaction->phys) {
// 			const shared_ptr<FrictMat>& sdec1 = YADE_PTR_CAST<FrictMat>(b1);
// 			const shared_ptr<FrictMat>& sdec2 = YADE_PTR_CAST<FrictMat>(b2);
// 
// 			if (!interaction->phys) interaction->phys = shared_ptr<CapillaryPhysDelaunay>(new CapillaryPhysDelaunay());
// 			const shared_ptr<CapillaryPhysDelaunay>& contactPhysics = YADE_PTR_CAST<CapillaryPhysDelaunay>(interaction->phys);
// 
// 			Real Ea = sdec1->young;
// 			Real Eb = sdec2->young;
// 			Real Va = sdec1->poisson;
// 			Real Vbb = sdec2->poisson;
// 			Real Da = geom->radius1; // FIXME - multiply by factor of sphere interaction distance (so sphere interacts at bigger range that its geometrical size)
// 			Real Db = geom->radius2; // FIXME - as above
// 			Real fa = sdec1->frictionAngle;
// 			Real fb = sdec2->frictionAngle;
// 			Real Kn = 2 * Ea * Da * Eb * Db / (Ea * Da + Eb * Db); //harmonic average of two stiffnesses
// 			Real Ks = 2 * Ea * Da * Va * Eb * Db * Vbb
// 			        / (Ea * Da * Va + Eb * Db * Va); //harmonic average of two stiffnesses with ks=V*kn for each sphere
// 
// 			contactPhysics->tangensOfFrictionAngle = math::tan(math::min(fa, fb));
// 			contactPhysics->kn                     = Kn;
// 			contactPhysics->ks                     = Ks;
// 			contactPhysics->computeBridge          = computeDefault;
// 		}
// 	}
// };


void Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay::go(
        const shared_ptr<Material>& b1 //FrictMat
        ,
        const shared_ptr<Material>& b2 // FrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;
	Ip2_FrictMat_FrictMat_FrictPhys::go(b1,b2,interaction);
	if (interaction->phys) {
		auto newPhys = shared_ptr<CapillaryPhysDelaunay>(new CapillaryPhysDelaunay(*YADE_PTR_CAST<FrictPhys>(interaction->phys)));
		newPhys->computeBridge          = computeDefault;
		interaction->phys = newPhys;
	}
};



void Ip2_FrictMat_FrictMat_CapillaryMindlinPhysDelaunay::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;
	Ip2_FrictMat_FrictMat_MindlinPhys::go(b1,b2,interaction);
	if (interaction->phys) {
		auto newPhys = shared_ptr<CapillaryMindlinPhysDelaunay>(new CapillaryMindlinPhysDelaunay(*YADE_PTR_CAST<MindlinPhys>(interaction->phys)));
		newPhys->computeBridge          = computeDefault;
		interaction->phys = newPhys;
	}
};

} // namespace yade

#endif //CAPILLARYPHYS1
#endif //YADE_CGAL
