//keep this #ifdef as long as you don't really want to realize a final version publicly, it will save compilation time for everyone else
//when you want it compiled, you can just uncomment the following line

#ifdef YADE_CGAL
#define CAPILLARYPHYS1
#ifdef CAPILLARYPHYS1

#include <core/Omega.hpp>
#include <core/Scene.hpp>
#include <pkg/dem/CapillaryPhysDelaunay.hpp>
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN(
        (CapillaryPhysDelaunay)(CapillaryMindlinPhysDelaunay)(Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay)(Ip2_FrictMat_FrictMat_CapillaryMindlinPhysDelaunay));

void Ip2_FrictMat_FrictMat_CapillaryPhysDelaunay::go(
        const shared_ptr<Material>& b1 //FrictMat
        ,
        const shared_ptr<Material>& b2 // FrictMat
        ,
        const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;
	Ip2_FrictMat_FrictMat_FrictPhys::go(b1, b2, interaction);
	if (interaction->phys) {
		auto newPhys           = shared_ptr<CapillaryPhysDelaunay>(new CapillaryPhysDelaunay(*YADE_PTR_CAST<FrictPhys>(interaction->phys)));
		newPhys->computeBridge = computeDefault;
		interaction->phys      = newPhys;
	}
};

} // namespace yade

#endif //CAPILLARYPHYS1
#endif //YADE_CGAL
