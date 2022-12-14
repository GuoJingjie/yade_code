#include "FrictPhys.hpp"
#include <pkg/dem/ScGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((FrictPhys)(ViscoFrictPhys)(RotStiffFrictPhys)(Ip2_FrictMat_FrictMat_ViscoFrictPhys)(Ip2_FrictMat_FrictMat_FrictPhys));

// The following code was moved from Ip2_FrictMat_FrictMat_FrictPhys.hpp

void Ip2_FrictMat_FrictMat_FrictPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;

	const shared_ptr<FrictMat>& mat1 = YADE_PTR_CAST<FrictMat>(b1);
	const shared_ptr<FrictMat>& mat2 = YADE_PTR_CAST<FrictMat>(b2);

	Real Ra, Rb;                                                           //Vector3r normal;
	assert(dynamic_cast<GenericSpheresContact*>(interaction->geom.get())); //only in debug mode
	GenericSpheresContact* sphCont = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	Ra                             = sphCont->refR1 > 0 ? sphCont->refR1 : sphCont->refR2;
	Rb                             = sphCont->refR2 > 0 ? sphCont->refR2 : sphCont->refR1;

	interaction->phys                           = shared_ptr<FrictPhys>(new FrictPhys());
	const shared_ptr<FrictPhys>& contactPhysics = YADE_PTR_CAST<FrictPhys>(interaction->phys);
	Real                         Ea             = mat1->young;
	Real                         Eb             = mat2->young;
	Real                         Va             = mat1->poisson;
	Real                         Vb             = mat2->poisson;
	Real                         kna            = Ea * Ra;
	Real                         knb            = Eb * Rb;
	Real                         ksa            = kna * Va;
	Real                         ksb            = knb * Vb;

	//match maker or half the harmonic average of the two stiffnesses, when (2*Ri*Ei=2*kni) is the stiffness of a contact point on sphere "i"
	Real Kn = (!kn) ? 2 * kna * knb / (kna + knb) : (*kn)(mat1->id, mat2->id, kna, knb);
	//same for shear stiffness
	Real Ks = (!ks) ? 2 * ksa * ksb / (ksa + ksb) : (*ks)(mat1->id, mat2->id, ksa, ksb);

	Real frictionAngle                     = (!frictAngle) ? math::min(mat1->frictionAngle, mat2->frictionAngle)
	                                                       : (*frictAngle)(mat1->id, mat2->id, mat1->frictionAngle, mat2->frictionAngle);
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	contactPhysics->kn                     = Kn;
	contactPhysics->ks                     = Ks;
};

void Ip2_FrictMat_FrictMat_ViscoFrictPhys::go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction)
{
	if (interaction->phys) return;
	const shared_ptr<FrictMat>& mat1                 = YADE_PTR_CAST<FrictMat>(b1);
	const shared_ptr<FrictMat>& mat2                 = YADE_PTR_CAST<FrictMat>(b2);
	interaction->phys                                = shared_ptr<ViscoFrictPhys>(new ViscoFrictPhys());
	const shared_ptr<ViscoFrictPhys>& contactPhysics = YADE_PTR_CAST<ViscoFrictPhys>(interaction->phys);
	Real                              Ea             = mat1->young;
	Real                              Eb             = mat2->young;
	Real                              Va             = mat1->poisson;
	Real                              Vb             = mat2->poisson;

	Real Ra, Rb;                                                           //Vector3r normal;
	assert(dynamic_cast<GenericSpheresContact*>(interaction->geom.get())); //only in debug mode
	GenericSpheresContact* sphCont = YADE_CAST<GenericSpheresContact*>(interaction->geom.get());
	Ra                             = sphCont->refR1 > 0 ? sphCont->refR1 : sphCont->refR2;
	Rb                             = sphCont->refR2 > 0 ? sphCont->refR2 : sphCont->refR1;

	//harmonic average of the two stiffnesses when (Ri.Ei/2) is the stiffness of a contact point on sphere "i"
	Real Kn = 2 * Ea * Ra * Eb * Rb / (Ea * Ra + Eb * Rb);
	//same for shear stiffness
	Real Ks = 2 * Ea * Ra * Va * Eb * Rb * Vb / (Ea * Ra * Va + Eb * Rb * Vb);

	Real frictionAngle                     = (!frictAngle) ? math::min(mat1->frictionAngle, mat2->frictionAngle)
	                                                       : (*frictAngle)(mat1->id, mat2->id, mat1->frictionAngle, mat2->frictionAngle);
	contactPhysics->tangensOfFrictionAngle = math::tan(frictionAngle);
	contactPhysics->kn                     = Kn;
	contactPhysics->ks                     = Ks;
};

Vector3r RotStiffFrictPhys::getRotStiffness() const { return Vector3r(ktw, kr, kr); };

} // namespace yade
