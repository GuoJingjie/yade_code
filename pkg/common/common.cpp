// =======================================================
// Some plugins from removed CPP-fiels

#include <core/Aabb.hpp>
#include <pkg/common/BoundaryController.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/CylScGeom6D.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/FieldApplier.hpp>
#include <pkg/common/ForceResetter.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/common/PyRunner.hpp>
#include <pkg/common/Recorder.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/common/StepDisplacer.hpp>
#include <pkg/common/TorqueEngine.hpp>
#include <pkg/dem/DemXDofGeom.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((ForceResetter)(TorqueEngine)(FieldApplier)(BoundaryController)(NormPhys)(NormShearPhys)(Recorder)(CylScGeom6D)(CylScGeom)(Box)(StepDisplacer)(
        GenericSpheresContact)(PeriodicEngine)(Sphere)(ElastMat)(FrictMat)(PyRunner));

} // namespace yade
