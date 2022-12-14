#pragma once

#include <core/GlobalEngine.hpp>
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

class Scene;
class ForceResetter : public GlobalEngine {
public:
	void action() override
	{
		scene->forces.reset(scene->iter);
		if (scene->trackEnergy) scene->energy->resetResettables();
	}
	// clang-format off
	YADE_CLASS_BASE_DOC(ForceResetter,GlobalEngine,"Reset all forces stored in Scene::forces (``O.forces`` in python). Typically, this is the first engine to be run at every step. In addition, reset those energies that should be reset, if energy tracing is enabled.");
	// clang-format on
};
REGISTER_SERIALIZABLE(ForceResetter);

} // namespace yade
