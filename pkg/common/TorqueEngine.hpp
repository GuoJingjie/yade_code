/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <core/PartialEngine.hpp>
#include <core/Scene.hpp>

namespace yade { // Cannot have #include directive inside.

class TorqueEngine : public PartialEngine {
public:
	void action() override
	{
		FOREACH(const Body::id_t id, ids)
		{
			// check that body really exists?
			scene->forces.addTorque(id, moment);
		}
	}
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(TorqueEngine,PartialEngine,"Apply given torque (momentum) value at every subscribed particle, at every step.",
		((Vector3r,moment,Vector3r::Zero(),,"Torque value to be applied."))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(TorqueEngine);

} // namespace yade
