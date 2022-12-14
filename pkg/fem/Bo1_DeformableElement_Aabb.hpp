/*************************************************************************
*  Copyright (C) 2013 by Burak ER                                 	 *
*									 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <core/Dispatching.hpp>
#include <pkg/fem/DeformableElement.hpp>

namespace yade { // Cannot have #include directive inside.

class Bo1_DeformableElement_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r&, const Body*) override;
	FUNCTOR1D(DeformableElement);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(Bo1_DeformableElement_Aabb,BoundFunctor,"Functor creating :yref:`Aabb` from :yref:`DeformableElement`.",
		((Real,aabbEnlargeFactor,((void)"deactivated",-1),,"Relative enlargement of the bounding box; deactivated if negative.\n\n.. note::\n\tThis attribute is used to create distant interaction, but is only meaningful with an :yref:`IGeomFunctor` which will not simply discard such interactions: :yref:`Ig2_Sphere_Sphere_ScGeom::interactionDetectionFactor` should have the same value as :yref:`aabbEnlargeFactor<Bo1_Sphere_Aabb::aabbEnlargeFactor>`."))
	);
	// clang-format on
};

REGISTER_SERIALIZABLE(Bo1_DeformableElement_Aabb);

} // namespace yade
