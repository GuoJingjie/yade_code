#pragma once
#include <lib/base/Math.hpp>
#include <core/Aabb.hpp>
#include <core/Dispatcher.hpp>
#include <core/Functor.hpp>
#include <core/IGeom.hpp>
#include <core/IPhys.hpp>
#include <core/Interaction.hpp>
#include <core/Scene.hpp>
#include <core/Shape.hpp>
#include <core/State.hpp>

namespace yade { // Cannot have #include directive inside.

/********
	functors
*********************/

class BoundFunctor : public Functor1D<
                             /*dispatch types*/ Shape,
                             /*return type*/ void,
                             /*argument types*/ TYPELIST_4(const shared_ptr<Shape>&, shared_ptr<Bound>&, const Se3r&, const Body*)> {
public:
	virtual ~BoundFunctor();
	// clang-format off
	YADE_CLASS_BASE_DOC(BoundFunctor,Functor,"Functor for creating/updating :yref:`Body::bound`.");
	// clang-format on
};
REGISTER_SERIALIZABLE(BoundFunctor);

class IGeomFunctor : public Functor2D<
                             /*dispatch types*/ Shape,
                             Shape,
                             /*return type*/ bool,
                             /*argument types*/
                             TYPELIST_7(
                                     const shared_ptr<Shape>&,
                                     const shared_ptr<Shape>&,
                                     const State&,
                                     const State&,
                                     const Vector3r&,
                                     const bool&,
                                     const shared_ptr<Interaction>&)> {
public:
	virtual ~IGeomFunctor();
	// called before every step once, from InteractionLoop (used to set Scene::flags & Scene::LOCAL_COORDS)
	virtual void preStep() {};
	// clang-format off
	YADE_CLASS_BASE_DOC(IGeomFunctor,Functor,"Functor for creating/updating :yref:`Interaction::geom` objects.");
	// clang-format on
};
REGISTER_SERIALIZABLE(IGeomFunctor);


class IPhysFunctor : public Functor2D<
                             /*dispatch types*/ Material,
                             Material,
                             /*retrun type*/ void,
                             /*argument types*/ TYPELIST_3(const shared_ptr<Material>&, const shared_ptr<Material>&, const shared_ptr<Interaction>&)> {
public:
	virtual ~IPhysFunctor();
	// clang-format off
	YADE_CLASS_BASE_DOC(IPhysFunctor,Functor,"Functor for creating/updating :yref:`Interaction::phys` objects from :yref:`bodies' material<Body::material>` properties.");
	// clang-format on
};
REGISTER_SERIALIZABLE(IPhysFunctor);


class LawFunctor : public Functor2D<
                           /*dispatch types*/ IGeom,
                           IPhys,
                           /*return type*/ bool,
                           /*argument types*/ TYPELIST_3(shared_ptr<IGeom>&, shared_ptr<IPhys>&, Interaction*)> {
public:
	virtual ~LawFunctor();
	/*! Convenience functions to get forces/torques quickly. */
	void addForce(const Body::id_t id, const Vector3r& f, Scene* rb) { rb->forces.addForce(id, f); }
	void addTorque(const Body::id_t id, const Vector3r& t, Scene* rb) { rb->forces.addTorque(id, t); }
	/*! Convenience function to apply force and torque from one force at contact point. Not sure if this is the right place for it. */
	void applyForceAtContactPoint(
	        const Vector3r& force, const Vector3r& contactPoint, const Body::id_t id1, const Vector3r& pos1, const Body::id_t id2, const Vector3r& pos2)
	{
		addForce(id1, force, scene);
		addTorque(id1, (contactPoint - pos1).cross(force), scene);
		addForce(id2, -force, scene);
		addTorque(id2, -(contactPoint - pos2).cross(force), scene);
	}
	// clang-format off
	YADE_CLASS_BASE_DOC(LawFunctor,Functor,"Functor for applying constitutive laws on :yref:`interactions<Interaction>`.");
	// clang-format on
};
REGISTER_SERIALIZABLE(LawFunctor);


/********
	dispatchers
*********************/

class BoundDispatcher : public Dispatcher1D<
                                /* functor type*/ BoundFunctor> {
public:
	void action() override;
	bool isActivated() override { return activated; }
	void processBody(const shared_ptr<Body>&);
	DECLARE_LOGGER;
	YADE_DISPATCHER1D_FUNCTOR_DOC_ATTRS_CTOR_PY(
	        BoundDispatcher,
	        BoundFunctor,
	        /*optional doc*/,
	        /*additional attrs*/
	        ((bool, activated, true, , "Whether the engine is activated (only should be changed by the collider)"))(
	                (Real,
	                 sweepDist,
	                 0,
	                 ,
	                 "Distance by which enlarge all bounding boxes, to prevent collider from being run at every step (only should be changed by the "
	                 "collider)."))(
	                (Real,
	                 minSweepDistFactor,
	                 0.2,
	                 ,
	                 "Minimal distance by which enlarge all bounding boxes; superseeds computed value of sweepDist when lower that (minSweepDistFactor x "
	                 "sweepDist). Updated by the collider. |yupdate|."))(
	                (Real, updatingDispFactor, -1, , "see :yref:`InsertionSortCollider::updatingDispFactor` |yupdate|"))(
	                (Real, targetInterv, -1, , "see :yref:`InsertionSortCollider::targetInterv` |yupdate|")),
	        /*ctor*/, /*py*/
	);
};
REGISTER_SERIALIZABLE(BoundDispatcher);


class IGeomDispatcher : public Dispatcher2D<
                                /* functor type*/ IGeomFunctor,
                                /* autosymmetry*/ false> {
	bool alreadyWarnedNoCollider;

public:
	void                    action() override;
	shared_ptr<Interaction> explicitAction(const shared_ptr<Body>& b1, const shared_ptr<Body>& b2, bool force);
	YADE_DISPATCHER2D_FUNCTOR_DOC_ATTRS_CTOR_PY(IGeomDispatcher, IGeomFunctor, /* doc is optional*/, /*attrs*/, /*ctor*/ alreadyWarnedNoCollider = false;,
	                                            /*py*/);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(IGeomDispatcher);


class IPhysDispatcher : public Dispatcher2D<
                                /*functor type*/ IPhysFunctor> {
public:
	void action() override;
	void explicitAction(shared_ptr<Material>& pp1, shared_ptr<Material>& pp2, shared_ptr<Interaction>& i);
	YADE_DISPATCHER2D_FUNCTOR_DOC_ATTRS_CTOR_PY(IPhysDispatcher, IPhysFunctor, /*doc is optional*/, /*attrs*/, /*ctor*/, /*py*/);
};
REGISTER_SERIALIZABLE(IPhysDispatcher);


class LawDispatcher : public Dispatcher2D<
                              /*functor type*/ LawFunctor,
                              /*autosymmetry*/ false> {
public:
	void action() override;
	YADE_DISPATCHER2D_FUNCTOR_DOC_ATTRS_CTOR_PY(LawDispatcher, LawFunctor, /*doc is optional*/, /*attrs*/, /*ctor*/, /*py*/);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(LawDispatcher);


} // namespace yade
