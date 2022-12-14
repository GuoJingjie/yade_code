// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2006 Janek Kozicki <cosurgi@berlios.de>

#pragma once

#include <lib/multimethods/FunctorWrapper.hpp>
#include <core/Body.hpp>
#include <core/Bound.hpp>
#include <core/Dispatcher.hpp>
#include <core/Functor.hpp>
#include <core/IGeom.hpp>
#include <core/IPhys.hpp>
#include <core/Interaction.hpp>
#include <core/Shape.hpp>
#include <core/State.hpp>

namespace forCtags { // for ctags, it can't find class inside these macros.
struct GlBoundFunctor {
};
struct GlShapeFunctor {
};
struct GlIGeomFunctor {
};
struct GlIPhysFunctor {
};
struct GlStateFunctor {
};
}

namespace yade { // Cannot have #include directive inside.

#define RENDERS(name)                                                                                                                                          \
public:                                                                                                                                                        \
	string renders() const override { return #name; };                                                                                                     \
	FUNCTOR1D(name);

struct GLViewInfo {
	GLViewInfo()
	        : sceneCenter(Vector3r::Zero())
	        , sceneRadius(1.)
	{
	}
	Vector3r sceneCenter;
	Real     sceneRadius;
};

class OpenGLRenderer;

#define GL_FUNCTOR(Klass, typelist, renderedType)                                                                                                              \
	class Klass : public Functor1D<renderedType, void, typelist> {                                                                                         \
	public:                                                                                                                                                \
		virtual ~Klass() {};                                                                                                                           \
		virtual string renders() const { throw std::runtime_error(#Klass ": unregistered gldraw class.\n"); };                                         \
		virtual void   initgl() { /*WARNING: it must deal with static members, because it is called from another instance!*/ };                        \
		YADE_CLASS_BASE_DOC(Klass, Functor, "Abstract functor for rendering :yref:`" #renderedType "` objects.");                                      \
	};                                                                                                                                                     \
	REGISTER_SERIALIZABLE(Klass);
#define GL_DISPATCHER(Klass, Functor)                                                                                                                          \
	class Klass : public Dispatcher1D<Functor> {                                                                                                           \
	public:                                                                                                                                                \
		YADE_DISPATCHER1D_FUNCTOR_DOC_ATTRS_CTOR_PY(Klass, Functor, /*optional doc*/, /*attrs*/, /*ctor*/, /*py*/);                                    \
	};                                                                                                                                                     \
	REGISTER_SERIALIZABLE(Klass);

GL_FUNCTOR(GlBoundFunctor, TYPELIST_2(const shared_ptr<Bound>&, Scene*), Bound);
GL_FUNCTOR(GlShapeFunctor, TYPELIST_4(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&), Shape);
GL_FUNCTOR(GlIGeomFunctor, TYPELIST_5(const shared_ptr<IGeom>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool), IGeom);
GL_FUNCTOR(GlIPhysFunctor, TYPELIST_5(const shared_ptr<IPhys>&, const shared_ptr<Interaction>&, const shared_ptr<Body>&, const shared_ptr<Body>&, bool), IPhys);
GL_FUNCTOR(GlStateFunctor, TYPELIST_1(const shared_ptr<State>&), State);

GL_DISPATCHER(GlBoundDispatcher, GlBoundFunctor);
GL_DISPATCHER(GlShapeDispatcher, GlShapeFunctor);
GL_DISPATCHER(GlIGeomDispatcher, GlIGeomFunctor);
GL_DISPATCHER(GlIPhysDispatcher, GlIPhysFunctor);
GL_DISPATCHER(GlStateDispatcher, GlStateFunctor);
#undef GL_FUNCTOR
#undef GL_DISPATCHER

} // namespace yade
