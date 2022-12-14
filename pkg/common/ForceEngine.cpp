// 2004 © Janek Kozicki <cosurgi@berlios.de>
// 2009 © Václav Šmilauer <eudoxos@arcig.cz>

#include "ForceEngine.hpp"
#include <lib/high-precision/Constants.hpp>
#include <lib/smoothing/LinearInterpolate.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <preprocessing/dem/Shop.hpp>

#include <core/IGeom.hpp>
#include <core/IPhys.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((ForceEngine)(InterpolatingDirectedForceEngine)(RadialForceEngine)(DragEngine)(LinearDragEngine)(HarmonicForceEngine));

void ForceEngine::action()
{
	FOREACH(Body::id_t id, ids)
	{
		if (!(scene->bodies->exists(id))) continue;
		scene->forces.addForce(id, force);
	}
}

void InterpolatingDirectedForceEngine::action()
{
	Real virtTime = wrap ? Shop::periodicWrap(scene->time, *times.begin(), *times.rbegin()) : scene->time;
	direction.normalize();
	force = linearInterpolate<Real, Real>(virtTime, times, magnitudes, _pos) * direction;
	ForceEngine::action();
}

void RadialForceEngine::postLoad(RadialForceEngine&) { axisDir.normalize(); }

void RadialForceEngine::action()
{
	FOREACH(Body::id_t id, ids)
	{
		if (!(scene->bodies->exists(id))) continue;
		const Vector3r& pos    = Body::byId(id, scene)->state->pos;
		Vector3r        radial = (pos - (axisPt + axisDir * /* t */ ((pos - axisPt).dot(axisDir)))).normalized();
		if (radial.squaredNorm() == 0) continue;
		scene->forces.addForce(id, fNorm * radial);
	}
}

void DragEngine::action()
{
	FOREACH(Body::id_t id, ids)
	{
		Body* b = Body::byId(id, scene).get();
		if (!b) continue;
		if (!(scene->bodies->exists(id))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere) {
			Real     A          = sphere->radius * sphere->radius * Mathr::PI; //Crossection of the sphere
			Vector3r velSphTemp = Vector3r::Zero();
			Vector3r dragForce  = Vector3r::Zero();

			if (scene->isPeriodic) {
				velSphTemp = scene->cell->bodyFluctuationVel(b->state->pos, b->state->vel, scene->cell->prevVelGrad);
			} else {
				velSphTemp = b->state->vel;
			}

			if (velSphTemp != Vector3r::Zero()) { dragForce = -0.5 * Rho * A * Cd * velSphTemp.squaredNorm() * velSphTemp.normalized(); }
			scene->forces.addForce(id, dragForce);
		}
	}
}

void LinearDragEngine::action()
{
	FOREACH(Body::id_t id, ids)
	{
		Body* b = Body::byId(id, scene).get();
		if (!b) continue;
		if (!(scene->bodies->exists(id))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere) {
			Vector3r velSphTemp = Vector3r::Zero();

			if (scene->isPeriodic) {
				velSphTemp = scene->cell->bodyFluctuationVel(b->state->pos, b->state->vel, scene->cell->prevVelGrad);
			} else {
				velSphTemp = b->state->vel;
			}

			Vector3r dragForce = Vector3r::Zero();

			Real b2 = 6. * Mathr::PI * nu * sphere->radius;

			if (velSphTemp != Vector3r::Zero()) { dragForce = -b2 * velSphTemp; }
			scene->forces.addForce(id, dragForce);
		}
	}
}

void HarmonicForceEngine::action()
{
	if (ids.size() > 0) {
		Vector3r w = f * 2.0 * Mathr::PI; //Angular frequency

		Vector3r force = (((w * scene->time + fi).array().sin()));

		force = force.cwiseProduct(A);

		FOREACH(Body::id_t id, ids)
		{
			assert(id < (Body::id_t)scene->bodies->size());
			Body* b = Body::byId(id, scene).get();
			if (!b) continue;
			if (!(scene->bodies->exists(id))) continue;
			scene->forces.addForce(id, force);
		}
	} else {
		LOG_WARN("The list of ids is empty! Can't apply any forces.");
	}
}

} // namespace yade
