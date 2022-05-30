// © 2021-2022 Anton Gladky
// © 2021-2022 Karol Brzeziński <brzezink.ubuntu@gmail.com>

#include <lib/high-precision/Constants.hpp>
#include <pkg/dem/ClumpBreakage.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade { // Cannot have #include directive inside.

YADE_PLUGIN((ClumpBreakage));

CREATE_LOGGER(ClumpBreakage);


ClumpBreakage::MapBody2StressTuple ClumpBreakage::getStressForEachBody()
{
	MapBody2StressTuple bStress;
	YADE_PARALLEL_FOREACH_BODY_BEGIN(const shared_ptr<Body>& b, scene->bodies)
	{
		// skip clumps, only apply forces on their constituents
		if (b->isClump()) continue;
		Vector3r f_u    = Vector3r::Zero(); // sum of all normal forces
		Vector3r t_u    = Vector3r::Zero(); // sum of all tangential forces
		Matrix3r stress = Matrix3r::Zero();
		for (const auto& i : b->intrs) {
			const auto& I = i.second;
			if (!I->isReal()) continue;
			const auto* geom = YADE_CAST<ScGeom*>(I->geom.get());
			const auto* phys = YADE_CAST<FrictPhys*>(I->phys.get());
			if (!geom || !phys) continue;

			if (b->isClumpMember()) {
				// for the clump member we check whether the partner body belongs to the same clump
				// if yes - we skip it
				Body::id_t otherId = -1;
				if (b->id == I->getId1()) {
					otherId = I->getId2();
				} else {
					otherId = I->getId1();
				}
				const auto& otherBody = Body::byId(otherId, scene);
				if (otherBody->isClump() and otherBody->clumpId == b->clumpId) continue;
			}
			const auto pp             = geom->contactPoint;
			const auto center         = b->state->pos;
			const auto ss             = phys->shearForce;
			const auto nn             = phys->normalForce;
			const auto ff             = nn + ss;
			const auto contact_vector = pp - center;
			f_u += ff;
			t_u += contact_vector.cross(ss);
			stress += ff * (contact_vector.transpose());
		}

		if (!b || !b->material || !b->shape || (b->shape->getClassName() != "Sphere")) continue;

		// Calculate particle volume
		const auto s = dynamic_cast<Sphere*>(b->shape.get());
		const auto v = 4.0 / 3.0 * M_PI * math::pow(s->radius, 3);

		if (stress_correction) {
			const auto     f_u_prim    = -f_u;
			const auto     normal_prim = f_u.normalized();
			const Vector3r f_t_prim    = (0.5 * t_u).cross(normal_prim) / s->radius;
			const Vector3r f_t_bis     = (0.5 * t_u).cross(-normal_prim) / s->radius;

			stress += (f_u_prim + f_t_prim) * ((s->radius * normal_prim).transpose());
			stress += (f_t_bis) * ((-1 * s->radius * normal_prim).transpose());
		}
		stress /= v;

		// std::cout << "CPP " << b->id << ": stress in getStress \n" << stress << std::endl;

#ifdef YADE_OPENMP
#pragma omp critical
#endif
		{
			bStress[b->id] = ClumpBreakage::bodyProp(stress, v);
		}
	}

	YADE_PARALLEL_FOREACH_BODY_END();
	return bStress;
}

void ClumpBreakage::checkFailure(MapBody2StressTuple& bStress)
{
	for (auto& i : bStress) {
		// const auto bId    = i.first;
		auto       stress = i.second.stress;
		const auto v      = i.second.v;

		if (weibull) {
			const auto coeff = math::pow((-Wei_V0 / v) * math::log(1 - Wei_P), (1.0 / Wei_m));
			stress /= coeff;
		}

		const auto        eivals = stress.eigenvalues();
		std::vector<Real> sigmas { std::real(eivals[0]), std::real(eivals[1]), std::real(eivals[2]) };
		std::sort(sigmas.begin(), sigmas.end());
		const auto sigma_1 = sigmas[2];
		// const auto sigma_2   = sigmas[1];
		const auto sigma_3   = sigmas[0];
		const auto sigma_tau = sigma_1 - sigma_3 * tension_strength / compressive_strength;

		Real sigma_eff = 0;
		Real effort    = 0;
		if (sigma_1 < 0 and sigma_3 < 0 and sigma_3 < -compressive_strength) {
			// case (a) # from GLadky and Kuna 2017
			sigma_eff = -sigma_3;
			effort    = sigma_eff / compressive_strength;
		} else if (sigma_1 > 0 and sigma_3 > 0 and sigma_1 > tension_strength) {
			// case (b)
			sigma_eff = sigma_1;
			effort    = sigma_eff / tension_strength;
		} else if (math::abs(sigma_tau) > math::abs(tension_strength)) {
			// case (c)
			sigma_eff = math::abs(sigma_tau);
			effort    = sigma_eff / math::abs(tension_strength);
		} else {
			effort = 0;
		}
		i.second.effort = effort;
		/*
		std::cout << " CPP: " << bId << ": sigma_1: " << sigma_1 << "; sigma_2: " << sigma_2 << "; sigma_3: " << sigma_3 << "; sigma_tau: " << sigma_tau
		          << ": effort " << effort << std::endl;
		*/
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


struct PredicateClumpBreakage {
public:
	virtual bool                           in(const Vector3r& pt, Real pad = 0.) const = 0;
	virtual std::tuple<Vector3r, Vector3r> aabb() const                                = 0;
	Vector3r                               dim() const
	{
		const auto mnmx = aabb();
		const auto mn   = std::get<0>(mnmx);
		const auto mx   = std::get<1>(mnmx);
		return (mx - mn).eval();
	}
	Vector3r center() const
	{
		const auto mnmx = aabb();
		const auto mn   = std::get<0>(mnmx);
		const auto mx   = std::get<1>(mnmx);
		return .5 * (mn + mx);
	}
	virtual std::string call() const  = 0;
	virtual ~PredicateClumpBreakage() = default;
};

class iSphere : public PredicateClumpBreakage {
public:
	Vector3r center = Vector3r::Zero();
	Real     radius = 0;
	iSphere(const Vector3r& _center, Real _radius)
	        : center(_center)
	        , radius(_radius) {};
	std::tuple<Vector3r, Vector3r> aabb() const override
	{
		return std::make_tuple(
		        Vector3r(center[0] - radius, center[1] - radius, center[2] - radius),
		        Vector3r(center[0] + radius, center[1] + radius, center[2] + radius));
	}

	bool in(const Vector3r& pt, Real pad = 0.) const override
	{
		if ((pt - center).norm() <= (radius - pad)) {
			return true;
		} else {
			return false;
		}
	}
	std::string call() const override { return "iSphere"; }
};


struct GeometriePredicateClumpBreakage {
	// add - true
	// sub - false
	Vector3r                                                              aabb_min = Vector3r::Zero();
	Vector3r                                                              aabb_max = Vector3r::Zero();
	std::vector<std::pair<bool, std::shared_ptr<PredicateClumpBreakage>>> predicates;

	void add(std::shared_ptr<PredicateClumpBreakage> p)
	{
		if (predicates.empty()) {
			aabb_min = std::get<0>(p->aabb());
			aabb_max = std::get<1>(p->aabb());
		} else {
			const auto aabb_min_new = std::get<0>(p->aabb());
			const auto aabb_max_new = std::get<1>(p->aabb());

			aabb_min = Vector3r(
			        std::min(aabb_min[0], aabb_min_new[0]), std::min(aabb_min[1], aabb_min_new[1]), std::min(aabb_min[2], aabb_min_new[2]));
			aabb_max = Vector3r(
			        std::max(aabb_max[0], aabb_max_new[0]), std::max(aabb_max[1], aabb_max_new[1]), std::max(aabb_max[2], aabb_max_new[2]));
		}
		predicates.push_back(std::make_pair(true, p));
	}

	void sub(std::shared_ptr<PredicateClumpBreakage> p)
	{
		if (predicates.empty()) {
			aabb_min = std::get<0>(p->aabb());
			aabb_max = std::get<1>(p->aabb());
		}

		predicates.push_back(std::make_pair(false, p));
	}

	void call() const
	{
		for (const auto& i : predicates) {
			std::cout << i.second->call() << std::endl;
		}
	}

	bool getPoint(const Vector3r& v) const
	{
		bool pointMark   = false;
		bool firstObject = true;
		for (const auto& i : predicates) {
			const auto resPoint = i.second->in(v, 0);

			if (firstObject) {
				pointMark = i.first && resPoint;
			} else {
				if (i.first == true && resPoint == true) {
					pointMark = true; // Add point
				} else if (i.first == false && resPoint == true) {
					pointMark = false; // Subtract point
				}
			}
			firstObject = false;
		}
		return pointMark;
	}

	void generateGeo() const
	{
		/*
		std::ofstream outfile("test.geo");

		Real        PointsNumber = 70;
		unsigned long i            = 0;

		for (Real x = aabb_min[0]; x <= aabb_max[0]; x += (aabb_max[0] - aabb_min[0]) / PointsNumber) {
			for (Real y = aabb_min[1]; y <= aabb_max[1]; y += (aabb_max[1] - aabb_min[1]) / PointsNumber) {
				for (Real z = aabb_min[2]; z <= aabb_max[2]; z += (aabb_max[2] - aabb_min[2]) / PointsNumber) {
					if (getPoint(Vector3r(x, y, z))) {
						outfile << "Point(" << i << ") = {" << x << ", " << y << ", " << z << "};" << std::endl;
					}
					i++;
				}
			}
		}
		outfile.close();
		*/
	}
	std::vector<shared_ptr<Body>> generateHexa(Real radius, Real gap, shared_ptr<Material> mat) const
	{
		std::vector<shared_ptr<Body>> ret;
		const auto                    a  = 2 * radius + gap;
		const auto                    hy = a * sqrt(3.) / 2.;
		const auto                    hz = a * sqrt(6.) / 3.;

		const auto dim = aabb_max - aabb_min;
		const auto ii  = math::ceil(dim[0] / a);
		const auto jj  = math::ceil(dim[1] / hy);
		const auto kk  = math::ceil(dim[0] / hz);

		for (unsigned int i = 0; i < ii; i++) {
			for (unsigned int j = 0; i < jj; j++) {
				for (unsigned int k = 0; i < kk; k++) {
					const auto coordSph
					        = Vector3r((2 * i + ((j + k) % 2)), (sqrt(3.) * (j + 1. / 3. * (k % 2))), (2. * sqrt(6.) / 3. * k)) * (a / 2.0)
					        + aabb_min;

					if (getPoint(coordSph)) { ret.push_back(Shop::sphere(coordSph, radius, mat)); }
				}
			}
		}
		return ret;
	}
};


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


bool ClumpBreakage::replaceSphere(
        MapBody2StressTuple& bStress,
        Body::id_t           sphere_id,
        Real                 subparticles_mass,
        Real                 radius_ratio,
        Real                 relative_gap,
        Real                 initial_packing_scale,
        Real                 max_scale,
        bool                 search_for_neighbours,
        Real                 grow_radius,
        Real                 max_grow_radius,
        Vector3r             shift_vector)
{
	bool res = true;
	if (grow_radius > max_grow_radius) {
		res = false;
		std::cout << "grow_radius > max_grow_radius, not possible to handle" << std::endl;
	} else {
		const auto bProp = bStress[sphere_id];
		if (bProp.effort >= 1) {
			std::cout << "1a: We will crash the particle with id: " << sphere_id << std::endl;
			auto&      b                  = Body::byId(sphere_id, scene);
			const auto sphere             = dynamic_cast<Sphere*>(b->shape.get());
			const auto sphere_radius      = sphere->radius;
			const auto sphere_center      = b->state->pos;
			const auto full_mass          = b->state->mass;
			auto       subparticles_mass_ = subparticles_mass;

			if (subparticles_mass_ <= 0) { subparticles_mass_ = full_mass; }

			const auto full_number_of_subparticles = math::pow(radius_ratio, 3);
			const auto req_number_of_spheres = static_cast<unsigned int>(math::ceil(full_number_of_subparticles * subparticles_mass / full_mass));
			const auto subparticle_mass      = subparticles_mass / req_number_of_spheres;
			const auto mat_label             = b->material->label;
			const auto density               = b->material->density;
			const auto exact_radius          = math::pow(((3 * subparticle_mass) / (4 * M_1_PI * density)), (1 / 3));
			//

			std::set<Body::id_t> neighbour_set;
			if (search_for_neighbours) {
				for (const auto& i : b->intrs) {
					const auto I       = std::get<1>(i);
					Body::id_t otherId = -1;
					if (b->id == I->getId1()) {
						otherId = I->getId2();
					} else {
						otherId = I->getId1();
					}
					const auto& otherBody = Body::byId(otherId, scene);
					std::cout << "2a: otherId: " << otherId << std::endl;

					if (otherBody && otherBody->material && otherBody->shape && (otherBody->shape->getClassName() == "Sphere")) {
						const auto other_sph = dynamic_cast<Sphere*>(otherBody->shape.get());
						const auto other_rad = other_sph->radius;
						const auto dist      = ((otherBody->state->pos - b->state->pos).norm()) - other_rad * 2;
						if (dist < 0) { neighbour_set.insert(otherId); }
					}
				}
			}
			unsigned int                  no_gen_subparticles = 0;
			std::vector<shared_ptr<Body>> sp;
			while (no_gen_subparticles < req_number_of_spheres) { // #make sure enought subparticles were create)
				yade::GeometriePredicateClumpBreakage master_predicate;
				master_predicate.add(std::make_shared<iSphere>(sphere_center, sphere_radius * initial_packing_scale));
				for (const auto i : neighbour_set) {
					const auto& body1 = Body::byId(i, scene);
					const auto  sph1  = dynamic_cast<Sphere*>(body1->shape.get());
					const auto  rad1  = sph1->radius;
					const auto  pos1  = body1->state->pos;
					master_predicate.sub(std::make_shared<iSphere>(pos1, rad1));
					sp = master_predicate.generateHexa(exact_radius / grow_radius, relative_gap * exact_radius, b->material);
					no_gen_subparticles = sp.size();
					if (no_gen_subparticles < req_number_of_spheres) {
						initial_packing_scale *= 1.05;
						if (initial_packing_scale > max_scale) {
							// if sphere packing cant be found, at least decrease the size of released sphere to conserve the mass
							const auto new_particle_ratio = math::pow((subparticles_mass / full_mass), (1. / 3.));
							Shop::growParticle(b->id, grow_radius, true);
							const auto shift_distance = sphere_radius * (1 - new_particle_ratio);
							if (shift_vector != Vector3r::Zero()) { b->state->pos += shift_distance * shift_vector; }
						}
					}
				}
			}
			// make sure that the number of the particles is exact, remove the furthest subparticles
			std::vector<std::pair<Real, shared_ptr<Body>>> sq_distances;
			for (const auto& i : sp) {
				const auto distance = (i->state->pos - sphere_center).norm();
				sq_distances.push_back(std::make_pair(distance, i));
			}
			std::sort(sq_distances.begin(), sq_distances.end(), [](auto& left, auto& right) { return left.first > right.first; });
			unsigned int cnt = 0;
			for (const auto& i : sq_distances) {
				if (cnt == req_number_of_spheres) break;
				const auto b_id = scene->bodies->insert(i.second);
				Shop::growParticle(b_id, grow_radius, true);
				cnt++;
			}

		} else {
			res = false;
			std::cout << "We will NOT crash the particle with id" << sphere_id << std::endl;
		}
	}
	return res;
}

void ClumpBreakage::action()
{
	auto bStress = getStressForEachBody();
	checkFailure(bStress);
	/*
	for (auto& i : bStress) {
		const auto bId    = i.first;
		const auto effort = i.second.effort;
	}
	*/
	replaceSphere(bStress, 0);
}

} // namespace yade
