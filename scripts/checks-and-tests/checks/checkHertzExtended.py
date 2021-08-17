#script with three contact pairs using Conical Damage Model as contact law:
#1st contact pair: sphere contacts facet
#2nd contact pair: sphere contacts wall
#3rd contact pair: sphere contacts sphere
#after 40 iterations, the following values are compared to precomputeted values: contact overlap, normal force, contact radius and stress
#------------------------------------------------------------------------------------------------------------
import numpy as np

if ('CGAL' in features):
	#-------------------------Define material
	young = 1e10
	poisson = 0.2
	sigmaMax = 6e7  #
	alpha = np.radians(85)  # in radians, 0<= alpha <= pi/2

	MatCDM = O.materials.append(
	        FrictMatCDM(young=young, density=3000, poisson=poisson, frictionAngle=np.arctan(0.1), sigmaMax=sigmaMax, alpha=alpha, c1=0.0, c2=0.0)
	)

	#-------------------------Define engines

	#geometric functors for Sphere-Wall Sphere-Facet contacts use Hertzian theory
	GeomFunctorList = [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(hertzian=True), Ig2_Wall_Sphere_ScGeom(hertzian=True)]

	O.engines = [
	        ForceResetter(),
	        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(),
	                               Bo1_Box_Aabb(), Bo1_Wall_Aabb(),
	                               Bo1_Polyhedra_Aabb()], allowBiggerThanPeriod=True),
	        InteractionLoop(
	                GeomFunctorList,
	                [Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM(),
	                 Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM(),
	                 Ip2_FrictMat_FrictMat_MindlinPhys()], [Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM(label="lawHMCDM")]
	        ),
	        NewtonIntegrator(gravity=(0, 0, 0.0), damping=0.0)
	]

	#------------------------- sphere radius and position
	rad_p = 0.1
	L = 1.0

	#------------------- Contact pair 1: Sphere contacts Facet box  ----------------------

	O.bodies.append([sphere((L / 5.0, L / 2.0, -rad_p), rad_p, material=MatCDM)])
	Sphere1 = O.bodies[-1]
	O.bodies.append(geom.facetBox(center=(L / 2.0, L / 2.0, L / 2.0), extents=(L / 2.0, L / 2.0, L / 2.0), material=MatCDM, fixed=True))
	FacetBox1 = O.bodies[10]

	#------------------- Contact pair 2: Sphere contacts Wall  ----------------------

	O.bodies.append([sphere((-L * 2, L / 2.0, -rad_p - L), rad_p, material=MatCDM)])
	Sphere2 = O.bodies[-1]
	O.bodies.append(utils.wall((-L * 2, L / 2.0, -L), axis=2, sense=0, material=MatCDM))
	Wall2 = O.bodies[-1]

	#------------------- Contact pair 3: Sphere contacts (fixed) Sphere  ----------------------

	O.bodies.append([sphere((2 * L, L / 2.0, -rad_p), rad_p, material=MatCDM)])
	Sphere3 = O.bodies[-1]
	O.bodies.append([sphere((2 * L, L / 2.0, rad_p), rad_p, material=MatCDM, fixed=True)])
	SphereFixed3 = O.bodies[-1]

	# ------------------------- apply velovity on spheres -------------------------
	O.bodies[Sphere1.id].state.vel = Vector3(0, 0, 0.05)
	O.bodies[Sphere2.id].state.vel = Vector3(0, 0, 0.05)
	O.bodies[Sphere3.id].state.vel = Vector3(0, 0, 0.05)

	#-------------------------  run 40 time steps
	O.dt = .8 * PWaveTimeStep()

	O.run(40, True)

	#-------------------------  check results
	if not (O.interactions[Sphere1.id, FacetBox1.id].geom):
		raise YadeCheckError("checkHertzExtended: sphere-facet contact not existing!")  #Test is failed
	if not (O.interactions[Sphere2.id, Wall2.id].geom):
		raise YadeCheckError("checkHertzExtended: sphere-wall contact not existing!")  #Test is failed
	if not (O.interactions[Sphere3.id, SphereFixed3.id].geom):
		raise YadeCheckError("checkHertzExtended: sphere-sphere contact not existing!")  #Test is failed

	Fz1 = np.abs(O.forces.f(Sphere1.id)[2])
	Da = O.interactions[Sphere1.id, FacetBox1.id].geom.refR1
	Db = O.interactions[Sphere1.id, FacetBox1.id].geom.refR2
	delta1 = O.interactions[Sphere1.id, FacetBox1.id].geom.penetrationDepth
	contactRadius1 = O.interactions[Sphere1.id, FacetBox1.id].phys.R
	Rorg = (Da * Db / (Da + Db))
	afac = (1 - np.sin(alpha)) / np.sin(alpha)
	delta_el = max((delta1 - (contactRadius1 - Rorg) * afac), 0)
	Es = young / 2.0 / (1 - poisson**2)
	stress1 = 2 * Es / np.pi * np.sqrt(delta_el / contactRadius1)

	Fz2 = np.abs(O.forces.f(Sphere2.id)[2])
	Da = O.interactions[Sphere2.id, Wall2.id].geom.refR1
	Db = O.interactions[Sphere2.id, Wall2.id].geom.refR2
	delta2 = O.interactions[Sphere2.id, Wall2.id].geom.penetrationDepth
	contactRadius2 = O.interactions[Sphere2.id, Wall2.id].phys.R
	Rorg = (Da * Db / (Da + Db))
	afac = (1 - np.sin(alpha)) / np.sin(alpha)
	delta_el = max((delta2 - (contactRadius2 - Rorg) * afac), 0)
	Es = young / 2.0 / (1 - poisson**2)
	stress2 = 2 * Es / np.pi * np.sqrt(delta_el / contactRadius2)
	Fz3 = np.abs(O.forces.f(Sphere3.id)[2])
	Da = O.interactions[Sphere3.id, SphereFixed3.id].geom.refR1
	Db = O.interactions[Sphere3.id, SphereFixed3.id].geom.refR2
	delta3 = O.interactions[Sphere3.id, SphereFixed3.id].geom.penetrationDepth
	contactRadius3 = O.interactions[Sphere3.id, SphereFixed3.id].phys.R
	Rorg = (Da * Db / (Da + Db))
	afac = (1 - np.sin(alpha)) / np.sin(alpha)
	delta_el = max((delta3 - (contactRadius3 - Rorg) * afac), 0)
	Es = young / 2.0 / (1 - poisson**2)
	stress3 = 2 * Es / np.pi * np.sqrt(delta_el / contactRadius3)

	tolerance = 1e-5

	if abs((delta1 - 5.61227e-05) / 5.61227e-05) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-facet contact: difference on overlap")  #Test is failed
	if abs((Fz1 - 459.1828) / 459.1828) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-facet contact: difference on normal force")  #Test is failed
	if abs((contactRadius1 - 0.1056368) / 0.1056368) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-facet contact: difference on contact radius")  #Test is failed
	if abs((stress1 - 60000000.0) / 60000000.0) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-facet contact: difference on stress")  #Test is failed

	if abs((delta2 - 5.61227e-05) / 5.61227e-05) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-wall contact: difference on overlap")  #Test is failed
	if abs((Fz2 - 459.1828) / 459.1828) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-wall contact: difference on normal force")  #Test is failed
	if abs((contactRadius2 - 0.1056368) / 0.1056368) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-wall contact: difference on contact radius")  #Test is failed
	if abs((stress2 - 60000000.0) / 60000000.0) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-wall contact: difference on stress")  #Test is failed

	if abs((delta3 - 7.43693e-05) / 7.43693e-05) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-sphere contact: difference on overlap")  #Test is failed
	if abs((Fz3 - 168.4617) / 168.4617) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-sphere contact: difference on normal force")  #Test is failed
	if abs((contactRadius3 - 0.0639842) / 0.0639842) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-sphere contact: difference on contact radius")  #Test is failed
	if abs((stress3 - 60000000.0) / 60000000.0) > tolerance:
		raise YadeCheckError("checkHertzExtended: sphere-sphere contact: difference on stress")  #Test is failed
else:
	print("This checkHertzExtended.py cannot be executed because CGAL is disabled")
