from yade import pack
import numpy as np
# required when functions are put in an external module
from yade import Sphere, math
from yade.utils import sphere, growParticle

"""
Functions required for clump breakage algorithm described in:
    Brzeziński, K., & Gladkyy, A. (2022), Clump breakage algorithm for DEM simulation of crushable aggregates.
Strength Criterion adopted from;
	Gladkyy, A., & Kuna, M. (2017). DEM simulation of polyhedral particle cracking using a combined Mohr–Coulomb–Weibull failure criterion. Granular Matter, 19(3), 41.
"""


def stressTensor(b, stress_correction=True):
    """
    Modification of Love-Weber stress tensor, that applied to the clump members gives results similar to standalone bodies.
    """
    center = b.state.pos
    radius = b.shape.radius
    volume = (np.pi*radius**3)*4/3
    sigma = np.zeros([3, 3]).astype(np.float32)
    interactions = b.intrs()

    if b.isStandalone:  # if not clumped
        f_u = O.forces.f(b.id)  # sum off forces for stress correction
        t_u = O.forces.t(b.id)  # sum off forces for stress correction
        for ii in interactions:
            pp = ii.geom.contactPoint
            contact_vector = pp-center
            cont_normal = ii.geom.normal
            ss = ii.phys.shearForce * \
                np.dot(contact_vector/np.linalg.norm(contact_vector), -cont_normal)
            nn = ii.phys.normalForce * \
                np.dot(contact_vector/np.linalg.norm(contact_vector), -cont_normal)
            ff = ss+nn
            for i in range(3):
                for j in range(3):
                    sigma[i, j] += ff[i]*contact_vector[j]
    else:
        c_id = b.clumpId
        # sum off forces for stress correction
        f_u = np.zeros(3).astype(np.float32)
        # sum of torques for stress correction
        t_u = np.zeros(3).astype(np.float32)
        for ii in interactions:
            # id of the other body (the one that b is contacting with)
            the_other_id = list({ii.id1, ii.id2}.difference({b.id}))[0]
            # ommit interaction if the contacting body belongs to the same clump
            if O.bodies[the_other_id].clumpId == c_id:
                continue
            pp = ii.geom.contactPoint
            contact_vector = pp-center
            cont_normal = ii.geom.normal
            ss = ii.phys.shearForce * \
                np.dot(contact_vector/np.linalg.norm(contact_vector), -cont_normal)
            nn = ii.phys.normalForce * \
                np.dot(contact_vector/np.linalg.norm(contact_vector), -cont_normal)
            ff = ss+nn
            f_u += ff
            t_u += np.cross(contact_vector, ss)
            for i in range(3):
                for j in range(3):
                    sigma[i, j] += ff[i]*contact_vector[j]
    # at least unbalanced force [Brzeziński and Gladkyy 2022] magnitude need to be greater than need to be greater than zero
    if stress_correction and np.linalg.norm(f_u) > 0:
        f_u_prim = -1*f_u
        normal_prim = f_u/np.linalg.norm(f_u)
        f_t_prim = np.cross(0.5*t_u, normal_prim)/radius
        f_t_bis = np.cross(0.5*t_u, -1*normal_prim)/radius
        false_interactions = [[radius*normal_prim, -1*radius*normal_prim],  # contact vectors in p_prim and p_bis respectively
                              [f_u_prim+f_t_prim, f_t_bis]]  # contact forces in p_prim and p_bis respectively
        for ii in range(len(false_interactions)):
            contact_vector = false_interactions[0][ii]
            ff = false_interactions[1][ii]
            for i in range(3):
                for j in range(3):
                    sigma[i, j] += ff[i]*contact_vector[j]

    sigma = sigma/volume
    return sigma


def checkFailure(b, tension_strength, compressive_strength, wei_V0=0.01, wei_P=0.63, wei_m=3, weibull=True,  stress_correction=True):
    """
    Strength criterion adopted from [Gladkyy and Kuna 2017]. Returns particles 'effort' (equivalent stress / strength) if the strength  is exceeded, and zero othervise.
    """
    sigma = stressTensor(b, stress_correction=stress_correction)

    if weibull:  # I just modify potential sigma_eff before checking conditions
        V = (np.pi*b.shape.radius**3)*4/3
        coeff = ((-wei_V0/V)*np.log(1-wei_P))**(1/wei_m)
        sigma /= coeff  # In gladky's paper there ia multiplication but sigma_eff_hat denotes strength, so I divide the acting stress

    # The next line does not work with the typ mpf! That is why we are skipping the check tests with mpf
    sigma_1, sigma_2, sigma_3 = np.sort(math.eig(sigma)[0])[::-1]
    sigma_tau = sigma_1 - sigma_3*tension_strength/compressive_strength

    # case (a) # from [Gladkyy and Kuna 2017]
    if sigma_1 < 0 and sigma_3 < 0 and sigma_3 < - compressive_strength:
        sigma_eff = -sigma_3
        effort = sigma_eff/compressive_strength
    # case (b)
    elif sigma_1 > 0 and sigma_3 > 0 and sigma_1 > tension_strength:
        sigma_eff = sigma_1
        effort = sigma_eff/tension_strength
    # case (c)
    elif abs(sigma_tau) > abs(tension_strength):
        sigma_eff = abs(sigma_tau)
        effort = sigma_eff/abs(tension_strength)
    else:
        effort = 0
    return effort


def replaceSphere(sphere_id, subparticles_mass=None, radius_ratio=2, relative_gap=0, neighbours_ids=[], initial_packing_scale=1.5, max_scale=3, scale_multiplier=None, search_for_neighbours=True, outer_predicate=None, grow_radius=1., max_grow_radius=2.0):
    """
    This function is intended to replace sphere with subparticles.
    It is dedicated for spheres replaced from clumps (but not only). Thus, two features are utilized:
    - subparticles_mass (mass of the subparticles after replacement), since in a clump only a fraction of original spheres mass is taken into account
    - neighbours_ids - list of ids of the neighbour bodies (e.g. other clump members, other bodies that sphere is contacting with) that we do not want to penetrate with new spheres (maybe it could be later use to avoid penetration of other bodies). However, passing neighbours_ids is not always necessary. By default (if search_for_neighbours==True), existing spheres are detected authomatically and added to neighbours_ids. Also, outer_predicate can beused to avoid penetrating other bodies with subparticles.
    Spheres will be initially populated in a hexagonal packing (predicate with dimension of sphere diameter multiplied by initial_packing_scale ). Initial packing scale is greater than one to make sure that sufficient number of spheres will be produced to obtain reguired mass (taking into account porosity).
    scale_multiplier - if a sufficient number of particles cannot be produced with initial packing scale, it is multiplied by scale multiplier. The procedure is repeated until initial_packing scale is reached. If scale_multiplier is None it will be changed to max_scale/initial_packing_scale, so the maximum range will be achieved in second iteration.
    max_scale - limits the initial_packing_scale which can be increased by the algorithm. If initial_packing_scale >max_scale, sphere will not be replaced (broken).
    outer_predicate - it is an additional constraint for subparticles generation. Can be used when non spherical bodies are in vicinity of the broken particle, particles are in box etc.
    search_for_neighbours - if True searches for additional neighbours (spheres whithin a range of initial_packing_scale *sphere_radius)

    Particles can be generated with smaller radius and than sligtly growed (by "grow_radius"). It allows for adding extra potential energy in the simulation, and increase the chances for successful packing.
    relative_gap - is the gap between packed subparticles (relalive to the radius of subparticle), note that if grow_radius > 1, during subparticles arangemenr their radius is temporarily decreased by 1/grow radius. It can be used to create special cases for overlapping (described in the paper).

    """
    assert grow_radius <= max_grow_radius
    sphere = O.bodies[sphere_id]
    sphere_radius = sphere.shape.radius
    sphere_center = sphere.state.pos
    full_mass = sphere.state.mass
    if subparticles_mass == None:
        subparticles_mass = full_mass
    if subparticles_mass == 0:
        print("Sphere (id={:d}) was erased without replacing with subparticles. The change of the mass is negligible (the updated mass of the clump has not changed). For grater accuracy increase the 'discretization' parameter.".format(sphere_id))
        O.bodies.erase(sphere_id)
        return None
    if scale_multiplier == None:
        scale_multiplier = max_scale/initial_packing_scale
    full_number_of_subparticles = radius_ratio**(3)
    req_number_of_spheres = int(
        np.ceil(full_number_of_subparticles*subparticles_mass/full_mass))
    subparticle_mass = subparticles_mass/req_number_of_spheres
    sphere_mat = sphere.material
    density = sphere_mat.density
    exact_radius = ((3*subparticle_mass)/(4*np.pi*density))**(1/3)
    # search for neighbours only once before the while loop
    # add to the list neighbours (spheres) whithin the range of max_scale (maximum_initial_packing_scale radius)
    if search_for_neighbours:
        neighbour_set = set(neighbours_ids)
        for bb in O.bodies:
            if bb.id != sphere_id and isinstance(bb.shape, Sphere):
                dist = ((np.array(sphere_center) -
                         np.array(bb.state.pos))**2).sum()**0.5
                if dist < sphere_radius*max_scale+bb.shape.radius:
                    neighbour_set.add(bb.id)
        neighbours_ids = list(neighbour_set)
    # populate space with new particles
    no_gen_subparticles = 0
    while no_gen_subparticles < req_number_of_spheres:  # make sure enought subparticles were create
        # make sure that there is enough of the particles
        master_predicate = pack.inSphere(
            sphere_center, sphere_radius*initial_packing_scale)
        for neighbour_id in neighbours_ids:  # only works for spheres
            nb = O.bodies[neighbour_id]
            master_predicate = master_predicate - \
                pack.inSphere(nb.state.pos, nb.shape.radius)
        if outer_predicate != None:
            master_predicate = master_predicate & outer_predicate
        sp = pack.regularHexa(master_predicate, radius=exact_radius /
                              grow_radius, gap=relative_gap*exact_radius, material=sphere_mat)
        no_gen_subparticles = len(sp)
        if no_gen_subparticles < req_number_of_spheres:
            initial_packing_scale *= scale_multiplier
            if initial_packing_scale > max_scale:
                print("We will NOT crash the particle with id id={:d}. ).".format(
                    sphere.id))
                return sphere.id
    # make sure that the number of the particles is exact, remove the furthest subparticles
    sq_distances = []
    for i in range(len(sp)):
        sq_distance = (np.array(sp[i].state.pos-sphere_center)**2).sum()
        # create a list that will store bodies and squared distances form the final master_sphere.
        sq_distances += [[sp[i], sq_distance]]
    # sort by distance - ascending
    sq_dist_sorted = sorted(sq_distances, key=lambda x: x[1])
    for i in range(req_number_of_spheres):
        b_id = O.bodies.append(sq_dist_sorted[i][0])
        growParticle(b_id, grow_radius)
    O.bodies.erase(sphere_id)
    return None


def evalClump(clump_id, radius_ratio, tension_strength, compressive_strength, relative_gap=0, wei_V0=0.001, wei_P=0.63, wei_m=3, weibull=True,  stress_correction=True, initial_packing_scale=1.5, max_scale=3, search_for_neighbours=True, outer_predicate=None, discretization=20, grow_radius=1., max_grow_radius=2.0):
    """
    Iterates over clump members with "checkFailure" function.
    Replaces the broken clump member with subparticles.
    Split new clump if necessary.
    If clump is not broken returns False, if broken True.
    """
    members_ids = O.bodies[clump_id].shape.members.keys()
    clump_mass = O.bodies[clump_id].state.mass
    old_clump_pos = O.bodies[clump_id].state.pos
    max_effort = 0
    member_to_break = -1
    for member_id in members_ids:
        effort = checkFailure(b=O.bodies[member_id], tension_strength=tension_strength, compressive_strength=compressive_strength,
                              wei_V0=wei_V0, wei_P=wei_P, wei_m=wei_m, weibull=weibull,  stress_correction=stress_correction)
        if effort > max_effort:
            max_effort = effort
            if effort >= 1:
                member_to_break = member_id
    if max_effort < 1:
        return False

    # breaking algorithm for two-sphere clump
    if len(members_ids) == 2:  # if only two spheres in the clump
        for member_id in members_ids:
            if member_id != member_to_break:
                conserved_mass = O.bodies[member_id].state.mass
        subparticles_mass = clump_mass - conserved_mass
        unbroken_sphere_id = replaceSphere(member_to_break, subparticles_mass=subparticles_mass, radius_ratio=radius_ratio, relative_gap=relative_gap, initial_packing_scale=initial_packing_scale,
                                           max_scale=max_scale, search_for_neighbours=search_for_neighbours, outer_predicate=outer_predicate, max_grow_radius=max_grow_radius, grow_radius=grow_radius)
        if not(unbroken_sphere_id is None):  # Finish if sphere cannot be broken
            return False
        # If the other sphere was broken, erase clump (but leave the sphere that wasn't broken).
        O.bodies.erase(clump_id)

    # breaking algorithm for more than two-sphere clump
    if len(members_ids) > 2:
        O.bodies.releaseFromClump(
            member_to_break, clump_id, discretization=discretization)
        conserved_mass = O.bodies[clump_id].state.mass
        subparticles_mass = clump_mass - conserved_mass
        # The line below is the same as in case of len(members_ids) == 2, but I can't do it outside condition because I need to check subparticles soundness in the same block of code
        unbroken_sphere_id = replaceSphere(member_to_break, subparticles_mass=subparticles_mass, radius_ratio=radius_ratio, relative_gap=relative_gap, initial_packing_scale=initial_packing_scale,
                                           max_scale=max_scale,  search_for_neighbours=search_for_neighbours, outer_predicate=outer_predicate, max_grow_radius=max_grow_radius, grow_radius=grow_radius)
        # if sphere cannot be broken first add it back to the clump, and then finish
        if not(unbroken_sphere_id is None):
            O.bodies.addToClump([member_to_break], clump_id,
                                discretization=discretization)
            return False
        # check if all the subparticles in the clump are connected otherwise create separate clumps/spheres
        list_of_sets = []  # sets to store members overlaping with the member 1 in each loop
        new_member_ids = O.bodies[clump_id].shape.members.keys()
        for member1 in new_member_ids:
            list_of_sets.append({member1})
            for member2 in new_member_ids:
                r1 = O.bodies[member1].shape.radius
                r2 = O.bodies[member2].shape.radius
                p1 = np.array(O.bodies[member1].state.pos)
                p2 = np.array(O.bodies[member2].state.pos)
                dist = ((p1-p2)**2).sum()**0.5
                if dist < r1 + r2:
                    list_of_sets[-1].add(member2)

        # Now sets are prepared, so verify whether there are any separete clumps.
        # The sets represent all the bodies that every member contacts with.
        # Union the intersecting sets until only distinct sets of elements remain.
        union_count = 1
        while union_count > 0:
            union_count = 0
            for i in range(len(list_of_sets)-1):
                for j in range(i+1, len(list_of_sets)):
                    s1 = list_of_sets[i]
                    s2 = list_of_sets[j]

                    if len(s1.intersection(s2)) > 0:
                        list_of_sets[i] = s1.union(s2)
                        list_of_sets[j].clear()
                        union_count += 1
                    else:
                        continue
        # This way I should have empty not important sets, sets of single values to replace by spheres, sets of multiple values to replace with clumps.
        # If any set before last one will not be empty, the whole clump needs to be erased. That means members will be erased to constitute other clumps and/or standalone spheres.
        len_of_list = len(list_of_sets)
        erase_clump = False
        for i in range(len(list_of_sets)):
            if len(list_of_sets[i]) == 0:
                continue
            else:
                if i < len_of_list:
                    erase_clump = True
            if len(list_of_sets[i]) == 1:
                member_id = list(list_of_sets[i])[0]
                O.bodies.append(sphere(
                    O.bodies[member_id].state.pos, O.bodies[member_id].shape.radius, material=O.bodies[member_id].material))
                O.bodies.erase(member_id)
            if len(list_of_sets[i]) > 1:
                members_ids = list(list_of_sets[i])
                list_of_spheres = []
                for member_id in members_ids:
                    list_of_spheres += [sphere(O.bodies[member_id].state.pos,
                                               O.bodies[member_id].shape.radius, material=O.bodies[member_id].material)]
                    O.bodies.erase(member_id)
                O.bodies.appendClumped(
                    list_of_spheres, discretization=discretization)

        if erase_clump:
            O.bodies.erase(clump_id)
    return True
