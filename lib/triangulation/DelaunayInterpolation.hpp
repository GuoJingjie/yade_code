/*************************************************************************
*  Copyright (C) 2013 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr> *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/*
This template class uses 3D triangulation to interpolate functions of three variables.
Exemple usage in CapillaryPhysDelaunay (see class "Meniscus") and CapillarityEngine.cpp (where the interpolation happens).
*/


#pragma once
#include <lib/base/AliasCGAL.hpp>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

namespace yade {

class DelaunayInterpolator {
    public:
    //helpful array for permutations
    static const unsigned  comb[]; // = { 1, 2, 3, 0, 1, 2 };
    
    // NOTE: make sure this static array is defined in some source file (currently CapillaryPhysDelaunay.cpp).
    // the folowing declararation/definition would be more convenient but not all gcc versions support this
    // static constexpr int  comb[6] = { 1, 2, 3, 0, 1, 2 };
    
    typedef CGAL::Exact_Real_predicates_inexact_constructions_kernel K;
    typedef CGAL::Delaunay_triangulation_3< K,
        CGAL::Triangulation_data_structure_3< CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>, CGAL::Triangulation_cell_base_3<K> > >
        Dt;
    typedef std::vector<std::pair<Dt::Vertex_handle, K::FT>>                    Vertex_weight_vector;
    typedef std::back_insert_iterator<Vertex_weight_vector>                     VtxWeightsIterator;
    typedef std::pair<VtxWeightsIterator , bool>                                    Weighted_vertices;   
    
    // helper class to convert an actual data class into an interpolator-friendly class which will add some internals of the interpolation
    // on the top of the actual solution, in order to speed-up things in repetead solve.
    // Exemple usage in CapillaryPhysDelaunay (see class "Meniscus") and CapillarityEngine.cpp (where the interpolation happens).
    template <class PhysicalData>
    class InterpolatorData : public PhysicalData {
        public:
        typedef PhysicalData Data;
        Data                         data;    // the interpolated solution 
        Dt::Cell_handle              cell;    // pointer to the last location in the triangulation, for faster locate()
        std::vector<K::Vector_3>     normals; // 4 normals relative to the current cell

        InterpolatorData()
                : data()
                , cell(Dt::Cell_handle())
                , normals(std::vector<K::Vector_3>())
        {
        }
    };

// The function returning vertices and their weights for an arbitrary point in R3 space.
// The returned triplet contains:
// 	- the output iterator pointing to the filled vector of vertices
// 	- the sum of weights for normalisation
// 	- a bool telling if the query point was located inside the convex hull of the data points
static Weighted_vertices getIncidentVtxWeights(
        const Dt&                            dt,
        const Dt::Geom_traits::Point_3&      Q,
        VtxWeightsIterator  nn_out,
        std::vector<Dt::Geom_traits::Vector_3>& normals,
        Dt::Cell_handle& start)
{

    typedef Dt::Geom_traits::FT FT;
	CGAL_triangulation_precondition(dt.dimension() == 3);
	Dt::Locate_type lt;
	int         li, lj;
	Dt::Cell_handle c             = dt.locate(Q, lt, li, lj, start);
	bool        updateNormals = (c != start || normals.size() < 4);
	if (updateNormals) normals.clear();
	if (lt == Dt::VERTEX) {
		*nn_out++ = std::make_pair(c->vertex(li), FT(1));
		return std::make_pair(nn_out, true);
	} else if (dt.is_infinite(c))
		return std::make_pair(nn_out, false); //point outside the convex-hull
	for (int k = 0; k < 4; k++) {
		if (updateNormals) {
			normals.push_back(cross_product(
			        c->vertex(comb[k])->point() - c->vertex(comb[k + 1])->point(), c->vertex(comb[k])->point() - c->vertex(comb[k + 2])->point()));
			normals[k] = normals[k] / ((c->vertex(k)->point() - c->vertex(comb[k])->point()) * normals[k]);
		}
        FT weight = ((Q - c->vertex(comb[k])->point()) * normals[k]);
		*nn_out++   = std::make_pair(c->vertex(k), weight);
	}
	start = c;
	return std::make_pair(nn_out, true);
}

template <class DataOwner>
static typename DataOwner::Data
interpolate1(const Dt& dt, const Dt::Geom_traits::Point_3& Q, DataOwner& owner, const std::vector<typename DataOwner::Data>& rawData, bool reset)
{
	Vertex_weight_vector weighted_vertices;
    weighted_vertices.reserve(4);
	if (reset) owner.cell = dt.cells_begin();
	auto result = getIncidentVtxWeights(dt, Q, std::back_inserter(weighted_vertices), owner.normals, owner.cell);

	typename DataOwner::Data data = typename DataOwner::Data(); //initialize null solution
	if (!result.second) return data;                             // out of the convex hull, we return the null solution
	//else, we compute the weighted sum
	for (unsigned int k = 0; k < weighted_vertices.size(); k++)
		data += (rawData[weighted_vertices[k].first->info()] * weighted_vertices[k].second);
	if (!data.ending) return data; //inside the convex hull but outside the concave hull of solutions
	else
		return typename DataOwner::Data();
}

};
    
} // namespace yade
