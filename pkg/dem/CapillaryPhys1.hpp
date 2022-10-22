/*************************************************************************
*  Copyright (C) 2013 by Bruno Chareyre    <bruno.chareyre@grenoble-inp.fr>  *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include <core/Dispatching.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/dem/FrictPhys.hpp>
#include <lib/triangulation/DelaunayInterpolation.hpp>

namespace yade { // Cannot have #include directive inside.

class MeniscusPhysicalData {
public:
	Real R;
	Real volume;
	Real distance;
	Real surface;
	Real energy;
	Real force;
	Real succion;
	Real delta1;
	Real delta2;
	Real arcLength;
	bool   ending;
	//default ctor
	MeniscusPhysicalData()
	        : R(0)
	        , volume(0)
	        , distance(0)
	        , surface(0)
	        , energy(0)
	        , force(0)
	        , succion(0)
	        , delta1(0)
	        , delta2(0)
	        , arcLength(0)
	        , ending(1)
	{
	}
	//ctor with input values
	MeniscusPhysicalData(
	        const Real& r,
	        const Real& v,
	        const Real& d,
	        const Real& s,
	        const Real& e,
	        const Real& f,
	        const Real& p,
	        const Real& a1,
	        const Real& a2,
	        const Real& arc,
	        bool          end)
	        : R(r)
	        , volume(v)
	        , distance(d)
	        , surface(s)
	        , energy(e)
	        , force(f)
	        , succion(p)
	        , delta1(a1)
	        , delta2(a2)
	        , arcLength(arc)
	        , ending(end)
	{
	}

	//a minimal list of operators for the interpolation
	//these operators are requirements for the DataType template parameter of interpolate()
	//FIXME: algebra operations must include energy, perimeter, and any other new variable for including them in the interpolation
	MeniscusPhysicalData& operator+=(const MeniscusPhysicalData& m2)
	{
		R += m2.R;
		volume += m2.volume;
		distance += m2.distance;
		surface += m2.surface;
		energy += m2.energy;
		force += m2.force;
		succion += m2.succion;
		delta1 += m2.delta1;
		delta2 += m2.delta2;
		arcLength += m2.arcLength;
		ending = (ending && m2.ending);
		return *this;
	}

	MeniscusPhysicalData operator*(const Real& fact) const
	{
		return MeniscusPhysicalData(
		        fact * R,
		        fact * volume,
		        fact * distance,
		        fact * surface,
		        fact * energy,
		        fact * force,
		        fact * succion,
		        fact * delta1,
		        fact * delta2,
		        fact * arcLength,
		        ending);
	}
};

//The structure for the meniscus data: physical properties + cached values for fast interpolation
typedef DelaunayInterpolator::InterpolatorData<MeniscusPhysicalData> Meniscus;

class CapillaryPhys1 : public FrictPhys {
public:
// 	int      currentIndexes[4]; // used for faster interpolation (stores previous positions in tables)
	Meniscus m;

	virtual ~CapillaryPhys1();

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CapillaryPhys1,FrictPhys,"Physics (of interaction) for Law2_ScGeom_CapillaryPhys_Capillarity1.",
				 ((bool,meniscus,false,,"Presence of a meniscus if true"))
				 ((bool,isBroken,false,,"If true, capillary force is zero and liquid bridge is inactive."))
				 ((bool,computeBridge,true,,"If true, capillary bridge will be computed if not it will be ignored."))
				 ((Real,capillaryPressure,0.,,"Value of the capillary pressure Uc defines as Ugas-Uliquid"))
				 ((Real,vMeniscus,0.,,"Volume of the menicus"))
				 ((Real,Delta1,0.,,"Defines the surface area wetted by the meniscus on the smallest grains of radius R1 (R1<R2)"))
				 ((Real,Delta2,0.,,"Defines the surface area wetted by the meniscus on the biggest grains of radius R2 (R1<R2)"))
				 ((Vector3r,fCap,Vector3r::Zero(),,"Capillary Force produces by the presence of the meniscus"))
				 ((Real,SInterface,0.,,"Fluid-Gaz Interfacial area"))
				 ((Real,arcLength,0.,,"Arc Length of the Fluid-Gaz Interface"))
				 ((short int,fusionNumber,0.,,"Indicates the number of meniscii that overlap with this one"))
				 ,createIndex();
				 );
	// clang-format on
	REGISTER_CLASS_INDEX(CapillaryPhys1, FrictPhys);
};
REGISTER_SERIALIZABLE(CapillaryPhys1);

class Ip2_FrictMat_FrictMat_CapillaryPhys1 : public IPhysFunctor {
public:
	virtual void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction) override;

	FUNCTOR2D(FrictMat, FrictMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_FrictMat_FrictMat_CapillaryPhys1,IPhysFunctor, "RelationShips to use with Law2_ScGeom_CapillaryPhys_Capillarity1\n\n In these RelationShips all the interaction attributes are computed. \n\n.. warning::\n\tas in the others :yref:`Ip2 functors<IPhysFunctor>`, most of the attributes are computed only once, when the interaction is new.",
				       ((bool,computeDefault,true,,"bool to assign the default value of computeBridge.")),;
	  
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_CapillaryPhys1);

} // namespace yade
