/*************************************************************************
*  Copyright (C) 2013 by Bruno Chareyre    <bruno.chareyre@grenoble-inp.fr>  *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once

#include <core/GlobalEngine.hpp>
#include <boost/tuple/tuple.hpp>
#include <set>

#include <core/Dispatching.hpp>
#include <pkg/dem/CapillaryPhys1.hpp>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace yade { // Cannot have #include directive inside.

class Interaction;

///This container class is used to check if meniscii overlap. Wet interactions are put in a series of lists, with one list per body.
class BodiesMenisciiList1 {
private:
	std::vector<std::list<shared_ptr<Interaction>>> interactionsOnBody;

	//shared_ptr<Interaction> empty;

public:
	BodiesMenisciiList1();
	BodiesMenisciiList1(Scene* body);
	bool                           prepare(Scene* scene);
	bool                           insert(const shared_ptr<Interaction>& interaction);
	bool                           remove(const shared_ptr<Interaction>& interaction);
	std::list<shared_ptr<Interaction>>& operator[](int index);
	int                            size();
	void                           display();
	bool initialized;
};

/// This is the constitutive law
class Law2_ScGeom_CapillaryPhys_Capillarity1 : public GlobalEngine {
public:
	void checkFusion();

	static DelaunayInterpolator::Dt                         dtVbased;
	static DelaunayInterpolator::Dt                         dtPbased;
	std::vector<MeniscusPhysicalData> solutions;
	int                               switchTriangulation; //to detect switches between P-based and V-based data


	BodiesMenisciiList1 bodiesMenisciiList;

	void action() override;
	Real intEnergy();
	Real swInterface();
	Real wnInterface();
	Real waterVolume();
	void solver(Real suction, bool reset);
	void triangulateData();
	shared_ptr<CapillaryPhys1> solveStandalone(Real R1, Real R2, Real pressure, Real gap, shared_ptr<CapillaryPhys1> bridge/*, bool pBased*/);

	// clang-format off
    YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_CapillaryPhys_Capillarity1,GlobalEngine,"This engine loops over interactions with physics :yref:`CapillaryPhys1` and it assign pendular bridges to them. It is a reimplementation of [Scholtes2009b]_, adding the option of imposing the bridge volume (instead of only capillary pressure) and enabling using unstructured input data by triangulation. This reimplementation also provides more geometrical quantities in order to define interfacial energy terms, it was used specifically in [Chalak2017].\n\nIf :yref:`Law2_ScGeom_CapillaryPhys_Capillarity1.imposePressure`==True, a uniform capillary pressure is assigned to all bridges, the liquid volume of each bridge is a result and it will change if the distance between the spheres changes. If :yref:`Law2_ScGeom_CapillaryPhys_Capillarity1.imposePressure`==False, then the volume of each bridge remains constant during motion, and capillary pressure is updated, instead. \n\nFor references, see [Scholtes2009b]_ and a couple papers by the same authors; [Scholtes2009d]_ (in french) is the most detailed.\n\nThe engine needs `an input data file <https://gitlab.com/yade-dev/yade-data/-/raw/887bfc12b6a8ad91024662fcf83efaf1dd01d968/capillaryFiles/capillaryfile.txt?inline=false>`__  available in yade-data package."
    "\n\nIn order to allow capillary forces between distant spheres, it is necessary to enlarge the bounding boxes using :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor` and make the Ig2 define define distant interactions via:yref:`interactionDetectionFactor<Ig2_Sphere_Sphere_ScGeom::interactionDetectionFactor>`. It is also necessary to disable interactions removal by the constitutive law (:yref:`Law2<Law2_ScGeom_FrictPhys_CundallStrack::neverErase>=True`). The only combinations of laws supported are currently capillary law + :yref:`Law2_ScGeom_FrictPhys_CundallStrack` and capillary law + :yref:`Law2_ScGeom_MindlinPhys_Mindlin` (and the other variants of Hertz-Mindlin).\n\nSee CapillaryPhys-example.py for an example script.",
		  ((Real,capillaryPressure,0.,,"Value of the capillary pressure Uc defines as Uc=Ugas-Uliquid"))
		  ((Real,totalVolumeofWater,-1.,,"Value of imposed water volume"))
		  ((Real,liquidTension,0.073,,"Value of the superficial water tension in N/m"))
		  ((Real,epsilonMean,0.,," Mean Value of the roughness"))//old value was epsilon
		  ((Real,disp,0.,," Dispersion from the mean Value of the roughness"))//added now
		  ((Real,interactionDetectionFactor,1.5,,"defines critical distance for deleting interactions. Must be consistent with the Ig2 value."))
                   ((bool,fusionDetection,false,,"If true potential menisci overlaps are checked"))
                   ((bool,initialized,false,," "))
                   ((bool,binaryFusion,true,,"If true, capillary forces are set to zero as soon as, at least, 1 overlap (menisci fusion) is detected"))
                   ((bool,hertzOn,false,,"|yupdate| true if hertz model is used"))
                   ((string,inputFilename,string("capillaryfile.txt"),,"the file with meniscus solutions, used for interpolation."))
                   ((bool,createDistantMeniscii,false,,"Generate meniscii between distant spheres? Else only maintain the existing one. For modeling a wetting path this flag should always be false. For a drying path it should be true for one step (initialization) then false, as in the logic of [Scholtes2009c]_"))
                   ((bool,imposePressure,true,," If True, suction is imposed and is constant if not Volume is imposed-Undrained test"))   
		  ((bool,totalVolumeConstant,true,," in undrained test there are 2 options, If True, the total volume of water is imposed,if false the volume of each meniscus is kept constant: in this case capillary pressure can be imposed for initial distribution of meniscus or it is the total volume that can be imposed initially")) 
		   ,switchTriangulation=-1,
		  .def("intEnergy",&Law2_ScGeom_CapillaryPhys_Capillarity1::intEnergy,"define the energy of interfaces in unsaturated pendular state")
// 		  .def("sninterface",&Law2_ScGeom_CapillaryPhys_Capillarity1::sninterface,"define the amount of solid-non-wetting interfaces in unsaturated pendular state")
		  .def("swInterface",&Law2_ScGeom_CapillaryPhys_Capillarity1::swInterface,"define the amount of solid-wetting interfaces in unsaturated pendular state")
		  .def("wnInterface",&Law2_ScGeom_CapillaryPhys_Capillarity1::wnInterface,"define the amount of wetting-non-wetiing interfaces in unsaturated pendular state")
		  .def("waterVolume",&Law2_ScGeom_CapillaryPhys_Capillarity1::waterVolume,"return the total value of water in the sample")		  
		  .def("solveStandalone",&Law2_ScGeom_CapillaryPhys_Capillarity1::solveStandalone, (py::args("R1"), py::args("R2"), py::args("pressure") , py::args("gap"),  py::args("bridge") = shared_ptr<CapillaryPhys1>()/*, py::args("pBased") = true*/), "Returns a :yref:`CapillaryPhys1` object representing a single bridge independently of the scene, using radii R1 and R2, capillary pressure, and gap between two spheres. The returned value contains internals of the interpolation process, it can be passed as an optional argument ('bridge'). If the resolution is repeated multiple times, re-using cached data will increase performance if the geometrical parameters are changing by small increments")
     );
	// clang-format on
};


REGISTER_SERIALIZABLE(Law2_ScGeom_CapillaryPhys_Capillarity1);

} // namespace yade
