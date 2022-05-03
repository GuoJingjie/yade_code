/*************************************************************************
*  Copyright (C) 2018 by Robert Caulk <rob.caulk@gmail.com>  		 *
*  Copyright (C) 2018 by Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>     *
*									 *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
/* This engine is under active development. Experimental only */
/* Thermal Engine was developed with funding provided by the Chateaubriand Fellowship */

/* Theoretical framework and experimental validation presented in:

Caulk, R. and Chareyre, B. (2019) An open framework for the simulation of thermal-hydraulic-mechanical processes in discrete element systems. Thermal Process Engineering: Proceedings of DEM8, Enschede Netherlands, July 2019.

*/

//#define THERMAL
#ifdef FLOW_ENGINE
#ifdef THERMAL
#ifdef YADE_OPENMP
#pragma once

#include <core/Dispatching.hpp>
#include <core/PartialEngine.hpp>
#include <core/State.hpp>
#include <pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include <pkg/dem/ScGeom.hpp>


//#include<pkg/pfv/FlowEngine.hpp>
#include "FlowEngine_FlowEngineT.hpp"
#include <pkg/pfv/FlowBoundingSphere.hpp>
#include <pkg/pfv/Network.hpp>
#include <lib/triangulation/Tesselation.h>
#include <pkg/dem/TesselationWrapper.hpp>

namespace yade { // Cannot have #include directive inside.

class ThermalEngine : public PartialEngine {
public:
	typedef TemplateFlowEngine_FlowEngineT<FlowCellInfo_FlowEngineT, FlowVertexInfo_FlowEngineT> FlowEngineT;
	typedef FlowEngineT::Tesselation                                                             Tesselation;
	typedef FlowEngineT::RTriangulation                                                          RTriangulation;
	typedef FlowEngineT::FiniteCellsIterator                                                     FiniteCellsIterator;
	typedef RTriangulation::Finite_facets_iterator                                               FiniteFacetsIterator;
	typedef FlowEngineT::CellHandle                                                              CellHandle;
	typedef FlowEngineT::VertexHandle                                                            VertexHandle;
	typedef FlowEngineT::Facet                                                                   Facet;
	typedef std::vector<CellHandle>                                                              VectorCell;
	typedef typename VectorCell::iterator                                                        VCellIterator;

public:
	Scene*       scene;
	bool         energySet; //initialize the internal energy of the particles
	FlowEngineT* flow;
	bool         timeStepEstimated;
	Real         thermalDT;
	bool         conductionIter;
	int          elapsedIters;
	int          conductionIterPeriod;
	Real         elapsedTime;
	bool         first;
	bool         runConduction;
	Real         maxTimeStep;
	Real         Pr;
	Real         Nu;
	Real         NutimesFluidK;
	Real         cavitySolidVolumeChange;
	Real         cavityVolume;
	Real         cavityDtemp;

	virtual ~ThermalEngine();
	void    action() override;
	void    setReynoldsNumbers();
	void    setInitialValues();
	void    applyTempDeltaToSolids(Real delT);
	void    resetFlowBoundaryTemps();
	void    resetBoundaryFluxSums();
	void    setConductionBoundary();
	void    thermalExpansion();
	void    initializeInternalEnergy();
	void    computeNewPoreTemperatures();
	void    computeNewParticleTemperatures();
	void    computeSolidFluidFluxes();
	void    computeFluidFluidConduction();
	void    updateForces();
	void    computeVertexSphericalArea();
	void    computeFlux(CellHandle& cell, const shared_ptr<Body>& b, const Real surfaceArea);
	void    computeSolidSolidFluxes();
	void    timeStepEstimate();
	CVector cellBarycenter(const CellHandle& cell);
	void    computeCellVolumeChangeFromDeltaTemp(CellHandle& cell, Real cavDens);
	void    accountForCavitySolidVolumeChange();
	void    accountForCavityThermalVolumeChange();
	void    unboundCavityParticles();
	void    computeCellVolumeChangeFromSolidVolumeChange(CellHandle& cell);
	Real    getThermalDT() const { return thermalDT; }
	int     getConductionIterPeriod() const { return conductionIterPeriod; }
	Real    getMaxTimeStep() const { return maxTimeStep; }
    void    makeThermal(); 
    bool    checkThermal(); // return false and print warning if some bodies don't have a thermal state
	//void         applyBoundaryHeatFluxes();
	// clang-format off
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalEngine,PartialEngine,"An engine typically used in combination with FlowEngine to simulate thermal-hydraulic-mechanical processes. Framework description and demonstration presented within the following paper [Caulk2019a]_ :Caulk, R.A. and Chareyre, B. (2019) An open framework for the simulation of thermal-hydraulic-mechanical processes in discrete element systems. Thermal Process Engineering: Proceedings of DEM8 International Conference for Discrete Element Methods, Enschede Netherlands, July 2019.",
		/*attributes*/
		((bool,advection,true,,"Activates advection"))
		((bool,fluidConduction,true,,"Activates conduction within fluid"))
		((bool,debug,false,,"debugging flags"))
		((bool,conduction,true,,"Activates conduction"))
		((bool,thermoMech,true,,"Activates thermoMech"))
        	((bool,fluidThermoMech,true,,"Activates thermoMech"))
        	((bool,solidThermoMech,true,,"Activates thermoMech"))
		((bool,ignoreFictiousConduction,false,,"Allows user to ignore conduction between fictious cells and particles. Mainly for debugging purposes."))
		((vector<bool>, bndCondIsTemperature, vector<bool>(6,false),,"defines the type of boundary condition for each side of particle packing. True if temperature is imposed, False for no heat-flux. Indices can be retrieved with :yref:`FlowEngine::xmin` and friends."))
		((vector<Real>, thermalBndCondValue, vector<Real>(6,0),,"Imposed temperature boundary condition for the particles."))
		((vector<Real>, thermalBndFlux, vector<Real>(6,0),,"Flux through thermal boundary."))
		((bool,boundarySet,false,,"set false to change boundary conditions"))
		((bool,useKernMethod,false,,"flag to use Kern method for thermal conductivity area calc"))
		((bool,useHertzMethod,false,,"flag to use hertzmethod for thermal conductivity area calc"))
        	((Real,fluidBeta,0.0002,,"volumetric temperature coefficient m^3/m^3C, default water, <= 0 deactivates"))
        	((Real,particleT0,0,,"Initial temperature of particles"))
		//((bool,useVolumeChange,false,,"Use volume change for thermal-mechanical-hydraulic coupling instead of pressure change. False by default."))
		((bool,letThermalRunFlowForceUpdates,false,,"If true, Thermal will run force updates according to new pressures instead of FlowEngine. only useful if useVolumeChange=false."))
		((bool,flowTempBoundarySet,true,,"set false to change boundary conditions"))
		((bool,unboundCavityBodies,true,,"automatically unbound bodies touching only cavity cells."))
		((Real,particleK,3.0,,"Particle thermal conductivity (W/(mK)"))
		((Real,particleCp,750.,,"Particle thermal heat capacity (J/(kgK)"))
		((Real,fluidConductionAreaFactor,1.,,"Factor for the porethroat area (used for fluid-fluid conduction model)"))
		((Real,particleAlpha,11.6e-6,,"Particle volumetric thermal expansion coeffcient"))
		((Real,particleDensity,0,,"If > 0, this value will override material density for thermodynamic calculations (useful for quasi-static simulations involving unphysical particle densities)"))
        ((Real,fluidK,0.580,,"Thermal conductivity of the fluid."))
		((Real,uniformReynolds,-1.,,"Control reynolds number in all cells (mostly debugging purposes). "))
		((Real,fluidBulkModulus,0,,"If > 0, thermalEngine uses this value instead of flow.fluidBulkModulus."))
		((Real, delT, 0,,"Allows user to apply a delT to solids and observe macro thermal expansion. Resets to 0 after one conduction step."))
        	((Real,tsSafetyFactor,0.8,,"Allow user to control the timstep estimate with a safety factor. Default 0.8. If <= 0, thermal timestep is equal to DEM"))
        	((Real,porosityFactor,0,,"If >0, factors the fluid thermal expansion. Useful for simulating low porosity matrices."))
        	((bool,tempDependentFluidBeta,false,,"If true, fluid volumetric thermal expansion coefficient, :yref:`ThermalEngine::fluidBeta`, is temperature dependent (linear model between 20-70 degC)"))
        	((Real,minimumFluidCondDist,0,,"Useful for maintaining stability despite poor external triangulations involving flat tetrahedrals. Consider setting to minimum particle diameter to keep scale."))
            ((unsigned,lenBodies,0,,"cache the number of thermal bodies to perform checks and raise warnings if newly inserted bodies are not thermal"))
		,
		/* extra initializers */
		,
		/* ctor */
		energySet=false;timeStepEstimated=false;thermalDT=0;elapsedTime=0;elapsedIters=0;conductionIterPeriod=1;first=true;runConduction=false;maxTimeStep=10000;Nu=0;NutimesFluidK=0;Pr=0;flow=NULL;
        makeThermal(); // automatically turn State's into ThermalState's, will not work if more spheres are instantiated after this engine
		,
		/* py */
        	.def("getThermalDT",&ThermalEngine::getThermalDT,"let user check estimated thermalDT .")
        	.def("getConductionIterPeriod",&ThermalEngine::getConductionIterPeriod,"let user check estimated conductionIterPeriod .")
        	.def("getMaxTimeStep",&ThermalEngine::getMaxTimeStep,"let user check estimated maxTimeStep.")
            .def("setReynoldsNumbers",&ThermalEngine::setReynoldsNumbers,"update the cell reynolds numbers manually (computationally expensive)")
            .def("makeThermal",&ThermalEngine::makeThermal,"Assign thermal states to all bodies.")
            .def("checkThermal",&ThermalEngine::checkThermal,"Check if all bodies have thermal states.")
	)
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ThermalEngine);

} // namespace yade

#endif //THERMAL
#endif //YADE_OPENMP
#endif //FLOW_ENGINE
