/*************************************************************************
*  Copyright (C) 2019 by Robert Caulk <rob.caulk@gmail.com>              *
*  Copyright (C) 2019 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//#ifndef FLOW_GUARD
//#define FLOW_GUARD
#pragma once
#ifdef FLOW_ENGINE
#ifdef LINSOLV
#ifdef PARTIALSAT
#include <core/Body.hpp>
#include <core/Omega.hpp>
#include <core/PartialEngine.hpp>
#include <core/Scene.hpp>
#include <core/State.hpp>
#include <core/Dispatching.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <Eigen/SparseLU>
#include <lib/triangulation/FlowBoundingSphere.hpp>
#include <lib/triangulation/Network.hpp>
#include <lib/triangulation/Tesselation.h>
#include <pkg/dem/TesselationWrapper.hpp>
#include "FlowEngine_PartialSatClayEngineT.hpp"
#include <cholmod.h>

namespace yade { // Cannot have #include directive inside.

class PartialSatCellInfo : public FlowCellInfo_PartialSatClayEngineT {
public:
	Real saturation; //the saturation of single pore (will be used in quasi-static imbibition and dynamic flow)
	Real porosity;
	//Real solidLine [4][4];//the length of intersecting line between sphere and facet. [i][j] is for facet "i" and sphere (facetVertices)"[i][j]". Last component [i][3] for 1/sumLines in the facet "i" (used by chao).
	Real              dsdp; // the change of saturation for given capillary pressure
	bool              crack;
	short int         crackNum;
	Real              crackArea;
	bool              isExposed;
	Real              initialPorosity;
	Real              initialSaturation;
	Real              Po;      // parameter for water retention curve, unique to cell
	Real              lambdao; // parameter for water retention curve
	Real              vSolids;
	Real              equivalentBulkModulus;
	Real              oldPressure;
	bool              clumped;
	bool              initiallyCracked; // flag to continue computing initial crack perm during simulation without needing displacements or broken bonds
	bool              entered;          // flag for a cracked cell that has become fully saturated
	std::vector<bool> opened;
	std::vector<bool> entry;
	std::vector<bool> visited;

	PartialSatCellInfo(void)
	{
		entered           = false;
		initialPorosity   = 0.25;
		initialSaturation = 0.01;
		oldPressure       = 0;
		crackArea         = 0;
		saturation        = 0.01;
		porosity          = 0.25;
		dsdp              = 0;
		vSolids           = 0;
		crack             = false;
		crackNum          = 0;
		isExposed         = false; // flag for determining if a pore is exposed to atmosphere, which controls the pressure force calculations
		Po                = 1.5;
		lambdao           = 0.2;
		clumped           = false;
		blocked           = false;
		initiallyCracked  = false;
		opened.resize(4, 0);
		entry.resize(4, 0);
		visited.resize(4, 0);
	}

	inline Real& sat(void) { return saturation; }
};

class PartialSatVertexInfo : public FlowVertexInfo_PartialSatClayEngineT {
public:
	//same here if needed
};

typedef CGT::_Tesselation<CGT::TriangulationTypes<PartialSatVertexInfo, PartialSatCellInfo> > PartialSatTesselation;
#ifdef LINSOLV
// #define PartialSatBoundingSphere CGT::PartialSatLinSolv<PartialSatTesselation>
class PartialSatBoundingSphere : public CGT::PartialSatLinSolv<PartialSatTesselation> {};
#endif

typedef TemplateFlowEngine_PartialSatClayEngineT<PartialSatCellInfo, PartialSatVertexInfo, PartialSatTesselation, PartialSatBoundingSphere>
        PartialSatClayEngineT;

REGISTER_SERIALIZABLE(PartialSatClayEngineT);

class PartialSatClayEngine : public PartialSatClayEngineT {
public:
	//typedef TemplateFlowEngine_FlowEngineT<FlowCellInfo_FlowEngineT,FlowVertexInfo_FlowEngineT> FlowEngineT;
	typedef PartialSatClayEngineT::Tesselation         Tesselation;
	typedef PartialSatClayEngineT::RTriangulation      RTriangulation;
	typedef PartialSatClayEngineT::FiniteCellsIterator FiniteCellsIterator;
	typedef PartialSatClayEngineT::CellHandle          CellHandle;
	typedef PartialSatClayEngineT::VertexHandle        VertexHandle;
	typedef std::vector<CellHandle>                    VectorCell;
	typedef typename VectorCell::iterator              VCellIterator;

public:
	//PartialSatClayEngineT* clayFlow;
	shared_ptr<FlowSolver> gas_solver;
	Real                   dsdp(CellHandle& cell);
	void                   initializeSaturations(FlowSolver& flow);
	void                   setSaturationFromPcS(CellHandle& cell);
	void                   setCellsDSDP(FlowSolver& flow);
	Real                   getTotalVolume();
	void                   updateSaturation(FlowSolver& flow);
	void                   updateBoundarySaturation(FlowSolver& flow);
	//	void triangulate(FlowSolver& flow);
	//Real diagonalSaturationContribution(CellHandle cell);
	//Real RHSSaturationContribution(CellHandle cell);
	void collectParticleSuction(FlowSolver& flow);
	void swellParticles();
	void setOriginalParticleValues();
	//void resetOriginalParticleValues() {this->setOriginalParticleValues();}
	void setInitialPorosity(FlowSolver& flow);
	void updatePorosity(FlowSolver& flow);
	void computeEquivalentBulkModuli(FlowSolver& flow);
	//	void updatePcSParams();
	Real laplaceDeviate(Real mu, Real b);
	Real exponentialDeviate(Real a, Real b);
	Real weibullDeviate(Real lambda, Real k);
	Real vanGenuchten(CellHandle& cell, Real pc);
	void exposureRecursion(CellHandle cell, Real bndPressure);
	void determineFracturePaths(FlowSolver& flow);
	void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius);
	void insertMicroPores(const Real fracMicroPore);
	bool findInscribedRadiusAndLocation(CellHandle& cell, std::vector<Real>& coordAndRad);
	bool checkSphereContainedInTet(CellHandle& cell, std::vector<Real>& coordAndRad);
	void setPorosityWithImageryGrid(string imageryFilePath, FlowSolver& flow);
	void resetPoresVolumeSolids(FlowSolver& flow);
	void getGasPermeability();
	void artificialParticleSwell(const Real volStrain);
	void setHomogeneousSuction(FlowSolver& flow);

	void       blockLowPoroRegions(FlowSolver& flow);
	void       blockMineralCellRecursion(CellHandle cell, std::vector<Body::id_t>& clumpIds);
	void       addIncidentParticleIdsToClumpList(CellHandle cell, std::vector<Body::id_t>& clumpIds);
	Body::id_t clump(vector<Body::id_t> ids, unsigned int discretization);
	void       printPorosityToFile(string file);
	void       addPermanentForces(FlowSolver& flow);

	// fracture network
	void crackCellsAbovePoroThreshold(FlowSolver& flow);
	void trickPermeability(FlowSolver& flow, bool gasPermFlag);
	//void trickPermOnCrackedCells(FlowSolver& flow);
	void interpolateCrack(Tesselation& Tes, Tesselation& NewTes);
	void computeFracturePerm(
	        RTriangulation::Facet_circulator&      facet,
	        Real                                   aperture,
	        RTriangulation::Finite_edges_iterator& edge,
	        const Real                             openingPressure,
	        bool                                   gasPermFlag,
	        FlowSolver&                            flow);
	void circulateFacets(RTriangulation::Finite_edges_iterator& edge, Real aperture, const Real openingPressure, bool gasPermFlag, FlowSolver& flow);
//void saveFractureNetworkVTK(string fileName);
#ifdef YADE_VTK
	void savePermeabilityNetworkVTK(string fileName);
	void saveUnsatVtk(const char* folder, bool withBoundaries);
#endif
	void blockCellsAbovePoroThreshold(FlowSolver& flow);
	void blockIsolatedCells(FlowSolver& flow);
	void simulateConfinement();
	void removeForceOnVertices(RTriangulation::Facet_circulator& facet, RTriangulation::Finite_edges_iterator& ed_it);
	void circulateFacetstoRemoveForces(RTriangulation::Finite_edges_iterator& edge);
	void removeForcesOnBrokenBonds();
	void timeStepControl();
	void computeVertexSphericalArea();
	void resetParticleSuctions();
	Real getEnteredRatio();

	//	void setPositionsBuffer(bool current);
	Real     leakOffRate                 = 0;
	Real     averageAperture             = 0;
	Real     sumOfApertures              = 0;
	Real     averageFracturePermeability = 0;
	Real     maxAperture                 = 0;
	Real     crackArea                   = 0;
	Real     crackVolume                 = 0;
	Real     waterVolume                 = 0;
	int      numCracks                   = 0;
	Real     totalFractureArea           = 0;
	Real     solverDT                    = 0;
	bool     emulatingAction             = false;
	Real     maxDSDPj                    = 0;
	Vector3r crack_fabric_vector;
	Real     crack_fabric_area;

	virtual void initializeVolumes(FlowSolver& flow);
	virtual void updateVolumes(FlowSolver& flow);
	virtual void buildTriangulation(Real pZero, Solver& flow, bool oneTes = false);
	virtual void initSolver(FlowSolver& flow);
	virtual void action();
	virtual void emulateAction()
	{
		scene           = Omega::instance().getScene().get();
		emulatingAction = true;
		action();
		emulatingAction = false;
	}

	void reloadSolver(FlowSolver& flow) { this->initSolver(flow); }

	virtual ~PartialSatClayEngine();
	Real getCrackArea() const { return crackArea; }
	Real getCrackVolume() const { return crackVolume; }
	Real getWaterVolume() const { return waterVolume; }
	Real getNumCracks() const { return numCracks; }
	Real getAverageAperture() const { return averageAperture; }
	void printPorosity(string file) { printPorosityToFile(file); }
	Real getCellSaturation(Vector3r pos) const { return solver->getCellSaturation(pos[0], pos[1], pos[2]); }
	//Real getCellVelocity(Vector3r pos){return solver->getCellVelocity(pos[0], pos[1], pos[2]);}
	std::vector<Real> getCellVelocity(Vector3r pos) const { return solver->getCellVelocity(pos[0], pos[1], pos[2]); }
	std::vector<Real> getCellGasVelocity(Vector3r pos) const { return gas_solver->getCellGasVelocity(pos[0], pos[1], pos[2]); }
	Real              getCellPorosity(Vector3r pos) const { return solver->getCellPorosity(pos[0], pos[1], pos[2]); }
	Real              getCellVolume(Vector3r pos) const { return solver->getCellVolume(pos[0], pos[1], pos[2]); }
	Real              getCellGasVolume(Vector3r pos) const { return gas_solver->getCellGasVolume(pos[0], pos[1], pos[2]); }
	Vector3r          getCellGasCenter(unsigned int id) const { return makeVector3r(gas_solver->T[gas_solver->currentTes].cellHandles[id]->info()); }
	bool              getCellGasPImposed(unsigned int id) const { return gas_solver->T[gas_solver->currentTes].cellHandles[id]->info().Pcondition; }
	bool              getCellCracked(Vector3r pos) const { return solver->getCellCracked(pos[0], pos[1], pos[2]); }
	Real              getAverageSaturation() const { return solver->getAverageSaturation(); }
	Real              getEnteredThroatRatio() { return getEnteredRatio(); }
	Real              getAverageSuction() const { return solver->getAverageSuction(); }
	Real              getTotalSpecimenVolume() { return getTotalVolume(); }
	Vector3r          getCrackFabricVector() const { return crack_fabric_vector / crack_fabric_area; }
	Real              getBoundaryGasFlux(unsigned int boundary) const { return gas_solver->boundaryFlux(boundary); }
	unsigned int      nGasCells() const { return gas_solver->T[gas_solver->currentTes].cellHandles.size(); }
	CELL_SCALAR_GETTER(Real, .sat(), cellSaturation);
	CELL_SCALAR_SETTER(Real, .sat(), setCellSaturation);

	//FlowEngineT* flow;

	//	PartialSatClayEngineT* flow;


	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(PartialSatClayEngine,PartialSatClayEngineT,"Engine designed to simulate the partial saturation of clay and associated swelling.",
	((Real,lmbda,0.2,,"Lambda parameter for Van Genuchten model. Free swelling 0.4. If porosity is distributed, this value becomes cell based."))
	((Real,pAir,0,,"Air pressure for calculation of capillary pressure (Pair - Pwater)"))
	((Real,Po,1.5e6,,"Po parameter for Van Genuchten model, Free swelling 0.04e6. If porosity is distributed, this value becomes cell based."))
	((Real,a,6.8,,"parameter a for evolution of Po as a function of porosity"))
	((Real,b,-1.5,,"parameter b for evolution of lambda as a function of porosity"))
	((Real,meanInitialPorosity,-1,, " if not negative, activates stochastic distribution for porosity. mean value of porosity for specimen"))
	((Real,maxPorosity,-1,, "max porosity found during stochastic poro distribution. used for evolution of porosity"))
	((Real,lambdaWeibullShape,6.,,"shape of weibull distribution of the correction factor used for porosity distribution."))
	((Real,kappaWeibullScale,1.,, "scale of weibull dist, this is the mean correction factor multiplied by :yref:`meanInitialPorosity<PartialSatClayEngine.meanInitialPorosity>`"))
	((Real,meanPoreSizeDiameter, 0.013e-6,, "mean pore size diameter, used for stochastic generation of porosity field"))
	((Real,phiZero,0.25,, "initial porosity of the specimen" ))
	((Real,nUnsatPerm,0,,"n parameter for empirical relative saturation based permeability relationship. Off by default. n=5 in Mokni2016b"))
	((Real,SrM,0,,"residual saturation for empirical relative saturation based permeability relationship"))
	((Real,SsM,1.,,"saturated saturation for empirical relative saturation based permeability relationship"))
	((Real,bIntrinsicPerm,0,,"b parameter for dependency of intrinsic permeability on macroporosity Gens 2011. Not active if 0 (default). Mokni2016 uses 8"))
	((bool,useKozeny,false,,"use Kozeny for determining the permeability based on porosity (off by default)"))
	((bool,partialSatEngine,1,,"Activates the partial sat clay engine"))
	((bool,allCellsFractured,0,,"use to simulate all pores fractured for debugging purposes only"))
	((Real,manualCrackPerm,0,,"If >0, it overrides the crack perm calculations (useful for setting cracked cells to extremely low perms to avoid fluid movement)"))
	((Real,permClamp,0,,"If >0, it prevents any permeabilities from increasing beyond this value (useful in case of very close cells"))
	((Real,crackAreaFactor,1,,"Factors the area used for crack geometry computations and capillary force removal inside cracks"))
	((Real,permAreaFactor,1,,"Factors the area used for representing roughness in cracks that still conduct flux."))
	//((bool,freezeRadii,0,,"use to stop swelling debugging purposes only"))
	((bool,crackModelActive,0,,"Activates the parallel plate approximation model for facets connected to cohesionBroken edges"))
	((Real,alpham,2.6048e-08,,"alpha parameter for particle volumetric strain model MPa^-1"))
	((Real,betam,2.10206e-08,,"beta parameter for particle volumetric strain model MPa^-1"))
	((bool,particleSwelling,1,,"set false to neglect particle swelling"))
	((bool,freeSwelling,1,,"if true, boundary forces are computed with pAir pressure only"))
	((bool,forceConfinement,false,,"If true, all the boundary particles are locked in place to simulate perfect constant volume."))
	((Real,totalVolChange,0,,"tracks the total volumetric strain that occured in each step"))
	((Real,matricSuctionRatio,1,,"The ratio of matric:osmotic suction. Facet forces computed for matricSuction fraction only."))
	((Real,residualAperture,0.,,"residual aperture of induced cracks"))
	((Real,apertureFactor,1.,,"factor to consider tortuosity"))
	((bool,calcCrackArea,true,,"The amount of crack per pore is updated if calcCrackArea=True"))
	((Real,microStructureE,1e6,,"The amount of crack per pore is updated if calcCrackArea=True"))
	((Real,microStructureNu,0.3,,"The amount of crack per pore is updated if calcCrackArea=True"))
	((Real,microStructurePhi,18.,,"The amount of crack per pore is updated if calcCrackArea=True"))
	((Real,microStructureRho,2600,,"The amount of crack per pore is updated if calcCrackArea=True"))
	((Real,microStructureAdh,6e6,,"Adhesion between microstructure particles"))
	((bool,swelling,true,,"turn just particle swelling off (for debug)"))
	((bool,suction,true,,"turn just particle suction off (for debug)"))
	((bool,volumes,true,,"turn just particle volumes off (for debug)"))
	((Real,minMicroRadFrac,0.1,,"Used during sphere insertion checks, if inscribed sphere contacts facet it cannot be reduced further than minMicroRadFrac*originalInscribedRadius"))
	((Real,minPoroClamp,0.1,,"min clamp for distribution of porosity"))
	((Real,maxPoroClamp,0.8,,"max clamp for distribution of porosity. Value over 0.8 messes with water retention curve"))
	((bool,freezePorosity,false,,"useful for freezing porosity values during stage for reaching initial conditions where volume changes should not impact porosity"))
	((bool,resetVolumeSolids,false,,"useful if genesis process was used to reach an initial condition. We don't want the volume changes that occured during geneis to affect porosity evolution."))
	((Real,Ka,101e3,,"bulk modulus of air used for equivalent compressibility model"))
	((Real,Kw,2.15e9,,"bulkmodulus of water used for equivalent compressibility model"))
	((Real,Ks,21.2e9,,"bulkmodulus of solid used for equivalent compressibility model"))
	((bool,useKeq,0,,"use the equivalent bulkmodulus for pressure field"))
	((bool,constantPorosity,0,,"use the meanInitialPorosity everywhere instead of random distribution"))
	((bool,resetOriginalParticleValues,0,,"use to reset initial volume and radii values for particles."))
	((short int,totalCracks,0,,"total discretely connected cracks."))
	((bool,onlyFractureExposedCracks,false,,"if true, only the exposed cracks have tricked permeability."))
	((bool,freezeSaturation,false,,"if true, saturation will not change in specimen."))
	((Real,fracPorosity,0.5,,"porosity value used for cracked cells"))
	((bool,displacementBasedCracks,false,,"fracture criteria will be based on displacement instead of broken bond status"))
	((bool,changeCrackSaturation,false,,"if cell becomes cracked, its saturation is reduced to residual saturation (warning this is not conservative). Useful for reducing partial sat permeability in these cells."))
	//((Real,alphaBoundSat,0,,"holds boundary saturation for transfering between triangulations"))
	// less useful items below
	((Real,betaLaplaceShape,0.23,,"shape of laplace distribution used for porosity distribution"))
	((Real,alphaExpRate,0.353,, "rate of exponential distribution for porosity distribution"))
	((Real,betaExpRate,0.353,, "rate of exponential distribution for porosity distribution"))
	((Real,totalSpecimenVolume,0,,"report the total specimen volume"))
	((int, crackedCellTotal,0,,"total number of cracked cells"))
	((string,imageryFilePath,"none",,"path to the porosity grid extracted from imagery"))
	((Real, mineralPoro,0,,"If >0, all cell with porosity below this threshold will be blocked from flow and any associated particles will be clumped together"))
	((Real, crackCellPoroThreshold,0,,"If >0, any cell above this porosity will follow crack logic from the beginning. (~0.35 for pellet imagery)"))
	((Real,blockCellPoroThreshold,0,,"If >0, any cell above this porosity will be blocked from the beginning (partially sat crack should not participate)."))
	((bool,blockIsoCells,true,,"search for cells that might be surrounded by blocked (minerals or cracks) and block them to avoid numerical instabilities."))
	((bool,brokenBondsRemoveCapillaryforces,false,,"if true, broken bonds will also remove any capillary forces associated with the area of the crack"))
	((bool,directlyModifySatFromPoro,false,,"if true, changes in porosity are used to directly change porosity. Normally, the water retention curve is taking care of this on its own."))
	((Real,partialSatDT,0,,"time step used for partial sat engine. If >0, the engine will only activate once every partialSatDT/scene->dt steps. Hydromechanical forces estimated and added as persistant forces to particles during non partial sat time steps. This value is not exact, see :yref:`PartialSatClayEngine.collectedDT`"))
	((int,elapsedIters,0,,"number of mechanical iters since last flow iter."))
	((Real,collectedDT,0,,"this is the exact time step that is computed, it enables the stiffness timestep estimate to change dynamically while maintaining an exact match for the flow timestep"))
	((bool,fracBasedPointSuctionCalc,0,,"if true, the suction per material point is computed based on fraction shared by incident cell."))
	((Real,minCellVol,0,,"Use for avoiding 0 volume cells that will interupt solution of linear system."))
	((Real,maxPo,2e6,,"Certain boundary situations where a low volume will develop and interpolate from a cell with high initial porosity leading to Po exponential estimate blowing up."))
	((Real,minLambdao,0,,"Maybe unnecessary since the lambdao function is a decay exponential for same situation described in maxPo above"))
	((Real,minParticleSwellFactor,0.1,,"If prevents particles from decreasing too far as their saturation decreases."))
	((Real,waterSurfaceTension,7.28e-2,,"Water surface tension at 20 degC used to determine entry pressure to cracks"))
	((bool,computeFracturePaths,1,,"if true, fracture paths connecting to boundary conditions will become pcondition cells and forces will be computed using atmospheric pressure."))
	((bool,useOpeningPressure,0,,"if true, cracks will be created based on local opening pressure criteria computed by waterSurfaceTension/aperture"))
	((bool,getGasPerm,0,,"If true, a gas permeability will be extracted during next timestep. This involves building another triangulation with a new conductivity matrix, factorizing the matrix, and solving i.e. this will double computational effort if performed every step."))
	((bool,gasPermFirst,0,,"Set true each time you want a new gas perm estimate."))
	((Real,airViscosity,1.8e-5,,"Used with :yref:`PartialSatClayEngine::getGasPerm` for crack permeability estimates."))
	((Real,timeDimension,0,,"Used to determine stability of system, partialSatEngine computes this value automatically."))
	((bool,useForceForCracks,0,,"Cracks are only considered if a normal force of 0 is encountered between two particles."))
	((Real,homogeneousSuctionValue,0,,"Will override the pressure solver and set all cells to the user provided value. Meant for testing non transient swelling conditions."))

	,/*PartialSatClayEngineT()*/,
	solver = shared_ptr<FlowSolver> (new FlowSolver);
	gas_solver = shared_ptr<FlowSolver> (new FlowSolver);
	,
	.def("getBoundaryGasFlux",&PartialSatClayEngine::getBoundaryGasFlux,(boost::python::arg("boundary")),"Get total Gas flux through boundary defined by its body id.\n\n.. Only useful if also useing :yref:`PartialSatClayEngine::getGasPerm` correctly. note:: The flux may be not zero even for no-flow condition. This artifact comes from cells which are incident to two or more boundaries (along the edges of the sample, typically). Such flux evaluation on impermeable boundary is just irrelevant, it does not imply that the boundary condition is not applied properly.")
	.def("setCellSaturation",&PartialSatClayEngine::setCellSaturation,(boost::python::arg("id"),boost::python::arg("temperature")),"set temperature in cell 'id'.")
	.def("getCellSaturation",&PartialSatClayEngine::getCellSaturation,(boost::python::arg("pos")),"Measure cell saturation in position pos[0],pos[1],pos[2]")
	.def("getCellVelocity", &PartialSatClayEngine::getCellVelocity,(boost::python::arg("pos")),"Get relative cell velocity at position pos[0] pos [1] pos[2].")
	.def("getCellGasVelocity", &PartialSatClayEngine::getCellGasVelocity,(boost::python::arg("pos")),"Get relative cell gas velocity at position pos[0] pos [1] pos[2]. Can only be used if :yref:`PartialSatEngine::getGasPerm`=True.")
	.def("getCellGasPImposed", &PartialSatClayEngine::getCellGasPImposed,(boost::python::arg("id")),"Get pressure condition of gas cell with id. Can only be used if :yref:`PartialSatEngine::getGasPerm`=True.")
	.def("getCellGasCenter", &PartialSatClayEngine::getCellGasCenter,(boost::python::arg("id")),"Get cell center of gas mesh with id. Can only be used if :yref:`PartialSatEngine::getGasPerm`=True.")
	.def("getCellGasVolume", &PartialSatClayEngine::getCellGasVolume,(boost::python::arg("id")),"Get volume of gas cell with id. Can only be used if :yref:`PartialSatEngine::getGasPerm`=True.")
	.def("nGasCells", &PartialSatClayEngine::nGasCells,"Get number of cells in gas mesh. Can only be used if :yref:`PartialSatEngine::getGasPerm`=True.")
	//.def("getBoundaryVel",&TemplateFlowEngine_@TEMPLATE_FLOW_NAME@::getBoundaryVel,(boost::python::arg("boundary")),"Get total avg cell velocity associated with boundary defined by its body id.")
	.def("getCrackArea",&PartialSatClayEngine::getCrackArea,"get the total cracked area.")
	.def("getCrackFabricVector",&PartialSatClayEngine::getCrackFabricVector,"get the crack fabric vector.")
	.def("getAverageAperture",&PartialSatClayEngine::getAverageAperture,"get the averageaperture.")
	.def("printPorosity",&PartialSatClayEngine::printPorosity,(boost::python::arg("file")="./porosity"),"save the porosity of the cell network.")
	.def("getCrackVolume",&PartialSatClayEngine::getCrackVolume,"get the total cracked volume.")
	.def("getWaterVolume",&PartialSatClayEngine::getWaterVolume,"get the total water volume (entered cracks only).")
	.def("getNumCracks",&PartialSatClayEngine::getNumCracks,"get the number of cracks.")
	.def("getCellPorosity",&PartialSatClayEngine::getCellPorosity,(boost::python::arg("pos")),"Measure cell porosity in position pos[0],pos[1],pos[2].")
	.def("getCellVolume",&PartialSatClayEngine::getCellVolume,(boost::python::arg("pos")),"Get cell volume in position pos[0],pos[1],pos[2].")
	.def("getCellCracked",&PartialSatClayEngine::getCellCracked,(boost::python::arg("pos")),"Get cell cracked in position pos[0],pos[1],pos[2].")
//	.def("getCellSaturation",&TwoPhaseFlowEngine::getCellSaturation,"Get saturation of cell")
	.def("getAverageSaturation",&PartialSatClayEngine::getAverageSaturation,"Get average saturation of entire specimen.")
	.def("getAverageSuction",&PartialSatClayEngine::getAverageSuction,"Get average suction of entire specimen.")
	.def("getEnteredThroatRatio",&PartialSatClayEngine::getEnteredThroatRatio,"Get ratio of entered to total cracked cells.")
	.def("artificialParticleSwell",&PartialSatClayEngine::artificialParticleSwell,(boost::python::arg("volStrain")),"Artificially swell all particles by the strain provided during next time step. Does not reactivate itself for next time step, user must call for each timestep they want to use it.")
	#ifdef YADE_VTK
	.def("saveUnsatVtk",&PartialSatClayEngine::saveUnsatVtk,(boost::python::arg("folder")="./VTK",boost::python::arg("withBoundaries")=false),"Save pressure and saturation field in vtk format. Specify a folder name for output. The cells adjacent to the bounding spheres are generated conditionally based on :yref:`withBoundaries` (not compatible with periodic boundaries)")
	//.def("saveFractureNetworkVTK",&PartialSatClayEngine::saveFractureNetworkVTK,(boost::python::arg("fileName")="./VTK"),"Save fracturenetwork as connections between cell centers")
	.def("savePermeabilityNetworkVTK",&PartialSatClayEngine::savePermeabilityNetworkVTK,(boost::python::arg("fileName")="./VTK"),"Save permeability network as connections between cell centers")
	#endif
	.def("insertMicroPores",&PartialSatClayEngine::insertMicroPores,(boost::python::arg("fracMicroPores")),"run to inscribe spheres in a desired fraction of existing pores.")
	.def("reloadSolver",&PartialSatClayEngine::reloadSolver,"use after reloading a partialSat simulation and before running next step")
	//.def("resetOriginalParticleValues",&PartialSatClayEngine::resetOriginalParticleValues,"reset the initial valume and radii values for particles.")
	.def("getTotalSpecimenVolume",&PartialSatClayEngine::getTotalSpecimenVolume,"get the total specimen volume")
	)
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(PartialSatClayEngine);


class PartialSatState : public State {
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(PartialSatState,State,"Hertz mindlin state information about each body. Only active if partially saturated clay model is active.",

		((Real,suctionSum,0,,"sum of suctions associated with incident cells"))
		((Real,suction,0,,"suction computed for particle (sum(sat of inc. cells)/num inc. cells)"))
		((Real,radiiChange,0,,"total change of particle radius due to swelling"))
		((Real,radiiOriginal,0,,"original particle radius prior to swelling"))
		((int,incidentCells,0,,"number of incident cells"))
		((int,lastIncidentCells,0,,"number of incident cells"))
		((Real,volumeOriginal,0,,"original particle volume stored for strain increments"))

		,
		createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(PartialSatState, State);
};
REGISTER_SERIALIZABLE(PartialSatState);

class PartialSatMat : public FrictMat {
public:
	virtual shared_ptr<State> newAssocState() const { return shared_ptr<State>(new PartialSatState); }
	virtual bool              stateTypeOk(State* s) const { return (bool)dynamic_cast<PartialSatState*>(s); }

	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(PartialSatMat,FrictMat,"Material used for :yref:`PartialSatClayEngine`. Necessary for the custom PartialSatState.",
		((int,num,0,,"Particle number"))
		,
		createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(PartialSatMat, FrictMat);
};
REGISTER_SERIALIZABLE(PartialSatMat);

class Ip2_PartialSatMat_PartialSatMat_MindlinPhys : public IPhysFunctor {
public:
	virtual void go(const shared_ptr<Material>& b1, const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(PartialSatMat, PartialSatMat);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(
			Ip2_PartialSatMat_PartialSatMat_MindlinPhys,IPhysFunctor,"PartialSat variant of HertzMindlin \n\n \
Calculate some physical parameters needed to obtain \
the normal and shear stiffnesses according to the Hertz-Mindlin formulation (as implemented in PFC).\n\n\
Viscous parameters can be specified either using coefficients of restitution ($e_n$, $e_s$) or viscous \
damping ratio ($\\beta_n$, $\\beta_s$). The following rules apply:\n#. If the $\\beta_n$ ($\\beta_s$) \
ratio is given, it is assigned to :yref:`MindlinPhys.betan` (:yref:`MindlinPhys.betas`) directly.\n#. \
If $e_n$ is given, :yref:`MindlinPhys.betan` is computed using $\\beta_n=-(\\log e_n)/\\sqrt{\\pi^2+(\\log e_n)^2}$. \
The same applies to $e_s$, :yref:`MindlinPhys.betas`.\n#. It is an error (exception) to specify both $e_n$ \
and $\\beta_n$ ($e_s$ and $\\beta_s$).\n#. If neither $e_n$ nor $\\beta_n$ is given, zero value \
for :yref:`MindlinPhys.betan` is used; there will be no viscous effects.\n#.If neither $e_s$ nor $\\beta_s$ \
is given, the value of :yref:`MindlinPhys.betan` is used for :yref:`MindlinPhys.betas` as well.\n\nThe \
$e_n$, $\\beta_n$, $e_s$, $\\beta_s$ are :yref:`MatchMaker` objects; they can be constructed from float \
values to always return constant value.\n\nSee :ysrc:`scripts/test/shots.py` for an example of specifying \
$e_n$ based on combination of parameters.",
			((Real,gamma,0.0,,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
			((Real,eta,0.0,,"Coefficient to determine the plastic bending moment"))
			((Real,krot,0.0,,"Rotational stiffness for moment contact law"))
			((Real,ktwist,0.0,,"Torsional stiffness for moment contact law"))
			((shared_ptr<MatchMaker>,en,,,"Normal coefficient of restitution $e_n$."))
			((shared_ptr<MatchMaker>,es,,,"Shear coefficient of restitution $e_s$."))
			((shared_ptr<MatchMaker>,betan,,,"Normal viscous damping ratio $\\beta_n$."))
			((shared_ptr<MatchMaker>,betas,,,"Shear viscous damping ratio $\\beta_s$."))
			((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute the friction angle of an interaction. If ``None``, minimum value is used."))
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_PartialSatMat_PartialSatMat_MindlinPhys);



} //namespaceyade
#endif
#endif //LINSOLV
#endif //FLOW_ENGINE
