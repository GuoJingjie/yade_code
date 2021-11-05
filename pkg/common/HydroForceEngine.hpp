// 2017 © Raphael Maurin <raphael.maurin@imft.fr>
// 2017 © Julien Chauchat <julien.chauchat@legi.grenoble-inp.fr>
// 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#pragma once

#include <core/PartialEngine.hpp>
#include <queue>

namespace yade { // Cannot have #include directive inside.

class HydroForceEngine : public PartialEngine {
public:
	void         initialization();
	void         averageProfile();
	void         turbulentFluctuation();
	void         turbulentFluctuationZDep();
	void         fluidResolution(Real tfin, Real dt);
	void         computeRadiusParts();
	Real         computePoreLength(int i);
	int          computeZbedIndex();
	Real         computeDiameter(int i);
	vector<Real> computePhiPartAverageOverTime();

public:
	void                     action() override;
	std::queue<vector<Real>> phiPartT;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(HydroForceEngine,PartialEngine,"Engine performing a coupling of the DEM with a volume-averaged 1D fluid resolution to simulate steady uniform unidirectional fluid flow. It has been developed and used to model steady uniform gravity-driven turbulent bedload transport [Maurin2015b]_ [Maurin2016]_ [Maurin2018]_, but can be also used in its current state for laminar or pressure-driven configurations. The fundamentals of the model can be found in [Maurin2015b]_ and [Maurin2015PhD]_, and in more details in [Maurin2018_VANSbasis]_, [Maurin2018_VANSfluidResol]_ and [Maurin2018_VANSvalidations]_. \n The engine can be decomposed in three different parts:\n (i) It applies the fluid force on the particles imposed by the fluid velocity profiles and fluid properties,\n (ii) It evaluates averaged solid depth profiles necessary for the fluid force application and for the fluid resolution,\n (iii) It solve the volume-averaged 1D fluid momentum balance. \nThe three different functions are detailed below: \n\n (i) Fluid force on particles \n Apply to each particles, buoyancy, drag and lift force due to a 1D fluid flow and can apply lubrication force between two particles. The applied drag force reads\n\n $F_{d}=\\frac{1}{2} C_d A\\rho^f|\\vec{v_f - v}| \\vec{v_f - v}$ \n\n where $\\rho$ is the fluid density (:yref:`densFluid<HydroForceEngine.densFluid>`), $v$ is particle's velocity,  $v_f$ is the velocity of the fluid at the particle center (taken from the fluid velocity profile :yref:`vxFluid<HydroForceEngine.vxFluid>`),  $A = \\pi d^2/4$ is particle projected area (disc), $C_d$ is the drag coefficient. The formulation of the drag coefficient depends on the local particle reynolds number and the solid volume fraction. The formulation of the drag is [Dallavalle1948]_ [RevilBaudard2013]_ with a correction of Richardson-Zaki [Richardson1954]_ to take into account the hindrance effect. This law is classical in sediment transport. The possibly activated lubrication force (with parameter:yref:`lubrication<HydroForceEngine.lubrication>` put to True) reads: $F_{lubrication} = \\frac{6 \\pi \\eta^f v_rel^n}{\\delta^n + \\epsilon_r}$, with $\\eta^f$ the fluid dynamic viscosity :yref:`viscoDyn<HydroForceEngine.viscoDyn>`, $v_rel^n$ the normal relative velocity of the two particles, $\\delta^n$ the distance between the two particles surface, and $\\epsilon_r$ the roughness scale of the particle (:yref:`roughnessPartScale<HydroForceEngine.roughnessPartScale>`).\n\nIt is possible to activate a fluctuation of the drag force for each particle which account for the turbulent fluctuation of the fluid velocity (:yref:`velFluct<HydroForceEngine.velFluct>`). Three simple discrete random walk model have been implemented for the turbulent velocity fluctuation. The main one (turbulentFluctuations) takes as input the Reynolds stress tensor $R^f_{xz}$ as a function of the depth, and allows to recover the main property of the fluctuations by imposing $<u_x'u_z'> (z) = <R^f_{xz}>(z)/\\rho^f$. It requires as input $<R^f_{xz}>(z)$ called :yref:`ReynoldStresses<HydroForceEngine.ReynoldStresses>` in the code. \n The formulation of the lift is taken from [Wiberg1985]_ and is such that : \n\n $F_{L}=\\frac{1}{2} C_L A\\rho^f((v_f - v)^2_{top} - (v_f - v)^2_{bottom})$ \n\n Where the subscript top and bottom means evaluated at the top (respectively the bottom) of the sphere considered. This formulation of the lift account for the difference of pressure at the top and the bottom of the particle inside a turbulent shear flow. As this formulation is controversial when approaching the threshold of motion [Schmeeckle2007]_ it is possible to desactivate it with the variable :yref:`lift<HydroForceEngine.lift>`.\n The buoyancy is taken into account through the buoyant weight : \n\n $F_{buoyancy}= - \\rho^f V^p g$ \n\n, where g is the gravity vector along the vertical, and $V^p$ is the volume of the particle. In the case where the fluid flow is steady and uniform, the buoyancy reduces to its wall-normal component (see [Maurin2018]_ for a full explanation), and one should put :yref:`steadyFlow<HydroForceEngine.steadyFlow>` to true in order to kill the streamwise component. \n\n (ii)  Averaged solid depth profiles\n The function averageProfile evaluates the volume averaged depth profiles (1D) of particle velocity, particle solid volume fraction and particle drag force. It uses a volume-weighting average following [Maurin2015PhD]_[Maurin2015b]_, i.e. the average of a variable $A^p$ associated to particles at a given discretized wall-normal position $z$ is given by: \n\n $\\left< A \\right>^s(z) = \\displaystyle \\frac{\\displaystyle \\sum_{p|z^p\\in[z-dz/2,z+dz/2]}  A^p(t) V^p_z}{\\displaystyle \\sum_{p|z^p\\in[z-dz/2,z+dz/2]}  V^p_z}$\n\n Where the sums are over the particles contained inside the slice between the wall-normal position $z-dz/2$ and $z+dz/2$, and $V^p$ represents the part of the volume of the given particle effectively  contained inside the slice. For more details, see [Maurin2015PhD]_. \n\n (iii) 1D volume-average fluid resolution\n The fluid resolution is based on the resolution of the 1D volume-averaged fluid momentum balance. It assumes by definition (unidirectional) that the fluid flow is steady and uniform. It is the same fluid resolution as [RevilBaudard2013]_. Details can be found in this paper and in [Maurin2015PhD]_ [Maurin2015b]_.\n\n The three different component can be used independently, e.g. applying a fluid force due to an imposed fluid profile or solving the fluid momentum balance for a given concentration of particles.",
		//// General parameters
                ((Vector3r,gravity,Vector3r(0,0,-9.81),,"Gravity vector"))
		((bool,convAccOption,false,,"To activate the convective acceleration option in order to account for a convective acceleration term inside the momentum balance. "))
		((vector<Real>,radiusParts,,,"Variables containing the number of different radius of particles in the simulation. Allow to perform class averaging by particle size."))
		((bool,compatibilityOldVersion,false,,"Option to make HydroForceEngine compatible with former scripts. Slow down slightly the calculation and will eventually be removed."))
		//// Mesh parameters
		((int,nCell,1,,"Number of cell in the depth"))
		((Real,zRef,0.,,"Position of the reference point which correspond to the first value of the fluid velocity, i.e. to the ground."))
		((Real,deltaZ,,,"Height of the discretization cell."))
		((Real,vCell,,,"Volume of averaging cell"))
		///  Fluid Resolution parameters
		((vector<Real>,vxFluid,,,"Discretized streamwise fluid velocity depth profile at t"))
		((Real,densFluid,1000,,"Density of the fluid, by default - density of water"))
		((Real,viscoDyn,1e-3,,"Dynamic viscosity of the fluid, by default - viscosity of water"))
		((Real,radiusPart,0.,,"Reference particle radius"))
		((Real,dpdx,0.,,"pressure gradient along streamwise direction"))
		((vector<Real>,turbulentViscosity,,,"Fluid Resolution: turbulent viscocity as a function of the depth"))
		((Real,phiMax,0.64,,"Fluid resolution: maximum solid volume fraction. "))
		((int,irheolf,0,,"Fluid resolution: effective fluid viscosity option: 0: pure fluid viscosity, 1: Einstein viscosity. "))
		((int,iturbu,1,,"Fluid resolution: activate the turbulence resolution, 1, or not, 0"))
		((int,ilm,2,,"Fluid resolution: type of mixing length resolution applied: 0: classical Prandtl mixing length, 1: Prandtl mixing length with free-surface effects, 2: Damp turbulence accounting for the presence of particles [Li1995]_, see [RevilBaudard2013]_ for more details."))
		((int,iusl,1,,"Fluid resolution: option to set the boundary condition at the top of the fluid, 0: Dirichlet, fixed ($u=uTop$ en $z=h$), 1:Neumann, free-surface ($du/dz=0$ en $z=h$)."))
		((Real,uTop,1.,,"Fluid resolution: fluid velocity at the top boundary when iusl = 0"))
		((Real,kappa,0.41,,"Fluid resolution: Von Karman constant. Can be tuned to account for the effect of particles on the fluid turbulence, see e.g. [RevilBaudard2015]_"))
		((int,viscousSubLayer,0,,"Fluid resolution: solve the viscous sublayer close to the bottom boundary if set to 1"))
		((bool,fluidWallFriction,false,,"Fluid resolution: if set to true, introduce a sink term to account for the fluid friction at the wall, see [Maurin2015]_ for details. Requires to set the width of the channel. It might slow down significantly the calculation."))
		((Real,fluidFrictionCoef,1.,,"Fluid resolution: fitting coefficient for the fluid wall friction"))
		((Real,channelWidth,1.,,"Fluid resolution: Channel width for the evaluation of the fluid wall friction inside the fluid resolution."))
		((Real,phiBed,0.08,,"Turbulence modelling parameter. Associated with mixing length modelling ilm = 5."))
		//// Fluid-particle interactions 
		((Real,expoRZ,3.1,,"Value of the Richardson-Zaki exponent, for the drag correction due to hindrance"))
		((bool,lift,false,,"Option to activate or not the evaluation of the lift"))
		((Real,Cl,0.2,,"Value of the lift coefficient taken from [Wiberg1985]_"))
		((vector<Real>,convAcc,,,"Convective acceleration, depth dependent"))
		((vector<Real>,taufsi,,,"Fluid Resolution: Create Taufsi/rhof = dragTerm/(rhof(vf-vxp)) to transmit to the fluid code"))
		((bool,steadyFlow,true,,"Condition to modify the buoyancy force according to the physical difference between a fluid at rest and a steady fluid flow. For more details see [Maurin2018]_"))
		((bool,lubrication,false,,"Condition to activate the calculation of the lubrication force."))
		((Real,roughnessPartScale,1e-3,,"Roughness length scale of the particle. In practice, the lubrication force is cut off when the two particles are at a distance roughnessPartScale."))
		//// Particle averaged depth profiles parameters Monodisperse: 
		((vector<Real>,phiPart,,,"Discretized solid volume fraction depth profile. Can be taken as input parameter or evaluated directly inside the engine, calling from python the averageProfile() function"))
		((vector<Vector3r>,vPart,,,"Discretized streamwise solid velocity depth profile, in x, y and z direction. Only the x direction measurement is taken into account in the 1D fluid coupling resolution. The two other can be used as output parameters. The x component can be taken as input parameter, or evaluated directly inside the engine, calling from python the averageProfile() function"))
		((unsigned int,nbAverageT,0,,"If >0, perform a time-averaging (in addition to the spatial averaging) over nbAverage steps."))
		((bool,pointParticleAverage,false,,"Evaluate the averaged with a point particle method. If False, consider the particle extent and weigth the averaged by the volume contained in each averaging cell."))
		// Polydisperse
		((bool,enableMultiClassAverage,false,,"Enables specific averaging for all the different particle size. Uses a lot of memory if using a lots of different particle size"))
		((vector<vector<Real>>,multiVxPart,,,"Spatial-averaged velocity in x direction for each class of particle."))
		((vector<vector<Real>>,multiVyPart,,,"Spatial-averaged velocity in y direction for each class of particle."))
		((vector<vector<Real>>,multiVzPart,,,"Spatial-averaged velocity in z direction for each class of particle."))
		((vector<vector<Real>>,multiPhiPart,,,"Spatial-averaged solid volume fraction for each class of particle."))
		//// Discrete Random Walk fluid velocity fluctuations model parameters
		((bool,velFluct,false,,"If true, activate the determination of turbulent fluid velocity fluctuation for the next time step only at the position of each particle, using a simple discrete random walk (DRW) model based on the Reynolds stresses profile (:yref:`ReynoldStresses<HydroForceEngine.ReynoldStresses>`)"))
		((vector<Real>,vFluctX,,,"Vector associating a streamwise fluid velocity fluctuation to each particle. Fluctuation calculated in the C++ code from the discrete random walk model"))
		((vector<Real>,vFluctY,,,"Vector associating a spanwise fluid velocity fluctuation to each particle. Fluctuation calculated in the C++ code from the discrete random walk model"))
		((vector<Real>,vFluctZ,,,"Vector associating a normal fluid velocity fluctuation to each particle. Fluctuation calculated in the C++ code from the discrete random walk model"))
		((vector<Real>,ReynoldStresses,,,"Vector of size equal to :yref:`nCell<HydroForceEngine.nCell>` containing the Reynolds stresses as a function of the depth. ReynoldStresses(z) $=  \\rho^f <u_x'u_z'>(z)^2$"))
		((Real,bedElevation,0.,,"Elevation of the bed above which the fluid flow is turbulent and the particles undergo turbulent velocity fluctuation."))
		((vector<Real>,fluctTime,,,"Vector containing the time of life of the fluctuations associated to each particles."))
		((Real,dtFluct,,,"Execution time step of the turbulent fluctuation model."))
		((bool,unCorrelatedFluctuations,false,,"Condition to generate uncorrelated fluid fluctuations. Default case represent in free-surface flows, for which the vertical and streamwise fluid velocity fluctuations are correlated (see e.g. reference book of Nezu & Nagakawa 1992, turbulence in open channel flows)."))

		// Old parameters, not necessary anymore, kept for compatibility with older scripts
		((bool,twoSize,false,,"Not maintained anymore. Option to activate when considering two particle size in the simulation. When activated evaluate the average solid volume fraction and drag force for the two type of particles of diameter diameterPart1 and diameterPart2 independently."))
		((Real,radiusPart1,0.,,"Radius of the particles of type 1. Useful only when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((Real,radiusPart2,0.,,"Radius of the particles of type 2. Useful only when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,phiPart1,,,"Discretized solid volume fraction depth profile of particles of type 1. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,phiPart2,,,"Discretized solid volume fraction depth profile of particles of type 2. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,averageDrag1,,,"Discretized average drag depth profile of particles of type 1. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,averageDrag2,,,"Discretized average drag depth profile of particles of type 2. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vxPart1,,,"Discretized solid streamwise velocity depth profile of particles of type 1. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vxPart2,,,"Discretized solid streamwise velocity depth profile of particles of type 2. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vyPart1,,,"Discretized solid spanwise velocity depth profile of particles of type 1. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vyPart2,,,"Discretized solid spanwise velocity depth profile of particles of type 2. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vzPart1,,,"Discretized solid wall-normal velocity depth profile of particles of type 1. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vzPart2,,,"Discretized solid wall-normal velocity depth profile of particles of type 2. Evaluated when :yref:`twoSize<HydroForceEngine.twoSize>` is set to True."))
		((vector<Real>,vxPart,,,"Discretized streamwise solid velocity depth profile. Can be taken as input parameter, or evaluated directly inside the engine, calling from python the averageProfile() function"))
		((vector<Real>,vyPart,,,"Discretized spanwise solid velocity depth profile. Can be taken as input parameter, or evaluated directly inside the engine, calling from python the averageProfile() function"))
		((vector<Real>,vzPart,,,"Discretized wall-normal solid velocity depth profile. Can be taken as input parameter, or evaluated directly inside the engine, calling from python the averageProfile() function"))
		((vector<Real>,averageDrag,,,"Discretized average drag depth profile. No role in the engine, output parameter. For practical reason, it can be evaluated directly inside the engine, calling from python the averageProfile() method of the engine"))
		((vector< vector<Real> >,multiDragPart,,,"Spatial-averaged mean drag force for each class of particle. Un-used ? Or just for debug."))

	,/*ctor*/
	,/*py*/
	 .def("initialization",&HydroForceEngine::initialization,"Initialize the necessary parameters to make HydroForceEngine run. Necessary to execute before any simulation run, otherwise it crashes")
	 .def("averageProfile",&HydroForceEngine::averageProfile,"Compute and store the particle velocity (:yref:`vxPart<HydroForceEngine.vxPart>`, :yref:`vyPart<HydroForceEngine.vyPart>`, :yref:`vzPart<HydroForceEngine.vzPart>`) and solid volume fraction (:yref:`phiPart<HydroForceEngine.phiPart>`) depth profile. For each defined cell z, the k component of the average particle velocity reads: \n\n $<v_k>^z= \\sum_p V^p v_k^p/\\sum_p V^p$,\n\n where the sum is made over the particles contained in the cell, $v_k^p$ is the k component of the velocity associated to particle p, and $V^p$ is the part of the volume of the particle p contained inside the cell. This definition allows to smooth the averaging, and is equivalent to taking into account the center of the particles only when there is a lot of particles in each cell. As for the solid volume fraction, it is evaluated in the same way:  for each defined cell z, it reads: \n\n $<\\phi>^z= \\frac{1}{V_{cell}}\\sum_p V^p$, where $V_{cell}$ is the volume of the cell considered, and $V^p$ is the volume of particle p contained in cell z.\n This function gives depth profiles of average velocity and solid volume fraction, returning the average quantities in each cell of height dz, from the reference horizontal plane at elevation :yref:`zRef<HydroForceEngine.zRef>` (input parameter) until the plane of elevation :yref:`zRef<HydroForceEngine.zRef>` plus :yref:`nCell<HydroForceEngine.nCell>` times :yref:`deltaZ<HydroForceEngine.deltaZ>` (input parameters). When the option :yref:`twoSize<HydroForceEngine.twoSize>` is set to True, evaluate in addition the average drag (:yref:`averageDrag1<HydroForceEngine.averageDrag1>` and :yref:`averageDrag2<HydroForceEngine.averageDrag2>`) and solid volume fraction (:yref:`phiPart1<HydroForceEngine.phiPart1>` and :yref:`phiPart2<HydroForceEngine.phiPart2>`) depth profiles considering only the particles of radius respectively :yref:`radiusPart1<HydroForceEngine.radiusPart1>` and :yref:`radiusPart2<HydroForceEngine.radiusPart2>` in the averaging.")
	 .def("fluidResolution",&HydroForceEngine::fluidResolution,"Solve the 1D volume-averaged fluid momentum balance on the defined mesh (:yref:`nCell<HydroForceEngine.nCell>`, :yref:`deltaZ<HydroForceEngine.deltaZ>`) from the volume-averaged solid profiles (:yref:`phiPart<HydroForceEngine.phiPart>`,:yref:`vxPart<HydroForceEngine.vxPart>`,:yref:`averageDrag<HydroForceEngine.averageDrag>`), which can be evaluated with the averageProfile function.")
	 .def("turbulentFluctuation",&HydroForceEngine::turbulentFluctuation,"Apply a discrete random walk model to the evaluation of the drag force to account for the fluid velocity turbulent fluctuations. Very simple model applying fluctuations from the values of the Reynolds stresses in order to recover the property $<u_x'u_z'> (z) = <R^f_{xz}>(z)/\\rho^f$. The random fluctuations are modified over a time scale given by the eddy turn over time.")
	 .def("turbulentFluctuationZDep",&HydroForceEngine::turbulentFluctuationZDep,"Apply turbulent fluctuation to the problem similarly to turbulentFluctuation but with an update of the fluctuation depending on the particle position.")
	 .def("computeRadiusParts", &HydroForceEngine::computeRadiusParts, "compute the different class of radius present in the simulation.")
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(HydroForceEngine);

} // namespace yade
