This folder illustrates the use of [Law2_ScGeom_CapillaryPhys_Capillarity](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Law2_ScGeom_CapillaryPhys_Capillarity) engine to describe capillary interaction between particle pairs in YADE, for simulating partial saturation in the pendular regime.


I. YADE simulation scripts
--------------------------

Two examples of simulations using Law2_ScGeom_CapillaryPhys_Capillarity are herein provided:

- CapillaryPhys-example.py is a (simple) packing-scale simulation considering capillary interaction and gravity
- capillaryBridge.py defines and makes evolve a capillary bridge between two particles only


II. Capillary scripts
---------------------

As explained in the doc of [Law2_ScGeom_CapillaryPhys_Capillarity](https://yade-dem.org/doc/yade.wrapper.html#yade.wrapper.Law2_ScGeom_CapillaryPhys_Capillarity), so-called capillary files are required in order to use that engine. These capillary files include capillary bridges configurations data (for positiv capillary pressure values). A version of capillary files can be downloaded from [yade-data/capillaryFiles](https://gitlab.com/yade-dev/yade-data/-/tree/main/capillaryFiles) (direct download [here](https://gitlab.com/yade-dev/yade-data/-/archive/main/yade-data-main.zip?path=capillaryFiles)), they consider a zero contact (wetting) angle and some given particle radius ratios, and disregard cases where a liquid bridge between contacting spheres would extend over less than 0.7 degrees.

Also, a suite of three files are herein provided in order to enable users to build their own capillary files for any contact angle or particle radius ratio, or to study capillary bridges as such (file `solveLiqBridge.m/py`). At the moment, the suite of files is fully coded in .m, for use with MATLAB, and partly in Python, with some .py counterparts of some .m files.
Capillary files generation requires launching `writesCapFile()` (in `writesCapFile.m`) only, after providing a couple of parameters therein (see description of `writesCapFile.m` below). The two other files (`solveLiqBridge.m/py` and `solveLaplace_uc.m/py`) define functions used by `writesCapFile()`DE and do not require any user input. Specifically, the roles of the three files are as follows:


- `solveLiqBridge.m/py` solves the Laplace-Young equation for one given bridge, defined in terms of the input attributes of the `solveLiqBridge` function (see therein). That `solveLiqBridge` function is usually called by other files (see below) during capillary files generation, however it can also be executed on its own in order to study (e.g. plot) capillary bridge profile.

- `solveLaplace_uc.m/py` calls several times `solveLiqBridge` in order to compute all possible bridges configurations for given contact angle, capillary pressure, and particle radius ratio. It is usually called by writesCapFile (see below) during capillary files generation, however it can also be executed on its own. In particular, it may output text files including bridges data for one given capillary pressure.


- `writesCapFile.m`, finally, is the conductor during the capillary files generation, and is the only one that actually requires the user's attention for such task. Parameters are to be defined by the user directly in the .m file (from l. 10 to 30) before launching.


Generating a complete set of 10 capillary files typically used to require few days in terms of computation time with MATLAB R2014 and much (~ 7 \*) more with Octave. Using MATLAB R2016b drastically increases the speed of execution (\* 10-30 with respect to R2014) and the incoveniency of Octave. Python 2.7.12 may be a good compromise, with a time cost ~ 6-7\* the one of MATLAB R2016b.

*References (in code comments)*

 * [Duriez2017](https://link.springer.com/article/10.1007/s11440-016-0500-6): J. Duriez and R. Wan, Contact angle mechanical influence for wet granular soils, Acta Geotechnica, 12, 2017 ([fulltext](https://hal.archives-ouvertes.fr/hal-01868741))
 * Lian1993: G. Lian and C. Thornton and M. J. Adams, A Theoretical Study of the Liquid Bridge Forces between Two Rigid Spherical Bodies, Journal of Colloid and Interface Science, 161(1), 1993
 * Scholtes2008 (french): L. Scholtes, Modelisation Micro-Mecanique des Milieux Granulaires Partiellement Satures, PhD Thesis from Institut polytechnique de Grenoble, 2008
 * Soulie2005 (french): F. Soulie, Cohesion par capillarite et comportement mecanique de milieux granulaires, PhD Thesis from Universite Montpellier II, 2005
 * Soulie2006: F. Soulie and F. Cherblanc and M. S. El Youssoufi and C. Saix, Influence of liquid bridges on the mechanical behaviour of polydisperse granular materials, International Journal for Numerical and Analytical Methods in Geomechanics, 30(3), 2006


