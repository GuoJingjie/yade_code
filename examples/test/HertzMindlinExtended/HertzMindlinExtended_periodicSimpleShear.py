# encoding: utf-8

# script for periodic simple shear test, with periodic boundary
# first compresses to attain some isotropic stress (checkStress),
# then loads in shear (checkDistorsion)
#
# the initial packing is either regular (hexagonal), with empty bands along the boundary,
# or periodic random cloud of spheres
#
#--------------------------------------------------------------------------------------------------
# standard example is adapted to use Hertz-Mindlin contact law OR Conical Damage Model
# Hertz-Minlin: purly elastic contact
# Conical Damage Model: elastic-plastic contact, plus pressure dependent inter-particle friction coefficient
#--------------------------------------------------------------------------------------------------
# for implementation of Conical Damage Model see:
# Suhr, B. & Six, K. Parametrisation of a DEM model for railway ballast under different load cases Granular Matter, 2017, 19, 64
# Suhr & Six 2016: On the effect of stress dependent interparticle friction in direct shear tests, Powder Technology , (2016) 294:211 - 220, https://arxiv.org/abs/1810.03308
#--------------------------------------------------------------------------------------------------
# Conical Damage Model details:
# normal contact:
#   sigmaMax=1e99 -> Hertzian contact, purely elastic
#   sigmaMax finite,  0<= alpha <= pi -> elastic-plastic contact
# tangential contact:
#   c1,c2 >0: usage of pressure dependent inter-particle friction coefficient
#   c1=c2=0.0: constant inter-particle friction coefficient

# setup the periodic boundary
O.periodic = True
O.cell.refSize = (2, 2, 2)

from yade import pack, plot
import numpy as np

#------------------------choose contact law:
#makeCDM=1 use Conical Damage Model
#makeCDM=0 use Hertz-Mindlin model
makeCDM = 1

#specify friction parameters to be used
#use constant interparticle friction coefficient:
#frict=0.5
#c1=0.0 # write zero here ...
#c2=0.0 # ... and here
##use pressure dependent inter-particle friction coefficient:
frict = .25
c1 = 1.3
c2 = 1.26e-5

if makeCDM == 0:
	O.materials.append(FrictMat(young=1e7, poisson=0.3, frictionAngle=np.arctan(frict), density=1000.0, label="Material1"))
else:
	#O.materials.append(RoughFrictMat(young=E, poisson=0.28, frictionAngle=0.01, density=2600.,jrc=0.35, jcs=1e6, phiR=0.52, A0=A0, label="ballast"))#ballast
	O.materials.append(
	        FrictMatCDM(
	                young=1e7,
	                poisson=0.3,
	                frictionAngle=np.arctan(frict),
	                density=1000.0,
	                sigmaMax=5e5,
	                alpha=np.radians(82.0),
	                c1=c1,
	                c2=c2,
	                label="Material1"
	        )
	)  #ballast

#for both contact laws:
#ATTENTION geometric functors for Sphere-Wall Sphere-Facet contacts use Hertzian theory
GeomFunctorList = [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(hertzian=True), Ig2_Wall_Sphere_ScGeom(hertzian=True)]

# the "if 0:" block will be never executed, therefore the "else:" block will be
# to use cloud instead of regular packing, change to "if 1:" or something similar
if 0:
	# create cloud of spheres and insert them into the simulation
	# we give corners, mean radius, radius variation
	sp = pack.SpherePack()
	sp.makeCloud((0, 0, 0), (2, 2, 2), rMean=.1, rRelFuzz=.1, periodic=True, seed=1)
	# insert the packing into the simulation
	sp.toSimulation(color=(0, 0, 1))  # pure blue
else:
	# in this case, add dense packing
	O.bodies.append(pack.regularHexa(pack.inAlignedBox((0, 0, 0), (2, 2, 2)), radius=.1, gap=0, color=(0, 0, 1)))

# simulation loop (will be run at every step)
O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb()]),
        InteractionLoop(
                #[Ig2_Sphere_Sphere_ScGeom()],
                #[Ip2_FrictMat_FrictMat_FrictPhys()],
                #[Law2_ScGeom_FrictPhys_CundallStrack()]
                GeomFunctorList,
                [Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM(),
                 Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM(),
                 Ip2_FrictMat_FrictMat_MindlinPhys()],
                [Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM(label="lawHMCDM"),
                 Law2_ScGeom_MindlinPhys_Mindlin(label="lawHM")]
        ),
        NewtonIntegrator(damping=.4),
        # run checkStress function (defined below) every second
        # the label is arbitrary, and is used later to refer to this engine
        PyRunner(command='checkStress()', iterPeriod=100, label='checker'),
        # record data for plotting every 100 steps; addData function is defined below
        PyRunner(command='addData()', iterPeriod=100, label='plotter', dead=True)
]

# set the integration timestep to be 1/2 of the "critical" timestep
O.dt = .5 * PWaveTimeStep()

# prescribe isotropic normal deformation (constant strain rate)
# of the periodic cell
O.cell.velGrad = Matrix3(-.1, 0, 0, 0, -.1, 0, 0, 0, -.1)

# when to stop the isotropic compression (used inside checkStress)
limitMeanStress = 5e5


# called every second by the PyRunner engine
def checkStress():
	# stress tensor as the sum of normal and shear contributions
	# Matrix3.Zero is the intial value for sum(...)
	stress = sum(normalShearStressTensors(), Matrix3.Zero)
	print('mean stress', stress.trace() / 3.)
	# if mean stress is below (bigger in absolute value) limitMeanStress, start shearing
	if stress.trace() / 3. > limitMeanStress:
		# apply constant-rate distorsion on the periodic cell
		O.cell.velGrad = Matrix3(0, 0, .1, 0, 0, 0, 0, 0, 0)
		# change the function called by the checker engine
		# (checkStress will not be called anymore)
		checker.command = 'checkDistorsion()'
		plotter.dead = False
		# block rotations of particles to increase tanPhi, if desired
		# disabled by default
		if 0:
			for b in O.bodies:
				# block X,Y,Z rotations, translations are free
				b.state.blockedDOFs = 'XYZ'
				# stop rotations if any, as blockedDOFs block accelerations really
				b.state.angVel = (0, 0, 0)


# called from the 'checker' engine periodically, during the shear phase
def checkDistorsion():
	# if the distorsion value is >.3, exit; otherwise do nothing
	if abs(O.cell.trsf[0, 2]) > .3:
		# save data from addData(...) before exiting into file
		# use O.tags['id'] to distinguish individual runs of the same simulation
		plot.saveDataTxt(O.tags['id'] + '.txt')
		# exit the program
		#import sys
		#sys.exit(0) # no error (0)
		O.pause()


# called periodically to store data history
def addData():
	# get the stress tensor (as 3x3 matrix)
	stress = sum(normalShearStressTensors(), Matrix3.Zero)
	Rmean = 0.0
	muMean = 0.0
	count = 0
	for i in O.interactions:
		Da = i.geom.refR1
		Db = i.geom.refR2
		count = count + 1
		if makeCDM:
			yieldCont = lawHMCDM.ratioYieldingContacts()
			Rmean = Rmean + i.phys.R
		else:
			yieldCont = 0
			Rmean = Rmean + (Da * Db / (Da + Db))
		muMean = muMean + i.phys.tangensOfFrictionAngle
	# give names to values we are interested in and save them
	plot.addData(
	        exz=O.cell.trsf[0, 2],
	        szz=stress[2, 2],
	        sxz=stress[0, 2],
	        tanPhi=stress[0, 2] / stress[2, 2],
	        i=O.iter,
	        averageContactRadius=Rmean / float(count),
	        PercentYieldCont=yieldCont * 100,
	        averageMu=muMean / float(count)
	)
	# color particles based on rotation amount
	for b in O.bodies:
		# rot() gives rotation vector between reference and current position
		b.shape.color = scalarOnColorScale(b.state.rot().norm(), 0, pi / 2.)


# define what to plot (3 plots in total)
## exz(i), [left y axis, separate by None:] szz(i), sxz(i)
## szz(exz), sxz(exz)
## tanPhi(i)
# note the space in 'i ' so that it does not overwrite the 'i' entry
#plot.plots={'i':('exz',None,'szz','sxz'),'exz':('szz','sxz'),'i ':('tanPhi',)}
#plot.plots={'exz':('szz','sxz'),'exz ':('R',None,'yieldCont',)}
plot.plots = {'exz': ('szz', 'sxz'), 'exz ': ('averageContactRadius'), 'exz  ': ('PercentYieldCont',), 'exz   ': ('averageMu',)}

# better show rotation of particles
Gl1_Sphere.stripes = True

# open the plot on the screen
plot.plot()

O.saveTmp()
