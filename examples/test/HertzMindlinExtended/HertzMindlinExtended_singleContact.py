#script with one single contact between a sphere and a wall/facet/sphere
#see the effect of contact laws:
#Hertz-Minlin: purly elastic contact
#Conical Damage Model: elastic-plastic contact, plus pressure dependent inter-particle friction coefficient
#--------------------------------------------------------------------------------------------------
#for implementation of Conical Damage Model see:
#Suhr, B. & Six, K. Parametrisation of a DEM model for railway ballast under different load cases Granular Matter, 2017, 19, 64
#--------------------------------------------------------------------------------------------------
# Conical Damage Model details:
# normal contact:
#   sigmaMax=1e99 -> Hertzian contact, purely elastic
#   sigmaMax finite,  0<= alpha <= pi/2 -> elastic-plastic contact
# tangential contact:
#   c1,c2 >0: usage of pressure dependent inter-particle friction coefficient
#   c1=c2=0.0: constant inter-particle friction coefficient

from yade import plot
import numpy as np

from yade import qt
qt.View()  # To visualize the system!

#------------------------choose contact law:
#makeCDM=1 use Conical Damage model
#makeCDM=0 use Hertz-Mindlin model
makeCDM = 1

#------------------------choose type of particle, which collides with the sphere:
TYPE = 'SPHERE'  #, 'FACET'# 'WALL'#

#------------------------- sphere radius and position
rad_p = 0.1
L = 1.0

#-------------------------Define materials
density = 3000
young = 1e10
poisson = 0.2

sigmaMax = 6e7  #
alpha = np.radians(85)  # in radians, 0<= alpha <= pi/2

Mat = O.materials.append(FrictMat(young=young, density=density, poisson=poisson, frictionAngle=np.arctan(0.1)))
MatCDM = O.materials.append(
        FrictMatCDM(young=young, density=density, poisson=poisson, frictionAngle=np.arctan(0.1), sigmaMax=sigmaMax, alpha=alpha, c1=0.0, c2=0.0)
)

if makeCDM == 1:
	usedMat = MatCDM
else:
	usedMat = Mat

#for both contact laws:
#ATTENTION geometric functors for Sphere-Wall Sphere-Facet contacts use Hertzian theory
GeomFunctorList = [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(hertzian=True), Ig2_Wall_Sphere_ScGeom(hertzian=True)]

#-------------------------Define engines

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb(), Bo1_Box_Aabb(),
                               Bo1_Wall_Aabb(), Bo1_Polyhedra_Aabb()], allowBiggerThanPeriod=True),
        InteractionLoop(
                GeomFunctorList, [Ip2_FrictMatCDM_FrictMatCDM_MindlinPhysCDM(),
                                  Ip2_FrictMat_FrictMatCDM_MindlinPhysCDM(),
                                  Ip2_FrictMat_FrictMat_MindlinPhys()],
                [Law2_ScGeom_MindlinPhysCDM_HertzMindlinCDM(label="lawHMCDM"),
                 Law2_ScGeom_MindlinPhys_Mindlin(label="lawHM")]
        ),
        PyRunner(command='pausing()', iterPeriod=1, label="checker"),
        PyRunner(command='addPlotData()', iterPeriod=1, label="plotter"),
        NewtonIntegrator(gravity=(0, 0, 0.0), damping=0.0)
]

#------------------- Sphere which will contact ...  ----------------------

O.bodies.append([sphere((L / 5.0, L / 2.0, -rad_p), rad_p, material=usedMat)])
Body1 = O.bodies[-1]

#------------------- ... one of the following 3 objects ----------------------

### ------------------- (I) Facetbox
if TYPE == 'FACET':
	O.bodies.append(geom.facetBox(center=(L / 2.0, L / 2.0, L / 2.0), extents=(L / 2.0, L / 2.0, L / 2.0), material=usedMat, fixed=True))
	Body2 = O.bodies[10]

### ------------------- (II) Plate
elif TYPE == 'WALL':
	O.bodies.append(utils.wall((L / 5.0, L / 2.0, 0.0), axis=2, sense=0, material=usedMat))
	Body2 = O.bodies[-1]

## ------------------- (III) Sphere
else:
	TYPE = 'SPHERE'  # otherwise choose sphere
	O.bodies.append([sphere((L / 5.0, L / 2.0, rad_p), rad_p, material=usedMat, fixed=True)])
	Body2 = O.bodies[-1]

# ------------------------- apply velovity on body1 -------------------------
vel = 0.05
O.bodies[Body1.id].state.vel = Vector3(0, 0, vel)


#-------------------------Define pyRunners
def pausing():
	if (O.iter >= stopAtIter - 1):
		print("pause at iter: ", O.iter)
		if makeCDM == 1:
			plot.saveDataTxt('CDMsingle' + TYPE + '.txt', vars=('i', 'overlap', 'Fz', 'yieldCont', 'contactRadius', 'stress'))
		else:
			plot.saveDataTxt('HMclassicSingle' + TYPE + '.txt', vars=('i', 'overlap', 'Fz', 'yieldCont', 'contactRadius', 'stress'))
		O.pause()


def addPlotData():
	if O.interactions[Body2.id, Body1.id].geom:
		Fz = O.forces.f(Body2.id)[2]
		Da = O.interactions[Body2.id, Body1.id].geom.refR1
		Db = O.interactions[Body2.id, Body1.id].geom.refR2
		delta = O.interactions[Body2.id, Body1.id].geom.penetrationDepth
		if makeCDM:
			yieldCont = lawHMCDM.ratioYieldingContacts()
			contactRadius = O.interactions[Body2.id, Body1.id].phys.R
			Rorg = (Da * Db / (Da + Db))
			afac = (1 - np.sin(alpha)) / np.sin(alpha)
			delta_el = max((delta - (contactRadius - Rorg) * afac), 0)
		else:
			yieldCont = 0
			contactRadius = (Da * Db / (Da + Db))
			delta_el = delta
		Es = young / 2.0 / (1 - poisson**2)
		stress = 2 * Es / np.pi * np.sqrt(delta_el / contactRadius)
		plot.addData(Fz=Fz, overlap=delta, contactRadius=contactRadius, yieldCont=yieldCont, stress=stress, i=O.iter)
	else:
		return


####### set timestep to a fraction of the critical timestep
O.dt = .8 * PWaveTimeStep()
stopAtIter = int(300)
O.run()

O.wait()

## plotting
plot.plots = {'overlap': ('stress',), 'overlap ': ('contactRadius',), ' overlap ': ('Fz',)}
plot.plot()
