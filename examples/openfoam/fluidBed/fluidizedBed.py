
from __future__ import print_function
import sys
from yadeimport import *
from yade.utils import *

initMPI() #Initialize the mpi environment, always required.
fluidCoupling = yade.FoamCoupling();
fluidCoupling.getRank();

class simulation():

  def __init__(self):

    young = 1.2e5; density = 1000;

    O.materials.append(FrictMat(young=young,poisson=0.3,frictionAngle=0.1,density=density,label='spheremat'))
    O.materials.append(FrictMat(young=young,poisson=0.5,frictionAngle=0,density=0,label='wallmat'))

    xmax = 0.044; ymax = 0.120; zmax = 0.01; ybed = 0.03; 
    #sphereIDs = []
    mn, mx= Vector3(0,0,0), Vector3(xmax, ybed, zmax)
    pred = pack.inAlignedBox(mn,mx)
    O.bodies.append(pack.regularHexa(pred, radius=0.60e-3, gap=0, material='spheremat'))
    sphereIDs = [b.id for b in O.bodies if type(b.shape)==Sphere]

    fluidCoupling.setNumParticles(len(sphereIDs))
    fluidCoupling.setIdList(sphereIDs)
    fluidCoupling.isGaussianInterp=True;  #use pimpleFoamYade for gaussianInterp
    
    #walls. 
    walls = aabbWalls([Vector3(0,0,0), Vector3(xmax, ymax, zmax)], thickness=0, material = 'wallmat')
    O.bodies.append(walls)
    newton=NewtonIntegrator(damping=0.0, gravity = (0.0 ,0.0, 0.0)) 

    O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.5, label = "ts"),
        fluidCoupling, #to be called after timestepper
        PyRunner(command='sim.printMessage()', iterPeriod= 1000, label='outputMessage'),
	newton,
        VTKRecorder(fileName='yadep/3d-vtk-',recorders=['spheres'],iterPeriod=1000)
    ]

  def printMessage(self):
     print("********************************YADE-ITER = " + str(O.iter) +" **********************************")



  def irun(self,num):
      O.run(num,1)


if __name__=="__main__":
    sim = simulation()
    sim.irun(1000000)
    fluidCoupling.killMPI()

import builtins
builtins.sim=sim
