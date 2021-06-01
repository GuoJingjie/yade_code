# -*- encoding=utf-8 -*-
#
# This is the same simulation as examples/test/performance/checkPerf.py, but does not perform a benchmark.
# It can be used to see the --stdbenchmark simulation in the 3D view.

from builtins import range
from yade import pack,export,geom,timing,bodiesHandling,math
from yade.math import Real
import time,numpy
if(math == sys.modules['math']): raise RuntimeError("To avoid precision problems avoid the native python math and use yade.math")

radRAD      = Real('30.77')
nbIter      = 7000
initIter    = 1 # number of initial iterations excluded from performance check (e.g. Collider is initialised in first step therefore exclude first step from performance check)
tStartAll   = time.time()
O.reset()
tc          = Real('0.001')
en          = Real('.003')
es          = Real('.003')
frictionAngle=math.radiansHP1(35)
density     = Real(2300)
defMat=O.materials.append(ViscElMat(density=density,frictionAngle=frictionAngle,tc=tc,en=en,et=es))

O.dt        = Real('.1')*tc  # time step
rad         = Real('0.5')    # particle radius
tolerance   = Real('0.0001')

SpheresID=[]
SpheresID+=O.bodies.append(pack.regularHexa(pack.inSphere((Vector3(0,0,0)),rad),radius=rad/radRAD,gap=(rad/radRAD)/2,material=defMat))

geometryParameters = bodiesHandling.spheresPackDimensions(SpheresID)
print(len(SpheresID))

floorId=[]
floorId+=O.bodies.append(geom.facetBox(geometryParameters['center'],geometryParameters['extends']/2*Real('1.05'),material=defMat)) # Floor

#Calculate the mass of spheres
sphMass = getSpheresVolume()*density

iterPyRunnerPeriod=100
# Create engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		[Ip2_ViscElMat_ViscElMat_ViscElPhys()],
		[Law2_ScGeom_ViscElPhys_Basic()],
	),
	NewtonIntegrator(damping=0,gravity=[0,Real('-9.81'),0]),
	PyRunner(dead=False,command='colorEnergyKinetic.go()'    ,initRun=True,iterPeriod=iterPyRunnerPeriod,firstIterRun=iterPyRunnerPeriod-1),
	PyRunner(dead=True ,command='colorEnergyNormalForce.go()',initRun=True,iterPeriod=iterPyRunnerPeriod,firstIterRun=iterPyRunnerPeriod-1),
]

print("number of bodies %d"%len(O.bodies))

if(opts.nogui==False):
    yade.qt.View()                # create a GUI view
    #Gl1_Sphere.stripes = True    # mark spheres with stripes
    rr = yade.qt.Renderer()       # get instance of the renderer
    #rr.intrAllWire = True        # draw wires
    #rr.intrPhys = True           # draw the normal forces between the spheres.
    rr.bgColor=Vector3(1,1,1)     # white background
    # If one needs to exactly copy camera position and settings between two different yade sessions, the following commands can be used:
    # https://yade-dem.org/doc/yade.qt.html#yade.qt._GLViewer.views
    v=yade.qt.views()[0]
    v.lookAt, v.viewDir, v.eyePosition, v.upVector = ((1.4550,1.4473,4.8060),(-0.2770,-0.2842,-0.9178),(1.7320,1.7315,5.7239),(-0.08516,0.9587,-0.2711))
    print("Press o/p keys to adjust camera's focal length. (press h for help)")

# FIXME: we should have some nice general way for this.
# Spheres coloring
def colorScale(maxVal,val):
	import math
	gray=val/maxVal if maxVal!=0 else 1.0
	return [   math.sqrt(gray), gray*gray, math.sin(gray*2*math.pi+0.3)  ]

class SimpleColor:
	def __init__(self,exclude,function,O,maxV):
		self.exclude  = exclude
		self.function = function
		self.O        = O
		self.maxV     = maxV # use this to have same color scale on whole simulation.
	def go(self):
		maxValue = self.maxV if self.maxV else 0
		if(not self.maxV):
			for body in self.O.bodies:
				if(body.id not in self.exclude):
					maxValue = max(maxValue,self.function(body,self.O))
		for body in self.O.bodies:
			if(body.id not in self.exclude):
				color = colorScale(maxVal = maxValue, val = self.function(body,self.O))
				body.shape.color = color
		# print("maxValue = ",maxValue) # find the optimal value for self.maxV

def kineticEnergy(body,O):
	return (body.state.mass*body.state.vel.norm()**2)/2

def totalForce(body,O):
	return O.forces.f(body.id).norm()

def totalMoment(body,O):
	return O.forces.m(body.id).norm()

excludeColoring        = floorId
maxEKinValue           = 0.2036816691496742
colorEnergyKinetic     = SimpleColor(excludeColoring,kineticEnergy   ,O,maxEKinValue)
colorTotalForce        = SimpleColor(excludeColoring,totalForce      ,O,None)
colorMoment            = SimpleColor(excludeColoring,totalMoment     ,O,None)

# color spheres on start
colorEnergyKinetic.go()

