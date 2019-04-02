
""" CAUTION:
Running this script can take very long!
"""
from __future__ import division

from past.utils import old_div
from numpy import arange
from yade import pack
import pylab
# define the section shape as polygon in 2d; repeat first point at the end to close the polygon
poly=((1e-2,5e-2),(5e-2,2e-2),(7e-2,-2e-2),(1e-2,-5e-2),(1e-2,5e-2))
# show us the meridian shape
#pylab.plot(*zip(*poly)); pylab.xlim(xmin=0); pylab.grid(); pylab.title('Meridian of the revolution surface\n(close to continue)'); pylab.gca().set_aspect(aspect='equal',adjustable='box'); pylab.show()
# angles at which we want this polygon to appear
thetas=arange(0,old_div(pi,2),old_div(pi,24))
# create 3d points from the 2d ones, turning the 2d meridian around the +y axis
# for each angle, put the poly a little bit higher (+2e-3*theta);
# this is just to demonstrate that you can do whatever here as long as the resulting
# meridian has the same number of points
#
# There is origin (translation) and orientation arguments, allowing to transform all the 3d points once computed.
#
# without these transformation, it would look a little simpler:
# 	pts=pack.revolutionSurfaceMeridians([[(pt[0],pt[1]+2e-3*theta) for pt in poly] for theta in thetas],thetas
#
pts=pack.revolutionSurfaceMeridians([[(pt[0],pt[1]+1e-2*theta) for pt in poly] for theta in thetas],thetas,origin=Vector3(0,-.05,.1),orientation=Quaternion((1,1,0),old_div(pi,4)))
# connect meridians to make surfaces
# caps will close it at the beginning and the end
# threshold will merge points closer than 1e-4; this is important: we want it to be closed for filling
surf=pack.sweptPolylines2gtsSurface(pts,capStart=True,capEnd=True,threshold=1e-4)
# add the surface as facets to the simulation, to make it visible
O.bodies.append(pack.gtsSurface2Facets(surf,color=(1,0,1)))
# now fill the inGtsSurface predicate constructed form the same surface with sphere packing generated by TriaxialTest
# with given radius and standard deviation (see documentation of pack.randomDensePack)
#
# The memoizeDb will save resulting packing into given file and next time, if you run with the same
# parameters (or parameters that can be scaled to the same one),
# it will load the packing instead of running the triaxial compaction again.
# Try running for the second time to see the speed difference!
memoizeDb='/tmp/gts-triax-packings.sqlite'
sp=SpherePack()
sp=pack.randomDensePack(pack.inGtsSurface(surf),radius=5e-3,rRelFuzz=1e-4,memoizeDb=memoizeDb,returnSpherePack=True)
sp.toSimulation()
# We could also fill the horse with triaxial packing, but have nice approximation, the triaxial would run terribly long,
# since horse discard most volume of its bounding box
# Here, we would use a very crude one, however
if 1:
	import gts
	horse=gts.read(open('horse.coarse.gts')) #; horse.scale(.25,.25,.25)
	O.bodies.append(pack.gtsSurface2Facets(horse))
	sp=pack.randomDensePack(pack.inGtsSurface(horse),radius=5e-3,memoizeDb=memoizeDb,returnSpherePack=True)
	sp.toSimulation()
	horse.translate(.07,0,0)
	O.bodies.append(pack.gtsSurface2Facets(horse))
	# specifying spheresInCell makes the packing periodic, with the given number of spheres, proportions being equal to that of the predicate
	sp=pack.randomDensePack(pack.inGtsSurface(horse),radius=1e-3,spheresInCell=2000,memoizeDb=memoizeDb,returnSpherePack=True)
	sp.toSimulation()
