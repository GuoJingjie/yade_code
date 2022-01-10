# -*- encoding=utf-8 -*-
# jerome.duriez@inrae.fr

# illustrating the use of LevelSet-shaped bodies and their export to Paraview

# import log
# log.setLevel("LevelSet",4)
# log.setLevel("ShopLS",4) # uncomment this group if you wish to have C++ debug messages of classes used

# I. Example of LevelSet Body definitions
#########################################

# first example of defining a level set unit sphere, as a predefined shape
sph1 = levelSetBody("sphere", radius=1, spacing=0.1)

# alternative definition of the same unit sphere, through direct assignement of distance field (which could be adapted to different cases)
grid = RegularGrid(
        -1.1, 1.1, 23
)  # the regular grid upon which the distance field will be defined. Syntax shortcut to RegularGrid(min=(-1.1,-1.1,-1.1),nGP=(23,23,23),spacing=(1.1+1.1)/(23-1)=0.1)
distField = []  # initial empty list
for xInd in range(grid.nGP[0]):
	field_x = []  # some x-cst data set (for some x-cst plane of gridpoints)
	for yInd in range(grid.nGP[1]):
		field_xy = []  # x and y being cst, z variable
		for zInd in range(grid.nGP[2]):
			field_xy.append(
			        (grid.gridPoint(xInd, yInd, zInd)[0]**2 + grid.gridPoint(xInd, yInd, zInd)[1]**2 + grid.gridPoint(xInd, yInd, zInd)[2]**2)**0.5
			        - 1
			)  # distance function to the unit sphere
		field_x.append(field_xy)
	distField.append(field_x)
sph2 = levelSetBody(grid=grid, distField=distField)

# print('sph1 vs sph2 volume comparison',sph1.shape.volume()/sph2.shape.volume()) # you can check eg here they're the same bodies

# yet another definition: direct assignment of distance field, but using numpy arrays
axis = numpy.linspace(-1.1, 1.1, 23)
X, Y, Z = numpy.meshgrid(axis, axis, axis, indexing='ij')
distField = (X**2 + Y**2 + Z**2)**0.5 - 1
sph3 = levelSetBody(grid=grid, distField=distField.tolist())

# II. Going further: interactions between two (non-spherical) LevelSet Bodies and Paraview export/visualisation
###############################################################################################################

lsb = levelSetBody('superellipsoid', extents=(0.07, 0.09, 0.14), epsilons=(1.1, 0.2), spacing=0.01, nSurfNodes=2502, nodesPath=1, dynamic=False)
O.bodies.append(lsb)
lsb = levelSetBody(
        'superellipsoid', center=(0, 0, 0.25), extents=(0.07, 0.09, 0.14), epsilons=(1.1, 0.2), spacing=0.01, nSurfNodes=2502, nodesPath=1, dynamic=False
)
O.bodies.append(lsb)

O.engines = [
        ForceResetter(),
        InsertionSortCollider([Bo1_LevelSet_Aabb()],
                              verletDist=0)  # absence of Sphere bodies prevents an automatic determination of verletDist. We use here 0 (non-optimal)
        ,
        InteractionLoop(
                [Ig2_LevelSet_LevelSet_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack(sphericalBodies=False)
                ]  #ScGeom and that Law2 also apply to LevelSet shapes, but sphericalBodies=False is required in the general case
        ),
        NewtonIntegrator(),
        VTKRecorder(recorders=["lsBodies"], iterPeriod=1)  # multiblockLS=True # if you prefer a unique output file
]
O.bodies[0].state.angVel = (1, 0, 0)
O.dt = 1
O.run(3, True)

# bodies can be displayed after simulation in Paraview / ipython calling e.g.
# pvVisu(idBodies = [0,1],itSnap = range(3))
# where pvVisu() is defined in pvVisu.py, same folder

O.save('lsBodyExample.yade.gz')  # saving the simulation: compressed .yade.gz is advised...
