# -*- encoding=utf-8 -*-
# jerome.duriez@inrae.fr

# illustrating the use of LevelSet-shaped bodies

# import log
# log.setLevel("LevelSet",4)
# log.setLevel("ShopLS",4) # uncomment this group if you wish to have C++ debug messages of classes used

# I. Example of LevelSet Body definitions
#########################################

# first example of defining a level set unit sphere, as a predefined shape
sph1 = levelSetBody("sphere",radius=1,spacing=0.1)

# alternative definition of the same unit sphere, through direct assignement of distance field (which could be adapted to different cases)
grid = RegularGrid(-1.1,1.1,23) # the regular grid upon which the distance field will be defined. Syntax shortcut to RegularGrid(min=(-1.1,-1.1,-1.1),nGP=(23,23,23),spacing=(1.1+1.1)/(23-1)=0.1)
distField = [] # initial empty list
for xInd in range(grid.nGP[0]):
    field_x = [] # some x-cst data set (for some x-cst plane of gridpoints)
    for yInd in range(grid.nGP[1]):
        field_xy=[] # x and y being cst, z variable
        for zInd in range(grid.nGP[2]):
            field_xy.append((grid.gridPoint(xInd,yInd,zInd)[0]**2 + grid.gridPoint(xInd,yInd,zInd)[1]**2 + grid.gridPoint(xInd,yInd,zInd)[2]**2)**0.5 - 1) # distance function to the unit sphere
        field_x.append(field_xy)
    distField.append(field_x)
sph2 = levelSetBody(grid=grid,distField=distField)

# yet another definition: direct assignment of distance field, but using numpy arrays
axis = numpy.linspace(-1.1,1.1,23)
X,Y,Z = numpy.meshgrid(axis,axis,axis,indexing='ij')
distField = (X**2+Y**2+Z**2)**0.5 - 1
sph3 = levelSetBody(grid=grid,distField=distField.tolist())

# print('sph1 vs sph2 volume comparison',sph1.shape.volume()/sph2.shape.volume()) # you can check eg here they're the same bodies


# II. Going further in subsequent YADE revisions
