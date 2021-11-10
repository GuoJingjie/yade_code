# -*- encoding=utf-8 -*-
# jerome.duriez@inrae.fr

from paraview.simple import * # import Paraview stuff (already done by Paraview itself when using that software, mandatory outside Paraview)

def pvVisu(prefix = '',itSnap = [0],thres = False,multiblockLS = False,idBodies = [0]):
    '''Displays, in Paraview typically, LevelSet-shaped bodies as exported by YADE's VTKRecorder. NB Paraview hint: use os.getcwd() to check Paraview's Python (2 in Paraview 5.4.1) interface working directory, and os.chdir() to change it.
    :param string prefix: the VTKRecoder.fileName path prefix to the files to read
    :param list itSnap: the snapshot iteration numbers
    :param bool thres: if True, bodies are sharply rendered from their inner gridpoints only. When False, rendering involves the zero-contour of the distance function
    :param bool multiblockLS: previously used value of VTKRecorder.multiblockLS
    :param list idBodies: the ids of LevelSet-shaped bodies to look at, necessary iff !multiblockLS
    :returns: zero-distance contours if not thres, or inner (thresholded) volume otherwise, in addition to displaying them. The returned item is a list of those if !multiblockLS. Display features can later on be controlled in Paraview, modifying eg display.DiffuseColor = RGB list with display = Show(a returned object from that function).
    '''
    #NB: some automatic handle of idBodies/itSnap could be pursued through
    # strF = subprocess.run("ls lsBody*.vts",shell=True,stdout=subprocess.PIPE)
    # bytcop = strF.stdout.read() # binary object to .decode()

    # Importing all LevelSet bodies data
    ####################################
    if multiblockLS: # either into one unique XMLMultiBlockDataReader
        lsFiles = [] # will include all (time passing) vtm files
        for it in itSnap:
            lsFiles.append(prefix+'lsBodies.'+str(it)+'.vtm')
        multiReader = XMLMultiBlockDataReader(FileName=lsFiles) # one reader item, for all bodies, and all time steps
    else: # or into a list of XMLStructuredGridReader
        allReaders = [] # a list of all readers, with each reader having all time steps for one body
        for bod in idBodies:
            lsBodFiles = [] # will include all (time passing) files for that body
            for it in itSnap:
                lsBodFiles.append(prefix+'lsBody'+str(bod)+'.'+str(it)+'.vts')
            # print 'Reading',lsBodFiles
            allReaders.append( XMLStructuredGridReader(FileName=lsBodFiles) )

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()

    # Extracting inner gridpoints / outer surface from the reader(s)
    ################################################################
    if multiblockLS:
        retItem = handleReader(multiReader,thres,multiblockLS)
    else:
        retItem = [] # returned item will be a list
        for bod in range(0,len(idBodies)):
            retItem.append(handleReader(allReaders[bod],thres,multiblockLS))

    ResetCamera() # useful for an appropriate zoom, probably calls Render() in background as well
    return retItem

def handleReader(reader,thres,multiblockLS):
    '''An auxiliary function for pvVisu. Can stay in the background'''
    if thres:
        threshold = Threshold(Input=reader)
        threshold.Scalars = ['POINTS', 'distField'] # threshold is on distance value, which is stored as distField
        (minDist,maxDist) = reader.PointData['distField'].GetRange()
        threshold.ThresholdRange = [minDist, 0]
        disp = Show(threshold)
        retObj = threshold
    else:
        contour = Contour(Input=reader,ContourBy='distField', Isosurfaces = [0.0])
        disp = Show(contour)
        disp.Representation = 'Surface'
        retObj = contour
    if multiblockLS: ColorBy(disp, ('FIELD', 'vtkBlockColors'))
    return retObj
