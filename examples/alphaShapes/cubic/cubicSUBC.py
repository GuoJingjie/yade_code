# -*- encoding=utf-8 -*-
################################################################################
# Load the saved confined cubic aggregate using a SUBC (Static Uniform 
# Boundary Condition) applied as a Cauchy stress tensor on the aggregate
# power diagram described by the alpha shape surface reconstruction algorithm
################################################################################
import os, sys, gts
from scipy.interpolate import griddata
from yade import pack, export
import numpy as np
#######################################
###       FUNCTION DEFINITIONS      ###
#######################################
def averageDefGrad(xq,Aq,vol0): # Can be used for Deformation/Velocity Gradient Tensors
  defGrad = Matrix3.Zero
  for i, el in enumerate(Aq):
    defGrad = defGrad + Vector3.outer(xq[i],Aq[i]) / vol0
  return defGrad

def triStressMicro(stressPt,stabThresh,triThresh,creepThresh):
  global stressM, currF
  deltaSigma = stressPt - stressM
  sigmaDot = 5.e3 * deltaSigma.normalized() / deltaSigma.normalized().maxAbsCoeff()
  count = int(deltaSigma.norm() / sigmaDot.norm())
  print('Incremental Loading Steps: ', count)
  stage = 1 # Loading stages 1) Increment Stress 2) Creep
  t0 = O.time
  while 1:
    for j in range(20):
      TW.triangulate()
      caps=TW.getAlphaCaps(alpha,alphaShrinked,True)
      velArr = [O.bodies[cap[0]].state.vel for cap in caps]
      areaArr = [cap[1] for cap in caps]
      lbar = averageDefGrad(velArr,areaArr,TW.alphaCapsVol)
      currF = currF + (lbar * currF) * (O.time - t0)
      currStrain = currF - Matrix3.Identity
      t0 = O.time
      meanBVel = 0.5 * lbar.spectralDecomposition()[1].maxAbsCoeff()
      integrator.damping=max(0.2,min(0.999,(meanBVel/boundVel)**4))
      stressM = (stressM + sigmaDot) if stage == 1 else stressPt
      count = (count - 1) if stage == 1 else 0
      if count <= 0 and stage == 1:
        stage = 2
      TW.applyAlphaForces(stressM,alpha,alphaShrinked,True)
      O.run(10, True) 
    TW.triangulate()
    caps=TW.getAlphaCaps(alpha,alphaShrinked,True)
    unb = unbalancedForce()
    meanStress = TW.calcAlphaStress(alpha,alphaShrinked,True)
    print('unb:',unb,' Boundary Velocity:',meanBVel)
    print('Stress Tensor\n',np.asarray(meanStress))
    print('Velocity Gradient Tensor\n',np.asarray(lbar))
    sflag = 0
    for i in range(0,3):
      for j in range(0,3):
        triDiff = abs(meanStress[i][j]-stressPt[i][j]) / max( abs(meanStress[i][j]) , abs(meanStress.trace()/3.) )
        if triDiff > triThresh:
          sflag += 1
    if unb < stabThresh and stage == 2 and sflag == 0 and meanBVel < abs(creepThresh * boundVel):
        break

def cart2sph(cartX, cartY, cartZ):
  point = Vector3(cartX,cartY,cartZ)
  sphR = np.sqrt(point.dot(point))
  sphTheta = np.arctan2(point[1],point[0])
  sphPhi = np.arccos(point[2]/sphR)
  return sphTheta, sphPhi, sphR

def sph2cart(theta, phi, r):
  x = r * np.cos(theta) * np.sin(phi)
  y = r * np.sin(theta) * np.sin(phi)
  z = r * np.cos(phi)
  return x, y, z
  
def ancillaryGrid(alpha,alphaShrinked,isFix,s,newBound):
  TW=TesselationWrapper()
  TW.triangulate()
  segments=TW.getAlphaGraph(alpha,alphaShrinked,isFix)
  caps=TW.getAlphaCaps(alpha,alphaShrinked,isFix)
  # aaBB Center and Size
  center, length = np.mean(aabbExtrema(),axis=0), np.mean(aabbDim(),axis=0)
  # Use alphaGraph and alphaCaps to build a boundary grid 
  polygons = np.asarray([cart2sph(*vertex-center) for vertex in segments])
  polygons = np.vstack([tuple(row) for row in polygons])
  boundGrid1 = np.vstack([polygons])
  # Mirror points to create a continuous surface in spherical coords
  boundGrid = np.vstack([boundGrid1+[revTh,revPh,0.]
                            for revTh in [-2*np.pi,0,2*np.pi] 
                            for revPh in [-np.pi,0,np.pi]])
  # Convert unit cube to spherical coordinates and interpolate surface
  gridGts = np.asarray([cart2sph(*vert) for vert in newBound])
  gridGts[:,2] = griddata(boundGrid[:,0:2], boundGrid[:,2], gridGts[:,0:2])
  for index, vertex in enumerate(s.vertices()):
      vertex.x, vertex.y, vertex.z = (sph2cart(*gridGts[index]) + center)
  return [[Vector3(np.mean([vertex.coords() for vertex in face.vertices()],axis=0)) 
                                                  for face in s.faces()],
          [Vector3(face.normal())/2. for face in s.faces()],
          s.volume()]

def sepath(alpha,alphaShrinked,isFix):
  global currF
  print('********** Saving State Variables ***************')
  newBound = [ogPt * currF for ogPt in ogBound]
  bSurf.append(ancillaryGrid(alpha,alphaShrinked,isFix,unitCube,newBound))
  print('currF\n',np.asarray(currF))
  Fbar.append(averageDefGrad(bSurf[-1][0],bSurf[0][1],bSurf[0][2]))
  print('Fbar\n',np.asarray(Fbar[-1]))
  # Note that currF and Fbar are both evaluations of the deformation gradient
  # currF tracks from velocity gradient, while Fbar obtains it after reinterpolating
  # the surface grid to the latest power diagram
  Sbar.append(TW.calcAlphaStress(alpha,alphaShrinked,True))
  FbarRef, FbarCurr, SbarRef, SbarCurr = Fbar[0], Fbar[-1], Sbar[0], Sbar[-1]
  eGreen = 0.5 * ( FbarCurr * FbarCurr.transpose() - Matrix3.Identity )
  hydroStress = SbarCurr.trace() / 3.e6
  dhydroStress = (SbarCurr.trace()-SbarRef.trace()) / 3.e6
  volStrain = FbarCurr.determinant() -1.
  pStress = SbarCurr.spectralDecomposition()[1] / 1.e6 - SbarRef.spectralDecomposition()[1] / 1.e6
  deviatorStress = pStress.maxCoeff() - pStress.minCoeff()
  axStrain = eGreen[0][0]
  pStrain = eGreen.spectralDecomposition()[1]
  deviatorStrain = pStrain.maxCoeff() - pStrain.minCoeff()
  with open('resultsSUBC/cauchygreen.txt', 'a') as f:
    f.write('%f %f %f %f %f\n'% (hydroStress,deviatorStress,volStrain,axStrain,deviatorStrain) )

def vtkOut(suffix,alpha,alphaShrinked,isFix):
  vtkExporter = yade.export.VTKExporter('resultsSUBC/spheres'+suffix)
  vtkExporter.exportSpheres()

  TW.triangulate()
  segments=TW.getAlphaGraph(alpha,alphaShrinked,True)
  caps=TW.getAlphaCaps(alpha,alphaShrinked,True)

  f=open('resultsSUBC/LVdiagram'+suffix+'.vtk','w')
  f.write("# vtk DataFile Version 2.0\nGenerated from YADE\n")
  f.write("ASCII\nDATASET POLYDATA\n")
  f.write("POINTS %d float \n" %(len(segments)))
  [f.write("%f %f %f\n" %(point[0],point[1],point[2])) for point in segments]
  f.write("LINES %d %d \n" %(len(segments)/2, 3*len(segments)/2))
  [f.write("2 %d %d\n" %(2*k,2*k+1)) for k in range(int(len(segments)/2))]
  f.close()

  f=open('resultsSUBC/gtsSurface'+suffix+'.vtk','w')
  unitCube.write_vtk(f)
  f.close()
############################################
###     LOAD REFERENCE CONFIGURATION     ###
############################################
if not(os.path.isdir('resultsSUBC')):
  os.mkdir('resultsSUBC')

O.load('confinedState.yade.gz')
O.resetTime()
############################################
###           FIND ALPHA GRAPH           ###
############################################
minRad = min(b.shape.radius for b in O.bodies)
meanDia = 2. * sum([b.shape.radius for b in O.bodies]) / len(O.bodies)
alpha = (22. * minRad) ** 2
alphaShrinked = (21. * minRad) ** 2
isFix = True

TW=TesselationWrapper()
TW.triangulate()
############################################
###         INITIALIZE VARIABLES         ###
############################################
stressM, rate_v = Matrix3.Identity * (-5.0e6), Matrix3.Zero
TW.applyAlphaForces(stressM,alpha,alphaShrinked,True)
O.step()
Sbar = [TW.calcAlphaStress(alpha,alphaShrinked,True)]
print("Initial Stress Tensor \n",np.asarray(Sbar[0]))
currF = Matrix3.Identity
Fbar = [currF]
with open('resultsSUBC/cauchygreen.txt', 'w') as f:
  f.write('pStress qStress vStrain aStrain dStrain')
###########################################
###      ANCILLARY GRID AT REF. STATE   ###
###########################################
# Build a unit gts cube and subdivide it 'subs' times
subs = 5
unitCube = gts.cube()
for i in range(subs):
  unitCube.tessellate()
  ogBound = [Vector3(vertex.coords()) for vertex in unitCube.vertices()]
  for i, vertex in enumerate(unitCube.vertices()):
    vertex.x, vertex.y, vertex.z = ( ogBound[i] / ogBound[i].maxAbsCoeff() )
# Use the unit gts cube to interpolate ancillary surface grid
ogBound = [Vector3(vertex.coords()) for vertex in unitCube.vertices()]
bSurf = []
sepath(alpha,alphaShrinked,True)
###########################################
###      TRIAXIAL STRESS/STRAIN GOALS   ###
###########################################
# Inertial Number
iNum = 7.9e-5
# Rate of Deformation Tensor Increment
boundVel = .45 * iNum * sqrt(5.e6 /O.bodies[0].mat.density) / meanDia
print('Quasistatic Strain Rate', boundVel)
###########################################
###             LOADING PATH            ###
###########################################
stressInc = -1.0e6
stressPath = [stressM]
# Load by an isotropic stress increment
stressPath.append(stressPath[-1] + Matrix3.Identity * stressInc)
# # Load by a uniaxial stress increment
# stressPath.append(stressPath[-1] + Matrix3(Vector3(stressInc,0.,0.),Vector3.Zero, Vector3.Zero))
# # Load by a deviatoric stress increment
# stressPath.append(stressPath[-1] + Matrix3(Vector3(0.,0.,0.5*stressInc),
#                                     Vector3.Zero,
#                                     Vector3(0.5*stressInc,0.,0.)))
###########################################
###            BOUNDARY CONDITION       ###
###########################################
print('APPLY SUBC')
for ptNum,stressPt in enumerate(stressPath[1:]):
  print('Load to \n', np.asarray(stressPt) / 1.e6,' MPa')
  stepper.timestepSafetyCoefficient = 0.6
  unbTol, triTol, creepTol = 0.01,0.01,0.1
  loader.iterPeriod = 0 # Loading is done in function triStressMicro()
  triStressMicro(stressPt,unbTol,triTol,creepTol)
  sepath(alpha,alphaShrinked,True)
######################################
###             VTK EXPORTS        ###
######################################
vtkOut("SUBC",alpha,alphaShrinked,isFix)
