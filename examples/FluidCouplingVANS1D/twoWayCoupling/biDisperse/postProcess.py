###############################################################################################################################################
# Authors: 	Remi Chassagne, remi.chassagne@univ-grenoble-alpes.fr
#
# 03/09/2021
#
# Script to post-process the data from script-segregation.py
# To use it: launch ipython and then type the command execfile('postProcess.py')
#
# To compare to a reference run: keep the figures open, put name = 'Ref' in postProcess.py and type again execfile('postProcess.py')
#
################################################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.ion()
name = ''
##############################
#Load time data from hdf5 file
##############################
loadFile = h5py.File("data" + name + ".hdf5")
#Keys of files are strings of int corresponding to time steps
timeSteps = list(loadFile.keys())
timeSteps = [int(timeSteps[i]) for i in range(len(timeSteps))]
timeSteps.sort()  #sort the array

#Define physical quantities
time = np.linspace(0.5, len(timeSteps) * 0.5, len(timeSteps))
grp = loadFile['0']
ndimz = len(grp['phiPart'].value)
fluidHeight = 0.117
zAxis = np.linspace(0, fluidHeight, ndimz)
diameterPart1 = 0.006
diameterPart2 = 0.003

########################
#Compute averaged values
########################
vxPart = np.zeros(ndimz)
vxPart1 = np.zeros(ndimz)
vxPart2 = np.zeros(ndimz)
phiPart = np.zeros(ndimz)
phiPart1 = np.zeros(ndimz)
phiPart2 = np.zeros(ndimz)
vxFluid = np.zeros(ndimz)
epsFluid = np.zeros(ndimz)
massCenter = np.zeros(len(time))
qs = np.zeros(len(time))
for timeStep in timeSteps:
	grp = loadFile[str(timeStep)]
	#Granular quantities
	phiPart += grp['phiPart'].value
	phiPart1 += grp['phiPart1'].value
	phiPart2 += grp['phiPart2'].value
	vxPart += grp['phiPart'].value * grp['vxPart'].value
	vxPart1 += grp['phiPart1'].value * grp['vxPart1'].value
	vxPart2 += grp['phiPart2'].value * grp['vxPart2'].value
	#Fluid quantities
	epsFluid += (1 - grp['phiPart'].value)
	vxFluid += (1 - grp['phiPart'].value) * grp['vxFluid'].value[1:]
	#Transport
	qs[timeStep] = grp['qsMean'].value
	#Small particle mass center position
	massCenter[timeStep] = np.sum(zAxis * grp['phiPart2'].value) / np.sum(grp['phiPart2'].value)

loadFile.close()

nonZeroIndices = np.where(phiPart > 0)
vxPart[nonZeroIndices] /= phiPart[nonZeroIndices]
phiPart /= len(timeSteps)

nonZeroIndices = np.where(phiPart1 > 0)
vxPart1[nonZeroIndices] /= phiPart1[nonZeroIndices]
phiPart1 /= len(timeSteps)

nonZeroIndices = np.where(phiPart2 > 0)
vxPart2[nonZeroIndices] /= phiPart2[nonZeroIndices]
phiPart2 /= len(timeSteps)

vxFluid /= epsFluid
epsFluid /= len(timeSteps)

##############
# Plot results
##############
plt.figure(1)
plt.plot(vxPart, zAxis / diameterPart1, label=r'$\left<v_x^p\right>^p$ ' + name)
plt.plot(vxFluid, zAxis / diameterPart1, label=r'$\left<v_x^f\right>^f$ ' + name)
plt.ylabel(r'$\frac{z}{d_1}$', rotation=True, horizontalalignment='right', fontsize=20)
plt.grid()
plt.legend()

plt.figure(2)
plt.plot(phiPart, zAxis / diameterPart1, label=r'$\phi$ ' + name)
plt.plot(epsFluid, zAxis / diameterPart1, label=r'$\epsilon$ ' + name)
plt.ylabel(r'$\frac{z}{d_1}$', rotation=True, horizontalalignment='right', fontsize=20)
plt.grid()
plt.legend()

plt.figure(3)
plt.plot(time, qs, label=name)
plt.xlabel(r'$t$ (s)', fontsize=25)
plt.ylabel(r'$q_s$', rotation=True, horizontalalignment='right', fontsize=25)
plt.grid()
plt.legend()

plt.figure(4)
plt.semilogx(time, massCenter / diameterPart1, label=name)
plt.xlabel(r'$t$ (s)', fontsize=25)
plt.ylabel(r'$\frac{z_c}{d_1}$', rotation=True, horizontalalignment='right', fontsize=25)
plt.grid()
plt.legend()
