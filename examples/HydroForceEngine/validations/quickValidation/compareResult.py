#########################################################################################################################################################################
# Author: Raphael Maurin, raphael.maurin@imft.fr
# 25/08/2021
#
# Post processing script to extract and compare the results from dataRef & data, coupled fluid-DEM simulations
#
############################################################################################################################################################################


import numpy as np
import os
from matplotlib.pyplot import *
import matplotlib.gridspec as gridspec

scriptPath = os.path.abspath(os.path.dirname(sys.argv[-1])) #Path where the script is stored


if os.path.exists(scriptPath +'/data/')==False:	#If the data folder does not exist, no data to extract, exit. 
	print('\n There is no data to extract in this folder ! Please first run sedimentTransportExample_1DRANSCoupling.py !\n')
	exit()
else:	#Else, extract the first file in order to get the size of the vectors, ndimz
	execfile(scriptPath +'/data/0.py')
	ndimz = len(phiPartPY)	#Size of the vectors, mesh parameter. 
if os.path.exists(scriptPath +'/dataRef/')==False:	#If the data folder does not exist, no data to extract, exit. 
	print('\n There is no data to extract in this folder ! Please first run sedimentTransportExample_1DRANSCoupling.py !\n')
	exit()

#Initilization of the variable to extract and plot
[qs,time,phiPart,vxPart,vxFluid] = [[],[],np.zeros(ndimz),np.zeros(ndimz),np.zeros(ndimz+1)]
[qsRef,phiPartRef,vxPartRef,vxFluidRef] = [[],np.zeros(ndimz),np.zeros(ndimz),np.zeros(ndimz+1)]
########
# LOOP over time to EXTRACT the data
########
nFiles = 10
for i in range(nFiles):
	#Data folder
	nameFile = scriptPath +'/data/' + str(i)+'.py'	#Name of the file at the considered time step
	#Extract the data from at the time (file) considered. Assign the vectors to qsMean, phiPartPY, vxPartPY, vxFluidPY, and zAxis see the structure of a saved file X.py. 
	execfile(nameFile)
	time+=[i*0.1]			#Recreate the time scale
	qs+=[qsMean]			#Extract the sediment transport rate as a function of time
	phiPart+=phiPartPY		#Average over time the solid volume fraction vertical profile
	vxPart+=vxPartPY*phiPartPY	#Average over time the solid velocity profile. The spatial averaging is weighted by the particle volume so that it is necessary to re-multiply the averaging by the solid volume fraction before re-normalizing. See the paper Maurin et al (2015) in the references for more details on the averaging. 
	vxFluid+=vxFluidPY		#Average over time the fluid velocity vertical profile

	#Ref folder
	nameFileRef = scriptPath +'/dataRef/' + str(i)+'.py'	#Name of the file at the considered time step
	#Extract the data from at the time (file) considered. Assign the vectors to qsMean, phiPartPY, vxPartPY, vxFluidPY, and zAxis see the structure of a saved file X.py. 
	execfile(nameFile)
	qsRef+=[qsMean]			#Extract the sediment transport rate as a function of time
	phiPartRef+=phiPartPY		#Average over time the solid volume fraction vertical profile
	vxPartRef+=vxPartPY*phiPartPY	#Average over time the solid velocity profile. The spatial averaging is weighted by the particle volume so that it is necessary to re-multiply the averaging by the solid volume fraction before re-normalizing. See the paper Maurin et al (2015) in the references for more details on the averaging. 
	vxFluidRef+=vxFluidPY		#Average over time the fluid velocity vertical profile

#Normalization
toto = np.where(phiPart>0)
vxPart[toto]/=phiPart[toto]	#Avoid division by zero
phiPart/=float(nFiles)
vxFluid/=float(nFiles)
#Normalization
toto = np.where(phiPartRef>0)
vxPartRef[toto]/=phiPartRef[toto]	#Avoid division by zero
phiPartRef/=float(nFiles)
vxFluidRef/=float(nFiles)

# Compare the results: 
if sum(vxPart - vxPartRef)==0. and sum(phiPart - phiPartRef)==0. and sum(vxFluid - vxFluidRef)==0.:
	print('\n The two simulations gives exactly the same results !\n VALIDATED!\n')
else: 
	[errVp,errVf,errPhi] = [0.,0.,0.]
	n=0
	for vp in vxPart: 
		errVp+=sqrt(pow(vxPart[n] - vxPartRef[n],2.))
		errVf+=sqrt(pow(phiPart[n] - phiPartRef[n],2.))
		errPhi+=sqrt(pow(vxFluid[n] - vxFluidRef[n],2.))
	errVp/=vxPartRef.mean()
	errVf/=vxFluidRef.mean()
	errPhi/=phiPartRef.mean()
	print('\nThe averaged particle velocity profile is different from {0}% between the two configurations'.format(errVp))
	print('\nThe averaged solid volume fraction profile is different from {0}% between the two configurations'.format(errPhi))
	print('\nThe averaged fluid velocity profile is different from {0}% between the two configurations\n'.format(errVf))


######
# PLOT
if True:
	# sediment transport rate as a function of time
	fig1 = figure()		#Create a figure
	gs = gridspec.GridSpec(1,1)	#Split the figure in 1 row and 1 column
	ax=subplot(gs[0,0])		#Assign the axis
	ax.plot(time,qsRef,'ok')		#Plot the sediment transport rate as a function of time
	ax.plot(time,qs,'+r')		#Plot the sediment transport rate as a function of time
	ax.grid()			#Put a grid on the figure
	ax.set_xlabel('t (s)')		#Assign the x axis label
	ax.set_ylabel(r'$Q_s^*$',rotation='horizontal')#Assign the y axis label

	# Solid and fluid velocity depth profiles, solid volume fraction and sediment transport rate density depth profiles
	fig2 = figure()
	gs = gridspec.GridSpec(1,3)	#Split the figure in 1 row and 3 column
	ax1=subplot(gs[0,0])		#ax1 = subfigure corresponding to 1st row 1st column
	ax2=subplot(gs[0,1])		#ax2 = subfigure corresponding to 1st row 2nd column
	ax3=subplot(gs[0,2])		#ax3 = subfigure corresponding to 1st row 3rd column

	ax1.plot(vxPartRef,zAxis,'ok')	#Plot the solid velocity profile on first subfigure
	ax1.plot(vxFluidRef[1:],zAxis,'-ob')	#Plot the fluid velocity profile on first subfigure
	ax2.plot(phiPartRef,zAxis,'ok')	#Plot the solid volume fraction on second subfigure
	ax3.plot(phiPartRef*vxPartRef,zAxis,'ok')	#Plot the sediment transport rate density on the third subfigure

	ax1.plot(vxPart,zAxis,'+r')	#Plot the solid velocity profile on first subfigure
	ax1.plot(vxFluid[1:],zAxis,'+c')	#Plot the fluid velocity profile on first subfigure
	ax2.plot(phiPart,zAxis,'+r')	#Plot the solid volume fraction on second subfigure
	ax3.plot(phiPart*vxPart,zAxis,'+r')	#Plot the sediment transport rate density on the third subfigure

	ax1.set_xlabel(r'$<v^k_x>$ (m/s)')	#Assign the x axis label of the first subfigure
	ax1.set_ylabel(r'$z/d$',rotation='horizontal')	#Assign the y axis label of the first subfigure
	ax2.set_xlabel(r'$\phi$') 	#Assign the x axis label of the second subfigure
	ax3.set_xlabel(r'$q (m/s)$')	#Assign the x axis label of the third subfigure
	ax2.set_yticklabels([])	#Remove the vertical scale for the second subfigure, same as first
	ax3.set_yticklabels([])	#Remove the vertical scale for the third subfigure, same as first
	ax1.grid()	#Put a grid on each subfigure
	ax2.grid()
	ax3.grid()
	fig2.subplots_adjust(left=0.08, bottom=0.11, right=0.98, top=0.98, wspace=.1)	#Adjust the spacing around the full figure (left, bottom, right, top) and between the subfigure (wspace), in percent of the figure size

	show()	#Show the two figures created
