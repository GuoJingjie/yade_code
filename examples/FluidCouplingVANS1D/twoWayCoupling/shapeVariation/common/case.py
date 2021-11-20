#########################################################################################################################################################################
# Author : Remi Monthiller, remi.monthiller@etu.enseeiht.fr
# Adapted from the code of Raphael Maurin, raphael.maurin@imft.fr
# 30/10/2018
#
# Incline plane simulations
#
#########################################################################################################################################################################

# import lib
import os
import matplotlib.pyplot as plt

# import params
execfile('params.py')

# Simulation
try:
	datas = os.listdir("data")
except:
	os.mkdir("data")
	datas = False
if datas:
	for i in range(len(datas)):
		datas[i] = float(datas[i].split(".yade")[0])
	datas.sort()
	# import PyRunners
	execfile('framework.py')
	execfile('../common/simulationPyRunners.py')
	O.load("data/" + str(datas[-1]) + ".yade")
	execfile('params.py')
	if pN.enable_new_engines:
		execfile('../common/simulationDefinition.py')
		for e in O.engines:
			del e
		sim.engineCreation()
		#sim.resetEngine()
		#sim.init()
		#O.resetTime()
		O.saveTmp()
	#O.run()
else:
	# import simulation
	execfile('framework.py')
	execfile('../common/simulationDefinition.py')
	sim.simulation()
