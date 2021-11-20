# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
execfile("../common/frameworkCreationUtils.py")

### Useful ids
ids = {}


def frameworkCreation():
	"""Creates the framework
	
	"""

	######################################################################################
	### Framework creation
	######################################################################################

	frCrea.defineMaterials()
	frCrea.createPeriodicCell(pM.l, 4.0 * pM.w, pM.h)
	frCrea.createWalls()
	#	frCrea.createGround()
	frCrea.createRugosity()
	frCrea.createParticles()

	f = open('.ids', 'w')
	f.write(str(ids))
	f.close()
