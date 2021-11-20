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
	frCrea.createPeriodicCell()
	frCrea.createRugosity()
	frCrea.createParticles()

	# Cosmetics
	#renderer = yade.qt._GLViewer.Renderer()
	#renderer.bgColor = Vector3(1.0, 1.0, 1.0)
	#renderer.cellColor = Vector3(1.0, 1.0, 1.0)
	#renderer.ghosts = False

	f = open('.ids', 'w')
	f.write(str(ids))
	f.close()
