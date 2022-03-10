#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
from testGuiHelper import TestGUIHelper

yade.log.setLevel("Default", yade.log.ERROR)
scr = TestGUIHelper("LevelSet")

if 'LS_DEM' in yade.config.features:

	guiIterPeriod = 10000

	O.engines = O.engines + [PyRunner(iterPeriod=guiIterPeriod, command='scr.screenshotEngine()')]

	# Create rectangular box from facets
	O.bodies.append(geom.facetBox((0.0, 0.0, 0.0), (5.1, 5.1, 5.1)))  #, wallMask=31

	# Create a complicated level set shape to render
	grid = RegularGrid(-4.99, 5.01, 41)
	fmm = FastMarchingMethod(grid=grid, phiIni=phiIniPy(lambda x, y, z: distApproxRose((x, y, z)), grid), speed=1)
	phiField = fmm.phi()
	O.bodies.append(levelSetBody(grid=grid, distField=phiField, center=(0, 0, 0), nSurfNodes=0))

	O.dt = O.dt * 0.0001
	O.run(guiIterPeriod * scr.getTestNum() + 1)
else:
	print("Skipping testGuiLevelSet, LS_DEM not available in features")
	scr.createEmptyFile("screenshots/testGui_LevelSet_OK_or_Skipped.txt")  # creating verification file so there is no error thrown
	os._exit(0)
