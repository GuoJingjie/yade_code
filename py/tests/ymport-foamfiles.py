# encoding: utf-8
##########################################################################
#  2022 Tóth János                                                       #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

import unittest
from pathlib import Path
from yade import ymport


# TODO: Add more conditions to the tests to check the validity of the data.
class TestYmportFoamFiles(unittest.TestCase):

	def testLoadBlockMeshDict(self):
		f = str(Path(__file__).parent) + "/ymport-files/blockMeshDict"
		facets = ymport.blockMeshDict(f)
		self.assertTrue(len(facets) == 40)

	def testLoadPolyMesh(self):
		f = str(Path(__file__).parent) + "/ymport-files/polyMesh"
		facets = ymport.polyMesh(f)
		self.assertTrue(len(facets) == 320)
