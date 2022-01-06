# encoding: utf-8
# 2021, János Tóth <toth-janos@outlook.com>

import unittest
from pathlib import Path
from yade import ymport


class TestYmportBlockMeshDict(unittest.TestCase):

	def testLoadFile(self):
		f = str(Path(__file__).parent) + "/ymport-files/blockMeshDict"
		faces = ymport.blockMeshDict(f)
		self.assertTrue(len(faces) == 40)
