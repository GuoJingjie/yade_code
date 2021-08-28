import unittest, yade


class TestEnum(unittest.TestCase):

	def setUp(self):
		if self.canTest():  # to test if assignment of enums works properly we need any kind enum.
			self.r = yade.OpenGLRenderingEngine()
		else:
			pass

	def tearDown(self):
		pass  # here finalize the tests, maybe close an opened file, or something else

	def canTest(self):
		return 'QT5' in yade.config.features

	def thisTestsException1(self):
		self.r.blinkHighlight = -100

	def thisTestsException2(self):
		self.r.blinkHighlight = 'incorrect string'

	def testException(self):
		if (not self.canTest()):
			return
		# make sure it's impossible to assign incorrect values.
		print("\n  Here it must throw two exceptions:")
		self.assertRaises(TypeError, self.thisTestsException1)
		self.assertRaises(TypeError, self.thisTestsException2)

	def testAssignment(self):
		if (not self.canTest()):
			return
		# This is not to test if useSolver by default equals CholeskyFactorization, if you need to change the default then change the test here. Or test a different eum.
		# Because this is to test function ArbitraryEnum_from_python<…> if assignment to any enum works as expected and if it throws when it should.
		self.assertTrue(self.r.blinkHighlight == yade.EnumClass_BlinkHighlight.NORMAL)
		self.r.blinkHighlight = 'NEVER'
		self.assertTrue(self.r.blinkHighlight == yade.EnumClass_BlinkHighlight.NEVER)
		self.r.blinkHighlight = 2
		self.assertTrue(self.r.blinkHighlight == yade.EnumClass_BlinkHighlight.WEAK)
