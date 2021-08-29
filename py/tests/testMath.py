# -*- coding: utf-8 -*-
# This is the test of all C++ lib/high-precision/MathFunctions.hpp exported to python via py/high-precision/_math.cpp
# (C) 2015 Anton Gladky <gladk@debian.org>
# (C) 2019 Janek Kozicki

import unittest, math, sys
import yade
from yade import math as mth

import testMathHelper

if (yade.config.highPrecisionMpmath):
	print('\n\033[92m' + "Using " + str(yade.math.getRealHPPythonDigits10()) + " decimal digits in python. Importing mpmath" + '\033[0m\n')
	import mpmath


class SimpleTests(unittest.TestCase):

	def needsMpmathAtN(self, N):
		return yade.math.needsMpmathAtN(N)

	def incompleteComplex(self):
		return ('COMPLEX_MP' not in yade.config.features)

	def hasMpfr(self):
		return ('MPFR' in yade.config.features)

	# flags: -Ofast -march=native -mtune=native -fno-associative-math -fno-finite-math-only -fsigned-zeros
	def isFastNative(self):
		return ('FAST_NATIVE' in yade.config.features)

	def nowUsesBoostBinFloat(self, N):
		return (not self.hasMpfr()) and ((yade.math.RealHPConfig.getDigits10(N) > 33) or (yade.math.RealHPConfig.getDigits10(N) in [24, 30]))

	def setUp(self):
		self.testRecordingMode = False  # if 'True' then it will record 'self.newTolerances' maximum errors encountered, to be put later in place of 'self.defaultTolerances'. See function tearDown() below.
		self.printedAlready = set()
		self.nonBoostMPFR = False  # I was testing non-boost MPFR before: /usr/include/eigen3/unsupported/Eigen/MPRealSupport. Might come handy later.
		# If failures appear and function is not broken then increase tolerance a little.
		# yapf: disable
		self.defaultTolerances={
		 #  function decimal places : tolerance factor. Each "10" corresponds to single wrong decimal place. But they are approximate and rounded up.
		 #
		 #                 float   double    long double float128        MPFR_100        MPFR_150     cpp_bin_float_100  cpp_bin_float_150
		 # Real C++ functions
		    "acos"      : {"6": 100 , "15": 100 , "18": 100  , "33": 1000   , "100": 1000  , "150" : 1000  , "100_b" : 1000    , "150_b" : 1000   }
		  , "atanh"     : {"6": 100 , "15": 100 , "18": 100  , "33": 1000   , "100": 1000  , "150" : 1000  , "100_b" : 1000    , "150_b" : 1000   }
		  , "atan"      : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "atan2"     : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "acosh"     : {"6": 100 , "15": 100 , "18": 100  , "33": 1000   , "100": 1000  , "150" : 1000  , "100_b" : 1000    , "150_b" : 1000   }
		  , "asin"      : {"6": 100 , "15": 100 , "18": 100  , "33": 1000   , "100": 1000  , "150" : 1000  , "100_b" : 1000    , "150_b" : 1000   }
		  , "asinh"     : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }

		 # maybe the error lies in  mpmath, because everything is compared with it.
		  , "sin"       : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "cos"       : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "tan"       : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "tanh"      : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }

		  , "exp"       : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }
		  , "exp2"      : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }
		  , "expm1"     : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }
		  , "cosh"      : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }
		  , "sinh"      : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }

		  , "log"       : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "log10"     : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "log1p"     : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "log2"      : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }

		  , "pow"       : {"6": 5   , "15": 10  , "18": 50   , "33": 50     , "100": 50    , "150" : 50    , "100_b" : 50      , "150_b" : 50     }
		  , "sqrt"      : {"6": 5   , "15": 10  , "18": 50   , "33": 50     , "100": 50    , "150" : 50    , "100_b" : 50      , "150_b" : 50     }

		  , "lgamma"    : {"6": 100 , "15": 500 , "18": 1000 , "33": 10000  , "100": 100000, "150" : 100000, "100_b" : 1000000 , "150_b" : 1000000}
		  , "tgamma"    : {"6": 100 , "15": 100 , "18": 1000 , "33": 10000  , "100": 100000, "150" : 100000, "100_b" : 1000000 , "150_b" : 1000000}
		  , "erfc"      : {"6": 100 , "15": 100 , "18": 2000 , "33": 20000  , "100": 200000, "150" : 200000, "100_b" : 4000000 , "150_b" : 8000000}
		  , "erf"       : {"6": 100 , "15": 100 , "18": 5    , "33": 5      , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }

		  , "modf"      : {"6": 10  , "15": 100 , "18": 5000 , "33": 300000 , "100": 10000 , "150" : 100000, "100_b" : 10000   , "150_b" : 10000  }
		  , "fmod"      : {"6": 10  , "15": 100 , "18": 5000 , "33": 10000  , "100": 10000 , "150" : 100000, "100_b" : 10000   , "150_b" : 10000  }
		  , "remainder" : {"6": 100 , "15": 5000, "18": 5000 , "33": 10000  , "100": 10000 , "150" : 100000, "100_b" : 10000   , "150_b" : 10000  }
		  , "remquo"    : {"6": 100 , "15": 5000, "18": 5000 , "33": 10000  , "100": 10000 , "150" : 100000, "100_b" : 10000   , "150_b" : 10000  }
		  , "fma"       : {"6": 10  , "15": 100 , "18": 10   , "33": 10     , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 1000   }

		# Same order of functions as in #include <lib/high-precision/MathComplexFunctions.hpp>
		# Complex C++ functions. Start names with "Complex " so that they can sit in the same defaultTolerances dictionary
		  , "Complex conj" : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 1      }
		  , "Complex real" : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 1      }
		  , "Complex imag" : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 1      }
		  , "Complex abs"  : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 1      }

		  , "Complex arg"  : {"6": 1   , "15": 1   , "18": 2    , "33": 4      , "100": 8     , "150" : 8     , "100_b" : 8       , "150_b" : 8      }
		  , "Complex norm" : {"6": 1   , "15": 1   , "18": 2    , "33": 4      , "100": 8     , "150" : 8     , "100_b" : 8       , "150_b" : 8      }
		  , "Complex proj" : {"6": 1   , "15": 1   , "18": 2    , "33": 4      , "100": 8     , "150" : 8     , "100_b" : 8       , "150_b" : 8      }
		  , "Complex polar": {"6": 1   , "15": 1   , "18": 2    , "33": 4      , "100": 8     , "150" : 8     , "100_b" : 8       , "150_b" : 8      }

		  , "Complex sin"  : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex sinh" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex cos"  : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex cosh" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex tan"  : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex tanh" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }

		  , "Complex asin" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex asinh": {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex acos" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex acosh": {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex atan" : {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "Complex atanh": {"6": 100 , "15": 100 , "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }

		  , "Complex exp"  : {"6": 10  , "15": 10  , "18": 10   , "33": 10     , "100": 10    , "150" : 10    , "100_b" : 100     , "150_b" : 100    }
		  , "Complex log"  : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "Complex log10": {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "Complex pow"  : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }
		  , "Complex sqrt" : {"6": 1000, "15": 1000, "18": 100  , "33": 100    , "100": 100   , "150" : 100   , "100_b" : 100     , "150_b" : 100    }

		# MathSpecialFunctions
		  , "cylBesselJ"        : {"6": 1000, "15": 1000, "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "factorial"         : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "laguerre"          : {"6": 1000, "15": 1000, "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }
		  , "sphericalHarmonic" : {"6": 1000, "15": 1000, "18": 20000, "33": 4000   , "100": 80000 , "150" : 80000 , "100_b" : 800000  , "150_b" : 800000 }

		 # these are not tolerances. These are EigenCostRealHP from lib/high-precision/EigenNumTraits.hpp
		  , "read"      : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 10000 , "150" : 10000 , "100_b" : 10000   , "150_b" : 10000  }
		  , "add"       : {"6": 1   , "15": 1   , "18": 1    , "33": 2      , "100": 10000 , "150" : 10000 , "100_b" : 10000   , "150_b" : 10000  }
		  , "mul"       : {"6": 1   , "15": 1   , "18": 1    , "33": 2      , "100": 10000 , "150" : 10000 , "100_b" : 10000   , "150_b" : 10000  }
		  , "cread"     : {"6": 2   , "15": 2   , "18": 2    , "33": 2      , "100": 20000 , "150" : 20000 , "100_b" : 20000   , "150_b" : 20000  }
		  , "cadd"      : {"6": 2   , "15": 2   , "18": 2    , "33": 4      , "100": 20000 , "150" : 20000 , "100_b" : 20000   , "150_b" : 20000  }
		  , "cmul"      : {"6": 6   , "15": 6   , "18": 6    , "33": 12     , "100": 60000 , "150" : 60000 , "100_b" : 60000   , "150_b" : 60000  }

		 # Euler–Mascheroni and Pi constants need higher tolerance in boost cpp_bin_float at very high precisions
		  , "euler"     : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "pi"        : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "logE2"     : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  , "catalan"   : {"6": 1   , "15": 1   , "18": 1    , "33": 1      , "100": 1     , "150" : 1     , "100_b" : 1       , "150_b" : 50     }
		  }
		# yapf: enable
		self.testLevelsHP = mth.RealHPConfig.getSupportedByMinieigen()
		self.baseDigits = mth.RealHPConfig.getDigits10(1)
		self.use33or30 = (33 if mth.RealHPConfig.isFloat128Present else 30)
		self.builtinHP = {
		        6: [6, 15, 18, 24, self.use33or30],
		        15: [15, self.use33or30]
		}  # higher precisions are multiplies of baseDigits, see NthLevelRealHP in lib/high-precision/RealHP.hpp
		if (self.testRecordingMode):
			self.startRecordingErrors()

	def tearDown(self):
		if (self.testRecordingMode):
			if (self.totalCount != 0):
				import pickle
				fname = "/tmp/" + str(
				        self.id().split('.')[-1] + "_dig" + str(self.baseDigits) + self.extraName + "_ex" +
				        str(mth.RealHPConfig.extraStringDigits10) + "__" + str(self.totalCount)
				) + ".pickle"
				print(str(self.newTolerances) + "\n\n saving: " + fname)
				pickle.dump(self.newTolerances, open(fname, "wb"))
				# this is how I sort by error, to find the worst performing functions:
				# d=pickle.load( open( "testMathFunctions_dig33_ex4__18000.pickle", "rb" ) );sorted([(d[key]['33'][0],key) for key in d], key=lambda tup: tup[0])

	def getDefaultTolerance(self, name, multiplyByTolerance=True):
		mult = self.tolerance
		key = str(self.digs0) + self.extraName
		if (self.testRecordingMode):
			if (name in self.newTolerances and key in self.newTolerances[name] and self.getMpmath().isfinite(self.newTolerances[name][key][0])):
				return self.newTolerances[name][key][0] * mult
		if (not multiplyByTolerance):
			mult = 1
		dictForThisFunc = self.defaultTolerances[name]
		if (key in dictForThisFunc):
			return dictForThisFunc[key] * mult
		## lower than 33 digits are all hardware precision: 6, 15, 18, 33 digits. But 4*float is 24 digits, and it can be achieved by MPFR only so add exception for 24 also.
		self.assertTrue(self.digs0 > 33 or self.digs0 in [24, 30])  ## 33 was here before
		low = dictForThisFunc["100" + self.extraName]
		high = dictForThisFunc["150" + self.extraName]
		import numpy
		return numpy.interp(self.digs0, [100, 150], [low, high]) * mult

	def storeArgs(self, args):
		self.storedArgs = args

	def storeDefaultTolerance(self, error, name):
		newFactor = (error / self.tolerance) * self.getMpmath().mpf(1.01)
		oldFactor = newFactor
		key = str(self.digs0) + self.extraName
		if (name in self.newTolerances):
			if (key in self.newTolerances[name] and self.getMpmath().isfinite(self.newTolerances[name][key][0])):
				oldFactor = self.newTolerances[name][key][0]
		else:
			self.newTolerances[name] = {}
		if (newFactor >= oldFactor):
			self.newTolerances[name][key] = (newFactor, self.storedArgs)

	def lgamma(self, r):
		mpmath.mp.dps = mth.RealHPConfig.getDigits10(self.maxN) + mth.RealHPConfig.extraStringDigits10
		return mpmath.log(abs(mpmath.gamma(r)))

	def startRecordingErrors(self):
		self.newTolerances = {}
		self.maxN = mth.RealHPConfig.getSupportedByMinieigen()[-1]
		if (mth.RealHPConfig.getDigits10(self.maxN) < 490):
			print(
			        "\n*****\nWarning: recording errors uses less than 490 digits precision. See commits ef1fed55f 015292c0a, they were removed afer this error search was finished.\n*****\n"
			)
		self.maxHPn = getattr(mth, "HP" + str(self.maxN))
		# make maxHPn very similar to mpmath - emulate it.
		self.maxHPn.mpf = mpmath.mpf
		self.maxHPn.mpc = mpmath.mpc
		self.maxHPn.power = self.maxHPn.pow
		self.maxHPn.pi = self.maxHPn.Pi()
		self.maxHPn.euler = self.maxHPn.Euler()
		self.maxHPn.catalan = self.maxHPn.Catalan()
		if (
		        self.maxN > mth.RealHPConfig.workaroundSlowBoostBinFloat
		):  # these functions are unavailable in C++ (because 'import minieigenHP' is too slow for cpp_bin_float), so emulate them
			self.maxHPn.lgamma = self.lgamma
			self.maxHPn.tgamma = mpmath.gamma
			self.maxHPn.gamma = mpmath.gamma
			self.maxHPn.erf = mpmath.erf
			self.maxHPn.erfc = mpmath.erfc
		else:
			self.maxHPn.gamma = self.maxHPn.tgamma
		self.storeArgs((mpmath.mpf('nan'),))
		self.totalCount = 0

	def testBasicHP(self):
		if (self.testRecordingMode):
			return  # skip this test if recording.
		if (mth.RealHPConfig.isEnabledRealHP):
			ec = (1, 2, 3, 4, 8, 10, 20)  # (1,2,3,4,5,6,7,8,9,10,20) #
			mn = (1, 2)  # ec                        # use these if changed something in lib/high-precision/RealHPConfig.hpp
			self.assertEqual(ec, mth.RealHPConfig.getSupportedByEigenCgal())
			self.assertEqual(mn, mth.RealHPConfig.getSupportedByMinieigen())
			if (not self.hasMpfr()):
				self.assertEqual(2, mth.RealHPConfig.workaroundSlowBoostBinFloat)
			else:
				self.assertEqual(ec[-1], mth.RealHPConfig.workaroundSlowBoostBinFloat)
		else:
			self.assertEqual((1,), mth.RealHPConfig.getSupportedByEigenCgal())
			self.assertEqual((1,), mth.RealHPConfig.getSupportedByMinieigen())
			self.assertEqual(1, mth.RealHPConfig.workaroundSlowBoostBinFloat)

	def getDigitsHP(self, N):
		ret = None
		if (self.baseDigits in self.builtinHP) and (N <= len(self.builtinHP[self.baseDigits])):
			ret = self.builtinHP[self.baseDigits][N - 1]
		else:
			ret = self.baseDigits * N
		self.assertEqual(ret, mth.RealHPConfig.getDigits10(N))
		return ret

	def adjustDigs0(self, N, HPn, MPn):
		self.HPnHelper = HPn
		self.digs0 = self.getDigitsHP(N)
		# tolerance = 1.2×10⁻ᵈ⁺¹, where ᵈ==self.digs0
		# so basically we store one more decimal digit, and expect one less decimal digit. That amounts to ignoring one (two, if the extra one is counted) least significant digits.
		self.tolerance = (MPn.mpf(10)**(-self.digs0 + 1)) * MPn.mpf("1.2")
		#self.bits       = MPn.ceil(MPn.mpf(self.digs0)/(MPn.log(2)/MPn.log(10)))+1 # Maybe a bug report against MPFR + cpp_bin_float? They don't use this formula for number of bits
		self.bits = MPn.ceil(
		        MPn.mpf(self.digs0) / (0.301)
		) + 1  # it is reproducing MPFR's formula for number of bits. Discovered by experiments. Adjustments are possible.
		mpmathVsMpfrBits = int(
		        self.bits / 2085
		)  # adjust discrepency between mpmath and MPFR due to incorrect log10/log2 value (above line). The 2085 was found empirically.
		# mpmath has 5 more internal bits, use its mechanisms to extract epsilon
		self.getMpmath().mp.dps = self.digs0 + 1
		self.expectedEpsilon = (2**5) * self.getMpmath().eps() / (2**mpmathVsMpfrBits)
		# now go back to using extraStringDigits10
		self.getMpmath().mp.dps = self.digs0 + mth.RealHPConfig.extraStringDigits10
		if (self.digs0 == 6):  # float case
			self.bits = 24
			self.expectedEpsilon = 1.1920928955078125e-07
		if (self.digs0 == 15):  # double case
			self.bits = 53
			self.expectedEpsilon = 2.220446049250313e-16
		if (self.digs0 == 18):  # long double case
			self.bits = 64
			self.expectedEpsilon = MPn.mpf('1.084202172485504433993e-19')
		if (self.digs0 == 33):  # float128 case
			self.bits = 113
			self.expectedEpsilon = MPn.mpf('1.925929944387235853055977942584926994e-34')
		if (self.needsMpmathAtN(N)):
			self.maxval = (MPn.mpf(1) - self.expectedEpsilon) * MPn.power(2, HPn.max_exp2)
		else:
			import sys
			self.maxval = sys.float_info.max
		if (self.nowUsesBoostBinFloat(N)):
			self.extraName = "_b"
		else:
			self.extraName = ""

	def getMpmath(self):
		if (self.needsMpmathAtN(self.currentN)):
			return mpmath
		else:
			return testMathHelper

	def runCheck(self, N, func):
		if (self.testRecordingMode and N == self.maxN):
			return  # no need to test maxN against itself. It is for testing lower precisions against it.
		self.currentN = N
		# the same as the line 'std::string name = "HP" + boost::lexical_cast<std::string>(N)' in function registerInScope in _math.cpp
		HPn = getattr(mth, "HP" + str(N))
		if (not self.testRecordingMode):
			MPn = self.getMpmath()
		else:
			MPn = self.maxHPn  # we are recording the errors, do all the tests against the max precision available
		if (N == 1):
			self.adjustDigs0(N, mth, MPn)
			func(N, mth, MPn)  # test global scope functions with RealHP<1>
		self.adjustDigs0(N, HPn, MPn)
		func(N, HPn, MPn)  # test scopes HP1, HP2, etc

	def printOnce(self, functionName, a):
		MPn = self.getMpmath()
		if (functionName and (functionName not in self.printedAlready) and (not MPn.isnan(abs(a)))):
			self.printedAlready.add(functionName)
			print(functionName.ljust(15) + " : " + a.__repr__())

	def checkRelativeError(self, a, b, tol=None, functionName=None, isComplex=False):
		if (functionName and self.incompleteComplex() and functionName[0:7] == "Complex"):
			# don't check complex functions
			return
		MPn = self.getMpmath()
		prevDps = MPn.mp.dps
		if (self.testRecordingMode):
			MPn.mp.dps = mth.RealHPConfig.getDigits10(self.maxN) + mth.RealHPConfig.extraStringDigits10
		denominator = max(abs(a), abs(b))  # avoid division by zero
		if (denominator == 0):  # they are both equal to zero
			error = 0
		else:
			if isComplex:
				error = abs((MPn.mpc(a) - MPn.mpc(b)) / MPn.mpc(denominator))
			else:
				error = abs((MPn.mpf(a) - MPn.mpf(b)) / MPn.mpf(denominator))
		if (abs(b) <= self.maxval and abs(b) >= self.HPnHelper.smallest_positive()):
			#print("a= ",a," b= ",b," smallest=",self.HPnHelper.smallest_positive(), " maxval=",self.maxval)
			self.printOnce(functionName, a)
			if ((not MPn.isfinite(a)) or (not MPn.isfinite(b))):
				if ((functionName != "lgamma") and (not self.testRecordingMode)):  # lgamma triggers this warning too often.
					print(
					        "\033[93m Warning: \033[0m got NaN or Inf, cannot verify if: ", a, " == ", b,
					        " that was for function: \033[93m ", functionName, " \033[0m"
					)
			else:
				if (tol != None):
					#print("a=",a," b=",b," tol=",tol)
					self.assertLessEqual(error, tol)
				else:
					if (functionName in self.defaultTolerances):
						if (self.testRecordingMode):
							self.storeDefaultTolerance(error, functionName)
						defaultToleranceForThisFunction = self.getDefaultTolerance(functionName)
						#print(defaultToleranceForThisFunction," ---- ",functionName)
						self.assertLessEqual(error, defaultToleranceForThisFunction)
					else:
						self.assertLessEqual(error, self.tolerance)
		elif (not self.testRecordingMode):
			print(
			        "Skipping ", functionName, " check, the builtin number: ", a, " cannot have value outside of its possible repesentation: ", b,
			        ", because it has only ", self.digs0, " digits."
			)
		MPn.mp.dps = prevDps

	def checkRelativeComplexError(self, a, b, tol=None, functionName=None):
		self.printOnce(functionName, a)
		self.checkRelativeError(abs(a), abs(b), tol, functionName, True)

	def oneArgMathCheck(self, N, HPn, MPn, r):
		# note: cos, tan, sin, lgamma, tgamma get wildly inaccurrate when |arg| > 20. Errors are in the range log₁₀(8000000)≈7 decimal places for most of RealHP<…> types.
		# these functions become more or less useless. So better to measure error in a usable range. I arbitrarily set it to 4*Pi and 20.
		# This strange behavior is explained by the error in the remainder(…) for which I do not restrict arguments (so you can look up its error in the table). These trig
		# functions try to remove periodicity by calculating remainder from division by Pi, but they can only be as good as the remainder calculation itself. And this calculation
		# cannot produce more precision than the number already has, after its first few digits are cut-off by the remainder calculation.
		cut1 = HPn.roundTrip(r % self.getMpmath().pi * 4)  # the HPn.identity(…) call is to cut the digits to those representible in HPn
		cut2 = HPn.roundTrip(r % 20)
		self.checkRelativeError(HPn.sin(cut1), MPn.sin(cut1), functionName="sin")
		self.checkRelativeError(HPn.sinh(r), MPn.sinh(r), functionName="sinh")
		self.checkRelativeError(HPn.cos(cut1), MPn.cos(cut1), functionName="cos")
		self.checkRelativeError(HPn.cosh(r), MPn.cosh(r), functionName="cosh")
		self.checkRelativeError(HPn.tan(cut1), MPn.tan(cut1), functionName="tan")
		self.checkRelativeError(HPn.tanh(r), MPn.tanh(r), functionName="tanh")
		# check math functions, but ensure that input arguments produce real (not complex) results
		self.checkRelativeError(HPn.abs(r), abs(r), functionName="abs")
		self.checkRelativeError(HPn.acos(r % 1), MPn.acos(r % 1), functionName="acos")
		self.checkRelativeError(HPn.acosh(abs(r) + 1), MPn.acosh(abs(r) + 1), functionName="acosh")
		self.checkRelativeError(HPn.asin(r % 1), MPn.asin(r % 1), functionName="asin")
		self.checkRelativeError(HPn.asinh(r), MPn.asinh(r), functionName="asinh")
		self.checkRelativeError(HPn.atan(r), MPn.atan(r), functionName="atan")
		self.checkRelativeError(HPn.atanh(r % 1), MPn.atanh(r % 1), functionName="atanh")
		self.checkRelativeError(HPn.cbrt(abs(r)), MPn.cbrt(abs(r)), functionName="cbrt")
		self.assertEqual(HPn.ceil(r), MPn.ceil(r))
		self.checkRelativeError(HPn.exp(r), MPn.exp(r), functionName="exp")
		self.checkRelativeError(HPn.sqrt(abs(r)), MPn.sqrt(abs(r)), functionName="sqrt")
		self.checkRelativeError(HPn.exp2(r), MPn.power(2, r), functionName="exp2")
		self.checkRelativeError(HPn.expm1(r), MPn.expm1(r), functionName="expm1")
		self.assertEqual(HPn.floor(r), MPn.floor(r))
		#print(HPn.ilogb(r).__repr__()) # ilogb is not present in mpmath
		if (N <= mth.RealHPConfig.workaroundSlowBoostBinFloat):
			#print(" N=",N , " digits10=", yade.math.RealHPConfig.getDigits10(N) ,"  self.expectedEpsilon = ",self.expectedEpsilon, " r=",r)
			if (self.testRecordingMode):
				self.checkRelativeError(HPn.lgamma(cut2), self.maxHPn.lgamma(cut2), functionName="lgamma")
				self.checkRelativeError(HPn.tgamma(cut2), self.maxHPn.tgamma(cut2), functionName="tgamma")
			else:
				self.checkRelativeError(HPn.lgamma(cut2), MPn.log(abs(MPn.gamma(cut2))), functionName="lgamma")
				self.checkRelativeError(HPn.tgamma(cut2), MPn.gamma(cut2), functionName="tgamma")
			self.checkRelativeError(HPn.erf(r), MPn.erf(r), functionName="erf")
			self.checkRelativeError(HPn.erfc(r), MPn.erfc(r), functionName="erfc")
		self.checkRelativeError(HPn.log(abs(r) + self.tolerance), MPn.log(abs(r) + self.tolerance), functionName="log")
		self.checkRelativeError(HPn.log10(abs(r) + self.tolerance), MPn.log10(abs(r) + self.tolerance), functionName="log10")
		self.checkRelativeError(HPn.log1p(abs(r) + self.tolerance), MPn.log(1 + abs(r) + self.tolerance), functionName="log1p")
		self.checkRelativeError(HPn.log2(abs(r) + self.tolerance), MPn.log(abs(r) + self.tolerance) / MPn.log(2), functionName="log2")
		#print(HPn.logb(r).__repr__()) # logb is not present in mpmath
		self.assertEqual(HPn.rint(r), round(r))
		self.assertTrue((HPn.round(r) == round(r)) or (r % 1 == 0.5))  # ignore rounding 0.5 up or down.
		self.assertEqual(HPn.trunc(abs(r)), int(abs(r)))

		self.checkRelativeError(HPn.fabs(r), abs(r), functionName="fabs")

		pair = HPn.frexp(abs(r))
		self.checkRelativeError(abs(r), pair[0] * MPn.power(2, pair[1]), functionName="frexp")

		pair = HPn.modf(abs(r))
		self.checkRelativeError(pair[0], (abs(r)) % 1, functionName="modf")
		self.assertEqual(pair[1], int(abs(r)))

		#self.assertEqual(HPn.frexp(abs(r)),HPn.frexp_c_test(abs(r)))
		#self.assertEqual(HPn.modf(abs(r)),HPn.modf_c_test(abs(r)))

		if (r == 0):
			self.assertEqual(HPn.sgn(r), 0)
			self.assertEqual(HPn.sign(r), 0)
		if (r > 0):
			self.assertEqual(HPn.sgn(r), 1)
			self.assertEqual(HPn.sign(r), 1)
		if (r < 0):
			self.assertEqual(HPn.sgn(r), -1)
			self.assertEqual(HPn.sign(r), -1)

		self.checkCgalNumTraits(HPn, MPn, r)

		# One arg MathSpecialFunctions
		self.checkRelativeError(HPn.factorial(int(abs(r))), MPn.factorial(int(abs(r))), functionName="factorial")

	def checkCgalNumTraits(self, HPn, MPn, r):
		if (HPn.testCgalNumTraits == False):
			print("Skipping test of CgalNumTraits")
			return
		self.assertEqual(HPn.CGAL_Is_valid(r), True)
		if (r != 0):
			self.checkRelativeError(HPn.CGAL_Square(r), MPn.power(r, 2), functionName="pow")
			self.checkRelativeError(HPn.CGAL_Sqrt(abs(r)), MPn.sqrt(abs(r)), functionName="sqrt")
			for kk in range(5):
				k = kk + 1
				self.checkRelativeError(HPn.CGAL_Kth_root(k, abs(r)), MPn.power(abs(r), 1 / MPn.mpf(k)), functionName="pow")
			# CGAL uses double for intervals
			interval = HPn.CGAL_To_interval(r)
			self.checkRelativeError(r, interval[0], 1e-14)
			self.checkRelativeError(r, interval[1], 1e-14)
		self.assertEqual(HPn.CGAL_Is_finite(r), True)
		if (r == 0):
			self.assertEqual(HPn.CGAL_Sgn(r), 0)
		if (r > 0):
			self.assertEqual(HPn.CGAL_Sgn(r), 1)
		if (r < 0):
			self.assertEqual(HPn.CGAL_Sgn(r), -1)
		self.assertEqual(HPn.CGAL_Sgn(0), 0)
		self.assertEqual(HPn.CGAL_Sgn(2.5), 1)
		self.assertEqual(HPn.CGAL_Sgn(-2.3), -1)

	def testInfinityNaN(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestInfinityNaN)

	def HPtestInfinityNaN(self, N, HPn, MPn):
		if (HPn.hasInfinityNan == False):
			print("Skipping inf,nan regular test\n")
			print("\033[91m *** Warning: usually YADE needs Inf and NaN for most of the calculations. *** \033[0m")
			return
		self.assertEqual(HPn.isinf(MPn.mpf(1)), False)
		self.assertEqual(HPn.isinf(MPn.mpf('nan')), False)
		self.assertEqual(HPn.isinf(MPn.mpf('inf')), True)
		self.assertEqual(HPn.isnan(MPn.mpf(1)), False)
		self.assertEqual(HPn.isnan(MPn.mpf('nan')), True)
		self.assertEqual(HPn.isnan(MPn.mpf('inf')), False)
		self.assertEqual(HPn.isfinite(MPn.mpf(1)), True)
		self.assertEqual(HPn.isfinite(MPn.mpf('nan')), False)
		self.assertEqual(HPn.isfinite(MPn.mpf('inf')), False)

	def testRealHPDiagnostics(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestRealHPDiagnostics)

	def HPtestRealHPDiagnostics(self, N, HPn, MPn):
		import random
		source = '1'
		for i in range(mth.getDigits2(N) - 1):
			source += random.choice(['0', '1'])
		for m in (
		        list(self.testLevelsHP) + list(reversed(self.testLevelsHP))
		):  # goes both directions, because that's how I discovered a small mistake in ToFromPythonConverter.hpp
			if (m >= N):
				HPm = getattr(mth, "HP" + str(m))
				toHPn = getattr(HPm, "toHP" + str(N))
				testNum = HPm.fromBits(source)
				if (m == N):
					source2 = source[:-1] + '1'
					testNum2 = HPm.fromBits(source2)
					self.assertLessEqual(HPm.getFloatDistanceULP(testNum, testNum2), 1)
					source3 = source[:-2] + '11'
					testNum3 = HPm.fromBits(source3)
					self.assertLessEqual(HPm.getFloatDistanceULP(testNum, testNum3), 3)
				if (mth.getDigits2(m) in [24, 53, 64, 113]):  # don't check MPFR, cpp_bin_float.
					self.assertTrue(
					        source in HPm.getRawBits(testNum)
					)  # checks only 'in', because sign+exponent (usually in front of raw bits) are architecture-specific.
				self.assertEqual(HPm.getDecomposedReal(testNum)['bits'], source.ljust(mth.getDigits2(m), '0'))
				self.assertEqual(HPn.fromBits(source), toHPn(testNum))
		if (mth.getDigits2(N) == 24):
			self.assertEqual(HPn.getDemangledName(), 'float')
		if (mth.getDigits2(N) == 53):
			self.assertEqual(HPn.getDemangledName(), 'double')
		if (mth.getDigits2(N) == 64):
			self.assertTrue('long double' in HPn.getDemangledName())
		if (mth.getDigits2(N) == 113):
			self.assertTrue('float128' in HPn.getDemangledName())

	def bitsToLevelHP(self, bits):
		N = -1
		for nn in self.testLevelsHP:
			if (mth.getDigits2(nn) == bits):
				N = nn
				break
		return N

	def testRealHPErrors(self):
		if (len(self.testLevelsHP) < 2):
			return
		testUlpRandom = yade.math.getRealHPErrors(
		        list(self.testLevelsHP), testCount=10000, minX=-100, maxX=100, useRandomArgs=True, printEveryNth=100000, extraChecks=False
		)
		testUlpLinear = yade.math.getRealHPErrors(
		        list(self.testLevelsHP), testCount=10000, minX=-100, maxX=100, useRandomArgs=False, printEveryNth=100000, extraChecks=False
		)
		self.showRealHPResults(testUlpRandom)
		self.showRealHPResults(testUlpLinear)
		self.processRealHPResults(testUlpRandom)
		self.processRealHPResults(testUlpLinear)

	def showRealHPResults(self, testULP):
		#print(testULP)
		for func in testULP:
			for bits in testULP[func]:
				ulp = testULP[func][bits][1]
				if (ulp > 4):
					N = self.bitsToLevelHP(bits)
					print(
					        "\033[93mWarning:\033[0m ULP error of\033[91m", func, "\033[0musing RealHP<", N, ">, ", bits, "bits, with arg:",
					        testULP[func][bits][0], "is ULP=\033[93m", ulp, "\033[0m"
					)

	def getBoostComplexTolerance(self, func, tol, bits):
		# g++ flags: -Ofast -march=native -mtune=native -fno-associative-math -fno-finite-math-only -fsigned-zeros
		isFast = self.isFastNative()
		mpfr = not self.nowUsesBoostBinFloat(self.bitsToLevelHP(bits))
		longDouble = (bits == 64)  # C++ double has 53
		# The self.defaultTolerances is about errors found on python side. This one is more precise about ULP errors found on C++ side.
		complexTolerancesInUnitsOfULP = {
		        "complex acos imag": 1e5,
		        "complex acos real": 3e5,
		        "complex asin imag": 3e5,
		        "complex asin real": 3e5,
		        "complex asinh imag": 3e5,
		        "complex asinh real": 1e5,
		        "complex atan imag": 5e6,
		        "complex atanh real": 5e5 if mpfr else 8e8,  # boost::cpp_bin_float has larger error
		        "complex cos imag": 8 if mpfr else 2e5,
		        "complex cos real": 8 if mpfr else 2e5,
		        "complex cosh imag": 8 if mpfr else 8e4,
		        "complex cosh real": 8 if mpfr else 2e5,
		        "complex exp imag": 8 if mpfr else 7e4,
		        "complex exp real": 8 if mpfr else 2e5,
		        "complex polar imag": 8 if mpfr else 9e4,
		        "complex polar real": 8 if mpfr else 2e5,
		        "complex pow imag": 3e7 if isFast else
		                            (9e6 if longDouble else
		                             (4e6 if mpfr else 7e8)),  # std::complex<double> -Ofast has 3e7 error then std::complex<long double> has 9e6.
		        "complex pow real": 5e7 if self.incompleteComplex() else (4e6 if mpfr else 2e7),
		        "complex sin imag": 8 if mpfr else 2e5,
		        "complex sin real": 8 if mpfr else 2e5,
		        "complex sinh imag": 8 if mpfr else 7e4,
		        "complex sinh real": 8 if mpfr else 2e5,
		        "complex tan imag": 8 if mpfr else 500,
		        "complex tan real": 8 if mpfr else 500,
		        "complex tanh imag": 8 if mpfr else 500,
		        "complex tanh real": 8 if mpfr else 500,
		}
		# FIXME: These tolerances ↑ also need to be reported to boost. But that's a smaller error anyway. the complex tan, tanh had error 1e38 ULP.
		#        The values here are chosen only for std::complex<long double>, complex_adaptor<cpp_bin_float_45>, mpc_complex_150.
		#        Other precisions may vary slightly.
		if ((func in complexTolerancesInUnitsOfULP)):
			tol = complexTolerancesInUnitsOfULP[func]
		return tol

	def getMathSpecialTolerance(self, func, tol, bits):
		# The self.defaultTolerances is about errors found on python side. This one is more precise about ULP errors found on C++ side.
		mpfr = not self.nowUsesBoostBinFloat(self.bitsToLevelHP(bits))
		tolerancesInUnitsOfULP = {
		        "cylBesselJ": 6e5,
		        "laguerre": 1e6,
		        "complex sphericalHarmonic imag": 3e7,
		        "complex sphericalHarmonic real": 3e7 if mpfr else 2e8,
		}
		if ((func in tolerancesInUnitsOfULP)):
			tol = tolerancesInUnitsOfULP[func]
		return tol

	def processRealHPResults(self, testULP):
		for func in testULP:
			for bits in testULP[func]:
				tolerateErrorULP = 8
				if (self.nowUsesBoostBinFloat(self.bitsToLevelHP(bits))):
					tolerateErrorULP = 127  # cpp_bin_float has larger errors
					if (func in ["tgamma", "acos", "erfc"]):
						tolerateErrorULP = 50000
					elif (func == "lgamma"):
						tolerateErrorULP = 1e10
					elif (func in ["sin", "cos", "tan", "fma"]):
						tolerateErrorULP = 2e8

				boostVer = yade.libVersions.getVersion('boost')
				tolerateErrorULP = self.getBoostComplexTolerance(func, tolerateErrorULP, bits)

				tolerateErrorULP = self.getMathSpecialTolerance(func, tolerateErrorULP, bits)

				# DONE: file a bug report about higher precision versions of these two functions. They have large error: log2(300000000)≈28.1 incorrect bits.
				#       https://github.com/boostorg/multiprecision/issues/262
				#       https://github.com/boostorg/multiprecision/issues/264
				# when it's fixed we can check boost version and skip this line below.
				if ((func in ["complex tan real", "complex tanh imag"]) and (boostVer < (1, 76, 0))):
					# for older boost there's nothing we can do. We have a large error: 1e38 ULP
					pass
				else:
					ulp = testULP[func][bits][1]
					if (ulp > tolerateErrorULP):
						print("\n\033[91mError with: ", func, "\033[0m having ULP=", ulp)
					self.assertLessEqual(ulp, tolerateErrorULP)

	def testCgalNumTraits(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestCgalNumTraits)

	def HPtestCgalNumTraits(self, N, HPn, MPn):
		if (HPn.testCgalNumTraits == False):
			print("Skipping test of CgalNumTraits")
			return
		self.checkCgalNumTraits(HPn, MPn, 0)
		self.checkCgalNumTraits(HPn, MPn, 0.5)
		self.checkCgalNumTraits(HPn, MPn, -1.5)
		self.checkCgalNumTraits(HPn, MPn, 55.5)
		self.assertEqual(HPn.CGAL_Is_valid(MPn.mpf(1)), True)
		self.assertEqual(HPn.CGAL_Is_valid(MPn.mpf('nan')), False)
		self.assertEqual(HPn.CGAL_Is_valid(MPn.mpf('inf')), True)
		self.assertEqual(HPn.CGAL_Is_finite(MPn.mpf(1)), True)
		self.assertEqual(HPn.CGAL_Is_finite(MPn.mpf('nan')), False)
		self.assertEqual(HPn.CGAL_Is_finite(MPn.mpf('inf')), False)
		self.assertEqual(HPn.CGAL_simpleTest(), MPn.mpf("3.0"))

	def twoArgMathCheck(self, HPn, MPn, r1, r2):
		# same order of complex functions as in lib/high-precision/MathComplexFunctions.hpp , py/high-precision/_math.cpp
		self.checkRelativeComplexError(HPn.conj(MPn.mpc(r1, r2)), MPn.conj(MPn.mpc(r1, r2)), functionName="Complex conj")
		self.checkRelativeComplexError(HPn.real(MPn.mpc(r1, r2)), r1, functionName="Complex real")
		self.checkRelativeComplexError(HPn.imag(MPn.mpc(r1, r2)), r2, functionName="Complex imag")
		self.checkRelativeComplexError(HPn.abs(MPn.mpc(r1, r2)), abs(MPn.mpc(r1, r2)), functionName="Complex abs")

		self.checkRelativeError(HPn.arg(MPn.mpc(r1, r2)), MPn.phase(MPn.mpc(r1, r2)), functionName="Complex arg")
		self.checkRelativeError(HPn.squaredNorm(MPn.mpc(r1, r2)), MPn.norm(MPn.mpc(r1, r2)) * MPn.norm(MPn.mpc(r1, r2)), functionName="Complex norm")
		# for now skip testing C++ std::proj, see note in py/tests/testMathHelper.py
		self.checkRelativeComplexError(HPn.polar(r1, r2), MPn.rect(r1, r2), functionName="Complex polar")

		self.checkRelativeComplexError(HPn.sin(MPn.mpc(r1, r2)), MPn.sin(MPn.mpc(r1, r2)), functionName="Complex sin")
		self.checkRelativeComplexError(HPn.sinh(MPn.mpc(r1, r2)), MPn.sinh(MPn.mpc(r1, r2)), functionName="Complex sinh")
		self.checkRelativeComplexError(HPn.cos(MPn.mpc(r1, r2)), MPn.cos(MPn.mpc(r1, r2)), functionName="Complex cos")
		self.checkRelativeComplexError(HPn.cosh(MPn.mpc(r1, r2)), MPn.cosh(MPn.mpc(r1, r2)), functionName="Complex cosh")
		self.checkRelativeComplexError(HPn.tan(MPn.mpc(r1, r2)), MPn.tan(MPn.mpc(r1, r2)), functionName="Complex tan")
		self.checkRelativeComplexError(HPn.tanh(MPn.mpc(r1, r2)), MPn.tanh(MPn.mpc(r1, r2)), functionName="Complex tanh")

		self.checkRelativeComplexError(HPn.asin(MPn.mpc(r1, r2)), MPn.asin(MPn.mpc(r1, r2)), functionName="Complex asin")
		self.checkRelativeComplexError(HPn.asinh(MPn.mpc(r1, r2)), MPn.asinh(MPn.mpc(r1, r2)), functionName="Complex asinh")
		self.checkRelativeComplexError(HPn.acos(MPn.mpc(r1, r2)), MPn.acos(MPn.mpc(r1, r2)), functionName="Complex acos")
		self.checkRelativeComplexError(HPn.acosh(MPn.mpc(r1, r2)), MPn.acosh(MPn.mpc(r1, r2)), functionName="Complex acosh")
		self.checkRelativeComplexError(HPn.atan(MPn.mpc(r1, r2)), MPn.atan(MPn.mpc(r1, r2)), functionName="Complex atan")
		self.checkRelativeComplexError(HPn.atanh(MPn.mpc(r1, r2)), MPn.atanh(MPn.mpc(r1, r2)), functionName="Complex atanh")

		self.checkRelativeComplexError(HPn.exp(MPn.mpc(r1, r2)), MPn.exp(MPn.mpc(r1, r2)), functionName="Complex exp")
		self.checkRelativeComplexError(HPn.log(MPn.mpc(r1, r2)), MPn.log(MPn.mpc(r1, r2)), functionName="Complex log")
		self.checkRelativeComplexError(HPn.log10(MPn.mpc(r1, r2)), MPn.log10(MPn.mpc(r1, r2)), functionName="Complex log10")
		self.checkRelativeComplexError(
		        HPn.pow(MPn.mpc(r1, r2), MPn.mpc(r1 + r2, r1 - r2)), (MPn.mpc(r1, r2)**MPn.mpc(r1 + r2, r1 - r2)), functionName="Complex pow"
		)
		self.checkRelativeComplexError(HPn.sqrt(MPn.mpc(r1, r2)), MPn.sqrt(MPn.mpc(r1, r2)), functionName="Complex sqrt")

		# Two argument MathSpecialFunctions
		if (self.needsMpmathAtN(self.currentN)):  # can test Bessel only if mpmath is available.
			# first is yade.math function here ↓ , next is mpmath function here ↓
			self.checkRelativeComplexError(HPn.cylBesselJ(int(abs(r1)), r2), MPn.besselj(int(abs(r1)), r2), functionName="cylBesselJ")

		# Other, non complex functions
		self.checkRelativeError(HPn.atan2(r1, r2), MPn.atan2(r1, r2), functionName="atan2")
		self.checkRelativeError(HPn.fmod(abs(r1), abs(r2)), MPn.fmod(abs(r1), abs(r2)), functionName="fmod")
		self.checkRelativeError(HPn.hypot(r1, r2), MPn.hypot(r1, r2), functionName="hypot")
		self.checkRelativeError(HPn.max(r1, r2), max(r1, r2), functionName="max")
		self.checkRelativeError(HPn.min(r1, r2), min(r1, r2), functionName="min")
		self.checkRelativeError(HPn.pow(abs(r1), r2), MPn.power(abs(r1), r2), functionName="pow")
		self.checkRelativeError(HPn.remainder(abs(r1), abs(r2)), abs(r1) - round(abs(r1) / abs(r2)) * abs(r2), functionName="remainder")
		pair = HPn.remquo(abs(r1), abs(r2))
		self.checkRelativeError(pair[0], abs(r1) - round(abs(r1) / abs(r2)) * abs(r2), functionName="remquo")
		self.assertEqual(pair[1] % 8, round(abs(r1 / r2)) % 8)

		self.checkRelativeError(HPn.ldexp(r1, int(r2)), MPn.mpf(r1) * MPn.power(2, int(r2)), functionName="ldexp")

	def threeArgMathCheck(self, HPn, MPn, r1, r2, r3):
		self.checkRelativeError(HPn.fma(r1, r2, r3), (MPn.mpf(r1) * r2) + r3, functionName="fma")

		# Three argument MathSpecialFunctions
		if (yade.config.highPrecisionMpmath == True):  # can test Bessel and sphericalHarmonic only if mpmath is available.
			# first is yade.math function here ↓ , next is mpmath function here ↓
			self.checkRelativeComplexError(
			        HPn.laguerre(int(abs(r1)), int(abs(r2)), r3), mpmath.laguerre(int(abs(r1)), int(abs(r2)), r3), functionName="laguerre"
			)
			theta = r2 % self.getMpmath().pi
			phi = r2 % (self.getMpmath().pi * 2)
			self.checkRelativeComplexError(
			        HPn.sphericalHarmonic(int(abs(r1)), int(abs(r1)), theta, phi),
			        mpmath.spherharm(int(abs(r1)), int(abs(r1)), theta, phi),
			        functionName="sphericalHarmonic"
			)

	def testMathFunctions(self):
		for N in self.testLevelsHP:
			self.printedAlready = set()
			self.runCheck(N, self.HPtestMathFunctions)

	def HPtestMathFunctions(self, N, HPn, MPn):
		self.assertEqual(HPn.defprec, self.bits)
		zz = MPn.acos(0)
		#print(zz.__repr__())
		#print("zz:",hex(id(zz)))
		#print("mpmath:",hex(id(mpmath)))
		a = HPn.Var()
		a.val = zz
		if (not self.testRecordingMode):
			self.assertEqual(MPn.mp.dps, self.digs0 + mth.RealHPConfig.extraStringDigits10)
		#print("---- a.val=",a.val.__repr__())
		#print("---- zz   =",zz   .__repr__())
		#print("---- DPS  =",mpmath.mp.dps)
		#print("---- abs  =",abs(MPn.mpf(a.val-zz)))
		#print("---- 10** =",self.tolerance)
		self.checkRelativeError(a.val, zz)
		self.assertEqual(HPn.IsInteger, 0)
		self.assertEqual(HPn.IsSigned, 1)
		self.assertEqual(HPn.IsComplex, 0)
		if (self.bits >= 64):
			self.assertEqual(HPn.RequireInitialization, 1)
		else:
			self.assertEqual(HPn.RequireInitialization, 0)
		self.assertGreaterEqual(HPn.ReadCost, 1)
		self.assertGreaterEqual(HPn.AddCost, 1)
		self.assertGreaterEqual(HPn.MulCost, 1)
		self.checkRelativeError(HPn.highest(), self.maxval, 2.1)
		self.checkRelativeError(-HPn.lowest(), self.maxval, 2.1)
		self.checkRelativeError(HPn.Pi(), MPn.pi, functionName="pi")
		print("HPn.Euler() ", HPn.Euler(), " N=", N, "  MPn.euler = ", str(MPn.euler))
		self.checkRelativeError(HPn.Euler(), MPn.euler, functionName="euler")
		self.checkRelativeError(HPn.Log2(), MPn.log(2), functionName="logE2")
		self.checkRelativeError(HPn.Catalan(), MPn.catalan, functionName="catalan")
		#print("HPn.epsilon() ",HPn.epsilon(), " N=",N ,"  self.expectedEpsilon = ",self.expectedEpsilon)
		self.checkRelativeError(HPn.epsilon(), self.expectedEpsilon, 10)
		if (self.digs0 == 6):  # exception for float
			self.assertLessEqual(HPn.dummy_precision(), 10e-6)
		else:
			self.checkRelativeError(MPn.log(HPn.dummy_precision() / HPn.epsilon()) / MPn.log(10), MPn.mpf(self.digs0) / 10, 1.5)
		for x in range(50):
			if (self.nonBoostMPFR):  # this looks like a bug in /usr/include/eigen3/unsupported/Eigen/MPRealSupport !
				self.assertLessEqual(abs(HPn.random() - 0.5), 0.5)
			else:
				self.assertLessEqual(abs(HPn.random()), 1.0)
		for aa in range(1000000 if self.testRecordingMode else 20):
			m = 100 if self.testRecordingMode else 20
			r = HPn.random(-m, m)
			r2 = HPn.random(-m, m)
			r3 = HPn.random(-m, m)
			if (self.testRecordingMode):
				self.totalCount = self.totalCount + 1
				self.storeArgs((r, r2, r3))
			#print("random=",r)
			self.assertLessEqual(r, m)
			self.assertGreaterEqual(r, -m)
			self.assertFalse(HPn.isMuchSmallerThan(r, 1, HPn.epsilon()))
			# NOTE: Below is a very sensitive test. If it starts failing, then see in function adjustDigs0, how expectedEpsilon is calculated.
			# Maybe MPFR or cpp_bin_float changed the number of bits or changed their internal approximation of log10/log2.
			self.assertTrue(HPn.isMuchSmallerThan(self.expectedEpsilon, 1 + abs(r), HPn.epsilon()))
			self.assertTrue(HPn.isEqualFuzzy(r + self.expectedEpsilon * 0.01, r, HPn.epsilon()))
			self.checkRelativeError(HPn.toLongDouble(r), float(r), 1e-14)  # FIXME - should be 1e-17, but python does not support that
			self.checkRelativeError(HPn.toDouble(r), float(r), 1e-14)
			self.checkRelativeError(HPn.toDouble(r), float(r), 1e-14)
			self.assertEqual(HPn.toLong(r), int(r))
			self.assertEqual(HPn.toInt(r), int(r))
			#
			#print(r.__repr__(),r2.__repr__(),r3.__repr__())
			self.oneArgMathCheck(N, HPn, MPn, r)
			self.oneArgMathCheck(N, HPn, MPn, r2)
			self.oneArgMathCheck(N, HPn, MPn, r3)
			self.twoArgMathCheck(HPn, MPn, r, r2)
			self.twoArgMathCheck(HPn, MPn, r, r3)
			self.twoArgMathCheck(HPn, MPn, r2, r3)
			self.threeArgMathCheck(HPn, MPn, r, r2, r3)
			if (self.testRecordingMode and (self.totalCount % 20000 == 0)):
				self.tearDown()  # save progress

	def testConstantsPythonSide(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestConstants)

	def HPtestConstants(self, N, HPn, MPn):
		HPn.testConstants()

	def testArray(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestArray)

	def HPtestArray(self, N, HPn, MPn):
		HPn.testArray()

	def testConstantsCppSide(self):
		mth.testLoopRealHP()

	def testBasicVariable(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestBasicVariable)

	def HPtestBasicVariable(self, N, HPn, MPn):
		a = HPn.Var()
		self.checkRelativeError(a.val, -71.23, 0.01)
		a.val = 10
		self.checkRelativeError(a.val, 10)
		self.checkRelativeComplexError(a.cpl, -71.23 + 33.23j, 0.01)
		a.cpl = MPn.mpc("1", "-1")
		self.checkRelativeComplexError(a.cpl, 1 - 1j, 1e-15)
		self.checkRelativeComplexError(a.cpl, MPn.mpc("1", "-1"))

	def thisTestsExceptionReal(self):
		a = self.HPnHelper.Var()
		a.val = "13123-123123*123"

	def thisTestsExceptionComplex(self):
		a = self.HPnHelper.Var()
		a.cpl = "13123-123123*123-50j"

	def testWrongInput(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestWrongInput)

	def HPtestWrongInput(self, N, HPn, MPn):
		if (self.nonBoostMPFR):  # this looks like another bug in /usr/include/mpreal.h
			print("skipping this test for non-boost /usr/include/mpreal.h")
			return
		# depending on backed Real use it throws TypeError or RuntimeError
		self.HPnHelper = HPn
		self.assertRaises(Exception, self.thisTestsExceptionReal)
		self.assertRaises(Exception, self.thisTestsExceptionComplex)

	def testEigenCost(self):
		for N in self.testLevelsHP:
			self.runCheck(N, self.HPtestEigenCost)

	def HPtestEigenCost(self, N, HPn, MPn):
		self.assertEqual(self.getDefaultTolerance("read", False), HPn.ReadCost)
		self.assertEqual(self.getDefaultTolerance("add", False), HPn.AddCost)
		self.assertEqual(self.getDefaultTolerance("mul", False), HPn.MulCost)
		self.assertEqual(self.getDefaultTolerance("cread", False), HPn.ComplexReadCost)
		self.assertEqual(self.getDefaultTolerance("cadd", False), HPn.ComplexAddCost)
		self.assertEqual(self.getDefaultTolerance("cmul", False), HPn.ComplexMulCost)

	def testDefReadonly(self):
		# boost::python def_readonly(…) has problems with higher precision Real
		# types: https://yade-dem.org/doc/prog.html#custom-converters
		# The solution (already applied almost everywhere in Serializable.hpp)
		# is to use ::boost::python::return_by_value in def_readonly(…) when
		# exposing Real member variable.
		t = yade.TriaxialStressController()
		p = t.porosity
