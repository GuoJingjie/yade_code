#!/usr/bin/python
# -*- coding: utf-8 -*-
##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

# When mpmath is not required implement a super-minimal version of mpmath, so that the tests will work and pretend that they use mpmath
# So this file is in fact to avoid python3-mpmath dependency when high-precision is not used.

import yade

# Without mpmath use a super-minimal version of mpmath, so that all the tests will work and pretend that they use mpmath
print('\n\033[92m' + "Using " + yade.config.highPrecisionName + " with " + str(yade.config.highPrecisionDecimalPlaces) + " digits." + '\033[0m')
#print('\033[92m'+"Will not 'import mpmath' to reduce dependency on mpmath. Using a local minimal replacement module instead."+'\033[0m\n')
from math import *
import math
import cmath


class MP:
	dps = None


mp = MP
mpf = float


# simulate names provided by Eigen NumTraits
def eps():
	return 0


def cbrt(a):
	return pow(a, 1. / 3.)


power = pow
euler = 0.57721566490153286060651209008240243104
catalan = 0.9159655941772190150546035149323841107


# mpc needs to accept two string arguments, apart from that it's just a complex number
class mpc(complex):

	def __new__(cls, *args, **kwargs):
		super_new = super(mpc, cls).__new__
		if super_new is object.__new__:
			return super_new(cls)
		if len(args) == 2:
			return super_new(cls, float(args[0]), float(args[1]))
		return super_new(cls, *args, **kwargs)


# simulate all complex functions, just like they work in mpmath
# these overrides serve only one purpose: so that math functions can have the same name for complex and real numbers.
# They have the same names as in mpmath
def conj(a):
	return complex(a).conjugate()


def re(a):
	return complex(a).real


def im(a):
	return complex(a).imag


def norm(a):
	return abs(complex(a))  # // Warning: C++ std::norm is squared length. In python/mpmath norm is C++ abs == length.


def phase(a):
	return cmath.phase(a)  # tests std::arg


# skip squaredNorm, for    testing std::norm I will use norm(a)*norm(a) here.
# skip proj       , I skip testing std::proj. Later maybe add, see https://en.cppreference.com/w/cpp/numeric/complex/proj is checking for infinities in either component.
def rect(a, b):
	return cmath.rect(a, b)  # tests std::polar


def sin(a):
	return (cmath.sin(a) if (a.__class__ == mpc) else math.sin(a))


def sinh(a):
	return (cmath.sinh(a) if (a.__class__ == mpc) else math.sinh(a))


def cos(a):
	return (cmath.cos(a) if (a.__class__ == mpc) else math.cos(a))


def cosh(a):
	return (cmath.cosh(a) if (a.__class__ == mpc) else math.cosh(a))


def tan(a):
	return (cmath.tan(a) if (a.__class__ == mpc) else math.tan(a))


def tanh(a):
	return (cmath.tanh(a) if (a.__class__ == mpc) else math.tanh(a))


def asin(a):
	return (cmath.asin(a) if (a.__class__ == mpc) else math.asin(a))


def asinh(a):
	return (cmath.asinh(a) if (a.__class__ == mpc) else math.asinh(a))


def acos(a):
	return (cmath.acos(a) if (a.__class__ == mpc) else math.acos(a))


def acosh(a):
	return (cmath.acosh(a) if (a.__class__ == mpc) else math.acosh(a))


def atan(a):
	return (cmath.atan(a) if (a.__class__ == mpc) else math.atan(a))


def atanh(a):
	return (cmath.atanh(a) if (a.__class__ == mpc) else math.atanh(a))


def exp(a):
	return (cmath.exp(a) if (a.__class__ == mpc) else math.exp(a))


def log(a):
	return (cmath.log(a) if (a.__class__ == mpc) else math.log(a))


def log10(a):
	return (cmath.log10(a) if (a.__class__ == mpc) else math.log10(a))


# skip pow        , for    testing std::pow I will use a**b
def sqrt(a):
	return (cmath.sqrt(a) if (a.__class__ == mpc) else math.sqrt(a))
