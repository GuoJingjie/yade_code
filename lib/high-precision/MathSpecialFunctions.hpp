/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// Add more functions as necessary, but remember to add them in:
// * py/high-precision/_math.cpp
// * py/tests/testMath.py
// * py/tests/testMathHelper.py

#ifndef YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP

#include <lib/high-precision/Real.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

/*
 * cover here:
 *

 boost::math::cyl_bessel_j
 boost::math::factorial;
 boost::math::laguerre;
 boost::math::spherical_harmonic;
*/

namespace forCtags {
struct MathSpecialFunctions {
}; // for ctags
}

/*************************************************************************/
/*********************** YADE_REAL_MATH_NAMESPACE ************************/
/*************************************************************************/
// The YADE_REAL_MATH_NAMESPACE macro makes sure that proper namespace is used for Real and eventual RealHP<…> types.
// On older systems RealHP<…> can't work, so we need to check if it was explicitly disabled during cmake autodetection.
#if ((YADE_REAL_BIT <= 80) and defined(YADE_DISABLE_REAL_MULTI_HP))
#define YADE_REAL_MATH_NAMESPACE ::std
#else
#define YADE_REAL_MATH_NAMESPACE ::boost::multiprecision
#endif

namespace yade {
namespace math {

	/********************************************************************************************/
	/**********************        bessel, factorial, laguerre, etc        **********************/
	/********************************************************************************************/
	// Add more functions as necessary, but remember to add them in:
	// * py/high-precision/_math.cpp
	// * py/tests/testMath.py
	// * py/tests/testMathHelper.py

} // namespace math
} // namespace yade

#undef YADE_REAL_MATH_NAMESPACE
#endif
