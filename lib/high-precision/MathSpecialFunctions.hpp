/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// NOTE: add more functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.
// This is the only way to make sure that ::std, ::boost::math, ::boost::multiprecision are all called correctly.

// TODO: Boost documentation recommends to link with tr1: -lboost_math_tr1 as it provides significant speedup. For example replace boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//     https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
//#include <boost/math/tr1.hpp>

#ifndef YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP

#include <lib/high-precision/Real.hpp>
#include <boost/config.hpp>
#include <boost/math/complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <limits>
#include <utility>

/*
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
 boost::math::cyl_bessel_j
 boost::math::factorial;
 boost::math::laguerre;
 boost::math::spherical_harmonic;
*/

namespace forCtags {
struct MathComplexFunctions {
}; // for ctags
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
	// use SFINAE to allow other non RealHP type, its level is treated as == 1
	//template <typename HP, typename Allow, typename boost::enable_if_c<(isComplexHP<HP> or std::is_same<HP, Allow>::value), int>::type = 0>
	//const constexpr int levelOfComplexHPAllow = ((isHP<HP>) ? (levelOrZero<HP>) : (1));

	/********************************************************************************************/
	/**********************   logarithm, exponential and power functions   **********************/
	/********************************************************************************************/
	// Add more functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	// They can be converted to accept complex by changing levelOfRealHP<> → levelOfHP<>, provided that a complex version exists.
	// But remember to add tests for complex versions in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	template <typename A, typename B, int Level = levelOfComplexHP<A>, typename Cc = PromoteHP<A>>
	inline typename boost::enable_if_c<(std::is_convertible<B, Cc>::value and isComplexHP<A::type>), Cc>::type pow(const A& a, const B& b)
	{
		using ::std::pow;
		using YADE_REAL_MATH_NAMESPACE::pow;
		return pow(static_cast<const UnderlyingHP<Cc>&>(a), static_cast<const UnderlyingHP<Cc>&>(b));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>>
	inline typename boost::enable_if_c<isComplexHP<typename std::decay<Cc>::type>, Cc>::type sqrt(const Cc& a)
	{
		using ::std::sqrt;
		using YADE_REAL_MATH_NAMESPACE::sqrt;
		return sqrt(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	// allow here Complex also
	//template <typename RC, int Level = levelOfHP<RC>> inline RealHP<Level> abs(const RC& a);

} // namespace math
} // namespace yade

#undef YADE_REAL_MATH_NAMESPACE
#endif
