/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file contains mathematical functions available in standard library and boost library.
//     https://en.cppreference.com/w/cpp/numeric/math
//     https://en.cppreference.com/w/cpp/numeric/special_functions
// They have to be provided here as inline redirections towards the correct implementation, depending on what precision type yade is being compiled with.
// This is the only way to make sure that ::std, ::boost::math, ::boost::multiprecision are all called correctly.

// TODO: Boost documentation recommends to link with tr1: -lboost_math_tr1 as it provides significant speedup. For example replace boost::math::acosh(x) ↔ boost::math::tr1::acosh(x)
//     https://www.boost.org/doc/libs/1_71_0/libs/math/doc/html/math_toolkit/overview_tr1.html
//#include <boost/math/tr1.hpp>

// All functions have to be tested in:
// * py/high-precision/_math.cpp
// * py/tests/testMath.py
// * py/tests/testMathHelper.py

#ifndef YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP
#define YADE_THIN_REAL_WRAPPER_MATH_COMPLEX_FUNCIONS_HPP

#include <lib/high-precision/Real.hpp>

namespace forCtags {
struct MathComplexFunctions {
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
	/**********************          Complex conj, real, imag, abs          *********************/
	/********************************************************************************************/
	// Add more complex functions as necessary, but remember to add them in py/high-precision/_math.cpp and py/tests/testMath.py
	// Note: most of the functions above can be converted to accept complex by changing levelOfRealHP<> → levelOfHP<>, provided that a complex version exists.
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade Complex supported HP types.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline Cc conj(const Cc& a)
	{
		using ::std::conj;
		using YADE_REAL_MATH_NAMESPACE::conj;
		return conj(static_cast<const UnderlyingHP<Cc>&>(a));
	}

	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> real(const Cc& a)
	{
		using ::std::real;
		using YADE_REAL_MATH_NAMESPACE::real;
		return real(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline RealHP<Level> imag(const Cc& a)
	{
		using ::std::imag;
		using YADE_REAL_MATH_NAMESPACE::imag;
		return imag(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, RealHP<Level>>::type abs(const Cc& a)
	{
		using ::std::abs;
		using YADE_REAL_MATH_NAMESPACE::abs;
		return abs(static_cast<const UnderlyingHP<Cc>&>(a));
	}

	/********************************************************************************************/
	/**********************     Real or Complex trigonometric functions    **********************/
	/********************************************************************************************/
	// typename RC is a type which can be Real or Complex, Rr → only Real, Cc → only Complex.
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade supported HP types.
	// int Level is the N in RealHP<N>.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sin(const Cc& a)
	{
		using ::std::sin;
		using YADE_REAL_MATH_NAMESPACE::sin;
		return sin(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sinh(const Cc& a)
	{
		using ::std::sinh;
		using YADE_REAL_MATH_NAMESPACE::sinh;
		return sinh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type cos(const Cc& a)
	{
		using ::std::cos;
		using YADE_REAL_MATH_NAMESPACE::cos;
		return cos(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type cosh(const Cc& a)
	{
		using ::std::cosh;
		using YADE_REAL_MATH_NAMESPACE::cosh;
		return cosh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type tan(const Cc& a)
	{
		using ::std::tan;
		using YADE_REAL_MATH_NAMESPACE::tan;
		return tan(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type tanh(const Cc& a)
	{
		using ::std::tanh;
		using YADE_REAL_MATH_NAMESPACE::tanh;
		return tanh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	// add more Real or Complex functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py

	/********************************************************************************************/
	/**********************      Real inverse trigonometric functions      **********************/
	/********************************************************************************************/
	// The check involving int Level = levelOfComplexHP<Cc> is necessary to make sure that function is called only with yade Real supported HP types.
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type asin(const Cc& a)
	{
		using ::std::asin;
		using YADE_REAL_MATH_NAMESPACE::asin;
		return asin(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type asinh(const Cc& a)
	{
		using ::std::asinh;
		using YADE_REAL_MATH_NAMESPACE::asinh;
		return asinh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type acos(const Cc& a)
	{
		using ::std::acos;
		using YADE_REAL_MATH_NAMESPACE::acos;
		return acos(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type acosh(const Cc& a)
	{
		using ::std::acosh;
		using YADE_REAL_MATH_NAMESPACE::acosh;
		return acosh(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type atan(const Cc& a)
	{
		using ::std::atan;
		using YADE_REAL_MATH_NAMESPACE::atan;
		return atan(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type atanh(const Cc& a)
	{
		using ::std::atanh;
		using YADE_REAL_MATH_NAMESPACE::atanh;
		return atanh(static_cast<const UnderlyingHP<Cc>&>(a));
	}

	/********************************************************************************************/
	/**********************   logarithm, exponential and power functions   **********************/
	/********************************************************************************************/
	// Add more functions as necessary, but remember to add them in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	// They can be converted to accept complex by changing levelOfRealHP<> → levelOfHP<>, provided that a complex version exists.
	// But remember to add tests for complex versions in py/high-precision/_math.cpp, py/tests/testMath.py and py/tests/testMathHelper.py
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type exp(const Cc& a)
	{
		using ::std::exp;
		using YADE_REAL_MATH_NAMESPACE::exp;
		return exp(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type log(const Cc& a)
	{
		using ::std::log;
		using YADE_REAL_MATH_NAMESPACE::log;
		return log(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type log10(const Cc& a)
	{
		using ::std::log10;
		using YADE_REAL_MATH_NAMESPACE::log10;
		return log10(static_cast<const UnderlyingHP<Cc>&>(a));
	}
	template <typename A, typename B, int Level = levelOfComplexHP<A>, typename Cc = PromoteHP<A>>
	inline typename boost::enable_if_c<(std::is_convertible<B, Cc>::value and isComplexHP<A>), Cc>::type pow(const A& a, const B& b)
	{
		using ::std::pow;
		using YADE_REAL_MATH_NAMESPACE::pow;
		return pow(static_cast<const UnderlyingHP<Cc>&>(a), static_cast<const UnderlyingHP<Cc>&>(b));
	}
	template <typename Cc, int Level = levelOfComplexHP<Cc>> inline typename boost::enable_if_c<isComplexHP<Cc>, Cc>::type sqrt(const Cc& a)
	{
		using ::std::sqrt;
		using YADE_REAL_MATH_NAMESPACE::sqrt;
		return sqrt(static_cast<const UnderlyingHP<Cc>&>(a));
	}


} // namespace math
} // namespace yade

#undef YADE_REAL_MATH_NAMESPACE
#endif
