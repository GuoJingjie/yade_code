/*************************************************************************
*  2019 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

/*
 * Add more functions as necessary, but remember to add them in:
 *

 lib/high-precision/MathSpecialFunctions.hpp  - add the function here

 py/high-precision/_math.cpp                  - export to python, under label:
                                                ********************************************************************************************
                                                **********************  Special Functions: bessel, factorial, etc ……  **********************
                                                ********************************************************************************************
                                                Remember: first export functions which take complex arguments.

 py/tests/testMath.py                         - test in python

 py/tests/testMathHelper.py                   - test helper. When mpmath isn't imported it provides functions with the same name exactly as in mpmath.
                                                if the function is not in mpmath, then it also is not in testMathHelper.
                                                It has to be calculated separately in a way which works for both mpmath and not mpmath.

 py/high-precision/_RealHPDiagnostics.cpp     - test on C++ side in checkSpecialFunctions()
                                                Then adjust tolerances in py/tests/testMath.py in getMathSpecialTolerance(…)

 */

#ifndef YADE_MATH_SPECIAL_FUNCIONS_HPP
#define YADE_MATH_SPECIAL_FUNCIONS_HPP

#include <lib/high-precision/Real.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

/*
 * covered here:
 *

 boost::math::cyl_bessel_j
 boost::math::factorial;
 boost::math::laguerre;
 boost::math::spherical_harmonic;

 * but use the YADE naming convention lowerUpperCase, e.g. cylBesselJ(…)
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
	/**********************  Special Functions: bessel, factorial, etc ……  **********************/
	/********************************************************************************************/
	// Add more functions as necessary, but remember to add them in:
	// * py/high-precision/_math.cpp
	// * py/tests/testMath.py
	// * py/tests/testMathHelper.py

	// https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/bessel/bessel_first.html
	template <typename Rr, int Level = levelOfRealHP<Rr>> inline Rr cylBesselJ(int k, const Rr& a)
	{
		using ::boost::math::cyl_bessel_j;
		return cyl_bessel_j(k, static_cast<const UnderlyingHP<Rr>&>(a));
	}

	// https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/factorials/sf_factorial.html
	template <typename Rr, int Level = levelOfRealHP<Rr>> inline Rr factorial(unsigned int a) // FIXME : templatyzacja po typie zwracanym
	{
		using ::boost::math::factorial;
		// careful:      ↓↓ calling this function always requires a templated return type.
		return factorial<Rr>(a);
	}

	// https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/sf_poly/laguerre.html
	template <typename Rr, int Level = levelOfRealHP<Rr>> inline Rr laguerre(unsigned int n, unsigned int m, const Rr& a)
	{
		using ::boost::math::laguerre;
		return laguerre(n, m, static_cast<const UnderlyingHP<Rr>&>(a));
	}

	// https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/math_toolkit/sf_poly/sph_harm.html
	template <typename Rr, int Level = levelOfRealHP<Rr>>
	inline ComplexHP<Level> sphericalHarmonic(unsigned int l, signed int m, const Rr& theta, const Rr& phi)
	{
		// TODO: report a bug to boost about forcing std::complex<cpp_bin_float<…>> instead of complex_adaptor<cpp_bin_float>. Same with MPFR mpc_complex_backend
		//       discussion started in https://github.com/BoostGSoC21/math/issues/16#issuecomment-907270124
		using ::boost::math::spherical_harmonic;
		// so I create tmp with improper complex type, then convert it.
		auto tmp = spherical_harmonic(l, m, static_cast<const UnderlyingHP<Rr>&>(theta), static_cast<const UnderlyingHP<Rr>&>(phi));
		return { tmp.real(), tmp.imag() };
	}

} // namespace math
} // namespace yade

#undef YADE_REAL_MATH_NAMESPACE
#endif
