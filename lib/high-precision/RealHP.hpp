/*************************************************************************
*  2020 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef YADE_REAL_HP_TYPE_HPP
#define YADE_REAL_HP_TYPE_HPP

namespace forCtags {
struct RealHP {
}; // for ctags, but also see RealHPEigenCgal
struct ComplexHP {
};
}

// In this file following things are declared for the users:
//
// (1) RealHP<N>, ComplexHP<N> types are declared in such a way that:
//     RealHP<1> == Real
//     RealHP<2> == has twice the precision of Real
//     RealHP<3> == has three times the precision of Real (this one is mainly for compatibility with the counting, usually this type is "slow")
//     RealHP<4> == has four  times the precision of Real (when https://github.com/boostorg/multiprecision/issues/184 is done it will be fast)
//     ...
//     RealHP<N> == has N     times the precision of Real (all the higher types are based on eiher MPFR or boost::cpp_bin_float)
//
// (2) constexpr bool isHP<Type>               - has value true, when Type belongs to RealHP<N> family
// (3) constexpr int  levelOfRealHP<Type>      - has value N in RealHP<N>   , provided that 'Type' is a valid RealHP<N>    type, otherwise it causes a compilation error.
// (4) constexpr int  levelOfComplexHP<Type>   - has value N in ComplexHP<N>, provided that 'Type' is a valid ComplexHP<N> type, otherwise it causes a compilation error.
// (5) constexpr int  levelOfHP<N>             - has value N for RealHP or ComplexHP

// note: int N in the RealHP<N> has type 'signed int N', so that compiler can catch errors when -1 is passed as argument.
// Older compilers could convert -1 to a large positive number which results in runtime segfault, instead of compilation error.

#ifndef YADE_REAL_MATH_NAMESPACE
#error "This file cannot be included alone, include Real.hpp instead"
#endif

#include <boost/mpl/at.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/utility/enable_if.hpp>
#ifdef YADE_MPFR
#include <boost/multiprecision/mpfr.hpp>
#if (BOOST_VERSION >= 107100)
#include <boost/multiprecision/mpc.hpp>
#endif
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Warray-bounds"
#include <boost/multiprecision/cpp_bin_float.hpp>
#pragma GCC diagnostic pop
#if (BOOST_VERSION >= 107100)
#include <boost/multiprecision/cpp_complex.hpp>
#endif
#endif
#if ((YADE_REAL_BIT <= 32) and (not defined(YADE_DISABLE_REAL_MULTI_HP)))
#include "ThinComplexWrapper.hpp"
#include "ThinRealWrapper.hpp"
// after the above two:
#include "NumericLimits.hpp"
#endif
#if (defined(YADE_THIN_REAL_WRAPPER_HPP) and defined(YADE_THIN_COMPLEX_WRAPPER_HPP))
namespace yade {
namespace math {
	template <typename A> const constexpr bool IsWrapped = std::is_same<ThinRealWrapper<long double>, typename std::decay<A>::type>::value;
	template <typename A> using MakeWrappedComplexHP     = ThinComplexWrapper<A>;
}
}
#else
namespace yade {
namespace math { // store info that ThinRealWrapper is not used.
	template <typename A> const constexpr bool IsWrapped = false;
	template <typename A> using MakeWrappedComplexHP     = A;
}
}
#endif
// clang currently does not support float128   https://github.com/boostorg/math/issues/181
// another similar include is in ThinRealWrapper.hpp, all other checks if we have float128 should be #ifdef BOOST_MP_FLOAT128_HPP or yade::math::isFloat128<T>
#if (((((YADE_REAL_BIT <= 64) and (not defined(YADE_DISABLE_REAL_MULTI_HP)))) or (YADE_REAL_BIT == 128)) and (not defined(YADE_FLOAT128_UNAVAILABLE)))
#include <boost/multiprecision/float128.hpp>
#if (BOOST_VERSION >= 107100)
#include <boost/multiprecision/complex128.hpp>
#endif
#endif

namespace yade {
namespace math {
	/*************************************************************************/
	/************************* MPFR / cpp_bin_float **************************/
	/*************************************************************************/
	// Here is declared the 'ultimate' precision type called NthLevelRealHP. Depending on settings it can be MPFR or cpp_bin_float.
	namespace detail {
		// Two approaches to complex numbers. Their skills:
		//
		// 1. ENABLE_COMPLEX_MP=ON (default)
		//    - uses best complex type available from boost::multiprecision , the (1) complex128 and (2) mpc_complex<…> and (3) complex_adaptor<cpp_bin_float<…>>
		//    - UpconversionOfBasicOperatorsHP.hpp            - with MPFR doesn't work for Real{+,-,*,/}Complex. The mpc_complex too greedily provides operator+-*/ which then we cannot overload.
		//                                                      this has been reported: https://github.com/boostorg/multiprecision/issues/363
		//    - MathComplexFunctions.hpp                      - always uses this better type, math functions have smaller error
		//    - doesn't work with boost < 1.71
		//
		// 2. ENABLE_COMPLEX_MP=OFF (optional, but default when boost < 1.71)
		//    - uses std::complex<for everything>             - generally not recommended because C++ standard guarantees std::complex to work only for fundamental types.
		//    - UpconversionOfBasicOperatorsHP.hpp            - works without any problems.
		//    - MathComplexFunctions.hpp                      - larger numerical error in the mathematical functions
		//    - always works on older compilers, linux distros
		//
		// Generally ENABLE_COMPLEX_MP=ON is numerically better and has only one drawback: not working UpconversionOfBasicOperatorsHP.hpp. It is the default selection.
		//
		// UpconversionOfBasicOperatorsHP.hpp works without problems only with ENABLE_COMPLEX_MP=OFF this is optional via cmake.
		// We closely watch if std::complex<User Defined Types> work. If they stop we will act accordingly.
		// There is a chance this problem will get some standardized solution: https://github.com/BoostGSoC21/math/issues/16#issuecomment-907625419
		// And then we are ready to switch to the standard solution.
		//
		template <typename T> struct MakeComplex {
			using type = std::complex<T>;
		};
#if (defined(YADE_THIN_REAL_WRAPPER_HPP) and defined(YADE_THIN_COMPLEX_WRAPPER_HPP))
		template <> struct MakeComplex<ThinRealWrapper<long double>> {
			using type = ThinComplexWrapper<std::complex<long double>>;
		};
#endif
#if ((BOOST_VERSION >= 107100) and defined(YADE_COMPLEX_MP))
// ENABLE_COMPLEX_MP=ON : The ComplexHP type is extended with (1) complex128 and (2) mpc_complex_backend<…> and (3) complex_adaptor<cpp_bin_float<…>>
#if (((((YADE_REAL_BIT <= 64) and (not defined(YADE_DISABLE_REAL_MULTI_HP)))) or (YADE_REAL_BIT == 128)) and (not defined(YADE_FLOAT128_UNAVAILABLE)))
		template <> struct MakeComplex<boost::multiprecision::float128> {
			using type = boost::multiprecision::complex128;
		};
#endif
		using ::boost::multiprecision::number;
#ifdef YADE_MPFR
		// here MakeComplex allows creation of mpc_complex_backend types which use MPFR
		using ::boost::multiprecision::expression_template_option;
		using ::boost::multiprecision::mpc_complex_backend;
		using ::boost::multiprecision::mpfr_allocation_type;
		using ::boost::multiprecision::backends::mpfr_float_backend;
		template <unsigned Digits, mpfr_allocation_type AllocationType, expression_template_option ExpressionTemplates>
		struct MakeComplex<number<mpfr_float_backend<Digits, AllocationType>, ExpressionTemplates>> {
			using type = number<mpc_complex_backend<Digits>, ExpressionTemplates>;
		};
#else
		// here MakeComplex allows creation of complex_adaptor types which use cpp_bin_float
		using ex_opt = ::boost::multiprecision::expression_template_option;
		using ::boost::multiprecision::backends::complex_adaptor;
		using ::boost::multiprecision::backends::cpp_bin_float;
		using ::boost::multiprecision::backends::digit_base_type;
		template <unsigned Digits, digit_base_type DigitBase, class Alloc, class Exp, Exp MinEx, Exp MaxEx, ex_opt ExpressionTemplates>
		struct MakeComplex<number<cpp_bin_float<Digits, DigitBase, Alloc, Exp, MinEx, MaxEx>, ExpressionTemplates>> {
			using type = number<complex_adaptor<cpp_bin_float<Digits, DigitBase, Alloc, Exp, MinEx, MaxEx>>, ExpressionTemplates>;
		};
#endif
#endif
#ifdef YADE_MPFR
		// if MPFR is available, then use it for higher types. YADE_MPFR is defined when yade links with MPFR after compilation. It is unrelated about how the Real type is defined by YADE_REAL_MPFR.
		template <int DecPlaces> using RealHPBackend = boost::multiprecision::mpfr_float_backend<DecPlaces>;
#else
		// otherwise use boost::cpp_bin_float
		template <int DecPlaces> using RealHPBackend = boost::multiprecision::cpp_bin_float<DecPlaces>;
#endif
		// the NthLevelRealHP uses MPFR or boost::cpp_bin_float to declare type with requested precision
		template <int Level>
		using NthLevelRealHP = boost::multiprecision::number<RealHPBackend<std::numeric_limits<Real>::digits10 * Level>, boost::multiprecision::et_off>;
	}

	/*************************************************************************/
	/*************************     RealHPLadder     **************************/
	/*************************************************************************/
	// Here is defined the progressing ladder of types 'RealHPLadder', each next type has higher precision, than the previous type.
	// depending on precision specified in compilation options, the types to use are a bit different. They "shift upwards" within RealHPLadder

#ifndef YADE_DISABLE_REAL_MULTI_HP
#if YADE_REAL_BIT <= 32
	// Real == float
	using RealHPLadder = boost::mpl::vector<Real, double, ThinRealWrapper<long double>, detail::NthLevelRealHP<4>, boost::multiprecision::float128>;
#elif YADE_REAL_BIT <= 64
	// Real == double
	// later quad-double will be added to this list, https://github.com/boostorg/multiprecision/issues/184, so that RealHP<2> and RealHP<4> will be quite fast.
#ifdef BOOST_MP_FLOAT128_HPP
	using RealHPLadder = boost::mpl::vector<Real, boost::multiprecision::float128>;
#else // clang does not support float128   https://github.com/boostorg/math/issues/181
	using RealHPLadder = boost::mpl::vector<Real>; // when float128 is not available the RealHPLadder does not include it.
#endif
#elif ((YADE_REAL_BIT <= 80) or defined(YADE_NON_386_LONG_DOUBLE))
	// Real == long double
	// Here it might become interesting in few years when templatized versions of doubling-* and quadding-* of types are implemented in boost::multiprecision.
	// Then it will be high speed, high precision. Because long double although not standarized is pretty fast, and then if it's doubled or quadded it might be good and fast.
	using RealHPLadder = boost::mpl::vector<Real>;
#elif YADE_REAL_BIT <= 128
	// Real == float128
	// later quad-double will be added here
	using RealHPLadder = boost::mpl::vector<Real>;
#elif (defined(YADE_REAL_MPFR) or defined(YADE_REAL_BBFLOAT))
	// Real == cpp_bin_float or MPFR
	using RealHPLadder = boost::mpl::vector<Real>;
#endif
#else
	using RealHPLadder = boost::mpl::vector<Real>;
#endif

	/*************************************************************************/
	/******   Find position on the RealHPLadder, find Underlying type   ******/
	/*************************************************************************/
	// first declare helper templates needed to find the position in RealHPLadder
	namespace detail {
		// posRealHP is position of HP in the RealHPLadder. If it's found then it equals to findPosRealHP + 1, otherwise it equals to 0.
		const constexpr int                        ladderSize    = boost::mpl::size<RealHPLadder>::value;
		template <typename HP> const constexpr int findPosRealHP = boost::mpl::find<RealHPLadder, HP>::type::pos::value;
		template <typename HP> const constexpr int posRealHP     = (findPosRealHP<HP> == ladderSize) ? (0) : (findPosRealHP<HP> + 1);

		// calculates the Level 'N' for types higher than those in RealHPLadder, they use MPFR or cpp_bin_float
		template <typename HP>
		const constexpr auto NthLevel = std::numeric_limits<typename std::decay<HP>::type>::digits10 / std::numeric_limits<Real>::digits10;
		// checks if HP is an MPFR or cpp_bin_float
		template <typename HP> const constexpr bool isNthLevel = std::is_same<NthLevelRealHP<NthLevel<HP>>, typename std::decay<HP>::type>::value;
	} // namespace detail

	/*************************************************************************/
	/*************************      levelOrZero     **************************/
	/*************************************************************************/
	// The level of RealHP. Equals to 0 if HP is not from RealHP<N> family
	template <typename HP> const constexpr int levelOrZero = (detail::isNthLevel<HP>) ? (detail::NthLevel<HP>) : (detail::posRealHP<HP>);

	/*************************************************************************/
	/*************************     Underlying*HP    **************************/
	/*************************************************************************/
	// Extract UnderlyingReal (strip ThinRealWrapper) from type 'HP' which may be any of the RealHP<N> types.
	// This implementation might be revised in the future if ThinRealWrapper or UnderlyingReal will acquire some new meaning. But interface shall remain the same.
	template <typename HP> using UnderlyingRealHP    = typename std::conditional<IsWrapped<HP>, long double, HP>::type;
	template <typename HP> using UnderlyingComplexHP = typename detail::MakeComplex<UnderlyingRealHP<HP>>::type;

	/*************************************************************************/
	/*************************      AND FINALLY     **************************/
	/*************************        declare       **************************/
	/*************************       RealHP<N>      **************************/
	/*************************************************************************/
	// OK. At this point we have the typelist 'RealHPLadder' which lists the consecutive types in the progressing precision ladder. Now we can use it to define RealHP<int> type.
	// Also we have some helpers: 'Underlying*HP<N>' and 'levelOrZero'

	// Note about Arb or flint library, and MPFI interval calculations or Automatic differentiation techniques:
	// they can be implemented here as a type trait for the Real type. Meaning that a Tag class can be used as a trait in definition of Real to select
	// which of these (regular one, MPFI, Arb, ……) is to be used.

#ifndef YADE_DISABLE_REAL_MULTI_HP
	template <int Level>
	using RealHP = typename std::conditional<
	        Level <= boost::mpl::size<RealHPLadder>::value,           // if Level <= than max size of RealHPLadder
	        typename boost::mpl::at_c<RealHPLadder, Level - 1>::type, // then use the type from the typelist
	        detail::NthLevelRealHP<Level>>::type;                     // otherwise use the type constructed from MPFR or from boost::cpp_bin_float

	template <int Level>
	using ComplexHP = typename std::conditional<
	        IsWrapped<RealHP<Level>>,                                 // construct ComplexHP type
	        MakeWrappedComplexHP<UnderlyingComplexHP<RealHP<Level>>>, // depending on whether it is wrapped or not
	        UnderlyingComplexHP<RealHP<Level>>>::type;

#else
	// RealHP<…> won't work on this system, cmake sets YADE_DISABLE_REAL_MULTI_HP to use RealHP<1> for all precisions RealHP<N>.
	template <int Level> using RealHP                     = Real;                      // ignore Level
	template <int Level> using ComplexHP                  = UnderlyingComplexHP<Real>; // ignore Level
#endif

	/*************************************************************************/
	/*************************   inspection of HP   **************************/
	/*************************************************************************/
	// boost::is_complex only checks for std::is_complex (confirmed by talking with boost devs). We need to check this ourselves.
	// Like in https://github.com/BoostGSoC21/math/blob/develop/include/boost/math/fft/multiprecision_complex.hpp
	//
	// So this detail is to properly recognize various complex types with struct IsHPComplex
	namespace detail {
		// Older compilers do not have std::void_t, this small helper type is used by HasPlus, HasMinus, etc, below.
		template <typename... Ts> struct make_void {
			typedef void type;
		};
		template <typename... Ts> using void_type = typename make_void<Ts...>::type;

		template <class, class = void> struct has_value_type : std::false_type {
		};
		template <class T> struct has_value_type<T, void_type<typename T::value_type>> : std::true_type {
		};
		template <typename T, bool = has_value_type<T>::value> struct IsHPComplex {
			static constexpr bool value = false;
		};
		template <typename T> struct IsHPComplex<T, true> {
			static constexpr bool value = std::is_same<typename MakeComplex<typename T::value_type>::type, typename std::decay<T>::type>::value;
		};
	} // namespace detail

	// Now that 'RealHP<N>' and 'ComplexHP<N>' are defined we can use the SFINAE techniques to define 'levelOfRealHP', 'levelOfComplexHP' and 'levelOfHP'
	// And make sure that they are  ⇒→ undefined ←⇐  when the argument is not from the RealHP<N> family.

	template <typename HP, bool = detail::IsHPComplex<HP>::value> struct InspectHP {
		using RT                         = typename std::decay<HP>::type;
		using Under                      = UnderlyingRealHP<RT>;
		static const constexpr bool isHP = (levelOrZero<HP> != 0 and std::is_same<RealHP<levelOrZero<HP>>, typename std::decay<HP>::type>::value);
	};
	template <typename HP> struct InspectHP<HP, true> {
		using RT                         = typename std::decay<typename HP::value_type>::type;
		using Under                      = UnderlyingComplexHP<RT>;
		static const constexpr bool isHP = (levelOrZero<RT> != 0 and std::is_same<ComplexHP<levelOrZero<RT>>, typename std::decay<HP>::type>::value);
	};

	template <typename HP> using RealOf                                                                               = typename InspectHP<HP>::RT;
	template <typename HP> using UnderlyingHP                                                                         = typename InspectHP<HP>::Under;
	template <typename HP> const constexpr bool                                                             isHP      = InspectHP<HP>::isHP;
	template <typename HP> const constexpr bool                                                             isReal    = not detail::IsHPComplex<HP>::value;
	template <typename HP> const constexpr bool                                                             isComplex = not isReal<HP>;
	template <typename HP> const constexpr bool                                                             isRealHP  = (isReal<HP> and isHP<HP>);
	template <typename HP> const constexpr bool                                                             isComplexHP      = (isComplex<HP> and isHP<HP>);
	template <typename HP, typename boost::enable_if_c<isHP<HP>, int>::type = 0> const constexpr int        levelOfHP        = levelOrZero<RealOf<HP>>;
	template <typename HP, typename boost::enable_if_c<isRealHP<HP>, int>::type = 0> const constexpr int    levelOfRealHP    = levelOfHP<HP>;
	template <typename HP, typename boost::enable_if_c<isComplexHP<HP>, int>::type = 0> const constexpr int levelOfComplexHP = levelOfHP<HP>;
	// check if it's float128
#ifdef BOOST_MP_FLOAT128_HPP
	template <typename A> const constexpr bool isFloat128 = std::is_same<boost::multiprecision::float128, typename std::decay<A>::type>::value;
#else
	template <typename A> const constexpr bool isFloat128 = false;
#endif
	// check if HP is an MPFR or cpp_bin_float
	template <typename HP> const constexpr bool isNthLevel = detail::isNthLevel<HP>;

	/*************************************************************************/
	/*************************  SelectHigherHP<A,B> **************************/
	/*************************************************************************/
	// SelectHigherHP helps selecting a maximum type encompassing both A and B.

	template <typename A, typename B, int LevelA = levelOfHP<A>, int LevelB = levelOfHP<B>, int MaxLevel = ((LevelA > LevelB) ? (LevelA) : (LevelB))>
	using SelectHigherHP = typename std::conditional<(isComplex<A> or isComplex<B>), ComplexHP<MaxLevel>, RealHP<MaxLevel>>::type;

	/*************************************************************************/
	/*************************   allow/exclude HP   **************************/
	/*************************************************************************/
	// sometimes MathFunctions.hpp needs to accept/exclude specific type arguments.

	// use SFINAE to allow other non HP type, its level is treated as == 1
	template <typename HP, typename Allow, typename boost::enable_if_c<(isHP<HP> or std::is_same<HP, Allow>::value), int>::type = 0>
	const constexpr int levelOfHPAllow = ((isHP<HP>) ? (levelOrZero<HP>) : (1));

	// use SFINAE to allow other non RealHP type, its level is treated as == 1
	template <typename HP, typename Allow, typename boost::enable_if_c<(isRealHP<HP> or std::is_same<HP, Allow>::value), int>::type = 0>
	const constexpr int levelOfRealHPAllow = ((isHP<HP>) ? (levelOrZero<HP>) : (1));

	// use SFINAE to filter out a specific type. Used by min/max in MathFunctions.hpp.
	template <typename HP, typename Except, typename boost::enable_if_c<(isRealHP<HP> and (not std::is_same<HP, Except>::value)), int>::type = 0>
	const constexpr int levelOfRealHPExcept = levelOfHP<HP>;

	// PromoteHP turns any non-HP type into RealHP<1>, otherwise it just keeps it as it is, be it Real or Complex.
	template <typename HP, int Level = levelOrZero<RealOf<HP>>> using PromoteHP = typename std::conditional<(Level == 0), RealHP<1>, HP>::type;

	/*************************************************************************/
	/*************************       Complex        **************************/
	/*************************************************************************/
	template <typename Rr>
	const constexpr bool isConstexprFloatingPoint = (std::numeric_limits<Rr>::digits10 <=
#ifdef BOOST_MP_FLOAT128_HPP
	                                                 33
#else
	                                                 15
#endif
	                                                 )
	        and
	        // 24 decimal places are used for RealHP<4> which cannot use constexpr when compiling yade with float, see lib/high-precision/RealHP.hpp:103
	        (std::numeric_limits<Rr>::digits10 != 24);

	using Complex = ComplexHP<1>;
} // namespace math

using math::ComplexHP;
using math::RealHP;

using Real    = math::Real;
using Complex = math::Complex;

static_assert(sizeof(Real) == sizeof(math::UnderlyingReal), "This compiler introduced padding, which breaks binary compatibility");
static_assert(
        (sizeof(Complex) == sizeof(std::complex<math::UnderlyingReal>))
                or ((not math::isConstexprFloatingPoint<math::UnderlyingReal>)and(
                        std::is_same<Complex, typename math::detail::MakeComplex<math::UnderlyingReal>::type>::value)),
        "This compiler introduced padding, which breaks binary compatibility");

// If necessary for debugging put here the static_assert tests from the end of py/high-precision/_math.cpp file

} // namespace yade

#endif
