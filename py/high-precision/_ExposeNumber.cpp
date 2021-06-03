/*************************************************************************
*  2021      Janek Kozicki                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This file is to workaround problem with exposing readonly Real via def_readonly(…), it is an alternative solution to
// commit 61df54d60862d0ee81f610f4e2275d78b53b6bf1

#include <lib/high-precision/RealIO.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
using namespace ::yade::MathEigenTypes;
namespace py = boost::python;

template <typename Number> class NumberVisitor : public py::def_visitor<NumberVisitor<Number>> {
	using Scalar     = Number;                              // could be complex number
	using RealScalar = typename yade::math::RealOf<Number>; // this is the "real" (math) scalar
public:
	template <class PyClass> void visit(PyClass& cl) const
	{
		cl.def(py::init<Number>(py::arg("other")))
		        .def("__str__", &NumberVisitor::__str__)
		        .def("__repr__", &NumberVisitor::__str__)
		        .def("__abs__", &NumberVisitor::__abs__)
		        .def("__neg__", &NumberVisitor::__neg__)
		        .def("__eq__", &NumberVisitor::__eq__)
		        .def("__ne__", &NumberVisitor::__ne__)
		        .def("__add__", &NumberVisitor::__add__)
		        .def("__iadd__", &NumberVisitor::__iadd__)
		        .def("__sub__", &NumberVisitor::__sub__)
		        .def("__isub__", &NumberVisitor::__isub__)
		        .def("__mul__", &NumberVisitor::__mul__scalar<long>)
		        .def("__imul__", &NumberVisitor::__imul__scalar<long>)
		        .def("__rmul__", &NumberVisitor::__rmul__scalar<long>);
		visit_if_float<Scalar, PyClass>(cl);
		visit_if_real<Scalar, PyClass>(cl);
		visit_if_complex<Scalar, RealScalar, PyClass>(cl);
	};

private:
	template <class PyClass> static std::string name(PyClass& cl) { return py::extract<std::string>(cl.attr("__name__"))(); }
	template <typename Scalar, class PyClass, typename boost::enable_if<boost::is_integral<Scalar>, int>::type = 0> static void visit_if_float(PyClass&)
	{ // do nothing
	}
	template <typename Scalar, class PyClass, typename boost::disable_if<boost::is_integral<Scalar>, int>::type = 0> static void visit_if_float(PyClass& cl)
	{
		cl.def("__mul__", &NumberVisitor::__mul__scalar<Scalar>)
		        .def("__rmul__", &NumberVisitor::__rmul__scalar<Scalar>)
		        .def("__imul__", &NumberVisitor::__imul__scalar<Scalar>)
		        .def("__div__", &NumberVisitor::__div__scalar<long>)
		        .def("__rtruediv__", &NumberVisitor::__rdiv__scalar<long>)
		        .def("__truediv__", &NumberVisitor::__div__scalar<long>)
		        .def("__idiv__", &NumberVisitor::__idiv__scalar<long>)
		        .def("__itruediv__", &NumberVisitor::__div__scalar<long>)
		        .def("__div__", &NumberVisitor::__div__scalar<Scalar>)
		        .def("__truediv__", &NumberVisitor::__div__scalar<Scalar>)
		        .def("__idiv__", &NumberVisitor::__idiv__scalar<Scalar>)
		        .def("__itruediv__", &NumberVisitor::__idiv__scalar<Scalar>);
	}
	template <typename Scalar, class PyClass, typename boost::disable_if_c<::yade::math::isRealHP<Scalar>, int>::type = 0>
	static void visit_if_real(PyClass&)
	{ // do nothing
	}
	template <typename Scalar, class PyClass, typename boost::enable_if_c<::yade::math::isRealHP<Scalar>, int>::type = 0>
	static void visit_if_real(PyClass& cl)
	{
		cl.def("__round__", &NumberVisitor::__round__scalar)
		        .def("__round__", &NumberVisitor::__round__places)
		        .def("__lt__", &NumberVisitor::__lt__)
		        .def("__le__", &NumberVisitor::__le__)
		        .def("__gt__", &NumberVisitor::__gt__)
		        .def("__ge__", &NumberVisitor::__ge__)
// not this way:        .def("isnan", &NumberVisitor::isnan)
/*
t=TriaxialStressController()
isnan(t.porosity)
TypeError: must be real number, not RealHP1

In [4]: yade.math.isnan(t.porosity)
Out[4]: False # OK.
So we have a working isnan, it needs to be imported in proper place.

mpmath.mpf(t.porosity)
TypeError: cannot create mpf from 1
So now mpmath compatibility must be added.
*/
			;
	}
	template <typename CNum, typename RScalar, class PyClass, typename boost::disable_if_c<::yade::math::isComplexHP<CNum>, int>::type = 0>
	static void visit_if_complex(PyClass&)
	{ // do nothing
	}
	template <typename CNum, typename RScalar, class PyClass, typename boost::enable_if_c<::yade::math::isComplexHP<CNum>, int>::type = 0>
	static void visit_if_complex(PyClass& cl)
	{
		cl.def("__rmul__", &NumberVisitor::__rmul__scalar<RScalar>)
		        .def("__imul__", &NumberVisitor::__imul__scalar<RScalar>)
		        .def("__div__", &NumberVisitor::__div__scalar<RScalar>)
		        .def("__truediv__", &NumberVisitor::__div__scalar<RScalar>)
		        .def("__idiv__", &NumberVisitor::__idiv__scalar<RScalar>)
		        .def("__itruediv__", &NumberVisitor::__idiv__scalar<RScalar>);
		;
	}
	static std::string __str__(const Number& a) { return ::yade::math::toStringHP(a); };
	static Number      __abs__(const Number& a) { return ::yade::math::abs(a); };
	static Number      __neg__(const Number& a) { return -a; };
	static Number      __round__scalar(const Number& a) { return ::yade::math::round(a); };
	static Number      __round__places(const Number& a, const long& places) { return ::yade::math::round(a * places) / static_cast<Number>(places); };
	static bool        __lt__(const Number& a, const Number& b) { return a < b; }
	static bool        __le__(const Number& a, const Number& b) { return a <= b; }
	static bool        __gt__(const Number& a, const Number& b) { return a > b; }
	static bool        __ge__(const Number& a, const Number& b) { return a >= b; }
	static bool        __eq__(const Number& a, const Number& b) { return a == b; }
	static bool        __ne__(const Number& a, const Number& b) { return !__eq__(a, b); }
//	static bool        isnan(const Number& a) { return ::yade::math::isnan(a); };

	static Number __add__(const Number& a, const Number& b) { return a + b; }
	static Number __iadd__(Number& a, const Number& b)
	{
		a += b;
		return a;
	};
	static Number __sub__(const Number& a, const Number& b) { return a - b; }
	static Number __isub__(Number& a, const Number& b)
	{
		a -= b;
		return a;
	};
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __mul__scalar(const Number& a, const Scalar2& scalar)
	{
		return a * scalar;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __mul__scalar(const Number& a, const Scalar2& scalar)
	{
		return a * static_cast<Scalar>(scalar);
	}
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __imul__scalar(Number& a, const Scalar2& scalar)
	{
		a *= scalar;
		return a;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __imul__scalar(Number& a, const Scalar2& scalar)
	{
		a *= static_cast<Scalar>(scalar);
		return a;
	}
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __rmul__scalar(const Number& a, const Scalar2& scalar)
	{
		return a * scalar;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __rmul__scalar(const Number& a, const Scalar2& scalar)
	{
		return a * static_cast<Scalar>(scalar);
	}
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __rdiv__scalar(const Number& a, const Scalar2& scalar)
	{
		return scalar / a;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __rdiv__scalar(const Number& a, const Scalar2& scalar)
	{
		return static_cast<Scalar>(scalar) / a;
	}
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __div__scalar(const Number& a, const Scalar2& scalar)
	{
		return a / scalar;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __div__scalar(const Number& a, const Scalar2& scalar)
	{
		return a / static_cast<Scalar>(scalar);
	}
	template <typename Scalar2, typename boost::enable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __idiv__scalar(Number& a, const Scalar2& scalar)
	{
		a /= scalar;
		return a;
	}
	template <typename Scalar2, typename boost::disable_if<std::is_convertible<Scalar2, const Scalar&>, int>::type = 0>
	static Number __idiv__scalar(Number& a, const Scalar2& scalar)
	{
		a /= static_cast<Scalar>(scalar);
		return a;
	}
};

template <int N, bool = (std::numeric_limits<RealHP<N>>::digits10 >= 18)> struct ExposeNum {
	//                                                         ↑ the lower Real types are covered by python by default.
	static void work(bool, const py::scope&) {};
};

template <int N> struct ExposeNum<N, true> {
	static void work(bool notDuplicate, const py::scope& topScope)
	{
		std::string numStr     = boost::lexical_cast<std::string>(N);
		std::string HPn        = "HP" + numStr;
		std::string realHPn    = "Real" + HPn;
		std::string complexHPn = "Complex" + HPn;
		if (notDuplicate) {
			// currently only for readonly access. The write access is currently through python mpmath package and ToFromPythonConverter
			py::class_<RealHP<N>>(realHPn.c_str(), ("The Real<" + numStr + "> type.").c_str(), py::init<>()).def(NumberVisitor<RealHP<N>>());
			py::class_<ComplexHP<N>>(complexHPn.c_str(), ("The Complex<" + numStr + "> type.").c_str(), py::init<>())
			        .def(NumberVisitor<ComplexHP<N>>());
		} else {
			py::scope().attr(realHPn.c_str())    = topScope.attr(realHPn.c_str());
			py::scope().attr(complexHPn.c_str()) = topScope.attr(complexHPn.c_str());
		}
	}
};

template <int N> void expose_number(bool notDuplicate, const py::scope& topScope) { ExposeNum<N>::work(notDuplicate, topScope); }

// explicit instantination - tell compiler to produce a compiled version of expose_converters (it is faster when done in parallel in .cpp files)
YADE_HP_PYTHON_REGISTER(expose_number)

