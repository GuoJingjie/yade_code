/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once


#include <lib/compatibility/LokiCompatibility.hpp>
#include <lib/serialization/Serializable.hpp>

#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>

#include <string>

namespace yade { // Cannot have #include directive inside.


///////////////////////////////////////////////////////////////////////////////////////////////////
/// base template for classes that provide virtual functions for multiple dispatch,		///
/// in other words for multivirtual function call						///
///
///
/// This is a base template for all classes that will provide functions for multiple dispatch.
///
/// To create a new category of virtual functions you must derive from this class, inside you will
/// have provided virtual functions that you must overload to get multivirtual behaviour.
///
/// Depending on the number of dimensions for which you need multiple dispatch you have to overload
/// different virtual functions:
///
/// for 1D, you overload function [*]:
///		public: virtual ResultType go( ArgumentTypeList );
///
/// for 2D, you overload functions [*]:
///		public: virtual ResultType go( ArgumentTypeList );
///		public: virtual ResultType goReverse( ArgumentTypeList );
///
///	function goReverse is called only when DynLibDispatcher has autoSymmetry set to false.
///	otherwise DynLibDispatcher automatically reverses first two arguments, and calls go(),
///	same applies for 3D.
///
/// for 3D, you overload functions (not implemented now, but easy to write) [*]:
///		public: virtual ResultType go( ArgumentTypeList );
///		public: virtual ResultType go012( ArgumentTypeList );  // forwards call to go()
///		public: virtual ResultType go120( ArgumentTypeList );
///		public: virtual ResultType go201( ArgumentTypeList );
///		public: virtual ResultType go021( ArgumentTypeList );
///		public: virtual ResultType go210( ArgumentTypeList );
///
///
/// Template parameters:
///
///	ResultType 		- is the type returned by multivirtual function
///
///	ArgumentTypeList	- is a TypeList of arguments accepted by multivirtual function,
///				  ATTENTION:
///					- for 1D first type in this list must be of type shared_ptr<YourClass>
///					  this first argument acts like *this, in C++ virtual functions
///
///					- for 2D first and second type in this list must be shared_ptr<YourClass>
///					  those argument act like *this1 and *this2
///
///					- for 3D first, second and third type in this list must be shared_ptr<YourClass>
///					  those argument act like *this1 , *this2 and *this3
///
///
///
/// [*]
/// Note about virtual function arguents ArgumentTypeList - all functions take arguments by value
/// only for fundametal types and pure pointers, all other types are passed by referece. For details look
/// into Loki::TypeTraits::ParameterType. For example if you your class is:
///
/// class ShapeDraw : public FunctorWrapper< std::string , TYPELIST_4(boost::shared_ptr<Shape>,double,char,const std::string) >
/// {}
///
/// then virtual function to overload is:
///		public: virtual std::string go(boost::shared_ptr<Shape>&,double,char,const std::string& );
///
/// references were added where necessary, to optimize call speed.
/// So pay attention when you overload this function.
///

template <
        class ResultType,      // type returned by multivirtual function
        class ArgumentTypeList // TypeList of arguments accepted by multivirtual function,
        >
class FunctorWrapper //: public Serializable // FIXME functor shouldn't be serializable
{
private:
	using Parm1 = typename boost::mpl::at_c<ArgumentTypeList, 0>::type;
	using Parm2 = typename boost::mpl::at_c<ArgumentTypeList, 1>::type;
	using Parm3 = typename boost::mpl::at_c<ArgumentTypeList, 2>::type;
	using Parm4 = typename boost::mpl::at_c<ArgumentTypeList, 3>::type;
	using Parm5 = typename boost::mpl::at_c<ArgumentTypeList, 4>::type;
	using Parm6 = typename boost::mpl::at_c<ArgumentTypeList, 5>::type;
	using Parm7 = typename boost::mpl::at_c<ArgumentTypeList, 6>::type;

	ResultType error(int n)
	{
		throw std::runtime_error(
		        ("Multimethods: bad virtual call (probably go/goReverse was not overridden with the same argument types; only fundamental types and "
		         "pure pointers are passed by value, all other types (including shared_ptr<>) are passed by reference); types in the call were:\n"
		         + string("1. ") + typeid(Parm1).name() + "\n" + "2. " + typeid(Parm2).name() + "\n" + "3. " + typeid(Parm3).name() + "\n" + "4. "
		         + typeid(Parm4).name() + "\n" + "5. " + typeid(Parm5).name() + "\n" + "6. " + typeid(Parm6).name() + "\n" + "7. "
		         + typeid(Parm7).name() + "\n" + "number of types used in the call: " + boost::lexical_cast<string>(n) + "\n")
		                .c_str());
	}

public:
	FunctorWrapper() {};
	virtual ~FunctorWrapper() {};
	virtual string checkOrder() const { return ""; };

	// in following functions a second throw was added - just to bypass compiler warnings - it will never be executed.

	virtual ResultType go(Parm1) { return error(1); };
	virtual ResultType go(Parm1, Parm2) { return error(2); };
	virtual ResultType go(Parm1, Parm2, Parm3) { return error(3); };
	virtual ResultType go(Parm1, Parm2, Parm3, Parm4) { return error(4); };
	virtual ResultType go(Parm1, Parm2, Parm3, Parm4, Parm5) { return error(5); };
	virtual ResultType go(Parm1, Parm2, Parm3, Parm4, Parm5, Parm6) { return error(6); };
	virtual ResultType go(Parm1, Parm2, Parm3, Parm4, Parm5, Parm6, Parm7) { return error(7); };

	virtual ResultType goReverse(Parm1) { return error(1); };
	virtual ResultType goReverse(Parm1, Parm2) { return error(2); };
	virtual ResultType goReverse(Parm1, Parm2, Parm3) { return error(3); };
	virtual ResultType goReverse(Parm1, Parm2, Parm3, Parm4) { return error(4); };
	virtual ResultType goReverse(Parm1, Parm2, Parm3, Parm4, Parm5) { return error(5); };
	virtual ResultType goReverse(Parm1, Parm2, Parm3, Parm4, Parm5, Parm6) { return error(6); };
	virtual ResultType goReverse(Parm1, Parm2, Parm3, Parm4, Parm5, Parm6, Parm7) { return error(7); };
};

#define DEFINE_FUNCTOR_ORDER_2D(class1, class2)                                                                                                                \
public:                                                                                                                                                        \
	std::string checkOrder() const override { return (string(#class1) + " " + string(#class2)); }


}; // namespace yade
