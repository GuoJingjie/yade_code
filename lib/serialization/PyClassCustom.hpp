/*************************************************************************
*  2021 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef PY_CLASS_CUSTOM_HPP
#define PY_CLASS_CUSTOM_HPP
#if YADE_REAL_BIT <= 64
// When Real == double, then there is nothing to do, we can use native boost::python::class_
// It has to stay inside namespace boost::python to avoid name lookup problems.
// This might change later, if someone will need to .def_readonly(…) for RealHP<2> member variables.
// Declaring them inside the YADE_CLASS_BASE_* macros works normally, without any issues.
namespace boost {
namespace python {
	template <typename... AllArguments> using PyClassCustom = ::boost::python::class_<AllArguments...>;
}
}
// Note: We have not (yet?) templatized entire Serialization.hpp file.
// (1) lots of work
// (2) it's possible that compilation will be slower after that
#else
/*
 * So this file is to solve this kind of runtime error:
 *
 *   TypeError: TypeError: No Python class registered for C++ class boost::multiprecision::number<………>
 *
 * The problem is that python tries to return internal reference to the Real type in def_readonly(…) calls, but is not possible,
 * because there is no python class corresponding to that type. The solution is to ::boost::python::return_by_value.
 * If this error occurs in more places, then it means that the solution for this error has to be improved here
 * (or to create an even larger wrapper for Real, like minieigenHP).
 *
 * See also:
 *
 *   https://yade-dem.org/doc/prog.html#custom-converters
 *   https://www.boost.org/doc/libs/1_68_0/libs/python/doc/html/faq/why_is_my_automatic_to_python_co.html
 *
 * So I will go the short route: override the def_readonly(…) function,
 * instead of modifying dozens of *.hpp files which call def_readonly(…)
 *
 */
//#include <boost/core/demangle.hpp>
//#include <iostream>

namespace boost {
namespace python {
	template <typename... AllArguments> struct PyClassCustom : public ::boost::python::class_<AllArguments...> {
		using Parent = ::boost::python::class_<AllArguments...>;
		template <typename... Args>
		PyClassCustom(Args... args) // use perfect forwarding on constructor arguments
		        : Parent(::std::forward<Args>(args)...)
		{
		}

		// safely detect if that is a high precision type.
		template <typename T, bool = (not std::is_function<T>::value)> struct NotClassMemberFunction {
			static const constexpr bool isHP = false;
		};
		template <typename T> struct NotClassMemberFunction<T, true> {    // true: means that this is not a function, but a member variable
			static const constexpr bool isHP = ::yade::math::isHP<T>; // so we can check isHP
		};
		template <typename T> struct ClassMember {
			static const constexpr bool isHP = false;
		};
		template <typename Class, typename TypeOfMember> struct ClassMember<TypeOfMember Class::*> {
			static const constexpr bool isHP = NotClassMemberFunction<TypeOfMember>::isHP;
		};
		//template <typename T, bool = (boost::is_floating_point<T>::value) and (not std::is_member_function_pointer<T>::value) and (not std::is_function<T>::value) and (not std::is_member_pointer<T>::value)> struct NotClassMemberFunction

		/********************************************************************************************/
		/**********************    The solution of the def_readonly problem    **********************/
		/********************************************************************************************/
		template <typename D, typename ::boost::enable_if_c<ClassMember<D>::isHP, int>::type = 0>
		auto def_readonly(char const* name, D d, char const* doc = 0)
		{
			// This whole file is to write only this single line with ::boost::python::return_by_value
			Parent::add_property(
			        name, ::boost::python::make_getter(d, ::boost::python::return_value_policy<::boost::python::return_by_value>()), doc);
			return *this;
		}
		template <typename D, typename ::boost::disable_if_c<ClassMember<D>::isHP, int>::type = 0>
		auto def_readonly(char const* name, D d, char const* doc = 0)
		{
			//std::cerr << "\nDEMANGLE 1 " << boost::core::demangle(typeid(D).name()) << "  \t\t name: " << name << "\n";
			Parent::def_readonly(name, d, doc);
			return *this;
		}

		/********************************************************************************************/
		/**********************   Make sure that chain calls will be working   **********************/
		/********************************************************************************************/
		// Redeclare all public member functions which return 'self&' in /usr/include/boost/python/class.hpp, this way we make sure that the chain calls will work.
		// This is necessary because static template overload cannot have virtual member functions (which is good btw, because we want maximum speed here).
		// These functions are optimized away by the compiler.
		template <typename... Args> auto def(Args... args)
		{
			Parent::def(::std::forward<Args>(args)...);
			return *this;
		}
		template <typename... Args> auto def_readwrite(Args... args)
		{
			Parent::def_readwrite(::std::forward<Args>(args)...);
			return *this;
		}
		template <typename... Args> auto add_property(Args... args)
		{
			Parent::add_property(::std::forward<Args>(args)...);
			return *this;
		}
		template <typename... Args> auto add_static_property(Args... args)
		{
			Parent::add_static_property(::std::forward<Args>(args)...);
			return *this;
		}
		template <typename... Args> auto setattr(Args... args)
		{
			Parent::setattr(::std::forward<Args>(args)...);
			return *this;
		}
		template <typename... Args> auto def_pickle(Args... args)
		{
			Parent::def_pickle(::std::forward<Args>(args)...);
			return *this;
		}
		auto enable_pickling()
		{
			Parent::enable_pickling();
			return *this;
		}
		auto staticmethod(char const* name)
		{
			Parent::staticmethod(name);
			return *this;
		}
	};
}
}
#endif
#endif
