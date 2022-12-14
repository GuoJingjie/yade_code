/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "DynLibManager.hpp"
#include <lib/base/Math.hpp>
#include <lib/base/Singleton.hpp>

/*! define the following macro to enable experimenta boost::serialization support
	Slows the compilation down about 2×.
	Python wrapper defines O.saveXML('file.xml') to try it out.
	Loading is not yet implemented (should be easy)
*/

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/serialization/export.hpp> // must come after all supported archive types

#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

namespace yade { // Cannot have #include directive inside.

#define REGISTER_FACTORABLE(name)                                                                                                                              \
	inline shared_ptr<Factorable> CreateShared##name() { return shared_ptr<name>(new name); }                                                              \
	inline Factorable*            Create##name() { return new name; }                                                                                      \
	inline void*                  CreatePureCustom##name() { return new name; }                                                                            \
	const bool                    registered##name __attribute__((unused))                                                                                 \
	        = ClassFactory::instance().registerFactorable(#name, Create##name, CreateShared##name, CreatePureCustom##name);


class Factorable;


/*! \brief The class factory of Yade used for serialization purpose and also as a dynamic library loader.
	All classes that call the macro REGISTER_FACTORABLE in their header are registered inside the factory
	so it is possible to ask the factory to create an instance of that class. This is automatic because the
	macro should be outside the class definition, so it is called automatically when the class is loaded by
	the program or when a dynamic library is loaded.
	
	This ClassFactory also acts as a dynamic library loader : when you ask for an instance, either the class
	already exists inside the factory and a new instance is created, or the class doesn't exist inside the
	factory and the ClassFactory will look on the hard drive to know if your class exists inside a dynamic library.
	If so the library is loaded and a new instance of the class can be created.

	\note ClassFactory is a singleton so you can't create an instance of it because its constructor is private.
	You should instead use ClassFactory::instance().createShared("Rigidbody") for example.
*/
class ClassFactory : public Singleton<ClassFactory> {
private:
	/// Pointer on a function that create an instance of a serializable class an return a shared pointer on it
	typedef shared_ptr<Factorable> (*CreateSharedFactorableFnPtr)();
	/// Pointer on a function that create an instance of a serializable class an return a C pointer on it
	typedef Factorable* (*CreateFactorableFnPtr)();
	/// Pointer on a function that create an instance of a custom class (i.e. not serializable) and return a void C pointer on it
	typedef void* (*CreatePureCustomFnPtr)();

	/// Description of a class that is stored inside the factory.
	struct FactorableCreators {
		CreateFactorableFnPtr       create;           /// Used to create a C pointer on the class (if serializable)
		CreateSharedFactorableFnPtr createShared;     /// Used to create a shared pointer on the class (if serializable)
		CreatePureCustomFnPtr       createPureCustom; /// Used to create a void C pointer on the class

		FactorableCreators() {};
		FactorableCreators(CreateFactorableFnPtr c, CreateSharedFactorableFnPtr cs, CreatePureCustomFnPtr cpc)
		{
			create           = c;
			createShared     = cs;
			createPureCustom = cpc;
		};
	};

	/// Type of a Stl map used to map the registered class name with their FactorableCreators
	typedef std::map<std::string, FactorableCreators> FactorableCreatorsMap;

	/// The internal dynamic library manager used to load dynamic libraries when an instance of a non loaded class is ask
	DynLibManager dlm;
	/// Map that contains the name of the registered class and their description
	FactorableCreatorsMap map;

	ClassFactory()
	{
		if (getenv("YADE_DEBUG")) fprintf(stderr, "Constructing ClassFactory.\n");
	}
	ClassFactory(const ClassFactory&);
	ClassFactory& operator=(const ClassFactory&);
	virtual ~ClassFactory() {};
	DECLARE_LOGGER;

public:
	/*! This method is used to register a Factorable class into the factory. It is called only from macro REGISTER_CLASS_TO_FACTORY
			\param name the name of the class
			\param create a pointer to a function that is able to return a C pointer on the given class
			\param createPureCustom a pointer to a function that is able to return a void C pointer on the given class
			\param verify a pointer to a function that is able to return the type_info of the given class
			\param type type of the class (SERIALIZABLE or CUSTOM)
			\param f is true is the class is a fundamental one (Vector3, Quaternion)
			\return true if registration is succesfull
		*/
	bool
	registerFactorable(std::string name, CreateFactorableFnPtr create, CreateSharedFactorableFnPtr createShared, CreatePureCustomFnPtr createPureCustom);

	/// Create a shared pointer on a serializable class of the given name
	shared_ptr<Factorable> createShared(std::string name);

	/// Create a C pointer on a serializable class of the given name
	Factorable* createPure(std::string name);

	/// Create a void C pointer on a class of the given name
	void* createPureCustom(std::string name);

	/*! Mainly used by the method findType for serialization purpose. Tells if a given type is a serilializable class
			\param tp type info of the type to test
			\param fundamental is true if the given type is fundamental (Vector3,Quaternion ...)
		*/

	bool        load(const string& fullFileName);
	std::string lastError();

	void              registerPluginClasses(const char* fileAndClasses[]);
	std::list<string> pluginClasses;

	virtual string getClassName() const { return "Factorable"; }; // FIXME - are they even inherited ????
	virtual string getBaseClassName(int) const { return ""; };

	FRIEND_SINGLETON(ClassFactory);
};


/*! Macro defining what classes can be found in this plugin -- must always be used in the respective .cpp file.
 * A function registerThisPluginClasses_FirstPluginName is generated at every occurence. The identifier should
 * be unique and avoids use of __COUNTER__ which didn't appear in gcc until 4.3.
 */

#define _YADE_PLUGIN_BOOST_REGISTER(x, y, z)                                                                                                                   \
	BOOST_CLASS_EXPORT_IMPLEMENT(yade::z);                                                                                                                 \
	BOOST_SERIALIZATION_FACTORY_0(yade::z);


#define _YADE_PLUGIN_REPEAT(x, y, z) BOOST_PP_STRINGIZE(z),

// If you get error "expected declaration before ‘}’ token" it means that you forgot put this macro inside yade namespace.
// To keep things consistent with REGISTER_SERIALIZABLE the YADE_PLUGIN is assumed to be inside yade namespace.
// However maco BOOST_SERIALIZATION_FACTORY_0 must be outside yade namespace, see also comment above REGISTER_SERIALIZABLE in lib/serialization/Serializable.hpp line 256
// https://github.com/boostorg/serialization/blob/develop/include/boost/serialization/factory.hpp#L91
// so here also to keep things inside single macro, there is closing and opening bracket for yade namespace.
#define YADE_PLUGIN(plugins)                                                                                                                                   \
	}                                                                                                                                                      \
	namespace {                                                                                                                                            \
		__attribute__((constructor)) void BOOST_PP_CAT(registerThisPluginClasses_, BOOST_PP_SEQ_HEAD(plugins))(void)                                   \
		{                                                                                                                                              \
			const char* info[] = { __FILE__, BOOST_PP_SEQ_FOR_EACH(_YADE_PLUGIN_REPEAT, ~, plugins) NULL };                                        \
			::yade::ClassFactory::instance().registerPluginClasses(info);                                                                          \
		}                                                                                                                                              \
	}                                                                                                                                                      \
	BOOST_PP_SEQ_FOR_EACH(_YADE_PLUGIN_BOOST_REGISTER, ~, plugins) namespace yade                                                                          \
	{
} // namespace yade
