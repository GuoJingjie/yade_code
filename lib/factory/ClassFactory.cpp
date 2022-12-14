/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "ClassFactory.hpp"

#include <lib/base/Logging.hpp>
#include <boost/algorithm/string/regex.hpp>

SINGLETON_SELF(yade::ClassFactory);

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(ClassFactory);

class Factorable;

bool ClassFactory::registerFactorable(
        std::string name, CreateFactorableFnPtr create, CreateSharedFactorableFnPtr createShared2, CreatePureCustomFnPtr createPureCustom2)
{
	bool tmp = map.insert(FactorableCreatorsMap::value_type(name, FactorableCreators(create, createShared2, createPureCustom2))).second;
	return tmp;
}

shared_ptr<Factorable> ClassFactory::createShared(std::string name)
{
	FactorableCreatorsMap::const_iterator i = map.find(name);
	if (i == map.end()) {
		dlm.load(name);
		if (dlm.isLoaded(name)) {
			if (map.find(name) == map.end()) {
				// Well, exception are also a way to return value, right?
				// Vaclav, this comment is so wrong I don't know where to start. / Janek.
				// This throws at startup for every .so that doesn't contain class named the same as the library.
				// I.e. almost everything in yade-libs and some others installed locally...
				// Don't log that, it would confuse users.
				//LOG_FATAL("Can't create class "<<name<<" - was never registered.");
				throw std::runtime_error(("Class " + name + " not registered in the ClassFactory.").c_str());
			}
			return createShared(name);
		}
		throw std::runtime_error(("Class " + name + " could not be factored in the ClassFactory.").c_str());
	}
	return (i->second.createShared)();
}

Factorable* ClassFactory::createPure(std::string name)
{
	FactorableCreatorsMap::const_iterator i = map.find(name);
	if (i == map.end()) {
		dlm.load(name);
		if (dlm.isLoaded(name)) {
			if (map.find(name) == map.end()) { throw std::runtime_error(("Class " + name + " not registered in the ClassFactory.").c_str()); }
			return createPure(name);
		}
		throw std::runtime_error(("Class " + name + " could not be factored in the ClassFactory.").c_str());
	}
	return (i->second.create)();
}

void* ClassFactory::createPureCustom(std::string name)
{
	FactorableCreatorsMap::const_iterator i = map.find(name);
	if (i == map.end()) throw std::runtime_error(("Class " + name + " could not be factored in the ClassFactory.").c_str());
	return (i->second.createPureCustom)();
}

bool ClassFactory::load(const string& name) { return dlm.load(name); }

string ClassFactory::lastError() { return dlm.lastError(); }

void ClassFactory::registerPluginClasses(const char* fileAndClasses[])
{
	assert(fileAndClasses[0] != NULL); // must be file name
	// only filename given, no classes names explicitly
	if (fileAndClasses[1] == NULL) {
		/* strip leading path (if any; using / as path separator) and strip one suffix (if any) to get the contained class name */
		string heldClass = boost::algorithm::replace_regex_copy(string(fileAndClasses[0]), boost::regex("^(.*/)?(.*?)(\\.[^.]*)?$"), string("\\2"));
#ifdef YADE_DEBUG
		if (getenv("YADE_DEBUG"))
			cerr << __FILE__ << ":" << __LINE__ << ": Plugin " << fileAndClasses[0] << ", class " << heldClass << " (deduced)" << endl;
#endif
		pluginClasses.push_back(heldClass); // last item with everything up to last / take off and .suffix strip
	} else {
		for (int i = 1; fileAndClasses[i] != NULL; i++) {
#ifdef YADE_DEBUG
			if (getenv("YADE_DEBUG"))
				cerr << __FILE__ << ":" << __LINE__ << ": Plugin " << fileAndClasses[0] << ", class " << fileAndClasses[i] << endl;
#endif
			pluginClasses.push_back(fileAndClasses[i]);
		}
	}
}

} // namespace yade
