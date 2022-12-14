/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Omega.hpp"
#include "Scene.hpp"
#include "ThreadRunner.hpp"
#include "TimeStepper.hpp"
#include <lib/base/Math.hpp>
#include <lib/multimethods/FunctorWrapper.hpp>
#include <lib/multimethods/Indexable.hpp>
#include <lib/serialization/ObjectIO.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/mutex.hpp>
#include <cxxabi.h>

SINGLETON_SELF(yade::Omega);

namespace yade { // Cannot have #include directive inside.

class RenderMutexLock : public std::lock_guard<std::mutex> {
public:
	RenderMutexLock()
	        : std::lock_guard<std::mutex>(Omega::instance().renderMutex)
	{
	}
	~RenderMutexLock() { }
};

CREATE_LOGGER(Omega);

const std::map<string, DynlibDescriptor>& Omega::getDynlibsDescriptor() const { return dynlibs; }

const shared_ptr<Scene>& Omega::getScene() const { return scenes.at(currentSceneNb); }
void                     Omega::resetCurrentScene()
{
	RenderMutexLock lock;
	scenes.at(currentSceneNb) = shared_ptr<Scene>(new Scene);
}
void Omega::resetScene() { resetCurrentScene(); }
void Omega::resetAllScenes()
{
	RenderMutexLock lock;
	scenes.resize(1);
	scenes[0]      = shared_ptr<Scene>(new Scene);
	currentSceneNb = 0;
}

int Omega::addScene()
{
	scenes.push_back(shared_ptr<Scene>(new Scene));
	return scenes.size() - 1;
}

void Omega::switchToScene(int i)
{
	if (i < 0 || i >= int(scenes.size())) {
		LOG_ERROR("Scene " << i << " has not been created yet, no switch.");
		return;
	}
	currentSceneNb = i;
}

Real Omega::getRealTime() const { return (boost::posix_time::microsec_clock::local_time() - startupLocalTime).total_milliseconds() / 1e3; }

boost::posix_time::time_duration Omega::getRealTime_duration() const { return boost::posix_time::microsec_clock::local_time() - startupLocalTime; }

void Omega::initTemps()
{
	char dirTemplate[] = "/tmp/yade-XXXXXX";
	tmpFileDir         = mkdtemp(dirTemplate);
	tmpFileCounter     = 0;
}

void Omega::cleanupTemps()
{
	stop();
	boost::filesystem::path tmpPath(tmpFileDir);
	boost::filesystem::remove_all(tmpPath);
}

std::string Omega::tmpFilename()
{
	if (tmpFileDir.empty()) throw runtime_error("tmpFileDir empty; Omega::initTemps not yet called()?");
	const std::lock_guard<std::mutex> lock(tmpFileCounterMutex);
	return tmpFileDir + "/tmp-" + boost::lexical_cast<string>(tmpFileCounter++);
}

void Omega::reset()
{
	stop();
	init();
}

void Omega::init()
{
	sceneFile = "";
	resetAllScenes();
	sceneAnother = shared_ptr<Scene>(new Scene);
	timeInit();
	createSimulationLoop();
}

void Omega::timeInit() { startupLocalTime = boost::posix_time::microsec_clock::local_time(); }

void Omega::createSimulationLoop() { simulationLoop = shared_ptr<ThreadRunner>(new ThreadRunner(&simulationFlow_)); }
void Omega::stop()
{
	LOG_DEBUG("");
	if (simulationLoop && simulationLoop->looping()) simulationLoop->stop();
	if (simulationLoop)
		// TODO: throwing from here causes segfault. But sometimes calling the destructor ~ThreadRunner() cannot be done - also causes segfault
		//       if(simulationLoop->permissionToDestroy()) ....
		simulationLoop = shared_ptr<ThreadRunner>();
}

/* WARNING: even a single simulation step is run asynchronously; the call will return before the iteration is finished. */
void Omega::step()
{
	if (simulationLoop) { simulationLoop->spawnSingleAction(); }
}

void Omega::run()
{
	if (!simulationLoop) {
		LOG_ERROR("No Omega::simulationLoop? Creating one (please report bug).");
		createSimulationLoop();
	}
	if (simulationLoop && !simulationLoop->looping()) { simulationLoop->start(); }
}

void Omega::pause()
{
	if (simulationLoop && simulationLoop->looping()) { simulationLoop->stop(); }
}

bool Omega::isRunning()
{
	if (simulationLoop) return simulationLoop->looping();
	else
		return false;
}

void Omega::buildDynlibDatabase(const vector<string>& dynlibsList)
{
	LOG_DEBUG("called with " << dynlibsList.size() << " plugins.");
	boost::python::object wrapperScope = boost::python::import("yade.wrapper");
	std::list<string>     pythonables;
	for (const auto& name : dynlibsList) {
		shared_ptr<Factorable> f;
		try {
			LOG_DEBUG("Factoring plugin " << name);
			f                            = ClassFactory::instance().createShared(name);
			dynlibs[name].isSerializable = ((YADE_PTR_DYN_CAST<Serializable>(f)).get() != 0);
			for (int i = 0; i < f->getBaseClassNumber(); i++) {
				dynlibs[name].baseClasses.insert(f->getBaseClassName(i));
			}
			if (dynlibs[name].isSerializable) pythonables.push_back(name);
		} catch (std::runtime_error& e) {
			/* FIXME: this catches all errors! Some of them are not harmful, however:
			 * when a class is not factorable, it is OK to skip it; */
		}
	}
	// handle Serializable specially
	//Serializable().pyRegisterClass(wrapperScope);
	/* python classes must be registered such that base classes come before derived ones;
	for now, just loop until we succeed; proper solution will be to build graphs of classes
	and traverse it from the top. It will be done once all classes are pythonable. */
	int numTries = 100;
	for (int i = 0; i <= numTries && pythonables.size() > 0; i++) {
		if (getenv("YADE_DEBUG")) cerr << endl << "[[[ Round " << i << " ]]]: ";
		for (std::list<string>::iterator I = pythonables.begin(); I != pythonables.end();) {
			shared_ptr<Serializable> s = boost::static_pointer_cast<Serializable>(ClassFactory::instance().createShared(*I));
			try {
				if (getenv("YADE_DEBUG")) cerr << "{{" << *I << "}}";
				s->pyRegisterClass(wrapperScope);
				std::list<string>::iterator prev = I++;
				pythonables.erase(prev);
			} catch (...) {
				if (getenv("YADE_DEBUG")) cerr << "[" << *I << "]";
				if (i == numTries)
					PyErr_Print(); // we want to see the actual error if it still fails after 100th attempt, else we hide useful errors
				boost::python::handle_exception();
				PyErr_Clear();
				I++;
			}
		}
	}

	std::map<string, DynlibDescriptor>::iterator dli    = dynlibs.begin();
	std::map<string, DynlibDescriptor>::iterator dliEnd = dynlibs.end();
	for (; dli != dliEnd; ++dli) {
		std::set<string>::iterator bci    = (*dli).second.baseClasses.begin();
		std::set<string>::iterator bciEnd = (*dli).second.baseClasses.end();
		for (; bci != bciEnd; ++bci) {
			string name = *bci;
			if (name == "Dispatcher1D" || name == "Dispatcher2D") (*dli).second.baseClasses.insert("Dispatcher");
			else if (name == "Functor1D" || name == "Functor2D")
				(*dli).second.baseClasses.insert("Functor");
			else if (name == "Serializable")
				(*dli).second.baseClasses.insert("Factorable");
			else if (name != "Factorable" && name != "Indexable") {
				shared_ptr<Factorable> f = ClassFactory::instance().createShared(name);
				for (int i = 0; i < f->getBaseClassNumber(); i++)
					dynlibs[name].baseClasses.insert(f->getBaseClassName(i));
			}
		}
	}
}

bool Omega::isInheritingFrom(const string& className, const string& baseClassName)
{
	return (dynlibs[className].baseClasses.find(baseClassName) != dynlibs[className].baseClasses.end());
}

bool Omega::isInheritingFrom_recursive(const string& className, const string& baseClassName)
{
	if (dynlibs[className].baseClasses.find(baseClassName) != dynlibs[className].baseClasses.end()) return true;
	for (const auto& parent : dynlibs[className].baseClasses) {
		if (isInheritingFrom_recursive(parent, baseClassName)) return true;
	}
	return false;
}

void Omega::loadPlugins(vector<string> pluginFiles)
{
	for (const auto& plugin : pluginFiles) {
		LOG_DEBUG("Loading plugin " << plugin);
		if (!ClassFactory::instance().load(plugin)) {
			string err = ClassFactory::instance().lastError();
			if (err.find(": undefined symbol: ") != std::string::npos) {
				size_t pos = err.rfind(":");
				assert(pos != std::string::npos);
				std::string sym(err, pos + 2); //2 removes ": " from the beginning
				int         status        = 0;
				char*       demangled_sym = abi::__cxa_demangle(sym.c_str(), 0, 0, &status);
				LOG_FATAL(plugin << ": undefined symbol `" << demangled_sym << "'");
				LOG_FATAL(plugin << ": " << err);
				LOG_FATAL("Bailing out.");
			} else {
				LOG_FATAL(plugin << ": " << err << " ."); /* leave space to not to confuse c++filt */
				LOG_FATAL("Bailing out.");
			}
			abort();
		}
	}
	std::list<string>& plugins(ClassFactory::instance().pluginClasses);
	plugins.sort();
	plugins.unique();
	buildDynlibDatabase(vector<string>(plugins.begin(), plugins.end()));
}

void Omega::loadSimulation(const string& f, bool quiet)
{
	bool isMem = boost::algorithm::starts_with(f, ":memory:");
	if (!isMem && !boost::filesystem::exists(f)) throw runtime_error("Simulation file to load doesn't exist: " + f);
	if (isMem && memSavedSimulations.count(f) == 0) throw runtime_error("Cannot load nonexistent memory-saved simulation " + f);

	if (!quiet) LOG_INFO("Loading file " + f);
	shared_ptr<Scene>& scene = scenes[currentSceneNb];
	{
		stop(); // stop current simulation if running
		resetScene();
		RenderMutexLock lock;
		if (isMem) {
			std::istringstream iss(memSavedSimulations[f]);
			yade::ObjectIO::load<decltype(scene), boost::archive::binary_iarchive>(iss, "scene", scene);
		} else {
			yade::ObjectIO::load(f, "scene", scene);
		}
	}
	if (scene->getClassName() != "Scene") throw logic_error("Wrong file format (scene is not a Scene!?) in " + f);
	sceneFile = f;
	timeInit();

	//ForceContainer is not serialized, better prepare it with correct size here
	//Assumption: maxId is size-1
	scene->forces.addMaxId(scene->bodies->size() - 1);

	if (!quiet) LOG_DEBUG("Simulation loaded");
}

void Omega::saveSimulation(const string& f, bool quiet)
{
	if (f.size() == 0) throw runtime_error("f of file to save has zero length.");
	if (!quiet) LOG_INFO("Saving file " << f);
	shared_ptr<Scene>& scene = scenes[currentSceneNb];
	if (boost::algorithm::starts_with(f, ":memory:")) {
		if (memSavedSimulations.count(f) > 0 && !quiet) LOG_INFO("Overwriting in-memory saved simulation " << f);
		std::ostringstream oss;
		yade::ObjectIO::save<decltype(scene), boost::archive::binary_oarchive>(oss, "scene", scene);
		memSavedSimulations[f] = oss.str();
	} else {
		yade::ObjectIO::save(f, "scene", scene);
	}
	sceneFile = f;
}

} // namespace yade
