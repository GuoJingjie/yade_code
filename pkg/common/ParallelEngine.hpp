#pragma once
#include <core/GlobalEngine.hpp>

namespace yade { // Cannot have #include directive inside.

class ParallelEngine;
shared_ptr<ParallelEngine> ParallelEngine_ctor_list(const boost::python::list& slaves);

class ParallelEngine : public Engine {
public:
	typedef vector<vector<shared_ptr<Engine>>> slaveContainer;
	void                                       action() override;
	bool                                       isActivated() override { return true; }
	// py access
	boost::python::list slaves_get();
	void                slaves_set(const boost::python::list& slaves);
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(ParallelEngine,Engine,"Engine for running other Engine in parallel.",
		((slaveContainer,slaves,,,"[will be overridden]"))
		,
		/*ctor*/
		ompThreads=2;
		,
		/*py*/
		.def("__init__",boost::python::make_constructor(ParallelEngine_ctor_list),"Construct from (possibly nested) list of slaves.")
		.add_property("slaves",&ParallelEngine::slaves_get,&ParallelEngine::slaves_set,"List of lists of Engines; each top-level group will be run in parallel with other groups, while Engines inside each group will be run sequentially, in given order.");
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(ParallelEngine);

} // namespace yade
