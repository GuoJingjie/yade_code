// 2008 © Václav Šmilauer <eudoxos@arcig.cz>
#pragma once
#include <lib/pyutil/gil.hpp>
#include <core/GlobalEngine.hpp>
#include <core/Scene.hpp>
#include <pkg/common/PeriodicEngines.hpp>

namespace yade { // Cannot have #include directive inside.

class PyRunner : public PeriodicEngine {
public:
	/* virtual bool isActivated: not overridden, PeriodicEngine handles that */
	void action() override
	{
		if (command.size() > 0) pyRunString(command, ignoreErrors, updateGlobals);
	}
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS(PyRunner,PeriodicEngine,
		"Execute a python command periodically, with defined (and adjustable) periodicity. See :yref:`PeriodicEngine` documentation for details.",
		((string,command,"",,"Command to be run by python interpreter. Not run if empty."))
		((bool,ignoreErrors,false,,"Debug only: set this value to true to tell PyRunner to ignore any errors encountered during command execution."))
		((bool,updateGlobals,true,,R"""(
Whether to workaround `ipython not recognizing local variables <https://github.com/ipython/ipython/issues/62>`__
by calling ``globals().update(locals())``. If ``true`` then PyRunner is able to call functions declared later locally in a running **live** yade session.
The ``PyRunner`` call is a bit slower because it updates ``globals()`` with recently declared python functions.

.. warning::
	When ``updateGlobals==False`` and a function was declared inside a *live* yade session (`ipython <http://ipython.org>`_)
	then an error ``NameError: name 'command' is not defined`` will occur unless python ``globals()`` are updated with command

	.. code-block:: python

		globals().update(locals())

)"""))
	);
	// clang-format on
};
REGISTER_SERIALIZABLE(PyRunner);

} // namespace yade
