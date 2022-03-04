/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// C++ backend of the ymport module.

#include <boost/python.hpp>

#include <lib/base/AliasNamespaces.hpp>
#include <lib/base/LoggingUtils.hpp>
#include <lib/pyutil/doc_opts.hpp>

#include <memory>

#include "ymport/foamfile/PolyMesh.hpp"

CREATE_CPP_LOCAL_LOGGER("_ymport.cpp");

namespace yade {
namespace ymport {
	namespace foamfile {

		py::list readPolyMesh(const std::string& directoryPath, bool patchAsWall, bool emptyAsWall)
		{
			std::unique_ptr<PolyMesh> polyMesh = std::make_unique<PolyMesh>(directoryPath, patchAsWall, emptyAsWall);

			return polyMesh->facetCoords();
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade

BOOST_PYTHON_MODULE(_ymport)
try {
	using namespace yade::ymport;
	namespace py = ::boost::python;
	YADE_SET_DOCSTRING_OPTS;
	py::def("readPolyMesh", foamfile::readPolyMesh, R"""(
	"""C++ backend of :yref:`yade.ymport.polyMesh`.

	:param str path: directory path. Typical value is: "constant/polyMesh".
	:param bool patchAsWall: load "patch"-es as walls.
	:param bool emptyAsWall: load "empty"-es as walls.
	:param \*\*kw: (unused keyword arguments) is passed to :yref:`yade.utils.facet`
	:returns: list of facets.
    )""");
} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}
