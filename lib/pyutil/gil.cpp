#include <lib/base/Logging.hpp>
#include <lib/pyutil/gil.hpp>
#include <boost/python.hpp>

CREATE_CPP_LOCAL_LOGGER(gil_cpp);
void pyRunString(const std::string& cmd, bool ignoreErrors, bool updateGlobals)
{
	namespace py = ::boost::python;
	gilLock    lock;
	py::object main    = py::import("__main__");
	py::object globals = main.attr("__dict__");
	py::scope  setScope(main);
	try {                                                         // https://docs.python.org/3/c-api/reflection.html#c.PyEval_GetFrame
		if (updateGlobals and PyEval_GetFrame() != nullptr) { // don't bother if there is no frame present. For example inside yade --check
			py::object ipython = py::import("IPython");   // access the running IPython session.
			// FIXED: when a new function is declared inside ipython session an extra python command must be called: globals().update(locals())
			//        https://stackoverflow.com/questions/43956636/internal-function-call-not-working-in-ipython
			globals.attr("update")(ipython.attr("__dict__")["__builtins__"]["locals"]());
		}
		py::exec(cmd.c_str(), globals);
	} catch (const py::error_already_set&) {
		// using approach similar to https://stackoverflow.com/questions/1418015/how-to-get-python-exception-text
		py::handle_exception();
		// prepare values of sys.last_type, sys.last_value, sys.last_traceback, see https://docs.python.org/3/c-api/exceptions.html#c.PyErr_PrintEx
		PyErr_Print();
		py::exec("import traceback, sys", globals);
		auto        err   = py::eval("str(sys.last_value)", globals);
		auto        trace = py::eval(R"""('\n'.join(traceback.format_exception(sys.last_type, sys.last_value, sys.last_traceback)))""", globals);
		std::string msg   = "PyRunner error.\n\nCOMMAND: '" + cmd;
		msg += "'\n\nERROR:\n" + py::extract<std::string>(err)() + "\n\nSTACK TRACE:\n" + py::extract<std::string>(trace)();
		if (ignoreErrors) {
			LOG_WARN(msg << "\nbut has ignoreErrors == true; not throwing exception.");
		} else {
			throw std::runtime_error(msg);
		}
	} catch (...) {
		py::handle_exception();
		if (ignoreErrors) {
			LOG_WARN("Error running command: '" << cmd << "', but has ignoreErrors == true; not throwing exception.");
		} else {
			LOG_ERROR("Error running command: '" << cmd << "'");
			throw std::runtime_error("PyRunner → pyRunString encountered error while executing command: '" + cmd + "'");
		}
	}
};
