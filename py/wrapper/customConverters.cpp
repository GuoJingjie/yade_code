/*************************************************************************
*  2009 Václav Šmilauer                                                  *
*  2010 Bruno Chareyre                                                   *
*  2010 Chiara Modenese                                                  *
*  2014 Anton Gladky                                                     *
*  2014 Jan Stránský                                                     *
*  2019 François Kneib                                                   *
*  2019 Janek Kozicki                                                    *
*  2021 Jérôme Duriez                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// for some reading about these conversions, see e.g. https://misspent.wordpress.com/2009/09/27/how-to-write-boost-python-converters/ (except the latter advocates using some "borrowed()" in the from-Python construc(), which is not done here)

#include <lib/base/Logging.hpp>
#include <lib/base/Math.hpp>
#include <lib/base/openmp-accu.hpp>
#include <lib/high-precision/ToFromPythonConverter.hpp>
#include <core/Callbacks.hpp>
#include <core/Dispatching.hpp>
#include <core/Engine.hpp>
#include <pkg/common/KinematicEngines.hpp>
#include <pkg/common/MatchMaker.hpp>
#include <preprocessing/dem/SpherePack.hpp>
#ifdef YADE_OPENGL
#include <pkg/common/GLDrawFunctors.hpp>
#include <pkg/common/OpenGLRenderer.hpp>
#endif

namespace forCtags {
struct customConverters {
}; // for ctags
}

CREATE_CPP_LOCAL_LOGGER("customConverters.cpp");

namespace yade { // Cannot have #include directive inside.

// move this to the miniEigen wrapper later

/* two-way se3 handling */
struct custom_se3_to_tuple {
	static PyObject* convert(const Se3r& se3)
	{
		boost::python::tuple ret = boost::python::make_tuple(se3.position, se3.orientation);
		return boost::python::incref(ret.ptr());
	}
};
struct custom_Se3r_from_seq {
	custom_Se3r_from_seq() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<Se3r>()); }
	static void* convertible(PyObject* obj_ptr)
	{
		if (!PySequence_Check(obj_ptr)) return 0;
		if (PySequence_Size(obj_ptr) != 2 && PySequence_Size(obj_ptr) != 7) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<Se3r>*)(data))->storage.bytes;
		new (storage) Se3r;
		Se3r* se3 = (Se3r*)storage;
		if (PySequence_Size(obj_ptr) == 2) { // from vector and quaternion
			se3->position    = boost::python::extract<Vector3r>(PySequence_GetItem(obj_ptr, 0));
			se3->orientation = boost::python::extract<Quaternionr>(PySequence_GetItem(obj_ptr, 1));
		} else if (PySequence_Size(obj_ptr) == 7) { // 3 vector components, 3 axis components, angle
			se3->position = Vector3r(
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 0)),
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 1)),
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 2)));
			Vector3r axis = Vector3r(
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 3)),
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 4)),
			        boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 5)));
			Real angle       = boost::python::extract<Real>(PySequence_GetItem(obj_ptr, 6));
			se3->orientation = Quaternionr(AngleAxisr(angle, axis));
		} else
			throw std::logic_error(__FILE__
			                       ": First, the sequence size for Se3r object was 2 or 7, but now is not? (programming error, please report!");
		data->convertible = storage;
	}
};


struct custom_OpenMPAccumulator_to_float {
#if (YADE_REAL_BIT <= 64)
	static PyObject* convert(const OpenMPAccumulator<Real>& acc) { return boost::python::incref(PyFloat_FromDouble(acc.get())); }
#else
	static PyObject* convert(const OpenMPAccumulator<Real>& acc) { return ArbitraryReal_to_python<Real>::convert(acc.get()); }
#endif
};
struct custom_OpenMPAccumulator_from_float {
	custom_OpenMPAccumulator_from_float()
	{
		boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<OpenMPAccumulator<Real>>());
	}
#if (YADE_REAL_BIT <= 64)
	static void* convertible(PyObject* obj_ptr) { return PyFloat_Check(obj_ptr) ? obj_ptr : 0; }
#else
	static void*     convertible(PyObject* obj_ptr) { return ArbitraryReal_from_python<Real>::convertible(obj_ptr); }
#endif
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<OpenMPAccumulator<Real>>*)(data))->storage.bytes;
		new (storage) OpenMPAccumulator<Real>;
		((OpenMPAccumulator<Real>*)storage)->set(boost::python::extract<Real>(obj_ptr));
		data->convertible = storage;
	}
};
struct custom_OpenMPAccumulator_to_int {
	static PyObject* convert(const OpenMPAccumulator<int>& acc)
	{
#if PY_MAJOR_VERSION >= 3
		return boost::python::incref(PyLong_FromLong((long)acc.get()));
#else
                return boost::python::incref(PyInt_FromLong((long)acc.get()));
#endif
	}
};

struct custom_OpenMPAccumulator_from_int {
	custom_OpenMPAccumulator_from_int()
	{
		boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<OpenMPAccumulator<int>>());
	}
	static void* convertible(PyObject* obj_ptr)
	{
#if PY_MAJOR_VERSION >= 3
		return PyLong_Check(obj_ptr) ? obj_ptr : 0;
#else
                return PyInt_Check(obj_ptr) ? obj_ptr : 0;
#endif
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<OpenMPAccumulator<int>>*)(data))->storage.bytes;
		new (storage) OpenMPAccumulator<int>;
		((OpenMPAccumulator<int>*)storage)->set(boost::python::extract<int>(obj_ptr));
		data->convertible = storage;
	}
};

template <typename T> struct custom_vvvector_to_list {
	static PyObject* convert(const std::vector<std::vector<std::vector<T>>>& vvv)
	{
		boost::python::list ret;
		for (const auto& vv : vvv) {
			boost::python::list ret2;
			for (const auto& v : vv) {
				boost::python::list ret3;
				for (const auto& e : v)
					ret3.append(e);
				ret2.append(ret3);
			}
			ret.append(ret2);
		}
		return boost::python::incref(ret.ptr());
	}
};

template <typename T> struct custom_vvector_to_list {
	static PyObject* convert(const std::vector<std::vector<T>>& vv)
	{
		boost::python::list ret;
		FOREACH(const std::vector<T>& v, vv)
		{
			boost::python::list ret2;
			FOREACH(const T& e, v) ret2.append(e);
			ret.append(ret2);
		}
		return boost::python::incref(ret.ptr());
	}
};

template <typename containedType> struct custom_list_to_list {
	static PyObject* convert(const std::list<containedType>& v)
	{
		boost::python::list ret;
		FOREACH(const containedType& e, v) ret.append(e);
		return boost::python::incref(ret.ptr());
	}
};
/*** c++-list to python-list */
template <typename containedType> struct custom_vector_to_list {
	static PyObject* convert(const std::vector<containedType>& v)
	{
		boost::python::list ret;
		FOREACH(const containedType& e, v) ret.append(e);
		return boost::python::incref(ret.ptr());
	}
};
template <typename containedType> struct custom_vector_from_seq {
	custom_vector_from_seq()
	{
		boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<std::vector<containedType>>());
	}
	static void* convertible(PyObject* obj_ptr)
	{
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if (!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr, "__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<std::vector<containedType>>*)(data))->storage.bytes;
		new (storage) std::vector<containedType>();
		std::vector<containedType>* v = (std::vector<containedType>*)(storage);
		int                         l = PySequence_Size(obj_ptr);
		if (l < 0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/
		v->reserve(l);
		for (int i = 0; i < l; i++) {
			v->push_back(boost::python::extract<containedType>(PySequence_GetItem(obj_ptr, i)));
		}
		data->convertible = storage;
	}
};

// custom_vvector_from_llist: will build a std::vector< std::vector<someType> > from a list (of list) Python input. Necessary for custom_vvvector_from_lllist below
template <typename someType> struct custom_vvector_from_llist {
	custom_vvector_from_llist()
	{
		boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<std::vector<std::vector<someType>>>());
	}
	static void* convertible(PyObject* obj_ptr)
	{
		if (!PyList_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr, "__len__"))
			return 0; // 1st check, with the 2nd condition kept as in custom_vector_from_seq (for no other reason ?)...
		PyObject* objItem = PyList_GetItem(obj_ptr, 0);
		if (!PyList_Check(objItem)) { // 2nd check
			LOG_ERROR("You did not provide a list of list, returning 0 (None ?)");
			return 0;
		} // now we're sure it is a list of list, at least.
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<std::vector<std::vector<someType>>>*)(data))->storage.bytes;
		new (storage) std::vector<std::vector<someType>>();
		std::vector<std::vector<someType>>* v = (std::vector<std::vector<someType>>*)(storage);
		int                                 l = PySequence_Size(obj_ptr);
		if (l < 0) abort();
		v->reserve(l);
		for (int i = 0; i < l; i++) {
			v->push_back(boost::python::extract<std::vector<someType>>(PySequence_GetItem(obj_ptr, i)));
		}
		data->convertible = storage;
	}
};

// custom_vvvector_from_lllist: will build a std::vector< std::vector< std::vector<someType> > > from a Python list (of list of list) input. E.g. for direct Python definition of LevelSet::distField
template <typename someType> struct custom_vvvector_from_lllist {
	custom_vvvector_from_lllist()
	{
		boost::python::converter::registry::push_back(
		        &convertible, &construct, boost::python::type_id<std::vector<std::vector<std::vector<someType>>>>());
	}
	static void* convertible(PyObject* obj_ptr)
	{
		if (!PyList_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr, "__len__"))
			return 0; // 1st check, with the 2nd condition kept as in custom_vector_from_seq (for no other reason ?)...
		PyObject* objItem = PyList_GetItem(obj_ptr, 0);
		if (!PyList_Check(objItem)) { // 2nd check
			LOG_ERROR("You did not provide a list of list (of list), returning 0 (None ?)");
			return 0;
		} // now we're sure it is a list of list, at least.
		PyObject* objItemItem = PyList_GetItem(objItem, 0);
		if (!PyList_Check(objItemItem)) { // 3rd check
			LOG_ERROR("You did not provide a list of list of list, returning 0 (None ?)");
			return 0;
		} // now we're sure it is a list of list of list
		//we could add a last test to check if someType is a number (using e.g. PyFloat_Check() from Python-C API)
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<std::vector<std::vector<std::vector<someType>>>>*)(data))->storage.bytes;
		new (storage) std::vector<std::vector<std::vector<someType>>>();
		std::vector<std::vector<std::vector<someType>>>* v = (std::vector<std::vector<std::vector<someType>>>*)(storage);
		int                                              l = PySequence_Size(obj_ptr);
		if (l < 0) abort();
		v->reserve(l);
		for (int i = 0; i < l; i++) {
			v->push_back(boost::python::extract<std::vector<std::vector<someType>>>(PySequence_GetItem(obj_ptr, i)));
		}
		data->convertible = storage;
	}
};

struct custom_ptrMatchMaker_from_float {
	custom_ptrMatchMaker_from_float()
	{
		boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<shared_ptr<MatchMaker>>());
	}
	static void* convertible(PyObject* obj_ptr)
	{
		if (!PyNumber_Check(obj_ptr)) {
			cerr << "Not convertible to MatchMaker" << endl;
			return 0;
		}
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<shared_ptr<MatchMaker>>*)(data))->storage.bytes;
		new (storage) shared_ptr<MatchMaker>(new MatchMaker);            // allocate the object at given address
		shared_ptr<MatchMaker>* mm = (shared_ptr<MatchMaker>*)(storage); // convert that address to our type
		(*mm)->algo                = "val";
		(*mm)->val                 = PyFloat_AsDouble(obj_ptr);
		(*mm)->postLoad(**mm);
		data->convertible = storage;
	}
};


#ifdef YADE_MASK_ARBITRARY
struct custom_mask_to_long {
	static PyObject* convert(const mask_t& mask) { return PyLong_FromString(const_cast<char*>(mask.to_string().c_str()), NULL, 2); }
};
struct custom_mask_from_long {
	custom_mask_from_long() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<mask_t>()); }
	static void* convertible(PyObject* obj_ptr)
	{
#if PY_MAJOR_VERSION >= 3
		return PyLong_Check(obj_ptr) ? obj_ptr : 0;
#else
		return (PyLong_Check(obj_ptr) || PyInt_Check(obj_ptr)) ? obj_ptr : 0;
#endif
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
	{
		void* storage = ((boost::python::converter::rvalue_from_python_storage<mask_t>*)(data))->storage.bytes;
		new (storage) mask_t;
		mask_t* mask = (mask_t*)storage;
#if PY_MAJOR_VERSION >= 3
		obj_ptr = _PyLong_Format(obj_ptr, 2);
		std::string s(PyUnicode_AsUTF8(obj_ptr));
#else
		obj_ptr = _PyLong_Format(obj_ptr, 2, 0, 0);
		std::string s(PyString_AsString(obj_ptr));
#endif
		//
		if (s.substr(0, 2).compare("0b") == 0) s = s.substr(2);
		if (s[s.length() - 1] == 'L') s = s.substr(0, s.length() - 1);
		// TODO?
		*mask             = mask_t(s);
		data->convertible = storage;
	}
};
#endif

} // namespace yade

// BOOST_PYTHON_MODULE cannot be inside yade namespace, it has 'extern "C"' keyword, which strips it out of any namespaces.
BOOST_PYTHON_MODULE(_customConverters)
try {
	namespace y = ::yade;

	y::custom_Se3r_from_seq();
	boost::python::to_python_converter<y::Se3r, y::custom_se3_to_tuple>();

	y::custom_OpenMPAccumulator_from_float();
	boost::python::to_python_converter<y::OpenMPAccumulator<y::Real>, y::custom_OpenMPAccumulator_to_float>();
	y::custom_OpenMPAccumulator_from_int();
	boost::python::to_python_converter<y::OpenMPAccumulator<int>, y::custom_OpenMPAccumulator_to_int>();
	// todo: OpenMPAccumulator<int>

	y::custom_ptrMatchMaker_from_float();

	// StrArrayMap (typedef for std::map<std::string,numpy_boost>) → python dictionary
	//custom_StrArrayMap_to_dict();
	// register from-python converter and to-python converter

	// register 1-way C++ -> Python conversion:
	boost::python::to_python_converter<std::vector<std::vector<std::string>>, y::custom_vvector_to_list<std::string>>();
	boost::python::to_python_converter<std::vector<std::vector<y::Matrix3r>>, y::custom_vvector_to_list<y::Matrix3r>>();
	boost::python::to_python_converter<std::vector<std::vector<int>>, y::custom_vvector_to_list<int>>();
	boost::python::to_python_converter<std::vector<std::vector<double>>, y::custom_vvector_to_list<double>>();
	//boost::python::to_python_converter<std::list<shared_ptr<Functor     > >, y::custom_list_to_list<shared_ptr<Functor> > >();
	//boost::python::to_python_converter<std::list<shared_ptr<Functor     > >, y::custom_list_to_list<shared_ptr<Functor> > >();
	boost::python::to_python_converter<std::vector<std::vector<std::vector<y::Real>>>, y::custom_vvvector_to_list<y::Real>>(); // eg for LevelSet.distField
	// register the reciprocal conversion Python -> C++
	y::custom_vvector_from_llist<y::Real>();
	y::custom_vvvector_from_lllist<y::Real>();

#ifdef YADE_MASK_ARBITRARY
	y::custom_mask_from_long();
	boost::python::to_python_converter<y::mask_t, y::custom_mask_to_long>();
#endif

// register 2-way conversion between c++ vector and python homogeneous sequence (list/tuple) of corresponding type
#define VECTOR_SEQ_CONV(Type)                                                                                                                                  \
	y::custom_vector_from_seq<Type>();                                                                                                                     \
	boost::python::to_python_converter<std::vector<Type>, y::custom_vector_to_list<Type>>();
	{
		using namespace yade; // but only inside this { … } block. Keep pollution under control. Make it convenient.
		VECTOR_SEQ_CONV(int);
		VECTOR_SEQ_CONV(bool);
		VECTOR_SEQ_CONV(Real);
		VECTOR_SEQ_CONV(Se3r);
		VECTOR_SEQ_CONV(Vector2r);
		VECTOR_SEQ_CONV(Vector2i);
		VECTOR_SEQ_CONV(Vector3r);
		VECTOR_SEQ_CONV(Vector3i);
		VECTOR_SEQ_CONV(Vector6r);
		VECTOR_SEQ_CONV(Vector6i);
		VECTOR_SEQ_CONV(Matrix3r);
		VECTOR_SEQ_CONV(Matrix6r);
		VECTOR_SEQ_CONV(std::string);
		VECTOR_SEQ_CONV(shared_ptr<Body>);
		VECTOR_SEQ_CONV(shared_ptr<Engine>);
		VECTOR_SEQ_CONV(shared_ptr<Material>);
		VECTOR_SEQ_CONV(shared_ptr<Serializable>);
		VECTOR_SEQ_CONV(shared_ptr<BoundFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IGeomFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IPhysFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<LawFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IntrCallback>);
		VECTOR_SEQ_CONV(shared_ptr<Interaction>);
#ifdef YADE_BODY_CALLBACK
		VECTOR_SEQ_CONV(shared_ptr<BodyCallback>);
#endif
		VECTOR_SEQ_CONV(shared_ptr<SpherePack>);
		VECTOR_SEQ_CONV(shared_ptr<KinematicEngine>);
#ifdef YADE_OPENGL
		VECTOR_SEQ_CONV(shared_ptr<GlBoundFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<GlStateFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<GlShapeFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<GlIGeomFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<GlIPhysFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<GlExtraDrawer>);
#endif
	}
#undef VECTOR_SEQ_CONV

} catch (...) {
	LOG_FATAL("Importing this module caused an exception and this module is in an inconsistent state now.");
	PyErr_Print();
	PyErr_SetString(PyExc_SystemError, __FILE__);
	boost::python::handle_exception();
	throw;
}
