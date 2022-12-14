/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once


#include <lib/multimethods/DynLibDispatcher.hpp>
#include <core/Engine.hpp>
#include <core/Functor.hpp>
#include <core/Omega.hpp>

namespace yade { // Cannot have #include directive inside.

// real base class for all dispatchers (the other one are templates)
class Dispatcher : public Engine {
public:
	// these functions look to be completely unused...?
	virtual string getFunctorType() const { throw; };
	virtual int    getDimension() const { throw; };
	virtual string getBaseClassType(unsigned int) const { throw; };
	//
	virtual ~Dispatcher() {};
	// clang-format off
	YADE_CLASS_BASE_DOC(Dispatcher,Engine,"Engine dispatching control to its associated functors, based on types of argument it receives. This abstract base class provides no functionality in itself.")
	// clang-format on
};
REGISTER_SERIALIZABLE(Dispatcher);

/* Each real dispatcher derives from Dispatcher1D or Dispatcher2D (both templates), which in turn derive from Dispatcher (an Engine) and DynLibDispatcher (the dispatch logic).
Because we need literal functor and class names for registration in python, we provide macro that creates the real dispatcher class with everything needed.
*/

#define _YADE_DISPATCHER1D_FUNCTOR_ADD(FunctorT, f)                                                                                                            \
	virtual void addFunctor(shared_ptr<FunctorT> f) { add1DEntry(f->get1DFunctorType1(), f); }
#define _YADE_DISPATCHER2D_FUNCTOR_ADD(FunctorT, f)                                                                                                            \
	virtual void addFunctor(shared_ptr<FunctorT> f) { add2DEntry(f->get2DFunctorType1(), f->get2DFunctorType2(), f); }

#define _YADE_DIM_DISPATCHER_FUNCTOR_DOC_ATTRS_CTOR_PY(Dim, DispatcherT, FunctorT, doc, attrs, ctor, py)                                                       \
	typedef FunctorT FunctorType;                                                                                                                          \
	void             updateScenePtr()                                                                                                                      \
	{                                                                                                                                                      \
		for (auto& f : functors) {                                                                                                                     \
			f->scene = scene;                                                                                                                      \
		}                                                                                                                                              \
	}                                                                                                                                                      \
	void postLoad(DispatcherT&)                                                                                                                            \
	{                                                                                                                                                      \
		clearMatrix();                                                                                                                                 \
		for (auto& f : functors)                                                                                                                       \
			add(YADE_PTR_CAST<FunctorT>(f));                                                                                                       \
	}                                                                                                                                                      \
	virtual void add(FunctorT* f) { add(shared_ptr<FunctorT>(f)); }                                                                                        \
	virtual void add(shared_ptr<FunctorT> f)                                                                                                               \
	{                                                                                                                                                      \
		bool   dupe = false;                                                                                                                           \
		string fn   = f->getClassName();                                                                                                               \
		for (const auto& fScan : functors) {                                                                                                           \
			if (fn == fScan->getClassName()) dupe = true;                                                                                          \
		}                                                                                                                                              \
		if (!dupe) functors.push_back(f);                                                                                                              \
		addFunctor(f);                                                                                                                                 \
	}                                                                                                                                                      \
	BOOST_PP_CAT(_YADE_DISPATCHER, BOOST_PP_CAT(Dim, D_FUNCTOR_ADD))(FunctorT, f) boost::python::list functors_get(void) const                             \
	{                                                                                                                                                      \
		boost::python::list ret;                                                                                                                       \
		for (const auto& f : functors) {                                                                                                               \
			ret.append(f);                                                                                                                         \
		}                                                                                                                                              \
		return ret;                                                                                                                                    \
	}                                                                                                                                                      \
	void functors_set(const vector<shared_ptr<FunctorT>>& ff)                                                                                              \
	{                                                                                                                                                      \
		functors.clear();                                                                                                                              \
		for (const auto& f : ff)                                                                                                                       \
			add(f);                                                                                                                                \
		postLoad(*this);                                                                                                                               \
	}                                                                                                                                                      \
	void pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict& /*d*/) override                                                              \
	{                                                                                                                                                      \
		if (boost::python::len(t) == 0) return;                                                                                                        \
		if (boost::python::len(t) != 1) throw invalid_argument("Exactly one list of " BOOST_PP_STRINGIZE(FunctorT) " must be given.");                 \
		typedef std::vector<shared_ptr<FunctorT>> vecF;                                                                                                \
		vecF                                      vf = boost::python::extract<vecF>(t[0])();                                                           \
		functors_set(vf);                                                                                                                              \
		t = boost::python::tuple();                                                                                                                    \
	}                                                                                                                                                      \
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(                                                                                                                     \
	        DispatcherT,                                                                                                                                   \
	        Dispatcher,                                                                                                                                    \
	        "Dispatcher calling :yref:`functors<" BOOST_PP_STRINGIZE(FunctorT) ">` based on received argument type(s).\n\n" doc,                                                                            \
	                ((vector<shared_ptr<FunctorT>>,                                                                                                        \
	                  functors,                                                                                                                            \
	                  ,                                                                                                                                    \
	                  ,                                                                                                                                    \
	                  "Functors active in the dispatch mechanism [overridden below].")) /*additional attrs*/ attrs,                                        \
	                /*ctor*/ ctor, /*py*/                                                                                                                  \
	                py.add_property("functors",                                                                                                            \
	                                &DispatcherT::functors_get,                                                                                            \
	                                &DispatcherT::functors_set,                                                                                            \
	                                "Functors associated with this dispatcher."                                                                            \
	                                " :yattrtype:`vector<shared_ptr<" BOOST_PP_STRINGIZE(FunctorT) "> >` ")                                                \
	                                        .def("dispMatrix",                                                                                             \
	                                             &DispatcherT::dump,                                                                                       \
	                                             boost::python::arg("names") = true,                                                                       \
	                                             "Return dictionary with contents of the dispatch matrix.")                                                \
	                                        .def("dispFunctor",                                                                                            \
	                                             &DispatcherT::getFunctor,                                                                                 \
	                                             "Return functor that would be dispatched for given argument(s); None if no dispatch; ambiguous dispatch " \
	                                             "throws.");)

#define YADE_DISPATCHER1D_FUNCTOR_DOC_ATTRS_CTOR_PY(DispatcherT, FunctorT, doc, attrs, ctor, py)                                                               \
	_YADE_DIM_DISPATCHER_FUNCTOR_DOC_ATTRS_CTOR_PY(1, DispatcherT, FunctorT, doc, attrs, ctor, py)
#define YADE_DISPATCHER2D_FUNCTOR_DOC_ATTRS_CTOR_PY(DispatcherT, FunctorT, doc, attrs, ctor, py)                                                               \
	_YADE_DIM_DISPATCHER_FUNCTOR_DOC_ATTRS_CTOR_PY(2, DispatcherT, FunctorT, doc, attrs, ctor, py)

// HELPER FUNCTIONS

/*! Function returning class name (as string) for given index and topIndexable (top-level indexable, such as Shape, Material and so on)
This function exists solely for debugging, is quite slow: it has to traverse all classes and ask for inheritance information.
It should be used primarily to convert indices to names in Dispatcher::dictDispatchMatrix?D; since it relies on Omega for RTTI,
this code could not be in Dispatcher itself.
s*/
template <class topIndexable> std::string Dispatcher_indexToClassName(int idx)
{
	boost::scoped_ptr<topIndexable> top(new topIndexable);
	std::string                     topName = top->getClassName();
	for (const auto& clss : Omega::instance().getDynlibsDescriptor()) {
		if (Omega::instance().isInheritingFrom_recursive(clss.first, topName) || clss.first == topName) {
			// create instance, to ask for index
			shared_ptr<topIndexable> inst = YADE_PTR_DYN_CAST<topIndexable>(ClassFactory::instance().createShared(clss.first));
			assert(inst);
			if (inst->getClassIndex() < 0 && inst->getClassName() != top->getClassName()) {
				throw logic_error(
				        "Class " + inst->getClassName() + " didn't use REGISTER_CLASS_INDEX(" + inst->getClassName() + "," + top->getClassName()
				        + ") and/or forgot to call createIndex() in the ctor. [[ Please fix that! ]]");
			}
			if (inst->getClassIndex() == idx) return clss.first;
		}
	}
	throw runtime_error("No class with index " + boost::lexical_cast<string>(idx) + " found (top-level indexable is " + topName + ")");
}

//! Return class index of given indexable
template <typename TopIndexable> int Indexable_getClassIndex(const shared_ptr<TopIndexable> i) { return i->getClassIndex(); }

//! Return sequence (hierarchy) of class indices of given indexable; optionally convert to names
template <typename TopIndexable> boost::python::list Indexable_getClassIndices(const shared_ptr<TopIndexable> i, bool convertToNames)
{
	int                 depth = 1;
	boost::python::list ret;
	int                 idx0 = i->getClassIndex();
	if (convertToNames) ret.append(Dispatcher_indexToClassName<TopIndexable>(idx0));
	else
		ret.append(idx0);

	if (idx0 < 0) return ret; // don't continue and call getBaseClassIndex(), since we are at the top already

	while (true) {
		int idx = i->getBaseClassIndex(depth++);
		if (convertToNames) ret.append(Dispatcher_indexToClassName<TopIndexable>(idx));
		else
			ret.append(idx);
		if (idx < 0) return ret;
	}
}

template <class FunctorType, bool autoSymmetry = true>
class Dispatcher1D : public Dispatcher,
                     public DynLibDispatcher<
                             TYPELIST_1(typename FunctorType::DispatchType1) // base classes for dispatch
                             ,
                             FunctorType // class that provides multivirtual call
                             ,
                             typename FunctorType::ReturnType // return type
                             ,
                             typename FunctorType::ArgumentTypes,
                             autoSymmetry> {
public:
	typedef typename FunctorType::DispatchType1 baseClass;
	typedef baseClass                           argType1;
	typedef FunctorType                         functorType;
	typedef DynLibDispatcher<TYPELIST_1(baseClass), FunctorType, typename FunctorType::ReturnType, typename FunctorType::ArgumentTypes, autoSymmetry>
	        dispatcherBase;

	shared_ptr<FunctorType> getFunctor(shared_ptr<baseClass> arg) { return dispatcherBase::getExecutor(arg); }
	boost::python::dict     dump(bool convertIndicesToNames)
	{
		boost::python::dict ret;
		for (const auto& item : dispatcherBase::dataDispatchMatrix1D()) {
			if (convertIndicesToNames) {
				string arg1                          = Dispatcher_indexToClassName<argType1>(item.ix1);
				ret[boost::python::make_tuple(arg1)] = item.functorName;
			} else
				ret[boost::python::make_tuple(item.ix1)] = item.functorName;
		}
		return ret;
	}

	int getDimension() const override { return 1; }

	string getFunctorType() const override
	{
		shared_ptr<FunctorType> eu(new FunctorType);
		return eu->getClassName();
	}

	string getBaseClassType(unsigned int i) const override
	{
		if (i == 0) {
			shared_ptr<baseClass> bc(new baseClass);
			return bc->getClassName();
		} else
			return "";
	}

public:
	REGISTER_ATTRIBUTES(Dispatcher, );
	REGISTER_CLASS_AND_BASE(Dispatcher1D, Dispatcher DynLibDispatcher);
};


template <class FunctorType, bool autoSymmetry = true>
class Dispatcher2D : public Dispatcher,
                     public DynLibDispatcher<
                             TYPELIST_2(typename FunctorType::DispatchType1, typename FunctorType::DispatchType2) // base classes for dispatch
                             ,
                             FunctorType // class that provides multivirtual call
                             ,
                             typename FunctorType::ReturnType // return type
                             ,
                             typename FunctorType::ArgumentTypes // argument of engine unit
                             ,
                             autoSymmetry> {
public:
	typedef typename FunctorType::DispatchType1 baseClass1;
	typedef typename FunctorType::DispatchType2 baseClass2;
	typedef baseClass1                          argType1;
	typedef baseClass2                          argType2;
	typedef FunctorType                         functorType;
	typedef DynLibDispatcher<
	        TYPELIST_2(baseClass1, baseClass2),
	        FunctorType,
	        typename FunctorType::ReturnType,
	        typename FunctorType::ArgumentTypes,
	        autoSymmetry>
	        dispatcherBase;

	shared_ptr<FunctorType> getFunctor(shared_ptr<baseClass1> arg1, shared_ptr<baseClass2> arg2) { return dispatcherBase::getExecutor(arg1, arg2); }

	boost::python::dict dump(bool convertIndicesToNames)
	{
		boost::python::dict ret;
		for (const auto& item : dispatcherBase::dataDispatchMatrix2D()) {
			if (convertIndicesToNames) {
				string arg1 = Dispatcher_indexToClassName<argType1>(item.ix1), arg2 = Dispatcher_indexToClassName<argType2>(item.ix2);
				ret[boost::python::make_tuple(arg1, arg2)] = item.functorName;
			} else
				ret[boost::python::make_tuple(item.ix1, item.ix2)] = item.functorName;
		}
		return ret;
	}

	int getDimension() const override { return 2; }

	string getFunctorType() const override
	{
		shared_ptr<FunctorType> eu(new FunctorType);
		return eu->getClassName();
	}
	string getBaseClassType(unsigned int i) const override
	{
		if (i == 0) {
			shared_ptr<baseClass1> bc(new baseClass1);
			return bc->getClassName();
		} else if (i == 1) {
			shared_ptr<baseClass2> bc(new baseClass2);
			return bc->getClassName();
		} else
			return "";
	}

public:
	REGISTER_ATTRIBUTES(Dispatcher, );
	REGISTER_CLASS_AND_BASE(Dispatcher2D, Dispatcher DynLibDispatcher);
};

} // namespace yade
