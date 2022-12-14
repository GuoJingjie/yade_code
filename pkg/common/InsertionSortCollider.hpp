// 2009 © Václav Šmilauer <eudoxos@arcig.cz>
// 2013 © Bruno Chareyre <bruno.chareyre@grenoble-inp.fr>
#pragma once
#include <core/Scene.hpp>
#include <pkg/common/Collider.hpp>
#include <pkg/dem/NewtonIntegrator.hpp>

namespace yade { // Cannot have #include directive inside.

class InteractionContainer;

/*! Periodic collider notes.

Architecture
============
Values from bounding boxes are added information about period in which they are.
Their container (VecBounds) holds position of where the space wraps.
The sorting algorithm is changed in such way that periods are changed when body crosses cell boundary.

Interaction::cellDist holds information about relative cell coordinates of the 2nd body
relative to the 1st one. Dispatchers (IGeomDispatcher and InteractionLoop)
use this information to pass modified position of the 2nd body to IGeomFunctors.
Since properly behaving IGeomFunctor's and LawFunctor's do not take positions
directly from Body::state, the interaction is computed with the periodic positions.

Positions of bodies (in the sense of Body::state) and their natural bboxes are not wrapped
to the periodic cell, they can be anywhere (but not "too far" in the sense of int overflow).

Since Interaction::cellDists holds cell coordinates, it is possible to change the cell boundaries
at runtime. This should make it easy to implement some stress control on the top.

Clumps do not interfere with periodicity in any way.

Rendering
---------
OpenGLRenderer renders Shape at all periodic positions that touch the
periodic cell (i.e. Bounds crosses its boundary).

It seems to affect body selection somehow, but that is perhaps not related at all.

Periodicity control
===================
c++:
	Scene::isPeriodic, Scene::cellSize
python:
	O.periodicCell=((0,0,0),(10,10,10)  # activates periodic boundary
	O.periodicCell=() # deactivates periodic boundary

Requirements
============
* By default, no body can have Aabb larger than about .499*cellSize. Exception is thrown if that is false.
	Large bodies are accepted if allowBiggerThanPeriod (experimental)
* Constitutive law must not get body positions from Body::state directly.
	If it does, it uses Interaction::cellDist to compute periodic position.
* No body can get further away than MAXINT periods. It will do horrible things if there is overflow. Not checked at the moment.
* Body cannot move over distance > cellSize in one step. Since body size is limited as well, that would mean simulation explosion.
	Exception is thrown if the movement is upwards. If downwards, it is not handled at all.

Possible performance improvements & bugs
========================================

* PeriodicInsertionSortCollider::{cellWrap,cellWrapRel} OpenGLRenderer::{wrapCell,wrapCellPt} Shop::PeriodicWrap
	are all very similar functions. They should be put into separate header and included from all those places.

*/


// #define this macro to enable timing within this engine
// #define ISC_TIMING

#ifdef ISC_TIMING
#define ISC_CHECKPOINT(cpt) timingDeltas->checkpoint(cpt)
#else
#define ISC_CHECKPOINT(cpt)
#endif

class GeneralIntegratorInsertionSortCollider; // Forward decleration of child to decleare it as friend

class InsertionSortCollider : public Collider {
	friend class GeneralIntegratorInsertionSortCollider;
	//! struct for storing bounds of bodies
	struct Bounds {
		//! coordinate along the given sortAxis
		Real coord;
		//! id of the body this bound belongs to
		Body::id_t id;
		//! periodic cell coordinate
		int period;
		//! is it the minimum (true) or maximum (false) bound?
		struct {
			bool hasBB : true, isMin : true;
		} flags;
		Bounds(Real coord_, Body::id_t id_, bool isMin)
		        : coord(coord_)
		        , id(id_)
		        , period(0)
		{
			flags.isMin = isMin;
		}
		bool operator<(const Bounds& b) const
		{
			/* handle special case of zero-width bodies, which could otherwise get min/max swapped in the unstable std::sort */
			if (id == b.id && coord == b.coord) return flags.isMin;
			return coord < b.coord;
		}
		bool operator>(const Bounds& b) const
		{
			if (id == b.id && coord == b.coord) return !flags.isMin;
			return coord > b.coord;
		}
	};

	// if False, no type of striding is used
	// if True, then either verletDist XOR nBins is set
	bool strideActive;
	struct VecBounds {
		// axis set in the ctor
		int  axis;
		Real cellDim;
		// index of the lowest coordinate element, before which the container wraps
		size_t loIdx;

		Bounds& operator[](long idx)
		{
			assert(idx < long(size()) && idx >= 0);
			return vec[idx];
		}
		const Bounds& operator[](long idx) const
		{
			assert(idx < long(size()) && idx >= 0);
			return vec[idx];
		}

		// update number of bodies, periodic properties and size from Scene
		void updatePeriodicity(Scene*);
		// normalize given index to the right range (wraps around)
		size_t norm(long i) const
		{
			if (i < 0) i += size();
			assert(i >= 0);
			size_t ret = i % size();
			assert(ret < size());
			return ret;
		}
		VecBounds()
		        : axis(-1)
		        , loIdx(0)
		{
		}
		void dump(std::ostream& os)
		{
			string ret;
			for (size_t i = 0; i < vec.size(); i++)
				os << (i == loIdx ? "@@ " : "") << vec[i].coord << "(id=" << vec[i].id << "," << (vec[i].flags.isMin ? "min" : "max") << ",p"
				   << vec[i].period << ") ";
			os << endl;
		}

		size_t size() const { return vec.size(); };

		void clear() { vec.clear(); }
		void reserve(size_t n) { vec.reserve(n); }
		void resize(size_t n)
		{
			if (n > vec.size()) LOG_ERROR("not supposed to increase size - shrink only");
			vec.resize(n, Bounds(0, 0, true));
		}
		void push_back(const Bounds& bb) { vec.push_back(bb); }
		// if the line below does not compile on older ubuntu 14.04, then I should add #ifdef guards to check compiler version. This line will make push_back faster when a newer compiler supports it.
		void                                push_back(Bounds&& bb) { vec.push_back(bb); }
		void                                sort() { std::sort(vec.begin(), vec.end()); }
		std::vector<Bounds>::const_iterator cbegin() const { return vec.cbegin(); }
		std::vector<Bounds>::const_iterator cend() const { return vec.cend(); }
		std::vector<Bounds>::iterator       begin() { return vec.begin(); }
		std::vector<Bounds>::iterator       end() { return vec.end(); }
		std::vector<Bounds>::iterator       insert(std::vector<Bounds>::iterator pos, Bounds& bb) { return vec.insert(pos, bb); }

	private:
		std::vector<Bounds> vec;
	};

private:
	//! storage for bounds
	VecBounds BB[3];
	//! storage for bb maxima and minima
	std::vector<Real> maxima, minima;
	//! Whether the Scene was periodic (to detect the change, which shouldn't happen, but shouldn't crash us either)
	bool periodic;
	//! Store inverse sizes to avoid repeated divisions within loops
	Vector3r invSizes;
	// return python representation of the BB struct, as ([...],[...],[...]).
	boost::python::tuple dumpBounds();

	/*! sorting routine; insertion sort is very fast for strongly pre-sorted lists, which is our case
  	    http://en.wikipedia.org/wiki/Insertion_sort has the algorithm and other details
	*/
	void insertionSort(VecBounds& v, InteractionContainer*, Scene*, bool doCollide = true);
	void findOverlaps(VecBounds& V);
#ifdef YADE_OPENMP
	void insertionSortParallel(VecBounds& v, InteractionContainer*, Scene*, bool doCollide = true);
#endif
	void handleBoundInversion(Body::id_t, Body::id_t, InteractionContainer*, Scene*);

	// periodic variants
	void insertionSortPeri(VecBounds& v, InteractionContainer*, Scene*, bool doCollide = true);
	void handleBoundInversionPeri(Body::id_t, Body::id_t, InteractionContainer*, Scene*);
	void handleBoundSplit(Body::id_t, Body::id_t, InteractionContainer*, Scene*);

	bool        spatialOverlapPeri(Body::id_t, Body::id_t, Scene*, Vector3i&) const;
	inline bool spatialOverlap(const Body::id_t& id1, const Body::id_t& id2) const
	{
		assert(!periodic);
		return (minima[3 * id1 + 0] <= maxima[3 * id2 + 0]) && (maxima[3 * id1 + 0] >= minima[3 * id2 + 0])
		        && (minima[3 * id1 + 1] <= maxima[3 * id2 + 1]) && (maxima[3 * id1 + 1] >= minima[3 * id2 + 1])
		        && (minima[3 * id1 + 2] <= maxima[3 * id2 + 2]) && (maxima[3 * id1 + 2] >= minima[3 * id2 + 2]);
	}

	static Real cellWrap(const Real, const Real, const Real, int&);
	static Real cellWrapRel(const Real, const Real, const Real);


public:
	//! Predicate called from loop within InteractionContainer::erasePending
	bool shouldBeErased(Body::id_t id1, Body::id_t id2, Scene* rb) const
	{
		if (!periodic) return !spatialOverlap(id1, id2);
		else {
			Vector3i periods;
			return !spatialOverlapPeri(id1, id2, rb, periods);
		}
	}
	bool isActivated() override;

	// force reinitialization at next run
	void invalidatePersistentData() override
	{
		for (int i = 0; i < 3; i++) {
			BB[i].clear();
		}
	}

	vector<Body::id_t> probeBoundingVolume(const Bound&) override;

	void action() override;
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(InsertionSortCollider,Collider,"\
		Collider with O(n log(n)) complexity, using :yref:`Aabb` for bounds.\
		\n\n\
		At the initial step, Bodies' bounds (along :yref:`sortAxis<InsertionSortCollider.sortAxis>`) are first std::sort'ed along this (sortAxis) axis, then collided. The initial sort has :math:`O(n^2)` complexity, see `Colliders' performance <https://yade-dem.org/wiki/Colliders_performace>`_ for some information (There are scripts in examples/collider-perf for measurements). \
		\n\n \
		Insertion sort is used for sorting the bound list that is already pre-sorted from last iteration, where each inversion	calls checkOverlap which then handles either overlap (by creating interaction if necessary) or its absence (by deleting interaction if it is only potential).	\
		\n\n \
		Bodies without bounding volume (such as clumps) are handled gracefully and never collide. Deleted bodies are handled gracefully as well.\
		\n\n \
		This collider handles periodic boundary conditions. There are some limitations, notably:\
		\n\n \
			#. No body can have Aabb larger than cell's half size in that respective dimension. You get exception if it does and gets in interaction. One way to explicitly by-pass this restriction is offered by ``allowBiggerThanPeriod``, which can be turned on to insert a floor in the form of a very large box for instance (see examples/periodicSandPile.py). \
			\n\n \
			#. No body can travel more than cell's distance in one step; this would mean that the simulation is numerically exploding, and it is only detected in some cases.\
		\n\n \
		**Stride** can be used to avoid running collider at every step by enlarging the particle's bounds, tracking their displacements and only re-run if they might have gone out of that bounds (see `Verlet list <http://en.wikipedia.org/wiki/Verlet_list>`_ for brief description and background) . This requires cooperation from :yref:`NewtonIntegrator` as well as :yref:`BoundDispatcher`, which will be found among engines automatically (exception is thrown if they are not found).\
		\n\n \
		If you wish to use strides, set ``verletDist`` (length by which bounds will be enlarged in all directions) to some value, e.g. 0.05 × typical particle radius. This parameter expresses the tradeoff between many potential interactions (running collider rarely, but with longer exact interaction resolution phase) and few potential interactions (running collider more frequently, but with less exact resolutions of interactions); it depends mainly on packing density and particle radius distribution.\
		\n\n \
		If ``targetInterv`` is >1, not all particles will have their bound enlarged by ``verletDist``; instead, they will have bounds increased by a length in order to trigger a new colliding after ``targetInterv`` iteration, assuming they move at almost constant velocity. Ideally in this method, all particles would reach their bounds at the sime iteration. This is of course not the case as soon as velocities fluctuate in time. :yref:`Bound::sweepLength` is tuned on the basis of the displacement recorded between the last two runs of the collider. In this situation, ``verletDist`` defines the maximum sweep length.\
		",
		((int,sortAxis,0,,"Axis for the initial contact detection."))
		((bool,allowBiggerThanPeriod,false,,"If true, tests on bodies sizes will be disabled, and the simulation will run normaly even if bodies larger than period are found. It can be useful when the periodic problem include e.g. a floor modelized with wall/box/facet.\nBe sure you know what you are doing if you touch this flag. The result is undefined if one large body moves out of the (0,0,0) period."))
		((bool,sortThenCollide,false,,"Separate sorting and colliding phase; it is MUCH slower, but all interactions are processed at every step; this effectively makes the collider non-persistent, not remembering last state. (The default behavior relies on the fact that inversions during insertion sort are overlaps of bounding boxes that just started/ceased to exist, and only processes those; this makes the collider much more efficient.)"))
		((int,targetInterv,100,,"(experimental) Target number of iterations between bound update, used to define a smaller sweep distance for slower grains if >0, else always use 1*verletDist. Useful in simulations with strong velocity contrasts between slow bodies and fast bodies."))
		((Real,overlapTolerance,1e-7,,"Tolerance on determining overlap. In rare cases different parts of the code can inconsistently lead to different results in terms of overlap, with false negative by spatialOverlapPeri possibly leading to nasty bugs in contact detection (false positive are harmless). This tolerance is to avoid false negative, the value can be understood as relative to 1 (i.e. independent of particle size or any other reference length). The default should be ok."))
		((Real,updatingDispFactor,-1,,"(experimental) Displacement factor used to trigger bound update: the bound is updated only if updatingDispFactor*disp>sweepDist when >0, else all bounds are updated."))
		((Real,verletDist,((void)"Automatically initialized",-.5),,"Length by which to enlarge particle bounds, to avoid running collider at every step. Stride disabled if zero. Negative value will trigger automatic computation, so that the real value will be *verletDist* × minimum spherical particle radius; if there are no spherical particles, it will be disabled. The actual length added to one bound can be only a fraction of verletDist when :yref:`InsertionSortCollider::targetInterv` is > 0."))
		((Real,minSweepDistFactor,0.1,,"Minimal distance by which enlarge all bounding boxes; superseeds computed value of verletDist when lower that (minSweepDistFactor x verletDist)."))
		((Real,fastestBodyMaxDist,0,,"Normalized maximum displacement of the fastest body since last run; if >= 1, we could get out of bboxes and will trigger full run. |yupdate|"))
		((int,numReinit,0,Attr::readonly,"Cummulative number of bound array re-initialization."))
		((int,numAction,0,,"Cummulative number of collision detection."))
		((bool,doSort,false,,"Do forced resorting of interactions."))
		((bool,keepListsShort,false,,"if true remove bounds of non-existent or unbounded bodies from the lists |yupdate|; turned true automatically in MPI mode and if bodies are erased with :yref:`BodyContainer.enableRedirection`=True."))
		((bool,smartInsertErase,false,,"Use an algorithm optimized for heavy insert/delete (avoid initSort) - experimental."))
		((shared_ptr<NewtonIntegrator>, newton,shared_ptr<NewtonIntegrator>(),,"reference to active :yref:`Newton integrator<NewtonIntegrator>`. |yupdate|"))
		, /* ctor */
			#ifdef ISC_TIMING
				timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
			#endif 
			for(int i=0; i<3; i++) BB[i].axis=i;
			periodic=false;
			strideActive=false;
			,
		/* py */
		.def_readonly("strideActive",&InsertionSortCollider::strideActive,"Whether striding is active (read-only; for debugging). |yupdate|")
		.def_readonly("periodic",&InsertionSortCollider::periodic,"Whether the collider is in periodic mode (read-only; for debugging) |yupdate|")
		.def("dumpBounds",&InsertionSortCollider::dumpBounds,"Return representation of the internal sort data. The format is ``([...],[...],[...])`` for 3 axes, where each ``...`` is a list of entries (bounds). The entry is a tuple with the fllowing items:\n\n* coordinate (float)\n* body id (int), but negated for negative bounds\n* period numer (int), if the collider is in the periodic regime.")
		.def("isActivated",&InsertionSortCollider::isActivated,"Return true if collider needs execution at next iteration.")
	);
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(InsertionSortCollider);

} // namespace yade
