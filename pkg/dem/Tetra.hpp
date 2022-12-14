// © 2007 Václav Šmilauer <eudoxos@arcig.cz>
// © 2013 Jan Stránský <jan.stransky@fsv.cvut.cz>

#pragma once

#include <core/GlobalEngine.hpp>
#include <core/IGeom.hpp>
#include <core/Shape.hpp>

#include <core/Aabb.hpp>
#include <core/Dispatching.hpp>
#include <pkg/common/NormShearPhys.hpp>

#ifdef YADE_CGAL
#include <CGAL/Cartesian.h>
//#include <CGAL/intersections.h>
#endif

#ifdef YADE_OPENGL
#include <pkg/common/GLDrawFunctors.hpp>
#endif

namespace yade { // Cannot have #include directive inside.

/* Our mold of tetrahedron: just 4 vertices.
 *
 * Self-contained. */
class Tetra : public Shape {
public:
	Tetra(Vector3r v0, Vector3r v1, Vector3r v2, Vector3r v3)
	{
		createIndex();
		v.resize(4);
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		v[3] = v3;
	}
	virtual ~Tetra();
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Tetra,Shape,"Tetrahedron geometry.",
		((std::vector<Vector3r>,v,std::vector<Vector3r>(4),,"Tetrahedron vertices (in local coordinate system).")),
		/*ctor*/createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(Tetra, Shape);
};
REGISTER_SERIALIZABLE(Tetra);


/*! Collision configuration for Tetra and something.
 * This is expressed as penetration volume properties: centroid, volume, orientation of principal axes, inertia.
 *
 * Self-contained. */
class TTetraGeom : public IGeom {
public:
	virtual ~TTetraGeom();
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(TTetraGeom,IGeom,"Geometry of interaction between 2 :yref:`tetrahedra<Tetra>`, including volumetric characteristics",
		((Real,penetrationVolume,NaN,,"Volume of overlap [m³]"))
		((Real,equivalentCrossSection,NaN,,"Cross-section of the overlap (perpendicular to the axis of least inertia"))
		((Real,maxPenetrationDepthA,NaN,,"??"))
		((Real,maxPenetrationDepthB,NaN,,"??"))
		((Real,equivalentPenetrationDepth,NaN,,"??"))
		((Vector3r,contactPoint,,,"Contact point (global coords)"))
		((Vector3r,normal,,,"Normal of the interaction, directed in the sense of least inertia of the overlap volume")),
		createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(TTetraGeom, IGeom);
};
REGISTER_SERIALIZABLE(TTetraGeom);


class TTetraSimpleGeom : public IGeom {
public:
	virtual ~TTetraSimpleGeom();
	// clang-format off
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(TTetraSimpleGeom,IGeom,"EXPERIMENTAL. Geometry of interaction between 2 :yref:`tetrahedra<Tetra>`",
		((Real,penetrationVolume,NaN,,"Volume of overlap [m³]"))
		((Vector3r,contactPoint,,,"Contact point (global coords)"))
		((Vector3r,normal,,,"Normal of the interaction TODO"))
		((int,flag,0,,"TODO")),
		createIndex();
	);
	// clang-format on
	REGISTER_CLASS_INDEX(TTetraSimpleGeom, IGeom);
};
REGISTER_SERIALIZABLE(TTetraSimpleGeom);


/*! Creates Aabb from Tetra. 
 *
 * Self-contained. */
class Bo1_Tetra_Aabb : public BoundFunctor {
public:
	void go(const shared_ptr<Shape>& ig, shared_ptr<Bound>& bv, const Se3r& se3, const Body*) override;
	FUNCTOR1D(Tetra);
	// clang-format off
	YADE_CLASS_BASE_DOC(Bo1_Tetra_Aabb,BoundFunctor,"Create/update :yref:`Aabb` of a :yref:`Tetra`");
	// clang-format on
};
REGISTER_SERIALIZABLE(Bo1_Tetra_Aabb);

#ifdef YADE_OPENGL
/*! Draw Tetra using OpenGL */
class Gl1_Tetra : public GlShapeFunctor {
public:
	void go(const shared_ptr<Shape>&, const shared_ptr<State>&, bool, const GLViewInfo&) override;
	// clang-format off
		YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_Tetra,GlShapeFunctor,"Renders :yref:`Tetra` object",
			((bool,wire,true,,"TODO"))
		);
	// clang-format on
	RENDERS(Tetra);
};
REGISTER_SERIALIZABLE(Gl1_Tetra);
#endif

/*! Calculate physical response based on penetration configuration given by TTetraGeom. */

class TetraVolumetricLaw : public GlobalEngine {
public:
	void action() override;
	DECLARE_LOGGER;
	// clang-format off
	YADE_CLASS_BASE_DOC(TetraVolumetricLaw,GlobalEngine,"Calculate physical response of 2 :yref:`tetrahedra<Tetra>` in interaction, based on penetration configuration given by :yref:`TTetraGeom`.");
	// clang-format on
};
REGISTER_SERIALIZABLE(TetraVolumetricLaw);


/*! @fixme implement Tetra2BoxBang by representing box as 6 tetrahedra. */

/*! Create TTetraGeom (collision geometry) from colliding Tetra's. */
class Ig2_Tetra_Tetra_TTetraGeom : public IGeomFunctor {
public:
	virtual bool
	             go(const shared_ptr<Shape>&       cm1,
	                const shared_ptr<Shape>&       cm2,
	                const State&                   state1,
	                const State&                   state2,
	                const Vector3r&                shift2,
	                const bool&                    force,
	                const shared_ptr<Interaction>& c) override;
	virtual bool goReverse(
	        const shared_ptr<Shape>& /*cm1*/,
	        const shared_ptr<Shape>& /*cm2*/,
	        const State& /*state1*/,
	        const State& /*state2*/,
	        const Vector3r& /*shift2*/,
	        const bool& /*force*/,
	        const shared_ptr<Interaction>& /*c*/) override
	{
		throw std::logic_error("Ig2_Tetra_Tetra_TTetraGeom::goReverse called, but the functor is symmetric.");
	}
	FUNCTOR2D(Tetra, Tetra);
	DEFINE_FUNCTOR_ORDER_2D(Tetra, Tetra);
	// clang-format off
		YADE_CLASS_BASE_DOC(Ig2_Tetra_Tetra_TTetraGeom,IGeomFunctor,"Create/update geometry of collision between 2 :yref:`tetrahedra<Tetra>` (:yref:`TTetraGeom` instance)");
	// clang-format on
	DECLARE_LOGGER;

private:
	std::list<Tetra> Tetra2TetraIntersection(const Tetra& A, const Tetra& B);
	std::list<Tetra> TetraClipByPlane(const Tetra& T, const Vector3r& P, const Vector3r& n);
	//! Intersection of line given by points A, B and plane given by P and its normal.
	Vector3r PtPtPlaneIntr(const Vector3r& A, const Vector3r& B, const Vector3r& P, const Vector3r& normal)
	{
		const Real t = (P - A).dot(normal) / (B - A).dot(normal);
		return A + t * (B - A);
	}
};

REGISTER_SERIALIZABLE(Ig2_Tetra_Tetra_TTetraGeom);


#ifdef YADE_CGAL
class Ig2_Tetra_Tetra_TTetraSimpleGeom : public IGeomFunctor {
protected:
	typedef CGAL::Cartesian<Real>  K;
	typedef K::Point_3             Point;
	typedef CGAL::Tetrahedron_3<K> CGALTetra;
	typedef CGAL::Segment_3<K>     Segment;
	typedef CGAL::Triangle_3<K>    Triangle;
	typedef CGAL::Vector_3<K>      Vector_3;
	bool                           checkVertexToTriangleCase(
	                                  const Triangle tA[4], const Point pB[4], const Segment sB[6], Vector3r& normal, Vector3r& contactPoint, Real& penetrationVolume);
	bool checkEdgeToEdgeCase(
	        const Segment  sA[6],
	        const Segment  sB[6],
	        const Triangle tA[4],
	        const Triangle tB[4],
	        Vector3r&      normal,
	        Vector3r&      contactPoint,
	        Real&          penetrationVolume);
	bool checkEdgeToTriangleCase1(
	        const Triangle tA[4], const Segment sB[6], const Point pB[4], Vector3r& normal, Vector3r& contactPoint, Real& penetrationVolume);
	bool checkEdgeToTriangleCase2(
	        const Triangle tA[4],
	        const Triangle tB[4],
	        const Segment  sA[6],
	        const Segment  sB[6],
	        Vector3r&      normal,
	        Vector3r&      contactPoint,
	        Real&          penetrationVolume);
	bool checkVertexToEdgeCase(
	        const Point    pA[4],
	        const Segment  sA[6],
	        const Triangle tA[4],
	        const Segment  sB[6],
	        const Triangle tB[4],
	        Vector3r&      normal,
	        Vector3r&      contactPoint,
	        Real&          penetrationVolume);
	bool checkVertexToVertexCase(
	        const Point    pA[4],
	        const Point    pB[4],
	        const Segment  sA[6],
	        const Triangle tA[4],
	        const Triangle tB[4],
	        Vector3r&      normal,
	        Vector3r&      contactPoint,
	        Real&          penetrationVolume);
	static const int psMap[4][3];
	static const int ptMap[4][3];
	static const int stMap[6][2];
	static const int tsMap[4][3];
	static const int ppsMap[4][4];
	static const int sstMap[6][6];

public:
	virtual bool
	go(const shared_ptr<Shape>&       cm1,
	   const shared_ptr<Shape>&       cm2,
	   const State&                   state1,
	   const State&                   state2,
	   const Vector3r&                shift2,
	   const bool&                    force,
	   const shared_ptr<Interaction>& c) override;
	virtual bool
	goReverse(const shared_ptr<Shape>&, const shared_ptr<Shape>&, const State&, const State&, const Vector3r&, const bool&, const shared_ptr<Interaction>&)
	        override
	{
		throw std::logic_error("Ig2_Tetra_Tetra_TTetraSimpleGeom::goReverse called, but the functor is symmetric.");
	}
	FUNCTOR2D(Tetra, Tetra);
	DEFINE_FUNCTOR_ORDER_2D(Tetra, Tetra);
	// clang-format off
	YADE_CLASS_BASE_DOC(Ig2_Tetra_Tetra_TTetraSimpleGeom,IGeomFunctor,"EXPERIMANTAL. Create/update geometry of collision between 2 :yref:`tetrahedra<Tetra>` (:yref:`TTetraSimpleGeom` instance)");
	// clang-format on
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ig2_Tetra_Tetra_TTetraSimpleGeom);


class Law2_TTetraSimpleGeom_NormPhys_Simple : public LawFunctor {
public:
	bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I) override;
	// clang-format off
	YADE_CLASS_BASE_DOC(Law2_TTetraSimpleGeom_NormPhys_Simple,LawFunctor,"EXPERIMENTAL. TODO");
	// clang-format on
	FUNCTOR2D(TTetraSimpleGeom, NormPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_TTetraSimpleGeom_NormPhys_Simple);
#endif


// Miscillaneous functions
//! Tetrahedron's volume.
/// http://en.wikipedia.org/wiki/Tetrahedron#Surface_area_and_volume
Real TetrahedronSignedVolume(const Vector3r v[4]);
Real TetrahedronVolume(const Vector3r v[4]);
Real TetrahedronSignedVolume(const vector<Vector3r>& v);
Real TetrahedronVolume(const vector<Vector3r>& v);
#ifdef YADE_CGAL
Real TetrahedronVolume(const CGAL::Point_3<CGAL::Cartesian<Real>>* v[4]);
Real TetrahedronVolume(const CGAL::Point_3<CGAL::Cartesian<Real>> v[4]);
#endif
Matrix3r TetrahedronInertiaTensor(const vector<Vector3r>& v);
//Matrix3r TetrahedronInertiaTensor(const Vector3r v[4]);
Matrix3r TetrahedronCentralInertiaTensor(const vector<Vector3r>& v);
//Matrix3r TetrahedronCentralInertiaTensor(const Vector3r v[4]);
Quaternionr TetrahedronWithLocalAxesPrincipal(shared_ptr<Body>& tetraBody);

} // namespace yade
