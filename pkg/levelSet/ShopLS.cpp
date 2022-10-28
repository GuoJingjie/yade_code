/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <lib/high-precision/Constants.hpp>
#include <core/Scene.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/levelSet/FastMarchingMethod.hpp>
#include <pkg/levelSet/ShopLS.hpp>

namespace yade { // Cannot have #include directive inside.

CREATE_LOGGER(ShopLS);

// Reminder: we re dealing with static functions here, there would be no point to const-qualify functions themselves

Real ShopLS::biInterpolate(std::array<Real, 2> pt, std::array<Real, 2> xExtr, std::array<Real, 2> yExtr, std::array<std::array<Real, 2>, 2> knownVal)
{
	// Performs interpolation in a 2D space that is denoted (x,y) just for the purpose of the present function, with
	//- pt the point where we want to know the value through interpolation
	//- knownVal, known values at (x0,y0), (x1,y0), (x0,y1), (x1,y1) with eg knownVal[0][1] = value at (x0,y1)
	//- xExtr = (x0,x1) and yExtr = (y0,y1)
	Real x0(xExtr[0]), y0(yExtr[0]);
	Real gx(xExtr[1] - x0), gy(yExtr[1] - y0); // spacings (equal in current implementation of RegularGrid, but let's make things general here)
	Real f00(knownVal[0][0]), f01(knownVal[0][1]), f10(knownVal[1][0]), f11(knownVal[1][1]);
	Real bracket((pt[1] - y0) / gy * (f11 - f10 - f01 + f00) + f10 - f00);
	return (pt[0] - x0) / gx * bracket + (pt[1] - y0) / gy * (f01 - f00) + f00;
}

Vector3r ShopLS::rigidMapping(const Vector3r& ptToMap, const Vector3r& baryIni, const Vector3r& baryEnd, const Quaternionr& rotIniToEnd)
{
	// describes a rigid body transformation of ptToMap, with baryIni/End the initial/final position of the body's center
	// and rotIniToEnd describing the body's rotation from that initial to this final configuration
	return baryEnd + rotIniToEnd * (ptToMap - baryIni);
}

int ShopLS::nGPr(Real min, Real max, Real step) { return nGPv(Vector3r(min, min, min), Vector3r(max, max, max), step)[0]; }

vector<vector<vector<Real>>> ShopLS::phiIni(int usage, Vector3r radii, Vector2r epsilons, shared_ptr<Clump> clump, shared_ptr<RegularGrid> grid)
{
	//- usage = 0: calling a ioFn Python Function => phiIniCppPy
	//- usage = 1: superellipsoids from radii and epsilons => distIniSE
	//- usage = 2: clump => distIniClump
	vector<vector<vector<Real>>> phiField;
	int                          nGPx(grid->nGP[0]), nGPy(grid->nGP[1]), nGPz(grid->nGP[2]); // how many grid pts we re working with
	// reserving memory for phiField and ingredients:
	phiField.reserve(nGPx);
	vector<vector<Real>> distanceInPlane;
	distanceInPlane.reserve(nGPy);
	vector<Real> distanceAlongZaxis;
	distanceAlongZaxis.reserve(nGPz);
	Real                  io(NaN); // return value of the appropriate (Python, SE, Clump) inside outside function at some gridpoint
	boost::python::object main = boost::python::import("__main__"); // useful only if usage == 0 but need to be declared here
	Vector3r              gp;
	for (int xInd = 0; xInd < nGPx; xInd++) {
		distanceInPlane.clear(); // https://stackoverflow.com/questions/19189699/clearing-a-vector-or-defining-a-new-vector-which-one-is-faster
		for (int yInd = 0; yInd < nGPy; yInd++) {
			distanceAlongZaxis.clear();
			for (int zInd = 0; zInd < nGPz; zInd++) {
				gp = grid->gridPoint(xInd, yInd, zInd);
				if (usage == 0) {
					io = boost::python::extract<Real>(main.attr("ioFn")(gp[0], gp[1], gp[2]));
					LOG_DEBUG("Getting " << io << " as return value from Python ioFn");
				} else if (usage == 1)
					io = insideOutsideSE(gp, radii, epsilons);
				else if (usage == 2)
					io = (insideClump(gp, clump) ? -1 : 1);
				if (io > 0) // outside
					distanceAlongZaxis.push_back(std::numeric_limits<Real>::infinity());
				else if (io < 0) // inside
					distanceAlongZaxis.push_back(-std::numeric_limits<Real>::infinity());
				else // io = 0
					distanceAlongZaxis.push_back(0);
			}
			distanceInPlane.push_back(distanceAlongZaxis);
		}
		phiField.push_back(distanceInPlane);
	}
	Real otherVal; // will be a neighboring (or the same, in edge cases) phi value in the grid
	for (int exterior = 0; exterior < 2; exterior++) {
		for (int xInd = 0; xInd < nGPx; xInd++) {
			for (int yInd = 0; yInd < nGPy; yInd++) {
				for (int zInd = 0; zInd < nGPz; zInd++) {
					// let s first check whether we are on the expected (exterior?) side, and continue if not:
					if ((phiField[xInd][yInd][zInd] < 0 && exterior) || (phiField[xInd][yInd][zInd] > 0 && !exterior))
						continue; // lower precedence of && vs < OK
					// taking the distToClumpSpheres value as granted for the exterior (and we do not build a corresponding narrow band for that side):
					if (usage == 2 && exterior) {
						phiField[xInd][yInd][zInd] = distToClumpSpheres(grid->gridPoint(xInd, yInd, zInd), clump);
						continue; // nothing more to be done for that gp, moving directly to the next one
					}
					// comparing now with neighbors, to detect the initial narrow band:
					for (unsigned int neigh = 0; neigh < 6; neigh++) {
						switch (neigh) {
							case 0: // neighbor towards x-, if possible
								otherVal = (xInd == 0 ? phiField[xInd][yInd][zInd] : phiField[xInd - 1][yInd][zInd]);
								break;
							case 1: // neighbor towards x+, if possible
								otherVal = (xInd == nGPx - 1 ? phiField[xInd][yInd][zInd] : phiField[xInd + 1][yInd][zInd]);
								break;
							case 2: // towards y-
								otherVal = (yInd == 0 ? phiField[xInd][yInd][zInd] : phiField[xInd][yInd - 1][zInd]);
								break;
							case 3: // towards y+
								otherVal = (yInd == nGPy - 1 ? phiField[xInd][yInd][zInd] : phiField[xInd][yInd + 1][zInd]);
								break;
							case 4: // towards z-
								otherVal = (zInd == 0 ? phiField[xInd][yInd][zInd] : phiField[xInd][yInd][zInd - 1]);
								break;
							case 5: // towards z+
								otherVal = (zInd == nGPz - 1 ? phiField[xInd][yInd][zInd] : phiField[xInd][yInd][zInd + 1]);
								break;
							default: break;
						}
						if ((otherVal < 0 && exterior)
						    || (otherVal > 0
						        && !exterior)) { // we are next to an inside gridpoint, ie in the initial front serving as BC
							gp = grid->gridPoint(xInd, yInd, zInd);
							if (usage == 0) {
								phiField[xInd][yInd][zInd]
								        = boost::python::extract<Real>(main.attr("ioFn")(gp[0], gp[1], gp[2]));
							} else if (usage == 1)
								phiField[xInd][yInd][zInd] = distApproxSE(gp, radii, epsilons);
							else if (
							        usage == 2
							        && !exterior) // inner part of a clump. Outer was done systematically (not only in the narrow band) above
								phiField[xInd][yInd][zInd] = distToClumpSpheres(gp, clump);
							if ((phiField[xInd][yInd][zInd] < 0 && exterior) || (phiField[xInd][yInd][zInd] > 0 && !exterior))
								LOG_ERROR(
								        "Not on the good side ! At "
								        << xInd << " " << yInd << " " << zInd << " since phi = " << phiField[xInd][yInd][zInd]
								        << " while we re looking at the " << (exterior ? "exterior" : "interior"));
							break; // no need to look at other neighbours along other axes, let s move on to the next gp
						}
					}
				}
			}
		} // closing the gp loop
	}         // closing the side loop
	return phiField;
}

vector<vector<vector<Real>>> ShopLS::phiIniCppPy(shared_ptr<RegularGrid> grid)
{
	// will give phiIni from an existing ioFn Python function, calling phiIni(usage=0, ..)
	// clang-format off
	shared_ptr<Clump> clump(new Clump); // we still need to pass an useless clump below. We can not use default arguments because they need to be trailing ones in C++
	// clang-format on
	return phiIni(0, Vector3r::Zero(), Vector2r::Zero(), clump, grid);
}

Vector3i ShopLS::nGPv(const Vector3r& min, const Vector3r& max, const Real& step)
{
	bool problem = false;
	for (int i = 0; i < 3; i++) {
		if (min[i] >= max[i]) problem = (problem || true);
	}
	if (problem) LOG_ERROR("min wrongly defined as >= max");
	int nGPx(int(ceil((max[0] - min[0]) / step)) + 1) // this number of grid points along the x-axis
	        ,
	        nGPy(int(ceil((max[1] - min[1]) / step)) + 1), nGPz(int(ceil((max[2] - min[2]) / step)) + 1);
	return Vector3i(nGPx, nGPy, nGPz);
}

vector<vector<vector<Real>>> ShopLS::distIniClump(shared_ptr<Clump> clump, shared_ptr<RegularGrid> grid)
{
	return phiIni(2, Vector3r::Zero(), Vector2r::Zero(), clump, grid);
}

vector<vector<vector<Real>>> ShopLS::distIniSE(const Vector3r& radii, const Vector2r& epsilons, shared_ptr<RegularGrid> grid)
{
	// clang-format off
	shared_ptr<Clump> clump(new Clump); // we still need to pass an useless clump below. We can not use default arguments because they need to be trailing ones in C++
	// clang-format on
	return phiIni(1, radii, epsilons, clump, grid);
}

shared_ptr<LevelSet> ShopLS::lsSimpleShape(int shape, const AlignedBox3r& aabb, const Real& step, const Real& smearCoeff, const Vector2r& epsilons, shared_ptr<Clump> clump)
{
	Vector3r minBod(aabb.min()), maxBod(aabb.max()), dimAabb(maxBod - minBod);
	if ((dimAabb[0] < 0) || (dimAabb[1] < 0) || (dimAabb[2] < 0)) LOG_ERROR("You specified negative extents for aabb, this is not expected.");

	shared_ptr<LevelSet> lsShape(new LevelSet);
	lsShape->smearCoeff = smearCoeff;
	Vector3r             maxGrid,
	        minGrid; // [minGrid,maxGrid] will be LevelSet->lsGrid. It may be different (bigger) than [minBod,maxBod] for the grid to go a little beyond the surface
	                 // clang-format off
	// ***** 1. We first define the grid (and twoD) depending upon the shape chosen *******
	if (shape == 0) { // disk (in x,y plane)
		         // clang-format on
		maxGrid       = Vector3r(maxBod[0] + step, maxBod[0] + step, step / 2.); // maxBod is Vector3(radius,radius,radius)
		minGrid       = -maxGrid;
		lsShape->twoD = true;
	} else if (shape <= 4) { // sphere, box, superellipsoid, clump of spherical particles, all combined
//		NB: pay attention to the following line, it has to enable the origin to belong to the grid (on all 3 axes)
		int nIntX(int(ceil(maxBod[0]/step)) + 1 ), nIntY( int(ceil(maxBod[1]/step)) + 1), nIntZ( int(ceil(maxBod[2]/step)) + 1);
		maxGrid       = Vector3r(nIntX*step,nIntY*step,nIntZ*step);
		minGrid       = -maxGrid;
		lsShape->twoD = false;
	} else
		LOG_FATAL("You asked for some shape value=" << shape << " which is not supported");
	Vector3i gpPerAxis = nGPv(minGrid, maxGrid, step);
	if (lsShape->twoD && gpPerAxis[2] != 2) LOG_WARN("Ending with " << gpPerAxis[2] << " gridpoints along z axis, instead of 2. This is not optimal.");
	lsShape->lsGrid->min     = minGrid;
	lsShape->lsGrid->spacing = step;
	lsShape->lsGrid->nGP     = gpPerAxis;
	int nGP_x(gpPerAxis[0]), nGP_y(gpPerAxis[1]), nGP_z(gpPerAxis[2]);

	// ***** 2. We now work on the distance field itself *******
	vector<vector<vector<Real>>> distanceVal;
	vector<vector<Real>>         distanceInPlane;
	vector<Real>                 distanceAlongZaxis;
	Real                         dist(std::numeric_limits<Real>::infinity());
	Vector3r                     gridPt;
	if (shape < 3) { // looping on all gridpoints to apply appropriate distance function when having an analytical function
		for (int xInd = 0; xInd < nGP_x; xInd++) {
			distanceInPlane.clear(); // https://stackoverflow.com/questions/19189699/clearing-a-vector-or-defining-a-new-vector-which-one-is-faster
			for (int yInd = 0; yInd < nGP_y; yInd++) {
				distanceAlongZaxis.clear();
				for (int zInd = 0; zInd < nGP_z; zInd++) {
					gridPt = lsShape->lsGrid->gridPoint(xInd, yInd, zInd);
					switch (shape) {
						case 0: // circles
							dist = distToCircle(gridPt, maxBod[0]);
							break;
						case 1: // sphere
							dist = distToSph(gridPt, maxBod[0]);
							break;
						case 2: // box
							dist = distToRecParallelepiped(gridPt, (maxBod - minBod) / 2.);
							break;
						default: LOG_FATAL("You asked for some shape value=" << shape << " which is not supported");
					}
					distanceAlongZaxis.push_back(dist);
				}
				distanceInPlane.push_back(distanceAlongZaxis);
			}
			distanceVal.push_back(distanceInPlane);
		}
	} else if (shape == 3) // FMM for se, with the whole distance field at once
	{
		FastMarchingMethod distFMM;
		distFMM.grid   = lsShape->lsGrid;
		distFMM.phiIni = distIniSE((maxBod - minBod) / 2., epsilons, lsShape->lsGrid);
		distanceVal    = distFMM.phi();
	} else if (shape == 4) {
		vector<int>        ids(clump->ids);
		FastMarchingMethod distFMM;
		distFMM.grid   = lsShape->lsGrid;
		distFMM.phiIni = distIniClump(clump, lsShape->lsGrid);
		distanceVal    = distFMM.phi();
	}
	if (!distanceVal.size()) // eg because there was no level set grid
		LOG_ERROR("We computed an empty level set... This is not expected and may crash e.g. when exposed to Python.")
	lsShape->distField = distanceVal;
	return lsShape;
}

Real ShopLS::distToSph(const Vector3r& point, const Real& rad, Vector3r center /*=Vector3r::Zero()*/)
{
	Real x(point[0]), y(point[1]), z(point[2]);
	Real xC(center[0]), yC(center[1]), zC(center[2]);
	return sqrt(pow(x - xC, 2) + pow(y - yC, 2) + pow(z - zC, 2)) - rad;
}

Real ShopLS::distToCircle(const Vector3r& point, const Real& rad)
{
	// for a circle in (x,y) plane
	Real x(point[0]), y(point[1]);
	return sqrt(pow(x, 2) + pow(y, 2)) - rad;
}

Real ShopLS::distToRecParallelepiped(Vector3r point, Vector3r extents)
{
	// Computing true outside distances is complicated: the homothetic squares, suggested by Fig. 2b Kawamoto2016, have to be rounded at their corners for an exact consideration of distances. We will neglect this rounding here, which is not a big deal since positive values of level set functions are anyway transparent for the DEM
	// in this case, the distance between the cube surface and the point considered is just the maximum of the positive values of distances, ie the maximum and there is no difference with inside case...
	Real x(point[0]), y(point[1]), z(point[2]);
	Real distX(distToInterval(x, -extents[0], extents[0]));
	Real distY(distToInterval(y, -extents[1], extents[1]));
	Real distZ(distToInterval(z, -extents[2], extents[2]));
	return math::max(math::max(distX, distY), distZ);
}

Real ShopLS::distToInterval(Real val, Real left, Real right)
{
	Real dist(-1.);
	if (val > left and val < right) dist = math::max(left - val, val - right);
	else if (val == right)
		dist = 0;
	else if (val > right)
		dist = val - right; // > 0 we're indeed outside
	else if (val == left)
		dist = 0;
	else if (val < left)
		dist = left - val; // > 0 we're indeed outside
	return dist;
}

Real ShopLS::fioRose(Vector3r gp)
{
	gp = ShopLS::cart2spher(gp);
	Real r(gp[0]), theta(gp[1]), phi(gp[2]);
	return r - 3 - 1.5 * sin(5 * theta) * sin(4 * phi); // value of inside/outside function
}

Vector3r ShopLS::grad_fioRose(Vector3r gp)
{
	gp = ShopLS::cart2spher(gp);
	Real r(gp[0]), theta(gp[1]), phi(gp[2]);
	if (sin(theta) == 0) LOG_ERROR("theta = 0 [pi], gradient of rose fction not defined for its z component");
	Vector3r grad_fio(1, -7.5 / r * cos(5 * theta) * sin(4 * phi), -6 / r * sin(5 * theta) / sin(theta) * cos(4 * phi));
	return grad_fio;
}
shared_ptr<ScGeom> ShopLS::geomPtrForLaterRemoval(const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c)
{
	// to use when we can not really compute anything, e.g. bodies lsGrid do not overlap anymore, but still need to have some geom data (while returning true as per general InteractionLoop workflow because it is an existing interaction. Otherwise we would need to update InteractionLoop itself to avoid LOG_WARN messages). Data mostly include an infinite tensile stretch to insure subsequent interaction removal (by Law2)
	return ShopLS::geomPtr(
	        Vector3r::Zero() /* inconsequential bullsh..*/,
	        -std::numeric_limits<Real>::infinity() /* arbitrary big tensile value to trigger interaction removal by Law2*/,
	        1, /* inconsequential bullsh..*/
	        1, /* inconsequential bullsh..*/
	        rbp1,
	        rbp2,
	        c,
	        Vector3r::UnitX() /* inconsequential bullsh..*/,
	        Vector3r::Zero() /* inconsequential bullsh..*/);
}

shared_ptr<ScGeom> ShopLS::geomPtr(
        Vector3r ctctPt, Real un, Real rad1, Real rad2, const State& rbp1, const State& rbp2, const shared_ptr<Interaction>& c, const Vector3r& currentNormal, const Vector3r& shift2)
{
	shared_ptr<ScGeom> geomPtr;
	bool               isNew = !c->geom;
	if (isNew) geomPtr = shared_ptr<ScGeom>(new ScGeom());
	else
		geomPtr = YADE_PTR_CAST<ScGeom>(c->geom);
	geomPtr->contactPoint     = ctctPt;
	geomPtr->penetrationDepth = un;
	// NB radius1, radius2: those are useful for
	// 1* contact kinematics description if and only if avoidGranularRatcheting, (not the case here)
	// 2* applying contact forces in C-S Law2, if sphericalBodies (not the case here)
	// 3* time step determination with respect to rotational stiffnesses
	// and also, as refR1, refR2, for
	// 4* contact stiffness expression in FrictPhys/FrictMat
	geomPtr->radius1 = rad1;
	geomPtr->radius2 = rad2;
	geomPtr->precompute(
	        rbp1,
	        rbp2,
	        Omega::instance().getScene().get(),
	        c,
	        currentNormal,
	        isNew,
	        shift2,
	        false); // the avoidGranularRatcheting=1 expression of relative velocity at contact does not make sense with radius1 and radius2 in all cases of LevelSet-sthg interaction (because the shapes are non-spherical): our present radius1+radius2 may have nothing in common with branch vector
	// precompute will take care of
	// * preparing the rotation of shearForce to the new tangent plane (done later, in Law2) defining these orthonormal_axis and twist_axis
	// * updating geomPtr->normal to normal
	// * computing the relative velocity at contact, through getIncidentVel(avoidGranularRatcheting=false), using now-defined contactPoint

	// Comparing with Ig2_Sphere_Sphere_ScGeom.cpp, I think everything is here
	return geomPtr;
}

Real ShopLS::distApproxRose(Vector3r gp)
{
	Vector3r grad_fio(grad_fioRose(gp));
	Real     normGrad(grad_fio.norm());
	if (normGrad == 0) LOG_ERROR("Zero gradient, approximate distance will be infinite");
	return fioRose(gp) / normGrad;
}

Real ShopLS::distToClumpSpheres(Vector3r pt, shared_ptr<Clump> clump)
{
	// gives min(dist to clump members) for spherical members
	// This is the true distance when outside, and just an approximation when inside
	// It will serve in phiIni(usage==2) and has some duplicated code with insideClump, let it be (insideClump has a break we do not have here for instance)
	Real retVal(std::numeric_limits<Real>::infinity());
	for (const auto& memberInfos : clump->members) { // clump->members is a std::map
		const Body::id_t& memberId(memberInfos.first);
		Vector3r          centerMember(memberInfos.second.position);
		Real              radiusMember(YADE_PTR_CAST<Sphere>(Body::byId(memberId)->shape)->radius);
		Real              dist(distToSph(pt, radiusMember, centerMember));
		if (dist < retVal) retVal = dist;
	}
	return retVal;
}

bool ShopLS::insideClump(const Vector3r& pt, shared_ptr<Clump> clump)
{ // checking whether we are strictly inside a Clump. Could fit into core/Clump ?
	bool retVal(false);
	for (const auto& memberInfos : clump->members) { // clump->members is a std::map
		const Body::id_t& memberId(memberInfos.first);
		Vector3r          centerMember(memberInfos.second.position);
		Real              radiusMember(YADE_PTR_CAST<Sphere>(Body::byId(memberId)->shape)->radius);
		LOG_INFO(
		        "Looking at Body " << memberId << " as a clump member. It has a " << radiusMember
		                           << " radius and is centered (in clump local coordinates) at " << centerMember);
		Real dist(distToSph(pt, radiusMember, centerMember));
		LOG_DEBUG("dist = " << dist);
		if (dist
		    < 0) { // dist = 0 actually tells nothing about the clump surface (unless we check we are outside every other member), and can not be used
			retVal = true;
			break; // pt is inside at least one clump member => it's inside the clump
		}
	}
	return retVal;
}

Real ShopLS::insideOutsideSE(const Vector3r& point, const Vector3r& radii, const Vector2r& epsilons)
{
	Real x(point[0]), y(point[1]), z(point[2]), rx(radii[0]), ry(radii[1]), rz(radii[2]), epsE = epsilons[0], epsN = epsilons[1];
	if (rx < 0 || ry < 0 || rz < 0) LOG_ERROR("You passed negative radii for a superellipsoid, this is not expected.");
	return pow(pow(math::abs(x / rx), 2. / epsE) + pow(math::abs(y / ry), 2. / epsE), epsE / epsN) + pow(math::abs(z / rz), 2. / epsN) - 1;
}

Real ShopLS::distApproxSE(const Vector3r& point, const Vector3r& radii, const Vector2r& epsilons)
{
	// Approximated distance for a superellipsoid, using the inside/outside fxyz divided by the norm of its gradient, grad_f.norm()
	// with a special treatment of the origin because grad_f.norm() -> 0 towards the origin (where it's not defined ?), whereas fxyz -> -1 and our approximated distance would -> infinity
	Real x(point[0]), y(point[1]), z(point[2]), rx(radii[0]), ry(radii[1]), rz(radii[2]), epsE = epsilons[0], epsN = epsilons[1];
	Real returnVal, minExtents(math::min(math::min(rx, ry), rz));
	if (point.norm() < 0.05 * minExtents) {
		returnVal = -minExtents; // using (signed) minExtents as a constant distance-approximation in that whole domain close to origin
		LOG_WARN("We used a crude distance approximation at point " << point << "! That could be a problem e.g. when using FastMarchingMethod");
	} else {
		Real fxyz(insideOutsideSE(point, radii, epsilons));
		if (math::abs(fxyz) < Mathr::EPSILON) // it s like we were exactly on the surface
			returnVal = 0.;
		else {
			if (x == 0 and y == 0) // precedence OK
				return z > 0 ? z - rz : -z - rz;
			Real mult(
			        2. / epsN
			        * pow(pow(math::abs(x / rx), 2. / epsE) + pow(math::abs(y / ry), 2. / epsE),
			              epsE / epsN - 1.)); //some number useful for grad f.
			if (!math::isfinite(mult))        // for x=y=0 and epsE < epsN typically, but this should have been taken care of in the above
				LOG_ERROR("Infinite (or C++ NaN) mult in distance to a super ellipsoid");
			Vector3r grad_f(
			        mult * pow(math::abs(x / rx), 2. / epsE - 1.) * (x > 0 ? 1. / rx : -1. / rx),
			        mult * pow(math::abs(y / ry), 2. / epsE - 1.) * (y > 0 ? 1. / ry : -1. / ry),
			        2. / epsN * pow(math::abs(z / rz), 2. / epsN - 1) * (z > 0 ? 1. / rz : -1. / rz));
			// special cases where gradient is actually not defined:
			// if(x==0 && epsE>=2) grad_f[0] = 0; // like when epsE < 2, and normal for a sphere
			// if(y==0 && epsE>=2) grad_f[1] = 0; // ditto
			// if(z==0 && epsN>=2) grad_f[2] = 0; // same, but for epsN
			Real denom(0);
			if (grad_f.norm()
			    < minExtents * Mathr::EPSILON) { // checking whether the gradient is not too small, with respect to a characteristic length.
				LOG_ERROR("Unexpected zero gradient in point = " << point);
				denom = 1;
			} // no real reason for this 1
			else
				denom = grad_f.norm();
			returnVal = fxyz / denom;
			LOG_DEBUG(
			        "Here returnVal = " << returnVal << " because fxyz = " << fxyz << " and denom = " << denom << " because grad_f = " << grad_f
			                            << " because mult = " << mult);
		}
	}
	return returnVal;
}

// The 2 following cart2spher and spher2cart functions assume:
// theta = (ez,er) angle with er(theta = 0,phi)=ez, in [0;pi]
// phi measured from ex axis: er(theta=pi/2,phi=0) = ex, in [0;2pi]
// same convention as Bazant1986, Eq. (4)

Vector3r ShopLS::cart2spher(const Vector3r& vec)
{ //	vec is here (x,y,z), we return (r,theta,phi)
	Real r(vec.norm());
	if (r == 0) return Vector3r::Zero();
	Real     theta(acos(vec[2] / r));
	Vector3r vecProj(vec[0], vec[1], 0); // projection of vec in x,y plane
	Real     normProj(vecProj.norm());
	Real     phi;
	if (normProj == 0) // vec is along z axis, phi=0 will be as good as any other value
		phi = 0;
	else {
		Real cosVal(vecProj[0] / normProj);
		phi = vec[1] > 0 ? acos(cosVal) : 2 * Mathr::PI - acos(cosVal);
	}
	return Vector3r(r, theta, phi);
}

Vector3r ShopLS::spher2cart(const Vector3r& vec)
{ //	vec is here (r,theta,phi), and we return (x,y,z)
	if (vec[0] < 0) LOG_ERROR("A negative r (" << vec[0] << " passed here) for spherical coordinates is impossible");
	if (vec[1] < 0 || vec[1] > Mathr::PI) LOG_ERROR("Spherical theta has to be between 0 and pi, passing " << vec[1] << " is impossible")
	return vec[0] * Vector3r(sin(vec[1]) * cos(vec[2]), sin(vec[1]) * sin(vec[2]), cos(vec[1]));
}

} // namespace yade
#endif // YADE_LS_DEM
