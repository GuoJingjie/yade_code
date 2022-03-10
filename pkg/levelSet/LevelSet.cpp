/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <lib/computational-geometry/MarchingCube.hpp>
#include <lib/high-precision/Constants.hpp>
#include <core/Scene.hpp>
#include <pkg/levelSet/LevelSet.hpp>
#include <pkg/levelSet/ShopLS.hpp>
#include <boost/math/tools/roots.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade {
YADE_PLUGIN((LevelSet));
CREATE_LOGGER(LevelSet);

Vector3r LevelSet::getCenter()
{
	if (!initDone) init();
	return center;
}

// // vector<int> LevelSet::getNodesInCellCube(int i, int j, int k){ Vector3i ?
// // 	return nodesIndicesInCell[Vector3i(i,j,k)];
// // }
void LevelSet::rayTrace(const Vector3r& ray)
{
	// we launch a ray from shape's center (= 0), and go through all grid cells until the end of the grid, looking for ray-surface intersection
	Vector3i indices(lsGrid->closestCorner(Vector3r::Zero())), // which grid cell do we start from ?
	        move;                                              // to control the move from one cell to another
	Vector3r pointP(Vector3r::Zero()),
	        point0; // pointP will always be the point where the ray starts from, while point0 is the lowest corner of the current grid cell
	std::array<Real, 3>
	        kVal; // points of the ray are of the form pointP + k * ray (k>= 0). Then, this kVal will include the 3 smallest k-values that make the ray exit of the cell, considering in an independent manner the three axis. E.G. for a ray = (1,0,0) starting from point0, kVal = (spacing, infinity, infinity)
	bool    diffSign // true if and only if at least one gp has a different sign than the others, ie iff at least one gp has a diff sign than the 1st one
	        ,
	        sign(false) // binary sign of current grid point
	        ,
	        touched(false) // true iff ray touches at least once the surface (at least one node found)
	        ,
	        touchedInCell(false) // true iff ray touches the surface in a given cell (one node found in that cell)
	        ;
	Real         minOtherAxes(-1.), gpDist;
	unsigned int xInd, yInd, zInd; // the grid indices of cell gridpoints, used below
	while (true) {
		LOG_DEBUG(
		        "\nConsidering cell " << indices[0] << " " << indices[1] << " " << indices[2] << " "
		                              << " while grid itself has " << lsGrid->nGP[0] << " " << lsGrid->nGP[1] << " " << lsGrid->nGP[2] << "gridpoints");
		point0   = lsGrid->gridPoint(indices[0], indices[1], indices[2]);
		diffSign = false;
		for (unsigned int gp = 0; gp < 8; gp++) { // passing through 8 cell gridpoints to check whether they all have the same distance sign
			xInd   = (gp % 2 ? indices[0] + 1
			                 : indices[0]); // better put parenthesis: "=" and ternary have same precedence, let s not look at associativity
			yInd   = ((gp & 2) / 2 ? indices[1] + 1 : indices[1]);
			zInd   = ((gp & 4) / 4 ? indices[2] + 1 : indices[2]);
			gpDist = distField[xInd][yInd][zInd]; // distance value for current gp
			                                      //			Same sign considerations:
			if (gp == 0) sign = gpDist > 0;       // precedence OK; the binary sign of the "1st" gridpoint (third case = 0 handled otherwise)
			if (gp > 0 and (gpDist > 0) != sign)  // "and" least precedence OK
				diffSign = diffSign || true;  // diffSign is now (and will remain) true because we found one gp with a different sign
				                              // including the case of a zero distance at one cell gp
			if (gpDist == 0) {
				LOG_DEBUG("We should find one gridpoint among the boundary nodes");
				diffSign = diffSign || true; // that s also a case where we need to trigger ray tracing (if not already done)..
			}
		}
		if (diffSign) { // zero distance surface can be detected as in the cell only with at least one change in sign between the gridpoints
			touchedInCell = rayTraceInCell(ray, pointP, point0, indices);
			LOG_INFO(
			        "Considering cell " << indices[0] << " " << indices[1] << " " << indices[2]
			                            << " which should include the surface. rayTraceInCell actually returned " << touchedInCell);
			touched = touched
			        || touchedInCell; // just touched || rayTraceInCell(..) would not work for multiple nodes along a ray because of short circuit !
		}
		// we check whether we still have to walk along ray, or whether we reached the grid boundary (whatever boundary it is):
		if ((indices[0] == lsGrid->nGP[0] - 2) or (indices[1] == lsGrid->nGP[1] - 2) or (!twoD and (indices[2] == lsGrid->nGP[2] - 2)) or !indices[0]
		    or !indices[1] or (!twoD and !indices[2])) {
			if (!touched) LOG_WARN("Ray (" << ray[0] << " " << ray[1] << " " << ray[2] << ") did not create a boundary node");
			LOG_INFO("Ray reached grid boundary because indices = " << indices << ", ray tracing for that one is over");
			break; // this should make sense as soon as center is not already in this boundary layer (which should not happen)
		}
		// we will now move to the next cell:
		// first, computing the path lengths (along ray) to exit our current cell, testing all three directions, and filling above-mentioned kVal
		for (int axis = 0; axis < 3; axis++) {
			if (ray[axis] > 0) kVal[axis] = (point0[axis] + lsGrid->spacing - pointP[axis]) / ray[axis]; // a positive value
			else if (ray[axis] < 0)
				kVal[axis] = (point0[axis] - pointP[axis]) / ray[axis]; // another positive value
			else
				kVal[axis] = std::numeric_limits<Real>::infinity();
		}
		//second, testing which exit is the closest and update move accordingly:
		move = Vector3i::Zero(); // after recalling the initial value here, setting -1,+1,0 in this Vector3i to update the indices below
		for (int axis = 0; axis < 3; axis++) {
			minOtherAxes = math::min(kVal[(axis + 1) % 3], kVal[(axis + 2) % 3]);
			// do we exit along axis ?
			if (kVal[axis] - minOtherAxes < Mathr::EPSILON * lsGrid->spacing) {
				//				yes: exit along axis is the quickest.
				//				We also take into account the possibility of kVal[axis] being very close from minOtherAxes (see below), be it by excess or by default
				move[axis] = int(math::sign(ray[axis]));
				for (int j = 1; j < 3; j++) { // still checking the two other axes because of finite precision
					if (math::abs(kVal[axis] - kVal[(axis + j) % 3]) < Mathr::EPSILON * lsGrid->spacing)
						move[(axis + j) % 3] = int(math::sign(ray[(axis + j) % 3])); // the exit along axis+j is actually just as close
				}
			}
		}
		// third, updating pointP, as well as cell indices for the next step in the loop (see in the above)
		pointP = pointP + *(std::min_element(kVal.begin(), kVal.end())) * ray;
		indices += move;
		// just checking finally whether we indeed changed cell:
		if (move.squaredNorm() == 0) {
			LOG_ERROR("We're stuck in the same cell !!!");
			break;
		}
	}
}

bool LevelSet::rayTraceInCell(const Vector3r& ray, const Vector3r& pointP, const Vector3r& point0, const Vector3i& indices)
{
	// Given the level set surface described by a trilinear distance function and an oriented ray starting from pointP in some cell,
	// where the cell is (redundantly, but this is to avoid computing things already computed in rayNode) defined by point0 (lowest corner) and its indices
	// it searches for the surface points intersected by the ray in that cell (=  surfNodes), and returns true if it found one
	int indX(indices[0]), indY(indices[1]), indZ(indices[2]);
	LOG_DEBUG("Ray tracing from pointP = " << pointP << " in cell " << indX << " " << indY << " " << indZ);
	Real     normNode(-1);
	Vector3r trialNode;
	Real     xP(pointP[0]), yP(pointP[1]), zP(pointP[2]);
	Real     ux(ray[0]), uy(ray[1]), uz(ray[2]); // the ray direction (oriented: vec is not the same as -vec)
	Real     spac = lsGrid->spacing;
	bool     touched(false);
	Real     x0(point0[0]), y0(point0[1]), z0(point0[2]);
	Real     f000(distField[indX][indY][indZ]), f111(distField[indX + 1][indY + 1][indZ + 1]), f100(distField[indX + 1][indY][indZ]),
	        f010(distField[indX][indY + 1][indZ]), f001(distField[indX][indY][indZ + 1]), f101(distField[indX + 1][indY][indZ + 1]),
	        f011(distField[indX][indY + 1][indZ + 1]), f110(distField[indX + 1][indY + 1][indZ]);
	// the coefficients of the trilinear interpolation, when written f(x,y,z) = A*x'y'z' + B*x'y' + C*y'z' + D*x'z' + E*x' + F*y' + G*z' + H (with x', .. the dimensionless coordinates, and the useless H = f000):
	Real A(f111 + f100 + f010 + f001 - f101 - f110 - f011 - f000), B(f110 - f100 - f010 + f000), C(f011 - f010 - f001 + f000), D(f101 - f100 - f001 + f000),
	        E(f100 - f000), F(f010 - f000), G(f001 - f000);
	// for the cubic polynom expression of that interpolation (when applied to a ray), we need:
	Real Ag3(A / pow(spac, 3)), Bg2(B / pow(spac, 2)), Cg2(C / pow(spac, 2)), Dg2(D / pow(spac, 2));
	// hence the polynom, in a dimensionless form:
	std::array<Real, 4> coeffs {
		{ distance(pointP) / spac,
		  Ag3 * (uz * (xP - x0) * (yP - y0) + uy * (xP - x0) * (zP - z0) + ux * (yP - y0) * (zP - z0)) + Bg2 * (ux * (yP - y0) + uy * (xP - x0))
		          + Cg2 * (uy * (zP - z0) + uz * (yP - y0)) + Dg2 * (ux * (zP - z0) + uz * (xP - x0)) + E / spac * ux + F / spac * uy + G / spac * uz,
		  (Ag3 * ((xP - x0) * uy * uz + (yP - y0) * ux * uz + (zP - z0) * ux * uy) + Bg2 * ux * uy + Cg2 * uy * uz + Dg2 * ux * uz) * spac,
		  (Ag3 * ux * uy * uz) * pow(spac, 2) }
	};
	LOG_INFO(
	        "We think finding surface-ray intersection is like finding the roots of this polynom (coeffs of increasing order):\n"
	        << coeffs[0] << " " << coeffs[1] << " " << coeffs[2] << " " << coeffs[3]);
	// Then, solving for the 1st root given by Boost root finding algorithm boost::math::tools::newton_raphson_iterate, in a -2 line syntax thanks to this C++11 lambda notation, that replaces going through e.g. an operator() of some struct (whose instance may depend on coeffs)
	Real root(boost::math::tools::newton_raphson_iterate(
	        [coeffs](Real k) {
		        return std::make_tuple( // (only ?) problematic part for High Precision (HP) compatibility and REAL_* cmake options different than 64
		                coeffs[0] + coeffs[1] * k + coeffs[2] * k * k + coeffs[3] * k * k * k,
		                coeffs[1] + 2 * coeffs[2] * k + 3 * coeffs[3] * k * k);
	        },
	        0.,
	        -sqrt(3.),
	        sqrt(3.),
	        std::numeric_limits<Real>::
	                digits)); // doc https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/roots_deriv.html suggests to use 0.6 * std::numeric_limits<Real>::digits for the last "digits" argument. With std::numeric_limits<double>::digits = 53, as per standard for the significand
	LOG_INFO(
	        "N-R root = " << root << " , leading to a dimensionless distance (through the ray cubic polynom) = "
	                      << coeffs[0] + coeffs[1] * root + coeffs[2] * root * root + coeffs[3] * root * root * root
	                      << " , using digits = " << std::numeric_limits<Real>::digits);
	Vector3i indicesPt; // in which cell do we find trialNode (see below) ? (quite useless since Feb 2021 actually, except for log messages)
	if (root >= 0) trialNode = pointP + spac * root * ray;
	else if (
	        math::abs(root) <= Mathr::
	                EPSILON) // tiny negative numbers count as well. See e.g. 'superellipsoid',extents=(0.07,0.09,0.14),epsilons=(1.1,0.2),spacing= 0.01 for the -z ray
		trialNode = pointP;
	else // truly negative roots do not count
		return false;
	// trialNode is not garanteed until now to be in the considered (indices) cell, to which the cubic polynom should restrict. In order to check for consistency, we no longer compare indices with indicesPt (that would be a discontinuous test too painful to make work with numeric precision), and just look at the distance value (which is continuous at the boundary of two cells)
	if (!Shop::isInBB(trialNode, lsGrid->min, lsGrid->max())) return false; // still starting to check whether trialNodes fits into lsGrid
	indicesPt = lsGrid->closestCorner(trialNode);
	LOG_INFO("We have a possible intersection point (" << trialNode << "): it would be in cell " << indicesPt << " whereas we considered cell " << indX
	                                                   << " " << indY << " " << indZ << endl;);
	LOG_INFO(
	        "Dimensionless distance, wrt to grid spacing, is here (from distance function) = " << distance(trialNode)
	                / spac << " ; while nodesTol*epsilon = " << nodesTol * Mathr::EPSILON);
	if (surfNodes.size() > 0)
		LOG_INFO(
		        "And previous node is " << surfNodes.back() << " whose comparison with gives (in absolute norm)"
		                                << (trialNode - surfNodes.back()).norm() / (surfNodes.back()).norm());
	if ((math::abs(distance(trialNode)) / lengthChar < nodesTol * Mathr::EPSILON)
	    && // looks like we found a node. Root finding algorithm normally had a Mathr::EPSILON precision wrt to spac but distance() computation is not the same, hence a different precision. Which is finally estimated here from lengthChar (which matters more than spac)
	    (surfNodes.size() == 0
	     || ((trialNode - surfNodes.back()).norm() / surfNodes.back().norm() > lengthChar * nodesTol
	                 * Mathr::
	                         EPSILON)) // but we still do not want a duplicate of last node, that would come from an adjacent cell, with phi=0 on the cell surface(s)
	) {
		LOG_INFO(
		        "ONE NODE FOUND ! (" << trialNode << "): it would be in cell " << indicesPt << " whereas we considered cell " << indX << " " << indY
		                             << " " << indZ);
		surfNodes.push_back(trialNode); // we insert it in shape.surfNodes
		touched  = true;
		normNode = trialNode.norm();
		if (normNode < minRad) minRad = normNode;
		if (normNode > maxRad) maxRad = normNode;
	}
	return touched;
}

Vector3r LevelSet::normal(const Vector3r& pt) const
{
	// returns the normal vector at pt, according to (2) Kawamoto2016
	Vector3i indices = lsGrid->closestCorner(pt);
	int      xInd(indices[0]), yInd(indices[1]), zInd(indices[2]);

	if (xInd < 0 || yInd < 0 || zInd < 0) { // operators precedence OK in || vs <
		LOG_ERROR("Can not compute the normal, returning a NaN vector");
		return Vector3r(NaN, NaN, NaN);
	}

	Real     spac   = lsGrid->spacing;
	Vector3r corner = lsGrid->gridPoint(xInd, yInd, zInd);
	Real     xRed((pt[0] - corner[0]) / spac), yRed((pt[1] - corner[1]) / spac),
	        zRed((pt[2] - corner[2]) / spac); // dimensionless x,y,z of 3., top p. 4 Kawamoto2016
	Real nx(0), ny(0), nz(0);                 // the x, y, z components of normal, computed below
	for (int indA = 0; indA < 2; indA++) {
		for (int indB = 0; indB < 2; indB++) {
			for (int indC = 0; indC < 2; indC++) {
				Real lsVal = distField[xInd + indA][yInd + indB][zInd + indC];
				nx += lsVal * (2 * indA - 1) * ((1 - indB) * (1 - yRed) + indB * yRed) * ((1 - indC) * (1 - zRed) + indC * zRed);
				ny += lsVal * (2 * indB - 1) * ((1 - indA) * (1 - xRed) + indA * xRed) * ((1 - indC) * (1 - zRed) + indC * zRed);
				nz += lsVal * (2 * indC - 1) * ((1 - indA) * (1 - xRed) + indA * xRed) * ((1 - indB) * (1 - yRed) + indB * yRed);
			}
		}
	}
	return Vector3r(nx, ny, nz).normalized();
}

void LevelSet::initSurfNodes()
{
	// Obtaining the boundary nodes through "ray tracing", starting from the center
	if (nSurfNodes <= 2 && nodesPath == 1 && !twoD)
		LOG_ERROR("You asked for a level set shape with no more than two boundary nodes, for contact detection purposes. This is too few and will lead "
		          "to square roots of negative numbers, then unexpected events.");
	int  nAngle(int(sqrt(nSurfNodes - 2)));    // useful only if(!twoD && case 1), but let s declare and assign it once for all
	Real phiMult(Mathr::PI * (3. - sqrt(5.))); // useful only if(!twoD && case 2), ditto
	for (int node = 0; node < nSurfNodes; node++) {
		if (twoD) // 2D analysis in (x,y) plane
			rayTrace(Vector3r(cos(float(node) / nSurfNodes * 2 * Mathr::PI), sin(float(node) / nSurfNodes * 2 * Mathr::PI), 0));
		else {
			Real theta(0), phi(0);
			switch (nodesPath) {
				case 1: // rectangular partition method
					if (node < 2) {
						theta = (node == 0 ? 0 : Mathr::PI);
						phi   = 0; // exact value does not really matter..
					} else {
						if ((node - 2) / nAngle >= nAngle) // possible when nSurfNodes != k^2 + 2...
							LOG_ERROR(
							        "Problems may come soon, please define nSurfNodes as a squared integer + 2 if not twoD. "
							        "Otherwise you will get phi = "
							        << (node - 2) / nAngle * 2 * Mathr::PI / nAngle);
						theta = ((node - 2) % nAngle + 1) * Mathr::PI / (nAngle + 1.); // won't include 0 neither pi, handled previously
						phi   = (node - 2) / nAngle * 2. * Mathr::PI
						        / nAngle; // from 0 included to 2 pi excluded. The integer division as first term is intended and prevents merging the 2 nAngle denominators. Also note the left-to-right associativity of * and / which make things work as expected
					}
					break;
				case 2: // spiral points: Elias' code in the spirit of Rakhmanov1994
					theta = acos(-1. + (2. * node + 1.) / nSurfNodes);
					phi   = node * phiMult;
					break;
				default:
					LOG_ERROR(
					        "You asked for locating boundary nodes according to a method #"
					        << nodesPath << " which is not implemented. Supported choices are 1 and 2.")
			}
			LOG_DEBUG(
			        "\nnode = " << node << " ; theta = " << theta << " and phi = " << phi << " ; launching ray " << Vector3r(1, theta, phi)
			                    << endl);
			rayTrace(ShopLS::spher2cart(Vector3r(1, theta, phi)));
		}
	}
	if (int(surfNodes.size()) != nSurfNodes)
		LOG_WARN(
		        "Ray tracing gave " << surfNodes.size() << " boundary nodes instead of " << nSurfNodes
		                            << " asked for: should not happen at least with convex bodies. You may want to adapt nodesTol.");
}

void LevelSet::init() // computes stuff (nVoxInside, center, volume, inertia, boundary nodes, ...) once distField exists
{
	if (initDone) LOG_WARN("How comes we run a second time init ?")
	if (!distField.size()) LOG_ERROR("You are interested into center/volume before that distField has been defined");
	Vector3i gpPerAxis = lsGrid->nGP;
	int      nGPx(gpPerAxis[0]), nGPy(gpPerAxis[1]), nGPz(gpPerAxis[2]);
	if ((int(distField.size()) != nGPx) or (int((distField[0][0]).size()) != nGPz))
		LOG_ERROR("There is a size-inconsistency between the current level set grid and shape.distField for this body! The level set grid has changed "
		          "since the creation of this body, this is not supported.");
	Real spac = lsGrid->spacing;
	Real dV   = pow(spac, 3);
	Real xMean(0), yMean(0), zMean(0); // the center of the level set volume, as determined below
	nVoxInside = 0;
	Vector3r gp;
	// We will base our level set volume description upon the voxellised description using "inside" voxels.
	// A voxel is said to be inside according to its minimum gridpoint only, boundary effects are thus unavoidable for now.
	// To have a consistent voxellised point of view, we will compute the mean coordinates (center's one) from the inside voxels gridpoints. Think to a 2*2 grid (from 0 to 1 along both axes, with just 1 pixel) with a positive level set value at the (0,0) corner only. The only pixel at hand will said to be inside, let's then define the center at the non gridpoint (0.5,0.5). This would not be possible defining the mean from the positive gridpoints (as opposed to the inside voxels considered here).
	for (int xIndex = 0; xIndex < nGPx - 1;
	     xIndex++) { // to nGPx-1 to avoid counting twice the last voxel (or the first if the test would consider the level set value at the max gridpoint)
		for (int yIndex = 0; yIndex < nGPy - 1; yIndex++) {
			for (int zIndex = 0; zIndex < nGPz - 1; zIndex++) {
				if (distField[xIndex][yIndex][zIndex] <= 0) {
					nVoxInside++;
					gp = lsGrid->gridPoint(xIndex, yIndex, zIndex);
					xMean += gp[0] + spac / 2.;
					yMean += gp[1] + spac / 2.;
					zMean += gp[2] + spac / 2.;
				}
			}
		}
	}
	if (nVoxInside == 0)
		LOG_ERROR(
		        "We have a level set body with nVoxInside = 0, this is not expected neither supported."); // may happen when positive values restrict to point, line, or surface just along some lsGrid boundaries. Note that a body with only one gridpoint with a positive level set value, strictly within lsGrid, leads to a non zero volume

	volume = nVoxInside * dV;
	xMean /= nVoxInside;
	yMean /= nVoxInside;
	zMean /= nVoxInside;
	center = Vector3r(xMean, yMean, zMean);
	if (center.norm()
	    > pow(3., 0.5) * spac) // possible offsets almost reach one voxel edge in 1D (e.g. an external but tangent voxel is considered to be inside)
		LOG_ERROR(
		        "Incorrect LevelSet description: shape center is equal to " << center << " in local axes, instead of 0 (modulo a " << spac
		                                                                    << " grid spacing).");
	// we had to first compute center before being able to compute inertia now (in a 2nd loop):
	Real Ixx(0), Ixy(0), Ixz(0), Iyy(0), Iyz(0), Izz(0);
	Real xV, yV, zV;
	for (int xIndex = 0; xIndex < nGPx - 1; xIndex++) { // we will stop before the last grid points to avoid counting twice the last (or first) voxels
		for (int yIndex = 0; yIndex < nGPy - 1; yIndex++) {
			for (int zIndex = 0; zIndex < nGPz - 1; zIndex++) {
				if (distField[xIndex][yIndex][zIndex] <= 0) {
					gp = lsGrid->gridPoint(xIndex, yIndex, zIndex);
					xV = gp[0] + spac / 2.; // again, we choose taking the middle of this inside voxel, instead of the gridpoint
					yV = gp[1] + spac / 2.;
					zV = gp[2] + spac / 2.;
					Ixx += (pow(yV - yMean, 2) + pow(zV - zMean, 2)) * dV;
					Iyy += (pow(xV - xMean, 2) + pow(zV - zMean, 2)) * dV;
					Izz += (pow(xV - xMean, 2) + pow(yV - yMean, 2)) * dV;
					Ixy -= (xV - xMean) * (yV - yMean) * dV;
					Ixz -= (xV - xMean) * (zV - zMean) * dV;
					Iyz -= (yV - yMean) * (zV - zMean) * dV;
				}
			}
		}
	}
	if (math::abs(Ixy) + math::abs(Ixz) + math::abs(Iyz) > 9.e-3 * (math::abs(Ixx) + math::abs(Iyy) + math::abs(Izz))) // TODO: sthg else than 9.e-3 ?
		LOG_ERROR(
		        "Incorrect LevelSet description: local frame is non_inertial with a out-of-diagonal / diagonal terms ratio = "
		        << (math::abs(Ixy) + math::abs(Ixz) + math::abs(Iyz)) / (math::abs(Ixx) + math::abs(Iyy) + math::abs(Izz)));
	inertia = Vector3r(Ixx, Iyy, Izz);
	if (!surfNodes
	             .size()) { // 0 boundary nodes as of now. Condition could be false after save/load, where we do not need to duplicate the nodes. Other choice to avoid this test would be to save everything else that's computed by init(), eg center, ..
		lengthChar = twoD ? sqrt(volume / lsGrid->spacing) : cbrt(volume); // for the node location tolerance
		initSurfNodes();
		sphericity = maxRad / minRad;
	} // note that sphericity is always saved (but maxRad, minRad not)
	initDone = true;
}

Real LevelSet::distance(const Vector3r& pt) const
{
	// We work here in the "reference configuration"
	Vector3i indices = lsGrid->closestCorner(pt);
	int      xInd(indices[0]), yInd(indices[1]), zInd(indices[2]);
	if (xInd < 0 || yInd < 0 || zInd < 0) { // operators precedence OK in || vs <
		LOG_ERROR("Can not compute the distance, returning NaN");
		return (NaN);
	}
	Real                               f0yz(NaN), f1yz(NaN); // distance values at the same y and z than pt and for a x-value just before (resp. after) pt
	std::array<Real, 2>                yzCoord = { pt[1], pt[2] };
	std::array<Real, 2>                yExtr   = { lsGrid->gridPoint(xInd, yInd, zInd)[1], lsGrid->gridPoint(xInd, yInd + 1, zInd)[1] };
	std::array<Real, 2>                zExtr   = { lsGrid->gridPoint(xInd, yInd, zInd)[2], lsGrid->gridPoint(xInd, yInd, zInd + 1)[2] };
	std::array<std::array<Real, 2>, 2> knownValx0;
	knownValx0[0][0] = distField[xInd][yInd][zInd];
	knownValx0[0][1] = distField[xInd][yInd][zInd + 1];
	knownValx0[1][0] = distField[xInd][yInd + 1][zInd];
	knownValx0[1][1] = distField[xInd][yInd + 1][zInd + 1];
	std::array<std::array<Real, 2>, 2> knownValx1;
	knownValx1[0][0] = distField[xInd + 1][yInd][zInd];
	knownValx1[0][1] = distField[xInd + 1][yInd][zInd + 1];
	knownValx1[1][0] = distField[xInd + 1][yInd + 1][zInd];
	knownValx1[1][1] = distField[xInd + 1][yInd + 1][zInd + 1];
	f0yz             = ShopLS::biInterpolate(yzCoord, yExtr, zExtr, knownValx0);
	f1yz             = ShopLS::biInterpolate(yzCoord, yExtr, zExtr, knownValx1);
	return (pt[0] - lsGrid->gridPoint(xInd, yInd, zInd)[0]) / lsGrid->spacing * (f1yz - f0yz) + f0yz;
}

Real LevelSet::getVolume()
{
	if (!initDone) init();
	return volume;
}

Vector3r LevelSet::getInertia()
{
	if (!initDone) init();
	return inertia;
}

void LevelSet::computeMarchingCubes()
{
	// Clear possible existing marching cube data
	// (might matter later if the level set grid changes size)
	marchingCubesData.triangles.clear();
	marchingCubesData.normals.clear();

	// Compute the Marching Cubes triangulation
	MarchingCube mc;
	Vector3i     gpPerAxis = lsGrid->nGP;
	int          nGPx(gpPerAxis[0]), nGPy(gpPerAxis[1]), nGPz(gpPerAxis[2]);
	mc.init(nGPx, nGPy, nGPz, lsGrid->min, lsGrid->max());
	mc.computeTriangulation(distField, 0);

	// Removing garbage values
	vector<Vector3r> obtTriangles(mc.getTriangles());
	obtTriangles.resize(3 * mc.getNbTriangles());
	vector<Vector3r> obtNormals(mc.getNormals());
	obtNormals.resize(3 * mc.getNbTriangles());

	// Store newly computed marching cube data
	marchingCubesData.triangles   = obtTriangles;
	marchingCubesData.normals     = obtNormals;
	marchingCubesData.nbTriangles = mc.getNbTriangles();

	initDoneMarchingCubes = true;
}

vector<Vector3r> LevelSet::getMarchingCubeTriangles()
{
	if (!initDone) init();
	if (!initDoneMarchingCubes) computeMarchingCubes();
	return marchingCubesData.triangles;
}

vector<Vector3r> LevelSet::getMarchingCubeNormals()
{
	if (!initDone) init();
	if (!initDoneMarchingCubes) computeMarchingCubes();
	return marchingCubesData.normals;
}

int LevelSet::getMarchingCubeNbTriangles()
{
	if (!initDone) init();
	if (!initDoneMarchingCubes) computeMarchingCubes();
	return marchingCubesData.nbTriangles;
}

} // namespace yade
#endif //YADE_LS_DEM
