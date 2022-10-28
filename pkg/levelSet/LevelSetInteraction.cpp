/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <pkg/levelSet/LevelSetInteraction.hpp>
#include <pkg/levelSet/ShopLS.hpp>
#include <preprocessing/dem/Shop.hpp>

namespace yade {
YADE_PLUGIN((Bo1_LevelSet_Aabb)(Ig2_Box_LevelSet_ScGeom)(Ig2_LevelSet_LevelSet_ScGeom)(Ig2_Wall_LevelSet_ScGeom));
CREATE_LOGGER(Bo1_LevelSet_Aabb);
CREATE_LOGGER(Ig2_Box_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_LevelSet_LevelSet_ScGeom);
CREATE_LOGGER(Ig2_Wall_LevelSet_ScGeom);

void Bo1_LevelSet_Aabb::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*)
{ //TODO: use Eigen Aligned Box.extend to avoid the 2*6 if below ?
	// NB: see BoundDispatcher::processBody() (called by BoundDispatcher::action(), called by InsertionSortCollider::action()) in pkg/common/Dispatching.cpp for the attributes used upon calling
	if (!bv) { bv = shared_ptr<Bound>(new Aabb); }
	Aabb*     aabb    = static_cast<Aabb*>(bv.get()); // no need to bother deleting that raw pointer: e.g. https://stackoverflow.com/q/53908753/9864634
	LevelSet* lsShape = static_cast<LevelSet*>(cm.get());
	Real      inf     = std::numeric_limits<Real>::infinity();
	// We compute the bounds from LevelSet->corners serving as an Aabb in local frame, and considering the transformation from that local frame
	// NB: it is useless to try to make something much similar to Bo1_Box_Aabb::go(), in case se3.position of our level set body would not be in the middle of the lsShape->corners box, contrary to Box bodies and their extents and "halfSize" in Bo1_Box_Aabb::go()

	if (!lsShape->corners.size()) { // kind of "first iteration"
		Vector3i nGP(lsShape->lsGrid->nGP);
		Vector3r gp;
		int      nGPx(nGP[0]), nGPy(nGP[1]), nGPz(nGP[2]);
		Real     xGP, yGP, zGP;
		Real     xMin(inf), xMax(-inf), yMin(inf), yMax(-inf), zMin(inf), zMax(-inf); // extrema values for the "inside" gridpoints
		// identifying the extrema coordinates within the surface using a nGridPoints^3 loop. Would be much better to do 3*nGridpoints^2 loops, if possible, searching for find_if(firstIterator to lsValAlongZAxis,lastIterator to lsValAlongZAxis, within surface) which is easy for the zAxis = lsGrid->gridPoint[xInd][yInd], but how accessing an yAxis ~ lsGrid->gridPoint[xInd][:][zInd] which does not exist in C++ ?
		for (int xInd = 0; xInd < nGPx; xInd++) {
			for (int yInd = 0; yInd < nGPy; yInd++) {
				for (int zInd = 0; zInd < nGPz; zInd++) {
					gp  = lsShape->lsGrid->gridPoint(xInd, yInd, zInd);
					xGP = gp[0];
					yGP = gp[1];
					zGP = gp[2];
					if (lsShape->distField[xInd][yInd][zInd] <= 0) {
						if (xGP < xMin) xMin = xGP;
						if (xGP > xMax) xMax = xGP;
						if (yGP < yMin) yMin = yGP;
						if (yGP > yMax) yMax = yGP;
						if (zGP < zMin) zMin = zGP;
						if (zGP > zMax) zMax = zGP;
					}
				}
			}
		}
		if ((xMin == xMax) or (yMin == yMax) or (zMin == zMax))
			LOG_WARN("One flat LevelSet body, as detected by shape.corners computation, was that expected ? (is the grid too coarse ?)");
		// right now, our *Min, *Max define a downwards-rounded Aabb (smaller than true surface), let's make it upwards-rounded below:
		Real g(lsShape->lsGrid->spacing);
		for (int iInd = 0; iInd < 2; iInd++) {
			for (int jInd = 0; jInd < 2; jInd++) {
				for (int kInd = 0; kInd < 2; kInd++)
					lsShape->corners.push_back(
					        Vector3r(iInd == 0 ? xMin - g : xMax + g, jInd == 0 ? yMin - g : yMax + g, kInd == 0 ? zMin - g : zMax + g));
			}
		}
	}
	if (lsShape->corners.size() != 8) LOG_ERROR("We have a LevelSet-shaped body with some shape.corners computed but not 8 of them !");
	std::array<Vector3r, 8> cornersCurrent; // current positions of corners (in global frame)
	for (int corner = 0; corner < 8; corner++)
		cornersCurrent[corner] = ShopLS::rigidMapping(lsShape->corners[corner], Vector3r::Zero(), se3.position, se3.orientation);
	// NB: corners should average to the origin. This is kind of tested through LevelSet.center in levelSetBody() Py function

	Real xCorner, yCorner, zCorner;                                           // the ones of cornersCurrent
	Real xMin(inf), xMax(-inf), yMin(inf), yMax(-inf), zMin(inf), zMax(-inf); // extrema values for the cornersCurrent
	for (int corner = 0; corner < 8; corner++) {
		xCorner = cornersCurrent[corner][0];
		yCorner = cornersCurrent[corner][1];
		zCorner = cornersCurrent[corner][2];
		if (xCorner < xMin) xMin = xCorner;
		if (xCorner > xMax) xMax = xCorner;
		if (yCorner < yMin) yMin = yCorner;
		if (yCorner > yMax) yMax = yCorner;
		if (zCorner < zMin) zMin = zCorner;
		if (zCorner > zMax) zMax = zCorner;
	}
	aabb->min = Vector3r(xMin, yMin, zMin);
	aabb->max = Vector3r(xMax, yMax, zMax);
}

bool Ig2_Box_LevelSet_ScGeom::go(
        const shared_ptr<Shape>& shape1,
        const shared_ptr<Shape>& shape2,
        const State&             state1,
        const State&             state2,
        const Vector3r&          shift2,
        const bool&              force,
        const shared_ptr<Interaction>& c)
{
	// 1.1 Preliminary declarations
	shared_ptr<Box>      boxSh = YADE_PTR_CAST<Box>(shape1);
	shared_ptr<LevelSet> lsSh  = YADE_PTR_CAST<LevelSet>(shape2);
	Vector3r             lsCenter(state2.pos+shift2), boxCenter(state1.pos);
	Vector3r             axisContact(Vector3r::Zero()); // I will need to .dot this vector with Vector3r => Vector3r (and not Vector3i) as well

	// 1.2 Checking the "contact direction", adopting the two following hypothesis:
	// - H1: lsCenter is outside of the box
	// - H2: projecting lsCenter along the contact direction towards the box will make lsCenter fall onto the  box, and not alongside. This is not always fulfilled, but will always be during triaxial box setups, for instance (as long as the box is closed).
	for (int axis = 0; axis < 3; axis++) {
		if (math::abs(lsCenter[axis] - boxCenter[axis]) > boxSh->extents[axis]) // TODO: account for a change in orientation of the box ?
			axisContact[axis] = 1; // TODO: once I'm sure this will happen only once, I could stop the loop here...
	}
	if (axisContact.norm() != 1) {
		LOG_ERROR(
		        "Problem while determining contact direction for "
		        << c->id1 << "-" << c->id2 << " : we got " << axisContact
		        << ". (0 0 0) means the LevelSet'd body has its center inside the box, which is not supported. Indeed, center = " << lsCenter
		        << " while boxCenter = " << boxCenter << " and extents = " << boxSh->extents);
	}
	LOG_DEBUG("axisContact = " << axisContact);
	Real boxC(axisContact.dot(boxCenter)), boxE(boxSh->extents.dot(axisContact)), lsC(lsCenter.dot(axisContact));

	// 2.1. Preliminary declarations for the surface nodes loop
	// clang-format off
	const int nNodes(lsSh->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");
	vector<Vector3r> lsNodes; // current positions of the boundary nodes (some of of those, at least) for the level set body. See for loop below
	lsNodes.reserve(nNodes); // nNodes will be a maximum size, reserve() is appropriate, not resize() (see also https://github.com/isocpp/CppCoreGuidelines/issues/493)
	vector<Real> distList; // will include all distance values from the node(s)On1 to shape2. It is a std::vector because we do not know beforehand the number of elements in this "list
	distList.reserve(nNodes); // nNodes might be a maximum size, reserve() is appropriate, not resize() (see also https://github.com/isocpp/CppCoreGuidelines/issues/493)
	vector<int> indicesNodes; // distList will include distance values corresponding to nodes e.g. 2, 5, 7 only out of 10 nodes among surfNodes. This indicesNodes vector will store these 2,5,7 indices
	indicesNodes.reserve(nNodes);
	// NB: it might actually be somewhat faster to not use these vectors and just compare new node with previous node, as done in Ig2_LevelSet_LevelSet*
	// I do not think it is critical for present Ig2_Box_LevelSet_ScGeom
	Real distToNode; // one distance value, for one node
	Real xNode;      // current position (on the axisContact of interest, can be something else than x-axis..) of the level set boundary node
	// clang-format on

	// 2.2. Actual loop over surface nodes
	for (int node = 0; node < nNodes; node++) {
		lsNodes[node] = ShopLS::rigidMapping(lsSh->surfNodes[node], Vector3r::Zero(), lsCenter, state2.ori);
		xNode         = lsNodes[node].dot(axisContact);
		if (xNode < boxC - boxE || xNode > boxC + boxE) continue;
		distToNode = (lsC > boxC ? boxC + boxE - xNode : xNode - (boxC - boxE));
		if (distToNode < 0) LOG_ERROR("Unexpected case ! We obtained " << distToNode << " while waiting for a positive quantity");
		distList.push_back(distToNode);
		indicesNodes.push_back(node);
	}

	// 2.3. Finishing the work when there is no contact
	if (!distList.size()) { // all boundary nodes are outside the bounds' overlap,
		if (!c->isReal()) return false;
		else {
			c->geom = ShopLS::geomPtrForLaterRemoval(state1, state2, c);
			return true;
		}
	}
	Real maxOverlap;
	maxOverlap = *std::max_element(distList.begin(), distList.end());
	if (maxOverlap < 0 && !c->isReal() && !force) // inspired by Ig2_Sphere_Sphere_ScGeom:
		return false;

	// 2.4. Finishing the work when there is a contact
	int indexCont = std::min_element(distList.begin(), distList.end())
	        - distList.begin();                             //this works: it seems min_element is returning here a random access iterator
	Vector3r normal((lsC > boxC ? 1. : -1.) * axisContact); // normal from the box body to the level set body, ie from 1 to 2, as expected.
	Real     rad((lsNodes[indicesNodes[indexCont]] - lsCenter).norm());
	c->geom = ShopLS::geomPtr(
	        lsNodes[indicesNodes[indexCont]] + maxOverlap / 2. * normal, // middle of overlapping volumes, as usual
	        maxOverlap,                                                  // does not work for very big/huge overlap
	        rad,
	        rad,
	        state1,
	        state2,
	        c,
	        normal,
	        shift2);
	return true;
}

bool Ig2_Wall_LevelSet_ScGeom::go(
        const shared_ptr<Shape>& shape1,
        const shared_ptr<Shape>& shape2,
        const State&             state1,
        const State&             state2,
        const Vector3r&          shift2,
        const bool&              force,
        const shared_ptr<Interaction>& c)
{
	shared_ptr<Wall>     wallSh = YADE_PTR_CAST<Wall>(shape1);
	shared_ptr<LevelSet> lsSh   = YADE_PTR_CAST<LevelSet>(shape2);
	Real                 lsPos(state2.pos[wallSh->axis]+shift2[wallSh->axis]), wallPos(state1.pos[wallSh->axis]);

	const int nNodes(lsSh->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");

	Real distToNode, // one wall-level set distance value (< 0 when contact), for one node
	        prevDistToNode(std::numeric_limits<Real>::infinity()),
	        nodePos, // current position along the Wall->axis of one given boundary node
	        maxOverlap(-1);
	Vector3r currNode,   // current position of one given boundary node
	        contactNode; // the boundary node which is the most inside the wall
	for (int node = 0; node < nNodes; node++) {
		currNode = ShopLS::rigidMapping(lsSh->surfNodes[node], Vector3r::Zero(), state2.pos+shift2, state2.ori);
		nodePos  = currNode[wallSh->axis];
		if (wallPos >= lsPos && wallPos <= nodePos) // first possibility for the wall to intersect the LevelSet body
			distToNode = wallPos - nodePos;
		else if (wallPos >= nodePos && wallPos <= lsPos) // second possibility for intersection
			distToNode = nodePos - wallPos;
		else
			continue; // go directly to next node
		if (distToNode > 0) LOG_ERROR("Unexpected case ! We obtained " << distToNode << " while waiting for a negative quantity");
		if (distToNode < prevDistToNode) {
			maxOverlap     = -distToNode;
			contactNode    = currNode;
			prevDistToNode = distToNode;
		}
	}
	if (maxOverlap < 0 && !c->isReal() && !force)
		return false; // we won't create the interaction in this case (but it is not our job here to delete it in case it already exists)
	Vector3r wallNormal(Vector3r::Zero()), normal(Vector3r::Zero());
	if (wallSh->axis == 0) wallNormal = Vector3r::UnitX();
	else if (wallSh->axis == 1)
		wallNormal = Vector3r::UnitY();
	else if (wallSh->axis == 2)
		wallNormal = Vector3r::UnitZ();
	normal = (wallPos - lsPos > 0 ? -1 : 1) * wallNormal; // Points from wall to particle center
	Real rad( (contactNode - state2.pos - shift2).norm() ); // Distance from surface to center of level-set body
	
	c->geom = ShopLS::geomPtr(
	        contactNode + maxOverlap / 2. * normal, // middle of overlapping volumes, as usual
	        maxOverlap,                             // does not work for very big/huge overlap
	        rad, // considering the 2* feature of radius* (see comments in ShopLS::geomPtr), this is what makes most sense ?
	        rad, // we keep in particular radius1/2 slightly greater than (contactPoint-center).norm. And we just use the same radii for the two particles, as  in Ig2_Box_Sphere_ScGeom
	        state1,
	        state2,
	        c,
	        normal,
	        shift2);
	return true;
}

bool Ig2_LevelSet_LevelSet_ScGeom::go(
        const shared_ptr<Shape>& shape1,
        const shared_ptr<Shape>& shape2,
        const State&             state1,
        const State&             state2,
        const Vector3r&          shift2,
        const bool&              force,
        const shared_ptr<Interaction>& c)
{
	// 1. We first determine the Aabb zone where bodies' bounds overlap. TODO: possible use of Eigen AlignedBox ?
	std::array<Real, 6>     overlap; // the xA,xB, yA,yB, zA,zB defining the Aabb where bounds overlap: for x in [xA,xB] ; y in [yA,yB] ; ...
	const shared_ptr<Bound> bound1 = Body::byId(c->id1, scene)->bound;
	const shared_ptr<Bound> bound2 = Body::byId(c->id2, scene)->bound;
	for (int axis = 0; axis < 3; axis++) {
	  overlap[2 * axis]     = math::max(bound1->min[axis], bound2->min[axis]+shift2[axis]);
	  overlap[2 * axis + 1] = math::min(bound1->max[axis], bound2->max[axis]+shift2[axis]);
		if (overlap[2 * axis + 1]
		    < overlap[2 * axis]) { // an empty bound overlap here (= is possible and OK: it happens when bodies' bounds are just tangent)
			if (!c->isReal()) return false;
			else {
				c->geom = ShopLS::geomPtrForLaterRemoval(state1, state2, c);
				return true;
			}
		}
	}
	Vector3r minBoOverlap(Vector3r(overlap[0], overlap[2], overlap[4])),
	        maxBoOverlap(Vector3r(overlap[1], overlap[3], overlap[5])); // current configuration obviously
	// End of overlap computations

	// 2. We go now for the master-slave contact algorithm with surface ie boundary nodes:
	// 2.1 Preliminary declarations:
	shared_ptr<LevelSet> sh1 = YADE_PTR_CAST<LevelSet>(shape1);
	shared_ptr<LevelSet> sh2 = YADE_PTR_CAST<LevelSet>(shape2);
	bool                 id1isBigger(
                sh1->getVolume() > sh2->getVolume()); // identifying the smallest particle (where the master-slave contact algorithm will locate boundary nodes)
	shared_ptr<LevelSet> shB(id1isBigger ? sh1 : sh2); // the "big" shape
	shared_ptr<LevelSet> shS(id1isBigger ? sh2 : sh1); // the "small" shape
	// centr*ini, centr*end, and rot* below define the bodies' changes in configuration since the beginning of the simulation:
	Vector3r centrSini(Vector3r::Zero()),
	        centrBini(Vector3r::Zero()); // with a correct ls body creation and use, does not matter whether we take sh1->getCenter or 0
	Vector3r    centrSend(id1isBigger ? (state2.pos+shift2) : state1.pos), centrBend(id1isBigger ? state1.pos : (state2.pos+shift2));
	Quaternionr rotS(id1isBigger ? state2.ori : state1.ori),                    // ori = rotation from reference configuration (local axes) to current one
	        rotB(id1isBigger ? state1.ori.conjugate() : state2.ori.conjugate()) // ori.conjugate() from the current configuration to the reference one
	        ;
	const int nNodes(shS->surfNodes.size());
	if (!nNodes) LOG_ERROR("We have one level-set body without boundary nodes for contact detection. Will probably crash");
	Real distToNode, // one distance value, for one node
	        prevDistToNode(std::numeric_limits<Real>::infinity()), maxOverlap(-1);
	Vector3r minLSgrid(shB->lsGrid->min), maxLSgrid(shB->lsGrid->max()), nodeOfS // some node of smaller Body, in current configuration
	        ,
	        nodeOfSinB // mapped into initial configuration of larger Body
	        ,
	        normal, contactNode;

	// 2.2 Actual loop over surface nodes:
	for (int node = 0; node < nNodes; node++) {
		nodeOfS = ShopLS::rigidMapping(shS->surfNodes[node], centrSini, centrSend, rotS); // current position of this boundary node
		if (!Shop::isInBB(nodeOfS, minBoOverlap, maxBoOverlap)) continue;
		nodeOfSinB = ShopLS::rigidMapping(nodeOfS, centrBend, centrBini, rotB); // mapping this node into the big Body's local frame
		if (!Shop::isInBB(
		            nodeOfSinB,
		            minLSgrid,
		            maxLSgrid)) // possible when bodies (and their shape.corners) rotate, leading their bounds to possibly "inflate" (think of a sphere)
			continue;
		distToNode = shB->distance(nodeOfSinB);
		if (distToNode < 0 and distToNode < prevDistToNode) {
			maxOverlap = -distToNode;
			normal     = rotS
			        * shS->normal(
			                shS->surfNodes
			                        [node]); // shS->surfNodes[..] (and normal() to itself) refers to the initial shape, current normal is obtained with rotS *. It is for now the outward normal to the small particle
			if (id1isBigger)                 // if necessary, we make the normal from 1 to 2, as expected
				normal *= -1;
			contactNode    = nodeOfS;
			prevDistToNode = distToNode;
		}
	}

	// 2.3 Finishing the work:
	if (!c->isReal() && !force && maxOverlap < 0)
		return false; // we won't create the interaction in this case (but it is not our job here to delete it -- outside returning false unless force, see Ig2_Sphere_Sphere_ScGeom -- in case it already exists)
	c->geom = ShopLS::geomPtr(
	        contactNode - maxOverlap / 2. * normal, // middle of overlapping volumes, as usual
	        maxOverlap,                             // does not work for very big/huge overlap, eg when one sphere's center gets into another.
	        (contactNode - state1.pos - maxOverlap * (id1isBigger ? normal : Vector3r::Zero()) ).norm(), // radius1 value: taking center to surface distance for id1, expressed from contactNode and possible offset (depending on which particle contactNode came from)
	        (contactNode - state2.pos - shift2 - maxOverlap * (id1isBigger ? Vector3r::Zero() : normal) ).norm(), // radius2: same thing for id2. Reminder: contactNode does belong to id2 if id1 is bigger
	        state1,
	        state2,
	        c,
	        normal,
	        shift2);
	return true;
}
} // namespace yade
#endif //YADE_LS_DEM
