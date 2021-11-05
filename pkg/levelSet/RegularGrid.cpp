/*************************************************************************
*  2021 jerome.duriez@inrae.fr                                           *
*  This program is free software, see file LICENSE for details.          *
*************************************************************************/

#ifdef YADE_LS_DEM
#include <lib/high-precision/Constants.hpp>
#include <pkg/levelSet/RegularGrid.hpp>

namespace yade {

YADE_PLUGIN((RegularGrid));
CREATE_LOGGER(RegularGrid);

void RegularGrid::pyHandleCustomCtorArgs(boost::python::tuple& t, boost::python::dict&)
{                                               // to enable the same syntax as numpy.linspace
	if (boost::python::len(t) == 0) return; // nothing to do
	if (boost::python::len(t) != 3) throw invalid_argument("Expecting nothing else than 3 numbers ");
	Real minR = boost::python::extract<Real>(t[0])();
	Real maxR = boost::python::extract<Real>(t[1])();
	if (maxR <= minR) LOG_ERROR("You asked for a grid with min = " << minR << " greater than max = " << maxR << ": this is impossible");
	int nGPi = boost::python::extract<int>(t[2])();
	if (nGPi < 2) LOG_ERROR("You asked for a grid with " << nGPi << " gridpoints per axis: this is impossible (2 minimum)");
	this->min     = Vector3r(minR, minR, minR);
	this->nGP     = Vector3i(nGPi, nGPi, nGPi);
	this->spacing = (maxR - minR) / (nGPi - 1);
	t             = boost::python::
	        tuple(); // empty the args, necessary to avoid "RuntimeError: Zero (not 1) non-keyword constructor arguments" (default constructor seems to take on finishing the job)
}

Vector3i RegularGrid::closestCorner(const Vector3r& pt) const
{
	//	returns the indices of the closest gridpoint which is smaller (for all components) than pt, for locating pt within lsGrid
	bool     prob(false);
	Vector3r gridMax(max()), gridMin(min);
	Real     epsilonScaled(Mathr::EPSILON * spacing);
	//	Checking first whether pt lies within the x-, y- and z-extents of the RegularGrid. Using Mathr::EPSILON (scaled) is necessary to face true-story situations where a pt with -0.1 y coordinate is detected to be outside a grid starting at y=-0.1 (because of the loss in precision in the C++ <-> Python crossing ?)
	prob = prob || ((pt[0] > gridMax[0] + epsilonScaled) or (pt[0] < gridMin[0] - epsilonScaled)) // ||= does not exist because of short-circuit feature
	        || ((pt[1] > gridMax[1] + epsilonScaled) or (pt[1] < gridMin[1] - epsilonScaled))     // of || see https://stackoverflow.com/a/11156602. We
	        || ((pt[2] > gridMax[2] + epsilonScaled)
	            or (pt[2]
	                < gridMin[2] - epsilonScaled)); // could use |= (there is no difference between || and | here) but short-circuiting is nice in this case
	if (prob) {
		LOG_ERROR(
		        "You're asking for the closest grid point to a pt "
		        << pt << " which is outside the grid.. Returning negative indices, it may crash in few seconds.");
		return Vector3i(-1, -1, -1);
	}
	//	the to-be-returned indices:
	Vector3i retIndices;
	for (unsigned int index = 0; index < 3; index++) {
		retIndices[index] = int((pt[index] - gridMin[index]) / spacing);
		//		in case pt is along one grid boundary, we actually return 2nd to last index, not the last one:
		if (retIndices[index] == nGP[index] - 1) // > nGP[*] - 1 would have been detected by prob above
			retIndices[index] = nGP[index] - 2;
		if (retIndices[index] < 0) {
			if (math::abs(pt[index] - gridMin[index]) < epsilonScaled) // then the problem is just a matter of numeric precision
				retIndices[index] = 0;                             // and we are right to change xInd to 0 (we are at the minimum boundary)
			else
				LOG_ERROR(
				        "The closest corner to " << pt << " is said to be at index = " << retIndices[index] << " along axis " << index << endl);
		}
	}
	return retIndices;
}

Vector3r RegularGrid::gridPoint(int xIndex, int yIndex, int zIndex) const
{
	if ((xIndex > nGP[0] - 1) || (yIndex > nGP[1] - 1) || (zIndex > nGP[2] - 1)) {
		// an assert here would be strong... Let's leave users the possibility asking for wrong indexes without YADE crashing right away
		LOG_ERROR("You asked for a grid point (" << xIndex << "," << yIndex << "," << zIndex << ") which does not exist");
		return Vector3r(NaN, NaN, NaN);
	}
	return Vector3r(min[0] + xIndex * spacing, min[1] + yIndex * spacing, min[2] + zIndex * spacing);
}

Vector3r RegularGrid::max() const { return min + Vector3r(nGP[0] - 1, nGP[1] - 1, nGP[2] - 1) * spacing; }

Vector3r RegularGrid::getDims() const
{
	if (!nGP.squaredNorm()) {
		LOG_WARN("Not logical to ask for RegularGrid.dims before that any grid points have been defined");
		return Vector3r::Zero();
	}
	return Vector3r(nGP[0] - 1, nGP[1] - 1, nGP[2] - 1) * spacing;
}
} // namespace yade
#endif // YADE_LS_DEM
