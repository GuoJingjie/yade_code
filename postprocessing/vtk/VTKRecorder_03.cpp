#ifdef YADE_VTK

#include "VTKRecorder.hpp"
// https://codeyarns.com/2014/03/11/how-to-selectively-ignore-a-gcc-warning/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wsuggest-override"
#include <lib/compatibility/VTKCompatibility.hpp> // fix InsertNextTupleValue → InsertNextTuple name change (and others in the future)

#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#pragma GCC diagnostic pop
#endif

// Code that generates this warning, Note: we cannot do this trick in yade. If we have a warning in yade, we have to fix it! See also https://gitlab.com/yade-dev/trunk/merge_requests/73
// This method will work once g++ bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431#c34 is fixed.
#include <vtkLine.h>
#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/DemXDofGeom.hpp>
#include <pkg/dem/HertzMindlin.hpp>
#include <pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include <pkg/dem/Lubrication.hpp>
#include <pkg/dem/Shop.hpp>
#include <pkg/dem/WirePM.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

void VTKRecorder::action_03()
{
	if (recActive[REC_INTR]) {
		// holds information about cell distance between spatial and displayed position of each particle
		vector<Vector3i> wrapCellDist;
		if (scene->isPeriodic) { wrapCellDist.resize(scene->bodies->size()); }
		// save body positions, referenced by ids by vtkLine

		// map to keep real body ids and their number in a vector (intrBodyPos)
		boost::unordered_map<Body::id_t, Body::id_t> bIdVector;
		Body::id_t                                   curId = 0;
		for (const auto& b : *scene->bodies) {
			if (b) {
				if (!scene->isPeriodic) {
					intrBodyPos->InsertNextPoint(b->state->pos);
				} else {
					Vector3r pos = scene->cell->wrapShearedPt(b->state->pos, wrapCellDist[b->id]);
					intrBodyPos->InsertNextPoint(pos);
				}
				bIdVector.insert(std::pair<Body::id_t, Body::id_t>(b->id, curId));
				curId++;
			}
		}
		FOREACH(const shared_ptr<Interaction>& I, *scene->interactions)
		{
			if (!I->isReal()) continue;
			if (skipFacetIntr) {
				if (!(Body::byId(I->getId1()))) continue;
				if (!(Body::byId(I->getId2()))) continue;
				if (!(dynamic_cast<Sphere*>(Body::byId(I->getId1())->shape.get()))) continue;
				if (!(dynamic_cast<Sphere*>(Body::byId(I->getId2())->shape.get()))) continue;
			}

			const auto iterId1 = bIdVector.find(I->getId1());
			const auto iterId2 = bIdVector.find(I->getId2());

			if (iterId2 == bIdVector.end() || iterId2 == bIdVector.end()) continue;

			const auto setId1Line = iterId1->second;
			const auto setId2Line = iterId2->second;

			/* For the periodic boundary conditions,
				find out whether the interaction crosses the boundary of the periodic cell;
				if it does, display the interaction on both sides of the cell, with one of the
				points sticking out in each case.
				Since vtkLines must connect points with an ID assigned, we will create a new additional
				point for each point outside the cell. It might create some data redundancy, but
				let us suppose that the number of interactions crossing the cell boundary is low compared
				to total numer of interactions
			*/
			// how many times to add values defined on interactions, depending on how many times the interaction is saved
			int numAddValues = 1;
			// aperiodic boundary, or interaction is inside the cell
			if (!scene->isPeriodic || (scene->isPeriodic && (I->cellDist == wrapCellDist[I->getId2()] - wrapCellDist[I->getId1()]))) {
				vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0, setId1Line);
				line->GetPointIds()->SetId(1, setId2Line);
				intrCells->InsertNextCell(line);
			} else {
				assert(scene->isPeriodic);
				// spatial positions of particles
				const Vector3r& p01(Body::byId(I->getId1())->state->pos);
				const Vector3r& p02(Body::byId(I->getId2())->state->pos);
				// create two line objects; each of them has one endpoint inside the cell and the other one sticks outside
				// A,B are the "fake" bodies outside the cell for id1 and id2 respectively, p1,p2 are the displayed points
				// distance in cell units for shifting A away from p1; negated value is shift of B away from p2
				Vector3r        ptA(p01 + scene->cell->hSize * (wrapCellDist[I->getId2()] - I->cellDist).cast<Real>());
				const vtkIdType idPtA = intrBodyPos->InsertNextPoint(ptA);

				Vector3r        ptB(p02 + scene->cell->hSize * (wrapCellDist[I->getId1()] - I->cellDist).cast<Real>());
				const vtkIdType idPtB = intrBodyPos->InsertNextPoint(ptB);

				vtkSmartPointer<vtkLine> line1B(vtkSmartPointer<vtkLine>::New());
				line1B->GetPointIds()->SetId(0, setId2Line);
				line1B->GetPointIds()->SetId(1, idPtB);

				vtkSmartPointer<vtkLine> lineA2(vtkSmartPointer<vtkLine>::New());
				lineA2->GetPointIds()->SetId(0, idPtA);
				lineA2->GetPointIds()->SetId(1, setId2Line);
				numAddValues = 2;
			}
			const NormShearPhys*         phys = YADE_CAST<NormShearPhys*>(I->phys.get());
			const GenericSpheresContact* geom = YADE_CAST<GenericSpheresContact*>(I->geom.get());
			// gives _signed_ scalar of normal force, following the convention used in the respective constitutive law
			Real     fn = phys->normalForce.dot(geom->normal);
			Vector3r fs((Real)math::abs(phys->shearForce[0]), (Real)math::abs(phys->shearForce[1]), (Real)math::abs(phys->shearForce[2]));
			// add the value once for each interaction object that we created (might be 2 for the periodic boundary)
			for (int i = 0; i < numAddValues; i++) {
				intrAbsForceT->InsertNextTuple(fs);
				if (recActive[REC_WPM]) {
					const WirePhys* wirephys = dynamic_cast<WirePhys*>(I->phys.get());
					if (wirephys != NULL && wirephys->isLinked) {
						wpmLimitFactor->InsertNextValue(wirephys->limitFactor);
						wpmNormalForce->InsertNextValue(fn);
						intrForceN->InsertNextValue(NaN);
					} else {
						intrForceN->InsertNextValue(fn);
						wpmNormalForce->InsertNextValue(NaN);
						wpmLimitFactor->InsertNextValue(NaN);
					}
				} else if (recActive[REC_JCFPM]) {
					const JCFpmPhys* jcfpmphys = YADE_CAST<JCFpmPhys*>(I->phys.get());
					intrIsCohesive->InsertNextValue(jcfpmphys->isCohesive);
					intrIsOnJoint->InsertNextValue(jcfpmphys->isOnJoint);
					intrForceN->InsertNextValue(fn);
					eventNumber->InsertNextValue(jcfpmphys->eventNumber);
				} else if (recActive[REC_HERTZMINDLIN]) {
#ifdef PARTIALSAT
					const auto mindlinphys = YADE_CAST<MindlinPhys*>(I->phys.get());
					const auto mindlingeom = YADE_CAST<ScGeom*>(I->geom.get());
					intrBrokenHertz->InsertNextValue(mindlinphys->isBroken);
					intrDisp->InsertNextValue(mindlingeom->penetrationDepth - mindlinphys->initD);
#endif
					intrForceN->InsertNextValue(fn);
				} else {
					intrForceN->InsertNextValue(fn);
				}
#ifdef YADE_LIQMIGRATION
				if (recActive[REC_LIQ]) {
					const ViscElCapPhys* capphys = YADE_CAST<ViscElCapPhys*>(I->phys.get());
					liqVol->InsertNextValue(capphys->Vb);
					liqVolNorm->InsertNextValue(capphys->Vb / capphys->Vmax);
				}
#endif
			}
		}
	}
	//Additional Vector for storing forces
	bodyStates.clear();
	if (recActive[REC_STRESS]) Shop::getStressForEachBody(bodyStates);

	bStresses.clear();
	if (recActive[REC_BSTRESS]) { Shop::getStressLWForEachBody(bStresses); }

	NCStresses.clear();
	SCStresses.clear();
	NLStresses.clear();
	SLStresses.clear();
	NPStresses.clear();
	if (recActive[REC_LUBRICATION]) {
		Law2_ScGeom_ImplicitLubricationPhys::getStressForEachBody(NCStresses, SCStresses, NLStresses, SLStresses, NPStresses);
	}
}

} // namespace yade

#endif /* YADE_VTK */
