#ifdef YADE_VTK

#include "VTKRecorder.hpp"
// https://codeyarns.com/2014/03/11/how-to-selectively-ignore-a-gcc-warning/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcomment"
#pragma GCC diagnostic ignored "-Wsuggest-override"
#include <lib/compatibility/VTKCompatibility.hpp> // fix InsertNextTupleValue → InsertNextTuple name change (and others in the future)

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#include <vtkXMLPMultiBlockDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#pragma GCC diagnostic pop
#endif

// Code that generates this warning, Note: we cannot do this trick in yade. If we have a warning in yade, we have to fix it! See also https://gitlab.com/yade-dev/trunk/merge_requests/73
// This method will work once g++ bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431#c34 is fixed.
#include <vtkHexahedron.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#ifdef YADE_VTK_MULTIBLOCK
#include <vtkMultiBlockDataSet.h>
#include <vtkXMLMultiBlockDataWriter.h>
#endif
#include <vtkTriangle.h>
#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/ConcretePM.hpp>
#include <pkg/dem/JointedCohesiveFrictionalPM.hpp>
#include <pkg/dem/Lubrication.hpp>
#include <pkg/dem/Shop.hpp>
#include <pkg/dem/WirePM.hpp>
#include <pkg/pfv/PartialSatClayEngine.hpp>
#include <pkg/pfv/Thermal.hpp>
#ifdef YADE_LIQMIGRATION
#include <pkg/dem/ViscoelasticCapillarPM.hpp>
#endif
#include <pkg/dem/HertzMindlin.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

#ifdef YADE_MASK_ARBITRARY
#define GET_MASK(b) b->groupMask.to_ulong()
#else
#define GET_MASK(b) b->groupMask
#endif

void VTKRecorder::action_02()
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
	vector<Shop::bodyState> bodyStates;
	if (recActive[REC_STRESS]) Shop::getStressForEachBody(bodyStates);

	vector<Matrix3r> bStresses;
	if (recActive[REC_BSTRESS]) { Shop::getStressLWForEachBody(bStresses); }

	vector<Matrix3r> NCStresses, SCStresses, NLStresses, SLStresses, NPStresses;
	if (recActive[REC_LUBRICATION]) {
		Law2_ScGeom_ImplicitLubricationPhys::getStressForEachBody(NCStresses, SCStresses, NLStresses, SLStresses, NPStresses);
	}

#ifdef YADE_MPI
	const auto& subD = YADE_PTR_CAST<Subdomain>(scene->subD);
	const auto& sz   = parallelMode ? subD->ids.size() : scene->bodies->size();
	for (unsigned bId = 0; bId != sz; ++bId) {
		const auto& b = parallelMode ? (*scene->bodies)[subD->ids[bId]] : (*scene->bodies)[bId];
#else
	for (const auto& b : *scene->bodies) {
#endif
		if (!b) continue;
		if (mask != 0 && !b->maskCompatible(mask)) continue;
		if (recActive[REC_SPHERES]) {
			const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
			if (sphere) {
				if (skipNondynamic && b->state->blockedDOFs == State::DOF_ALL) continue;
				vtkIdType pid[1];
				Vector3r  pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				pid[0] = spheresPos->InsertNextPoint(pos);
				spheresCells->InsertNextCell(1, pid);
				radii->InsertNextValue(sphere->radius);
				if (recActive[REC_BSTRESS]) {
					const Matrix3r&                         bStress = bStresses[b->getId()];
					Eigen::SelfAdjointEigenSolver<Matrix3r> solver(
					        bStress); // bStress is probably not symmetric (= self-adjoint for real matrices), but the solver still works, considering only one half of bStress. Which is good since existence of (real) eigenvalues is not sure for not symmetric bStress..
					Matrix3r dirAll = solver.eigenvectors();
					Vector3r eigenVal
					        = solver.eigenvalues(); // cf http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html#a30caf3c3884a7f4a46b8ec94efd23c5e to be sure that eigenVal[i] * dirAll.col(i) = bStress * dirAll.col(i) and that eigenVal[0] <= eigenVal[1] <= eigenVal[2]
					spheresSigI->InsertNextValue(eigenVal[2]);
					spheresSigII->InsertNextValue(eigenVal[1]);
					spheresSigIII->InsertNextValue(eigenVal[0]);
					Vector3r dirI((Real)dirAll(0, 2), (Real)dirAll(1, 2), (Real)dirAll(2, 2));
					Vector3r dirII((Real)dirAll(0, 1), (Real)dirAll(1, 1), (Real)dirAll(2, 1));
					Vector3r dirIII((Real)dirAll(0, 0), (Real)dirAll(1, 0), (Real)dirAll(2, 0));
					spheresDirI->InsertNextTuple(dirI);
					spheresDirII->InsertNextTuple(dirII);
					spheresDirIII->InsertNextTuple(dirIII);
				}
				if (recActive[REC_ID]) spheresId->InsertNextValue(b->getId());
				if (recActive[REC_MASK]) spheresMask->InsertNextValue(GET_MASK(b));
				if (recActive[REC_MASS]) spheresMass->InsertNextValue(b->state->mass);
#ifdef YADE_MPI
				if (recActive[REC_SUBDOMAIN]) spheresSubdomain->InsertNextValue(b->subdomain);
#endif

#ifdef THERMAL
				if (recActive[REC_TEMP]) {
					auto* thState = b->state.get();
					spheresTemp->InsertNextValue(thState->temp);
				}
#endif
#ifdef PARTIALSAT
				if (recActive[REC_PARTIALSAT]) {
					PartialSatState* state = dynamic_cast<PartialSatState*>(b->state.get());
					spheresRadiiChange->InsertNextValue(state->radiiChange);
					spheresSuction->InsertNextValue(state->suction);
					spheresIncidentCells->InsertNextValue(state->lastIncidentCells);
				}
#endif
				if (recActive[REC_CLUMPID]) clumpId->InsertNextValue(b->clumpId);
				if (recActive[REC_COLORS]) {
					const Vector3r& color = sphere->color;
					spheresColors->InsertNextTuple(color);
				}
				if (recActive[REC_VELOCITY]) {
					Vector3r vel = Vector3r::Zero();
					if (scene->isPeriodic) { // Take care of cell deformation
						vel = scene->cell->bodyFluctuationVel(b->state->pos, b->state->vel, scene->cell->prevVelGrad)
						        + scene->cell->prevVelGrad * scene->cell->wrapShearedPt(b->state->pos);
					} else {
						vel = b->state->vel;
					}
					spheresLinVelVec->InsertNextTuple(vel);
					spheresLinVelLen->InsertNextValue(vel.norm());
					const Vector3r& angVel = b->state->angVel;
					spheresAngVelVec->InsertNextTuple(angVel);
					spheresAngVelLen->InsertNextValue(angVel.norm());
				}
				if (recActive[REC_STRESS]) {
					const Vector3r& stress = bodyStates[b->getId()].normStress;
					const Vector3r& shear  = bodyStates[b->getId()].shearStress;
					spheresNormalStressVec->InsertNextTuple(stress);
					spheresShearStressVec->InsertNextTuple(shear);
					spheresNormalStressNorm->InsertNextValue(stress.norm());
				}
				if (recActive[REC_LUBRICATION]) {
					const Matrix3r& ncs = NCStresses[b->getId()];
					const Matrix3r& scs = SCStresses[b->getId()];
					const Matrix3r& nls = NLStresses[b->getId()];
					const Matrix3r& sls = SLStresses[b->getId()];
					const Matrix3r& nps = NPStresses[b->getId()];

					spheresLubricationNormalContactStress->InsertNextTuple(ncs);
					spheresLubricationShearContactStress->InsertNextTuple(scs);
					spheresLubricationNormalLubricationStress->InsertNextTuple(nls);
					spheresLubricationShearLubricationStress->InsertNextTuple(sls);
					spheresLubricationNormalPotentialStress->InsertNextTuple(nps);
				}
				if (recActive[REC_FORCE]) {
					scene->forces.sync();
					const Vector3r& f  = scene->forces.getForce(b->getId());
					const Vector3r& t  = scene->forces.getTorque(b->getId());
					Real            fn = f.norm();
					Real            tn = t.norm();
					spheresForceLen->InsertNextValue(fn);
					spheresTorqueLen->InsertNextValue(tn);
					spheresForceVec->InsertNextTuple(f);
					spheresTorqueVec->InsertNextTuple(t);
				}

				if (recActive[REC_CPM]) {
					cpmDamage->InsertNextValue(YADE_PTR_CAST<CpmState>(b->state)->normDmg);
					const Matrix3r& ss = YADE_PTR_CAST<CpmState>(b->state)->stress;
					//Real s[3]={ss[0],ss[1],ss[2]};
					cpmStress->InsertNextTuple(ss);
				}

				if (recActive[REC_JCFPM]) {
					nbCracks->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->nbBrokenBonds);
					jcfpmDamage->InsertNextValue(YADE_PTR_CAST<JCFpmState>(b->state)->damageIndex);
				}

				if (recActive[REC_COORDNUMBER]) { spheresCoordNumb->InsertNextValue(b->coordNumber()); }
#ifdef YADE_SPH
				if (recActive[REC_SPH]) {
					spheresRhoSPH->InsertNextValue(b->state->rho);
					spheresPressSPH->InsertNextValue(b->state->press);
					spheresCoordNumbSPH->InsertNextValue(b->coordNumber());
				}
#endif

#ifdef YADE_DEFORM
				if (recActive[REC_DEFORM]) {
					const Sphere* sphereDef = dynamic_cast<Sphere*>(b->shape.get());
					spheresRealRad->InsertNextValue(b->state->dR + sphereDef->radius);
				}
#endif

#ifdef YADE_LIQMIGRATION
				if (recActive[REC_LIQ]) {
					spheresLiqVol->InsertNextValue(b->state->Vf);
					const Real tmpVolIter = liqVolIterBody(b);
					spheresLiqVolIter->InsertNextValue(tmpVolIter);
					spheresLiqVolTotal->InsertNextValue(tmpVolIter + b->state->Vf);
				}
#endif
				if (recActive[REC_MATERIALID]) spheresMaterialId->InsertNextValue(b->material->id);
				continue;
			}
		} // end rec sphere.
		if (recActive[REC_FACETS]) {
			const Facet* facet = dynamic_cast<Facet*>(b->shape.get());
			if (facet) {
				Vector3r                     pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				const vector<Vector3r>&      localPos   = facet->vertices;
				Matrix3r                     facetAxisT = b->state->ori.toRotationMatrix();
				vtkSmartPointer<vtkTriangle> tri        = vtkSmartPointer<vtkTriangle>::New();
				vtkIdType                    nbPoints   = facetsPos->GetNumberOfPoints();
				for (int i = 0; i < 3; ++i) {
					Vector3r globalPos = pos + facetAxisT * localPos[i];
					facetsPos->InsertNextPoint(globalPos);
					tri->GetPointIds()->SetId(i, nbPoints + i);
				}
				facetsCells->InsertNextCell(tri);
				if (recActive[REC_COLORS]) {
					const Vector3r& color = facet->color;
					facetsColors->InsertNextTuple(color);
				}
				if (recActive[REC_STRESS]) {
					const Vector3r& stress = bodyStates[b->getId()].normStress + bodyStates[b->getId()].shearStress;
					facetsStressVec->InsertNextTuple(stress);
					facetsStressLen->InsertNextValue(stress.norm());
				}
				if (recActive[REC_FORCE]) {
					scene->forces.sync();
					const Vector3r& f  = scene->forces.getForce(b->getId());
					const Vector3r& t  = scene->forces.getTorque(b->getId());
					Real            fn = f.norm();
					Real            tn = t.norm();
					facetsForceLen->InsertNextValue(fn);
					facetsTorqueLen->InsertNextValue(tn);
					facetsForceVec->InsertNextTuple(f);
					facetsTorqueVec->InsertNextTuple(t);
				}

				if (recActive[REC_MATERIALID]) facetsMaterialId->InsertNextValue(b->material->id);
				if (recActive[REC_MASK]) facetsMask->InsertNextValue(GET_MASK(b));
				if (recActive[REC_COORDNUMBER]) { facetsCoordNumb->InsertNextValue(b->coordNumber()); }
				continue;
			}
		}
		if (recActive[REC_BOXES]) {
			const Box* box = dynamic_cast<Box*>(b->shape.get());
			if (box) {
				Vector3r                 pos(scene->isPeriodic ? scene->cell->wrapShearedPt(b->state->pos) : b->state->pos);
				Quaternionr              ori(b->state->ori);
				Vector3r                 ext(box->extents);
				vtkSmartPointer<vtkQuad> boxes = vtkSmartPointer<vtkQuad>::New();
				Vector3r                 A     = Vector3r(-ext[0], -ext[1], -ext[2]);
				Vector3r                 B     = Vector3r(-ext[0], +ext[1], -ext[2]);
				Vector3r                 C     = Vector3r(+ext[0], +ext[1], -ext[2]);
				Vector3r                 D     = Vector3r(+ext[0], -ext[1], -ext[2]);

				Vector3r E = Vector3r(-ext[0], -ext[1], +ext[2]);
				Vector3r F = Vector3r(-ext[0], +ext[1], +ext[2]);
				Vector3r G = Vector3r(+ext[0], +ext[1], +ext[2]);
				Vector3r H = Vector3r(+ext[0], -ext[1], +ext[2]);

				A = pos + ori * A;
				B = pos + ori * B;
				C = pos + ori * C;
				D = pos + ori * D;
				E = pos + ori * E;
				F = pos + ori * F;
				G = pos + ori * G;
				H = pos + ori * H;

				addWallVTK(boxes, boxesPos, A, B, C, D);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, E, H, G, F);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, A, E, F, B);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, G, H, D, C);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, F, G, C, B);
				boxesCells->InsertNextCell(boxes);

				addWallVTK(boxes, boxesPos, D, H, E, A);
				boxesCells->InsertNextCell(boxes);

				for (int i = 0; i < 6; i++) {
					if (recActive[REC_COLORS]) {
						const Vector3r& color = box->color;
						boxesColors->InsertNextTuple(color);
					}
					if (recActive[REC_STRESS]) {
						const Vector3r& stress = bodyStates[b->getId()].normStress + bodyStates[b->getId()].shearStress;
						boxesStressVec->InsertNextTuple(stress);
						boxesStressLen->InsertNextValue(stress.norm());
					}
					if (recActive[REC_FORCE]) {
						scene->forces.sync();
						const Vector3r& f  = scene->forces.getForce(b->getId());
						const Vector3r& t  = scene->forces.getTorque(b->getId());
						Real            fn = f.norm();
						Real            tn = t.norm();
						boxesForceVec->InsertNextTuple(f);
						boxesTorqueVec->InsertNextTuple(t);
						boxesForceLen->InsertNextValue(fn);
						boxesTorqueLen->InsertNextValue(tn);
					}
					if (recActive[REC_MATERIALID]) boxesMaterialId->InsertNextValue(b->material->id);
					if (recActive[REC_MASK]) boxesMask->InsertNextValue(GET_MASK(b));
				}
				continue;
			}
		}
	} // end bodies loop.

#ifdef YADE_MPI
	if ((!parallelMode and recActive[REC_PERICELL]) or (scene->subdomain == 0 and recActive[REC_PERICELL]))
#else
	if (recActive[REC_PERICELL])
#endif
	{
		const Matrix3r& hSize = scene->cell->hSize;
		Vector3r        v0    = hSize * Vector3r(0, 0, 1);
		Vector3r        v1    = hSize * Vector3r(0, 1, 1);
		Vector3r        v2    = hSize * Vector3r(1, 1, 1);
		Vector3r        v3    = hSize * Vector3r(1, 0, 1);
		Vector3r        v4    = hSize * Vector3r(0, 0, 0);
		Vector3r        v5    = hSize * Vector3r(0, 1, 0);
		Vector3r        v6    = hSize * Vector3r(1, 1, 0);
		Vector3r        v7    = hSize * Vector3r(1, 0, 0);
		pericellPoints->InsertNextPoint(v0);
		pericellPoints->InsertNextPoint(v1);
		pericellPoints->InsertNextPoint(v2);
		pericellPoints->InsertNextPoint(v3);
		pericellPoints->InsertNextPoint(v4);
		pericellPoints->InsertNextPoint(v5);
		pericellPoints->InsertNextPoint(v6);
		pericellPoints->InsertNextPoint(v7);
		vtkSmartPointer<vtkHexahedron> h = vtkSmartPointer<vtkHexahedron>::New();
		vtkIdList*                     l = h->GetPointIds();
		for (int i = 0; i < 8; i++) {
			l->SetId(i, i);
		}
		pericellHexa->InsertNextCell(h);
	}
}
#undef GET_MASK

} // namespace yade

#endif /* YADE_VTK */
