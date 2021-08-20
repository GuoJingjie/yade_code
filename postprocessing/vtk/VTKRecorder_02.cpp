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

#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <pkg/dem/ConcretePM.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

void VTKRecorder::action_02()
{
	spheresNormalStressVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresNormalStressVec->SetNumberOfComponents(3);
	spheresNormalStressVec->SetName("normalStress");

	spheresShearStressVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresShearStressVec->SetNumberOfComponents(3);
	spheresShearStressVec->SetName("shearStress");

	spheresNormalStressNorm = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresNormalStressNorm->SetNumberOfComponents(1);
	spheresNormalStressNorm->SetName("normalStressNorm");

	spheresLubricationNormalContactStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLubricationNormalContactStress->SetNumberOfComponents(9);
	spheresLubricationNormalContactStress->SetName("lubrication_NormalContactStress");

	spheresLubricationShearContactStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLubricationShearContactStress->SetNumberOfComponents(9);
	spheresLubricationShearContactStress->SetName("lubrication_ShearContactStress");

	spheresLubricationNormalLubricationStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLubricationNormalLubricationStress->SetNumberOfComponents(9);
	spheresLubricationNormalLubricationStress->SetName("lubrication_NormalLubricationStress");

	spheresLubricationShearLubricationStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLubricationShearLubricationStress->SetNumberOfComponents(9);
	spheresLubricationShearLubricationStress->SetName("lubrication_ShearLubricationStress");

	spheresLubricationNormalPotentialStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLubricationNormalPotentialStress->SetNumberOfComponents(9);
	spheresLubricationNormalPotentialStress->SetName("lubrication_NormalPotentialStress");

	spheresMaterialId = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresMaterialId->SetNumberOfComponents(1);
	spheresMaterialId->SetName("materialId");

	spheresForceVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresForceVec->SetNumberOfComponents(3);
	spheresForceVec->SetName("forceVec");

	spheresForceLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresForceLen->SetNumberOfComponents(1);
	spheresForceLen->SetName("forceLen");

	spheresTorqueVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresTorqueVec->SetNumberOfComponents(3);
	spheresTorqueVec->SetName("torqueVec");

	spheresTorqueLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresTorqueLen->SetNumberOfComponents(1);
	spheresTorqueLen->SetName("torqueLen");

	// facets
	facetsPos    = vtkSmartPointer<vtkPointsReal>::New();
	facetsCells  = vtkSmartPointer<vtkCellArray>::New();
	facetsColors = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsColors->SetNumberOfComponents(3);
	facetsColors->SetName("color");

	facetsStressVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsStressVec->SetNumberOfComponents(3);
	facetsStressVec->SetName("stressVec");

	facetsStressLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsStressLen->SetNumberOfComponents(1);
	facetsStressLen->SetName("stressLen");

	facetsMaterialId = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsMaterialId->SetNumberOfComponents(1);
	facetsMaterialId->SetName("materialId");

	facetsMask = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsMask->SetNumberOfComponents(1);
	facetsMask->SetName("mask");

	facetsForceVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsForceVec->SetNumberOfComponents(3);
	facetsForceVec->SetName("forceVec");

	facetsForceLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsForceLen->SetNumberOfComponents(1);
	facetsForceLen->SetName("forceLen");

	facetsTorqueVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsTorqueVec->SetNumberOfComponents(3);
	facetsTorqueVec->SetName("torqueVec");

	facetsTorqueLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsTorqueLen->SetNumberOfComponents(1);
	facetsTorqueLen->SetName("torqueLen");

	facetsCoordNumb = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	facetsCoordNumb->SetNumberOfComponents(1);
	facetsCoordNumb->SetName("coordNumber");

#ifdef YADE_MPI
	spheresSubdomain = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresSubdomain->SetNumberOfComponents(1);
	spheresSubdomain->SetName("subdomain");
#endif


	// boxes
	boxesPos    = vtkSmartPointer<vtkPointsReal>::New();
	boxesCells  = vtkSmartPointer<vtkCellArray>::New();
	boxesColors = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesColors->SetNumberOfComponents(3);
	boxesColors->SetName("color");

	boxesStressVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesStressVec->SetNumberOfComponents(3);
	boxesStressVec->SetName("stressVec");

	boxesStressLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesStressLen->SetNumberOfComponents(1);
	boxesStressLen->SetName("stressLen");

	boxesMaterialId = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesMaterialId->SetNumberOfComponents(1);
	boxesMaterialId->SetName("materialId");

	boxesMask = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesMask->SetNumberOfComponents(1);
	boxesMask->SetName("mask");

	boxesForceVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesForceVec->SetNumberOfComponents(3);
	boxesForceVec->SetName("forceVec");

	boxesForceLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesForceLen->SetNumberOfComponents(1);
	boxesForceLen->SetName("forceLen");

	boxesTorqueVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesTorqueVec->SetNumberOfComponents(3);
	boxesTorqueVec->SetName("torqueVec");

	boxesTorqueLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	boxesTorqueLen->SetNumberOfComponents(1);
	boxesTorqueLen->SetName("torqueLen");

	// interactions
	intrBodyPos = vtkSmartPointer<vtkPointsReal>::New();
	intrCells   = vtkSmartPointer<vtkCellArray>::New();
	intrForceN  = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrForceN->SetNumberOfComponents(1);
	intrForceN->SetName("forceN");
	intrAbsForceT = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrAbsForceT->SetNumberOfComponents(3);
	intrAbsForceT->SetName("absForceT");

	// pericell
	pericellPoints = vtkSmartPointer<vtkPointsReal>::New();
	pericellHexa   = vtkSmartPointer<vtkCellArray>::New();

	// extras for CPM
	if (recActive[REC_CPM]) {
		CpmStateUpdater csu;
		csu.update(scene);
	}
	cpmDamage = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	cpmDamage->SetNumberOfComponents(1);
	cpmDamage->SetName("cpmDamage");
	cpmStress = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	cpmStress->SetNumberOfComponents(9);
	cpmStress->SetName("cpmStress");

	// extras for JCFpm
	nbCracks = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	nbCracks->SetNumberOfComponents(1);
	nbCracks->SetName("nbCracks");
	jcfpmDamage = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	jcfpmDamage->SetNumberOfComponents(1);
	jcfpmDamage->SetName("damage");
	intrIsCohesive = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrIsCohesive->SetNumberOfComponents(1);
	intrIsCohesive->SetName("isCohesive");
	intrIsOnJoint = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrIsOnJoint->SetNumberOfComponents(1);
	intrIsOnJoint->SetName("isOnJoint");
	eventNumber = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	eventNumber->SetNumberOfComponents(1);
	eventNumber->SetName("eventNumber");

	// extras for cracks
	crackPos   = vtkSmartPointer<vtkPointsReal>::New();
	crackCells = vtkSmartPointer<vtkCellArray>::New();
	crackIter  = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackIter->SetNumberOfComponents(1);
	crackIter->SetName("iter");
	crackTime = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackTime->SetNumberOfComponents(1);
	crackTime->SetName("time");
	crackType = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackType->SetNumberOfComponents(1);
	crackType->SetName("type");
	crackSize = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackSize->SetNumberOfComponents(1);
	crackSize->SetName("size");
	crackNorm = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackNorm->SetNumberOfComponents(3);
	crackNorm->SetName("norm");
	crackNrg = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackNrg->SetNumberOfComponents(1);
	crackNrg->SetName("nrg");
	crackOnJnt = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	crackOnJnt->SetNumberOfComponents(1);
	crackOnJnt->SetName("onJnt");

	// extras for moments
	momentPos   = vtkSmartPointer<vtkPointsReal>::New();
	momentCells = vtkSmartPointer<vtkCellArray>::New();
	momentiter  = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	momentiter->SetNumberOfComponents(1);
	momentiter->SetName("momentiter");
	momenttime = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	momenttime->SetNumberOfComponents(1);
	momenttime->SetName("momenttime");
	momentSize = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	momentSize->SetNumberOfComponents(1);
	momentSize->SetName("momentSize");
	momentEventNum = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	momentEventNum->SetNumberOfComponents(1);
	momentEventNum->SetName("momentEventNum");
	momentNumInts = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	momentNumInts->SetNumberOfComponents(1);
	momentNumInts->SetName("momentNumInts");


#ifdef YADE_LIQMIGRATION
	liqVol = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	liqVol->SetNumberOfComponents(1);
	liqVol->SetName("liqVol");

	liqVolNorm = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	liqVolNorm->SetNumberOfComponents(1);
	liqVolNorm->SetName("liqVolNorm");
#endif

	// extras for WireMatPM
	wpmNormalForce = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	wpmNormalForce->SetNumberOfComponents(1);
	wpmNormalForce->SetName("wpmNormalForce");
	wpmLimitFactor = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	wpmLimitFactor->SetNumberOfComponents(1);
	wpmLimitFactor->SetName("wpmLimitFactor");

#ifdef PARTIALSAT
	// extras for hertzmindlin
	intrBrokenHertz = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrBrokenHertz->SetNumberOfComponents(1);
	intrBrokenHertz->SetName("broken");
	intrDisp = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	intrDisp->SetNumberOfComponents(1);
	intrDisp->SetName("disp");


	spheresRadiiChange = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresRadiiChange->SetNumberOfComponents(1);
	spheresRadiiChange->SetName("radiiChange");

	spheresSuction = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresSuction->SetNumberOfComponents(1);
	spheresSuction->SetName("suction");

	spheresIncidentCells = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresIncidentCells->SetNumberOfComponents(1);
	spheresIncidentCells->SetName("incidentCells");
#endif
}

} // namespace yade

#endif /* YADE_VTK */
