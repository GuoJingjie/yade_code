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

YADE_PLUGIN((VTKRecorder));
CREATE_LOGGER(VTKRecorder);

void VTKRecorder::action()
{
	action_01();
	action_02();
	action_03();
	//action_04(); // remember to uncomment this if you create postprocessing/vtk/VTKRecorder_04.cpp
};

void VTKRecorder::action_01()
{
#ifdef YADE_MPI
	if (parallelMode and !sceneRefreshed) {
		// update scene pointer (to avoid issues after scene broadcast in mpi case init)
		scene = Omega::instance().getScene().get();
		MPI_Comm_size(scene->getComm(), &commSize);
		MPI_Comm_rank(scene->getComm(), &rank);
		sceneRefreshed = true;
	}
#endif
	recActive = vector<bool>(REC_SENTINEL, false);
	FOREACH(string & rec, recorders)
	{
		if (rec == "all") {
			recActive[REC_SPHERES]     = true;
			recActive[REC_VELOCITY]    = true;
			recActive[REC_FACETS]      = true;
			recActive[REC_BOXES]       = true;
			recActive[REC_COLORS]      = true;
			recActive[REC_MASS]        = true;
			recActive[REC_INTR]        = true;
			recActive[REC_ID]          = true;
			recActive[REC_MASK]        = true;
			recActive[REC_CLUMPID]     = true;
			recActive[REC_MATERIALID]  = true;
			recActive[REC_STRESS]      = true;
			recActive[REC_FORCE]       = true;
			recActive[REC_COORDNUMBER] = true;
			if (scene->isPeriodic) { recActive[REC_PERICELL] = true; }
			recActive[REC_BSTRESS] = true;
		} else if (rec == "spheres") {
			recActive[REC_SPHERES] = true;
		} else if (rec == "velocity")
			recActive[REC_VELOCITY] = true;
		else if (rec == "facets")
			recActive[REC_FACETS] = true;
		else if (rec == "boxes")
			recActive[REC_BOXES] = true;
		else if (rec == "mass")
			recActive[REC_MASS] = true;
		else if (rec == "thermal")
			recActive[REC_TEMP] = true;
		else if ((rec == "colors") || (rec == "color"))
			recActive[REC_COLORS] = true;
		else if (rec == "cpm")
			recActive[REC_CPM] = true;
		else if (rec == "wpm")
			recActive[REC_WPM] = true;
		else if (rec == "intr")
			recActive[REC_INTR] = true;
		else if ((rec == "ids") || (rec == "id"))
			recActive[REC_ID] = true;
		else if (rec == "mask")
			recActive[REC_MASK] = true;
		else if ((rec == "clumpids") || (rec == "clumpId"))
			recActive[REC_CLUMPID] = true;
		else if (rec == "materialId")
			recActive[REC_MATERIALID] = true;
		else if (rec == "stress")
			recActive[REC_STRESS] = true;
		else if (rec == "force")
			recActive[REC_FORCE] = true;
		else if (rec == "jcfpm")
			recActive[REC_JCFPM] = true;
		else if (rec == "cracks")
			recActive[REC_CRACKS] = true;
		else if (rec == "moments")
			recActive[REC_MOMENTS] = true;
		else if (rec == "pericell" && scene->isPeriodic)
			recActive[REC_PERICELL] = true;
		else if (rec == "liquidcontrol")
			recActive[REC_LIQ] = true;
		else if (rec == "bstresses")
			recActive[REC_BSTRESS] = true;
		else if (rec == "coordNumber")
			recActive[REC_COORDNUMBER] = true;

		else if (rec == "SPH")
			recActive[REC_SPH] = true;
		else if (rec == "deform")
			recActive[REC_DEFORM] = true;
		else if (rec == "lubrication")
			recActive[REC_LUBRICATION] = true;
		else if (rec == "hertz")
			recActive[REC_HERTZMINDLIN] = true;
		else if (rec == "partialSat")
			recActive[REC_PARTIALSAT] = true;
		else
			LOG_ERROR(
			        "Unknown recorder named `" << rec
			                                   << "' (supported are: all, spheres, velocity, facets, boxes, color, stress, cpm, wpm, intr, id, "
			                                      "clumpId, materialId, jcfpm, cracks, moments, pericell, liquidcontrol, bstresses). Ignored.");
	}
#ifdef YADE_MPI
	if (parallelMode) { recActive[REC_SUBDOMAIN] = true; }
#endif
	// cpm needs interactions
	if (recActive[REC_CPM]) recActive[REC_INTR] = true;

	// jcfpm needs interactions
	if (recActive[REC_JCFPM]) recActive[REC_INTR] = true;

	// wpm needs interactions
	if (recActive[REC_WPM]) recActive[REC_INTR] = true;

	// liquid control needs interactions
	if (recActive[REC_LIQ]) recActive[REC_INTR] = true;

	// spheres
	spheresPos   = vtkSmartPointer<vtkPointsReal>::New();
	spheresCells = vtkSmartPointer<vtkCellArray>::New();

	radii = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	radii->SetNumberOfComponents(1);
	radii->SetName("radii");

	spheresSigI = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresSigI->SetNumberOfComponents(1);
	spheresSigI->SetName("sigI");

	spheresSigII = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresSigII->SetNumberOfComponents(1);
	spheresSigII->SetName("sigII");

	spheresSigIII = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresSigIII->SetNumberOfComponents(1);
	spheresSigIII->SetName("sigIII");

	spheresDirI = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresDirI->SetNumberOfComponents(3);
	spheresDirI->SetName("dirI");

	spheresDirII = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresDirII->SetNumberOfComponents(3);
	spheresDirII->SetName("dirII");

	spheresDirIII = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresDirIII->SetNumberOfComponents(3);
	spheresDirIII->SetName("dirIII");

	spheresMass = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresMass->SetNumberOfComponents(1);
	spheresMass->SetName("mass");

	spheresId = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresId->SetNumberOfComponents(1);
	spheresId->SetName("id");

	spheresTemp = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresTemp->SetNumberOfComponents(1);
	spheresTemp->SetName("temp");

#ifdef YADE_SPH
	spheresRhoSPH = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresRhoSPH->SetNumberOfComponents(1);
	spheresRhoSPH->SetName("SPH_Rho");

	spheresPressSPH = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresPressSPH->SetNumberOfComponents(1);
	spheresPressSPH->SetName("SPH_Press");

	spheresCoordNumbSPH = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresCoordNumbSPH->SetNumberOfComponents(1);
	spheresCoordNumbSPH->SetName("SPH_Neigh");
#endif

#ifdef YADE_DEFORM
	spheresRealRad = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresRealRad->SetNumberOfComponents(1);
	spheresRealRad->SetName("RealRad");
#endif

#ifdef YADE_LIQMIGRATION
	spheresLiqVol = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLiqVol->SetNumberOfComponents(1);
	spheresLiqVol->SetName("Liq_Vol");

	spheresLiqVolIter = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLiqVolIter->SetNumberOfComponents(1);
	spheresLiqVolIter->SetName("Liq_VolIter");

	spheresLiqVolTotal = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLiqVolTotal->SetNumberOfComponents(1);
	spheresLiqVolTotal->SetName("Liq_VolTotal");
#endif

	spheresMask = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresMask->SetNumberOfComponents(1);
	spheresMask->SetName("mask");

	spheresCoordNumb = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresCoordNumb->SetNumberOfComponents(1);
	spheresCoordNumb->SetName("coordNumber");

	clumpId = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	clumpId->SetNumberOfComponents(1);
	clumpId->SetName("clumpId");

	spheresColors = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresColors->SetNumberOfComponents(3);
	spheresColors->SetName("color");

	spheresLinVelVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLinVelVec->SetNumberOfComponents(3);
	spheresLinVelVec->SetName("linVelVec"); //Linear velocity in Vector3 form

	spheresLinVelLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresLinVelLen->SetNumberOfComponents(1);
	spheresLinVelLen->SetName("linVelLen"); //Length (magnitude) of linear velocity

	spheresAngVelVec = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresAngVelVec->SetNumberOfComponents(3);
	spheresAngVelVec->SetName("angVelVec"); //Angular velocity in Vector3 form

	spheresAngVelLen = vtkSmartPointer<vtkDoubleArrayFromReal>::New();
	spheresAngVelLen->SetNumberOfComponents(1);
	spheresAngVelLen->SetName("angVelLen"); //Length (magnitude) of angular velocity

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

void VTKRecorder::addWallVTK(vtkSmartPointer<vtkQuad>& boxes, vtkSmartPointer<vtkPointsReal>& boxesPos2, Vector3r& W1, Vector3r& W2, Vector3r& W3, Vector3r& W4)
{
	//Function for exporting walls of boxes
	vtkIdType nbPoints = boxesPos2->GetNumberOfPoints();

	boxesPos2->InsertNextPoint(W1);
	boxes->GetPointIds()->SetId(0, nbPoints + 0);

	boxesPos2->InsertNextPoint(W2);
	boxes->GetPointIds()->SetId(1, nbPoints + 1);

	boxesPos2->InsertNextPoint(W3);
	boxes->GetPointIds()->SetId(2, nbPoints + 2);

	boxesPos2->InsertNextPoint(W4);
	boxes->GetPointIds()->SetId(3, nbPoints + 3);
};

} // namespace yade

#endif /* YADE_VTK */
