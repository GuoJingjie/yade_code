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
#include <vtkQuad.h>
#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
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
//	action_04();
	action_05();
//	action_06();
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
