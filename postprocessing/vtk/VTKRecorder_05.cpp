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
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkZLibDataCompressor.h>

#ifdef YADE_MPI
#include <core/Subdomain.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#include <mpi.h>
#include <vtkXMLPUnstructuredGridWriter.h>
#pragma GCC diagnostic pop
#endif

// Code that generates this warning, Note: we cannot do this trick in yade. If we have a warning in yade, we have to fix it! See also https://gitlab.com/yade-dev/trunk/merge_requests/73
// This method will work once g++ bug https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431#c34 is fixed.
#pragma GCC diagnostic pop

#include <core/Scene.hpp>
#include <pkg/common/Box.hpp>
#include <pkg/common/Facet.hpp>
#include <pkg/common/Sphere.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/unordered_map.hpp>

namespace yade { // Cannot have #include directive inside.

void VTKRecorder::action_05()
{
	if (compress) compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();

	/* In mpiruns, spheres and other types of bodies are held by workers */
	spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifdef YADE_MPI
	if ((recActive[REC_SPHERES] and scene->subdomain != 0) or (!parallelMode and recActive[REC_SPHERES]))
#else
	if (recActive[REC_SPHERES])
#endif
	{
		spheresUg->SetPoints(spheresPos);
		spheresUg->SetCells(VTK_VERTEX, spheresCells);
		spheresUg->GetPointData()->AddArray(radii);
		if (recActive[REC_ID]) spheresUg->GetPointData()->AddArray(spheresId);
		if (recActive[REC_MASK]) spheresUg->GetPointData()->AddArray(spheresMask);
		if (recActive[REC_MASS]) spheresUg->GetPointData()->AddArray(spheresMass);
		if (recActive[REC_TEMP]) spheresUg->GetPointData()->AddArray(spheresTemp);
		if (recActive[REC_CLUMPID]) spheresUg->GetPointData()->AddArray(clumpId);
		if (recActive[REC_COLORS]) spheresUg->GetPointData()->AddArray(spheresColors);
		if (recActive[REC_VELOCITY]) {
			spheresUg->GetPointData()->AddArray(spheresLinVelVec);
			spheresUg->GetPointData()->AddArray(spheresAngVelVec);
			spheresUg->GetPointData()->AddArray(spheresLinVelLen);
			spheresUg->GetPointData()->AddArray(spheresAngVelLen);
		}
#ifdef YADE_MPI
		if (recActive[REC_SUBDOMAIN]) spheresUg->GetPointData()->AddArray(spheresSubdomain);
#endif
#ifdef YADE_SPH
		if (recActive[REC_SPH]) {
			spheresUg->GetPointData()->AddArray(spheresRhoSPH);
			spheresUg->GetPointData()->AddArray(spheresPressSPH);
			spheresUg->GetPointData()->AddArray(spheresCoordNumbSPH);
		}
#endif

#ifdef YADE_DEFORM
		if (recActive[REC_DEFORM]) { spheresUg->GetPointData()->AddArray(spheresRealRad); }
#endif

#ifdef YADE_LIQMIGRATION
		if (recActive[REC_LIQ]) {
			spheresUg->GetPointData()->AddArray(spheresLiqVol);
			spheresUg->GetPointData()->AddArray(spheresLiqVolIter);
			spheresUg->GetPointData()->AddArray(spheresLiqVolTotal);
		}
#endif

#ifdef PARTIALSAT
		if (recActive[REC_PARTIALSAT]) {
			spheresUg->GetPointData()->AddArray(spheresRadiiChange);
			spheresUg->GetPointData()->AddArray(spheresSuction);
			spheresUg->GetPointData()->AddArray(spheresIncidentCells);
		}
#endif
		if (recActive[REC_STRESS]) {
			spheresUg->GetPointData()->AddArray(spheresNormalStressVec);
			spheresUg->GetPointData()->AddArray(spheresShearStressVec);
			spheresUg->GetPointData()->AddArray(spheresNormalStressNorm);
		}
		if (recActive[REC_LUBRICATION]) {
			spheresUg->GetPointData()->AddArray(spheresLubricationNormalContactStress);
			spheresUg->GetPointData()->AddArray(spheresLubricationShearContactStress);
			spheresUg->GetPointData()->AddArray(spheresLubricationNormalLubricationStress);
			spheresUg->GetPointData()->AddArray(spheresLubricationShearLubricationStress);
			spheresUg->GetPointData()->AddArray(spheresLubricationNormalPotentialStress);
		}
		if (recActive[REC_FORCE]) {
			spheresUg->GetPointData()->AddArray(spheresForceVec);
			spheresUg->GetPointData()->AddArray(spheresForceLen);
			spheresUg->GetPointData()->AddArray(spheresTorqueVec);
			spheresUg->GetPointData()->AddArray(spheresTorqueLen);
		}
		if (recActive[REC_CPM]) {
			spheresUg->GetPointData()->AddArray(cpmDamage);
			spheresUg->GetPointData()->AddArray(cpmStress);
		}

		if (recActive[REC_JCFPM]) {
			spheresUg->GetPointData()->AddArray(nbCracks);
			spheresUg->GetPointData()->AddArray(jcfpmDamage);
		}
		if (recActive[REC_BSTRESS]) {
			spheresUg->GetPointData()->AddArray(spheresSigI);
			spheresUg->GetPointData()->AddArray(spheresSigII);
			spheresUg->GetPointData()->AddArray(spheresSigIII);
			spheresUg->GetPointData()->AddArray(spheresDirI);
			spheresUg->GetPointData()->AddArray(spheresDirII);
			spheresUg->GetPointData()->AddArray(spheresDirIII);
		}
		if (recActive[REC_MATERIALID]) spheresUg->GetPointData()->AddArray(spheresMaterialId);
		if (recActive[REC_COORDNUMBER]) spheresUg->GetCellData()->AddArray(spheresCoordNumb);

#ifdef YADE_VTK_MULTIBLOCK
		if (!multiblock)
#endif

#ifdef YADE_MPI
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			std::string fn;
			if (!parallelMode) {
				fn = fileName + "spheres." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			} else {
				fn = fileName + "spheres_" + boost::lexical_cast<string>(scene->iter) + "_" + boost::lexical_cast<string>(rank - 1) + ".vtu";
			}
			writer->SetFileName(fn.c_str());
			writer->SetInputData(spheresUg);
			writer->Write();

			if (rank == 1 && parallelMode) {
				vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
				string                                         pfn = fileName + "spheres_" + boost::lexical_cast<string>(scene->iter) + ".pvtu";
				pwriter->EncodeAppendedDataOff();
				pwriter->SetFileName(pfn.c_str());
				pwriter->SetNumberOfPieces(commSize - 1);
				pwriter->SetStartPiece(0);

				pwriter->SetInputData(spheresUg);
				pwriter->Update();
				pwriter->Write();
			}
		}
#else
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			string fn = fileName + "spheres." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			writer->SetFileName(fn.c_str());
			writer->SetInputData(spheresUg);
			writer->Write();
		}
#endif
	}


	/* facet bodies can be held by workers or master (like a wall body) */
	facetsUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if (recActive[REC_FACETS]) {
		facetsUg->SetPoints(facetsPos);
		facetsUg->SetCells(VTK_TRIANGLE, facetsCells);
		if (recActive[REC_COLORS]) facetsUg->GetCellData()->AddArray(facetsColors);
		if (recActive[REC_STRESS]) {
			facetsUg->GetCellData()->AddArray(facetsStressVec);
			facetsUg->GetCellData()->AddArray(facetsStressLen);
		}
		if (recActive[REC_FORCE]) {
			facetsUg->GetCellData()->AddArray(facetsForceVec);
			facetsUg->GetCellData()->AddArray(facetsForceLen);
			facetsUg->GetCellData()->AddArray(facetsTorqueVec);
			facetsUg->GetCellData()->AddArray(facetsTorqueLen);
		}
		if (recActive[REC_MATERIALID]) facetsUg->GetCellData()->AddArray(facetsMaterialId);
		if (recActive[REC_MASK]) facetsUg->GetCellData()->AddArray(facetsMask);
		if (recActive[REC_COORDNUMBER]) facetsUg->GetCellData()->AddArray(facetsCoordNumb);
#ifdef YADE_VTK_MULTIBLOCK
		if (!multiblock)
#endif
#ifdef YADE_MPI
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			std::string fn;
			if (!parallelMode) {
				fn = fileName + "facets." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			} else {
				fn = fileName + "facets_" + boost::lexical_cast<string>(scene->iter) + "_" + boost::lexical_cast<string>(rank) + ".vtu";
			}
			writer->SetFileName(fn.c_str());
			writer->SetInputData(facetsUg);
			writer->Write();

			if (rank == 0 && parallelMode) {
				vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
				string                                         pfn = fileName + "facets_" + boost::lexical_cast<string>(scene->iter) + ".pvtu";
				pwriter->EncodeAppendedDataOff();
				pwriter->SetFileName(pfn.c_str());
				pwriter->SetNumberOfPieces(commSize);
				pwriter->SetStartPiece(0);

				pwriter->SetInputData(facetsUg);
				pwriter->Update();
				pwriter->Write();
			}
		}
#else
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			string fn = fileName + "facets." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			writer->SetFileName(fn.c_str());
			writer->SetInputData(facetsUg);
			writer->Write();
		}
#endif
	}

	/* mpi case : assuming box bodies  are held by master  */
	vtkSmartPointer<vtkUnstructuredGrid> boxesUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifdef YADE_MPI
	if ((!parallelMode and recActive[REC_BOXES]) or (scene->subdomain == 0 and recActive[REC_BOXES]))
#else
	if (recActive[REC_BOXES])
#endif
	{
		boxesUg->SetPoints(boxesPos);
		boxesUg->SetCells(VTK_QUAD, boxesCells);
		if (recActive[REC_COLORS]) boxesUg->GetCellData()->AddArray(boxesColors);
		if (recActive[REC_STRESS]) {
			boxesUg->GetCellData()->AddArray(boxesStressVec);
			boxesUg->GetCellData()->AddArray(boxesStressLen);
		}
		if (recActive[REC_FORCE]) {
			boxesUg->GetCellData()->AddArray(boxesForceVec);
			boxesUg->GetCellData()->AddArray(boxesForceLen);
			boxesUg->GetCellData()->AddArray(boxesTorqueVec);
			boxesUg->GetCellData()->AddArray(boxesTorqueLen);
		}
		if (recActive[REC_MATERIALID]) boxesUg->GetCellData()->AddArray(boxesMaterialId);
		if (recActive[REC_MASK]) boxesUg->GetCellData()->AddArray(boxesMask);
#ifdef YADE_VTK_MULTIBLOCK
		if (!multiblock)
#endif
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			string fn = fileName + "boxes." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			writer->SetFileName(fn.c_str());
			writer->SetInputData(boxesUg);
			writer->Write();
		}
	}
}

} // namespace yade

#endif /* YADE_VTK */
