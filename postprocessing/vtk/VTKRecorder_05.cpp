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

void VTKRecorder::action_05()
{

	vtkSmartPointer<vtkDataCompressor> compressor;
	if (compress) compressor = vtkSmartPointer<vtkZLibDataCompressor>::New();


/* In mpiruns, spheres and other types of bodies are held by workers */
#ifdef YADE_MPI
	vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
	if ((recActive[REC_SPHERES] and scene->subdomain != 0) or (!parallelMode and recActive[REC_SPHERES]))
#else
	vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
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
	vtkSmartPointer<vtkUnstructuredGrid> facetsUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
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


	/* all subdomains have interactions in mpi case */
	vtkSmartPointer<vtkPolyData> intrPd = vtkSmartPointer<vtkPolyData>::New();
	if (recActive[REC_INTR]) {
		intrPd->SetPoints(intrBodyPos);
		intrPd->SetLines(intrCells);
		intrPd->GetCellData()->AddArray(intrForceN);
		intrPd->GetCellData()->AddArray(intrAbsForceT);
#ifdef YADE_LIQMIGRATION
		if (recActive[REC_LIQ]) {
			intrPd->GetCellData()->AddArray(liqVol);
			intrPd->GetCellData()->AddArray(liqVolNorm);
		}
#endif
		if (recActive[REC_JCFPM]) {
			intrPd->GetCellData()->AddArray(intrIsCohesive);
			intrPd->GetCellData()->AddArray(intrIsOnJoint);
			intrPd->GetCellData()->AddArray(eventNumber);
		}
		if (recActive[REC_WPM]) {
			intrPd->GetCellData()->AddArray(wpmNormalForce);
			intrPd->GetCellData()->AddArray(wpmLimitFactor);
		}
		if (recActive[REC_HERTZMINDLIN]) {
#ifdef PARTIALSAT
			intrPd->GetCellData()->AddArray(intrBrokenHertz);
			intrPd->GetCellData()->AddArray(intrDisp);
#endif
		}
#ifdef YADE_VTK_MULTIBLOCK
		if (!multiblock)
#endif

#ifdef YADE_MPI
		{
			vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			std::string fn;
			if (!parallelMode) {
				fn = fileName + "intrs." + boost::lexical_cast<string>(scene->iter) + ".vtp";
			} else {
				fn = fileName + "intrs_" + boost::lexical_cast<string>(scene->iter) + "_" + boost::lexical_cast<string>(rank) + ".vtp";
			}
			writer->SetFileName(fn.c_str());
			writer->SetInputData(intrPd);
			writer->Write();

			if (rank == 0) {
				vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
				string                                 pfn     = fileName + "intrs_" + boost::lexical_cast<string>(scene->iter) + ".pvtp";
				pwriter->EncodeAppendedDataOff();
				pwriter->SetFileName(pfn.c_str());
				pwriter->SetNumberOfPieces(commSize);
				pwriter->SetStartPiece(0);

				pwriter->SetInputData(intrPd);
				pwriter->Update();
				pwriter->Write();
			}
		}
#else
		{
			vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			string fn = fileName + "intrs." + boost::lexical_cast<string>(scene->iter) + ".vtp";
			writer->SetFileName(fn.c_str());
			writer->SetInputData(intrPd);
			writer->Write();
		}
#endif
	}


	vtkSmartPointer<vtkUnstructuredGrid> pericellUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
#ifdef YADE_MPI
	if ((!parallelMode and recActive[REC_PERICELL]) or (scene->subdomain == 0 and recActive[REC_PERICELL]))
#else
	if (recActive[REC_PERICELL])
#endif
	{
		pericellUg->SetPoints(pericellPoints);
		pericellUg->SetCells(12, pericellHexa);
#ifdef YADE_VTK_MULTIBLOCK
		if (!multiblock)
#endif
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			if (compress) writer->SetCompressor(compressor);
			if (ascii) writer->SetDataModeToAscii();
			string fn = fileName + "pericell." + boost::lexical_cast<string>(scene->iter) + ".vtu";
			writer->SetFileName(fn.c_str());
			writer->SetInputData(pericellUg);
			writer->Write();
		}
	}


	/* in mpi case, just use master to do this */
#ifdef YADE_MPI
	if ((!parallelMode && recActive[REC_CRACKS]) || (scene->subdomain == 0 && recActive[REC_CRACKS]))
#else
	if (recActive[REC_CRACKS])
#endif
	{
		string                               fileCracks = "cracks_" + Key + ".txt";
		std::ifstream                        file(fileCracks.c_str(), std::ios::in);
		vtkSmartPointer<vtkUnstructuredGrid> crackUg = vtkSmartPointer<vtkUnstructuredGrid>::New();

		if (file) {
			while (!file.eof()) {
				std::string line;
				Real        iter, time, p0, p1, p2, type, size, n0, n1, n2, nrg, onJnt;
				while (std::getline(file, line)) { /* writes into string "line", a line of file "file". To go along diff. lines*/
					file >> iter >> time >> p0 >> p1 >> p2 >> type >> size >> n0 >> n1 >> n2 >> nrg >> onJnt;
					vtkIdType pid[1];
					pid[0] = crackPos->InsertNextPoint(Vector3r(p0, p1, p2));
					crackCells->InsertNextCell(1, pid);
					crackIter->InsertNextValue(iter);
					crackTime->InsertNextValue(time);
					crackType->InsertNextValue(type);
					crackSize->InsertNextValue(size);
					Vector3r n(n0, n1, n2);
					crackNorm->InsertNextTuple(n);
					crackNrg->InsertNextValue(nrg);
					crackOnJnt->InsertNextValue(onJnt);
				}
			}
			file.close();
		}
		//
		crackUg->SetPoints(crackPos);
		crackUg->SetCells(VTK_VERTEX, crackCells);
		crackUg->GetPointData()->AddArray(crackIter);
		crackUg->GetPointData()->AddArray(crackTime);
		crackUg->GetPointData()->AddArray(crackType);
		crackUg->GetPointData()->AddArray(crackSize);
		crackUg->GetPointData()->AddArray(
		        crackNorm); //see https://www.mail-archive.com/paraview@paraview.org/msg08166.html to obtain Paraview 2D glyphs conforming to this normal
		crackUg->GetPointData()->AddArray(crackNrg);
		crackUg->GetPointData()->AddArray(crackOnJnt);

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		if (compress) writer->SetCompressor(compressor);
		if (ascii) writer->SetDataModeToAscii();
		string fn = fileName + "cracks." + boost::lexical_cast<string>(scene->iter) + ".vtu";
		writer->SetFileName(fn.c_str());
		writer->SetInputData(crackUg);
		writer->Write();
	}

// doing same thing for moments that we did for cracks:
#ifdef YADE_MPI
	if ((!parallelMode && recActive[REC_MOMENTS]) || (scene->subdomain == 0 && recActive[REC_MOMENTS]))
#else
	if (recActive[REC_MOMENTS])
#endif
	{
		string                               fileMoments = "moments_" + Key + ".txt";
		std::ifstream                        file(fileMoments.c_str(), std::ios::in);
		vtkSmartPointer<vtkUnstructuredGrid> momentUg = vtkSmartPointer<vtkUnstructuredGrid>::New();

		if (file) {
			while (!file.eof()) {
				std::string line;
				Real        i, p0, p1, p2, moment, numInts, eventNum, time;
				while (std::getline(file, line)) { /* writes into string "line", a line of file "file". To go along diff. lines*/
					file >> i >> p0 >> p1 >> p2 >> moment >> numInts >> eventNum >> time;
					vtkIdType pid[1];
					pid[0] = momentPos->InsertNextPoint(Vector3r(p0, p1, p2));
					momentCells->InsertNextCell(1, pid);
					momentSize->InsertNextValue(moment);
					momentiter->InsertNextValue(i);
					momenttime->InsertNextValue(time);
					momentNumInts->InsertNextValue(numInts);
					momentEventNum->InsertNextValue(eventNum);
				}
			}
			file.close();
		}
		//
		momentUg->SetPoints(momentPos);
		momentUg->SetCells(VTK_VERTEX, momentCells);
		momentUg->GetPointData()->AddArray(momentiter);
		momentUg->GetPointData()->AddArray(momenttime);
		momentUg->GetPointData()->AddArray(momentSize);
		momentUg->GetPointData()->AddArray(momentNumInts);
		momentUg->GetPointData()->AddArray(momentEventNum);

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		if (compress) writer->SetCompressor(compressor);
		if (ascii) writer->SetDataModeToAscii();
		string fn = fileName + "moments." + boost::lexical_cast<string>(scene->iter) + ".vtu";
		writer->SetFileName(fn.c_str());
		writer->SetInputData(momentUg);
		writer->Write();
	}

#ifdef YADE_VTK_MULTIBLOCK
#ifdef YADE_MPI
	if (multiblock) {
		if (parallelMode) { LOG_WARN("Multiblock feature in MPI case untested."); }
		vtkSmartPointer<vtkMultiBlockDataSet> multiblockDataset = vtkSmartPointer<vtkMultiBlockDataSet>::New();
		int                                   i                 = 0;
		if (recActive[REC_SPHERES]) multiblockDataset->SetBlock(i++, spheresUg);
		if (recActive[REC_FACETS]) multiblockDataset->SetBlock(i++, facetsUg);
		if (recActive[REC_INTR]) multiblockDataset->SetBlock(i++, intrPd);
		if (recActive[REC_PERICELL]) multiblockDataset->SetBlock(i++, pericellUg);
		vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
		if (ascii) writer->SetDataModeToAscii();
		std::string fn;
		if (!parallelMode) {
			fn = fileName + boost::lexical_cast<string>(scene->iter) + ".vtm";
		} else {
			fn = fileName + boost::lexical_cast<string>(scene->iter) + "_" + boost::lexical_cast<string>(rank) + ".vtm";
		}
		writer->SetFileName(fn.c_str());
		writer->SetInputData(multiblockDataset);
		writer->Write();

		if (parallelMode and rank == 0) {
			// use vtkController ; controller->GetProcessCount();
			// 					vtkSmartPointer<vtkXMLPMultiBlockDataWriter> parWriter = vtkXMLPMultiBlockDataWriter::New()
			// 					std::string pfn = fileName+boost::lexical_cast<string>(scene->iter+d)+".pvtm";
			// 					parWriter->EncodeAppendedDataOff();
			// 					parWriter->SetFileName(pfn.c_str());
			// 						parWriter->SetInputData(multiblockDataset);
			// 					parWriter->SetWriteMetaFile(commSize);
			// 					parWriter->Update();
			LOG_ERROR("PARALLEL MULTIBLOCKDATA NOT IMPLEMENTED. ");
		}
	}
#else
	if (multiblock) {
		vtkSmartPointer<vtkMultiBlockDataSet> multiblockDataset = vtkSmartPointer<vtkMultiBlockDataSet>::New();
		int                                   i                 = 0;
		if (recActive[REC_SPHERES]) multiblockDataset->SetBlock(i++, spheresUg);
		if (recActive[REC_FACETS]) multiblockDataset->SetBlock(i++, facetsUg);
		if (recActive[REC_INTR]) multiblockDataset->SetBlock(i++, intrPd);
		if (recActive[REC_PERICELL]) multiblockDataset->SetBlock(i++, pericellUg);
		vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
		if (ascii) writer->SetDataModeToAscii();
		string fn = fileName + boost::lexical_cast<string>(scene->iter) + ".vtm";
		writer->SetFileName(fn.c_str());
		writer->SetInputData(multiblockDataset);
		writer->Write();
	}
#endif

#endif
}

} // namespace yade

#endif /* YADE_VTK */
