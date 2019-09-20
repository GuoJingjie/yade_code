#if defined (YADE_MPI) & defined (YADE_VTK)
#pragma once
#include<lib/compatibility/VTKCompatibility.hpp> // fix InsertNextTupleValue → InsertNextTuple name change (and others in the future)
#include<pkg/common/PeriodicEngines.hpp>
#include<vtkQuad.h>
#include<vtkSmartPointer.h>
#include <core/Subdomain.hpp> 
#include<vtkXMLPUnstructuredGridWriter.h>
#include<mpi.h>

// multiblock features don't seem to exist prioor to 5.2 
#if (VTK_MAJOR_VERSION==5 && VTK_MINOR_VERSION>=2) || (VTK_MAJOR_VERSION > 5)
	#define YADE_VTK_MULTIBLOCK
#endif


class VTKRecorderParallel: public PeriodicEngine {
	public:
  enum {REC_SPHERES=0,REC_FACETS,REC_BOXES,REC_COLORS,REC_MASS,REC_TEMP,REC_CPM,REC_INTR,REC_VELOCITY,REC_ID,REC_CLUMPID,REC_SENTINEL,REC_MATERIALID,REC_STRESS,REC_MASK,REC_RPM,REC_JCFPM,REC_CRACKS,REC_MOMENTS,REC_WPM,REC_PERICELL,REC_LIQ,REC_BSTRESS,REC_FORCE,REC_COORDNUMBER, REC_SUBDOMAIN, REC_SPH, REC_DEFORM};
		virtual void action();
		int procRank, commSize; 
		bool sceneRefreshed = false; 
		
		void addWallVTK (vtkSmartPointer<vtkQuad>& boxes, vtkSmartPointer<vtkPoints>& boxesPos, Vector3r& W1, Vector3r& W2, Vector3r& W3, Vector3r& W4);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(VTKRecorderParallel,PeriodicEngine,"Engine recording snapshots of simulation into series of \\*.vtu files, readable by VTK-based postprocessing programs such as Paraview. Both bodies (spheres and facets) and interactions can be recorded, with various vector/scalar quantities that are defined on them.\n\n:yref:`PeriodicEngine.initRun` is initialized to ``True`` automatically.",
		((bool,compress,false,,"Compress output XML files [experimental]."))
		((bool,ascii,false,,"Store data as readable text in the XML file (sets `vtkXMLWriter <http://www.vtk.org/doc/nightly/html/classvtkXMLWriter.html>`__ data mode to ``vtkXMLWriter::Ascii``, while the default is ``Appended``)"))
		((bool,skipFacetIntr,false,,"Skip interactions that are not of sphere-sphere type (e.g. sphere-facet, sphere-box...), when saving interactions")) //FIXME : turning this true leads to segfaults (seg fault 11, possibly accessing some null ptr
		((bool,skipNondynamic,false,,"Skip non-dynamic spheres (but not facets)."))
		#ifdef YADE_VTK_MULTIBLOCK
			((bool,multiblock,false,,"Use multi-block (``.vtm``) files to store data, rather than separate ``.vtu`` files."))
		#endif
		((string,fileName,"",,"Base file name; it will be appended with {spheres,intrs,facets}-243100.vtu (unless *multiblock* is ``True``) depending on active recorders and step number (243100 in this case). It can contain slashes, but the directory must exist already."))
		((vector<string>,recorders,vector<string>(1,string("all")),,"List of active recorders (as strings). ``all`` (the default value) enables all base and generic recorders.\n\n.. admonition:: Base recorders\n\n\tBase recorders save the geometry (unstructured grids) on which other data is defined. They are implicitly activated by many of the other recorders. Each of them creates a new file (or a block, if :yref:`multiblock <VTKRecorder.multiblock>` is set).\n\n\t``spheres``\n\t\tSaves positions and radii (``radii``) of :yref:`spherical<Sphere>` particles.\n\t``facets``\n\t\tSave :yref:`facets<Facet>` positions (vertices).\n\t``boxes``\n\t\tSave :yref:`boxes<Box>` positions (edges).\n\t``intr``\n\t\tStore interactions as lines between nodes at respective particles positions. Additionally stores magnitude of normal (``forceN``) and shear (``absForceT``) forces on interactions (the :yref:`geom<Interaction.geom> must be of type :yref:`NormShearPhys`). \n\n.. admonition:: Generic recorders\n\n\tGeneric recorders do not depend on specific model being used and save commonly useful data.\n\n\t``id``\n\t\tSaves id's (field ``id``) of spheres; active only if ``spheres`` is active.\n\t``mass``\n\t\tSaves masses (field ``mass``) of spheres; active only if ``spheres`` is active.\n\t``clumpId``\n\t\tSaves id's of clumps to which each sphere belongs (field ``clumpId``); active only if ``spheres`` is active.\n\t``colors``\n\t\tSaves colors of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>` (field ``color``); only active if ``spheres`` or ``facets`` are activated.\n\t``mask``\n\t\tSaves groupMasks of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>` (field ``mask``); only active if ``spheres`` or ``facets`` are activated.\n\t``materialId``\n\t\tSaves materialID of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>`; only active if ``spheres`` or ``facets`` are activated.\n\t``coordNumber``\n\t\tSaves coordination number (number of neighbours) of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>`; only active if ``spheres`` or ``facets`` are activated.\n\t``velocity``\n\t\tSaves linear and angular velocities of spherical particles as Vector3 and length(fields ``linVelVec``, ``linVelLen`` and ``angVelVec``, ``angVelLen`` respectively``); only effective with ``spheres``.\n\t``stress``\n\t\tSaves stresses of :yref:`spheres<Sphere>` and of :yref:`facets<Facet>`  as Vector3 and length; only active if ``spheres`` or ``facets`` are activated.\n\t``force``\n\t\tSaves force and torque of :yref:`spheres<Sphere>`, :yref:`facets<Facet>` and :yref:`boxes<Box>` as Vector3 and length (norm); only active if ``spheres``, ``facets`` or ``boxes`` are activated.\n\t``pericell``\n\t\tSaves the shape of the cell (simulation has to be periodic).\n\t``bstresses``\n\t\tConsidering the per-particle stress tensors as given by :yref:`bodyStressTensors<yade.utils.bodyStressTensors>`, saves the per-particle principal stresses, sigI (most tensile) $\\geq$ sigII $\\geq$ sigIII (most compressive), and the associated principal directions dirI, dirII, dirIII.\n\n.. admonition:: Specific recorders\n\n\tThe following should only be activated in appropriate cases, otherwise crashes can occur due to violation of type presuppositions.\n\n\t``cpm``\n\t\tSaves data pertaining to the :yref:`concrete model<Law2_ScGeom_CpmPhys_Cpm>`: ``cpmDamage`` (normalized residual strength averaged on particle), ``cpmStress`` (stress on particle); ``intr`` is activated automatically by ``cpm``\n\t``wpm``\n\t\tSaves data pertaining to the :yref:`wire particle model<Law2_ScGeom_WirePhys_WirePM>`: ``wpmForceNFactor`` shows the loading factor for the wire, e.g. normal force divided by threshold normal force.\n\t``jcfpm``\n\t\tSaves data pertaining to the :yref:`rock (smooth)-jointed model<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM>`: ``damage`` is defined by :yref:`JCFpmState.tensBreak` + :yref:`JCFpmState.shearBreak`; ``intr`` is activated automatically by ``jcfpm``, and :yref:`on joint<JCFpmPhys.isOnJoint>` or :yref:`cohesive<JCFpmPhys.isCohesive>` interactions can be vizualized.\n\
t``cracks``\n\t\tSaves other data pertaining to the :yref:`rock model<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM>`: ``cracks`` shows locations where cohesive bonds failed during the simulation, with their types (0/1  for tensile/shear breakages), their sizes (0.5*(R1+R2)), and their normal directions. The :yref:`corresponding attribute<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.recordCracks>` has to be activated, and Key attributes have to be consistent.\n\n"))
		((string,Key,"",,"Necessary if :yref:`recorders<VTKRecorder.recorders>` contains 'cracks'. A string specifying the name of file 'cracks___.txt' that is considered in this case (see :yref:`corresponding attribute<Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM.Key>`)."))
		((int,mask,0,,"If mask defined, only bodies with corresponding groupMask will be exported. If 0, all bodies will be exported.")),
		/*ctor*/
		initRun=true;
	);
	DECLARE_LOGGER;
};


REGISTER_SERIALIZABLE(VTKRecorderParallel);
#endif 