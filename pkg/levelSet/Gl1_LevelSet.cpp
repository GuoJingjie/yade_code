/*************************************************************************
*  2022 DLH van der Haven, dannyvdhaven@gmail.com                        *
*  This program is free software, see file LICENSE for details.          *
*  Code here is based on the Gl1_PotentialParticle.*pp by CWBoon 2015    *
*  for the implementation of potential particles in YADE.  				 *
*************************************************************************/

#ifdef YADE_LS_DEM
#include "Gl1_LevelSet.hpp"
#ifdef YADE_OPENGL
#include <lib/opengl/OpenGLWrapper.hpp>
#endif
#include <core/Aabb.hpp>
#include <pkg/dem/ScGeom.hpp>


namespace yade { // Cannot have #include directive inside.

#ifdef YADE_OPENGL

bool Gl1_LevelSet::recompute;
bool Gl1_LevelSet::wire;

void Gl1_LevelSet::go(const shared_ptr<Shape>& cm, const shared_ptr<State>& /*state*/, bool wire2, const GLViewInfo&)
{
	// Get the level set shape
	LevelSet* LS = static_cast<LevelSet*>(cm.get());

	if (recompute) { // Update the triangulation if we want that (necessary for fracturing or deforming particles)
		LS->computeMarchingCubes();
	}

	// Retrieve the triangulation (initialisation is done automatically within the level set shape)
	const vector<Vector3r>& triangles   = LS->getMarchingCubeTriangles();   //mc.getTriangles();
	const vector<Vector3r>& normals     = LS->getMarchingCubeNormals();     //mc.getNormals();
	int                     nbTriangles = LS->getMarchingCubeNbTriangles(); //mc.getNbTriangles();

	glColor3v(cm->color); //glColor3v is used when lighting is not enabled

	if (wire || wire2) {
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // Turn on wireframe mode. Render front and back faces of the wireframe
	} else {
		//TODO: Using glMaterialv(GL_FRONT,...) in conjunction with: glCullFace(GL_BACK); glEnable(GL_CULL_FACE); is the most cost-effective approach, since culling the back faces makes the visualisation lighter. An example why I don't activate this for now, is that in cubePPscaled.py we visualise the faces of the box as hollow, even with wire=False and culling the back faces makes the visualisation of the hollow particles confusing. Thus, for the time being I chose to keep and color the back faces; to be revisited @vsangelidakis
		//			glEnable(GL_NORMALIZE); //Not needed for vertex-based shading. The normals have been normalised inside the Marching Cubes script
		glMaterialv(
		        GL_FRONT_AND_BACK,
		        GL_AMBIENT_AND_DIFFUSE,
		        Vector3r(cm->color[0], cm->color[1], cm->color[2])); //glMaterialv is used when lighting is enabled
		glDisable(GL_CULL_FACE);
		//			glCullFace(GL_BACK); glEnable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);            // 2D
		glPolygonMode(GL_FRONT, GL_FILL); // Turn off wireframe mode
	}

	//			// VERTEX-BASED SHADING: Use the normal vector of each vertex of each triangle (makes the shading of each face look smoother)
	glBegin(GL_TRIANGLES)
		;
		for (int i = 0; i < 3 * nbTriangles; i += 3) {
			glNormal3v(normals[i + 2]);
			glVertex3v(triangles[i + 2]); //vertex #2 The sequence of the vertices specifies which side of the faces is front and which is back
			glNormal3v(normals[i + 1]);
			glVertex3v(triangles[i + 1]); //vertex #1
			glNormal3v(normals[i + 0]);
			glVertex3v(triangles[i + 0]); //vertex #0
		}
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	return;
}

YADE_PLUGIN((Gl1_LevelSet));
#endif // YADE_OPENGL

} // namespace yade
#endif //YADE_LS_DEM
