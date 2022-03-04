/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "PolyMesh.hpp"

#include <lib/base/Math.hpp>

#include <boost/python/tuple.hpp>

#include <stdarg.h>
#include <stdio.h>

namespace yade {
namespace ymport {
	namespace foamfile {

		PolyMesh::PolyMesh(const std::string& baseDirPath, bool patchAsWall, bool emptyAsWall)
		        : m_baseDirPath(baseDirPath)
		{
			if (boost::filesystem::is_directory(m_baseDirPath) == false) { error("Not a directory."); }

			// NOTE: Support for .gz files?
			m_pointsParser = std::make_unique<PointsParser>(m_baseDirPath + "/points");

			m_pointsParser->parse();

			m_points = m_pointsParser->points();

			if (m_points.size() == 0) { error("No points were found."); }

			m_facesParser = std::make_unique<FacesParser>(m_baseDirPath + "/faces", m_points.size());

			m_facesParser->parse();

			m_faces = m_facesParser->faces();

			if (m_faces.size() == 0) { error("No faces were found."); }

			m_boundaryParser = std::make_unique<BoundaryParser>(m_baseDirPath + "/boundary");

			m_boundaryParser->parse();

			m_boundaries = m_boundaryParser->boundaries();

			if (m_boundaries.size() == 0) { error("No boundaries were found."); }

			buildFacets(patchAsWall, emptyAsWall);
		}

		void PolyMesh::buildFacets(bool patchAsWall, bool emptyAsWall)
		{
			for (auto& b : m_boundaries) {
				if (b.type == BoundaryParser::Boundary::Type::PATCH && patchAsWall == false) { continue; }

				if (b.type == BoundaryParser::Boundary::Type::EMPTY && emptyAsWall == false) { continue; }

				std::vector<FacesParser::Face> faces(b.nFaces);
				for (int i = 0; i < b.nFaces; ++i) {
					faces.emplace_back(m_faces[b.startFace + i]);
				}

				for (FacesParser::Face f : faces) {
					int v0i = f.v0;
					int v1i = f.v1;
					int v2i = f.v2;
					int v3i = f.v3;

					Vector3r f0p0(m_points[v0i].x, m_points[v0i].y, m_points[v0i].z);
					Vector3r f0p1(m_points[v1i].x, m_points[v1i].y, m_points[v1i].z);
					Vector3r f0p2(m_points[v2i].x, m_points[v2i].y, m_points[v2i].z);

					m_facetCoords.append(py::make_tuple(f0p0, f0p1, f0p2));

					Vector3r f1p0(m_points[v2i].x, m_points[v2i].y, m_points[v2i].z);
					Vector3r f1p1(m_points[v3i].x, m_points[v3i].y, m_points[v3i].z);
					Vector3r f1p2(m_points[v0i].x, m_points[v0i].y, m_points[v0i].z);

					m_facetCoords.append(py::make_tuple(f1p0, f1p1, f1p2));
				}
			}
		}

		void PolyMesh::error(const char* fmt, ...)
		{
			va_list args;

			va_start(args, fmt);
			size_t size = vsnprintf(NULL, 0, fmt, args);
			va_end(args);

			std::vector<char> data(size + 1);

			va_start(args, fmt);
			vsnprintf(data.data(), data.size(), fmt, args);
			va_end(args);

			std::string message = m_baseDirPath + ": ";
			message.append(std::string(data.begin(), data.end()));

			throw std::runtime_error(message);
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
