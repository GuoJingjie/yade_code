/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "BoundaryParser.hpp"
#include "FacesParser.hpp"
#include "PointsParser.hpp"

#include <lib/base/AliasNamespaces.hpp>

#include <boost/filesystem.hpp>
#include <boost/python/list.hpp>

#include <memory>
#include <string>

namespace yade {
namespace ymport {
	namespace foamfile {
		class PolyMesh {
		public:
			PolyMesh(const std::string& baseDirPath, bool patchAsWall, bool emptyAsWall);

			const py::list& facetCoords() const { return m_facetCoords; }

		private:
			void buildFacets(bool patchAsWall, bool emptyAsWall);

			void error(const char* fmt, ...);

			std::string m_baseDirPath;

			std::unique_ptr<PointsParser>   m_pointsParser { nullptr };
			std::unique_ptr<FacesParser>    m_facesParser { nullptr };
			std::unique_ptr<BoundaryParser> m_boundaryParser { nullptr };

			std::vector<PointsParser::Point>      m_points;
			std::vector<FacesParser::Face>        m_faces;
			std::vector<BoundaryParser::Boundary> m_boundaries;
			py::list                              m_facetCoords;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
