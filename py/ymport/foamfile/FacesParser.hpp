/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include "Parser.hpp"

#include <vector>

namespace yade {
namespace ymport {
	namespace foamfile {

		class FacesParser final : public Parser {
		public:
			struct Face {
				std::size_t v0, v1, v2, v3;
			};

			FacesParser(const boost::filesystem::path& path, std::size_t numberOfPoints)
			        : Parser(path, "faceList", "faces")
			        , m_numberOfPoints(numberOfPoints)
			{
			}

			void parse() override;

			const std::vector<Face>& faces() const { return m_faces; }

		private:
			int               m_facesToRead { 0 };
			std::size_t       m_numberOfPoints;
			std::vector<Face> m_faces;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
