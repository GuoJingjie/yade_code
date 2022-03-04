/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "BoundaryParser.hpp"

namespace yade {
namespace ymport {
	namespace foamfile {

		void BoundaryParser::parse()
		{
			m_boundariesToRead = getInt();

			expect({ '(' });

			m_boundaries.reserve(m_boundariesToRead);

			for (int i = 0; i < m_boundariesToRead; ++i) {
				readBoundary();
			}

			expect({ ')' });
		}

		void BoundaryParser::readBoundary()
		{
			Token       token;
			Boundary    boundary;
			std::string type;
			while (!token.asCharEquals('}')) {
				token = m_lexer->getNextToken();

				if (token.isWord() && boundary.name.empty()) {
					boundary.name = token.asWord();

					expect({ '{' });
				} else if (token.asWordEquals("type")) {
					type = getWord();
					expect({ ';' });

					if (type == "patch") {
						boundary.type = Boundary::Type::PATCH;
					} else if (type == "wall") {
						boundary.type = Boundary::Type::WALL;
					} else if (type == "empty") {
						boundary.type = Boundary::Type::EMPTY;
					} else {
						error("Unknown type: %s", type.c_str());
					}
				} else if (token.asWordEquals("nFaces")) {
					boundary.nFaces = getInt();
					expect({ ';' });
				} else if (token.asWordEquals("startFace")) {
					boundary.startFace = getInt();
					expect({ ';' });
				} else if (token.asWordEquals("inGroups")) {
					skipUntil({ ';' });
				} else if (token.asCharEquals('}')) {
					VERIFY(boundary.name != "");
					VERIFY(boundary.type != Boundary::Type::UNKNOWN);
					m_boundaries.emplace_back(boundary);
					boundary.name = "";
				} else {
					error("Unexpected token: %s", token.debugString().c_str());
				}
			}
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
