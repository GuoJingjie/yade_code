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

		class BoundaryParser final : public Parser {
		public:
			struct Boundary {
				enum Type { UNKNOWN = 0, PATCH, WALL, EMPTY };

				std::string name;
				Type        type;
				int         nFaces;
				int         startFace;
			};

			BoundaryParser(const boost::filesystem::path& path)
			        : Parser(path, "polyBoundaryMesh", "boundary")
			{
			}

			void parse() override;

			const std::vector<Boundary>& boundaries() const { return m_boundaries; }

		private:
			void readBoundary();

			int                   m_boundariesToRead { 0 };
			std::vector<Boundary> m_boundaries;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
