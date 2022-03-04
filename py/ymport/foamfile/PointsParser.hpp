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

		class PointsParser final : public Parser {
		public:
			struct Point {
				double x, y, z;
			};

			PointsParser(const boost::filesystem::path& path)
			        : Parser(path, "vectorField", "points")
			{
			}

			void parse() override;

			const std::vector<Point>& points() const { return m_points; }

		private:
			int                m_pointsToRead { 0 };
			std::vector<Point> m_points;
		};

	} // namespace foamfile
} // namespace ymport
} // namespace yade
