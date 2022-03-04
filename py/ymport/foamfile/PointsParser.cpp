/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "PointsParser.hpp"
#include "VerifyMacros.hpp"

namespace yade {
namespace ymport {
	namespace foamfile {

		void PointsParser::parse()
		{
			m_pointsToRead = getInt();

			expect({ '(' });

			double x, y, z;

			m_points.reserve(m_pointsToRead);

			for (int i = 0; i < m_pointsToRead; i++) {
				expect({ '(' });
				x = getNumber();
				y = getNumber();
				z = getNumber();

				m_points.emplace_back(Point { x, y, z });

				expect({ ')' });
			}

			expect({ ')' });
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
