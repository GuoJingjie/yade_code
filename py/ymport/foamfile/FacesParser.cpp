/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "FacesParser.hpp"
#include "VerifyMacros.hpp"

namespace yade {
namespace ymport {
	namespace foamfile {

		void FacesParser::parse()
		{
			m_facesToRead = getInt();

			expect({ '(' });

			int n;
			int v0, v1, v2, v3;

			m_faces.reserve(m_facesToRead);

			for (int i = 0; i < m_facesToRead; i++) {
				n = getInt();

				if (n != 4) { error("Face must contain 4 points, got: %d", n); }

				expect({ '(' });

#define CHECK_INDEX(ind)                                                                                                                                       \
	if ((ind) < 0) {                                                                                                                                       \
		error("Face index cannot be less than 0.");                                                                                                    \
	} else if ((std::size_t)(ind) >= m_numberOfPoints) {                                                                                                   \
		error("Face index is out of bounds (%d >= %d).", (ind), m_numberOfPoints);                                                                     \
	}

				v0 = getInt();
				CHECK_INDEX(v0);

				v1 = getInt();
				CHECK_INDEX(v1);

				v2 = getInt();
				CHECK_INDEX(v2);

				v3 = getInt();
				CHECK_INDEX(v3);

#undef CHECK_INDEX

				m_faces.emplace_back(Face { (std::size_t)v0, (std::size_t)v1, (std::size_t)v2, (std::size_t)v3 });

				expect({ ')' });
			}

			expect({ ')' });
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
