/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "VerifyMacros.hpp"

#include <boost/filesystem.hpp>
#include <stdexcept>

namespace yade {
namespace ymport {
	namespace foamfile {

		void _verify(bool cond, const char* file, int line, const char* msg)
		{
			if (cond == false) {
				boost::filesystem::path p(file);

				std::string message = p.filename().string() + ":" + std::to_string(line) + ": " + std::string(msg);
				throw std::runtime_error(message);
			}
		}

	} // namespace foamfile
} // namespace ymport
} // namespace yade
