/*************************************************************************
*  2022 Tóth János                                                       *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//	Useful macros for development. Will not survive optimization.

#pragma once

#ifndef NDEBUG
#define TODO() _verify(false, __FILE__, __LINE__, "TODO")
#define UNREACHABLE() _verify(false, __FILE__, __LINE__, "UNREACHABLE")
#define VERIFY(x) _verify((x), __FILE__, __LINE__, "VERIFY FAILED: " #x)
#else
#define TODO()
#define UNREACHABLE()
#define VERIFY(x)
#endif

namespace yade {
namespace ymport {
	namespace foamfile {

		void _verify(bool cond, const char* file, int line, const char* msg);

	} // namespace foamfile
} // namespace ymport
} // namespace yade
