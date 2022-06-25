/*************************************************************************
*  2022 Janek Kozicki                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "SupportStatus.hpp"

namespace yade { // Cannot have #include directive inside.

YADE_ENUM(yade::Status, SupportStatus, (MATURE)(DEVELOPMENT)(ABANDONED)(UNKNOWN));
YADE_PLUGIN((Support));

}
