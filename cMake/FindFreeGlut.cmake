# - Find FREEGLUT library
#
# This module defines
#  FREEGLUT_INCLUDE_DIR, where to find loki/Typelist.h, etc.
#  FREEGLUT_LIBRARY, libraries to link against to use GL2PS.
#  FREEGLUT_FOUND, If false, do not try to use GL2PS.

FIND_PATH(FREEGLUT_INCLUDE_DIR NAMES freeglut.h freeglut_std.h glut.h HINTS /usr/include/GL)
FIND_LIBRARY(FREEGLUT_LIBRARY NAMES glut)

# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FREEGLUT  DEFAULT_MSG  FREEGLUT_INCLUDE_DIR FREEGLUT_LIBRARY)

MARK_AS_ADVANCED(FREEGLUT_INCLUDE_DIR FREEGLUT_LIBRARY)
