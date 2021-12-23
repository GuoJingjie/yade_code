# This is a supplementary file for libVersions module. It forces finding versions of those
# libraries that do not have a detectable version number inside them
# The only way to determine version is by reading the actual source files.
# I do this by checking their md5sum.

# This doesn't work inside docker. Because uname -m returns parent architectuure, not the docker one.
# if(CMAKE_VERSION VERSION_LESS "3.17.0") # https://stackoverflow.com/questions/11944060/how-to-detect-target-architecture-using-cmake
#     EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
# else()                                  # https://cmake.org/cmake/help/latest/variable/CMAKE_HOST_SYSTEM_PROCESSOR.html  - new cmake version calls uname -m natively
#     set(ARCHITECTURE CMAKE_HOST_SYSTEM_PROCESSOR)
# endif()

# So let's use dpkg. A very simple and efficient method....
EXECUTE_PROCESS( COMMAND /usr/bin/dpkg --print-architecture COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
if (NOT ARCHITECTURE) # .... which doesn't work on other linux distributions. So someone might want to fix this later.
    # Warning: uname -m doesn't work properly inside i386 docker (but seems to work with ppc64le and aarch64)
    EXECUTE_PROCESS( COMMAND /usr/sbin/uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
    if (NOT ARCHITECTURE)
        set(ARCHITECTURE "unknown")
    endif()
endif()

message( STATUS "Architecture: ${ARCHITECTURE}" )