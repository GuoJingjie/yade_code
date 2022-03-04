#==================pygts=========================================
IF(ENABLE_GTS)
  find_python_module(gts)
  IF (PY_gts)
    MESSAGE(STATUS "Use system gts version")
  ELSE (PY_gts)
    MESSAGE(STATUS "Use embedded version of gts. Please, consider installing the corresponding package")
    FILE(GLOB SRC_PYGTS "${CMAKE_CURRENT_SOURCE_DIR}/3rd-party/pygts-0.3.1/*.cpp")
    ADD_LIBRARY(_gts SHARED ${SRC_PYGTS})
    SET_TARGET_PROPERTIES(_gts PROPERTIES PREFIX "")
    SET_TARGET_PROPERTIES(_gts PROPERTIES  UNITY_BUILD OFF) # a bunch of multiple declarations there
    TARGET_LINK_LIBRARIES(_gts gts)
    INSTALL(FILES 3rd-party/pygts-0.3.1/__init__.py DESTINATION ${YADE_LIB_PATH}/py/gts)
    INSTALL(FILES 3rd-party/pygts-0.3.1/pygts.py DESTINATION ${YADE_LIB_PATH}/py/gts)
    INSTALL(TARGETS _gts DESTINATION ${YADE_LIB_PATH}/py/gts)
  ENDIF (PY_gts)
ENDIF(ENABLE_GTS)
#==================pygts=========================================

#==================pyModules=========================================

FILE(GLOB filesPY "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
SET(filesPY "${filesPY};${CMAKE_CURRENT_SOURCE_DIR}/pack/pack.py")
INSTALL(FILES ${filesPY} DESTINATION ${YADE_PY_PATH}/yade)
INSTALL(FILES 3rd-party/mtTkinter-0.4/mtTkinter.py DESTINATION ${YADE_LIB_PATH}/py)
FILE(GLOB filesPYTests "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.py")
INSTALL(FILES ${filesPYTests} DESTINATION ${YADE_PY_PATH}/yade/tests)
FILE(GLOB filesPYPerf "${CMAKE_CURRENT_SOURCE_DIR}/../examples/test/performance/*")
INSTALL(FILES ${filesPYPerf} DESTINATION ${YADE_PY_PATH}/yade/tests/checks/performance)
FILE(GLOB filesFEMxDEM "${CMAKE_CURRENT_SOURCE_DIR}/FEMxDEM/*.py")
INSTALL(FILES ${filesFEMxDEM} DESTINATION ${YADE_PY_PATH}/yade/FEMxDEM)

ADD_LIBRARY(_libVersions SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_libVersions.cpp")
SET_TARGET_PROPERTIES(_libVersions PROPERTIES PREFIX "")
TARGET_LINK_LIBRARIES(_libVersions ${VTK_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(_libVersions PUBLIC ${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS})
INSTALL(TARGETS _libVersions DESTINATION "${YADE_PY_PATH}/yade/")

# hint: follow changes in d067b0696a8 to add new modules.
ADD_LIBRARY(_log SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_log.cpp")
SET_TARGET_PROPERTIES(_log PROPERTIES PREFIX "")
INSTALL(TARGETS _log DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_math SHARED "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_math.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeStorageOrdering.cpp")
SET_TARGET_PROPERTIES(_math PROPERTIES PREFIX "")
INSTALL(TARGETS _math DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_minieigenHP SHARED "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_minieigenHP.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeBoxes.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeComplex1.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeComplex2.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeConverters.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeMatrices1.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeMatrices2.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeQuaternion.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeStorageOrdering.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeVectors1.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_ExposeVectors2.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/high-precision/_RealHPDiagnostics.cpp")
SET_TARGET_PROPERTIES(_minieigenHP PROPERTIES PREFIX "")
INSTALL(TARGETS _minieigenHP DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(WeightedAverage2d SHARED "${CMAKE_CURRENT_SOURCE_DIR}/WeightedAverage2d.cpp")
SET_TARGET_PROPERTIES(WeightedAverage2d PROPERTIES PREFIX "")
INSTALL(TARGETS WeightedAverage2d DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_utils.cpp")
SET_TARGET_PROPERTIES(_utils PROPERTIES PREFIX "")
INSTALL(TARGETS _utils DESTINATION "${YADE_PY_PATH}/yade/")


ADD_LIBRARY(_polyhedra_utils SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_polyhedra_utils.cpp")
SET_TARGET_PROPERTIES(_polyhedra_utils PROPERTIES PREFIX "")
INSTALL(TARGETS _polyhedra_utils DESTINATION "${YADE_PY_PATH}/yade/")


ADD_LIBRARY(_packPredicates SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packPredicates.cpp")
SET_TARGET_PROPERTIES(_packPredicates PROPERTIES PREFIX "")
IF(ENABLE_GTS AND NOT(PY_gts))
  TARGET_LINK_LIBRARIES(_packPredicates _gts)
ENDIF(ENABLE_GTS AND NOT(PY_gts))


INSTALL(TARGETS _packPredicates DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_packSpheres SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packSpheres.cpp")
SET_TARGET_PROPERTIES(_packSpheres PROPERTIES PREFIX "")
INSTALL(TARGETS _packSpheres DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_packObb SHARED "${CMAKE_CURRENT_SOURCE_DIR}/pack/_packObb.cpp")
SET_TARGET_PROPERTIES(_packObb PROPERTIES PREFIX "")
INSTALL(TARGETS _packObb DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(wrapper SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/yadeWrapper.cpp")
SET_TARGET_PROPERTIES(wrapper PROPERTIES PREFIX "")
INSTALL(TARGETS wrapper DESTINATION "${YADE_PY_PATH}/yade/")

ADD_LIBRARY(_customConverters SHARED "${CMAKE_CURRENT_SOURCE_DIR}/wrapper/customConverters.cpp")
SET_TARGET_PROPERTIES(_customConverters PROPERTIES PREFIX "")
INSTALL(TARGETS _customConverters DESTINATION "${YADE_PY_PATH}/yade/")

FILE(GLOB_RECURSE YMPORT_CPP_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/ymport/*.cpp)
ADD_LIBRARY(_ymport SHARED "${CMAKE_CURRENT_SOURCE_DIR}/_ymport.cpp" ${YMPORT_CPP_SRCS})
SET_TARGET_PROPERTIES(_ymport PROPERTIES PREFIX "")
INSTALL(TARGETS _ymport DESTINATION "${YADE_PY_PATH}/yade/")

IF(ENABLE_LOGGER AND CMAKE_UNITY_BUILD)
	MESSAGE(STATUS "${Esc}[93mENABLE_LOGGER needs a workaround with CMAKE_UNITY_BUILD, adding .cpp files containing CREATE_CPP_LOCAL_LOGGER to SKIP_UNITY_BUILD_INCLUSION.${Esc}[0m")
	SET(FIND_CPP_LOGGER_ARGS -rZl CREATE_CPP_LOCAL_LOGGER ${CMAKE_SOURCE_DIR}/py --include=*.cpp)
	EXECUTE_PROCESS( COMMAND grep ${FIND_CPP_LOGGER_ARGS} COMMAND sed -e "s/\\x0/;/g" OUTPUT_VARIABLE FILES_WITH_CPP_LOGGER OUTPUT_STRIP_TRAILING_WHITESPACE)
	FOREACH(CPP_LOGGER_USED_HERE ${FILES_WITH_CPP_LOGGER})
		MESSAGE(STATUS "Adding ${CPP_LOGGER_USED_HERE} to SKIP_UNITY_BUILD_INCLUSION")
		SET_PROPERTY(SOURCE ${CPP_LOGGER_USED_HERE} PROPERTY SKIP_UNITY_BUILD_INCLUSION ON)
	ENDFOREACH(CPP_LOGGER_USED_HERE)
ENDIF(ENABLE_LOGGER AND CMAKE_UNITY_BUILD)

