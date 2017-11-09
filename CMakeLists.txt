cmake_minimum_required(VERSION 3.1.0)
project(Subdiv2D)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-local")

###
# CMake Modules
###
#include(EnableExtraCompilerWarnings)

#include(SetDefaultBuildType)
#set_default_build_type(RelWithDebInfo)

#include(MSVCMultipleProcessCompile)

###
# Main library
###
add_subdirectory(src)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	# We're a root project - OK to test.
	include(CTest)
	add_subdirectory(vendor)
	add_subdirectory(tests)
endif()

