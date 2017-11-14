cmake_minimum_required(VERSION 3.1.0)
project(Subdiv2D)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	set(IS_SUBPROJECT FALSE)
else()
	set(IS_SUBPROJECT TRUE)
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake-local")

###
# CMake Modules
###
#include(EnableExtraCompilerWarnings)

#include(SetDefaultBuildType)
#set_default_build_type(RelWithDebInfo)

#include(MSVCMultipleProcessCompile)

set(SUBDIV2D_USE_BOOST FALSE)
find_package(Boost QUIET)

if(Boost_FOUND)
	if(${Boost_VERSION} GREATER 105699)
		# static_vector introduced in boost 1.54, with initializer list support added in 1.57
		# this is new enough for stable vector with initializer lists.
		option(SUBDIV2D_USE_BOOST_STATIC_VECTOR "Allow use of Boost static_vector for simpler, faster passing of collections?" ON)

		if(NOT IS_SUBPROJECT)
			set(VERB "NOT use")
			if(SUBDIV2D_USE_BOOST_STATIC_VECTOR)
				set(VERB "use")
				set(SUBDIV2D_USE_BOOST TRUE)
			endif()

			message(STATUS "Subdivision2d: Configured to ${VERB} boost::container::static_vector")
		endif()
	endif()

	if(${Boost_VERSION} GREATER 105799)
		# small_vector introduced in boost 1.58
		option(SUBDIV2D_USE_BOOST_SMALL_VECTOR "Allow use of Boost small_vector to reduce runtime allocations?" ON)
		if(NOT IS_SUBPROJECT)
			set(VERB "NOT use")
			if(SUBDIV2D_USE_BOOST_SMALL_VECTOR)
				set(VERB "use")
				set(SUBDIV2D_USE_BOOST TRUE)
			endif()

			message(STATUS "Subdivision2d: Configured to ${VERB} boost::container::small_vector")
		endif()
	endif()
	if(IS_SUBPROJECT)
		mark_as_advanced(SUBDIV2D_USE_BOOST_STATIC_VECTOR SUBDIV2D_USE_BOOST_SMALL_VECTOR)
	endif()
endif()


###
# Main library
###
add_subdirectory(include/subdiv2d)
add_subdirectory(src)

if(NOT IS_SUBPROJECT)
	# We're a root project - OK to test.
	list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/catch2")
	include(CTest)
	add_subdirectory(vendor)
	add_subdirectory(tests)
	add_subdirectory(samples)
endif()

