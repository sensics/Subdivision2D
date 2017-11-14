list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(Eigen3 QUIET)

if(EIGEN3_FOUND)
	add_executable(AngleMesh
		AngleMesh.cpp
		AngleMeshUtils.cpp
		AngleMeshUtils.h
		EigenStdArrayInterop.h
		GenericExtremaFinder.h)
	target_link_libraries(AngleMesh PRIVATE Subdivision2D)
	target_include_directories(AngleMesh PRIVATE ${EIGEN3_INCLUDE_DIR})
	set_property(TARGET AngleMesh PROPERTY FOLDER Samples)
endif()
