# Files in include/subdiv2d: API/public headers
set(API
	AssertAndError.h
	SubdivContainer.h
	Subdivision2D.h
	Types.h)

# Files in src: implementation and private headers
set(SOURCES
	AssertAndError.cpp
	SubdivContainer.cpp
	Subdivision2D.cpp)

# Produce full paths to API headers.
set(FULLAPI)
foreach(APISRC ${API})
	set(APIFULLPATH "${CMAKE_CURRENT_SOURCE_DIR}/../include/subdiv2d/${APISRC}")
	list(APPEND FULLAPI "${APIFULLPATH}")
endforeach()

source_group(API FILES ${FULLAPI})

add_library(Subdivision2D STATIC
	${SOURCES}
	${FULLAPI})

target_include_directories(Subdivision2D
	PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include"
	INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>")
