# Files in include/subdiv2d: API/public headers
set(API
	AssertAndError.h
	FixedMaxSizeArray.h
	IdTypes.h
	SubdivContainer.h
	Subdivision2D.h
	Types.h
	TypeSafeIndex.h
	TypeSafeIndexOutput.h)

# Files in src: implementation and private headers
set(SOURCES
	AssertAndError.cpp
	SubdivContainer.cpp
	Subdivision2D.cpp
	TypeSafeIndexIterable.h)

# Produce full paths to API headers.
set(FULLAPI "${CONFIG_HEADER}")
foreach(APISRC ${API})
	set(APIFULLPATH "${CMAKE_CURRENT_SOURCE_DIR}/../include/subdiv2d/${APISRC}")
	list(APPEND FULLAPI "${APIFULLPATH}")
endforeach()

source_group(API FILES ${FULLAPI})

add_library(Subdivision2D STATIC
	${SOURCES}
	${FULLAPI})

target_include_directories(Subdivision2D
	PRIVATE
	"${CMAKE_CURRENT_SOURCE_DIR}/../include"
	"${CONFIG_HEADER_DIR}"
	INTERFACE
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>"
	"$<BUILD_INTERFACE:${CONFIG_HEADER_DIR}>")
if(SUBDIV2D_USE_BOOST)
	target_include_directories(Subdivision2D
		PRIVATE "${Boost_INCLUDE_DIR}"

		# Only hard-code this path for boost in a build tree.
		INTERFACE "$<BUILD_INTERFACE:${Boost_INCLUDE_DIR}>")
endif()
