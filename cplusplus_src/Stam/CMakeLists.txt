cmake_minimum_required(VERSION 2.8)
project( SubMatrix )

set( CMAKE_BUILD_TYPE Debug )
# set( CMAKE_BUILD_TYPE Release )

include_directories( "ext/libigl/include" )
add_definitions( -DIGL_HEADER_ONLY )

set(SUBMATRIX_SOURCES
	SubMatrix.cpp
	)

add_executable( SubMatrix ${SUBMATRIX_SOURCES} )

set_target_properties( SubMatrix PROPERTIES COMPILE_FLAGS "-Wall -Woverloaded-virtual -Winit-self -Wno-sign-compare" )
