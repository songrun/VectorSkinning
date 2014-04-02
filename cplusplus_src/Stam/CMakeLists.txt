cmake_minimum_required(VERSION 2.8)
project( STAM )

set( CMAKE_BUILD_TYPE Debug )
# set( CMAKE_BUILD_TYPE Release )

INCLUDE_DIRECTORIES(/opt/local/include/eigen3)
LINK_DIRECTORIES(/opt/local/lib)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
## Use "export CXX=/opt/local/bin/g++-mp-4.8" on the command line first.
## or "export CXX=/usr/bin/clang++"
# SET(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-4.8")
# SET(CMAKE_CXX_COMPILER "/opt/local/bin/g++-mp-4.8")

# IF(APPLE)
#    INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
#    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
#    FIND_LIBRARY(GLUT_LIBRARY GLUT )
#    FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
#    MARK_AS_ADVANCED (COCOA_LIBRARY
#                      GLUT_LIBRARY
#                      OpenGL_LIBRARY)
#    SET(EXTRA_LIBS ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})
# ENDIF (APPLE)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

set(EXECUTABLE_OUTPUT_PATH 
	build
	)

include_directories( "ext/libigl/include" )
add_definitions( -DIGL_HEADER_ONLY )

set(STAM_SOURCES
	src/main.cpp src/catmull.c src/stam_eigen_reader.c src/catmull_obj.cpp
	)

add_executable( main ${STAM_SOURCES} )

target_link_libraries(main ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )


set_target_properties( main PROPERTIES COMPILE_FLAGS "-Wall -Woverloaded-virtual -Winit-self -Wno-sign-compare" )
