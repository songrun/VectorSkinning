import os
from cffi import FFI

## Compile the library with:
'''
# OSX
g++ -fPIC \
    bbw.cpp \
    -I/Users/yotam/Work/ext/libigl/include \
    -I/Installers/eigen/eigen-3.2-ffa86ffb5570/eigen-eigen-ffa86ffb5570 \
    -I/Users/yotam/Work/ext/mosektoolsosx64x86/mosek/7/tools/platform/osx64x86/h \
    -L/Users/yotam/Work/ext/mosektoolsosx64x86/mosek/7/tools/platform/osx64x86/bin \
    -lmosek \
    -dynamiclib -o bbw.dylib \
    -g -O2 -Wall -Wshadow -Wno-sign-compare


# Linux
g++ -fPIC \
    bbw.cpp \
    -Ipath/to/igl????
    -shared -o bbw.so \
    -g -O2 -Wall -Wshadow -Wno-sign-compare

# Cygwin?
g++ -fPIC \
    bbw.cpp \
    -Ipath/to/igl????
    -shared -o bbw.dll \
    -g -O2 -Wall -Wshadow -Wno-sign-compare
'''


ffi = FFI()
ffi.cdef("""
typedef double real_t;

// Returns 0 for success, anything else is an error.
int bbw2D(
    /// Input Parameters
    // 'vertices' is a pointer to num_vertices*2 floating point values,
    // packed: x0, y0, x1, y1, ...
    int num_vertices, real_t* vertices,
    // 'faces' is a pointer to num_faces*3 integers,
    // where each face is three vertex indices: f0.v0, f0.v1, f0.v2, f1.v0, f1.v1, f1.v2, ...
    // Face i's vertices are: vertices[ faces[3*i]*2 ], vertices[ faces[3*i+1]*2 ], vertices[ faces[3*i+2]*2 ]
    int num_faces, int* faces,
    // 'skeleton_vertices' is a pointer to num_skeleton_vertices*2 floating point values,
    // packed the same way as 'vertices'
    int num_skeleton_vertices, real_t* skeleton_vertices,
    // 'skeleton_point_handles' is a pointer to num_skeleton_point_handles integers,
    // where each element "i" in skeleton_point_handles references the vertex whose data
    // is located at skeleton_vertices[ skeleton_point_handles[i]*2 ].
    int num_skeleton_point_handles, int* skeleton_point_handles,
    // TODO: Take skeleton bone edges and cage edges
    
    /// Output Parameters
    // 'Wout' is a pointer to num_vertices*num_skeleton_vertices values.
    // Upon return, W will be filled with each vertex in 'num_vertices' weight for
    // each skeleton vertex in 'num_skeleton_vertices'.
    // The data layout is that all num_vertices values for skeleton vertex 0
    // appear before all num_vertices values for skeleton vertex 1, and so on.
    // TODO: This data layout may be transposed, please check.
    real_t* Wout
    );
""")

def platform_shared_library_suffix():
    import sys
    result = '.so'
    if 'win' in sys.platform.lower(): result = '.dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = '.dylib'
    return result

libbbw = ffi.dlopen( os.path.join( os.path.dirname( __file__ ), 'bbw' + platform_shared_library_suffix() ) )

def bbw( vertices, faces, skeleton_handle_vertices, skeleton_point_handles ):
    '''
    Given an N-by-2 numpy array 'vertices' of 2d vertices,
    an M-by-3 numpy array 'faces' of indices into 'vertices',
    an H-by-2 numpy.array 'skeleton_handle_vertices' of 2d vertices,
    a numpy array 'skeleton_point_handles' of indices into 'skeleton_handle_vertices'
    which are the point handles,
    returns a N-by-H numpy.array of weights per vertex per handle.
    '''
    
    import numpy
    
    ## Make sure the input values have their data in a way easy to access from C.
    vertices = numpy.ascontiguousarray( asarray( vertices, dtype = numpy.float64 ) )
    faces = numpy.ascontiguousarray( asarray( faces, dtype = int ) )
    skeleton_handle_vertices = numpy.ascontiguousarray( asarray( skeleton_handle_vertices, dtype = numpy.float64 ) )
    skeleton_point_handles = numpy.ascontiguousarray( asarray( skeleton_point_handles, dtype = int ) )
    
    Wout = numpy.empty( ( len( vertices ), len( skeleton_handle_vertices ) ), dtype = float64 )
    result = libbbw.bbw2D(
        len( vertices ), vertices.ctypes.data,
        len( faces ), faces.ctypes.data,
        len( skeleton_handle_vertices ), skeleton_handle_vertices.ctypes.data,
        len( skeleton_point_handles ), skeleton_point_handles.ctypes.data,
        Wout.ctypes.data
        )
    if result != 0:
        raise RuntimeError
    
    return Wout
