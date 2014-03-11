import os
from cffi import FFI
from numpy import *

try:
   from pydb import debugger

   ## Also add an exception hook.
   import pydb, sys
   sys.excepthook = pydb.exception_hook

except ImportError:
   import pdb
   def debugger():
       pdb.set_trace()

## Compile the library with:
'''
# OSX
gcc -fPIC \
    FitCurves.c GGVecLib.c \
    -dynamiclib -o FitCurves.dylib \
    -g -O3 -Wall -Wshadow -Wno-sign-compare

# Linux
gcc -fPIC \
    FitCurves.c GGVecLib.c \
    -shared -o FitCurves.so \
    -g -O2 -Wall -Wshadow -Wno-sign-compare

# Cygwin?
gcc -fPIC \
    FirCurve.c GGVecLib.c \
    -shared -o FitCurves.dll \
    -g -O2 -Wall -Wshadow -Wno-sign-compare
'''


ffi = FFI()
ffi.cdef("""
typedef struct Point2Struct {	/* 2d point */
	double x, y;
	} Point2;
typedef Point2 Vector2;

typedef Point2 *BezierCurve;

void (*DrawBezierCurve)(int n, BezierCurve curve);

void FitCurve(Point2 *d, int nPts, double error);
""")

import ctypes

def platform_shared_library_suffix():
    import sys
    result = '.so'
    if 'win' in sys.platform.lower(): result = '.dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = '.dylib'
    return result

libFitCurves = ffi.dlopen( os.path.join( os.path.dirname( __file__ ), 'FitCurves' + platform_shared_library_suffix() ) )

@ffi.callback("void(int,BezierCurve)")
def DrawBezierCurve( n, curve ):
    assert curve_in_progress is not None
    
    bezier = []
    for i in xrange( n+1 ):
        bezier.append( ( curve[i].x, curve[i].y ) )
    
    curve_in_progress.append( bezier )

libFitCurves.DrawBezierCurve = DrawBezierCurve

curve_in_progress = None
def FitCurve( vertices ):
    '''
    Given an N-by-2 numpy array 'vertices' of 2D vertices representing a line strip,
    returns an N-by-4-by-2 numpy.array of N cubic bezier curves approximating 'vertices'.
    '''
    
    import numpy
    
    global curve_in_progress
    assert curve_in_progress is None
    
    ## Make sure the input values have their data in a way easy to access from C.
    vertices = numpy.ascontiguousarray( numpy.asarray( vertices, dtype = ctypes.c_double ) )
    
    ## 'vertices' must be 2D
    assert vertices.shape[1] == 2
    
    ## This calls a callback function that appends to the global variable 'curve_in_progress'.
    curve_in_progress = []
    libFitCurves.FitCurve(
        ffi.cast( 'Point2*',  vertices.ctypes.data ),
        len( vertices ),
        1e-1
        )
    result = asarray( curve_in_progress )
    curve_in_progress = None
    
    return result

def test_simple( N = 10 ):
    print 'test_simple( %d )' % N
    
    from pprint import pprint
    
    assert N > 1
    
    line_strip = zeros( ( N, 2 ) )
    line_strip[:,0] = linspace( 0, 1, N )
    line_strip[:,1] = linspace( -1, 1, N )
    pprint( line_strip )
    
    beziers = FitCurve( line_strip )
    pprint( beziers )

def main():
    import sys
    
    N = 10
    if len( sys.argv ) > 1: N = int( sys.argv[1] )
    
    test_simple( N )

if __name__ == '__main__': main()
