from math import *
from numpy import *

'''
gcc raytri.c raytri_wrapper.cpp -lstdc++ -fkeep-inline-functions \
    -dynamiclib -o libraytri.dylib \
    -g -O3 -Wall -Wshadow -Woverloaded-virtual -Winit-self -Wno-sign-compare
'''

try:
    import ctypes
except:
    print '''
ERROR: ctypes not installed properly. (ctypes needed for laplacian_editing_1d())
        '''
    import sys
    sys.exit()

def platform_shared_library_suffix():
    import sys
    result = 'so'
    if 'win' in sys.platform.lower(): result = 'dll'
    ## No else if, because we want darwin to override win (which is a substring of darwin)
    if 'darwin' in sys.platform.lower(): result = 'dylib'
    return result

libraytri = ctypes.cdll.LoadLibrary( 'libraytri.' + platform_shared_library_suffix() )
real_t = ctypes.c_float
#real_t = ctypes.c_double

def to_ctypes_array( sequence, ctypes_type ):
    '''
    Returns a contiguous ctypes array suitable for extracting
    a pointer (with .ctypes.data_as() or ctypes_array_ref()) to pass to a C function.
    '''
    return ascontiguousarray( asarray( sequence, dtype = ctypes_type ) )
def ctypes_array_ref( ctypes_array, ctypes_type ):
    '''
    Returns a ctypes.POINTER to the internal data of 'ctypes_array' of type 'ctypes_type'.
    
    WARNING: Do not chain the calls to to_ctypes_array() and ctypes_array_ref() as in
               ctypes_array_ref( to_ctypes_array( sequence, ctypes_type ), ctypes_type ).
             ctypes_array_ref returns a pointer to the internal numpy storage of the array and
             cannot keep a reference to the array, so it will be garbage collected if these
             functions are chained together.
             Instead, keep a temporary reference to the result of to_ctypes_array():
               contig = to_ctypes_array( sequence, ctypes_type )
               ctypes_array_ref( contig, ctypes_type )
    '''
    return ctypes_array.ctypes.data_as( ctypes.POINTER( ctypes_type ) )

__intersect_triangle1 = libraytri.intersect_triangle1
__intersect_triangle1.argtypes = \
    [
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double )
    ]
__intersect_triangle1.restype = ctypes.c_int
def ray_triangle_intersection( ray, tri ):
    '''
    Given a ray (3d point, 3d direction) and 'triangle' (a list of three 3d points),
    returns a vector (t u v), where t is the distance to the plane in which the triangle lies
    and (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the triangle, returns None.
    NOTE: t may be negative, indicating that the intersection is in the negative ray direction.
    
    NOTE: The exact intersection point can be reproduced in two ways:
          ray[0] + t * ray[1]
          or
          (1-u-v)*tri[0] + u*tri[1] + v*tri[2]
    
    used
    '''
    
    ## Convert parameters into a ctypes-friendly format.
    rorig = ascontiguousarray( asarray( ray[0], dtype = ctypes.c_double ) )
    assert len( rorig ) == 3
    rdir = ascontiguousarray( asarray( ray[1], dtype = ctypes.c_double ) )
    assert len( rdir ) == 3
    v0 = ascontiguousarray( asarray( tri[0], dtype = ctypes.c_double ) )
    v1 = ascontiguousarray( asarray( tri[1], dtype = ctypes.c_double ) )
    v2 = ascontiguousarray( asarray( tri[2], dtype = ctypes.c_double ) )
    assert len( v0 ) == 3
    assert len( v1 ) == 3
    assert len( v2 ) == 3
    
    ## Make the output parameters in a ctypes-friendly way.
    t = ctypes.c_double(0)
    u = ctypes.c_double(0)
    v = ctypes.c_double(0)
    
    intersect = __intersect_triangle1(
        rorig.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        rdir.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v0.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v1.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        v2.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) ),
        ctypes.byref( t ),
        ctypes.byref( u ),
        ctypes.byref( v )
        )
    
    ## TODO Q: Why is the return value getting converted to an 'int' instead of remaining a c_int?
    # if intersect.value:
    if intersect:
        return t.value, u.value, v.value
    else:
        return None

def ray_mesh_intersections_slow( ray, mesh ):
    '''
    Given a ray (3d point, 3d direction) and 'mesh' (a structure containing a list of
    3d vertices 'vs' and a list of triangles (triplets of indices into 'vs') named 'faces'),
    returns a list of vectors ( t, fi, (u,v) ), where
        t is the distance to the plane in which the triangle lies,
        fi is the index of the face in mesh.faces, and
        (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the mesh, returns an empty list.
    
    NOTE: The exact intersection point can be reproduced in two ways:
          ray[0] + t * ray[1]
          or
          (1-u-v)*p[0] + u*p[1] + v*p[2], where the p's are the three vertices of triangle fi.
    
    used
    '''
    
    result = []
    for fi, face in enumerate( mesh.faces ):
        isect = ray_triangle_intersection(
            ray,
            [ mesh.vs[ vi ] for vi in face ]
            )
        if isect is None: continue
        result.append( ( isect[0], fi, ( isect[1], isect[2] ) ) )
    
    return result

class raymesh_intersection_t( ctypes.Structure ):
    _fields_ = [
        ( "t", ctypes.c_double ),
        ( "fi", ctypes.c_int ),
        ( "u", ctypes.c_double ),
        ( "v", ctypes.c_double )
        ]
__raytri_wrapper_ray_mesh_intersections = libraytri.ray_mesh_intersections
__raytri_wrapper_ray_mesh_intersections.argtypes = \
    [
        ctypes.POINTER( ctypes.c_double ),
        ctypes.POINTER( ctypes.c_double ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_double ),
        ctypes.c_int,
        ctypes.POINTER( ctypes.c_int ),
        ctypes.POINTER( ctypes.c_int )
    ]
__raytri_wrapper_ray_mesh_intersections.restype = ctypes.POINTER( raymesh_intersection_t )
__raytri_wrapper_ray_mesh_intersections_free = libraytri.ray_mesh_intersections_free
__raytri_wrapper_ray_mesh_intersections_free.argtypes = [ ctypes.POINTER( raymesh_intersection_t ) ]
__raytri_wrapper_ray_mesh_intersections_free.restype = None
def _raytri_wrapper_ray_mesh_intersections( ray, mesh ):
    '''
    Given a ray (3d point, 3d direction) and 'mesh' (a structure containing a list of
    3d vertices 'vs' and a list of triangles (triplets of indices into 'vs') named 'faces'),
    returns a list of vectors ( t, fi, (u,v) ), where
        t is the distance to the plane in which the triangle lies,
        fi is the index of the face in mesh.faces, and
        (u, v) represents the coordinates inside the triangle.
    If the ray does not intersect the mesh, returns an empty list.
    
    NOTE: The exact intersection point can be reproduced in two ways:
          ray[0] + t * ray[1]
          or
          (1-u-v)*p[0] + u*p[1] + v*p[2], where the p's are the three vertices of triangle fi.
    
    used
    '''
    
    ray = (
        to_ctypes_array( ray[0], ctypes.c_double ),
        to_ctypes_array( ray[1], ctypes.c_double )
        )
    vs = to_ctypes_array( mesh.vs, ctypes.c_double )
    faces = to_ctypes_array( mesh.faces, ctypes.c_int )
    
    assert ray[0].shape == (3,)
    assert ray[1].shape == (3,)
    assert vs.shape[1] == 3
    assert faces.shape[1] == 3
    
    num_intersections_out = ctypes.c_int()
    intersections = __raytri_wrapper_ray_mesh_intersections(
        ctypes_array_ref( ray[0], ctypes.c_double ),
        ctypes_array_ref( ray[1], ctypes.c_double ),
        len( vs ),
        ctypes_array_ref( vs, ctypes.c_double ),
        len( faces ),
        ctypes_array_ref( faces, ctypes.c_int ),
        ctypes.byref( num_intersections_out )
        )
    num_intersections_out = num_intersections_out.value
    result = [
        ( isect.t, isect.fi, ( isect.u, isect.v ) )
        for isect in intersections[:num_intersections_out]
        ]
    __raytri_wrapper_ray_mesh_intersections_free( intersections )
    
    ## These generate the same thing:
    #print 'ray_mesh_intersections fast', result
    #print 'ray_mesh_intersections slow', ray_mesh_intersections_slow( ray, mesh )
    #print 'ray_mesh_intersections difference:', set( result ).symmetric_difference( set( ray_mesh_intersections_slow( ray, mesh ) ) )
    
    return result

ray_mesh_intersections = _raytri_wrapper_ray_mesh_intersections

def point2d_in_mesh2d_barycentric( point2d, vertices2d, faces ):
    '''
    Given a 2d point 'point2d' and mesh in the form of a list of
    2d vertices 'vertices2d' and a list of triangles (triplets
    of indices into 'vertices2d') named 'faces',
    returns a tuple ( fi, (b0, b1, b2) ), where
        fi is the index of the face in mesh.faces, and
        (b0, b1, b2) represents the barycentric coordinates inside the triangle.
    If the point is not inside the mesh, returns None.
    
    tested (see test_one_triangle(), below)
    '''
    
    assert len( point2d ) == 2
    assert len( vertices2d[0] ) == 2
    
    ## Make a ray above the mesh pointing down
    point = ( point2d[0], point2d[1], 1 )
    dir = ( 0, 0, -1 )
    
    vertices = zeros( ( len( vertices2d ), 3 ) )
    vertices[:,:2] = vertices2d
    
    class Object( object ): pass
    
    mesh = Object()
    mesh.vs = vertices
    mesh.faces = faces
    
    intersections = ray_mesh_intersections( ( point, dir ), mesh )
    if len( intersections ) == 0: return None
    ## Otherwise, return the intersection with the largest smallest barycentric coordinate,
    ## so the "most" inside the triangle. This should give a reasonable result
    ## with regards to negative epsilon issues.
    intersections = [ ( fi, ( 1.-u-v, u, v ) ) for ( t, fi, ( u,v ) ) in intersections ]
    smallest = asarray([ bary for fi, bary in intersections ]).min( axis = 1 ).argmax()
    return intersections[ smallest ]

def closest_distsqr_and_edge_index_and_t_on_line_strip_to_point( line_strip, pt ):
    '''
    Given a line strip (sequence of n-dimensional points) 'linestrip' and n-dimensional point 'pt',
    returns the tuple (
        squared distance to the closest point on 'line_strip',
        index of point in 'line_strip' where the edge (index, index+1) contains the closest point ),
        t along the line strip such that the closest point is (1-t)*line_strip[index] + t*(line_strip[index+1])
        ).
    
    tested (see test_line_strip_and_loop_distances())
    '''
    
    assert len( line_strip ) > 0
    
    from edge_distances import min_distanceSqr_edge_t_to_edges
    ## This function takes a sequence of points, so we have to pack 'pt' into a length-1 list.
    ## It also takes edges as pairs of points, not as a line loop, so we have to duplicate points.
    result = min_distanceSqr_edge_t_to_edges( [ pt ], list( zip( line_strip[:-1], line_strip[1:] ) ) )
    ## Unpack the result, since we passed a length-1 list containing 'pt'.
    return result[0][0], result[1][0], result[2][0]

def closest_distsqr_and_edge_index_and_t_on_line_loop_to_point( line_loop, pt ):
    '''
    Same as closest_distsqr_and_point_and_edge_index_on_line_strip_to_point(), but
    takes a line loop (closed path) instead of a line strip (open path).
    
    NOTE: The index of the closing edge is len( line_loop )-1.
    
    tested (see test_line_strip_and_loop_distances())
    '''
    ## Append the first point to the end to turn the line loop into a line strip.
    return closest_distsqr_and_edge_index_and_t_on_line_strip_to_point(
        list( line_loop ) + [line_loop[0]], pt
        )

def test_one_triangle( sx = None, sy = None ):
    if sx is None: sx = 1.
    if sy is None: sy = 1.
    
    vertices = [ (0,0), (0,sy), (sx,0) ]
    faces = [ (0,1,2) ]
    
    print 'vertices:', vertices
    print 'faces:', faces
    
    pts = vertices + [ ( .1,.1 ), ( -.1, -.1 ), ( .5, 0 ), ( .5*sx, 0 ), ( 0, .5 ), ( 0, .5*sy ) ]
    
    for pt in pts:
        print 'point:', pt
        bary = point2d_in_mesh2d_barycentric( pt, vertices, faces )
        print 'bary:', bary

def test_line_strip_and_loop_distances():
    pt = (0,0)
    line_strip = [ ( 0, -.1 ), ( 1, -.1 ), ( 1, .9 ), ( 0, .9 ) ]
    print 'pt:', pt
    print 'line_strip:', line_strip
    print 'closest_distsqr_and_edge_index_and_t_on_line_strip_to_point():'
    print closest_distsqr_and_edge_index_and_t_on_line_strip_to_point( line_strip, pt )
    
    ## Interpreted as a line loop, 'line_strip' should pass directly through 'pt'.
    print 'closest_distsqr_and_edge_index_and_t_on_line_loop_to_point():'
    print closest_distsqr_and_edge_index_and_t_on_line_loop_to_point( line_strip, pt )

def main():
    import sys
    
    sx = None
    sy = None
    if len( sys.argv ) == 3:
        sx = float( sys.argv[1] )
        sy = float( sys.argv[2] )
    
    test_one_triangle( sx, sy )
    # test_line_strip_and_loop_distances()

if __name__ == '__main__': main()
