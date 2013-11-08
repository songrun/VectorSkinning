import sys
import subprocess
import os
from numpy import asarray

#triangle_path = os.path.join( "C:\\Users\\Mai\\Dropbox\\Research\\Deformation\\src\\py\\triangle", "triangle.exe")
triangle_path = os.path.join( os.getcwd(), "triangle/triangle")

if not os.path.exists( triangle_path ):
    raise ImportError, "Triangle not found: " + triangle_path


def triangles_for_points( points, boundary_edges = None ):
    '''
    Given a sequence of 2D points 'points' and
    optional sequence of 2-tuples of indices into 'points' 'boundary_edges',
    returns a triangulation of the points as a sequence
    of length-three tuples ( i, j, k ) where i,j,k are
    the indices of the triangle's vertices in 'points'.
    
    If 'boundary_edges' is not specified or is an empty sequence,
    a convex triangulation will be returned.
    Otherwise, 'boundary_edges' indicates the boundaries of the desired mesh.
    '''
    
    import os, subprocess
    
    options = [ '-q', '-a100', '-g' ]
    
    if boundary_edges is None: boundary_edges = []
    
    if len( boundary_edges ) == 0:
        input_path = write_node_file( points )
        print triangle_path, input_path
        subprocess.call( [ triangle_path ] + options + [ input_path ] )
    else:
        input_path = write_poly_file( points, boundary_edges )
        subprocess.call( [ triangle_path ] + options + [ '-p', input_path ] )
    
    ele_path = os.path.splitext( input_path )[0] + '.1.ele'
    triangles = read_ele_file( ele_path )

    node_path = os.path.splitext( input_path )[0] + '.1.node'
    points = read_node_file( node_path)
    
    #os.remove( poly_path )
    #os.remove( ele_path )
    
    return points, triangles

def __write_node_portion_of_file_to_object( obj, points, boundary_indices = set() ):
    '''
    Given an object 'obj' that can be passed as a parameter to
    print >> 'obj', "Something to print",
    a sequence of 2D points 'points', and
    an optional set of indices in 'points' that are to be considered 'boundary_indices',
    writes the '.node' portion of the file suitable for passing to 'triangle'
    ( http://www.cs.cmu.edu/~quake/triangle.node.html ).
    Does not return a value.
    '''
    
    ## 'points' must be a non-empty sequence of x,y positions.
    points = asarray( points )
    assert points.shape == ( len( points ), 2 )
    assert points.shape[0] > 0
    ## The elements in 'boundary_indices' must be a subset of indices into 'points'.
    ## NOTE: set.issuperset() returns True if the sets are the same.
    assert set(range(len(points))).issuperset( boundary_indices )
    
    print >> obj, '## The vertices'
    print >> obj, len( points ), 2, 0, len( boundary_indices )
    for i, ( x, y ) in enumerate( points ):
        print >> obj, i, x, y, ( 1 if i in boundary_indices else 0 )

def write_poly_file( points, boundary_edges ):
    '''
    Given a sequence of 2D points 'points'
    and a potentially empty sequence 'boundary_edges' of
    2-tuples of indices into 'points',
    writes a '.poly' file suitable for passing to 'triangle'
    ( http://www.cs.cmu.edu/~quake/triangle.poly.html )
    and returns the path to the '.poly' file.
    '''
    
    ## Each of the two elements of each 2-tuple in 'boundary_edges'
    ## must be indices into 'points'.
    assert all([ i >= 0 and i < len( points ) and j >= 0 and j < len( points ) and i != j for i,j in boundary_edges ])
    ## They must be unique and undirected.
    assert len( boundary_edges ) == len( set([ frozenset( edge ) for edge in boundary_edges ]) )
    
    
    ## Create 'boundary_indices', the set of all indices that appear
    ## in 'boundary_edges'.
    boundary_indices = frozenset( asarray( boundary_edges ).ravel() )
    
    
    import tempfile
    ## This only works on Python 2.6+
    #poly_file = tempfile.NamedTemporaryFile( suffix = '.poly', delete = False )
    #poly_file_name = poly_file.name
    poly_file, poly_file_name = tempfile.mkstemp( suffix = '.poly' )
    poly_file = os.fdopen( poly_file, 'w' )
    
    print >> poly_file, '## Written by triangle.py'
    
    __write_node_portion_of_file_to_object( poly_file, points, boundary_indices )
    
    print >> poly_file, ''
    print >> poly_file, '## The segments'
    print >> poly_file, len( boundary_edges ), len( boundary_edges )
    for i, ( e0, e1 ) in enumerate( boundary_edges ):
        print >> poly_file, i, e0, e1, 1
    
    print >> poly_file, ''
    print >> poly_file, '## The holes'
    print >> poly_file, 0
    
    poly_file.close()
    return poly_file_name

def write_node_file( points ):
    '''
    Given a sequence of 2D points 'points',
    writes a '.node' file suitable for passing to 'triangle'
    ( http://www.cs.cmu.edu/~quake/triangle.node.html )
    and returns the path to the '.node' file.
    '''
    
    import tempfile
    ## This only works on Python 2.6+
    #node_file = tempfile.NamedTemporaryFile( suffix = '.node', delete = False )
    #node_file_name = node_file.name
    node_file, node_file_name = tempfile.mkstemp( suffix = '.node' )
    node_file = os.fdopen( node_file, 'w' )
    
    print >> node_file, '## Written by triangle.py'
    
    __write_node_portion_of_file_to_object( node_file, points )
    
    node_file.close()
    return node_file_name

def read_ele_file( ele_path ):
    '''
    Reads a '.ele' file generated by 'triangle'.
    Returns the list of triangles as indices into the
    corresponding '.node' file.
    '''
    
    ele_file = open( ele_path )
    
    ## Ignore top line.
    ele_file.readline()
    
    triangles = []
    for line in ele_file:
        sline = line.strip().split()
        if len( sline ) == 0: continue
        if sline[0][0] == '#': continue
        
        triangles.append( tuple([ int( index ) for index in sline[1:4] ]) )
        assert len( triangles[-1] ) == 3
    
    ele_file.close()
    
    return triangles

def read_node_file( node_path ):
    '''
    Reads a '.node' file generated by 'triangle'.
    Returns the list of points as tuples.
    '''
    
    node_file = open( node_path )
    
    ## Ignore top line.
    node_file.readline()
    
    triangles = []
    for line in node_file:
        sline = line.strip().split()
        if len( sline ) == 0: continue
        if sline[0][0] == '#': continue
        
        triangles.append( tuple([ float( index ) for index in sline[1:4] ]) )
        #assert len( triangles[-1] ) == 3
    
    node_file.close()
    
    return triangles

# def main():
#     pts = [ ( -1,-1 ), ( 1, -1 ), ( 1, 1 ), ( -1, 1 ), ( 0, 0 ) ]
#     edges = [ ( 0, 1 ), ( 1, 2 ), ( 2, 3 ), ( 3, 0 ) ]
#     
#     ## This isn't very good, because 4 random points may be self-intersecting
#     ## when viewed as a polyline loop.
#     #import random
#     #pts = [ ( random.uniform( -1, 1 ), random.uniform( -1, 1 ) ) for i in xrange(4) ]
#     
#     print 'pts:', pts
#     
#     points, triangles = triangles_for_points( pts )
#     print 'points (no boundary edges):', points
#     print 'triangles (no boundary edges):', triangles
#     
#     print 'width edges:', edges
#     points, triangles = triangles_for_points( pts, edges )
#     print 'points (with edges):', points
#     print 'triangles (with edges):', triangles
# 
# if __name__ == '__main__': main()
