from numpy import *

def distancesSqr_and_t_to_edges( pts, edges ):
    '''
    Input parameter 'pts' has dimensions #pts x (x,y,...).
    Input parameter 'edges' has dimensions ... x #edges x 2 endpoints x N coordinates (x,y,...).
    Returns an array of distances squared with dimensions #edges x #pts.
    '''
    
    pts = asfarray( pts ).T
    ## pts has dimensions N (x,y,...) x #pts
    edges = asfarray( edges )
    ## edges has dimensions ... x #edges x 2 endpoints x N coordinates (x,y,...)
    
    N = pts.shape[0]
    
    assert len( pts.shape ) == 2 and pts.shape[0] == N
    assert edges.shape[-2] == 2 and edges.shape[-1] == N
    #print 'pts.shape:', pts.shape
    #print 'edges.shape:', edges.shape
    
    
    ## get distance squared to each edge:
    ##   let p = black_pixel_pos, a = endpoint0, b = endpoint1, d = ( b-a ) / dot( b-a,b-a )
    ##   dot( p-a, d ) < 0 => dot( p-a, p-a )
    ##   dot( p-a, d ) > 1 => dot( p-b, p-b )
    ##   else              => dot( dot( p-a, d ) * (b-a) - p, same )
    p_a = pts[newaxis,...] - edges[...,:,0,:,newaxis]
    p_b = pts[newaxis,...] - edges[...,:,1,:,newaxis]
    ## p_a and p_b have dimensions ... x #edges x N coordinates (x,y,...) x #pts
    b_a = edges[...,:,1,:] - edges[...,:,0,:]
    ## b_a has dimensions ... x #edges x N coordinates (x,y,...)
    d = b_a / ( b_a**2 ).sum( -1 )[...,newaxis]
    ## d has same dimensions as b_a
    assert b_a.shape == d.shape
    cond = ( p_a * d[...,newaxis] ).sum( -2 )
    ## cond has dimensions ... x #edges x #pts
    assert cond.shape[-2:] == (edges.shape[-3], pts.shape[-1])
    
    ## clip cond so that values are between 0 and 1
    cond.clip( 0, 1, cond )
    
    #distancesSqr = empty( cond.shape, Real )
    ## distancesSqr has dimensions ... x #edges x #pts
    #assert distancesSqr.shape[-2:] == (edges.shape[-3], pts.shape[-1])
    
    # distancesSqr = p_a - cond[:,newaxis,:] * b_a[...,newaxis]
    # distancesSqr = ( distancesSqr**2 ).sum( 1 )
    # <=>
    distancesSqr = ( ( p_a - cond[:,newaxis,:] * b_a[...,newaxis] )**2 ).sum( -2 )
    
    #print 'distancesSqr:', distancesSqr
    #print 'distances:', sqrt( distancesSqr.min(0) )
    
    return distancesSqr, cond

def min_distanceSqr_edge_t_to_edges( pts, edges ):
    '''
    Input parameter 'pts' has dimensions #pts x 2 (x,y).
    Input parameter 'edges' has dimensions ... x #edges x 2 endpoints x 2 coordinates (x,y).
    Returns an array of the minimum distance to edges with dimensions #pts.
    '''
    
    distancesSqr, cond = distancesSqr_and_t_to_edges( pts, edges )
    edge_index = distancesSqr.argmin( -2 )
    return distancesSqr[ edge_index ], edge_index, cond[ edge_index ]

def test_timing():
    import timeit, time
    
    N = 10
    dim = 3
    npts = 300
    nedges = 100
    
    #random.uniform( 0, 1 )
    
    setup2d_basic = \
'''
from numpy import asarray, array
from edge_distances import distancesSqr_and_t_to_edges
pts = [(.5, 0.), (.5,1.), (0.,.25), (1.,.25), (0.,0.), (1.,0.), (.3,.3)]
pts = asarray( pts )
edges = [[(.25, .25), (.75, .25)], [(.25, .25), (.25, .75)], [(.75, .25), (.75, .75)], [(.25, .75), (.75, .75)]]
'''
    
    setupNd = \
'''
from numpy import asarray, array
from edge_distances import distancesSqr_and_t_to_edges
import random
r = random.Random( 7 )
u = r.uniform
dim = %s
npts = %s
nedges = %s
pts = [ [ u( 0,1 ) for d in range(dim) ] for i in xrange( npts ) ]
pts = asarray( pts )
edges = [ [ [ u( 0,1 ) for d in range(dim) ], [ u( 0,1 ) for d in range(dim) ] ] for i in xrange( nedges ) ]
''' % (dim, npts, nedges)
    
    print 'repititions:', N
    print 'dimension:', dim
    print 'num points:', npts
    print 'num edges:', nedges
    
    #timeit.Timer( 'print abs( distancesSqr_to_edges( pts, edges ) - distancesSqr_to_edges_Nd( pts, edges ) ).sum(),', setup2d, time.clock ).timeit(100)
    #print 
    
    ## I have to interleave calls to the two functions I'm comparing, since whichever executes
    ## second has an advantage.  I verified this by calling identical functions.
    
    
    print 'distancesSqr_to_edges:',
    print timeit.Timer( 'distancesSqr_and_t_to_edges( pts, edges )', setupNd ).repeat(4,N)

def test1():
    pt = [ 0,0 ]
    edges = [
        ## distance should be 0 at t = 0
        [ (0,0), (1,0) ],
        [ (0,0), (0,1) ],
        ## distance should be 0 at t = 1
        [ (1,0), (0,0) ],
        [ (0,1), (0,0) ],
        ## distance should be 1 at t = 0
        [ (0,1), (1,1) ],
        ## distance should be 1 at t = 1
        [ (1,1), (0,1) ],
        ## distance should be .3 at t = .75
        [ (-.75,.3), (.25,.3) ],
        ## distance should be 0 at t = .5
        [ (-1,0), (1,0) ],
        ## distance should be 0 at t = .25
        [ (-.25,0), (.75,0) ],
        ]
    distSqrs, ts = distancesSqr_and_t_to_edges( [ pt ], edges )
    print 'distSqrs:'
    print distSqrs
    print 'ts:'
    print ts
    minDistSqr, min_edge_index, min_t = min_distanceSqr_edge_t_to_edges( [ pt ], edges )
    print 'min_distanceSqr_edge_t_to_edges():', minDistSqr, min_edge_index, min_t

def test2():
    pts = [(.5, 0.), (.5,1.), (0.,.25), (1.,.25), (0.,0.), (1.,0.), (.3,.3)]
    pts = asarray( pts )
    edges = [[(.25, .25), (.75, .25)], [(.25, .25), (.25, .75)], [(.75, .25), (.75, .75)], [(.25, .75), (.75, .75)]]
    distSqrs, ts = distancesSqr_and_t_to_edges( pts, edges )
    print 'distSqrs:'
    print distSqrs
    print 'ts:'
    print ts
    print 'min_distanceSqr_edge_t_to_edges():', min_distanceSqr_edge_t_to_edges( pts, edges )

def main():
    test1()
    #test_timing()

if __name__ == '__main__':
    main()
