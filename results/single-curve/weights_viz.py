from numpy import *

def shepard_fast( vs, skeleton_handle_vertices ):
	'''
	Given an N-by-(2 or 3) sequence 'vs' of 2D or 3D vertices and
	an H-by-(2 or 3) sequence 'skeleton_handle_vertices' of 2D or 3D vertices,
	returns a N-by-H numpy.array of weights per vertex per handle.
	
	tested
	
	>>> skeleton_handle_vertices = [ (0,0), (1,0), (1,1), (0,1) ]
	>>> vs = [ (0,0), (1,0), (1,1), (0,1), (.2,.1), (.9,.8), (.8,.9), ( -1, -1 ), ( -1, 1 ) ]
	>>> shepard_fast( vs, skeleton_handle_vertices ) - shepard( vs, skeleton_handle_vertices )
	>>> abs( shepard_fast( vs, skeleton_handle_vertices ) - shepard( vs, skeleton_handle_vertices ) ).max()
	'''
	
	vs = asarray( vs )
	skeleton_handle_vertices = asarray( skeleton_handle_vertices )
	
	assert len( vs.shape ) == 2
	assert len( skeleton_handle_vertices.shape ) == 2
	assert vs.shape[1] == skeleton_handle_vertices.shape[1]
	
	diffs = (skeleton_handle_vertices[newaxis,...] - vs[:,newaxis,:])
	diffs = ( diffs**2 ).sum( -1 )
	## 'diffs' is N-by-H
	
	small = diffs < 1e-8
	
	## NOTE: Shepard weights could allow one to take this to a power other than 2,
	##       which is the implied power since 'diffs' stores the squared distance.
	diffs = 1./diffs
	
	mask = small.any( 1 )[:,newaxis].repeat( diffs.shape[1], 1 )
	diffs[ mask ] = small[ mask ]
	
	## Normalize rows.
	diffs /= diffs.sum( 1 )[:,newaxis]
	
	return diffs

xdpi = 9.
artboard = zeros( ( 30*xdpi, 100*xdpi, 3 ) )
w = where( ones( artboard.shape[:2], dtype = bool ) )
weights = shepard_fast( array( w ).T, xdpi*array([ ( 20, 10 ), ( 20, 90 ), ( 12.585, 50 ) ]) )
artboard[:,:,0][ w ] = weights[:,2]
## Other handles green
#artboard[:,:,1:2] = ( 1-artboard[:,:,0] )[...,newaxis]
## Other handles cyan
artboard[:,:,1:] = ( 1-artboard[:,:,0] )[...,newaxis]

import Image

Image.fromarray( ( artboard * 255. ).round().clip( 0, 255 ).astype( uint8 ) ).save( 'weights.png' )
