from bezier_utility import *
from triangle import *
import bbw_wrapper.bbw as bbw
from itertools import izip as zip

def uniquify_points_and_return_input_index_to_unique_index_map( pts, boundary_pts ):
	'''
	Given a sequence of N points 'pts',
	returns two items:
	   a sequence of all the unique elements in 'pts'
	   and
	   a list of length N where the i-th item in the list tells you where
	   pts[i] can be found in the unique elements.
	'''
	
	kRound = 5
	from collections import OrderedDict
	unique_pts = OrderedDict()
	pts_map = []
	## Add rounded points to a dictionary and set the key to the index into the ordered dictionary.
	for i, pt in enumerate( map( tuple, asarray( pts ).round( kRound ) ) ):
		index = unique_pts.setdefault( pt, len( unique_pts ) )
		pts_map.append( index )
	
	boundary_map = []
	for pt in map( tuple, asarray( boundary_pts ).round( kRound ) ):
		boundary_map.append( unique_pts[ pt ] )
		## The above will raise a KeyError if boundary_pts contains points not in all_pts.
		# raise RuntimeError('boundary_pts contains points not in all_pts')
	
	return unique_pts.keys(), pts_map, boundary_map

def triangulate_and_compute_weights(boundary_pts, skeleton_handle_vertices, all_pts=None):
	'''
	trianglue a region closed by a bunch of bezier curves, precompute the vertices at each sample point.
	'''	
	if all_pts is None:
		all_pts = [boundary_pts]
	
	boundary_pts = asarray( boundary_pts )
	all_pts = asarray( all_pts )
	all_shapes = [ asarray( pts ).shape[:-1] for pts in all_pts ]
	
	boundary_pts = concatenate( boundary_pts )
	all_pts = concatenate( [ concatenate( curve_pts ) for curve_pts in all_pts ] )
	
	print 'Removing duplicate points...'
	all_clean_pts, pts_maps, boundary_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts, boundary_pts )
	print '...finished.'
	
	boundary_edges = [ ( i, (i+1) % len( set(boundary_maps) ) ) for i in xrange(len( set(boundary_maps) )) ]

# 	boundary_edges = [ ( boundary_maps[i], boundary_maps[ (i+1) % len(boundary_maps) ] ) for i in xrange(len( boundary_maps )) ]
	
	all_maps = []
	pos = 0
	for shape in all_shapes:
		maps = []
		for j in range( shape[0] ):
			maps.append( pts_maps[ pos : pos + shape[1] ] )
			pos += shape[1]
		all_maps.append( maps )
	
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	
	if len( skeleton_handle_vertices ) > 0:
		skeleton_handle_vertices = asarray( skeleton_handle_vertices )[:, :2]
	skeleton_point_handles = list( range( len(skeleton_handle_vertices) ) )
	
	registered_pts = concatenate( ( all_clean_pts, skeleton_handle_vertices ), axis = 0 )
	print 'Computing triangulation...'
	vs, faces = triangles_for_points( registered_pts, boundary_edges )
	print '...finished.'
	
	vs = asarray(vs)[:, :2] 
	faces = asarray(faces)

	print 'Computing BBW...'
	all_weights = bbw.bbw(vs, faces, skeleton_handle_vertices, skeleton_point_handles)
	print '...finished.'
	
	return vs, all_weights, all_maps

def precompute_W_i_bbw( vs, weights, i, sampling_index2vs_index, sampling, ts, dts = None ):
	'''
	Given an N-by-k numpy.array 'vs' of all points represented in 'weights',
	an N-by-num-handles numpy.array 'weights' of all the weights for each sample vertex,
	an index 'i' specifying which handle,
	a map 'sampling_index2vs_index' from the index of a point in 'sampling' to the index of its closest vertex in 'vs',
	an M-by-k numpy.array 'sampling' containing sampled positions,
	a length-M numpy.array of t values corresponding to each sample in 'sampling',
	an optional length-M numpy.array of dt values corresponding to each sample in 'sampling',
	returns W, a 4-by-4 numpy.array defined as:
		\int weights_i( sample ) \overbar{t}^T \overbar{t}^T dt
	where sample is drawn from 'sampling', t is drawn from 'ts', and dt is drawn from 'dts'.
	
	If 'dts' is not given, it defaults to 1/len(sampling).
	'''
	
	if dts is None: dts = ones( len( sampling )-1 ) * (1./(len(sampling)-1) )
	
	### Asserts
	## Ensure our inputs are numpy.arrays:
	weights = asarray( weights )
	vs = asarray( vs )
	sampling = asarray( sampling )
	ts = asarray( ts )
	dts = asarray( dts )
	
	assert len( weights.shape ) == 2
	assert len( vs.shape ) == 2
	assert len( sampling.shape ) == 2
	assert len( ts.shape ) == 1
	assert len( dts.shape ) == 1

	result = zeros((4,4 ))

	## Vertices and sampling must have the same dimension for each point.
#	debugger()
#	sampling = sampling[:,:-1]
	assert vs.shape[1] == sampling.shape[1]
	
	## The index 'i' must be valid.
	assert i >= 0 and i < weights.shape[1]
	
	def weight_function( pi ):
		## Find the closest vertex in 'vs' to 'p'
		#vi = argmin( ( ( vs - p )**2 ).sum( axis = 1 ) )
		vi = sampling_index2vs_index[ pi ]
		# assert allclose( vs[vi], p, 1e-5 )
		return weights[ vi, i ]
	
	result[:] = precompute_W_i_with_weight_function_and_sampling( weight_function, sampling, ts, dts )		
				
	return result

# def precompute_W_i_with_weight_function_and_sampling( weight_function, sampling, ts, dts ):
#	'''
#	Given a function 'weight_function' that takes a point and returns its weight,
#	a N-by-k numpy.array 'sampling' containing the positions of the control points as the rows,
#	corresponding t values 'ts' for each point in 'sampling',
#	and an optional corresponding 'dt' for each point in sampling (default is 1/len(sampling)),
#	returns W, a 4-by-4 numpy.array defined as:
#		\int_i weight_function( sample ) \overbar{t}^T \overbar{t}^T dt
#	where sample, t, and dt are drawn from the corresponding input arrays.
#	
#	The optional parameter 'num_samples' determines how many samples to use to compute
#	the integral.
#	'''
#	
#	### Asserts
#	## Ensure our inputs are the same lengths:
#	assert len( sampling ) == len( ts )
#	assert len( sampling ) == len( dts ) + 1
#	
#	## Ensure our inputs are numpy.arrays:
#	sampling = asarray( sampling )
#	ts = asarray( ts )
#	dts = asarray( dts )
#	
#	## sampling must be N-by-k.
#	assert len( sampling.shape ) == 2
#	assert len( ts.shape ) == 1
#	assert len( dts.shape ) == 1
#	
#	### Compute the integral.
#	W_i = zeros( ( 4,4 ) )
#	tbar = ones( 4 )
#	
#	for i in range(len(dts)):
#		t = (ts[i] + ts[i+1])/2
#		dt = dts[i]
#		## For arc-length parameterization:
#		# dt = magnitude( sampling[i] - sampling[i+1] )
#				
#		tbar[0] = t**3
#		tbar[1] = t**2
#		tbar[2] = t
#		tbar = tbar.reshape( (4,1) )
#		
#		w = (weight_function( sampling[i] ) + weight_function( sampling[i+1] ))/2
#		
#		W_i += dot( dt * w * tbar, tbar.T )
#	
#	return W_i

def precompute_W_i_with_weight_function_and_sampling( weight_function, sampling, ts, dts ):
	'''
	R = sum( T_i * P.T * M * partofR ) 
	compute integral of w * tbar * (M * tbar1)
			integral of w * tbar * (M * tbar2)
			integral of w * tbar * (M * tbar3)
			integral of w * tbar * (M * tbar4)
	'''
	### Asserts
	## Ensure our inputs are the same lengths:
	assert len( sampling ) == len( ts )
	assert len( sampling ) == len( dts ) + 1
	
	## Ensure our inputs are numpy.arrays:
	sampling = asarray( sampling )
	ts = asarray( ts )
	dts = asarray( dts )
	
	## Compute the integral.	
	R = zeros( ( 4, 4 ) )
	tbar = ones( 4 ).reshape( (4,1) )
	
	for i in range(len(dts)):
		t = (ts[i] + ts[i+1])/2
		dt = dts[i]
		
		tbar[0] = t*t*t
		tbar[1] = t*t
		tbar[2] = t
		
		w = (weight_function( i ) + weight_function( i+1 ))/2
		
		## M * tbar
		C_P = dot( M, tbar )
		
		R += dot( ((w*dt)*tbar), C_P.T )
		#R[:, 0] += asarray(w * tbar * C_P[0] *dt).reshape(-1)
		#R[:, 1] += asarray(w * tbar * C_P[1] *dt).reshape(-1)
		#R[:, 2] += asarray(w * tbar * C_P[2] *dt).reshape(-1)
		#R[:, 3] += asarray(w * tbar * C_P[3] *dt).reshape(-1)
	
	return R

def ts_and_dts_for_num_samples( a, b, num_samples ):
	'''
	Given two endpoints of integration 'a' and 'b',
	and positive integer 'num_samples' determining how many samples to use to compute
	the integral,
	returns two same-length arrays for numerical integration 'ts' and 'dts',
	where 'ts' contains the points at which to integrate and 'dts' contains
	the weight of the corresponding sample.
	'''
	dts = ( float(b-a)/num_samples ) * ones( len( num_samples ) )
	ts = [ a + ( ti + .5 ) * dt for ti in xrange( num_samples ) ]
	return ts, dts

def default_w_i( handle_positions, i, p ):	  
	'''
	Given an N-by-(2 or 3) numpy array 'vertices' of 2D or 3D vertices,
	an M-by-3 numpy array 'faces' of indices into 'vertices',
	an H-by-(2 or 3) numpy.array 'skeleton_handle_vertices' of 2D or 3D vertices,
	a numpy array 'skeleton_point_handles' of indices into 'skeleton_handle_vertices'
	which are the point handles,
	returns a N-by-H numpy.array of weights per vertex per handle.
	'''
	### Without error checking:
	'''
	inverse_distances = 1.0 / ( ( handle_positions - p )**2 ).sum( axis = 1 )
	return inverse_distances[i] / sum( inverse_distances )
	'''
	### With error checking:
	
	## Ensure our inputs are numpy.arrays:
	handle_positions = asarray( handle_positions )
	p = asarray( p ).reshape(3)
	
	## 'handle_positions' must be a num-handles-by-k array.
	assert len( handle_positions.shape ) == 2
	assert handle_positions.shape[1] == p.shape[0]
	## 'p' must be a k-vector.
	assert len( p.shape ) == 1
	
	## Compute the distance squared from each handle position.
	diff = handle_positions - p
	diff = diff**2
	diff = diff.sum( axis = 1 )
	
	## If any handle is within epsilon of p, we will get an almost divide-by-zero.
	## In this case, the weight should be 1 if the point is near that handle,
	## and 0 otherwise.
	eps = 1e-7
	## where() gives us the indices where the condition is true.
	wh = where( abs( diff ) < eps )
	## Actually, where() returns an array of arrays, one for each dimension of the input.
	## Since we have a one-dimensional input vector, where() should return an array of
	## length one containing the thing we actually care about, the indices into 'diff'.
	assert len( wh ) == 1
	## If where() returned any indices, then the weight function should be 1 if
	## 'p' is near handle 'i', and 0 otherwise.
	if len( wh[0] ) > 0: return 1.0 if i in wh[0] else 0.0
	
	## We want inverse distance:
	diff = 1.0 / diff
	
	return diff[i] / diff.sum()

def compute_error_metric( bbw_curve, spline_skin_curve, dts ):
	
	bbw_curve = asarray(bbw_curve).reshape(-1,2)
	spline_skin_curve = asarray(spline_skin_curve).reshape(-1,2)
	
	assert bbw_curve.shape == spline_skin_curve.shape
	
	diffs = ((spline_skin_curve - bbw_curve)**2 ).sum( axis = 1 )
	diffs = (diffs[:-1] + diffs[1:])/2
	
	diffs = dot(diffs, dts)
	
	bbw_lengths = [mag(bbw_curve[i]-bbw_curve[i+1]) for i in xrange( len( bbw_curve )-1 )]
	spline_lengths = [mag(spline_skin_curve[i]-spline_skin_curve[i+1]) for i in xrange( len( spline_skin_curve )-1 )]
	
	scale = sum( spline_lengths )
	
#	if diffs*scale > 100: debugger()
	
	return diffs*scale
