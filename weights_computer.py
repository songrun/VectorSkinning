from bezier_utility import *
from triangle import *
from bbw_wrapper import bbw
from itertools import izip as zip

# kEnableBBW = True

def uniquify_points_and_return_input_index_to_unique_index_map( pts, threshold = 0 ):
	'''
	Given a sequence of N points 'pts',
	and an optional 'threshold' indicating how many decimal places of accuracy (default: 0)
	returns two items:
	   a sequence of all the unique elements in 'pts'
	   and
	   a list of length N where the i-th item in the list tells you where
	   pts[i] can be found in the unique elements.
	'''
	
	from collections import OrderedDict
	unique_pts = OrderedDict()
	pts_map = []
	## Add rounded points to a dictionary and set the key to
	## ( the index into the ordered dictionary, the non-rounded point )
	for i, ( pt, rounded_pt ) in enumerate( zip( pts, map( tuple, asarray( pts ).round( threshold ) ) ) ):
		index = unique_pts.setdefault( rounded_pt, ( len( unique_pts ), pt ) )[0]
		pts_map.append( index )
	
	## Return the original resolution points.
	return [ tuple( pt ) for i, pt in unique_pts.itervalues() ], pts_map

def flatten_paths( all_pts ):
	'''
	Given a sequence of paths, each of which is a sequence of chains, each of which is a sequence of points,
	return the points as a flat sequence with a list of chain shapes (length of chain, number of points in each piece of chain).
	'''
	
	all_pts = asarray( all_pts )
	all_shapes = [ asarray( pts ).shape[:-1] for pts in all_pts ]
	
	all_pts = concatenate( [ concatenate( curve_pts ) for curve_pts in all_pts ] )
	return all_pts, all_shapes

def unflatten_data( flattened, all_shapes ):
	'''
	Given a flat sequence of data 'flattened' and
	an 'all_shapes' list of shapes as returned by flatten_paths(),
	returns the data from flattened after "unflattening"
	to have the same shape as 'all_shapes'.
	'''
	all_maps = []
	pos = 0
	for shape in all_shapes:
		maps = []
		for j in range( shape[0] ):
			maps.append( flattened[ pos : pos + shape[1] ] )
			pos += shape[1]
		all_maps.append( maps )
	
	return all_maps
	
	all_clean_pts = asarray( all_clean_pts )[:, :2]

def compute_all_weights( all_pts, skeleton_handle_vertices, boundary_index, which = None ):
	'''
	triangulate a region closed by a bunch of bezier curves if needed, and precompute the vertices at each sample point.
	
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	a sequence of M skeleton handle vertices, and
	the index into 'all_pts' of the boundary_curve (may be -1 for no boundary),
	a parameter 'which' specifying the style of weights ('bbw' or 'shepherd'),
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	## Debugging shepherd.
	# which = 'shepherd'
	
	if which is None: which = 'bbw'
	
	if 'bbw' == which:
		try:
			return compute_all_weights_bbw( all_pts, skeleton_handle_vertices, boundary_index )
		except bbw.BBWError as e:
			print 'BBW Computation failed:', e
			print 'Falling back to Shepherd weights.'
			which = 'shepherd'
	
	if 'shepherd' == which:
		return compute_all_weights_shepherd( all_pts, skeleton_handle_vertices )
	
	raise RuntimeError( "Unknown weight type" )

def compute_all_weights_shepherd( all_pts, skeleton_handle_vertices ):
	'''
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	and a sequence of M skeleton handle vertices
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	all_pts, all_shapes = flatten_paths( all_pts )
	
	print 'Removing duplicate points...'
	## Use 7 digits of accuracy. We're really only looking to remove actual duplicate
	## points.
	all_clean_pts, pts_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts, threshold = 7 )
	print '...finished.'
	
	all_maps = unflatten_data( pts_maps, all_shapes )
	
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	print 'Computing Shepherd weights...'
	all_weights = shepherd( all_clean_pts, skeleton_handle_vertices )
	print '...finished.'
	
	return all_clean_pts, all_weights, all_maps

def compute_all_weights_bbw( all_pts, skeleton_handle_vertices, boundary_index ):
	'''
	triangulate a region closed by a bunch of bezier curves if needed, and precompute the vertices at each sample point.
	
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	a sequence of M skeleton handle vertices, and
	the index into 'all_pts' of the boundary_curve (may be -1 for no boundary),
	a parameter 'which' specifying the style of weights ('bbw' or 'shepherd'),
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	if boundary_index < 0 or boundary_index >= len( all_pts ):
		raise RuntimeError( "compute_all_weights_bbw() got an invalid boundary curve" )
	
	all_pts, all_shapes = flatten_paths( all_pts )
	
	print 'Removing duplicate points...'
	all_clean_pts, pts_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts )
	print '...finished.'
	
	all_maps = unflatten_data( pts_maps, all_shapes )
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	
	## This will store a sequence of tuples ( edge_start_index, edge_end_index ).
	boundary_edges = []
	for curve in all_maps[ boundary_index ]:
		for vi in curve:
			## An edge in progress isn't a tuple, it's directly edge_start_index.
			if len( boundary_edges ) == 0:
				boundary_edges.append( vi )
			## Skip repeated points
			elif boundary_edges[-1] != vi:
				## Replace the edge in progress with a proper tuple.
				boundary_edges[-1] = ( boundary_edges[-1], vi )
				boundary_edges.append( vi )
	## We don't need to do anything to close the curve, because a closed curve
	## will have its last and first points overlapping.
	assert boundary_edges[-1] == boundary_edges[0][0]
	del boundary_edges[-1]
	
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

def shepherd( vs, skeleton_handle_vertices ):
	'''
	Given an N-by-(2 or 3) sequence 'vs' of 2D or 3D vertices and
	an H-by-(2 or 3) sequence 'skeleton_handle_vertices' of 2D or 3D vertices,
	returns a N-by-H numpy.array of weights per vertex per handle.
	'''
	
	vs = asarray( vs )
	
	assert len( vs.shape ) == 2
	assert vs.shape[1] == 2
	
	weights = ones( ( len( vs ), len( skeleton_handle_vertices ) ) )
	for vi, p in enumerate( vs ):
		for hi, p in enumerate( skeleton_handle_vertices ):
			weights[ vi, hi ] = shepherd_w_i( skeleton_handle_vertices, hi, p )
	
	return weights


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
	
	return R

# def ts_and_dts_for_num_samples( a, b, num_samples ):
#	'''
#	Given two endpoints of integration 'a' and 'b',
#	and positive integer 'num_samples' determining how many samples to use to compute
#	the integral,
#	returns two same-length arrays for numerical integration 'ts' and 'dts',
#	where 'ts' contains the points at which to integrate and 'dts' contains
#	the weight of the corresponding sample.
#	'''
#	dts = ( float(b-a)/num_samples ) * ones( len( num_samples ) )
#	ts = [ a + ( ti + .5 ) * dt for ti in xrange( num_samples ) ]
#	return ts, dts

def shepherd_w_i( handle_positions, i, p ):	  
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
	p = asarray( p ).reshape(-1)
	
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

def compute_error_metric( bbw_curve, skin_spline_curve, path_dts, lengths ):
	'''
	Total energy is the sum of each pair of points' square distance
	'''
	path_dts = asarray( path_dts )
	
	assert len( bbw_curve ) == len( skin_spline_curve ) == len( path_dts ) == len( lengths )

	energy = []
	for bbw_samplings, spline_samplings, segment_dts, length in zip( bbw_curve, skin_spline_curve, path_dts, lengths ):
		bbw_samplings = asarray( bbw_samplings ).reshape( -1, 2 )
		spline_samplings = asarray( spline_samplings ).reshape( -1, 2)
		dists = ( ( spline_samplings - bbw_samplings )**2 ).sum( axis = 1 )
		dists = (dists[:-1] + dists[1:])/2
	
		energy.append( dot( dists, segment_dts )*length )
	
	return energy
	
	
def compute_arc_length_error_metric( bbw_curve, skin_spline_curve, path_dss ):
	'''
	Total energy is the sum of each pair of points' square distance, compute with arc lengths
	'''
	path_dss = asarray( path_dss )
	
	assert len( bbw_curve ) == len( skin_spline_curve ) == len( path_dss )
	
	energy = []
	for bbw_samplings, spline_samplings, segment_dss in zip( bbw_curve, skin_spline_curve, path_dss ):
		bbw_samplings = asarray( bbw_samplings ).reshape( -1, 2 )
		spline_samplings = asarray( spline_samplings ).reshape( -1, 2)
		dists = ( ( spline_samplings - bbw_samplings )**2 ).sum( axis = 1 )
		dists = (dists[:-1] + dists[1:])/2
	
		energy.append( dot( dists, segment_dss ) )
	
	return energy	
