from bezier_utility import *
from triangle import *
from bbw_wrapper import bbw
from itertools import izip as zip
from tictoc import tic, toc

# kEnableBBW = True
kBarycentricProjection = False

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
		## For fancier schemes:
		# index = unique_pts.setdefault( rounded_pt, ( len( unique_pts ), [] ) )[0]
		# unique_pts[ rounded_pt ][1].append( pt )
		pts_map.append( index )
	
	## Return the original resolution points.
	## The average of all points that round:
	# return [ tuple( average( pt, axis = 0 ) ) for i, pt in unique_pts.itervalues() ], pts_map
	## The closest point to the rounded point:
	# return [ tuple( pt[ abs( asarray( pt ).round( threshold ) - asarray( pt ) ).sum(axis=1).argmin() ] ) for i, pt in unique_pts.itervalues() ], pts_map
	## Simplest, the first rounded point:
	return [ tuple( pt ) for i, pt in unique_pts.itervalues() ], pts_map

def barycentric_projection( vs, faces, boundary_edges, weights, pts ):
	'''
	Given a sequence 'vertices' and 'faces' representing a 2D triangle mesh,
	a sequence of pairs of indices into 'vertices' corresponding to the
	boundary edges of the mesh,
	a sequence of (not necessarily scalar-valued) values 'weights', one for each vertex in 'vs',
	and a sequence of points 'pts'
	returns
		a sequence of uniqified points from 'pts',
		a corresponding interpolated weight for the uniqified points,
		and map from each element of 'pts' to the uniqified sequence.
	
	
	tested:
	vs = [ (0,0), (1,0), (1,1), (0,1) ]
	faces = [ ( 0,1,2 ), ( 2, 3, 0 ) ]
	boundary_edges = [ ( 0,1 ), ( 1,2 ), ( 2,3 ), ( 3, 0 ) ]
	weights = asarray([ [ 1,0,0,0 ], [ 0,1,0,0 ], [ 0,0,1,0 ], [ 0,0,0,1 ] ])
	pts = [ (0,0), (1,0), (1,1), (0,1), (.2,.1), (.9,.8), (.8,.9), ( -1, -1 ), ( -1, 1 ) ]
	unique_pts, unique_weights, pts_map = barycentric_projection( vs, faces, boundary_edges, weights, pts )
	out: [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.20000000000000001, 0.10000000000000001), (0.90000000000000002, 0.80000000000000004), (0.80000000000000004, 0.90000000000000002), (-1.0, -1.0), (-1.0, 1.0)]
	out: array([[ 1. ,	0. ,  0. ,	0. ],
	   [ 0. ,  1. ,	 0. ,  0. ],
	   [ 0. ,  0. ,	 1. ,  0. ],
	   [ 0. ,  0. ,	 0. ,  1. ],
	   [ 0.8,  0.1,	 0.1,  0. ],
	   [ 0.1,  0.1,	 0.8,  0. ],
	   [ 0.1,  0. ,	 0.8,  0.1],
	   [ 1. ,  0. ,	 0. ,  0. ],
	   [ 0. ,  0. ,	 0. ,  1. ]])
	out: [0, 1, 2, 3, 4, 5, 6, 7, 8]
	'''
	
	tic( 'Barycentric projection...' )
	
	from raytri import raytri
	
	pts = asarray( pts )
	
	## TODO Q: Should we uniquify points even though we don't have to?
	kRemoveDuplicates = True
	## A1: Yes, because our point2d_in_mesh2d_barycentric() function is slow.
	if kRemoveDuplicates:
		tic( 'Removing duplicate points...' )
		## Use 7 digits of accuracy. We're really only looking to remove actual duplicate
		## points.
		unique_pts, unique_map = uniquify_points_and_return_input_index_to_unique_index_map( pts, threshold = 7 )
		unique_pts = asarray( unique_pts )
		toc()
	## A2: No, because we don't have to.
	else:
		unique_pts = pts
		unique_map = range(len( pts ))
	
	
	edges = zeros( ( len( boundary_edges ), 2, len( vs[0] ) ) )
	for bi, ( e0, e1 ) in enumerate( boundary_edges ):
		edges[ bi ] = vs[ e0 ], vs[ e1 ]
	
	## Using vertex positions as weights should lead to the
	## identity transformation. (See comment d987dsa98d7h below.)
	# weights = array( vs )
	
	misses = 0
	misses_total_distance = 0.
	misses_max_distance = -31337.
	unique_weights = zeros( ( len( unique_pts ), len( weights[0] ) ) )
	for pi, pt in enumerate( unique_pts ):
		bary = raytri.point2d_in_mesh2d_barycentric( pt, vs, faces )
		## Did we hit the mesh?
		if bary is not None:
			fi, ( b0, b1, b2 ) = bary
			#assert abs( b0 + b1 + b2 - 1 ) < 1e-5
			#assert b0 > -1e-5
			#assert b1 > -1e-5
			#assert b2 > -1e-5
			#assert b0 < 1+1e-5
			#assert b1 < 1+1e-5
			#assert b2 < 1+1e-5
			unique_weights[pi] = b0*weights[ faces[ fi ][0] ] + b1*weights[ faces[ fi ][1] ] + b2*weights[ faces[ fi ][2] ]
		else:
			#print 'pi outside:', pi
			dist, ei, t = raytri.closest_distsqr_and_edge_index_and_t_on_edges_to_point( edges, pt )
			#assert t > -1e-5
			#assert t < 1+1e-5
			dist = sqrt( dist )
			misses += 1
			misses_total_distance += dist
			misses_max_distance = max( misses_max_distance, dist )
			unique_weights[pi] = (1-t)*weights[ boundary_edges[ ei ][0] ] + t*weights[ boundary_edges[ ei ][1] ]
	
	## And indeed it does come out nearly identical. (See comment d987dsa98d7h above.)
	# assert ( unique_weights - pts ).allclose()
	
	#assert unique_weights.min() > -1e-4
	#assert unique_weights.max() < 1 + 1e-4
	## Clip the weights?
	# unique_weights = unique_weights.clip( 0, 1 )
	
	## Re-normalize the weights?
	unique_weights *= 1./unique_weights.sum( axis = 1 )[...,newaxis]
	
	if misses == 0:
		print 'Barycentric projection: No one missed the mesh.'
	else:
		print 'Barycentric projection:', misses, 'points missed the mesh. Average distance was', misses_total_distance/misses, ' and maximum distance was', misses_max_distance
	
	toc()
	
	return unique_pts, unique_weights, unique_map

def flatten_paths( all_pts ):
	'''
	Given a sequence of paths, each of which is a sequence of chains, each of which is a sequence of points,
	return the points as a flat sequence with a list of chain shapes (length of chain, number of points in each piece of chain).
	'''
	
	all_pts = asarray( all_pts )
	all_shapes = [ asarray( pts ).shape[:-1] for pts in all_pts ]
	
	all_new_pts = concatenate( [ concatenate( path_pts ) for path_pts in all_pts ] )

	return all_new_pts, all_shapes

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

def compute_all_weights( all_pts, skeleton_handle_vertices, boundary_index, which = None ):
	'''
	triangulate a region closed by a bunch of bezier curves if needed, and precompute the vertices at each sample point.
	
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	a sequence of M skeleton handle vertices, and
	the index into 'all_pts' of the boundary_curve (may be -1 for no boundary),
	a parameter 'which' specifying the style of weights ('bbw' or 'shepard' or 'mvc'),
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	## To try shepard no matter what:
	# which = 'shepard'
	
	if which is None: which = 'bbw'
	
	if 'bbw' == which:
		try:
			return compute_all_weights_bbw( all_pts, skeleton_handle_vertices, boundary_index )
		except bbw.BBWError as e:
			print 'BBW Computation failed:', e
			print 'Falling back to Shepard weights.'
			which = 'shepard'
	
	if 'harmonic' == which:
		try:
			return compute_all_weights_harmonic( all_pts, skeleton_handle_vertices )
		except bbw.BBWError as e:
			print 'Harmonic Computation failed:', e
			print 'Falling back to Shepard weights.'
			which = 'shepard'
	
	if 'mvc' == which:
		return compute_all_weights_mvc( all_pts, skeleton_handle_vertices )
	
	if 'shepard' == which:
		return compute_all_weights_shepard( all_pts, skeleton_handle_vertices )
	
	raise RuntimeError( "Unknown weight type" )

def compute_all_weights_shepard( all_pts, skeleton_handle_vertices ):
	'''
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	and a sequence of M skeleton handle vertices
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	all_pts, all_shapes = flatten_paths( all_pts )
	
	#tic( 'Removing duplicate points...' )
	## Use 7 digits of accuracy. We're really only looking to remove actual duplicate
	## points.
	#all_clean_pts, pts_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts, threshold = 7 )
	#toc()
	## UPDATE: There's no need to remove duplicates.
	all_clean_pts = asarray( all_pts )
	pts_maps = range( len( all_clean_pts ) )
	
	all_maps = unflatten_data( pts_maps, all_shapes )
	
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	tic( 'Computing Shepard weights...' )
	all_weights = shepard( all_clean_pts, skeleton_handle_vertices )
	toc()
	
	return all_clean_pts, all_weights, all_maps

def compute_all_weights_mvc( all_pts, cage_loop ):
	'''
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	and a sequence of M cage loop vertices 'cage_loop'
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	all_pts, all_shapes = flatten_paths( all_pts )
	all_maps = unflatten_data( range(len( all_pts )), all_shapes )
	
	tic( 'Computing Mean Value Coordinate weights...' )
	all_weights = bbw.mvc( all_pts, cage_loop )
	toc()
	
	return all_pts, all_weights, all_maps

def compute_all_weights_bbw( all_pts, skeleton_handle_vertices, boundary_index, customized = False ):
	'''
	triangulate a region closed by a bunch of bezier curves if needed, and precompute the vertices at each sample point.
	
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	a sequence of M skeleton handle vertices, and
	the index into 'all_pts' of the boundary_curve (may be -1 for no boundary),
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	if boundary_index < 0 or boundary_index >= len( all_pts ):
		raise RuntimeError( "compute_all_weights_bbw() got an invalid boundary curve" )
	
	all_pts, all_shapes = flatten_paths( all_pts )
	
	tic( 'Removing duplicate points...' )
	all_clean_pts, pts_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts, threshold = 0 )
	toc()
	
	all_maps = unflatten_data( pts_maps, all_shapes )
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	
	## This will store a sequence of tuples ( edge_start_index, edge_end_index ).
	## UPDATE: We need to make sure that this boundary loop stays manifold.
	##		   That means: no vertex index should be the start index more than once,
	##		   and no vertex index should be the end index more than once.
	boundary_edges = []
	for curve in all_maps[ boundary_index ]:
		for vi in curve:
			## An edge in progress isn't a tuple, it's directly edge_start_index.
			if len( boundary_edges ) == 0:
				boundary_edges.append( vi )
			## Skip repeated points
			elif boundary_edges[-1] == vi:
				## This happens a lot.
				# print 'Skipping a collapsed boundary edge'
				pass
			## UPDATE: And skip a point that would make us fold back on ourselves!
			##		   It's possible due to rounding that the two points on the
			##		   bottom of a ^ sticking out of the mesh collapses, leading
			##		   us with a | shape and a duplicate (undirected) edge.
			##		   It's easy to check for that here.
			##		   If it comes up again, though, just wait until the
			##		   end of the loop and filter boundary_edges like so:
			##			   boundary_edges = [ tuple( edge ) for edge in set([ frozenset( edge ) for edge in boundary_edges ]) ]
			##		   For now, just check for immediately collapsed edges.
			##		   NOTE: This came up with the alligator.
			## UPDATE 2: Filtering at the end still might not work in weird cases.
			##			 If it comes up again, we'll need to make sure that
			##			 the boundary loop stays manifold, which means:
			##			 no vertex index should be the start index more than once,
			##			 and no vertex index should be the end index more than once.
			elif len( boundary_edges ) > 1 and boundary_edges[-2][0] == vi:
				print 'Skipping a boundary foldback'
				pass
			else:
				## Replace the edge-in-progress with a proper tuple.
				boundary_edges[-1] = ( boundary_edges[-1], vi )
				boundary_edges.append( vi )
	## We don't need to do anything to close the curve, because a closed curve
	## will have its last and first points overlapping.
	assert boundary_edges[-1] == boundary_edges[0][0]
	del boundary_edges[-1]
	
	## The list of handles.
	if len( skeleton_handle_vertices ) > 0:
		skeleton_handle_vertices = asarray( skeleton_handle_vertices )[:, :2]
	skeleton_point_handles = list( range( len(skeleton_handle_vertices) ) )
	
	registered_pts = concatenate( ( all_clean_pts, skeleton_handle_vertices ), axis = 0 )
	tic( 'Computing triangulation...' )
	vs, faces = triangles_for_points( registered_pts, boundary_edges )
	toc()
	
	vs = asarray(vs)[:, :2] 
	faces = asarray(faces)
	
	tic( 'Computing BBW...' )
	all_weights = bbw.bbw(vs, faces, skeleton_handle_vertices, skeleton_point_handles)
	toc()
	
	if kBarycentricProjection:
		if __debug__: old_weights = asarray([ all_weights[i] for i in pts_maps ])
		
		vs, all_weights, pts_maps = barycentric_projection( vs, faces, boundary_edges, all_weights, all_pts )
		all_maps = unflatten_data( pts_maps, all_shapes )
		
		if __debug__:
			new_weights = asarray([ all_weights[i] for i in pts_maps ])
			total_weight_change = abs(old_weights-new_weights).sum()
			print 'Barycentric projection led to an average change in weights of', total_weight_change/prod( new_weights.shape ), 'and a total change of', total_weight_change
	
	if customized == False:
		return vs, all_weights, all_maps
	## for the test of naive approaches.
	else:
		return vs, faces, boundary_edges, all_weights, all_maps

def compute_all_weights_harmonic( all_pts, skeleton_handle_vertices, customized = False ):
	'''
	triangulate a region closed the handles as a cage, and precompute the vertices at each sample point.
	
	Given a sequence of sequences of sequences of points 'all_pts' (paths of chains of sampled bezier curves),
	a sequence of M skeleton handle vertices, and
	the index into 'all_pts' of the boundary_curve (may be -1 for no boundary),
	returns
		a sequence of vertices,
		a M-dimensional weight for each vertex,
		and a sequence of sequences mapping the index of a point in 'all_pts' to a vertex index.
	'''
	
	all_pts, all_shapes = flatten_paths( all_pts )
	
	tic( 'Removing duplicate points...' )
	all_clean_pts, pts_maps = uniquify_points_and_return_input_index_to_unique_index_map( all_pts, threshold = 0 )
	toc()
	
	all_maps = unflatten_data( pts_maps, all_shapes )
	all_clean_pts = asarray( all_clean_pts )[:, :2]
	
	## The list of handles.
	if len( skeleton_handle_vertices ) > 0:
		skeleton_handle_vertices = asarray( skeleton_handle_vertices )[:, :2]
	skeleton_point_handles = list( range( len(skeleton_handle_vertices) ) )
	
	registered_pts = concatenate( ( all_clean_pts, skeleton_handle_vertices ), axis = 0 )
	## The boundary edges are the handle vertices as a loop.
	off = len(all_clean_pts)
	boundary_edges = [ ( off + i, off + ( (i+1) % len( skeleton_handle_vertices ) ) ) for i in xrange(len( skeleton_handle_vertices )) ]
	
	tic( 'Computing triangulation...' )
	vs, faces = triangles_for_points( registered_pts, boundary_edges )
	toc()
	
	vs = asarray(vs)[:, :2] 
	faces = asarray(faces)
	
	tic( 'Computing Harmonic Coordinates...' )
	all_weights = bbw.harmonic( vs, faces, [ i for i,j in boundary_edges ], 1 )
	toc()
	
	if kBarycentricProjection:
		if __debug__: old_weights = asarray([ all_weights[i] for i in pts_maps ])
		
		vs, all_weights, pts_maps = barycentric_projection( vs, faces, boundary_edges, all_weights, all_pts )
		all_maps = unflatten_data( pts_maps, all_shapes )
		
		if __debug__:
			new_weights = asarray([ all_weights[i] for i in pts_maps ])
			total_weight_change = abs(old_weights-new_weights).sum()
			print 'Barycentric projection led to an average change in weights of', total_weight_change/prod( new_weights.shape ), 'and a total change of', total_weight_change
	
	if customized == False:
		return vs, all_weights, all_maps
	## for the test of naive approaches.
	else:
		return vs, faces, boundary_edges, all_weights, all_maps

def shepard( vs, skeleton_handle_vertices ):
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
		for hi in range( len( skeleton_handle_vertices ) ):
			weights[ vi, hi ] = shepard_w_i( skeleton_handle_vertices, hi, p )
	
	# print 'Shepard fast - slow:', abs( weights - shepard_fast( vs, skeleton_handle_vertices ) ).max()
	
	return weights

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

shepard = shepard_fast

def precompute_W_i( vs, weights, i, sampling_index2vs_index, sampling, ts, dts ):
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
	
	#if dts is None: dts = ones( len( sampling )-1 ) * (1./(len(sampling)-1) )
#	dts = ones( len( sampling )-1 ) * (1./(len(sampling)-1) ) * sum( dts )
	
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

def precompute_W_i_fast( vs, weights, i, sampling_index2vs_index, sampling, ts, dts ):
	sampling_weights = weights[ sampling_index2vs_index, i ]
	
	wdt = .5*(sampling_weights[:-1] + sampling_weights[1:])*dts
	midts = .5*(ts[:-1] + ts[1:])
	tbars = array([midts**3, midts**2, midts, ones(midts.shape).squeeze()]).T
	
	C_P = dot( asarray( M ), tbars.T )
	R = dot( ( wdt[:,newaxis]*tbars ).T, C_P.T )
	return R

precompute_W_i = precompute_W_i_fast

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
	
	#dtts = zeros(len(ts))
	#dtts[1:] += dts*.5
	#dtts[:-1] += dts*.5
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

def shepard_w_i( handle_positions, i, p ):	  
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

def compute_error_metric( target_path, deformed_path, path_dts, lengths ):
	'''
	Total energy is the sum of each pair of points' square distance
	'''
	path_dts = asarray( path_dts )
	
	assert len( target_path ) == len( deformed_path ) == len( path_dts ) == len( lengths )
	target_path, deformed_path = asarray ( target_path ), asarray( deformed_path )
	assert target_path.shape == deformed_path.shape

	energy = []
	for target_curve, deformed_curve, segment_dts, length in zip( target_path, deformed_path, path_dts, lengths ):
		target_curve = asarray( target_curve ).reshape( -1, 2 )
		deformed_curve = asarray( deformed_curve ).reshape( -1, 2)
		dists = ( ( deformed_curve - target_curve )**2 ).sum( axis = 1 )
		dists = (dists[:-1] + dists[1:])/2
	
		energy.append( dot( dists, segment_dts )*length )
	
	return energy
	
def compute_maximum_distances( target_path, spline_path ):
	'''
	Find the approximate largest distance between a spline curve and its target curve.
	
	untested
	'''
	assert len( target_path ) == len( spline_path )
	target_path, spline_path = asarray ( target_path ), asarray( spline_path )
	assert target_path.shape == spline_path.shape
	
	distances = []
	for target_curve, spline_curve in zip( target_path, spline_path ):
		## allDistSqrs[i][j] is the distance squared from target_curve[i] to spline_curve[j].
		allDistSqrs = ( (spline_curve[newaxis,...] - target_curve[:,newaxis,:])**2 ).sum(-1)
		## Hausdorff distance is the longest shortest distance from either to either.
		dist2 = max( allDistSqrs.min(0).max(), allDistSqrs.min(1).max() )
		indices = where( allDistSqrs == dist2 )
		assert len( indices ) > 0
		
		target_index = indices[0][0]
		spline_index = indices[1][0]
		dist = sqrt( dist2 )
		
		distances.append({
			'spline_pos': spline_curve[ spline_index ].tolist(),
			'target_pos': target_curve[ target_index ].tolist(),
			'maximum_distance': dist
			})
	
	return distances		
