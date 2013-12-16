from bezier_utility import *
from triangle import *
import bbw_wrapper.bbw as bbw
from itertools import izip as zip


def triangulate_and_compute_weights(pts, skeleton_handle_vertices):
	'''
	trianglue a region closed by a bunch of bezier curves, precompute the vertices at each sample point.
	'''
	boundary_edges = [ ( i, (i+1) % len(pts) ) for i in xrange(len( pts )) ]
	
	pts = asarray( pts )[:, :2]
	
	if len( skeleton_handle_vertices ) > 0:
		skeleton_handle_vertices = asarray( skeleton_handle_vertices )[:, :2]
	skeleton_point_handles = list( range( len(skeleton_handle_vertices) ) )
	
	pts = concatenate( ( pts, skeleton_handle_vertices ), axis = 0 )
	vs, faces = triangles_for_points( pts, boundary_edges )
	vs = asarray(vs)[:, :2].tolist()
	
	faces = asarray(faces).tolist()
	all_weights = bbw.bbw(vs, faces, skeleton_handle_vertices, skeleton_point_handles)
	
	return vs, faces, all_weights

def precompute_W_i_bbw( vs, weights, i, sampling, ts, dts = None ):
	'''
	Given an N-by-k numpy.array 'vs' of all points represented in 'weights',
	an N-by-num-handles numpy.array 'weights' of all the weights for each sample vertex,
	an index 'i' specifying which handle,
	an M-by-k numpy.array 'sampling' containing sampled positions,
	a length-M numpy.array of t values corresponding to each sample in 'sampling',
	an optional length-M numpy.array of dt values corresponding to each sample in 'sampling',
	returns W, a 4-by-4 numpy.array defined as:
		\int weights_i( sample ) \overbar{t}^T \overbar{t}^T dt
	where sample is drawn from 'sampling', t is drawn from 'ts', and dt is drawn from 'dts'.
	
	If 'dts' is not given, it defaults to 1/len(sampling).
	'''
	
	if dts is None: dts = ones( len( sampling ) ) * (1./len(sampling))
	
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

	result = zeros((2,4,4 ))

	## Vertices and sampling must have the same dimension for each point.
	sampling = sampling[:,:-1]
	assert vs.shape[1] == sampling.shape[1]
	
	## The index 'i' must be valid.
	assert i >= 0 and i < weights.shape[1]
	
	def weight_function( p ):
		## Find the closest vertex in 'vs' to 'p'
		vi = argmin( ( ( vs - p )**2 ).sum( axis = 1 ) )
		assert allclose( vs[vi], p, 1e-5 )
		return weights[ vi, i ]
	
	result[0] = precompute_W_i_with_weight_function_and_sampling( weight_function, 
				sampling, ts, dts )
	result[1] = precompute_partOfR_with_weight_function_and_sampling( weight_function, 
				sampling, ts, dts )
				
	return result

def precompute_W_i_with_weight_function_and_sampling( weight_function, sampling, ts, dts ):
	'''
	Given a function 'weight_function' that takes a point and returns its weight,
	a N-by-k numpy.array 'sampling' containing the positions of the control points as the rows,
	corresponding t values 'ts' for each point in 'sampling',
	and an optional corresponding 'dt' for each point in sampling (default is 1/len(sampling)),
	returns W, a 4-by-4 numpy.array defined as:
		\int_i weight_function( sample ) \overbar{t}^T \overbar{t}^T dt
	where sample, t, and dt are drawn from the corresponding input arrays.
	
	The optional parameter 'num_samples' determines how many samples to use to compute
	the integral.
	'''
	
	### Asserts
	## Ensure our inputs are the same lengths:
	assert len( sampling ) == len( ts )
	assert len( sampling ) == len( dts )
	
	## Ensure our inputs are numpy.arrays:
	sampling = asarray( sampling )
	ts = asarray( ts )
	dts = asarray( dts )
	
	## sampling must be N-by-k.
	assert len( sampling.shape ) == 2
	assert len( ts.shape ) == 1
	assert len( dts.shape ) == 1
	
	### Compute the integral.
	W_i = zeros( ( 4,4 ) )
	tbar = ones( 4 )
	
	for sample, t, dt in zip( sampling, ts, dts ):
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )
		
		w = weight_function( sample )
		
		W_i += dot( dt * w * tbar, tbar.T )
	
	return W_i

def precompute_partOfR_with_weight_function_and_sampling( weight_function, sampling, ts, dts ):
	'''
	compute integral of w * tbar * (M * tbar11)
			integral of w * tbar * (M * tbar21)
			integral of w * tbar * (M * tbar31)
			integral of w * tbar * (M * tbar41)
	'''
	## Compute the integral.	
	R = zeros( ( 4, 4 ) )
	tbar = ones( 4 )
	
	for sample, t, dt in zip( sampling, ts, dts ):
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )
		
		w = weight_function( sample )
		
		## M * tbar
		C_P = dot( M, tbar )
		
		R[:, 0] += asarray(w * tbar * C_P[0] *dt).reshape(-1) 
		R[:, 1] += asarray(w * tbar * C_P[1] *dt).reshape(-1)
		R[:, 2] += asarray(w * tbar * C_P[2] *dt).reshape(-1)
		R[:, 3] += asarray(w * tbar * C_P[3] *dt).reshape(-1)	
	
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

def precompute_W_i_default( handle_positions, i, P, M, a, b, num_samples = 100 ):
	'''
	Given a sequence of k-dimensional handle positions,
	an index 'i' specifying which handle,
	a 4-by-k numpy.array P containing the positions of the control points as the rows,
	a 4-by-4 numpy.array M containing the Bezier curve weights,
	and the interval to integrate from 'a' to 'b',
	returns W, a 4-by-4 numpy.array defined as:
		\int_a^b w_i( P^T M^T \overbar{t}^T ) \overbar{t}^T \overbar{t}^T dt
	
	The optional parameter 'num_samples' determines how many samples to use to compute
	the integral.
	'''
	
	### Asserts
	## Ensure our inputs are numpy.arrays:
	handle_positions = asarray( handle_positions )
	P = asarray( P )
	M = asarray( M )
	
	## We must have at least one handle position.
	assert len( handle_positions ) > 0
	## The index 'i' must be valid.
	assert i >= 0 and i < len( handle_positions )
	## 'handle_positions' must be a num-handles-by-k array.
	assert len( handle_positions.shape ) == 2
	assert handle_positions.shape[1] == P.shape[1]
	
	## P must be 4-by-k.
	assert len( P.shape ) == 2
	assert P.shape[0] == 4
	
	## M must be 4-by-4.
	assert M.shape == (4,4)
	
	## Use a default number of samples of 100.
	assert num_samples > 0
	
	ts_and_dts_for_num_samples( a, b, num_samples )
	samplings = []
	for t, dt in zip( ts, dts ):
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )
		
		samplings.append( dot( P.T, dot( M.T, tbar ) ) )
	
	def weight_function( p ):
		return default_w_i( handle_positions, i, p )
	
	return precompute_W_i_with_weight_function_and_sampling( weight_function, sampling, ts, dts )

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

def compute_error_metric( transforms, handle_positions, P_primes, P, M, num_samples = 100, dim = 3 ):
	
	P_primes = asarray( P_primes )
	P = asarray( P )
	M = asarray( M )
	
	assert P_primes.shape == (4, dim)
	assert P.shape == (4, dim)
	assert M.shape == (4, 4)
	
	Error = 0.0
	tbar = ones( 4 )
	dt = 1./num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
		t = ( ti + .5 ) * dt
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )
		
		tmp = 0
		for i in range( len( transforms ) ):
			w = w_i( handle_positions, i, dot( P.T, dot( M.T, tbar ) ) )
			tmp += w * dot( transforms[i].reshape(3,3),	 dot( dot(P.T, M), tbar ) )
			
		diff = dot( dot(P_primes.T, M), tbar ) - tmp 
		Error += dot( diff.T, diff ) * dt
		
	return Error
