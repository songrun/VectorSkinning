from numpy import *

def precompute_W_i( handle_positions, i, P, M, a, b, num_samples = 100 ):
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
	
	### Compute the integral.
	W_i = zeros( ( 4,4 ) )
	tbar = ones( 4 )
	dt = (b-a)/num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
		t = (ti+.5) / num_samples
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
 		tbar = tbar.reshape( (4,1) )
		
		w = w_i( handle_positions, i, dot( P.T, dot( M.T, tbar ) ) )
		
		W_i += dot( dt * w * tbar, tbar.T )
	
	return W_i

def w_i( handle_positions, i, p ):
	'''
	Given a sequence of k-dimensional handle positions,
	an index 'i' specifying which handle,
	and a k-dimensional position 'p',
	returns	the i-th handle's weight function evaluated at p.
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

def precompute_inv_A( M, a, b, num_samples = 100 ):

	M = asarray( M )
	assert M.shape == (4,4)
	a = float(a)
	b = float(b)
	
	A = zeros( (4, 4) )
	tbar = ones( 4 )
	dt = (b-a)/num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
		t = (ti+.5) / num_samples
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )

		A = A + dot( dt * tbar, tbar.T )
	
	return linalg.inv( dot(M, A) )	
	
def control_points_after_split( P, S ):
	
	assert len( P.shape ) == 2
	assert P.shape[0] == 4
	P = asarray( P )
	
	S = asarray( S )
	
	if len( S ) < 2: 
		return P
	
	result = []
	for i, k in enumerate( S[1:-1] ):
		assert k > 0 and k < 1
		
		r1 = P[:-1]*(1.-k) + P[1:]*k
		r2 = r1[:-1]*(1.-k) + r1[1:]*k
		r3 = r2[:-1]*(1.-k) + r2[1:]*k
		
# 		print 'r1 ', r1, 'r2 ', r2, 'r3 ', r3, 'P ', P
		result = result + [P[0].tolist(), r1[0].tolist(), r2[0].tolist(), r3[0].tolist()]
# 		print 'result ', result
		P = array( [r3[-1], r2[-1], r1[-1], P[-1]] )
	
	result = result + P.tolist()
	
	return result
	