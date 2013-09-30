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

# 	dt = (b-a)/num_samples
# 	#for t in linspace( a, b, num_samples ):
# 	for ti in xrange( num_samples ):
# 		t = a + ( ti + .5 ) * dt
# 		
# 		tbar[0] = t**3
# 		tbar[1] = t**2
# 		tbar[2] = t
#  		tbar = tbar.reshape( (4,1) )
# 		
# 		w = w_i( handle_positions, i, dot( P.T, dot( M.T, tbar ) ) )
# 		
# 		W_i += dot( dt * w * tbar, tbar.T )
	
	dt = (b-a)/num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
		t = a + ( ti + .5 ) * dt
		
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

def precompute_A( M, a, b, num_samples = 100 ):

	M = asarray( M )
	assert M.shape == (4,4)
	a = float(a)
	b = float(b)
	'''
	t_prime = (1-t)*a + t*b
	'''
	a_prime = (1-a)*a + a*b
	b_prime = (1-b)*a + b*b
	
	A = zeros( (4, 4) )
	tbar = ones( 4 )
	dt = (1.0 - 0.)/num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
	
		t = ( ti + 0.5 ) * dt
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		tbar = tbar.reshape( (4,1) )

		A = A + dot( dt * tbar, tbar.T )
	
	A = A*(b-a)
	
	return dot(M, A)

def precompute_partOfR( handle_positions, i, P, M ,a, b, num_samples = 100 ):
	### Compute the integral.
	R = zeros( ( 4, 4 ) )
	tbar = ones( 4 )

	dt = (b-a)/num_samples
	#for t in linspace( a, b, num_samples ):
	for ti in xrange( num_samples ):
		t = a + ( ti + .5 ) * dt
			
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
 		tbar = tbar.reshape( (4,1) )
		
		w = w_i( handle_positions, i, dot( P.T, dot( M.T, tbar ) ) )
		
		C_P = dot( M, tbar )
		
		R[:, 0] += asarray(w * tbar * C_P[0] *dt).reshape(-1) 
		R[:, 1] += asarray(w * tbar * C_P[1] *dt).reshape(-1)
		R[:, 2] += asarray(w * tbar * C_P[2] *dt).reshape(-1)
		R[:, 3] += asarray(w * tbar * C_P[3] *dt).reshape(-1)	
	
	return R

## Compute error computed by numerical approach
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
			tmp += w * dot( transforms[i].reshape(3,3),  dot( dot(P.T, M), tbar ) )
			
		diff = dot( dot(P_primes.T, M), tbar ) - tmp 
		Error += dot( diff.T, diff ) * dt
		
	return Error
