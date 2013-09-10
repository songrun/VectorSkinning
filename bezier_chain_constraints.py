from numpy import *
from weight_inverse_distance import *

M = matrix('-1. 3. -3. 1.; 3. -6. 3. 0.; -3. 3. 0. 0.; 1. 0. 0. 0.')

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
		
		result = result + [P[0].tolist(), r1[0].tolist(), r2[0].tolist(), r3[0].tolist()]

		P = array( [r3[-1], r2[-1], r1[-1], P[-1]] )
	
	result = result + P.tolist()
	
	return result

def compute_control_points_chain_with_constraint( partition, original_controls, handles, transforms, constraint_level ):

	result = []
	
	Cset =  control_points_after_split( original_controls, partition ) 
	Cset = Cset + original_controls.tolist()
	Cset = asarray( Cset )
	
	if constraint_level == 0:
		
 		result = compute_control_points_chain_without_constraint( Cset, handles, transforms )
#  		result = compute_control_points_chain_without_constraint_2( original_controls, partition, handles, transforms )	
	elif constraint_level == 1:
		
		result = compute_control_points_chain_with_continuity( Cset, handles, transforms )
	
	elif constraint_level == 2:
	
		result = compute_control_points_chain_with_derivative_continuity( Cset, handles, transforms )
	
	elif constraint_level == 3:
		
		result = compute_control_points_chain_with_second_derivative_continuity( Cset, handles, transforms )

	return result
	
def compute_control_points_chain_without_constraint( controls, handles, transforms ): 

	result = []
	
	inv_A = linalg.inv( precompute_A( M, 0., 1., 50 ) )					
	for k in range( len( controls ) /4  ):
		
		P_prime = zeros(12).reshape(3,4)
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(3,3) )
			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )		
			
			#for num_samples in xrange(1,501):
			#	writer.writerow( [ num_samples, i ] + list( precompute_W_i( handles, i, cps, M, 0., 1., num_samples ).flat ) )
			
			P_prime = P_prime + T_i * (controls[k*4:k*4+4].T) * M * mat(W_i) * inv_A	

		result.append( asarray(P_prime.T) )
	
	return result

'''
	has bugs
'''
def compute_control_points_chain_without_constraint_2( controls, partition, handles, transforms ):

	result = []
	
	for k in range( len( partition ) - 1  ):
		
		a, b = partition[k], partition[k+1]
		a, b = 0.0, 0.5
		inv_A = linalg.inv( precompute_A( M, a, b, 50 )	)
		P_prime = zeros(12).reshape(3,4)
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(3,3) )
			W_i = precompute_W_i( handles, i, controls, M, a, b, 50 )
			
			P_prime = P_prime + T_i * mat(controls.T) * M * mat(W_i) * inv_A	

		result.append( P_prime.T )
	
	return result

'''
	Compute the control points for each splited curve with endpoints connected.
	Boundary Conditions are as follows:
		lambda1 * ( P41' - Q11' ) = 0
		lambda2 * ( P42' - Q12' ) = 0
		lambda3 * ( P43' - Q13' ) = 0
'''	
def compute_control_points_chain_with_continuity( controls, handles, transforms ):

	temps = []
	result = []
	
	const_k = 3
	num = len( controls ) /4 -1 # num is the number of splited curves
	
	A = precompute_A( M, 0., 1., 50 )
	inv_M = linalg.inv( M )
	
	for k in range( num ):

		C = zeros(const_k*4).reshape(const_k,4)
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(const_k,const_k) )
			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )	
	
			C = C + T_i * (controls[k*4:k*4+4].T) * M * mat(W_i) * M
	
		temps.append( asarray(C.T) )
	
	Right = []

	for R in temps:
		Right = Right + (R.T.reshape(const_k*4)).tolist()
	Right = array( Right + zeros( const_k * (num-1) ).tolist() )
	
	dim = len( Right )
	AA = M.T * A.T
	assert AA.shape == (4,4)
	
	Left =  zeros( (dim, dim) )		
	for i in range( num * const_k ):
		Left[ i*4:i*4+4, i*4:i*4+4 ] = AA[:,:]
	
	R = zeros( ( 8*const_k, const_k ) )
	for i in range( const_k ):
		R[i*4+3, i] = 1
		R[4*const_k + i*4, i] = -1
		
	for i in range( num-1 ): 
		Left[ 4*const_k*i:(4*i+8)*const_k, num*const_k*4+const_k*i:num*const_k*4+const_k*(i+1) ] = R
		Left[ num*const_k*4+const_k*i:num*const_k*4+const_k*(i+1), 4*const_k*i:(4*i+8)*const_k ] = R.T
		
	X = linalg.solve(Left, Right)	
	X = array( X[:num*const_k*4] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*const_k:(i+1)*const_k ] )		
	
	return result

'''
	Compute the control points for each splited curve with endpoints connected and the derivative at the endpoints equaled.
	Boundary Conditions are as follows:
		lambda1 * ( P41' - Q11' ) = 0
		lambda2 * ( P42' - Q12' ) = 0
		lambda3 * ( P43' - Q13' ) = 0
		lambda4 * ( P41' - P31' + Q11' - Q21') = 0
		lambda5 * ( P42' - P32' + Q12' - Q22') = 0
		lambda6 * ( P43' - P33' + Q13' - Q23') = 0
'''		
def compute_control_points_chain_with_derivative_continuity( controls, handles, transforms ):

	temps = []
	result = []
	
	const_k = 3
	num = len( controls ) /4 -1 # num is the number of splited curves
	
	A = precompute_A( M, 0., 1., 50 )
	inv_M = linalg.inv( M )
	
	for k in range( num ):

		C = zeros(const_k*4).reshape(const_k,4)
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(const_k,const_k) )
			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )	
	
			C = C + T_i * (controls[k*4:k*4+4].T) * M * mat(W_i) * M
	
		temps.append( asarray(C.T) )
	
	Right = []

	for R in temps:
		Right = Right + (R.T.reshape(const_k*4)).tolist()
	Right = array( Right + zeros( 2 * const_k * (num-1) ).tolist() )
	
	dim = len( Right )
	AA = M.T * A.T
	assert AA.shape == (4,4)
	
	Left =  zeros( (dim, dim) )		
	for i in range( num * const_k ):
		Left[ i*4:i*4+4, i*4:i*4+4 ] = AA[:,:]
	
	R = zeros( ( 8*const_k, 2*const_k ) )
	for i in range( const_k ):
		R[i*4+3, i] = R[i*4+3, i+const_k] = 1
		R[i*4+2, i+const_k] = -1
		R[4*const_k+i*4, i] = R[4*const_k+i*4+1, i+const_k] = -1
		R[4*const_k+i*4, i+const_k] = 1
		
	for i in range( num-1 ): 
		Left[ 4*const_k*i:(4*i+8)*const_k, num*const_k*4+2*const_k*i:num*const_k*4+2*const_k*(i+1) ] = R
		Left[ num*const_k*4+2*const_k*i:num*const_k*4+2*const_k*(i+1), 4*const_k*i:(4*i+8)*const_k ] = R.T
		
	X = linalg.solve(Left, Right)	
	X = array( X[:num*const_k*4] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*const_k:(i+1)*const_k ] )		
	
	return result

'''
'''	
def compute_control_points_chain_with_second_derivative_continuity( controls, handles, transforms ):

	result = []
		
	return result