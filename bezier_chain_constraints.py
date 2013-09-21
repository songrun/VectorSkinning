from numpy import *
from weight_inverse_distance import *
from copy import copy, deepcopy
from myarray import *

try:
   from pydb import debugger

   ## Also add an exception hook.
   import pydb, sys
   sys.excepthook = pydb.exception_hook

except ImportError:
   import pdb
   def debugger():
       pdb.set_trace()


M = matrix('-1. 3. -3. 1.; 3. -6. 3. 0.; -3. 3. 0. 0.; 1. 0. 0. 0.')

def control_points_after_split( P, S ):
	
# 	fake = asarray([[200., 200., 1.], [200., 400., 1.], [600., 400., 1.], [600., 200., 1.]])
# 	P = fake
	
	assert len( P.shape ) == 2
	assert P.shape[0] == 4
	P = asarray( P )
	
	S = asarray( S )
	
	ref = 1.0
	for i in range( len(S) ):
		S[i] = S[i] / ref
		ref = ( 1 - S[i] ) * ref
	
	result = []
	for i, k in enumerate( S[:-1] ):
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
# 	Cset = Cset + original_controls.tolist()
	Cset = asarray( Cset )
	
	if constraint_level == 0:
		
 		result = compute_control_points_chain_without_constraint( Cset, handles, transforms, partition )
#  		result = compute_control_points_chain_without_constraint_2( original_controls, partition, handles, transforms )	
	elif constraint_level == 1:
		
		result = compute_control_points_chain_with_C0_continuity( Cset, handles, transforms, partition )
	
	elif constraint_level == 2:
	
		result = compute_control_points_chain_with_C1_continuity( Cset, handles, transforms, partition )
	
	elif constraint_level == 3:
		
		result = compute_control_points_chain_with_G1_continuity( Cset, handles, transforms, partition )

	return result
	
def compute_control_points_chain_without_constraint( controls, handles, transforms, partitions ): 

	result = []
	num = len( partitions )
	
	inv_A = linalg.inv( precompute_A( M, 0., 1., 50 ) )					
	for k in range( num ):
		
		P_prime = zeros(12).reshape(3,4)
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(3,3) )
			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )		
			
			P_prime = P_prime + T_i * (controls[k*4:k*4+4].T) * M * mat(W_i) * inv_A	

		result.append( asarray(P_prime.T) )
	
	return result

'''
	Compute the control points for each splited curve with endpoints connected.
	Boundary Conditions are as follows:
		lambda1 * ( P41' - Q11' ) = 0
		lambda2 * ( P42' - Q12' ) = 0
		lambda3 * ( P43' - Q13' ) = 0
'''	
def compute_control_points_chain_with_C0_continuity( controls, handles, transforms, partitions ):

	temps = []
	result = []
	
	const_k = 3
	num = len( controls ) /4 # num is the number of splited curves
	
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
	
	'''
	add weights
	'''
	for i in range( num ):
		Left[4*const_k*i:4*(i+1)*const_k, 4*const_k*i:4*(i+1)*const_k] *= partitions[i]
		Right[4*const_k*i:4*(i+1)*const_k] *= partitions[i]
		
	X = linalg.solve(Left, Right)	
	X = array( X[:num*const_k*4] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*const_k:(i+1)*const_k ] )		
	
	'''
	for test
	'''
	
	
	
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
def compute_control_points_chain_with_C1_continuity( controls, handles, transforms, partitions ):
	
	num = len( partitions )
	mags = ones( (num, 2) )
	for i in range( len( partitions ) ):
		mags[i] *= partitions[i]
	print 'mags: ', mags
	result = compute_control_points_chain_with_derivative_continuity_with_weight( controls, handles, transforms, partitions, mags )
	
	return result

def compute_control_points_chain_with_derivative_continuity_with_weight( controls, handles, transforms, partitions, mags ):

	temps = []
	result = []
	
	const_k = 3
	num = len( controls ) /4 # num is the number of splited curves
	
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
		
	'''
	add weights to lambda
	'''
	assert len(partitions) == num	
	for i in range( num-1 ): 
		R_copy = deepcopy( R )	
		R_copy[ :4*const_k, const_k: ] /= mags[i][1]
		R_copy[ 4*const_k:, const_k: ] /= mags[i+1][0]
		print 'R', R, 'R_copy ', R_copy
		Left[ 4*const_k*i:(4*i+8)*const_k, num*const_k*4+2*const_k*i:num*const_k*4+2*const_k*(i+1) ] = R_copy
		Left[ num*const_k*4+2*const_k*i:num*const_k*4+2*const_k*(i+1), 4*const_k*i:(4*i+8)*const_k ] = R_copy.T
	
	'''
	add weights
	'''
	for i in range( num ):
		Left[4*const_k*i:4*(i+1)*const_k, 4*const_k*i:4*(i+1)*const_k] *= partitions[i]
		Right[4*const_k*i:4*(i+1)*const_k] *= partitions[i]
		
	X = linalg.solve(Left, Right)	
	X = array( X[:num*const_k*4] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*const_k:(i+1)*const_k ] )		
	
	return result

'''
	Compute the control points for each splited curve with endpoints connected and the direction of the derivative at the endpoints collineared.
'''	
def compute_control_points_chain_with_G1_continuity( controls, handles, transforms, partitions ):
	
	result = []
	const_k = 3
	num = len( partitions ) # num is the number of splited curves
	
	base = (const_k+1)*2	
	solution = zeros( base*num )
	mags = ones( (num, 2) )
	for i in range( len( partitions ) ):
		mags[i] *= partitions[i]
	v = []
	r = []
	
# 	debugger()
	count = 0
	while True:
		count += 1
		if count == 20:
			print 'stop at count: ', count
			result = []
			
			for i in range( num ):
				P1 = solution[base*i: base*i+const_k*2 : 2]
				P4 = solution[base*i+1: base*i+const_k*2 : 2]
				P2 = P1 + solution[base*(i+1)-2]*v[i]
				P3 = P4 + solution[base*(i+1)-1]*r[i]
				#print solution, P1, P2, P3, P4
				result.append( array([P1, P2, P3, P4]) )
			break
			
		del v[:]
		del r[:]
		
		result = compute_control_points_chain_with_derivative_continuity_with_weight( controls, handles, transforms, partitions, mags )
		
		old_solution = deepcopy( solution )
		for i in range( num ):
			solution[base*i: base*i+const_k*2 : 2] = result[i][0].reshape(-1) 
			solution[base*i+1: base*i+const_k*2 : 2] = result[i][3].reshape(-1)
			solution[base*(i+1)-2] = mag( result[i][1]-result[i][0] ) 		
			solution[base*(i+1)-1] = mag( result[i][3]-result[i][2] )
		## check if the new solution is close enough to previous solution
		## if not, update the direction vector for each split
		if allclose( old_solution, solution, atol = 0.01 ):
			print 'stop at count: ', count
			break
		else:
			for i in range( num ):
				v.append( dir( result[i][1]-result[i][0] ) )
				r.append( dir( result[i][3]-result[i][2] ) )		
	
		Right = zeros( base*num + (num-1)*const_k )
# 		Right = zeros( base*num )
		for k in range( num ):	
		
			temp = zeros( (const_k, 4) )
			for i in range( len( handles ) ):

				T_i = mat( asarray(transforms[i]).reshape(const_k, const_k) )
				partOfR = precompute_partOfR( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )		
				
 				temp = temp + T_i * (controls[k*4:k*4+4].T) * M * mat(partOfR)
			
			Right[k*base : k*base + const_k*2 : 2] = asarray(temp[:,0] +  temp[:,1]).reshape(-1)
			Right[k*base + 1 : k*base + const_k*2 : 2 ] = asarray(temp[:,2] +  temp[:,3]).reshape(-1)
			Right[(k+1)*base-2] = dot( asarray(temp[:,1]).reshape(-1), v[k] )
			Right[(k+1)*base-1] = dot( asarray(temp[:,2]).reshape(-1), r[k] )	
	
# 		print 'G1 Right: ', Right.shape, Right
		
		AA1 = array([[(13./35.), (9./70.)], [(9./70.), (13./35.)]])
		AA2 = array([[(11./70.), (13./140.)], [(13./140.), (11./70.)]])
		
		dim = len( Right )
		Left = zeros( ( dim, dim ) )
		
		for k in range( num ):
		
			for i in range( const_k ):
				
				Left[ base*k+i*2:base*k+i*2+2, base*k+i*2:base*k+i*2+2 ] = AA1[:,:]
				AA2_copy = deepcopy( AA2 )
				AA2_copy[:, 0] *= v[k][i]
				AA2_copy[:, 1] *= r[k][i]
				Left[ base*k+i*2:base*k+i*2+2, base*(k+1)-2:base*(k+1) ] = AA2_copy[:,:]
				Left[ base*(k+1)-2:base*(k+1), base*k+i*2:base*k+i*2+2 ] = AA2_copy[:,:].T
			
			tmp = (9./140.)*dot(v[k], r[k])	
			Left[ base*(k+1)-2:base*(k+1), base*(k+1)-2:base*(k+1) ] = array([[(3./35.), tmp], [tmp, (3./35.)]])
			
		R = zeros( ( 2*base, const_k ) )
		for i in range( const_k ):
			R[2*i+1, i] = 1
			R[base+2*i, i] = -1	
			
		for i in range( num-1 ): 
			Left[ base*i:base*(i+2), num*base+const_k*i:num*base+const_k*(i+1) ] = R
			Left[ num*base+const_k*i:num*base+const_k*(i+1), base*i:base*(i+2) ] = R.T
		'''
		add weights
		'''
		for i in range( num ):
			Left[base*i:base*(i+1), base*i:base*(i+1)] *= partitions[i]
			Right[base*i:base*(i+1)] *= partitions[i]
	
		X = linalg.solve(Left, Right)	
		
		if allclose( asarray(X[:base*num]), solution, atol = 0.01 ):
			print 'stop at count: ', count
			result = []
			for i in range( num ):
				P1 = solution[base*i: base*i+const_k*2 : 2]
				P4 = solution[base*i+1: base*i+const_k*2 : 2]
				P2 = P1 + solution[base*(i+1)-2]*v[i]
				P3 = P4 + solution[base*(i+1)-1]*r[i]
				result.append( array([P1, P2, P3, P4]) )
			break
			
		else:
			solution = asarray(X[:base*num])
			for i in range( num ):
				mags[i, :] = X[base*(i+1)-2:base*(i+1)]  
			
	print 'result: ', result
	return result