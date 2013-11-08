from weight_inverse_distance import *
from copy import copy, deepcopy
from bezier_utility import *

## The dimensions of a point represented in the homogeneous coordinates
dim = 3
## Two set of precomputed parameters
precomputes_1 = []
precomputes_2 = []

## update the two set of precomputed parameters
def update_precomputation_of_controls_or_handles(cps, handles):
	
	del precomputes_1[:]
	del precomputes_2[:]
	
	for i in range( len(handles ) ):
	
		precomputes_1.append( precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 100 ) )
		precomputes_2.append( precompute_partOfR( handles, i, controls[k*4:k*4+4], M, 0., 1., 100 ) )
		

## Compute the control points for each split curve, and return them in a list
## input: 	The 4-by-dim control points of the original bezier curve
## 			A list containing the fractions for each split
## output:	A list containing the same number of control points as the splits

def control_points_after_split( P, S ):
	
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


## A method the dispatch different specific methods to do the approximation 
## according to the constraint requirement indicated by constraint_level

def compute_control_points_chain_with_constraint( partition, original_controls, handles, transforms, constraint_level ):

	result = []
	
	Cset =  control_points_after_split( original_controls, partition ) 
	Cset = asarray( Cset )
	
# 	s = [0.] + cumsum(partition).tolist()
# 	print s
# 	for k in range( len( partition ) ):
# 		for i in range( len( handles ) ):
# 			print 'haha1: ', precompute_W_i( handles, i, original_controls, M, s[k], s[k+1], 100 )
# 			print 'haha2: ', precompute_W_i( handles, i, Cset[k*4:k*4+4], M, 0., 1., 100 )
	
	if constraint_level == 0:
		
 		result = compute_control_points_chain_without_constraint( Cset, handles, transforms, partition )

	elif constraint_level == 1:
		
		result = compute_control_points_chain_with_C0_continuity( Cset, handles, transforms, partition )
	
	elif constraint_level == 2:
	
		result = compute_control_points_chain_with_C1_continuity( Cset, handles, transforms, partition )
	
	elif constraint_level == 3:
			
#  		debugger()
		## precompute the Right part for G1 fixing directions
		right_fixing_dirs = precompute_rightside_for_fixing_directions( Cset, handles, transforms )
		
		right_fixing_mags = precompute_rightside_for_fixing_magnitudes( Cset, handles, transforms )
		
		result = compute_control_points_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, partition )

	return result


## Optimize the approximation without any constraint on joints and return the control points of each optimized split bezier in a list	
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


##	Optimize the approximation with all splits connected and return the control points of each optimized split bezier in a list.
##	Boundary Conditions are as follows:
## 		lambda1 * ( P41' - Q11' ) = 0
## 		lambda2 * ( P42' - Q12' ) = 0
## 		lambda3 * ( P43' - Q13' ) = 0
	
def compute_control_points_chain_with_C0_continuity( controls, handles, transforms, partitions ):

	temps = []
	result = []

	num = len( partitions ) # num is the number of splited curves
	base = dim * 4
	
# 	A = precompute_A( M, 0., 1., 50 )
	
	for k in range( num ):

		C = zeros( (dim, 4) )
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(dim,dim) )
			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 50 )	
	
			C = C + T_i * (controls[k*4:k*4+4].T) * M * mat(W_i) * M
	
		temps.append( asarray(C.T) )
	
	Right = []

	for R in temps:
		Right = Right + (R.T.reshape(-1)).tolist()
	Right = array( Right + zeros( dim * (num-1) ).tolist() )
	
	rank = len( Right )
# 	AA = M.T * A.T
	## AA is computed using Sage.
	AA = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,  3./35, 9./140,  1./35], [ 1./35, 9./140,  3./35,  1./14], [1./140,  1./35,  1./14,   1./7]] )

	assert AA.shape == (4,4)
	
	Left =  zeros( (rank, rank) )		
	for i in range( num * dim ):
		Left[ i*4:i*4+4, i*4:i*4+4 ] = AA[:,:]
	
	R = zeros( ( 2*base, dim ) )
	for i in range( dim ):
		R[i*4+3, i] = 1
		R[base + i*4, i] = -1
		
	for i in range( num-1 ): 
		Left[ base*i:(i+2)*base, num*base+dim*i:num*base+dim*(i+1) ] = R
		Left[ num*base+dim*i:num*base+dim*(i+1), base*i:(i+2)*base ] = R.T
	
	'''
	add weights
	'''
	for i in range( num ):
		Left[base*i:base*(i+1), base*i:base*(i+1)] *= partitions[i]
		Right[base*i:base*(i+1)] *= partitions[i]
		
	X = linalg.solve(Left, Right)	
	X = array( X[:num*base] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*dim:(i+1)*dim ] )		
		
	return result


## 	Optimize the approximation with constraints that all splits connected, 
##	the magnitudes of derivatives at joints are identical with the corresponding fraction, 
##	and the directions of derivatives at joints are collinear.
##	Then return the control points of each optimized split bezier in a list.
## 	Boundary Conditions are as follows:
## 		lambda1 * ( P41' - Q11' ) = 0
## 		lambda2 * ( P42' - Q12' ) = 0
## 		lambda3 * ( P43' - Q13' ) = 0
## 		lambda4 * ( P41' - P31' + Q11' - Q21') = 0
## 		lambda5 * ( P42' - P32' + Q12' - Q22') = 0
## 		lambda6 * ( P43' - P33' + Q13' - Q23') = 0
		
def compute_control_points_chain_with_C1_continuity( controls, handles, transforms, partitions ):
	
	num = len( partitions )
	mags = ones( (num, 2) )
	for i in range( len( partitions ) ):
		mags[i] *= partitions[i]
		
	right = precompute_rightside_for_fixing_magnitudes( controls, handles, transforms )	
	result = compute_control_points_chain_with_derivative_continuity_with_weight( right, partitions, mags )

	return result


## 	Optimize the approximation with constraints that all splits connected, 
##	a given magnitudes proportion of derivatives at joints, 
##	and the directions of derivatives at joints are collinear.
##	then return the control points of each optimized split bezier in a list.

def compute_control_points_chain_with_derivative_continuity_with_weight( right, partitions, mags ):

	result = []
	
	num = len( partitions ) # num is the number of splited curves
	base = dim * 4
	
	Right = deepcopy( right )
	rank = len( Right )
	
	## AA is computed using Sage.
	AA = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,  3./35, 9./140,  1./35], [ 1./35, 9./140,  3./35,  1./14], [1./140,  1./35,  1./14,   1./7]] )
	
	assert AA.shape == (4,4)
	
	Left =  zeros( (rank, rank) )		
	for i in range( num * dim ):
		Left[ i*4:i*4+4, i*4:i*4+4 ] = AA[:,:]
	
	R = zeros( ( 2*base, 2*dim ) )
	for i in range( dim ):
		R[i*4+3, i] = R[i*4+3, i+dim] = 1
		R[i*4+2, i+dim] = -1
		R[base+i*4, i] = R[base+i*4+1, i+dim] = -1
		R[base+i*4, i+dim] = 1
		
	'''
	add weights to lambda
	'''
	assert len(partitions) == num	
	for i in range( num-1 ): 
		R_copy = deepcopy( R )
		R_copy[ :base, dim: ] *= mags[i+1][0]
		R_copy[ base:, dim: ] *= mags[i][1]
		Left[ base*i:(i+2)*base, num*base+2*dim*i:num*base+2*dim*(i+1) ] = R_copy
		Left[ num*base+2*dim*i:num*base+2*dim*(i+1), base*i:(i+2)*base ] = R_copy.T
	
	'''
	add weights
	'''
	for i in range( num ):
		Left[base*i:base*(i+1), base*i:base*(i+1)] *= partitions[i]
		Right[base*i:base*(i+1)] *= partitions[i]

	X = linalg.solve(Left, Right)	
	X = array( X[:num*base] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*dim:(i+1)*dim ] )		
	
	return result


## 	Optimize the approximation with all splits connected and the directions of derivatives at joints collinear. 
##	then return the control points of each optimized split bezier in a list.
	
def compute_control_points_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, partitions, old_solution = None, mags = None, dirs = None, index = 0):
	
	assert mags is None or dirs is None
	
	if index >= 20:
		print 'stop at index: ', index
 		return old_solution
#		return [old_solution, index]
		
	rs = []
	num = len( partitions ) # num is the number of splited curves		
	
	if mags is not None:
		
		assert mags.shape == (num, 2)
 		dirs = ones( (num, 2, dim-1) )
			
		solution = compute_control_points_chain_with_derivative_continuity_with_weight( right_fixing_mags, partitions, mags )
		## check if the new solution is close enough to previous solution
		## if not, update the direction vector for each split
# 		error = 0.0
# 		for k in range( len( solution ) ):
# 	 		error += compute_error_metric(transforms, handles, solution[k], controls[k*4:k*4+4], M )*partitions[k]
# 	 	print 'mags error: ', error
	 	
		if old_solution is not None and allclose( old_solution, solution, atol = 0.1 ):
			print 'stop at index: ', index
#			return [solution, index]
 			return solution
		else:
			for i in range( num ):
				dirs[i][0] = dir_allow_zero( (solution[i][1]-solution[i][0])[:2] ) 
				dirs[i][1] = dir_allow_zero( (solution[i][2]-solution[i][3])[:2] ) 
	
#		return [solution, index+1]
 		return compute_control_points_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, partitions, old_solution = solution, dirs = dirs, index = index+1 )
		
	elif dirs is not None:	
	
		assert dirs.shape == (num, 2, dim-1 )
		mags = ones( (num, 2) )		
# 		debugger()
 		solution = compute_control_points_chain_fixing_directions( right_fixing_dirs, partitions, dirs[:, 0, :], dirs[:, 1, :])
		
		## get result of the iteration, and check if the result is stable
		
# 		error = 0.0	
# 		for k in range( len( solution ) ):
# 	 		error += compute_error_metric(transforms, handles, solution[k], controls[k*4:k*4+4], M )*partitions[k]
# 	 	print 'dirs error: ', error
		
		if old_solution is not None and allclose( old_solution, solution, atol = 0.1 ):
			print 'stop at index: ', index
#			return [solution, index]
  			return solution
			
		else:
			for i in range( num ):
				mags[i][0] = mag( (solution[i][1]-solution[i][0])[:2] )  
				mags[i][1] = mag( (solution[i][2]-solution[i][3])[:2] ) 

#		return [solution, index+1]
 		return compute_control_points_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, partitions, old_solution = solution, mags = mags, index = index+1 )

	else:	
	
		mags = ones( (num, 2) )
		for i in range( num ):
			mags[i] *= partitions[i]
			
 		return compute_control_points_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, partitions, mags = mags )
 		

## 	Optimize the approximation with contraints that all splits connected,
## 	A given magnitude proportion of derivatives at joints,
##	and the directions of derivatives at joints are collinear. 
##	then return the control points of each optimized split bezier in a list.
	
def compute_control_points_chain_fixing_directions( raw_rights, partitions, dir1, dir2): 

	dir1 = append( dir1, zeros((len( dir1 ),1)), axis = 1 )
	dir2 = append( dir2, zeros((len( dir2 ),1)), axis = 1 )
	
	num = len( partitions ) # num is the number of splited curves
	result = []
	base = (dim+1)*2
	
	Right = zeros( base*num + (num-1)*dim )
	for k in range( num ):
		
		temp = asarray( raw_rights[k] )
		Right[k*base : k*base + dim*2 : 2] = (temp[:,0] +  temp[:,1]).reshape(-1)
		Right[k*base + 1 : k*base + dim*2 : 2 ] = (temp[:,2] +  temp[:,3]).reshape(-1)
		Right[(k+1)*base-2] = dot( temp[:,1], dir1[k] )
		Right[(k+1)*base-1] = dot( temp[:,2], dir2[k] )		
	
	## AA1 and AA2 is computed using Sage.
	AA1 = array([[(13./35.), (9./70.)], [(9./70.), (13./35.)]])
	AA2 = array([[(11./70.), (13./140.)], [(13./140.), (11./70.)]])
	
	rank = len( Right )
	Left = zeros( ( rank, rank ) )

	for k in range( num ):
	
		for i in range( dim ):
			
			Left[ base*k+i*2:base*k+i*2+2, base*k+i*2:base*k+i*2+2 ] = AA1[:,:]
			AA2_copy = deepcopy( AA2 )
			AA2_copy[:, 0] *= dir1[k][i]
			AA2_copy[:, 1] *= dir2[k][i]
			Left[ base*k+i*2:base*k+i*2+2, base*(k+1)-2:base*(k+1) ] = AA2_copy[:,:]
			Left[ base*(k+1)-2:base*(k+1), base*k+i*2:base*k+i*2+2 ] = AA2_copy[:,:].T
		
		## The coefficients are computed using Sage.
		tmp = (9./140.)*dot(dir1[k], dir2[k])	
		Left[ base*(k+1)-2:base*(k+1), base*(k+1)-2:base*(k+1) ] = array([[(3./35.)*mag2(dir1[k]), tmp], [tmp, (3./35.)*mag2(dir2[k])]])
		
		## deal with the case when direction is zero
		if mag(dir1[k]) == 0:
			Left[ base*(k+1)-2, base*(k+1)-2 ] = 1.	
		if mag(dir2[k]) == 0:
			Left[ base*(k+1)-1, base*(k+1)-1 ] = 1.
			
	R = zeros( ( 2*base, dim ) )
	for i in range( dim ):
		R[2*i+1, i] = 1 # should be 0.5, but it doesn't affect the result
		R[base+2*i, i] = -1
		
	for i in range( num-1 ): 
		Left[ base*i:base*(i+2), num*base+dim*i:num*base+dim*(i+1) ] = R
		Left[ num*base+dim*i:num*base+dim*(i+1), base*i:base*(i+2) ] = R.T
	'''
	add weights
	'''
	for i in range( num ):
		Left[base*i:base*(i+1), base*i:base*(i+1)] *= partitions[i]
		Right[base*i:base*(i+1)] *= partitions[i]

	X = linalg.solve(Left, Right)
	solution = asarray(X[:base*num])
	
	for i in range( num ):
		P1, P4 = solution[base*i: base*i+dim*2 : 2], solution[base*i+1: base*i+dim*2 : 2]	
		P2, P3 = P1 + solution[base*(i+1)-2]*dir1[i], P4 + solution[base*(i+1)-1]*dir2[i]
		
		result.append(asarray([P1, P2, P3, P4]))
		
	return result	

## Precompute the right side for the G1 continuity for fixing directions		
def precompute_rightside_for_fixing_directions( controls, handles, transforms ):
	
	num = len( controls ) / 4
	result = []
	for k in range( num ):	
	
		temp = zeros( (dim, 4) )
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(dim, dim) )
 			partOfR = precompute_partOfR( handles, i, controls[k*4:k*4+4], M, 0., 1., 100 )		
			
			temp = temp + T_i * (controls[k*4:k*4+4].T) * M * mat(partOfR) 
		
		result.append( temp )
	
	result = asarray( result )
		
	return result

## Precompute the left side for the G1 continuity for fixing magnitudes
def precompute_rightside_for_fixing_magnitudes( controls, handles, transforms ):	
		
	temps = []
	num = len( controls ) / 4
	for k in range( num ):

		C = zeros( (dim, 4) )
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(dim,dim) )
 			W_i = precompute_W_i( handles, i, controls[k*4:k*4+4], M, 0., 1., 100 )	
	
			C = C + T_i * (controls[k*4:k*4+4].T) * M * mat( W_i ) * M
	
		temps.append( asarray(C.T) )
	
	Right = []

	for R in temps:
		Right = Right + (R.T.reshape(-1)).tolist()
	Right = array( Right + zeros( 2 * dim * (num-1) ).tolist() )
	
	return Right	