from weight_inverse_distance import *
from copy import copy, deepcopy
from bezier_utility import *

## The dimensions of a point represented in the homogeneous coordinates
dim = 3
	
def approximate_beziers(W_matrices, Cset, handles, transforms, constraints, if_closed):
	
	Left = []
	Right = []
	## compute the weight of each segment according to its length
	num = len(Cset)
	partition = asarray([length_of_cubic_bezier_curve(P) for P in Cset])
	partition = partition/sum(partition)
	
	## select the joint points' constraints
	constraints = [constraints[key] for key in sorted(constraints.iterkeys())[: :3]]
	base = dim * 4
	rank = base * num
	
	for i in range(num):
		Left = construct_coefficients_for_general_constraint(Left)
		Right = construct_solution_for_general_constraint(Right, W_matrices, i, Cset[i], 
				handles, transforms, scale=1.0)
	
# 	elif constraint_level == 3:			
# 		## precompute the Right part for G1 fixing directions
# 		right_fixing_dirs = precompute_rightside_for_fixing_directions(W_matrices, Cset, handles, transforms )	
# 		
# 		result = approximate_bezier_chain_with_G1_continuity( right_fixing_dirs, right, partition )

	X = linalg.solve(Left, Right)	
	X = array( X[:rank] ).reshape(-1,4).T
	
	result = []
	for i in range(num):
		result.append( X[:, i*dim:(i+1)*dim ] )		
		
	return result		

## Construct the coefficient matrix for one plain curve
def construct_coefficients_for_general_constraint(Left, scale=1.0):
	'''
	## A is computed using Sage, integral of (tbar.T * tbar) with respect to t.
	# 	A = asarray( [[  1./7,  1./6,  1./5, 1./4], [ 1./6,  1./5, 1./4,  1./3], 
	# 		[ 1./5, 1./4,  1./3,  1./2], [1./4,  1./3,  1./2,   1.]] )
	## MAM is computed using Sage. MAM = M * A * M
	'''
	MAM = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,  3./35, 9./140,  1./35], 
		[ 1./35, 9./140,  3./35,  1./14], [1./140,  1./35,  1./14,   1./7]] )
		
	Left = asarray(Left)
	if len(Left) != 0:
		assert len(Left.shape) == 2
		assert Left.shape[0] == Left.shape[1]
	old_rank = Left.shape[0]
	
	rank = old_rank + 4*dim
	result = zeros((rank, rank))
	result[:old_rank, :old_rank] = Left[:]
	
	for i in range(dim):		
		result[ old_rank+i*4:old_rank+(i+1)*4, old_rank+i*4:old_rank+(i+1)*4 ] = MAM[:,:]*scale
		
	return result	

## construct the solution for the general case including G1 continuity for fixing magnitudes
def construct_solution_for_general_constraint(Right, W_matrices, k, controls, handles, 
											transforms, scale=1.0):	
	Right = asarray(Right)
	assert len(Right.shape) == 1
	old_rank = Right.shape[0]
	
	rank = old_rank + 4*dim
	result = zeros(rank)
	result[:old_rank] = Right[:]
	
	P = zeros( (dim, 4) )
	for i in range( len( handles ) ):

		T_i = mat( asarray(transforms[i]).reshape(dim,dim) )
		W_i = W_matrices[k,i,0]	

		P = P + T_i * (controls.T) * M * mat( W_i ) * M
	
	result[old_rank:] = asarray(P).reshape(-1)*scale	
	
	return result	


# 	for P in temps:
# 		Right = Right + (P.T.reshape(-1)).tolist()	
# 	if(	closed ): 
# 		Right = array( Right + zeros( 2 * dim * num ).tolist() )
# 	else:
# 		Right = array( Right + zeros( 2 * dim * (num-1) ).tolist() )

## Optimize the approximation without any constraint on joints and return the control points of each optimized split bezier in a list	
def approximate_bezier_chain_without_constraint(left, right): 
	'''
	left is a square matrix of existed equation coefficients n-by-n
	right is the solutions of a linear system with rank n
	'''
	old_rank = len( right )
	rank = old_rank + 4*dim + dim 
	## make a new coefficient matrix and a solution array and copy the existed values.
	Left =  zeros((rank, rank))
	Right = zeros(rank)
	Left[:old_rank, :old_rank] = left[:]
	Right[:old_rank] = right[:]
	
	for i in range(dim):		
		result[ old_rank+i*4:old_rank+(i+1)*4, old_rank+i*4:old_rank+(i+1)*4 ] = MAM[:,:]
	
	return Left, Right


##	Optimize the approximation with all splits connected and return the control points of each optimized split bezier in a list.
##	Boundary Conditions are as follows:
## 		lambda1 * ( P41' - Q11' ) = 0
## 		lambda2 * ( P42' - Q12' ) = 0
## 		lambda3 * ( P43' - Q13' ) = 0
	
def approximate_bezier_chain_with_C0_continuity(left, right, portion ):
	'''
	left is a square matrix of existed equation coefficients n-by-n
	right is the solutions of a linear system with rank n
	'''
	old_rank = len( right )
	rank = old_rank + 4*dim + dim 
	base = dim * 4
	## make a new coefficient matrix and a solution array and copy the existed values.
	Left =  zeros((rank, rank))
	Right = zeros(rank)
	Left[:old_rank, :old_rank] = left[:]
	Right[:old_rank] = right[:]
	
	for i in range(dim):		
		result[ old_rank+i*4:old_rank+(i+1)*4, old_rank+i*4:old_rank+(i+1)*4 ] = MAM[:,:]*scale
	
	## add lambdas
	R = zeros( ( 2*base, dim ) )
	for i in range( dim ):
		R[i*4+3, i] = 1
		R[base + i*4, i] = -1		 
	Left[ -dim-2*base:-dim, -dim: ] = R
	Left[ -dim:, -dim-2*base:-dim ] = R.T
	
# 	## For close loop
# 	if(	array_equal(controls[0][0], controls[-1][-1]) ):
# 		Left[ :base, -dim:] = R[base:]
# 		Left[ (num-1)*base:num*base, -dim: ] = R[:base]
# 		Left[ -dim:, :base ] = R[base:].T
# 		Left[ -dim:, (num-1)*base:num*base ] = R[:base].T
		
	return Left, Right


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
		
def approximate_bezier_chain_with_C1_continuity(W_matrices, controls, handles, transforms, partitions ):
	
	num = len( partitions )
	## mags is the magnitudes of the derivatives of each curve at endpoints. 
	## The mags are proportional to the lengths of each curve for C1 continuity.
	mags = ones( (num, 2) )
	for i in range( len( partitions ) ):
		mags[i] *= partitions[i]
	closed = False
	if(	array_equal(controls[0][0], controls[-1][-1]) ): closed = True
	
	
	right = compute_rightside_for_general(W_matrices, controls, handles, transforms, closed )	
	result = approximate_bezier_chain_with_G1_magnitude_fixed( right, partitions, mags, closed )

	return result


## 	Optimize the approximation with constraints that all splits connected, 
##	a given magnitudes proportion of derivatives at joints, 
##	and the directions of derivatives at joints are collinear.
##	then return the control points of each optimized split bezier in a list.

def approximate_bezier_chain_with_G1_magnitude_fixed(right, partitions, mags, closed = True ):

	result = []
	
	num = len( partitions ) # num is the number of splited curves
	base = dim * 4
	
	Right = deepcopy( right )
	rank = len( Right )
	
	Left =  zeros( (rank, rank) )		
	for i in range( num * dim ):
		Left[ i*4:i*4+4, i*4:i*4+4 ] = MAM[:,:]
	
	## Apply Langrange Multiplier to add constraints.
	R = zeros( ( 2*base, 2*dim ) )
	for i in range( dim ):
		R[i*4+3, i] = R[i*4+3, i+dim] = 1
		R[i*4+2, i+dim] = -1
		R[base+i*4, i] = R[base+i*4+1, i+dim] = -1
		R[base+i*4, i+dim] = 1
		
	## add weights to lambda
	assert len(partitions) == num	
	for i in range( num-1 ): 
		R_copy = deepcopy( R )
		R_copy[ :base, dim: ] *= mags[i+1][0]
		R_copy[ base:, dim: ] *= mags[i][1]
		Left[ base*i:(i+2)*base, num*base+2*dim*i:num*base+2*dim*(i+1) ] = R_copy
		Left[ num*base+2*dim*i:num*base+2*dim*(i+1), base*i:(i+2)*base ] = R_copy.T
		
	## For close loop
	if(	closed ):
		R[ :base, dim: ] *= mags[0,0]
		R[ base:, dim: ] *= mags[-1,1]
		Left[ :base, -2*dim:] = R[base:]
		Left[ (num-1)*base:num*base, -2*dim: ] = R[:base]
		Left[ -2*dim:, :base ] = R[base:].T
		Left[ -2*dim:, (num-1)*base:num*base ] = R[:base].T
	
	## add weights
	for i in range( num ):
		Left[base*i:base*(i+1), base*i:base*(i+1)] *= partitions[i]
		Right[base*i:base*(i+1)] *= partitions[i]

	X = linalg.solve(Left, Right)	
	X = array( X[:num*base] ).reshape(-1,4).T
	
	for i in range(num):
		result.append( X[:, i*dim:(i+1)*dim ] )		
	
	return result

## 	Optimize the approximation with contraints that all splits connected,
## 	A given magnitude proportion of derivatives at joints,
##	and the directions of derivatives at joints are collinear. 
##	then return the control points of each optimized split bezier in a list.
	
def approximate_chain_with_G1_direction_fixed( raw_rights, partitions, dir1, dir2): 

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

## 	Optimize the approximation with all splits connected and the directions of derivatives at joints collinear. 
##	then return the control points of each optimized split bezier in a list.
	
def approximate_bezier_chain_with_G1_continuity( right_fixing_dirs, right_fixing_mags, 
	partitions, old_solution = None, mags = None, dirs = None, index = 0):
	
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
			
		solution = approximate_bezier_chain_with_G1_magnitude_fixed( 
					right_fixing_mags, partitions, mags )
		## check if the new solution is close enough to previous solution
		## if not, update the direction vector for each split
	 	
		if old_solution is not None and allclose( old_solution, solution, atol = 0.1 ):
			print 'stop at index: ', index
#			return [solution, index]
 			return solution
		else:
			for i in range( num ):
				dirs[i][0] = dir_allow_zero( (solution[i][1]-solution[i][0])[:2] ) 
				dirs[i][1] = dir_allow_zero( (solution[i][2]-solution[i][3])[:2] ) 
	
#		return [solution, index+1]
 		return approximate_bezier_chain_with_G1_continuity( right_fixing_dirs, 
 		right_fixing_mags, partitions, old_solution = solution, dirs = dirs, index = index+1 )
		
	elif dirs is not None:	
	
		assert dirs.shape == (num, 2, dim-1 )
		mags = ones( (num, 2) )		
# 		debugger()
 		solution = approximate_chain_with_G1_direction_fixed( right_fixing_dirs, 
 					partitions, dirs[:, 0, :], dirs[:, 1, :])
		
		## get result of the iteration, and check if the result is stable
		
		if old_solution is not None and allclose( old_solution, solution, atol = 0.1 ):
			print 'stop at index: ', index
#			return [solution, index]
  			return solution
			
		else:
			for i in range( num ):
				mags[i][0] = mag( (solution[i][1]-solution[i][0])[:2] )  
				mags[i][1] = mag( (solution[i][2]-solution[i][3])[:2] ) 

#		return [solution, index+1]
 		return approximate_bezier_chain_with_G1_continuity( right_fixing_dirs, 
 		right_fixing_mags, partitions, old_solution = solution, mags = mags, index = index+1 )

	else:	
	
		mags = ones( (num, 2) )
		for i in range( num ):
			mags[i] *= partitions[i]
			
 		return approximate_bezier_chain_with_G1_continuity( right_fixing_dirs, 
 				right_fixing_mags, partitions, mags = mags )
 			

## Precompute the right side for the G1 continuity for fixing directions		
def precompute_rightside_for_fixing_directions(W_matrices, controls, handles, transforms ):
	
	num = len( controls )
	result = []
	for k in range( num ):	
	
		temp = zeros( (dim, 4) )
		for i in range( len( handles ) ):

			T_i = mat( asarray(transforms[i]).reshape(dim, dim) )
 			partOfR = W_matrices[k,i,1]		
			
			temp = temp + T_i * (controls[k].T) * M * mat(partOfR) 
		
		result.append( temp )
	
	result = asarray( result )
		
	return result

