from weight_inverse_distance import *
from copy import copy, deepcopy
from bezier_utility import *
from bezier_constraint_odd_solver import *
from bezier_constraint_even_solver import *

## The dimensions of a point represented in the homogeneous coordinates
dim = 2
	
def approximate_beziers(W_matrices, Cset, transforms, constraints, if_closed):
	
	solution = None
	
	## select the joint points' constraints
	joint_constraints = [constraints[key] for key in sorted(constraints.iterkeys())]
	joint_constraints = asarray(joint_constraints)

	odd = BezierConstraintSolverOdd(W_matrices, Cset, joint_constraints, transforms, if_closed)
# 	odd.update_rhs_for_handles( transforms )
	solution = odd.solve()
	
	if 2 not in joint_constraints[:,0] and 4 not in joint_constraints[:,0]: 
		return solution
	
	else: 
		even = BezierConstraintSolverEven(W_matrices, Cset, joint_constraints, transforms, if_closed)	
# 		even.update_rhs_for_handles( transforms )
	
		for iter in xrange( 10 ):
			even.update_system_with_result_of_previous_iteration( solution )
			solution = even.solve()

			## Check if error is low enough and terminate
			odd.update_system_with_result_of_previous_iteration( solution )
			solution = odd.solve()
    
    	return solution
    	
def main():
	cps = array([[100, 300,   1],
	   [200, 400,   1],
	   [300, 400,   1],
	   [400, 300,   1],
	   [300, 200,   1],
	   [200, 200,   1]])
	skeleton_handle_vertices = [[200.0, 300.0, 1.0], [300.0, 300.0, 1.0]] 
	trans = [array([ 1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.]), array([ 1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.])]	  
	
	Cset = make_control_points_chain( cps, True )  

	all_pts = sample_cubic_bezier_curve_chain( Cset, 100 )
	from itertools import chain
	loop = list( chain( *[ samples for samples, ts in asarray(all_pts)[:,:,:-1] ] ) )

	all_vertices, facets, all_weights = (triangulate_and_compute_weights
													(loop, skeleton_handle_vertices))
													
# 	debugger()
	## boundaries is a table of all the indices of the points on the boundary 
	## in all_vertices
	boundaries = [range(len(pts)) for pts, ts in all_pts]
	last = 0
	for i in range(len(boundaries)):
		boundaries[i] = asarray(boundaries[i])+last
		last = boundaries[i][-1]
	boundaries[-1][-1] = boundaries[0][0]	
	
	W_matrices = zeros( ( len( Cset ), len( skeleton_handle_vertices ), 2, 4, 4 ) )
	for k in xrange(len( Cset )):	
		for i in xrange(len( skeleton_handle_vertices )):
			## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
			## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
			W_matrices[k,i] = precompute_W_i_bbw( all_vertices, 
									all_weights, i, all_pts[k][0], all_pts[k][1])
	
							       
	constraints = {1: [4.0, 1.0], 4: [4.0, 0.0]}
	
	P_primes = approximate_beziers(W_matrices, Cset, trans, constraints, True )
	
	print 'HAHA ~ '
	print P_primes
	
if __name__ == '__main__': main()		





# 	for P in temps:
# 		Right = Right + (P.T.reshape(-1)).tolist()	
# 	if(	closed ): 
# 		Right = array( Right + zeros( 2 * dim * num ).tolist() )
# 	else:
# 		Right = array( Right + zeros( 2 * dim * (num-1) ).tolist() )
	
	# 	## For close loop
# 	if(	array_equal(controls[0][0], controls[-1][-1]) ):
# 		Left[ :base, -dim:] = R[base:]
# 		Left[ (num-1)*base:num*base, -dim: ] = R[:base]
# 		Left[ -dim:, :base ] = R[base:].T
# 		Left[ -dim:, (num-1)*base:num*base ] = R[:base].T
	
	
	
	
# 	num = len( partitions )
# 	## mags is the magnitudes of the derivatives of each curve at endpoints. 
# 	## The mags are proportional to the lengths of each curve for C1 continuity.
# 	mags = ones( (num, 2) )
# 	for i in range( len( partitions ) ):
# 		mags[i] *= partitions[i]
# 	closed = False
# 	if(	array_equal(controls[0][0], controls[-1][-1]) ): closed = True
# 	
# 	
# 	right = compute_rightside_for_general(W_matrices, controls, handles, transforms, closed )	
# 	result = approximate_bezier_chain_with_G1_magnitude_fixed( right, partitions, mags, closed )
# 
# 	return result


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

