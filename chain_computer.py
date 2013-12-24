from weights_computer import *
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
	last_solution = solution = odd.solve()
	
	if 2 not in joint_constraints[:,0] and 4 not in joint_constraints[:,0]: 
		return solution
	
	else: 
		even = BezierConstraintSolverEven(W_matrices, Cset, joint_constraints, transforms, if_closed)	
# 		even.update_rhs_for_handles( transforms )
	
		for iter in range( 10 ):
			even.update_system_with_result_of_previous_iteration( solution )
			last_solution = solution
			solution = even.solve()
			
			if allclose(last_solution, solution, atol=1):
				print 'odd even iter num: ', iter
				return solution
				
			## Check if error is low enough and terminate
			odd.update_system_with_result_of_previous_iteration( solution )
			last_solution = solution
			solution = odd.solve()
			
			if allclose(last_solution, solution, atol=1):
				print 'even odd iter num: ', iter
				return solution
				
    	print 'exceed: ', iter
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

