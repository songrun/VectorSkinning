from copy import copy, deepcopy
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
	
	even = BezierConstraintSolverEven(W_matrices, Cset, joint_constraints, transforms, if_closed)	
# 		even.update_rhs_for_handles( transforms )

	for iter in xrange( 1 ):
		even.update_system_with_result_of_previous_iteration( solution )
		last_solution = solution
		solution = even.solve()
		
# 		debugger()
		if allclose(last_solution, solution, atol=1.0, rtol=1e-03):
			return solution
		
		## Check if error is low enough and terminate
		odd.update_system_with_result_of_previous_iteration( solution )
		last_solution = solution
		solution = odd.solve()
		
		if allclose(last_solution, solution, atol=0.5, rtol=1e-03):
			return solution

	return solution
    	
def main():
	Cset = array([[[100, 300,   1],
        [200, 400,   1],
        [300, 400,   1],
        [400, 300,   1]],

       [[400, 300,   1],
        [300, 200,   1],
        [200, 200,   1],
        [100, 300,   1]]])
        

	skeleton_handle_vertices = [[200.0, 300.0, 1.0], [300.0, 300.0, 1.0]] 
	trans = [array([ 1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.]), array([ 1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.])]	  

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
	
	W_matrices = zeros( ( len( Cset ), len( skeleton_handle_vertices ), 4, 4 ) )
	for k in xrange(len( Cset )):	
		for i in xrange(len( skeleton_handle_vertices )):
			## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
			## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
			W_matrices[k,i] = precompute_W_i_bbw( all_vertices, all_weights, i, all_pts[k][0], all_pts[k][1])
	
							       
	constraints = {1: [1.0, 0.0], 4: [2.0, 0.0]}
	
	P_primes = approximate_beziers(W_matrices, Cset, trans, constraints, True )
	
	print 'HAHA ~ '
	print P_primes
	
if __name__ == '__main__': main()		
