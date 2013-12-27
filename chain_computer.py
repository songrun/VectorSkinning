from copy import copy, deepcopy
from bezier_constraint_odd_solver import *
from bezier_constraint_even_solver import *

## The dimensions of a point represented in the homogeneous coordinates
dim = 2
	
def approximate_beziers(W_matrices, Cset, transforms, constraints, all_weights, all_vertices, all_indices, all_pts, if_closed=True):
	
	'''
	### 1 construct and solve the linear system for the odd iteration. if the constraints don't contain fixed angle and G1, skip ### 2.
	### 2 If the constraints contain any fixed angle and G1, iterate between the even and odd system until two consecutive solutions are close enough.
	### 3 refine the solutions based on the error of each curve. If it is larger than a threshold, split the curve into two.
	### 4 compute the bbw curves
	### 5 compute all the points along the curves.
	'''
	
	solutions = None
	## select the joint points' constraints
	joint_constraints = [constraints[key] for key in sorted(constraints.iterkeys())]
	joint_constraints = asarray(joint_constraints)

	### 1
	odd = BezierConstraintSolverOdd(W_matrices, Cset, joint_constraints, transforms, if_closed)
# 	odd.update_rhs_for_handles( transforms )
	last_solutions = solutions = odd.solve()

	### 2	
	if 2 in joint_constraints[:,0] or 4 in joint_constraints[:,0]: 

		even = BezierConstraintSolverEven(W_matrices, Cset, joint_constraints, transforms, if_closed)	
	# 		even.update_rhs_for_handles( transforms )

		for iter in xrange( 1 ):
			even.update_system_with_result_of_previous_iteration( solutions )
			last_solutions = solutions
			solutions = even.solve()
		
			if allclose(last_solutions, solutions, atol=1.0, rtol=1e-03):
				break
		
			## Check if error is low enough and terminate
			odd.update_system_with_result_of_previous_iteration( solutions )
			last_solutions = solutions
			solutions = odd.solve()
		
			if allclose(last_solutions, solutions, atol=0.5, rtol=1e-03):
				break

	### 3
	
	
	### 4 
	bbw_curves = []
	for indices in all_indices:
		tps = []	
		for i in indices:
			m = zeros(9)
			p = asarray(all_vertices[i] + [1.0])
			for h in range(len(transforms)):
				m = m + transforms[h]*all_weights[i,h]
		
			p = dot( m.reshape(3, 3), p.reshape(3,-1) ).reshape(-1)
			tps = tps + [p[0], p[1]] 	
		bbw_curves.append(tps)
	
	### 5
	spline_skin_curves = []
	for k, solution in enumerate(solutions):
		tps = []
		for t in asarray(all_pts)[k, 1]:
			tbar = asarray([t**3, t**2, t, 1.])
			p = dot(tbar, asarray( M * solution ) )
			tps = tps + [p[0], p[1]]
		spline_skin_curves.append(tps)
		
	return solutions, bbw_curves, spline_skin_curves
    	
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

	all_pts, all_dts = sample_cubic_bezier_curve_chain( Cset, 100 )
	from itertools import chain
	loop = list( chain( *[ samples for samples, ts in asarray(all_pts)[:,:,:-1] ] ) )

	all_vertices, facets, all_weights = (triangulate_and_compute_weights
													(loop, skeleton_handle_vertices))
	
	## all_indices is a table of all the indices of the points on the boundary 
	## in all_vertices
	all_indices = [range(len(pts)) for pts, ts in all_pts]
	last = 0
	for i in range(len(all_indices)):
		all_indices[i] = asarray(all_indices[i])+last
		last = all_indices[i][-1]
	all_indices[-1][-1] = all_indices[0][0]	
		
		
	W_matrices = zeros( ( len( Cset ), len( skeleton_handle_vertices ), 4, 4 ) )
	for k in xrange(len( Cset )):	
		for i in xrange(len( skeleton_handle_vertices )):
			## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
			## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
			W_matrices[k,i] = precompute_W_i_bbw( all_vertices, all_weights, i, all_pts[k][0], all_pts[k][1], all_dts[k])
	
							       
	constraints = {1: [1.0, 0.0], 4: [2.0, 0.0]}
	
	P_primes, bbw_curves, spline_skin_curves = approximate_beziers(W_matrices, Cset, trans, constraints, all_weights, all_vertices, all_indices, all_pts )
	
	print 'HAHA ~ '
	print P_primes
	
if __name__ == '__main__': main()		
