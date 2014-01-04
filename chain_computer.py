from copy import copy, deepcopy
from bezier_constraint_odd_solver import *
from bezier_constraint_even_solver import *


class Control_point:
	'''
	Class of control point.
	id is the control point's id.
	position is its relative position with respect to the canvas
	is_joint tells if it is an joint point. If so, it has constraint whose first element corresponds to one of the four types of C^0, fixed angle, C^1 and G^1, and second element indicates if its position is fixed.
	'''
	id = -1
	position = zeros(2)
	is_joint = False
	constraint = None
	
	def __init__(self, id, pos, is_joint=False, constraint = [1,0] ):
		pos = asarray( pos )
		constraint = asarray( constraint )
		assert len( pos ) == 2
		assert len( constraint ) == 2
	
		self.id = id
		self.position = pos
		self.is_joint = is_joint
		if is_joint:
			self.constraint = constraint
			
	def compare_shape( compared_control ):
		assert isinstance( compared_control, Handle )
		if compared_control.id == self.id and array_equal( compared_control.position, self.position) and compared_control.is_joint == self.is_joint and array_equal( compared_control.constraint, self.constraint ) :
			return True
		else: 
			return False

def get_controls( controls ):
	'''
	given a list of Control_point classes, return each control point's position and the joint's constraint.
	'''
	control_pos = []
	constraints = []
	for control in controls:
		control_pos.append( control.position )
		if control.is_joint:
			assert control.constraint is not None
			constraints.append( control.constraint )
		
	control_pos = make_control_points_chain( control_pos )

	return asarray(control_pos), asarray(constraints)

## The dimensions of a point represented in the homogeneous coordinates
dim = 2
	
def approximate_beziers(W_matrices, controls, handles, transforms, all_weights, all_vertices, all_indices, all_pts, all_dts, enable_refinement=True):
	
	'''
	### 1 construct and solve the linear system for the odd iteration. if the constraints don't contain fixed angle and G1, skip ### 2.
	### 2 If the constraints contain any fixed angle and G1, iterate between the even and odd system until two consecutive solutions are close enough.
	### 3 compute the bbw curves
	### 4 compute all the points along the curves.
	### 5 refine the solutions based on the error of each curve. If it is larger than a threshold, split the curve into two.
	'''
	
	solutions = None
	control_pos, constraints = get_controls( controls )
	control_pos = concatenate((control_pos, ones((control_pos.shape[0],4,1))), axis=2)

	### 1
	odd = BezierConstraintSolverOdd(W_matrices, control_pos, constraints, transforms )
#	odd.update_rhs_for_handles( transforms )
	last_solutions = solutions = odd.solve()

	### 2	
	if 2 in constraints[:,0] or 4 in constraints[:,0]: 

		even = BezierConstraintSolverEven(W_matrices, control_pos, constraints, transforms )	
	#		even.update_rhs_for_handles( transforms )

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
	
	### 4
	spline_skin_curves = []
	for k, solution in enumerate(solutions):
		tps = []
		for t in asarray(all_pts)[k, 1]:
			tbar = asarray([t**3, t**2, t, 1.])
			p = dot(tbar, asarray( M * solution ) )
			tps = tps + [p[0], p[1]]
		spline_skin_curves.append(tps)
	
	
	### 5
	new_controls = adapt_configuration_based_on_diffs( controls, bbw_curves, spline_skin_curves, all_dts )
	
	if enable_refinement and len( new_controls ) > len( controls ):	
# 		debugger()
		new_control_pos = get_controls( new_controls )[0]

		W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts = precompute_all_when_configuration_change( new_control_pos, handles  )
	
		solutions, bbw_curves, spline_skin_curves = approximate_beziers(W_matrices, new_controls, handles, transforms, all_weights, all_vertices, all_indices, all_pts, all_dts, False)	
	
	return solutions, bbw_curves, spline_skin_curves


def adapt_configuration_based_on_diffs( controls, bbw_curves, spline_skin_curves, all_dts ):
	'''
	 sample the bezier curve solution from optimization at the same "t" locations as bbw-affected curves. Find the squared distance between each corresponding point, multiply by the corresponding "dt", and sum that up. That's the energy. Then scale it by the arc length of each curve.
	'''
	assert len( bbw_curves ) == len( spline_skin_curves )
	diffs = [compute_error_metric(bbw_curve, spline_skine_curve, dts) for bbw_curve, spline_skine_curve, dts in zip(bbw_curves, spline_skin_curves, all_dts) ]
	print 'differences: ', diffs
	
	new_controls = []
	partition = [0.5, 0.5]
	threshold = 100 
	
	all_pos = asarray([x.position for x in controls])
	
	for k, diff in enumerate( diffs ):
		control_pos = all_pos[ k*3 : k*3+4 ]
		if len(control_pos) == 3:	
			control_pos = concatenate((control_pos, all_pos[0].reshape(1,2)))
		
		if diff > threshold*length_of_cubic_bezier_curve(control_pos):
			splitted = split_cublic_beizer_curve( control_pos, partition )
			splitted = asarray( splitted ).astype(int)
# 			debugger()
			
			new_controls.append( controls[ k*3 ] )
			for j, each in enumerate(splitted):
				new_controls += [ Control_point(-1, each[1], False), Control_point(-1, each[2], False) ]
				if j != len(splitted)-1:
					new_controls.append( Control_point(-1, each[-1], True, [4,0]) )	
			
		else:
			new_controls += [ controls[i] for i in range( k*3, k*3+3 ) ]
			
	'''
	if is not closed, add the last control at the end.
	'''
	
	return new_controls


def precompute_all_when_configuration_change( control_pos, skeleton_handle_vertices  ):
	'''
	precompute everything when the configuration changes, in other words, when the number of control points and handles change.
	W_matrices is the table contains all integral result corresponding to each sample point on the boundaries.
	all_weights is an array of num_samples-by-num_handles
	all_vertices is an array of positions of all sampling points. It contains no duplicated points, and matches to all_weights one-on-one
	all_indices is an array of all indices in all_vertices of those sampling points on the boundaries(the curves we need to compute).
	all_pts is an array containing all sampling points and ts for each curve.(boundaries)
	all_dts contains all dts for each curve. It is in the shape of num_curve-by-(num_samples-1)
	'''
	all_pts, all_dts = sample_cubic_bezier_curve_chain( control_pos, 100 )
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
		
		
	W_matrices = zeros( ( len( control_pos ), len( skeleton_handle_vertices ), 4, 4 ) )
	for k in xrange(len( control_pos )):	
		for i in xrange(len( skeleton_handle_vertices )):
			## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
			## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
			W_matrices[k,i] = precompute_W_i_bbw( all_vertices, all_weights, i, all_pts[k][0], all_pts[k][1], all_dts[k])
			
	return W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts

	
		
def main():
	'''
	a console test.
	'''
	control_pos = array([[[100, 300,	1],
		[200, 400,	 1],
		[300, 400,	 1],
		[400, 300,	 1]],

	   [[400, 300,	 1],
		[300, 200,	 1],
		[200, 200,	 1],
		[100, 300,	 1]]])
		
	skeleton_handle_vertices = [[200.0, 300.0, 1.0], [300.0, 300.0, 1.0]] 
	
	W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts = precompute_all_when_configuration_change( control_pos, skeleton_handle_vertices  )
	
	trans = [array([ 1.,  0.,  0.,	0.,	 1.,  0.,  0.,	0.,	 1.]), array([ 1.,	0.,	 0.,  0.,  1., 0., 0., 0., 1.])]	  
						   
	constraints = [[1, 0], [2, 0]]
	
	P_primes, bbw_curves, spline_skin_curves = approximate_beziers(W_matrices, control_pos, skeleton_handle_vertices, trans, constraints, all_weights, all_vertices, all_indices, all_pts, all_dts )
	
	print 'HAHA ~ '
	print P_primes
	
if __name__ == '__main__': main()		
