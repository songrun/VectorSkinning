from copy import copy, deepcopy
from bezier_constraint_odd_solver import *
from bezier_constraint_even_solver import *

class Engine:
	'''
	A data persistant that have all information needed to precompute and the system matrix of the previous state.
	'''
	
	def set_control_positions( self, paths_info, boundary_index ):
		'''
		initialize the control points for multiple paths and make the default constraints
		boundary_index tells which path is the outside boundary
		'''
		self.boundary_index = boundary_index
		
		all_controls = [ make_control_points_chain( path[u'cubic_bezier_chain'], path[u'closed'] ) for path in paths_info]
		all_constraints = [ make_constraints_from_control_points( controls, path[u'closed'] ) for controls, path in zip( all_controls, paths_info ) ]
		
		self.all_controls = all_controls
		self.all_constraints = all_constraints	
		self.all_lengths = [ [length_of_cubic_bezier_curve(P) for P in control_points] for control_points in all_controls ]
		
		self.transforms = []
		self.handle_positions = []	
		self.precomputed_parameter_table = []
			
	def constraint_change( self, path_index, joint_index, constraint ):
		'''
		change the constraint at a joint of a path.
		path_index tells which path, joint_index tells which joint
		'''
		constraint = asarray( constraint )
		assert constraint.shape == (2, )
		
		self.all_constraints[ path_index ][ joint_index ] = constraint

	def transform_change( self, i, transform ):
		'''
		change the transform at the index i
		'''
		assert i in range( len( self.transforms ) )
		transform = asarray( transform )
		if len( transform ) == 2:
			transform = concatenate( ( transform, array( [[0, 0, 1]] ) ), axis=0 )
		assert transform.shape == (3,3)
		
		self.transforms[i] = transform
	
	def set_handle_positions( self, handles ):
		'''
		set new handles with identity transforms and keep old handles and transforms unchanged.
		'''
		handles = asarray( handles )
		handle_positions = self.handle_positions
		handle_positions = asarray( handle_positions )
		num_adding = len( handles ) - len( handle_positions )
		
		assert num_adding >= 0
		if len( handle_positions ) != 0:
			assert array_equal( handle_positions, handles[ :len( handle_positions ) ] )
		
		self.handle_positions = handles.tolist()
		
		for i in range( num_adding ):
			self.transforms.append( identity(3) )
	
	def precompute_configuration( self ):
		'''
		precompute W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts
		'''
		handles = self.handle_positions
		all_controls = self.all_controls
		boundary = all_controls[ self.boundary_index ]
		
		layer1 = precompute_all_when_configuration_change( boundary, all_controls, handles )		
		self.precomputed_parameter_table = []
		self.precomputed_parameter_table.append( layer1 )
		
	def solve( self ):
		'''
		send back all groups of controls
		'''
		all_controls = self.all_controls
		all_constraints = self.all_constraints
		all_lengths = self.all_lengths
		
		handles = self.handle_positions
		transforms = self.transforms
		if len( all_controls ) == 0:
			return 'Error Message: No control points'
		elif len( handles ) == 0:
			return 'Error Message: No handles'
		elif len( self.precomputed_parameter_table ) == 0:
			self.precompute_configuration()			
		precomputed_parameters = self.precomputed_parameter_table[0]
		
		result = []
		self.fast_update_functions = []
		for i, controls, constraints, lengths in zip( range( len( all_controls ) ), all_controls, all_constraints, all_lengths ):
			W_matrices = precomputed_parameters[0][i]
			all_weights = precomputed_parameters[1]
			all_vertices = precomputed_parameters[2]
			all_indices = precomputed_parameters[3][i]
			all_pts = precomputed_parameters[4][i]
			all_dts = precomputed_parameters[5][i]
			
			fast_update = prepare_approximate_beziers( controls, constraints, lengths, handles, transforms, W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts )
			self.fast_update_functions.append( fast_update )
			
			result.append(	fast_update( transforms ) )
			
		return result	
		
	
	def solve_transform_change( self ):
		'''
		solve for the new control points when only transform changes
		'''
		result = []
		for fast_update in self.fast_update_functions:
			result.append(	fast_update( self.transforms ) )
			
		return result
	
			
	def compute_tkinter_bbw_affected_curve_per_path( self, all_indices, all_vertices, transforms, all_weights ):
		'''
		compute the bbw curves
		'''
		### 3 
		bbw_curves = []
		for indices in all_indices:
			tps = []	
			for i in indices:
				m = zeros((3,3))
				p = concatenate( ( all_vertices[i], [1.0] ) )
				for h in range(len(transforms)):
					m = m + transforms[h]*all_weights[i,h]
		
				p = dot( m.reshape(3, 3), p.reshape(3,-1) ).reshape(-1)
				tps = tps + [p[0], p[1]]	
			bbw_curves.append(tps)
			
		return bbw_curves

	def compute_tkinter_curve_per_path_solutions( self, solutions ):
		'''
		compute all the points along the curves.
		'''
		spline_skin_curves = []
		for k, solution in enumerate(solutions):
			tps = []
			for t in asarray(all_pts)[k, 1]:
				tbar = asarray([t**3, t**2, t, 1.])
				p = dot(tbar, asarray( M * solution ) )
				tps = tps + [p[0], p[1]]
			spline_skin_curves.append(tps)	
			
		return spline_skin_curves		
		

## The dimensions of a point represented in the homogeneous coordinates
dim = 2

def prepare_approximate_beziers( controls, constraints, lengths, handles, transforms, W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts ):
	
	'''
	### 1 construct and solve the linear system for the odd iteration. if the constraints don't contain fixed angle and G1, skip ### 2.
	### 2 If the constraints contain any fixed angle and G1, iterate between the even and odd system until two consecutive solutions are close enough.
	### 3 refine the solutions based on the error of each curve. If it is larger than a threshold, split the curve into two.
	'''
	solutions = None
	controls = concatenate((controls, ones((controls.shape[0],4,1))), axis=2)
	is_closed = array_equal( controls[0,0], controls[-1,-1])
	### 1
	odd = BezierConstraintSolverOdd(W_matrices, controls, constraints, lengths, transforms, is_closed )
#	odd.update_rhs_for_handles( transforms )
	last_solutions = solutions = odd.solve()

	if 2 in constraints[:,0] or 4 in constraints[:,0]: 
		even = BezierConstraintSolverEven(W_matrices, controls, constraints, lengths, transforms, is_closed )	
	
	def update_with_transforms( transforms ):
		odd.update_rhs_for_handles( transforms )
		last_solutions = solutions = odd.solve()
		
		### 2	
		if 2 in constraints[:,0] or 4 in constraints[:,0]: 
	
			even.update_rhs_for_handles( transforms )
	
			for iter in xrange( 10 ):
				even.update_system_with_result_of_previous_iteration( solutions )
				last_solutions = solutions
				solutions = even.solve()
			
				if allclose(last_solutions, solutions, atol=1.0, rtol=1e-03):
					break
			
				## Check if error is low enough and terminate
				odd.update_system_with_result_of_previous_iteration( solutions )
				last_solutions = solutions
				solutions = odd.solve()
			
				if allclose(last_solutions, solutions, atol=1.0, rtol=1e-03):
					break
		
		return solutions
	
	return update_with_transforms					
	
	
	### 3
#	new_controls = adapt_configuration_based_on_diffs( controls, bbw_curves, spline_skin_curves, all_dts )
#	
#	if enable_refinement and len( new_controls ) > len( controls ): 
# #			debugger()
#		new_control_pos = get_controls( new_controls )[0]
# 
#		W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts = precompute_all_when_configuration_change( new_control_pos, handles  )
#	
#		solutions, bbw_curves, spline_skin_curves = approximate_beziers(W_matrices, new_controls, handles, transforms, all_weights, all_vertices, all_indices, all_pts, all_dts, False) 



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
#			debugger()
			
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


def precompute_all_when_configuration_change( controls_on_boundary, all_control_positions, skeleton_handle_vertices	 ):
	'''
	precompute everything when the configuration changes, in other words, when the number of control points and handles change.
	W_matrices is the table contains all integral result corresponding to each sample point on the boundaries.
	all_weights is an array of num_samples-by-num_handles
	all_vertices is an array of positions of all sampling points. It contains no duplicated points, and matches to all_weights one-on-one
	all_indices is an array of all indices in all_vertices of those sampling points on the boundaries(the curves we need to compute).
	all_pts is an array containing all sampling points and ts for each curve.(boundaries)
	all_dts contains all dts for each curve. It is in the shape of num_curve-by-(num_samples-1)
	'''
	num_samples = 100
	all_pts = []
	all_dts = []
	for control_pos in all_control_positions:
		pts, dts = sample_cubic_bezier_curve_chain( control_pos, num_samples ) 
		all_pts.append( pts )
		all_dts.append( dts )	
	
	boundary_pts, boundary_dts = sample_cubic_bezier_curve_chain( controls_on_boundary, num_samples )
	
	boundary_pos = [ curve[0] for curve in boundary_pts ]
	all_pos = [ [ pts[0] for pts in curve_pts ] for curve_pts in all_pts ]
	
	all_vertices, all_weights, all_indices= triangulate_and_compute_weights( boundary_pos, skeleton_handle_vertices, all_pos )
	
	print 'Precomputing W_i...'
	W_matrices = []
	for j, control_pos in enumerate( all_control_positions ):
		W_matrices.append( zeros( ( len( control_pos ), len( skeleton_handle_vertices ), 4, 4 ) ) )		
		for k in xrange(len( control_pos )):	
			for i in xrange(len( skeleton_handle_vertices )):
				## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
				## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
				W_matrices[j][k,i] = precompute_W_i_bbw( all_vertices, all_weights, i, all_indices[j][k], all_pts[j][k][0], all_pts[j][k][1], all_dts[j][k])
				
	W_matrices = asarray( W_matrices )
	print '...finished.'

	return [W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts]

	

def get_test1():
	
	paths_info =  [
	{u'bbox_area': 81583.4111926838,
  u'closed': True,
  u'cubic_bezier_chain': [[46.95399856567383, 114.95899963378906],
						  [35.944000244140625, 177.95700073242188],
						  [96.1259994506836, 266.40399169921875],
						  [198.39999389648438, 266.40399169921875],
						  [300.67401123046875, 266.40399169921875],
						  [342.614990234375, 182.7259979248047],
						  [342.614990234375, 122.19000244140625],
						  [342.614990234375, 61.65399932861328],
						  [366.375, 19.503999710083008],
						  [241.58200073242188, 21.156999588012695],
						  [116.78900146484375, 22.809999465942383],
						  [61.83000183105469, 29.834999084472656],
						  [46.95399856567383, 114.95899963378906],
						  [46.95399856567383, 114.95899963378906],
						  [46.95399856567383, 114.95899963378906],
						  [46.95399856567383, 114.95899963378906]]}
#	{u'bbox_area': 55.29089948625202,
#	u'closed': False,
#	u'cubic_bezier_chain': [[-255.1510009765625, 5.1479997634887695],
#							[-255.76300048828125, 9.116000175476074],
#							[-263.0260009765625, 8.20199966430664],
#							[-263.8190002441406, 5.1479997634887695],
#							[-263.51300048828125, -0.24000000953674316],
#							[-255.78399658203125, 0.5950000286102295],
#							[-255.1510009765625, 5.1479997634887695],
#							[-260.4859924316406, 5.1479997634887695],
#							[-259.3039855957031, 4.995999813079834],
#							[-257.14300537109375, 5.821000099182129],
#							[-257.8190002441406, 3.815000057220459],
#							[-259.3370056152344, 3.628000020980835],
#							[-260.32598876953125, 3.9749999046325684],
#							[-260.4859924316406, 5.1479997634887695]]},						   
#							
#  {u'bbox_area': 4.065760665711228,
#	u'closed': True,
#	u'cubic_bezier_chain': [[155.34100341796875, 86.31900024414062],
#							[156.80299377441406, 86.19999694824219],
#							[157.47000122070312, 86.86699676513672],
#							[157.34300231933594, 88.32099914550781],
#							[157.34300231933594, 88.32099914550781],
#							[155.34100341796875, 88.32099914550781],
#							[155.34100341796875, 88.32099914550781],
#							[155.34100341796875, 88.32099914550781],
#							[155.34100341796875, 86.31900024414062],
#							[155.34100341796875, 86.31900024414062]]},							
#  {u'bbox_area': 6.86434952491282,
#	u'closed': False,
#	u'cubic_bezier_chain': [[-272.48699951171875, -4.85099983215332],
#							[-270.177001953125, -5.317999839782715],
#							[-270.0920104980469, -2.513000011444092],
#							[-271.1549987792969, -1.5190000534057617],
#							[-272.614990234375, -1.61899995803833],
#							[-272.5870056152344, -3.197000026702881],
#							[-272.48699951171875, -4.85099983215332]]}
					]
	
	skeleton_handle_vertices = [[176, 126]] 
#	skeleton_handle_vertices = [[200.0, 300.0, 1.0], [300.0, 300.0, 1.0]] 
	
	constraint = [0, 3, (2,1) ]
	
	return paths_info, skeleton_handle_vertices, constraint

def get_test2():
	paths_info = [{u'bbox_area': 81583.4111926838,
  u'closed': True,
  u'cubic_bezier_chain': [[46.95399856567383, 114.95899963378906],
						  [35.944000244140625, 177.95700073242188],
						  [96.1259994506836, 266.40399169921875],
						  [198.39999389648438, 266.40399169921875],
						  [300.67401123046875, 266.40399169921875],
						  [342.614990234375, 182.7259979248047],
						  [342.614990234375, 122.19000244140625],
						  [342.614990234375, 61.65399932861328],
						  [366.375, 19.503999710083008],
						  [241.58200073242188, 21.156999588012695],
						  [116.78900146484375, 22.809999465942383],
						  [61.83000183105469, 29.834999084472656],
						  [46.95399856567383, 114.95899963378906]]},
 {u'bbox_area': 15.526111421524547,
  u'closed': False,
  u'cubic_bezier_chain': [[134.5, 128.5],
						  [132.61599731445312, 127.67900085449219],
						  [130.8800048828125, 126.66699981689453],
						  [129.3159942626953, 125.50499725341797]]},
 {u'bbox_area': 9399.832030713733,
  u'closed': False,
  u'cubic_bezier_chain': [[121.60099792480469, 115.66500091552734],
						  [115.31800079345703, 99.75399780273438],
						  [130.28199768066406, 78.03600311279297],
						  [188.5, 85.5],
						  [256.49200439453125, 94.21700286865234],
						  [272.8139953613281, 111.29199981689453],
						  [272.5719909667969, 137.718994140625]]},
 {u'bbox_area': 4.357566170394421,
  u'closed': False,
  u'cubic_bezier_chain': [[272.23199462890625, 144.0469970703125],
						  [272.04998779296875, 145.98500061035156],
						  [271.802001953125, 147.968994140625],
						  [271.5, 150]]}]
	
	skeleton_handle_vertices = [[176, 126]] 
#	skeleton_handle_vertices = [[200.0, 300.0, 1.0], [300.0, 300.0, 1.0]] 
	
	constraint = [0, 3, (2,1) ]
	
	return paths_info, skeleton_handle_vertices, constraint 
	
def main():
	'''
	a console test.
	'''
	
	# paths_info, skeleton_handle_vertices, constraint = get_test1()
	paths_info, skeleton_handle_vertices, constraint = get_test2()
	
	engine = Engine()
	boundary_path = max(paths_info, key=lambda e : e[u'bbox_area']) 
	boundary_index = paths_info.index( boundary_path )
	
	engine.set_control_positions( paths_info, boundary_index )
	
	engine.constraint_change( constraint[0], constraint[1], constraint[2] )

	engine.set_handle_positions( skeleton_handle_vertices )
	
#	engine.set_transforms()
	engine.precompute_configuration()
	all_paths = engine.solve()
	
	for path in all_paths:
		if len( path ) > 1:
			chain = concatenate( asarray(path)[:-1, :-1] )
			chain = concatenate( ( chain, path[-1] ) )
		else:
			chain = path[0]
		print chain
#	debugger()
#	parameters = precompute_all_when_configuration_change( control_pos, skeleton_handle_vertices  )
#	
#	trans = [array([ 1.,  0.,  0.,	0.,	 1.,  0.,  0.,	0.,	 1.]), array([ 1.,	0.,	 0.,  0.,  1., 0., 0., 0., 1.])]	  
#						   
	
	print 'HAHA ~ '
	
if __name__ == '__main__': main()		
