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
		
		self.transforms = []
		self.handle_positions = []	
		self.precomputed_parameter_table = []
			
	def constraint_change( self, path_index, joint_index, constraint ):
		'''
		change the constraint at a joint of a path.
		path_index tells which path, joint_index tells which joint
		'''
		constraint = list( constraint )
		assert len( constraint ) == 2
		
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
		for i, controls, constraints in zip( range( len( all_controls ) ), all_controls, all_constraints ):
			W_matrices = precomputed_parameters[0][i]
			all_weights = precomputed_parameters[1]
			all_vertices = precomputed_parameters[2]
			all_indices = precomputed_parameters[3][i]
			all_pts = precomputed_parameters[4][i]
			all_dts = precomputed_parameters[5][i]
			
			fast_update = prepare_approximate_beziers( controls, constraints, handles, transforms, W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts )
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

def prepare_approximate_beziers( controls, constraints, handles, transforms, W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts ):
	
	'''
	### 1 construct and solve the linear system for the odd iteration. if the constraints don't contain fixed angle and G1, skip ### 2.
	### 2 If the constraints contain any fixed angle and G1, iterate between the even and odd system until two consecutive solutions are close enough.
	### 3 refine the solutions based on the error of each curve. If it is larger than a threshold, split the curve into two.
	'''
	solutions = None
	controls = concatenate((controls, ones((controls.shape[0],4,1))), axis=2)
	is_closed = array_equal( controls[0,0], controls[-1,-1])
	### 1
	odd = BezierConstraintSolverOdd(W_matrices, controls, constraints, transforms, is_closed )
	last_solutions = solutions = odd.solve()

	smoothness = [ constraint[0] for constraint in constraints ]
	if 'A' in smoothness or 'G1' in smoothness: 
		even = BezierConstraintSolverEven(W_matrices, controls, constraints, transforms, is_closed )	
	
	def update_with_transforms( transforms ):
		odd.update_rhs_for_handles( transforms )
		last_solutions = solutions = odd.solve()
		
		### 2
		smoothness = [ constraint[0] for constraint in constraints ]
		if 'A' in smoothness or 'G1' in smoothness: 
	
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
	
	constraint = constraint = [0, 3, ('A',True) ]
	
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
	
	constraint = [0, 3, ('A',True) ]
	
	return paths_info, skeleton_handle_vertices, constraint 

def get_test_pebble():
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
 {u'bbox_area': 1990.207536521426,
  u'closed': True,
  u'cubic_bezier_chain': [[83.2509994506836, 124.2969970703125],
                          [83.2509994506836, 119.40899658203125],
                          [83.09300231933594, 115.46700286865234],
                          [82.93599700927734, 111.83999633789062],
                          [82.93599700927734, 111.83999633789062],
                          [89.16500091552734, 111.83999633789062],
                          [89.16500091552734, 111.83999633789062],
                          [89.16500091552734, 111.83999633789062],
                          [89.4800033569336, 118.38400268554688],
                          [89.4800033569336, 118.38400268554688],
                          [89.4800033569336, 118.38400268554688],
                          [89.63800048828125, 118.38400268554688],
                          [89.63800048828125, 118.38400268554688],
                          [92.47599792480469, 113.73200225830078],
                          [96.97100067138672, 110.9729995727539],
                          [103.1989974975586, 110.9729995727539],
                          [112.42400360107422, 110.9729995727539],
                          [119.36199951171875, 118.77799987792969],
                          [119.36199951171875, 130.3679962158203],
                          [119.36199951171875, 144.08700561523438],
                          [111.00499725341797, 150.86700439453125],
                          [102.01699829101562, 150.86700439453125],
                          [96.97100067138672, 150.86700439453125],
                          [92.55599975585938, 148.66000366210938],
                          [90.26899719238281, 144.875],
                          [90.26899719238281, 144.875],
                          [90.11100006103516, 144.875],
                          [90.11100006103516, 144.875],
                          [90.11100006103516, 144.875],
                          [90.11100006103516, 165.61000061035156],
                          [90.11100006103516, 165.61000061035156],
                          [90.11100006103516, 165.61000061035156],
                          [83.25199890136719, 165.61000061035156],
                          [83.25199890136719, 165.61000061035156],
                          [83.25199890136719, 165.61000061035156],
                          [83.25199890136719, 124.2969970703125],
                          [83.25199890136719, 124.2969970703125],
                          [83.25199890136719, 124.2969970703125],
                          [83.2509994506836, 124.2969970703125],
                          [83.2509994506836, 124.2969970703125]]},
 {u'bbox_area': 645.6041814603959,
  u'closed': True,
  u'cubic_bezier_chain': [[90.11100006103516, 134.46800231933594],
                          [90.11100006103516, 135.4929962158203],
                          [90.26899719238281, 136.43899536132812],
                          [90.4260025024414, 137.30599975585938],
                          [91.68699645996094, 142.11599731445312],
                          [95.86599731445312, 145.42599487304688],
                          [100.83300018310547, 145.42599487304688],
                          [108.16500091552734, 145.42599487304688],
                          [112.4229965209961, 139.4340057373047],
                          [112.4229965209961, 130.68299865722656],
                          [112.4229965209961, 123.03600311279297],
                          [108.4020004272461, 116.49199676513672],
                          [101.06900024414062, 116.49199676513672],
                          [96.33899688720703, 116.49199676513672],
                          [91.9229965209961, 119.88200378417969],
                          [90.58300018310547, 125.08599853515625],
                          [90.34600067138672, 125.9530029296875],
                          [90.11000061035156, 126.97799682617188],
                          [90.11000061035156, 127.92400360107422],
                          [90.11000061035156, 127.92400360107422],
                          [90.11000061035156, 134.46800231933594],
                          [90.11000061035156, 134.46800231933594],
                          [90.11000061035156, 134.46800231933594],
                          [90.11100006103516, 134.46800231933594],
                          [90.11100006103516, 134.46800231933594]]},
 {u'bbox_area': 1343.071580939577,
  u'closed': True,
  u'cubic_bezier_chain': [[131.97500610351562, 132.1820068359375],
                          [132.13299560546875, 141.56399536132812],
                          [138.125, 145.427001953125],
                          [145.06300354003906, 145.427001953125],
                          [150.02999877929688, 145.427001953125],
                          [153.0260009765625, 144.55999755859375],
                          [155.6280059814453, 143.45599365234375],
                          [155.6280059814453, 143.45599365234375],
                          [156.81100463867188, 148.42300415039062],
                          [156.81100463867188, 148.42300415039062],
                          [154.36700439453125, 149.52699279785156],
                          [150.18800354003906, 150.86700439453125],
                          [144.11700439453125, 150.86700439453125],
                          [132.3699951171875, 150.86700439453125],
                          [125.35299682617188, 143.06100463867188],
                          [125.35299682617188, 131.55099487304688],
                          [125.35299682617188, 120.04100036621094],
                          [132.13299560546875, 110.9729995727539],
                          [143.25, 110.9729995727539],
                          [155.70700073242188, 110.9729995727539],
                          [159.0189971923828, 121.93199920654297],
                          [159.0189971923828, 128.94900512695312],
                          [159.0189971923828, 130.3679962158203],
                          [158.86099243164062, 131.4720001220703],
                          [158.7830047607422, 132.18099975585938],
                          [158.7830047607422, 132.18099975585938],
                          [131.97500610351562, 132.18099975585938],
                          [131.97500610351562, 132.18099975585938],
                          [131.97500610351562, 132.18099975585938],
                          [131.97500610351562, 132.1820068359375],
                          [131.97500610351562, 132.1820068359375]]},
 {u'bbox_area': 229.33571337140165,
  u'closed': True,
  u'cubic_bezier_chain': [[152.3159942626953, 127.21499633789062],
                          [152.39500427246094, 122.79900360107422],
                          [150.5030059814453, 115.94100189208984],
                          [142.69700622558594, 115.94100189208984],
                          [135.67999267578125, 115.94100189208984],
                          [132.60499572753906, 122.40599822998047],
                          [132.05299377441406, 127.21499633789062],
                          [132.05299377441406, 127.21499633789062],
                          [152.3159942626953, 127.21499633789062],
                          [152.3159942626953, 127.21499633789062]]},
 {u'bbox_area': 2075.0125521887094,
  u'closed': True,
  u'cubic_bezier_chain': [[167.61099243164062, 94.02200317382812],
                          [167.61099243164062, 94.02200317382812],
                          [174.47000122070312, 94.02200317382812],
                          [174.47000122070312, 94.02200317382812],
                          [174.47000122070312, 94.02200317382812],
                          [174.47000122070312, 117.98999786376953],
                          [174.47000122070312, 117.98999786376953],
                          [174.47000122070312, 117.98999786376953],
                          [174.6280059814453, 117.98999786376953],
                          [174.6280059814453, 117.98999786376953],
                          [177.07200622558594, 113.73200225830078],
                          [181.48699951171875, 110.9729995727539],
                          [187.63699340820312, 110.9729995727539],
                          [197.09800720214844, 110.9729995727539],
                          [203.7989959716797, 118.85700225830078],
                          [203.72000122070312, 130.44700622558594],
                          [203.72000122070312, 144.08700561523438],
                          [195.12600708007812, 150.86700439453125],
                          [186.61199951171875, 150.86700439453125],
                          [181.09300231933594, 150.86700439453125],
                          [176.67799377441406, 148.73800659179688],
                          [173.83999633789062, 143.69200134277344],
                          [173.83999633789062, 143.69200134277344],
                          [173.60400390625, 143.69200134277344],
                          [173.60400390625, 143.69200134277344],
                          [173.60400390625, 143.69200134277344],
                          [173.28799438476562, 150],
                          [173.28799438476562, 150],
                          [173.28799438476562, 150],
                          [167.29600524902344, 150],
                          [167.29600524902344, 150],
                          [167.45399475097656, 147.3979949951172],
                          [167.61099243164062, 143.53500366210938],
                          [167.61099243164062, 140.14500427246094],
                          [167.61099243164062, 140.14500427246094],
                          [167.61099243164062, 94.02200317382812],
                          [167.61099243164062, 94.02200317382812]]},
 {u'bbox_area': 643.8640927977394,
  u'closed': True,
  u'cubic_bezier_chain': [[174.4709930419922, 134.7050018310547],
                          [174.4709930419922, 135.57200622558594],
                          [174.62899780273438, 136.43899536132812],
                          [174.78599548339844, 137.22799682617188],
                          [176.12600708007812, 142.03799438476562],
                          [180.14700317382812, 145.3489990234375],
                          [185.1929931640625, 145.3489990234375],
                          [192.4459991455078, 145.3489990234375],
                          [196.7830047607422, 139.43499755859375],
                          [196.7830047607422, 130.6840057373047],
                          [196.7830047607422, 123.03700256347656],
                          [192.84100341796875, 116.49299621582031],
                          [185.4290008544922, 116.49299621582031],
                          [180.69900512695312, 116.49299621582031],
                          [176.2830047607422, 119.7249984741211],
                          [174.86399841308594, 125.00800323486328],
                          [174.70599365234375, 125.7969970703125],
                          [174.47000122070312, 126.74299621582031],
                          [174.47000122070312, 127.84600067138672],
                          [174.47000122070312, 127.84600067138672],
                          [174.47000122070312, 134.7050018310547],
                          [174.47000122070312, 134.7050018310547],
                          [174.47000122070312, 134.7050018310547],
                          [174.4709930419922, 134.7050018310547],
                          [174.4709930419922, 134.7050018310547]]},
 {u'bbox_area': 2075.0697996569797,
  u'closed': True,
  u'cubic_bezier_chain': [[212.4720001220703, 94.02200317382812],
                          [212.4720001220703, 94.02200317382812],
                          [219.3300018310547, 94.02200317382812],
                          [219.3300018310547, 94.02200317382812],
                          [219.3300018310547, 94.02200317382812],
                          [219.3300018310547, 117.98999786376953],
                          [219.3300018310547, 117.98999786376953],
                          [219.3300018310547, 117.98999786376953],
                          [219.48800659179688, 117.98999786376953],
                          [219.48800659179688, 117.98999786376953],
                          [221.9320068359375, 113.73200225830078],
                          [226.3470001220703, 110.9729995727539],
                          [232.4969940185547, 110.9729995727539],
                          [241.95799255371094, 110.9729995727539],
                          [248.66000366210938, 118.85700225830078],
                          [248.58099365234375, 130.44700622558594],
                          [248.58099365234375, 144.08700561523438],
                          [239.98699951171875, 150.86700439453125],
                          [231.47300720214844, 150.86700439453125],
                          [225.9530029296875, 150.86700439453125],
                          [221.53799438476562, 148.73800659179688],
                          [218.7010040283203, 143.69200134277344],
                          [218.7010040283203, 143.69200134277344],
                          [218.46400451660156, 143.69200134277344],
                          [218.46400451660156, 143.69200134277344],
                          [218.46400451660156, 143.69200134277344],
                          [218.14700317382812, 150],
                          [218.14700317382812, 150],
                          [218.14700317382812, 150],
                          [212.156005859375, 150],
                          [212.156005859375, 150],
                          [212.31300354003906, 147.3979949951172],
                          [212.4709930419922, 143.53500366210938],
                          [212.4709930419922, 140.14500427246094],
                          [212.4709930419922, 140.14500427246094],
                          [212.4709930419922, 94.02200317382812],
                          [212.4709930419922, 94.02200317382812],
                          [212.4709930419922, 94.02200317382812],
                          [212.4720001220703, 94.02200317382812],
                          [212.4720001220703, 94.02200317382812]]},
 {u'bbox_area': 643.8640927977394,
  u'closed': True,
  u'cubic_bezier_chain': [[219.3300018310547, 134.7050018310547],
                          [219.3300018310547, 135.57200622558594],
                          [219.48800659179688, 136.43899536132812],
                          [219.64500427246094, 137.22799682617188],
                          [220.98599243164062, 142.03799438476562],
                          [225.00599670410156, 145.3489990234375],
                          [230.052001953125, 145.3489990234375],
                          [237.30599975585938, 145.3489990234375],
                          [241.64199829101562, 139.43499755859375],
                          [241.64199829101562, 130.6840057373047],
                          [241.64199829101562, 123.03700256347656],
                          [237.7010040283203, 116.49299621582031],
                          [230.28900146484375, 116.49299621582031],
                          [225.5590057373047, 116.49299621582031],
                          [221.14300537109375, 119.7249984741211],
                          [219.7239990234375, 125.00800323486328],
                          [219.56700134277344, 125.7969970703125],
                          [219.32899475097656, 126.74299621582031],
                          [219.32899475097656, 127.84600067138672],
                          [219.32899475097656, 127.84600067138672],
                          [219.32899475097656, 134.7050018310547],
                          [219.32899475097656, 134.7050018310547],
                          [219.32899475097656, 134.7050018310547],
                          [219.3300018310547, 134.7050018310547],
                          [219.3300018310547, 134.7050018310547]]},
 {u'bbox_area': 388.37468598783016,
  u'closed': True,
  u'cubic_bezier_chain': [[257.3320007324219, 94.02200317382812],
                          [257.3320007324219, 94.02200317382812],
                          [264.2699890136719, 94.02200317382812],
                          [264.2699890136719, 94.02200317382812],
                          [264.2699890136719, 94.02200317382812],
                          [264.2699890136719, 150],
                          [264.2699890136719, 150],
                          [264.2699890136719, 150],
                          [257.3320007324219, 150],
                          [257.3320007324219, 150],
                          [257.3320007324219, 150],
                          [257.3320007324219, 94.02200317382812],
                          [257.3320007324219, 94.02200317382812]]},
 {u'bbox_area': 1343.0320132162888,
  u'closed': True,
  u'cubic_bezier_chain': [[279.802001953125, 132.1820068359375],
                          [279.9590148925781, 141.56399536132812],
                          [285.95098876953125, 145.427001953125],
                          [292.8890075683594, 145.427001953125],
                          [297.85699462890625, 145.427001953125],
                          [300.8529968261719, 144.55999755859375],
                          [303.4540100097656, 143.45599365234375],
                          [303.4540100097656, 143.45599365234375],
                          [304.6369934082031, 148.42300415039062],
                          [304.6369934082031, 148.42300415039062],
                          [302.1929931640625, 149.52699279785156],
                          [298.0140075683594, 150.86700439453125],
                          [291.9440002441406, 150.86700439453125],
                          [280.1960144042969, 150.86700439453125],
                          [273.1789855957031, 143.06100463867188],
                          [273.1789855957031, 131.55099487304688],
                          [273.1789855957031, 120.04100036621094],
                          [279.9590148925781, 110.9729995727539],
                          [291.07598876953125, 110.9729995727539],
                          [303.5329895019531, 110.9729995727539],
                          [306.843994140625, 121.93199920654297],
                          [306.843994140625, 128.94900512695312],
                          [306.843994140625, 130.3679962158203],
                          [306.68701171875, 131.4720001220703],
                          [306.6080017089844, 132.18099975585938],
                          [306.6080017089844, 132.18099975585938],
                          [279.802001953125, 132.18099975585938],
                          [279.802001953125, 132.18099975585938],
                          [279.802001953125, 132.18099975585938],
                          [279.802001953125, 132.1820068359375],
                          [279.802001953125, 132.1820068359375]]},
 {u'bbox_area': 229.33536931639537,
  u'closed': True,
  u'cubic_bezier_chain': [[300.14300537109375, 127.21499633789062],
                          [300.22198486328125, 122.79900360107422],
                          [298.3299865722656, 115.94100189208984],
                          [290.52398681640625, 115.94100189208984],
                          [283.5069885253906, 115.94100189208984],
                          [280.4320068359375, 122.40599822998047],
                          [279.8800048828125, 127.21499633789062],
                          [279.8800048828125, 127.21499633789062],
                          [300.14300537109375, 127.21499633789062],
                          [300.14300537109375, 127.21499633789062]]}]
	handle_positions = [[130,158]]
	constraint = None
	return paths_info, handle_positions, constraint

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
	
	if constraint is not None:
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
