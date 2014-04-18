from copy import copy, deepcopy
from bezier_constraint_odd_solver import *
from bezier_constraint_even_solver import *

from tictoc import tic, toc

class EngineError( Exception ): pass
class NoControlPointsError( EngineError ): pass
class NoHandlesError( EngineError ): pass

systems = []
############################### Basic Engine #################################
class Engine:
	def init_engine( self, paths_info, boundary_index ):
		'''
		initialize the control points for multiple paths and make the default constraints
		boundary_index tells which path is the outside boundary
		'''
		self.boundary_index = boundary_index	
		self.all_controls = [ make_control_points_chain( path[u'cubic_bezier_chain'], path[u'closed'] ) for path in paths_info]
		self.all_constraints = [ make_constraints_from_control_points( controls, path[u'closed'] ) for controls, path in zip( self.all_controls, paths_info ) ]
		self.num_of_paths = len( self.all_controls )
		
		self.transforms = []
		self.handle_positions = []
		self.weight_function = 'shepard'
	
	def copy_engine( self, engine ):
		self.boundary_index = engine.boundary_index
		self.all_controls = engine.all_controls
		self.all_constraints = engine.all_constraints
		self.num_of_paths = engine.num_of_paths
		
		self.transforms = engine.transforms
		self.handle_positions = engine.handle_positions
		self.weight_function = engine.weight_function
	
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
			
	def set_handle_positions( self, new_handle_positions, new_transforms = None ):
		'''
		Replaces the set of known handles with the positions and transforms
		in new_handle_positions, new_transforms.
		If 'new_transforms' is not specified or is None,
		all handles will be set to the identity transformation.
		'''
		
		assert new_transforms is None or len( new_transforms ) == len( new_handle_positions )
		
		## Set the new handle positions.
		self.handle_positions = asarray( new_handle_positions ).tolist()
		## Initialize all the transforms to the identity matrix.
		## NOTE: We don't directly use 'new_transforms' here, because they're
		##		 in a different format.
		self.transforms = [ identity(3) for i in range(len( new_handle_positions )) ]
		
		## Set the corresponding transformations.
		if new_transforms is not None:
			for i, transform in enumerate( new_transforms ):
				self.transform_change( i, transform )

	def set_weight_function( self, weight_function ):
		'''
		set weight_function
		'''
		if self.weight_function != weight_function:
			self.weight_function = weight_function
			self.precompute_configuration( ) 
	
	def get_weight_function( self ):
		'''
		gets weight_function
		'''
		return self.weight_function

	def get_engine_type( self ):
		'''
		gets engine_type
		'''
		return self.engine_type

	def compute_target_path( self, path_indices, all_vertices, transforms, all_weights ):
		'''
		compute the target curves for just one path.
		'''
		### 3 
		target_curves_per_path = []
		for indices in path_indices:
				
			## http://jameshensman.wordpress.com/2010/06/14/multiple-matrix-multiplication-in-numpy/
			As = dot( asarray(transforms).T, all_weights[ indices ].T ).T
			Bs = append( all_vertices[ indices ], ones( ( len( indices ), 1 ) ), axis = 1 )
			tps = sum(As*Bs[:,newaxis,:],-1)[:,:2]
			
			target_curves_per_path.append(tps)
		
		return target_curves_per_path

	def compute_deformed_path( self, solutions, all_ts ):
		'''
		compute all the points along the curves for just one path.
		'''
		deformed_curves_per_path = []
		for k, solution in enumerate(solutions):

			tps = dot( array([all_ts[k]**3, all_ts[k]**2, all_ts[k], ones(all_ts[k].shape)]).T, asarray( dot( M, solution ) ) )
			deformed_curves_per_path.append(tps)
			
		return deformed_curves_per_path	
	
	def constraint_change( self, path_index, joint_index, constraint ): pass
	def precompute_configuration( self ): pass
	def set_enable_arc_length( self, is_arc_enabled ): pass
	def get_enable_arc_length( self ): return False
	def set_iterations( self, whether ): pass
	def prepare_to_solve( self ): 
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def solve_transform_change( self ): 
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def compute_energy_and_maximum_distance( self ):
		'''
		compute the error between the skinning spline and the bbw_affected position.
		'''	
		all_controls, handle_positions, transforms = self.all_controls, self.handle_positions, self.transforms
		solutions = self.solutions
		boundary_index = self.boundary_index
		weight_function = self.weight_function
		
		num_samples = 100
		all_pts, all_ts, all_dts = [], [], []
		for control_pos in all_controls:
			pts, ts, dts = sample_cubic_bezier_curve_chain( control_pos, num_samples )
			all_pts.append( pts )
			all_ts.append( ts )
			all_dts.append( dts )
			
		all_vertices, all_weights, all_indices = compute_all_weights( all_pts, handle_positions, boundary_index, weight_function )
		
		all_pts, all_ts, all_dts = asarray( all_pts ), asarray( all_ts ), asarray( all_dts )
		
		energy, target_paths, distances = [], [], []
		for i in range( self.num_of_paths ):

			path_indices = all_indices[i]		
			path_pts, path_ts, path_dts = all_pts[i], all_ts[i], all_dts[i]
			path_lengths = asarray( [ map( mag, ( curve_pts[1:] - curve_pts[:-1] ) ) for curve_pts in path_pts ]).sum( axis=1 )
			
			target_path = self.compute_target_path( path_indices, all_vertices, transforms, all_weights )
			target_paths.append( target_path )
			
			deformed_path = self.compute_deformed_path( solutions[i], path_ts )			
 			energy.append( compute_error_metric( target_path, deformed_path, path_dts, path_lengths ) )
 			distances.append( compute_maximum_distances( target_path, deformed_path ) )
		
 		return energy, target_paths, distances

			
############################### Basic Engine End #################################

class FourControlsEngine(Engine) :
	'''
	- apply deformation to all control points
	'''
	
	def prepare_to_solve( self ): 
		all_controls, handle_positions, transforms = self.all_controls, self.handle_positions, self.transforms
		boundary_index, weight_function = self.boundary_index, self.weight_function

  		all_vertices, all_weights, all_indices = compute_all_weights_shepard( all_controls, handle_positions )
		if 'bbw' == weight_function:
 			
			all_endpoints = [ asarray([ [curve[0], curve[3]] for curve in path ]) for path in all_controls ] 			
 			vs, faces, boundary_edges, all_weights, all_maps = compute_all_weights_bbw( all_endpoints, handle_positions, boundary_index, customized = True )
			
			flatten_ctrls = flatten_paths( all_controls )[0]
			ctrl_vertices, ctrl_weights, ctrl_maps = barycentric_projection( vs, faces, boundary_edges, all_weights, flatten_ctrls )
			
			all_weights = asarray([ ctrl_weights[index] for index in ctrl_maps])
		
		elif 'mvc' == weight_function:
			all_vertices, all_weights, all_indices = compute_all_weights_mvc( all_controls, handle_positions )
	
			
		self.all_vertices, self.all_weights, self.all_indices = all_vertices, all_weights, all_indices
	
	def solve_transform_change( self ):
		
		transforms = self.transforms
		all_vertices, all_weights, all_indices = self.all_vertices, self.all_weights, self.all_indices
		result = []
		for path_indices in all_indices:
			path_controls = []
			for indices in path_indices:

				As = dot( asarray(transforms).T, all_weights[ indices ].T ).T
				Bs = append( all_vertices[ indices ], ones( ( len( indices ), 1 ) ), axis = 1 )
				#tps = sum(As*Bs[:,newaxis,:],-1)[:,:2]
				tps = asarray([ dot( A, B ) for A, B in zip( As, Bs ) ])
		
				path_controls.append(tps)
			result.append( path_controls )
		
		result = [ asarray(path)[:,:,:2] for path in result ]
			
		self.solutions = result
		return result
	
class TwoEndpointsEngine(Engine):
	'''
	- apply endpoint deformation to closest interior points
	return new control points
	'''
	def prepare_to_solve( self ):
		all_controls, handle_positions, transforms = self.all_controls, self.handle_positions, self.transforms
		boundary_index, weight_function = self.boundary_index, self.weight_function
	
  		all_vertices, all_weights, all_indices = compute_all_weights_shepard( all_controls, handle_positions )
  		all_weights[ 1: : 4] = all_weights[ 0: : 4]
  		all_weights[ 2: : 4] = all_weights[ 3: : 4]
  		
 		if 'bbw' == weight_function:
 			
			all_endpoints = [ asarray([ [curve[0], curve[3]] for curve in path ]) for path in all_controls ] 			
 			vs, faces, boundary_edges, all_weights, all_maps = compute_all_weights_bbw( all_endpoints, handle_positions, boundary_index, customized = True )
	
		
			flatten_maps = concatenate( [ concatenate( path_map ) for path_map in all_maps ] )
			end_weights = asarray( [all_weights[index] for index in flatten_maps] )
			
			flatten_ctrls = flatten_paths( all_controls )[0]
			ctrl_vertices, ctrl_weights, ctrl_maps = barycentric_projection( vs, faces, boundary_edges, all_weights, flatten_ctrls )
			
			all_weights = asarray([ ctrl_weights[index] for index in ctrl_maps])
				
			all_weights[ 0: : 4] = all_weights[ 1: : 4] = end_weights[ 0 : : 2 ]
			all_weights[ 2: : 4] = all_weights[ 3: : 4] = end_weights[ 1 : : 2 ]
	
		elif 'mvc' == weight_function:
			all_vertices, all_weights, all_indices = compute_all_weights_mvc( all_controls, handle_positions )
			all_weights[ 1: : 4] = all_weights[ 0: : 4]
  			all_weights[ 2: : 4] = all_weights[ 3: : 4]
			
		self.all_vertices, self.all_weights, self.all_indices = all_vertices, all_weights, all_indices
			
	def solve_transform_change( self ):

		transforms = self.transforms
		all_vertices, all_weights, all_indices = self.all_vertices, self.all_weights, self.all_indices		
		result = []
		for path_indices in all_indices:
			path_controls = []
			for indices in path_indices:

				As = dot( asarray(transforms).T, all_weights[ indices ].T ).T
				Bs = append( all_vertices[ indices ], ones( ( len( indices ), 1 ) ), axis = 1 )
				#tps = sum(As*Bs[:,newaxis,:],-1)[:,:2]
				tps = asarray([ dot( A, B ) for A, B in zip( As, Bs ) ])
		
				path_controls.append(tps)
			result.append( path_controls )
		
		result = [ asarray(path)[:,:,:2] for path in result ]
		
		self.solutions = result	
		return result


class JacobianEngine(Engine):
	'''
	apply jacobian of endpoint deformation to vector towards closest interior points 
	(n.b. not the inverse transpose of it; these are tangents, not normals)
	return new control points 
	'''
	def prepare_to_solve( self ):
		pass

	def solve_transform_change( self ):
		all_controls, handle_positions, transforms = self.all_controls, self.handle_positions, self.transforms

		all_vertices, all_weights, all_indices = compute_all_weights_shepard( all_controls, handle_positions )
		
		def shepard_jacobian( handle_positions, transforms, T_p, p ):
			'''
			return a 2-by-2 jacobian matrix
			'''
			jac = zeros( (2,2) )
			 
			def derivative_ws( hs ): 
				'''
				compute an array of derivative of w_i with respect to x or y, k = 0 means x, k = 1 means y
				'''
				hs = asarray( hs )
				assert len( hs.shape ) == 2			
			 	assert hs.shape[1] == p.shape[0]
			 	
			 	## Compute the distance squared from each handle position.
				diffs = p - hs
				diffs = diffs**2
				diffs = diffs.sum( axis = 1 )
				
				eps = 1e-7
				## where() gives us the indices where the condition is true.
				wh = where( abs( diffs ) < eps )
				assert len( wh ) == 1
				if len( wh[0] ) > 0: return zeros( len( hs ) )

				us = asarray( 1. / diffs )
				dus = asarray ([-2.*sep_diff/square_sum**2 for sep_diff, square_sum in zip( (p-hs), diffs )])
				
				derivatives = [ dus[i] / us.sum() - us[i]*dus.sum(axis=0) / us.sum()**2 for i in range( len(hs) ) ]
				
				return derivatives
				
			derivatives = asarray( derivative_ws( handle_positions ) )
			transforms = asarray( transforms )

			T_x = asarray([d*transform for d, transform in zip(derivatives[:, 0], transforms)]).sum(axis=0).squeeze()
			T_y = asarray([d*transform for d, transform in zip(derivatives[:, 1], transforms)]).sum(axis=0).squeeze()
			
			jac[0,0] = T_p[0,0] + T_x[0,0]*p[0] + T_x[0,1]*p[1] + T_x[0,2]
 			jac[0,1] = T_p[0,1] + T_y[0,0]*p[0] + T_y[0,1]*p[1] + T_y[0,2]
			jac[1,0] = T_p[1,0] + T_x[1,0]*p[0] + T_x[1,1]*p[1] + T_x[1,2]
 			jac[1,1] = T_p[1,1] + T_y[1,0]*p[0] + T_y[1,1]*p[1] + T_y[1,2]
 			
 			return jac
 			
 			
		result = []
		for path_indices in all_indices:
			path_controls = []
			for indices in path_indices:
				assert len( indices ) == 4

				As = dot( asarray(transforms).T, all_weights[ indices ].T ).T
				Bs = append( all_vertices[ indices ], ones( ( len( indices ), 1 ) ), axis = 1 )
				
				tps = asarray([ dot( A, B ) for A, B in zip( As, Bs ) ])
				## replace the two interior control points' position to endpoint + Jacobian * derivative vector
				vectors = asarray( [all_vertices[indices[1]] - all_vertices[indices[0]],
									all_vertices[indices[2]] - all_vertices[indices[3]] ] ) 
				
				tps[1,:2] = tps[0,:2] + dot( shepard_jacobian( handle_positions, transforms, As[0], all_vertices[indices[0]] ), vectors[0] )
				tps[2,:2] = tps[3,:2] + dot( shepard_jacobian( handle_positions, transforms, As[3], all_vertices[indices[3]] ), vectors[1] )

				path_controls.append(tps)
			result.append( path_controls )
		
		result = [ asarray(path)[:,:,:2] for path in result ]
		
		self.solutions = result	
		return result
		

	
class YSEngine(Engine):
	'''
	A data persistant that have all information needed to precompute and the system matrix of the previous state.
	'''
	
	def init_engine( self, paths_info, boundary_index ):
		'''
		initialize the control points for multiple paths and make the default constraints
		boundary_index tells which path is the outside boundary
		'''
		Engine.init_engine( self, paths_info, boundary_index )
	
		self.weight_function = 'bbw'
		self.is_arc_enabled = parameters.kArcLengthDefault
		self.perform_multiple_iterations = True
		self.precomputed_parameter_table = []

	def copy_engine( self, engine ):
		Engine.copy_engine( self, engine )

		self.is_arc_enabled = parameters.kArcLengthDefault
		self.perform_multiple_iterations = True
		self.precomputed_parameter_table = []
		 		
	def constraint_change( self, path_index, joint_index, constraint ):
		'''
		change the constraint at a joint of a path.
		path_index tells which path, joint_index tells which joint
		'''
		constraint = list( constraint )
		assert len( constraint ) == 2
		
		self.all_constraints[ path_index ][ joint_index ] = constraint
	
	def precompute_configuration( self ):
		'''
		precompute W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts
		'''
		handles = self.handle_positions
		all_controls = self.all_controls
		
		if len( handles ) == 0: return
		
		weight_function = self.weight_function
		is_arc_enabled = self.is_arc_enabled
		layer1 = precompute_all_when_configuration_change( self.boundary_index, all_controls, handles, weight_function, is_arc_enabled )
		self.precomputed_parameter_table = []
		self.precomputed_parameter_table.append( layer1 )
		
	def prepare_to_solve( self ):
		'''
		call this and then call solve_transform_change() to get back all groups of controls
		'''
		if len( self.all_controls ) == 0:
			raise NoControlPointsError()
		elif len( self.handle_positions ) == 0:
			raise NoHandlesError()
		elif len( self.precomputed_parameter_table ) == 0:
			self.precompute_configuration()			
		
		all_controls = self.all_controls
		all_constraints = self.all_constraints
	
		handles = self.handle_positions
		transforms = self.transforms
		precomputed_parameters = self.precomputed_parameter_table[0]
		
		is_arc_enabled = self.is_arc_enabled
		
		tic( 'Generating system matrices...' )
		self.fast_update_functions = []
		for i, controls, constraints in zip( range( len( all_controls ) ), all_controls, all_constraints ):
			W_matrices = precomputed_parameters.W_matrices[i]
			ts = precomputed_parameters.all_ts[i]
			dts = precomputed_parameters.all_dts[i]
			lengths = precomputed_parameters.all_lengths[i]
			
			fast_update = prepare_approximate_beziers( controls, constraints, handles, transforms, lengths, W_matrices, ts, dts, is_arc_enabled )
			self.fast_update_functions.append( fast_update )
		toc()
		
	
	def solve_transform_change( self ):
		'''
		solve for the new control points when only transform changes
		'''
		result = []
		for fast_update in self.fast_update_functions:
			result.append(	fast_update( self.transforms, self.perform_multiple_iterations ) )
		
		self.solutions = result

		return result
	
	def set_enable_arc_length( self, is_arc_enabled ):
		'''
		set is_arc_enabled flag on/off
		'''
		if self.is_arc_enabled != is_arc_enabled:
			self.is_arc_enabled = is_arc_enabled
			self.precompute_configuration() 
	
	def get_enable_arc_length( self ):
		'''
		get is_arc_enabled flag
		'''
		return self.is_arc_enabled
	
	def set_iterations( self, whether ):
		self.perform_multiple_iterations = whether		
			
	def compute_energy_and_maximum_distance( self ):
		'''
		compute the error between the skinning spline and the bbw_affected position.
		'''
		if len( self.precomputed_parameter_table ) == 0:
			raise EngineError( "compute_energy() can't be called before solve()" )
		
		all_controls = self.all_controls
		transforms = self.transforms
		precomputed_parameters = self.precomputed_parameter_table[0]
		solutions = self.solutions
		
		all_weights = precomputed_parameters.all_weights
		all_vertices = precomputed_parameters.all_vertices
			
		energy = []
		target_paths = []
		distances = []
		for i in range( self.num_of_paths ):

			path_indices = precomputed_parameters.all_indices[i]		
			path_pts = precomputed_parameters.all_pts[i]
			path_dts = precomputed_parameters.all_dts[i]
			path_ts = precomputed_parameters.all_ts[i]
			path_lengths = precomputed_parameters.all_lengths[i]
			
			target_path = self.compute_target_path( path_indices, all_vertices, transforms, all_weights )

			deformed_path = self.compute_deformed_path( solutions[i], path_ts )
			
			energy.append( compute_error_metric( target_path, deformed_path, path_dts, path_lengths ) )
			
			target_paths.append( target_path )

			distances.append( compute_maximum_distances( target_path, deformed_path ) )
		
		return energy, target_paths, distances
			
				
## The dimensions of a point represented in the homogeneous coordinates
# dim = 2

def prepare_approximate_beziers( controls, constraints, handles, transforms, lengths, W_matrices, ts, dts, kArcLength=False ):
	
	'''
	### 1 construct and solve the linear system for the odd iteration. if the constraints don't contain fixed angle and G1, skip ### 2.
	### 2 If the constraints contain any fixed angle and G1, iterate between the even and odd system until two consecutive solutions are close enough.
	### 3 refine the solutions based on the error of each curve. If it is larger than a threshold, split the curve into two.
	'''
	
	kPickleDebug = False
	if kPickleDebug:
		import cPickle as pickle
		debug_out = 'solutions-arc' + str(kArcLength) + '.pickle'
		all_solutions = []
	
	solutions = None
	controls = concatenate((controls, ones((controls.shape[0],4,1))), axis=2)
	is_closed = array_equal( controls[0,0], controls[-1,-1])
	### 1
	odd = BezierConstraintSolverOdd(W_matrices, controls, constraints, transforms, lengths, ts, dts, is_closed, kArcLength )

	smoothness = [ constraint[0] for constraint in constraints ]
	if 'A' in smoothness or 'G1' in smoothness: 
		even = BezierConstraintSolverEven(W_matrices, controls, constraints, transforms, lengths, ts, dts, is_closed, kArcLength )
	
	def update_with_transforms( transforms, multiple_iterations = True ):
		iteration = 1
		odd.update_rhs_for_handles( transforms )
		last_odd_solutions = solutions = odd.solve()
		if not multiple_iterations: return solutions
		## Keep track of even and odd solutions separately, because they may converge
		## separately.
		last_even_solutions = None
		
		systems.append( odd.system )
		
		if kPickleDebug:
			all_solutions.append( solutions )
			pickle.dump( all_solutions, open( debug_out, "wb" ) )
		
		### 2
		smoothness = [ constraint[0] for constraint in constraints ]
		if 'A' in smoothness or 'G1' in smoothness: 
			## TODO Q: Why does is sometimes seem like this code only runs if there
			##         a print statement inside? It seems haunted.
			even.update_rhs_for_handles( transforms )
			
			for i in xrange( 10 ):
				iteration += 1
				even.update_system_with_result_of_previous_iteration( solutions )
				solutions = even.solve()
				
				if kPickleDebug:
					all_solutions.append( solutions )
					pickle.dump( all_solutions, open( debug_out, "wb" ) )
				
				#print 'max |last solutions - solutions|:', abs( asarray( last_solutions ) - asarray( solutions ) ).max()
				#from pprint import pprint
				#pprint( solutions )
				if allclose(last_odd_solutions, solutions, atol=1.0, rtol=1e-03):
				    break
				if last_even_solutions is not None and allclose(last_even_solutions, solutions, atol=1.0, rtol=1e-03):
				    break
				
				last_even_solutions = solutions
				
				## For debugging, randomly don't perform the last even iteration.
				#import random
				#if i == 9 and random.randint(0,1): break
				
				## Check if error is low enough and terminate
				iteration += 1
				odd.update_system_with_result_of_previous_iteration( solutions )
				solutions = odd.solve()
				
				#print 'max |last solutions - solutions|:', abs( asarray( last_solutions ) - asarray( solutions ) ).max()
				#pprint( solutions )
				if allclose(last_even_solutions, solutions, atol=1.0, rtol=1e-03):
					break
				if allclose(last_odd_solutions, solutions, atol=1.0, rtol=1e-03):
					break
				
				last_odd_solutions = solutions
		
		print 'iterations:', iteration
		return solutions
	
	return update_with_transforms					
	

def precompute_all_when_configuration_change( boundary_index, all_control_positions, skeleton_handle_vertices, weight_function = 'bbw', kArcLength=False ):
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
	all_ts = []
	all_lengths = []
	for control_pos in all_control_positions:
		path_pts, path_ts, path_dts = sample_cubic_bezier_curve_chain( control_pos, num_samples )
		all_pts.append( path_pts )
		all_ts.append( path_ts )
		
		## Compute all_lengths
		path_dss = [ map( mag, ( curve_pts[1:] - curve_pts[:-1] ) ) for curve_pts in path_pts ]
		path_dss = asarray( path_dss )
		path_lengths = [ sum( path_dss[i] ) for i in range( len( path_dss ) ) ]
		all_lengths.append( path_lengths )
		## Then normalize dss
		dss = [ ds / length for ds, length in zip( path_dss, path_lengths ) ]
		
		
		if kArcLength:
			all_dts.append( path_dss )
		else:
			all_dts.append( path_dts )

	
	all_vertices, all_weights, all_indices = compute_all_weights( all_pts, skeleton_handle_vertices, boundary_index, weight_function )
	
	tic( 'Precomputing W_i...' )
	W_matrices = []
	for j, control_pos in enumerate( all_control_positions ):
		W_matrices.append( zeros( ( len( control_pos ), len( skeleton_handle_vertices ), 4, 4 ) ) )		
		for k in xrange(len( control_pos )):	
			for i in xrange(len( skeleton_handle_vertices )):
				## indices k, i, 0 is integral of w*tbar*tbar.T, used for C0, C1, G1,
				## indices k, i, 1 is integral of w*tbar*(M*tbar), used for G1
				W_matrices[j][k,i] = precompute_W_i( all_vertices, all_weights, i, all_indices[j][k], all_pts[j][k], all_ts[j][k], all_dts[j][k])
				
	W_matrices = asarray( W_matrices )
	toc()
	
	class Layer( object ): pass
	layer = Layer()
	layer.W_matrices = W_matrices
	layer.all_weights = all_weights
	layer.all_vertices = all_vertices
	layer.all_indices = all_indices
	layer.all_pts = all_pts
	layer.all_dts = all_dts
	layer.all_ts = all_ts
	layer.all_lengths = all_lengths
	return layer

	

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

from chain_computer_tests import *

def test_fancy():
	
	import sys
	argv = list( sys.argv )
	## Remove the first item in argv, which is always the program name itself.
	argv.pop(0)
	
	if len( argv ) == 1:
		if argv[0].isdigit():
			paths_info, skeleton_handle_vertices, constraint = eval( 'get_test_infinite(' + argv[0] + ')' )
		else:
			paths_info, skeleton_handle_vertices, constraint = eval( 'get_test_' + argv[0] + '()' )
	else:
		#paths_info, skeleton_handle_vertices, constraint = get_test_steep_closed_curve()
		# paths_info, skeleton_handle_vertices, constraint = get_test1()
		# paths_info, skeleton_handle_vertices, constraint = get_test2()
		paths_info, skeleton_handle_vertices, constraint = get_test_simple_closed()
		#paths_info, skeleton_handle_vertices, constraint = get_test_pebble()
		#paths_info, skeleton_handle_vertices, constraint = get_test_alligator()
		#paths_info, skeleton_handle_vertices, constraint = get_test_box()
	
	engine = YSEngine()
	
	try:
		boundary_index = argmax([ info['bbox_area'] for info in paths_info if info['closed'] ])
	except ValueError:
		boundary_index = -1
	
	engine.init_engine( paths_info, boundary_index )
	
	if constraint is not None:
		engine.constraint_change( constraint[0], constraint[1], constraint[2] )

	engine.set_handle_positions( skeleton_handle_vertices )
	
#	engine.precompute_configuration()
#	engine.prepare_to_solve()
	direct = compute_transformed_by_control_points( engine.all_controls,engine.handle_positions, engine.transforms )
	
	
	## Transform a handle
#	engine.transform_change( 0, [[1,0,-20],[0,1,20]] )
#	
#	all_paths = engine.solve_transform_change()
#	
#	for path in all_paths:
#		if len( path ) > 1:
#			chain = concatenate( asarray(path)[:-1, :-1] )
#			chain = concatenate( ( chain, path[-1] ) )
#		else:
#			chain = path[0]
#		print chain
#		
#	print engine.compute_energy_and_maximum_distance()

def test_actually_solve():
	
	import sys
	argv = list( sys.argv )
	## Remove the first item in argv, which is always the program name itself.
	argv.pop(0)
	
	if len( argv ) == 1:
		if argv[0].isdigit():
			paths_info, skeleton_handle_vertices, constraint = eval( 'get_test_infinite(' + argv[0] + ')' )
		else:
			paths_info, skeleton_handle_vertices, constraint = eval( 'get_test_' + argv[0] + '()' )
	else:
		#paths_info, skeleton_handle_vertices, constraint = get_test_steep_closed_curve()
		#paths_info, skeleton_handle_vertices, constraint = get_test1()
		#paths_info, skeleton_handle_vertices, constraint = get_test2()
		#paths_info, skeleton_handle_vertices, constraint = get_test_simple_closed()
		#paths_info, skeleton_handle_vertices, constraint = get_test_pebble()
		#paths_info, skeleton_handle_vertices, constraint = get_test_alligator()
		#paths_info, skeleton_handle_vertices, constraint = get_test_box()
		paths_info, skeleton_handle_vertices, constraint = get_test_turtle_glasses()
	
	engine = YSEngine()
	
	try:
		boundary_index = argmax([ info['bbox_area'] for info in paths_info if info['closed'] ])
	except ValueError:
		boundary_index = -1
	
	engine.init_engine( paths_info, boundary_index )
	
	if constraint is not None:
		engine.constraint_change( constraint[0], constraint[1], constraint[2] )

	engine.set_handle_positions( skeleton_handle_vertices )
	
	engine.precompute_configuration()
# 	engine.set_weight_function( 'shepard' )
	engine.prepare_to_solve()
	
	
	## Transform a handle
#  	engine.transform_change( 0, [[1,0,-2],[0,1,2]] )
	
	all_paths = engine.solve_transform_change()
	
	print 'diff: ', (asarray( all_paths ) - asarray( engine.all_controls ))[-2:]
	print 'error range: ', max( [ (asarray( path ) - asarray( controls )).max() for path, controls in zip( all_paths, engine.all_controls ) ] )
	
	debugger()
	for path in all_paths:
		if len( path ) > 1:
			chain = concatenate( asarray(path)[:-1, :-1] )
			chain = concatenate( ( chain, path[-1] ) )
		else:
			chain = path[0]
			
# 		print chain
		
# 	print engine.compute_energy_and_maximum_distance()

def test_simple():
	bbw_curve, spline_curve = get_test_distances()					   
	distances = compute_maximum_distances( bbw_curve, spline_curve )
	
	print distances
	print 'HAHA ~ '
		

def main():
	'''
	a console test.
	'''
	
	#test_fancy()
	#test_simple()
	test_actually_solve()


if __name__ == '__main__': main()		
