from bezier_utility import *
from weights_computer import *

import systems_and_solvers

class Bundle( object ):
	def __init__( self, W_matrices, control_points, constraints, length, ts, dts,  mags = None, dirs = None ):
		self.W_matrices = W_matrices
		self.control_points = control_points
		self.constraints = constraints
		self.length = length
		self.ts = ts
		self.dts = dts
		controls = asarray(self.control_points)
		if mags is None:
			self.magnitudes = [mag(controls[1] - controls[0]), mag(controls[2] - controls[3])]
		else:
			self.magnitudes = mags
		
		if dirs is None:
					
 			self.directions = [ dir_allow_zero((controls[1] - controls[0])[:2]), dir_allow_zero((controls[2] - controls[3])[:2]) ]
		else:
			self.directions = dirs
				

def compute_angle( bundle0, bundle1 ):
		## search out for the first control point making non-zero directions
		vec0, vec1 = zeros( 2 ), zeros( 2 )
		for j in range( 2, 5 ):
			if not array_equal( bundle0.control_points[-j], bundle0.control_points[3] ):
				vec0 = (bundle0.control_points[-j]-bundle0.control_points[3])[:2]
				break
				
		for j in range( 1, 4 ):
			if not array_equal( bundle1.control_points[j], bundle1.control_points[0] ):
				vec1 = (bundle1.control_points[j]-bundle1.control_points[0])[:2]
				break
		
		if  mag(vec0)*mag(vec1) != 0:
			cos_theta = dot(vec0, vec1)/( mag(vec0)*mag(vec1) )
			sin_theta = cross(vec0, vec1)/( mag(vec0)*mag(vec1) )
		else:
			cos_theta = 1.0
			sin_theta = 0.0
			
		
		return [ cos_theta, sin_theta ]

class BezierConstraintSolver( object ):
	def __init__( self, W_matrices, control_points, constraints, transforms, lengths, ts, dts, is_closed, kArcLength ):
		## compute the weight of each segment according to its length
		# num = len(control_points)
		control_points = asarray(control_points)
		self.build_system( W_matrices, control_points, constraints, transforms, lengths, ts, dts, is_closed, kArcLength )

	def build_system( self, W_matrices, control_points, constraints, transforms, lengths, ts, dts, is_closed, kArcLength ):
		
		### 1 Bundle all data for each bezier curve together
		### 2 Allocate space for the system matrix
		### 3 Gather the pieces of the system for each curve
		### 4 Insert them into the system matrix and right-hand-side
		### 5 Gather the lagrange equations between adjacent curves
		### 6 Insert them into the system matrix and right-hand-side

		### 1
		self.bundles = [ Bundle( W_matrices[i], control_points[i], [constraints[i], constraints[(i+1)%len(control_points)]], lengths[i], ts[i], dts[i] ) for i in xrange(len( control_points )) ]
		
		if is_closed:
			self.angles = [ compute_angle( self.bundles[i], self.bundles[(i+1)%len( self.bundles ) ] ) for i in xrange( len( self.bundles ) ) ]
		else:
			self.angles = [ compute_angle( self.bundles[i], self.bundles[i+1] ) for i in xrange( len( self.bundles )-1 ) ]
						
		self.dofs_per_bundle = [ self.compute_dofs_per_curve( bundle ) for bundle in self.bundles ]
						
		self.lambdas_per_joint = [ self.constraint_number_per_joint( constraint ) for constraint in constraints]

		self.MAM = self.get_default_MAM()
		self.coefficient_matrices = self.get_default_left_matrices_for_even_iterations()
		
		### 2
		self.total_dofs = sum( self.dofs_per_bundle ) 
		self.system_size = self.total_dofs + sum( self.lambdas_per_joint )
		# print 'system_size:', self.system_size
		'''
		test
		'''
		
		# self.zeros_system_build_t, self.to_system_solve_t, self.compute_symbolic_factorization = systems_and_solvers.get_system_and_factor_funcs()
		smoothness = [ constraint[0] for constraint in constraints ]
		G1orA = 'A' in smoothness or 'G1' in smoothness
		del smoothness
		self.zeros_system_build_t, self.to_system_solve_t, self.compute_symbolic_factorization = systems_and_solvers.get_system_and_factor_funcs_for_system_size( self.system_size, G1orA )
		
		# self.system = zeros( ( self.system_size, self.system_size ) )
		self.system	 = self.zeros_system_build_t( ( self.system_size, self.system_size ) )
		self.system_symbolic_factored = None
		self.system_factored = None
		self.rhs = zeros( self.system_size )
		self.transforms = transforms
		self.is_closed = is_closed
		self.kArcLength = kArcLength
		
		self._update_bundles( )


	def _update_bundles( self, lagrange_only = False ):
		## For convenience, set local variables from instance variables.
		## WARNING: If you re-assign one of these, the instance variable will not be updated!
		bundles = self.bundles
		dofs_per_bundle = self.dofs_per_bundle
		lambdas_per_joint = self.lambdas_per_joint
		total_dofs = self.total_dofs
		system_size = self.system_size
		system = self.system
		angles = self.angles
		kArcLength = self.kArcLength
		rhs = self.rhs
		transforms = self.transforms
		is_closed = self.is_closed

		### 3
		dof_offset = 0
		for i in range(len( bundles )):
			bundle = bundles[i]
			dofs = sum(dofs_per_bundle[i])
			
			if not lagrange_only:
# 				small_system = self.system_for_curve_with_arc_length( bundle )
				
				if kArcLength:
					small_system = self.system_for_curve_with_arc_length( bundle )
				else:
					small_system = self.system_for_curve( bundle )
				
				small_rhs = self.rhs_for_curve( bundle, transforms)
				### 4
				system[ dof_offset : dof_offset + dofs, dof_offset : dof_offset + dofs ] = small_system
				rhs[ dof_offset : dof_offset + dofs ] = small_rhs
	
			dof_offset += dofs

		assert dof_offset == total_dofs
		
				
		### 5
		dof_offset = 0
		constraint_equation_offset = total_dofs
		for i in range( len( bundles ) - 1 ):
			dofs = sum(dofs_per_bundle[i])
			dofs_next = sum(dofs_per_bundle[i+1])
			constraint_eqs = lambdas_per_joint[i+1]
	
			small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_curve_constraints( bundles[i], bundles[i+1], angles[i] )
			
			### 4
			system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, dof_offset : dof_offset + dofs + dofs_next ] = small_lagrange_system
			rhs[ constraint_equation_offset : constraint_equation_offset + constraint_eqs ] = small_lagrange_rhs
	
			dof_offset += dofs
			constraint_equation_offset += constraint_eqs

		## Handle the connection between the last and first bezier curves if it is a closed curve.
		if is_closed:
			dofs = sum(dofs_per_bundle[-1])
			dofs_next = sum(dofs_per_bundle[0])
			constraint_eqs = lambdas_per_joint[0]
	
			small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_curve_constraints( bundles[-1], bundles[0], angles[-1] )
	
			### 4
			system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, dof_offset : dof_offset + dofs  ] = small_lagrange_system[ :, :dofs ]
			system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, : dofs_next ] = small_lagrange_system[ :, dofs: ]
			rhs[ constraint_equation_offset : constraint_equation_offset + constraint_eqs ] = small_lagrange_rhs
		
# 		if len ( self.bundles ) == 1:
# 			small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_single_bundle( bundles[0] )	
				

		## Set the upper-right portion of the system matrix, too
		system[ : total_dofs, total_dofs : ] = system.T[ : total_dofs, total_dofs : ]
		
		# self.system	 = self.to_system_solve_t( system )
		## Reset system_factored, but leave 'self.system_symbolic_factorization' alone,
		## because we didn't change the sparsity pattern of the matrix.
		## UPDATE: Actually, if constrained directions align with coordinate axes
		##         or have zero magnitude, then the matrix may gain
		##		   or lose zeros.
		##         So, whoever calls this function should make sure to set whichever
		##         ones should be set to None.
		# self.system_factored = None
		
	
	def update_rhs_for_handles( self, transforms ):
		dof_offset = 0
		for i in range(len( self.bundles )):
			bundle = self.bundles[i]
			dofs = sum(self.dofs_per_bundle[i])
	
			small_rhs = self.rhs_for_curve( bundle, transforms )
			### 4
			self.rhs[ dof_offset : dof_offset + dofs ] = small_rhs
	
			dof_offset += dofs

		assert dof_offset == self.total_dofs

	def lagrange_equations_for_single_bundle( self, bundle ):
		
		assert len( self.bundles ) == 1
		dof = sum( self.compute_dofs_per_curve( bundle ) )
		dim = 2
		R = zeros( ( dof, dim*2 ) )
		rhs = zeros( dim * 2 )
		assert type( bundle.constraints[0][1] ) == bool
		if bundle.constraints[0][1]:

			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			
			for i in range( dim ):
				R[i*4, i] = 1
	
			rhs[ :dim ] = bundle.control_points[0][ :dim ]
				
		assert type( bundle.constraints[1][1] ) == bool
		if bundle.constraints[1][1]:

			fixed = asarray(fixed)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			for i in range( dim ):
				R[i*4+3, dim+i] = 1
				
			rhs[ -dim: ] = bundle.control_points[-1][:dim]

		return R.T, rhs

	def get_default_MAM( self, num_samples = 100 ):
		
		ts = linspace( 0, 1, num_samples )
		dts = ones( num_samples-1 ) * (1./(num_samples-1) )
		
		tbar = ones( ( 4, 1 ) )
		MAM = zeros( ( 4, 4 ) )
# 		MAM = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,	3./35, 9./140,	1./35], [ 1./35, 9./140,  3./35,  1./14], [1./140,	1./35,	1./14,	 1./7]] )
		for i in range(len(dts)):
			t = (ts[i] + ts[i+1])/2
			dt = dts[i]
			
			tbar[0] = t**3
			tbar[1] = t**2
			tbar[2] = t
			
			Mtbar = dot( M.T, tbar )

			MAM += dot( Mtbar, Mtbar.T )*dt

		return MAM
		
	def get_default_left_matrices_for_even_iterations( self ):
	
		MAM = self.MAM	
		if MAM == None:	
			MAM = self.get_default_MAM()
			
		Left1 = zeros( ( 8, 8 ) )	
		Left1[ : : 2, : : 2 ] = MAM
		Left1[ 1: : 2, 1: : 2 ] = MAM
		
		Left2 = zeros( ( 7,  7 ) )
		## right bottom
		Left2[ 3: : 2, 3: : 2 ] = MAM[-2:,-2:]
		Left2[ 4: : 2, 4: : 2 ] = MAM[-2:,-2:]
		## left botton
		Left2[ 3: : 2, 0 ] = Left2[ 4: : 2, 1 ] = MAM[-2:,0] + MAM[-2:,1]
		Left2[ 3: : 2, 2 ] = Left2[ 4: : 2, 2 ] = MAM[-2:, 1]
		## right top
		Left2[ :3, 3: ] = Left2[ 3:, :3 ].T
		## left top
		Left2[0,0] = Left2[1,1] = MAM[0,0] + MAM[0,1] + MAM[1,0] + MAM[1,1]
		Left2[2,0] = Left2[0,2] = Left2[2,1] = Left2[1,2] = ( MAM[1,0] + MAM[1,1] )
		Left2[2,2] = MAM[1,1]
		
		Left3 = zeros( ( 7,  7 ) )
		## left top
		Left3[ :4 : 2, :4 : 2 ] = Left3[ 1:4 : 2, 1:4 : 2 ] = MAM[:2, :2]
		## right top
		Left3[ :4 : 2, -3 ] = Left3[ 1:4 : 2, -2 ] = MAM[:2,-2] + MAM[:2,-1]
		Left3[ :4 : 2, -1 ] = Left3[ 1:4 : 2, -1 ] = MAM[:2, -2]
		## left bottom
		Left3[ 4:, :4 ] = Left3[ :4, 4: ].T
		## right bottom
		Left3[-3,-3] = Left3[-2,-2] = MAM[-1,-1] + MAM[-1,-2] + MAM[-2,-1] + MAM[-2,-2]
		Left3[-1,-3] = Left3[-3,-1] = Left3[-1,-2] = Left3[-2,-1] = ( MAM[-2,-2] + MAM[-2,-1] )
		Left3[-1,-1] = MAM[-2,-2]
		
		Left4 = zeros( ( 6,  6 ) )
		## left top
		Left4[0,0] = Left4[1,1] = MAM[0,0] + MAM[0,1] + MAM[1,0] + MAM[1,1]
		Left4[2,0] = Left4[0,2] = Left4[2,1] = Left4[1,2] = ( MAM[1,0] + MAM[1,1] )
		Left4[2,2] = MAM[1,1]
		## right top
		Left4[0,-3] = Left4[1,-2] = MAM[0,-2] + MAM[0,-1] + MAM[1,-2] + MAM[1,-1]
		Left4[2,-3] = Left4[2,-2] = ( MAM[0,-2] + MAM[1,-2] ) 
		Left4[0,-1] = Left4[1,-1] = ( MAM[1,-2] + MAM[1,-1] )
		Left4[2,-1] = MAM[1,-2]
		## left bottom
		Left4[ 3:, :3 ] = Left4[ :3, 3: ].T
		## right bottom
		Left4[-3,-3] = Left4[-2,-2] = MAM[-1,-1] + MAM[-1,-2] + MAM[-2,-1] + MAM[-2,-2]
		Left4[-1,-3] = Left4[-3,-1] = Left4[-1,-2] = Left4[-2,-1] = ( MAM[-2,-2] + MAM[-2,-1] )
		Left4[-1,-1] = MAM[-2,-2]
		'''
		standard1 = asarray( [[1./7, 0., 1./14, 0., 1./35, 0., 1./140, 0.],
					[0., 1./7, 0., 1./14, 0., 1./35, 0., 1./140],
					[1./14, 0., 3./35, 0., 9./140, 0., 1./35, 0.],
					[0., 1./14, 0., 3./35, 0., 9./140, 0., 1./35],
					[1./35, 0., 9./140, 0., 3./35, 0., 1./14, 0.],
					[0., 1./35, 0., 9./140, 0., 3./35, 0., 1./14],
					[1./140, 0., 1./35, 0., 1./14, 0., 1./7, 0.],
					[0., 1./140, 0., 1./35, 0., 1./14, 0., 1./7]] )
					
		standard2 = asarray( [[13./35, 0., 11./70, 13./140, 0., 1./28, 0.],
					[0., 13./35, 11./70, 0., 13./140, 0., 1./28],
					[11./70, 11./70, 3./35, 9./140, 9./140, 1./35, 1./35],
					[13./140, 0., 9./140, 3./35, 0., 1./14, 0.],
					[0., 13./140, 9./140, 0., 3./35, 0., 1./14],
					[1./28, 0., 1./35, 1./14, 0., 1./7, 0.],
					[0., 1./28, 1./35, 0., 1./14, 0., 1./7]] )			
					
		standard3 = asarray(
					[[1./7, 0., 1./14, 0., 1./28, 0., 1./35],
					[0., 1./7, 0., 1./14, 0., 1./28, 1./35],
					[1./14, 0., 3./35, 0., 13./140, 0., 9./140],
					[0., 1./14, 0., 3./35, 0., 13./140, 9./140],
					[1./28, 0., 13./140, 0., 13./35, 0., 11./70],
					[0., 1./28, 0., 13./140, 0., 13./35, 11./70],
					[1./35, 1./35, 9./140, 9./140, 11./70, 11./70, 3./35]]
					)	
		standard4 = asarray(
					[[13./35, 0., 11./70, 9./70, 0., 13./140],
					[0., 13./35, 11./70, 0., 9./70, 13./140],
					[11./70, 11./70, 3./35, 13./140, 13./140, 9./140],
					[9./70, 0., 13./140, 13./35, 0., 11./70],
					[0., 9./70, 13./140, 0., 13./35, 11./70],
					[13./140, 13./140, 9./140, 22./140, 22./140, 3./35]]
					)
					
		'''						
		'''
		## p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y
		if array_equal(dofs, (4,4)):
			Left = [[1./7, 0., 1./14, 0., 1./35, 0., 1./140, 0.],
					[0., 1./7, 0., 1./14, 0., 1./35, 0., 1./140],
					[1./14, 0., 3./35, 0., 9./140, 0., 1./35, 0.],
					[0., 1./14, 0., 3./35, 0., 9./140, 0., 1./35],
					[1./35, 0., 9./140, 0., 3./35, 0., 1./14, 0.],
					[0., 1./35, 0., 9./140, 0., 3./35, 0., 1./14],
					[1./140, 0., 1./35, 0., 1./14, 0., 1./7, 0.],
					[0., 1./140, 0., 1./35, 0., 1./14, 0., 1./7]]
		## p1x, p1y, s, p3x, p3y, p4x, p4y
		elif array_equal(dofs, (3,4)):
			Left = [[13./35, 0., 11./70*dirs[0,0], 13./140, 0., 1./28, 0.],
					[0., 13./35, 11./70*dirs[0,1], 0., 13./140, 0., 1./28],
					[11./70*dirs[0,0], 11./70*dirs[0,1], 3./35*mag2(dirs[0]), 9./140*dirs[0,0], 9./140*dirs[0,1], 1./35*dirs[0,0], 1./35*dirs[0,1]],
					[13./140, 0., 9./140*dirs[0,0], 3./35, 0., 1./14, 0.],
					[0., 13./140, 9./140*dirs[0,1], 0., 3./35, 0., 1./14],
					[1./28, 0., 1./35*dirs[0,0], 1./14, 0., 1./7, 0.],
					[0., 1./28, 1./35*dirs[0,1], 0., 1./14, 0., 1./7]]	
		## p1x, p1y, p2x, p2y, p4x, p4y, u
		elif array_equal(dofs, (4,3)):
			Left = [[1./7, 0., 1./14, 0., 1./28, 0., 1./35*dirs[1,0]],
					[0., 1./7, 0., 1./14, 0., 1./28, 1./35*dirs[1,1]],
					[1./14, 0., 3./35, 0., 13./140, 0., 9./140*dirs[1,0]],
					[0., 1./14, 0., 3./35, 0., 13./140, 9./140*dirs[1,1]],
					[1./28, 0., 13./140, 0., 13./35, 0., 11./70*dirs[1,0]],
					[0., 1./28, 0., 13./140, 0., 13./35, 11./70*dirs[1,1]],
					[1./35*dirs[1,0], 1./35*dirs[1,1], 9./140*dirs[1,0], 9./140*dirs[1,1], 11./70*dirs[1,0], 11./70*dirs[1,1], 3./35*mag2(dirs[1])]]	
		## p1x, p1y, s, p4x, p4y, u
		elif array_equal(dofs, (3,3)):
			Left = [[13./35, 0., 11./70*dirs[0,0], 9./70, 0., 13./140*dirs[1,0]],
					[0., 13./35, 11./70*dirs[0,1], 0., 9./70, 13./140*dirs[1,1]],
					[11./70*dirs[0,0], 11./70*dirs[0,1], 3./35*mag2(dirs[0]), 13./140*dirs[0,0], 13./140*dirs[0,1], 9./140*dot(dirs[0], dirs[1])],
					[9./70, 0., 13./140*dirs[0,0], 13./35, 0., 11./70*dirs[1,0]],
					[0., 9./70, 13./140*dirs[0,1], 0., 13./35, 11./70*dirs[1,1]],
					[13./140*dirs[1,0], 13./140*dirs[1,1], 9./140*dot(dirs[0], dirs[1]), 22./140*dirs[1,0], 22./140*dirs[1,1], 3./35*mag2(dirs[1])]]	
		'''	
		
		return [ Left1, Left2, Left3, Left4 ]
		

	### For subclasses to implement:
	def update_system_with_result_of_previous_iteration( self, previous_solution ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )

	def solve( self ):
		### Return a nicely formatted chain of bezier curves,
		### even if some of the variables were substituted in the actual system matrix.
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )

	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1, angle):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def system_for_curve( self, bundle ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def compute_dofs_per_curve( self, bundle ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def constraint_number_per_joint( self, constraint ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def rhs_for_curve( self, bundle, transforms ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	### solve by arc length parameterization
	def system_for_curve_with_arc_length( self, bundle ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )




