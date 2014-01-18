from bezier_utility import *
from weights_computer import *

import systems_and_solvers

class Bundle( object ):
	def __init__( self, W_matrices, control_points, constraints, length, ts, dss,  mags = None, dirs = None ):
		self.W_matrices = W_matrices
		self.control_points = control_points
		self.constraints = constraints
		self.length = length
		self.ts = ts
		self.dss = dss
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
			sin_theta = (1.-cos_theta**2) ** 0.5
		else:
			cos_theta = 1.0
			sin_theta = 0.0
			
		if cross(vec0, vec1) < 0:	
			sin_theta = -sin_theta
			
		return [ cos_theta, sin_theta ]

class BezierConstraintSolver( object ):
	def __init__( self, W_matrices, control_points, constraints, transforms, lengths, ts, dss, is_closed, kArcLength ):
		## compute the weight of each segment according to its length
		# num = len(control_points)
		control_points = asarray(control_points)
		self.build_system( W_matrices, control_points, constraints, transforms, lengths, ts, dss, is_closed, kArcLength )

	def build_system( self, W_matrices, control_points, constraints, transforms, lengths, ts, dss, is_closed, kArcLength ):
		
		### 1 Bundle all data for each bezier curve together
		### 2 Allocate space for the system matrix
		### 3 Gather the pieces of the system for each curve
		### 4 Insert them into the system matrix and right-hand-side
		### 5 Gather the lagrange equations between adjacent curves
		### 6 Insert them into the system matrix and right-hand-side

		### 1
		self.bundles = [ Bundle( W_matrices[i], control_points[i], [constraints[i], constraints[(i+1)%len(control_points)]], lengths[i], ts[i], dss[i] ) for i in xrange(len( control_points )) ]
		
		if is_closed:
			self.angles = [ compute_angle( self.bundles[i], self.bundles[(i+1)%len( self.bundles ) ] ) for i in xrange( len( self.bundles ) ) ]
		else:
			self.angles = [ compute_angle( self.bundles[i], self.bundles[i+1] ) for i in xrange( len( self.bundles )-1 ) ]
						
		self.dofs_per_bundle = [ self.compute_dofs_per_curve( bundle ) for bundle in self.bundles ]
						
		self.lambdas_per_joint = [ self.constraint_number_per_joint( constraint ) for constraint in constraints]
		
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




