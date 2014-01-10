from bezier_utility import *
from weights_computer import *

kSystemSolvePackage = 'numpy-solve'
kBuildDense = True
if 'numpy-inv' == kSystemSolvePackage:
	if kBuildDense:
		zeros_system_build_t = zeros
		as_system_build_t = lambda x: x
		to_system_solve_t = lambda x: x
	else:
		zeros_system_build_t = scipy.sparse.lil_matrix
		as_system_build_t = scipy.sparse.lil_matrix
		to_system_solve_t = scipy.sparse.csr_matrix
	
	def compute_symbolic_factorization( system ):
		'''
		Given a scipy.sparse system matrix 'system',
		return a function that can be used to solve for 'x' in
			system * x = b
		as follows:
			x = compute_symbolic_factorization( system )( system )( b )
		'''
		def compute_numeric_factorization( system ):
			try:
				inverse = linalg.inv( system )
			except numpy.linalg.linalg.LinAlgError as e:
				print e
				inverse = eye( len( system ) )
			def solve( rhs ):
				return dot( inverse, rhs )
			return solve
		return compute_numeric_factorization

elif 'numpy-solve' == kSystemSolvePackage:
	if kBuildDense:
		zeros_system_build_t = zeros
		as_system_build_t = lambda x: x
		to_system_solve_t = lambda x: x
	else:
		zeros_system_build_t = scipy.sparse.lil_matrix
		as_system_build_t = scipy.sparse.lil_matrix
		to_system_solve_t = scipy.sparse.csr_matrix
	
	def compute_symbolic_factorization( system ):
		'''
		Given a scipy.sparse system matrix 'system',
		return a function that can be used to solve for 'x' in
			system * x = b
		as follows:
			x = compute_symbolic_factorization( system )( system )( b )
		'''
		def compute_numeric_factorization( system ):
			def solve( rhs ):
				save( 'broken.npy', system )
				return linalg.solve( system, rhs )
			return solve
		return compute_numeric_factorization

elif 'scipy' == kSystemSolvePackage:
	import scipy.sparse.linalg
	## Sparse matrix types.
	zeros_system_build_t = scipy.sparse.lil_matrix
	to_system_solve_t = scipy.sparse.csr_matrix
	## UPDATE: building as a dense matrix is much faster at least for a small system.
	if kBuildDense:
		as_system_build_t = lambda x: x.todense()
	else:
		as_system_build_t = scipy.sparse.lil_matrix
	
	def compute_symbolic_factorization( system ):
		'''
		Given a scipy.sparse system matrix 'system',
		return a function that can be used to solve for 'x' in
			system * x = b
		as follows:
			x = compute_symbolic_factorization( system )( system )( b )
		'''
		def compute_numeric_factorization( system ):
			factorized = scipy.sparse.linalg.factorized( system )
			def solve( rhs ):
				## spsolve() is slower.
				# return scipy.sparse.linalg.spsolve( system, rhs )
				return factorized( rhs )
			return solve
		return compute_numeric_factorization

elif 'cvxopt' == kSystemSolvePackage:
	import cvxopt
	## UPDATE: cholmod dies with our system for some reason.
	#import cvxopt.cholmod as cvxopt_solver
	import cvxopt.umfpack as cvxopt_solver
	## Sparse matrix types.
	zeros_system_build_t = lambda shape: cvxopt.spmatrix( [], [], [], shape )
	
	if kBuildDense:
		## Build dense
		as_system_build_t = lambda x: asarray( cvxopt.matrix( x ) )
		to_system_solve_t = lambda x: cvxopt.sparse( cvxopt.matrix( x ) )
	else:
		## Build sparse
		as_system_build_t = lambda x: x
		to_system_solve_t = lambda x: x
	
	def compute_symbolic_factorization( system ):
		'''
		Given a scipy.sparse system matrix 'system',
		return a function that can be used to solve for 'x' in
			system * x = b
		as follows:
			x = compute_symbolic_factorization( system )( system )( b )
		'''
		symbolic_factorization = cvxopt_solver.symbolic( system )
		def compute_numeric_factorization( system ):
			assert abs(asarray( cvxopt.matrix( (system - system.T) ) )).max() < 1e-8
			full_factorization = cvxopt_solver.numeric( system, symbolic_factorization )
			def solve( rhs ):
				x = cvxopt.matrix( rhs )
				## If cvxopt_solver is cvxopt.umfpack
				cvxopt_solver.solve( system, full_factorization, x )
				## If cvxopt_solver is cvxopt.cholmod
				## UPDATE: Neither of these work.
				# cvxopt_solver.solve( full_factorization, x, sys = 1 )
				# cvxopt_solver.spsolve( system, x, sys = 0 )
				return asarray( x ).reshape( rhs.shape )
			return solve
		return compute_numeric_factorization

else:
	raise RuntimeError( "Unknown system solve package " + str( kSystemSolvePackage ) )


class Bundle( object ):
	def __init__( self, W_matrices, control_points, constraints, mags = None, dirs = None ):
		self.W_matrices = W_matrices
		self.control_points = control_points
		self.constraints = constraints
		controls = asarray(self.control_points)
		if mags is None:
			self.magnitudes = [mag(controls[1] - controls[0]), mag(controls[2] - controls[3])]
		else:
			self.magnitudes = mags
		
		if dirs is None:
			self.directions = [ dir_allow_zero((controls[1] - controls[0])[:2]), dir_allow_zero((controls[2] - controls[3])[:2]) ]
		else:
			self.directions = dirs


class BezierConstraintSolver( object ):
	def __init__( self, W_matrices, control_points, constraints, transforms, is_closed ):
		## compute the weight of each segment according to its length
		num = len(control_points)
		control_points = asarray(control_points)
		self.build_system( W_matrices, control_points, constraints, transforms, is_closed )

	def build_system( self, W_matrices, control_points, constraints, transforms, is_closed ):
		
		### 1 Bundle all data for each bezier curve together
		### 2 Allocate space for the system matrix
		### 3 Gather the pieces of the system for each curve
		### 4 Insert them into the system matrix and right-hand-side
		### 5 Gather the lagrange equations between adjacent curves
		### 6 Insert them into the system matrix and right-hand-side

		### 1
		self.bundles = [ Bundle( W_matrices[i], control_points[i], [constraints[i], constraints[(i+1)%len(control_points)]] ) for i in xrange(len( control_points )) ]
						
		self.dofs_per_bundle = [ self.compute_dofs_per_curve( bundle ) for bundle in self.bundles ]
						
		self.lambdas_per_joint = [ self.constraint_number_per_joint( constraint ) for constraint in constraints]
						
		### 2
		self.total_dofs = sum( self.dofs_per_bundle ) 
		self.system_size = self.total_dofs + sum( self.lambdas_per_joint )
		'''
		test
		'''
		# self.system = zeros( ( self.system_size, self.system_size ) )
		self.system	 = zeros_system_build_t( ( self.system_size, self.system_size ) )
		self.system_symbolic_factorization = None
		self.system_factored = None
		self.rhs = zeros( self.system_size )
		self.transforms = transforms
		self.is_closed = is_closed
		
		self._update_bundles()


	def _update_bundles( self, lagrange_only = False ):
		## For convenience, set local variables from instance variables.
		## WARNING: If you re-assign one of these, the instance variable will not be updated!
		bundles = self.bundles
		dofs_per_bundle = self.dofs_per_bundle
		lambdas_per_joint = self.lambdas_per_joint
		total_dofs = self.total_dofs
		system_size = self.system_size
		## TODO: Figure out the right sparse matrix type.
		system = as_system_build_t( self.system )
		# system = self.system
		rhs = self.rhs
		transforms = self.transforms
		is_closed = self.is_closed

		### 3
		dof_offset = 0
		for i in range(len( bundles )):
			bundle = bundles[i]
			dofs = sum(dofs_per_bundle[i])
			
			if not lagrange_only:
				small_system = self.system_for_curve( bundle )
				small_rhs = self.rhs_for_curve(bundle, transforms)
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
	
			small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_curve_constraints( bundles[i], bundles[i+1] )
			
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
	
			small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_curve_constraints( bundles[-1], bundles[0] )
	
			### 4
			system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, dof_offset : dof_offset + dofs  ] = small_lagrange_system[ :, :dofs ]
			system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, : dofs_next ] = small_lagrange_system[ :, dofs: ]
			rhs[ constraint_equation_offset : constraint_equation_offset + constraint_eqs ] = small_lagrange_rhs

		## Set the upper-right portion of the system matrix, too
		system[ : total_dofs, total_dofs : ] = system.T[ : total_dofs, total_dofs : ]
		
		self.system	 = to_system_solve_t( system )
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

	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def system_for_curve( self, bundle ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def compute_dofs_per_curve( self, bundle ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def constraint_number_per_joint( self, constraint ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )
	def rhs_for_curve( self, bundle, transforms ):
		raise NotImplementedError( "This is an abstract base class. Only call this on a subclass." )




