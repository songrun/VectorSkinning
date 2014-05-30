from generate_chain_system import *

dim = 3

class BezierConstraintSolverOddFast( BezierConstraintSolver ):
	'''
	Free direction, magnitude fixed (for G1 or A).
	'''
	
	def _update_bundles( self, lagrange_only = False ):
		if not lagrange_only:
			## The default from the superclass doesn't work for us.
			self.rhs = zeros( ( len( self.transforms ), self.system_size, 10 ) )
		
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
				rhs[ :, dof_offset : dof_offset + dofs, :9 ] = small_rhs
	
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
			rhs[ :, constraint_equation_offset : constraint_equation_offset + constraint_eqs, 9 ] = small_lagrange_rhs
	
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
			rhs[ :, constraint_equation_offset : constraint_equation_offset + constraint_eqs, 9 ] = small_lagrange_rhs
			constraint_equation_offset += constraint_eqs
			
		else:
			dofs_head = sum(dofs_per_bundle[0])
			dofs_tail = sum(dofs_per_bundle[-1])
			constraint_eqs = 2
			
			if lambdas_per_joint[-1] == 2:
				small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_fixed_opening( bundles[-1], is_head = False )
				system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, dof_offset : dof_offset + dofs_tail  ] = small_lagrange_system
				rhs[ :, constraint_equation_offset : constraint_equation_offset + constraint_eqs, 9 ] = small_lagrange_rhs
				constraint_equation_offset += constraint_eqs
				
			if lambdas_per_joint[0] == 2:
				small_lagrange_system, small_lagrange_rhs = self.lagrange_equations_for_fixed_opening( bundles[0], is_head = True )							
				system[ constraint_equation_offset : constraint_equation_offset + constraint_eqs, dof_offset : dof_offset + dofs_tail  ] = small_lagrange_system
				rhs[ :, constraint_equation_offset : constraint_equation_offset + constraint_eqs, 9 ] = small_lagrange_rhs
				constraint_equation_offset += constraint_eqs
				
				
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
			self.rhs[ :, dof_offset : dof_offset + dofs, :9 ] = small_rhs
	
			dof_offset += dofs

		assert dof_offset == self.total_dofs
	
	def update_system_with_result_of_previous_iteration( self, solution ):
		### Iterate only over the parts of the matrix that will change,
		### such as the lagrange multipliers across G1 or A edges and the right-hand-side.
		solution = asarray(solution)
		num = len(self.bundles)
		assert solution.shape == (num, 4, 2)
		
		for i in range(num):
			dir1 = dir_allow_zero( solution[i][1]-solution[i][0] )
			dir2 = dir_allow_zero( solution[i][2]-solution[i][3] )
			
			self.bundles[i].directions[0] = dir1
			self.bundles[i].directions[1] = dir2
			
			mag1 = mag( solution[i][1]-solution[i][0] )
			mag2 = mag( solution[i][2]-solution[i][3] )
			
			self.bundles[i].magnitudes[0] = mag1
			self.bundles[i].magnitudes[1] = mag2			
		
		
		## The lagrange multipliers changed, but not the locations of the zeros.
		self._update_bundles( lagrange_only = True )
		self.system_factored = None
		## UPDATE: Actually, if fixed angles are parallel or perpendicular,
		##		   then the lagrange multiplier systems may gain
		##		   or lose zeros. So, reset the symbolic factorization.
		## UPDATE 2: If we could update_bundles once with all fixed angles
		##			 not parallel or perpendicular, and then compute the symbolic
		##			 factorization, we could keep it.
		## UPDATE 3: Let's try it assuming that the first time through there are no zeros.
		## UPDATE 4: I tried it and it makes no difference to performance at all
		##			 up to alec's alligator. So, we'll reset the symbolic factorization
		##			 in case the initial configuration has zeros.
		self.system_symbolic_factored = None
		self.Os = None
		
	
	def solve( self ):
		num = len(self.bundles)
		
		#print 'rhs:'
		#print self.rhs.tolist()
		
		if self.system_symbolic_factored is None:
			#print 'odd symbolic factoring'
			system = self.to_system_solve_t( self.system )
			self.system_symbolic_factored = self.compute_symbolic_factorization( system )
			self.system_factored = self.system_symbolic_factored( system )
			self.Os = None
		
		elif self.system_factored is None:
			#print 'odd numeric factoring'
			system = self.to_system_solve_t( self.system )
			self.system_factored = self.system_symbolic_factored( system )
			self.Os = None
		
		if self.Os is None:
			self.Os = []
			for i in xrange(len( self.Ts )):
				## Doesn't take a matrix:
				# self.Os.append( self.system_factored( self.rhs[i] ) )
				
				solved = zeros( self.rhs[i].shape )
				for j in xrange( self.rhs[i].shape[1] ):
					solved[:,j] = self.system_factored( self.rhs[i,:,j] )
				self.Os.append( solved )
		
		x = zeros( self.rhs.shape[1] )
		for i in xrange(len( self.Ts )):
			T = self.Ts[i]
			O = self.Os[i]
			x += dot( O, append( T.ravel(), 1. ) )
		#x = x[:2,:]
		
		#print 'odd solve'
		# x = self.system_factored( self.rhs )
		# x = linalg.solve( self.system, self.rhs )
		# x = scipy.sparse.linalg.spsolve( self.system, self.rhs )
		### Return a nicely formatted chain of bezier curves.
		x = array( x[:self.total_dofs] ).reshape(-1,4).T
		
		solution = []
		for i in range(num):
			P = x[:, i*dim:(i+1)*dim ]
			solution.append( P[:,:2] )
		
		if parameters.kClampOn == True: solution = clamp_solution( self.bundles, solution )
		
		return solution 
			
	def lagrange_equations_for_fixed_opening( self, bundle, is_head ):
		## handle the case of open end path.
		dofs = self.compute_dofs_per_curve(bundle)
		
		R = zeros( ( sum(dofs), dim ) )
		rhs = zeros(R.shape[1])
		
		if is_head: 
#			assert bundle.constraints[0][1] == True
			fixed_positions = bundle.control_points[0][:2]
			fixed_positions = asarray(fixed_positions)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P1x' - constraint_X' ) = 0
			lambda2 * ( P1y' - constraint_Y' ) = 0
			'''
			for i in range( dim ):
				R[i*4, i] = 1
		
			rhs = fixed_positions
		else:
#			assert bundle.constraints[-1][1] == True
			fixed_positions = bundle.control_points[-1][:2]
			fixed_positions = asarray(fixed_positions)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			for i in range( dim ):
				R[i*4+3, i] = 1
		
			rhs = fixed_positions
				
		return R.T, rhs
	
	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1, angle ):
		mag0, mag1 = bundle0.magnitudes[1], bundle1.magnitudes[0]
		cos_theta = angle[0]
		sin_theta = angle[1]
		
		dofs0 = self.compute_dofs_per_curve(bundle0)
		dofs1 = self.compute_dofs_per_curve(bundle1)
		dofs = sum(dofs0) + sum(dofs1)

		smoothness = bundle0.constraints[1][0]
		if smoothness == 'C0':			## C0
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			'''
			R = zeros( ( dofs, dim ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0) + i*4, i] = -1

		elif smoothness == 'A':		 ## fixed angle
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x - Q1x ) = 0
			lambda2 * ( P4y - Q1y ) = 0
			lambda3 * ( mag1(P4x-P3x) + mag0[cos_theta(Q2x-Q1x)-sin_theta(Q2y-Q1y)] ) = 0
			lambda4 * ( mag1(P4y-P3y) + mag0[sin_theta(Q2x-Q1x)+cos_theta(Q2y-Q1y)] ) = 0
			'''
			R = zeros( ( dofs, dim+2 ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0)+i*4, i] = -1
			
			## Rotations only happen like this in the first two dimensions.
			for i in range( 2 ):
				R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
			
			##          Qn x   eq
			R[sum(dofs0)+1+0*4, 0+dim] = cos_theta
			R[sum(dofs0)+0+0*4, 0+dim] = -cos_theta
			R[sum(dofs0)+1+1*4, 0+dim] = -sin_theta
			R[sum(dofs0)+0+1*4, 0+dim] = sin_theta
			
			R[sum(dofs0)+1+0*4, 1+dim] = sin_theta
			R[sum(dofs0)+0+0*4, 1+dim] = -sin_theta
			R[sum(dofs0)+1+1*4, 1+dim] = cos_theta
			R[sum(dofs0)+0+1*4, 1+dim] = -cos_theta
			
			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= mag1
			R[ sum(dofs0):, dim: ] *= -mag0
			
		elif smoothness == 'C1':		 ## C1
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			lambda3 * ( w_q(P4x' - P3x') + w_p(Q1x' - Q2x')) = 0
			lambda4 * ( w_q(P4y' - P3y') + w_p(Q1y' - Q2y')) = 0
			'''
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0)+i*4, i] = -1
				
				R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				R[sum(dofs0)+i*4+1, i+dim] = -1
				R[sum(dofs0)+i*4, i+dim] = 1

			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= mag1
			R[ sum(dofs0):, dim: ] *= mag0
		elif smoothness == 'G1':		 ## G1
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0)+i*4, i] = -1
				
				R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				R[sum(dofs0)+i*4+1, i+dim] = -1
				R[sum(dofs0)+i*4, i+dim] = 1

			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= mag1
			R[ sum(dofs0):, dim: ] *= mag0
		
		else:
			R = zeros( ( dofs, 0 ) )
		
		rhs = zeros(R.shape[1])
		
		fixed_positions = bundle0.control_points[-1][:2]
		is_fixed = bundle0.constraints[1][1]
		assert type( is_fixed ) == bool
		if is_fixed:

			fixed_positions = asarray(fixed_positions)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			R2 = zeros( ( dofs, dim ) )
			for i in range( dim ):
				R2[i*4+3, i] = 1
		
			R = concatenate((R, R2), axis=1)
			rhs = concatenate((rhs, fixed_positions))
	
		return R.T, rhs
		
		
	def system_for_curve( self, bundle ):
		'''
		## A is computed using Sage, integral of (tbar.T * tbar) with respect to t.
		#	A = asarray( [[	 1./7,	1./6,  1./5, 1./4], [ 1./6,	 1./5, 1./4,  1./3], 
		#		[ 1./5, 1./4,  1./3,  1./2], [1./4,	 1./3,	1./2,	1.]] )
		## MAM is computed using Sage. MAM = M * A * M
		'''
		length = bundle.length
		MAM = asarray( self.MAM )
		
		Left = zeros((4*dim, 4*dim))

		for i in range(dim):		
			Left[ i*4:(i+1)*4, i*4:(i+1)*4 ] = MAM[:,:]
		
		return Left*length
		
	def system_for_curve_with_arc_length( self, bundle ):
		'''
		## Solve the same integral as system__for_curve only with dt replaced by ds
		'''
		length = bundle.length
		ts = bundle.ts
		dts = bundle.dts
		Left = zeros( ( 4*dim, 4*dim ) )
		tbar = ones( ( 4, 1 ) )
		MAM = zeros( ( 4, 4 ) )
		
		for i in range(len(dts)):
			t = (ts[i] + ts[i+1])/2
			ds = dts[i]
			
			tbar[0] = t*t*t
			tbar[1] = t*t
			tbar[2] = t
			
			Mtbar = dot( M.T, tbar )

			MAM += dot( Mtbar, Mtbar.T )*ds
	
		for i in range( dim ):		
			Left[ i*4:( i+1 )*4, i*4:( i+1 )*4 ] = MAM[:,:]
		
		return Left*length
			
			
	def compute_dofs_per_curve( self, bundle ):
	
		dofs = zeros( 2, dtype = int )
		'''
		assume open end points can only emerge at the endpoints
		'''
		for i, (smoothness, is_fixed) in enumerate(bundle.constraints):
			if smoothness == 'C0': dofs[i] += 4			## C0
			elif smoothness == 'A': dofs[i] += 4		## fixed angle
			elif smoothness == 'C1': dofs[i] += 4		## C1
			elif smoothness == 'G1': dofs[i] += 4		## G1
			elif smoothness == 'None': dofs[i] += 4		## Free of constraint
			
		return (dofs/2)*dim
		
	
	def constraint_number_per_joint(self, constraint ):	   
	
		assert len(constraint) == 2
		smoothness = constraint[0]	
		is_fixed = constraint[1]	 
		
		num = 0
		if smoothness == 'C0': num = dim			## C0
		elif smoothness == 'A': num = dim+2		## fixed angle
		elif smoothness == 'C1': num = 2*dim		## C1
		elif smoothness == 'G1': num = 2*dim		## G1
		
		assert type( is_fixed ) == bool
		if is_fixed:
			num += dim
			
		return num
		
		
	def rhs_for_curve( self, bundle, transforms ):
		'''
		The rhs is computed according to the formula:
			rhs = sum(Ti * P.T * M.T * W_i * M)
		'''
		length = bundle.length
		
#		W_matrices = bundle.W_matrices
#		controls = bundle.control_points
#  
#		Right = zeros( (3, 4) )
#		for i in range( len( transforms ) ):
# 
#			T_i = mat( asarray(transforms[i]).reshape(3,3) )
#			W_i = W_matrices[i,0]	
# 
#			Right = Right + T_i * (controls.T) * M * mat( W_i ) * M
# 
#		Right = asarray(Right).reshape(-1)
#		Right = Right[:8]	

		W_matrices = bundle.W_matrices
		controls = bundle.control_points
		
		Right = zeros( ( len( transforms ), 4*dim, 9 ) )
		self.Os = None
		self.Ts = []
		for i in range( len( transforms ) ):
		
			T_i = mat( asarray(transforms[i]).reshape(3, 3) )
			W_i = asarray(W_matrices[i])

			## 3x3 * 3xN * NxN * NxN
			# temp = temp + dot(asarray(T_i*(controls.T)*M), W_i)
			
			#vec(C)A = vec(TCMW)
			## matrix calculus says:
			#vec(C)A = (I kronecker T) vec(CMW)
			
			self.Ts.append( asarray(transforms[i]).reshape(3, 3) )
			Right[ i ] = kron( identity(3), dot( dot( controls.T, M ), W_i ).T )

		return Right*length
	
	def update_rhs_for_handles( self, transforms ):
		self.Ts = []
		for i in range( len( transforms ) ):
			self.Ts.append( asarray(transforms[i]).reshape(3, 3) )
