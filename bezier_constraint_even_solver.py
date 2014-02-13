from generate_chain_system import *

class BezierConstraintSolverEven( BezierConstraintSolver ):
	'''
	Fixed direction, magnitude free (for G1 or A).
	'''
	
	def update_system_with_result_of_previous_iteration( self, solution ):
		### Iterate only over the parts of the matrix that will change,
		### such as the system matrix involving fixed directions,
		### lagrange multipliers across G1 or A edges, and the right-hand-side.	
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
			
		
		self._update_bundles( )
		self.system_factored = None
		## UPDATE: Actually, if constrained directions align with coordinate axes
		##         or have zero magnitude, then the systems may gain
		##		   or lose zeros. So, reset the symbolic factorization.
		## UPDATE 2: If we could update_bundles once with all directions zero-free,
		##           and then compute the symbolic factorization, we could keep it.
		## UPDATE 3: Let's try it assuming that the first time through there are no zeros.
		## UPDATE 4: I tried it and it makes no difference to performance at all
		##           up to alec's alligator. So, we'll reset the symbolic factorization
		##           in case the initial configuration has zeros.
		self.system_symbolic_factored = None
		
	
	def solve( self ):
		
		### Return a nicely formatted chain of bezier curves, un-substituting the fixed
		### directions.
		dofs_per_bundle = self.dofs_per_bundle
		dirs_per_bundle = [bundle.directions for bundle in self.bundles]
		# num = len( dofs_per_bundle )
		if self.system_symbolic_factored is None:
			#print 'even symbolic factoring'
			system = self.to_system_solve_t( self.system )
			self.system_symbolic_factored = self.compute_symbolic_factorization( system )
			self.system_factored = self.system_symbolic_factored( system )
		
		elif self.system_factored is None:
			#print 'even numeric factoring'
			system = self.to_system_solve_t( self.system )
			self.system_factored = self.system_symbolic_factored( system )
		
		#print 'even solve'
		x = self.system_factored( self.rhs )
		# x = linalg.solve( self.system, self.rhs )
		# x = scipy.sparse.linalg.spsolve( self.system, self.rhs )
	
		### Return a nicely formatted chain of bezier curves.
		result = []
		dofs_offset = 0
		for dofs, dirs in zip(dofs_per_bundle, dirs_per_bundle):
			solution = zeros( (4, 2) )
			if dofs[0] == 4:
				solution[ :2 ] = x[dofs_offset : dofs_offset + 4].reshape(2,2)
				
			elif dofs[0] == 3:
				solution[0] = x[dofs_offset : dofs_offset + 2]
				solution[1] = solution[0] + asarray(dirs[0]) * x[dofs_offset + 2]
			
			dofs_offset += dofs[0]
			
			if dofs[1] == 4:
				solution[ 2:4 ] = x[dofs_offset : dofs_offset + 4].reshape(2,2)
				
			elif dofs[1] == 3:
				solution[3] = x[dofs_offset : dofs_offset + 2]
				solution[2] = solution[3] + asarray(dirs[1]) * x[dofs_offset + 2]
				

			dofs_offset += dofs[1]
				
			result.append(solution)
		
		return result	
	
	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1, angle ):
		mag0, mag1 = bundle0.magnitudes[1], bundle1.magnitudes[0]
		dim = 2
		dofs0 = self.compute_dofs_per_curve(bundle0)
		dofs1 = self.compute_dofs_per_curve(bundle1)
		dofs = sum(dofs0) + sum(dofs1)
		dirs0 = asarray(bundle0.directions)
		dirs1 = asarray(bundle1.directions)
		
		assert bundle0.constraints[1][0] == bundle1.constraints[0][0]
		smoothness = bundle0.constraints[1][0]
		if smoothness == 'C0':         ## C0
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			'''
			R = zeros( ( dofs, dim ) )
			if dofs0[1] == 3:
				R[dofs0[0] : dofs0[0]+dim, :] = identity(dim)
			else:
				R[sum(dofs0)-dim : sum(dofs0), :] = identity(dim)
			R[sum(dofs0) : sum(dofs0)+dim, :] = identity(dim) * -1

		elif smoothness == 'A':        ## fixed angle
			R = zeros( ( dofs, dim ) )
			if dofs0[1] == 3:
				R[dofs0[0] : dofs0[0]+dim, :] = identity(dim)
			else:
				R[sum(dofs0)-dim : sum(dofs0), :] = identity(dim)
			R[sum(dofs0) : sum(dofs0)+dim, :] = identity(dim) * -1
			
		elif smoothness == 'C1':        ## C1
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			lambda3 * ( w_q(P4x' - P3x') + w_p(Q1x' - Q2x')) = 0
			lambda4 * ( w_q(P4y' - P3y') + w_p(Q1y' - Q2y')) = 0
			'''
			R = zeros( ( dofs, 2*dim ) )
			## about bundle0
			if dofs0[1] == 3:
				R[ dofs0[0] : dofs0[0]+dim, :dim ] = identity( dim )
				R[ sum(dofs0)-1, -dim: ] = -mag1*dirs0[1]
			else:
				R[ sum(dofs0)-dim : sum(dofs0), :dim ] = identity( dim )
				R[ sum(dofs0)-dim : sum(dofs0), -dim: ] = identity( dim ) * mag1
				R[ dofs0[0] : dofs0[0]+dim, -dim: ] = -identity( dim ) * mag1
			## about bundle1
			R[ sum(dofs0) : sum(dofs0)+dim, :dim ] = -identity( dim )
			R[ sum(dofs0) : sum(dofs0)+dim, -dim: ] = identity( dim ) * mag0
			if dofs1[0] == 4:
				R[ sum(dofs0)+dofs1[0]-dim : sum(dofs0)+dofs1[0], -dim:] = -mag0 * identity( dim )
			elif dofs1[0] == 3:
				R[ sum(dofs0)+dofs1[0]-1, -dim: ] = -mag0 * dirs1[0]
				
				
		elif smoothness == 'G1':        ## G1
			R = zeros( ( dofs, dim ) )
			if dofs0[1] == 3:
				R[dofs0[0] : dofs0[0]+dim, :] = identity(dim)
			else:
				R[sum(dofs0)-dim : sum(dofs0), :] = identity(dim)
			R[sum(dofs0) : sum(dofs0)+dim, :] = identity(dim) * -1
		else:
			R = zeros( ( dofs, 0 ) )
		
		rhs = zeros(R.shape[1])
		
		fixed = bundle0.control_points[-1][:dim]
		is_fixed = bundle0.constraints[1][1]
		assert type( is_fixed ) == bool
		if is_fixed:

			fixed = asarray(fixed)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			R2 = zeros( ( dofs, dim ) )
			R2[sum(dofs0)-dim : sum(dofs0), :] = identity(dim)
		
			R = concatenate((R, R2), axis=1)
			rhs = concatenate((rhs, fixed))
	
		return R.T, rhs
		
		
	def system_for_curve( self, bundle ):
		dofs = self.compute_dofs_per_curve(bundle)
		dirs = asarray(bundle.directions)
		length = bundle.length

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
					
		else:
			raise RuntimeError('bundle return wrong dofs.')
		
		return asarray( Left )*length
		
	
	def system_for_curve_with_arc_length( self, bundle ):
		'''
		## Solve the same integral as system__for_curve only with dt replaced by ds
		'''
		ts = bundle.ts
  		dts = bundle.dts
		dim = 2
		
		dofs = self.compute_dofs_per_curve(bundle)
		dirs = asarray(bundle.directions)
		length = bundle.length
		
		if array_equal(dofs, (4,4)):
			
			tbar = ones( ( 4, 1 ) )
			MAM = zeros( ( 4, 4 ) )
			for i in range(len(dts)):
				t = (ts[i] + ts[i+1])/2
				ds = dts[i]
		
				tbar[0] = t**3
				tbar[1] = t**2
				tbar[2] = t
		
				Mtbar = dot( M.T, tbar )
				MAM += dot( Mtbar, Mtbar.T )*ds
			
			Left = zeros( ( 8, 8 ) )	
			Left[ : : 2, : : 2 ] = MAM
			Left[ 1: : 2, 1: : 2 ] = MAM
			print 'arc length 4*4', Left	

		## p1x, p1y, s, p3x, p3y, p4x, p4y
		elif array_equal(dofs, (3,4)):
			Left = zeros( ( 7,  7 ) )
			
			coefs = zeros( 4 )
			for i in range(len(dts)):
				temp = zeros( ( 7,  7 ) )
				t = (ts[i] + ts[i+1])/2
				ds = dts[i]
		
				coefs[0] = 2*t**3 - 3*t**2 + 1
				coefs[1] = 3*( t**3 - 2*t**2 + t )
				coefs[2] = -3*( t**3 - t**2 )
				coefs[3] = t**3
			
				## top half
				temp[0,0] = temp[1,1] = coefs[0]
				temp[0,2] = coefs[1] * dirs[0,0]
				temp[1,2] = coefs[1] * dirs[0,1]
				temp[0,3] = temp[1,4] = coefs[2]
				temp[0,5] = temp[1,6] = coefs[3]			
				temp[2] = temp[0] * coefs[1] * dirs[0,0] + temp[1] * coefs[1] * dirs[0,1]
				temp[0] = temp[0] * coefs[0]
				temp[1] = temp[1] * coefs[0]
				## bottom half
				temp[3,0] = temp[4,1] = temp[5,0] = temp[6,1] = coefs[0]
				temp[3,2] = temp[5,2] = coefs[1] * dirs[0,0]
				temp[4,2] = temp[6,2] = coefs[1] * dirs[0,1]
				temp[3,3] = temp[4,4] = temp[5,3] = temp[6,4] = coefs[2]				
				temp[3,5] = temp[4,6] = temp[5,5] = temp[6,6] = coefs[3]
				temp[3:5] *= coefs[2]
				temp[5:7] *= coefs[3]
				
				Left = Left + temp * ds

						
		## p1x, p1y, p2x, p2y, p4x, p4y, u
		elif array_equal(dofs, (4,3)):
			Left = zeros( ( 7,  7 ) )
			
			coefs = zeros( 4 )
			for i in range(len(dts)):
				temp = zeros( ( 7,  7 ) )
				t = (ts[i] + ts[i+1])/2
				ds = dts[i]
		
				coefs[0] = t**3 - 3*t**2 + 3*t - 1
				coefs[1] = -3*( t**3 - 2*t**2 + t )
				coefs[2] = 2*t**3 - 3*t**2
				coefs[3] = 3*( t**3 - t**2 )
			
				## top half
				temp[0,0] = temp[1,1] = temp[2,0] = temp[3,1] = coefs[0]
				temp[0,2] = temp[1,3] = temp[2,2] = temp[3,3] = coefs[1]				
				temp[0,4] = temp[1,5] = temp[2,4] = temp[3,5] = coefs[2]
				temp[0,6] = temp[2,6] = coefs[3] * dirs[1,0]
				temp[1,6] = temp[3,6] = coefs[3] * dirs[1,1]
				temp[0:2] *= coefs[0]
				temp[2:4] *= coefs[1]
				## bottom half
				temp[4,0] = temp[5,1] = coefs[0]
				temp[4,2] = temp[5,3] = coefs[1]
				temp[4,4] = temp[5,5] = coefs[2]
				temp[4,6] = coefs[3] * dirs[1,0]
				temp[5,6] = coefs[3] * dirs[1,1]		
				temp[6] = temp[4] * coefs[3] * dirs[1,0] + temp[5] * coefs[3] * dirs[1,1]
				temp[4] = temp[4] * coefs[2]
				temp[5] = temp[5] * coefs[2]
				
				Left = Left + temp * ds
					
		## p1x, p1y, s, p4x, p4y, u
		elif array_equal(dofs, (3,3)):
			Left = zeros( ( 6,  6 ) )
			
			coefs = zeros( 4 )
			for i in range(len(dts)):
				temp = zeros( ( 6,  6 ) )
				t = (ts[i] + ts[i+1])/2
				ds = dts[i]
		
				coefs[0] = -( 2*t**3 - 3*t**2 + 1 )
				coefs[1] = 2*t**3 - 3*t**2
				coefs[2] = -3*( t**3 - 2*t**2 + t )
				coefs[3] = 3*( t**3 - t**2 )
			
				## top half
				temp[0,0] = temp[1,1] = coefs[0]
				temp[0,2] = coefs[2] * dirs[0,0]
				temp[1,2] = coefs[2] * dirs[0,1]
				temp[0,3] = temp[1,4] = coefs[1]
				temp[0,5] = coefs[3] * dirs[1,0]
				temp[1,5] = coefs[3] * dirs[1,1]
				temp[2] = temp[0] * coefs[2] * dirs[0,0] + temp[1] * coefs[2] * dirs[0,1]
				temp[0] = temp[0] * coefs[0]
				temp[1] = temp[1] * coefs[0]
				## bottom half
				temp[3,0] = temp[4,1] = coefs[0]
				temp[3,2] = coefs[2] * dirs[0,0]
				temp[4,2] = coefs[2] * dirs[0,1]
				temp[3,3] = temp[4,4] = coefs[1]
				temp[3,5] = coefs[3] * dirs[1,0]
				temp[4,5] = coefs[3] * dirs[1,1]
				temp[5] = temp[3] * coefs[3] * dirs[1,0] + temp[4] * coefs[3] * dirs[1,1]
				temp[3] = temp[3] * coefs[1]
				temp[4] = temp[4] * coefs[1]
				
				Left = Left + temp * ds
					
		else:
			raise RuntimeError('bundle return wrong dofs.')
		
		return Left*length
		
	def compute_dofs_per_curve( self, bundle ):
		dofs = zeros( 2, dtype = int )
		'''
		assume open end points can only emerge at the endpoints
		'''
		for i, (smoothness, fixed) in enumerate(bundle.constraints):
			
			if smoothness == 'C0': dofs[i] += 4			## C0
			elif smoothness == 'A': dofs[i] += 3		## fixed angle
			elif smoothness == 'C1': dofs[i] += 4		## C1
			elif smoothness == 'G1': dofs[i] += 3		## G1
			elif smoothness == 'None': dofs[i] += 4  	## Free of constraint
			
		return dofs

	
	def constraint_number_per_joint(self, constraint ):	   
		assert len(constraint) == 2
		smoothness = constraint[0]  
		is_fixed = constraint[1]    

		num = 0
		if smoothness == 'C0': num = 2         ## C0
		elif smoothness == 'A': num = 2       ## fixed angle
		elif smoothness == 'C1': num = 4       ## C1
		elif smoothness == 'G1': num = 2       ## G1
		
		assert type( is_fixed ) == bool
		if is_fixed:
			num += 2
	
		return num	
		  
	def rhs_for_curve( self, bundle, transforms ):
		'''
		The rhs is computed according to the formula:
			R = sum(wi * T_i * P.T * M * tbar )
			rhs(P1) = integral(R * M_hat, dt)
		'''
		length = bundle.length
		W_matrices = bundle.W_matrices
		dofs = self.compute_dofs_per_curve(bundle)
		Right = zeros( sum(dofs) )
		
		controls = asarray(bundle.control_points)
		dirs = asarray(bundle.directions)
		temp = zeros( (3, 4) )
		
		for i in range( len( transforms ) ):

			T_i = mat( asarray(transforms[i]).reshape(3, 3) )
 			W_i = asarray(W_matrices[i])
			
			temp = temp + dot(asarray(T_i*(controls.T)*M), W_i)
		R = temp[:2,:]
		
		## p1x, p1y, p2x, p2y, p3x, p3y, p4x, p4y
		if array_equal(dofs, (4,4)):
			Right[:] = concatenate((R[:,0], R[:,1], R[:,2], R[:,3])) 
		## p1x, p1y, s, p3x, p3y, p4x, p4y
		elif array_equal(dofs, (3,4)):
			Right[:2] = R[:,0] + R[:,1]
			Right[2] = dot( R[:,1], dirs[0] )
			Right[-4:] = concatenate((R[:,2], R[:,3]))
		## p1x, p1y, p2x, p2y, p4x, p4y, u
		elif array_equal(dofs, (4,3)):
			Right[:4] = concatenate((R[:,0], R[:,1]))
			Right[-3:-1] = R[:,2] + R[:,3]
			Right[-1] = dot( R[:,2], dirs[1] )
		## p1x, p1y, s, p4x, p4y, u
		elif array_equal(dofs, (3,3)):
			Right[:2] = R[:,0] + R[:,1]
			Right[2] = dot( R[:,1], dirs[0] )
			Right[-3:-1] = R[:,2] + R[:,3]
			Right[-1] = dot( R[:,2], dirs[1] )
		else:
			raise RuntimeError('bundle return wrong dofs.')
		
		return Right*length
