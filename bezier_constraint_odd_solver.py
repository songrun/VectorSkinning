from generate_chain_system import *

import scipy.sparse.linalg

class BezierConstraintSolverOdd( BezierConstraintSolver ):
	'''
	Free direction, magnitude fixed (for G1 or A).
	'''
	
	def update_system_with_result_of_previous_iteration( self, solution ):
		### Iterate only over the parts of the matrix that will change,
		### such as the lagrange multipliers across G1 or A edges and the right-hand-side.
		solution = asarray(solution)
		num = len(self.bundles)
		assert solution.shape == (num, 4, 2)
		magnitudes = [[mag( solution[i][1]-solution[i][0] ), mag( solution[i][2]-solution[i][3] )] for i in range(num) ]
		
		for i in range(num):
			self.bundles[i].magnitudes = magnitudes[i]
		
		self._update_bundles( lagrange_only = False )
		## The lagrange multipliers changed, but not the locations of the zeros.
		self.system_factorization = None
		## UPDATE: Actually, if constrained directions align with coordinate axes
		##         or have zero magnitude, then the lagrange multiplier systems may gain
		##		   or lose zeros.
		## UPDATE 2: If we could update_bundles once with all fixed angles
		##           not parallel or perpendicular, and then compute the symbolic
		##           factorization, we could keep it.
		self.system_symbolic_factorization = None
		
		
	
	def solve( self ):
		dim = 2
		num = len(self.bundles)
#		debugger()
		
		if self.system_symbolic_factorization is None:
			self.system_symbolic_factorization = compute_symbolic_factorization( self.system )
		if self.system_factored is None:
			self.system_factored = self.system_symbolic_factorization( self.system )
		x = self.system_factored( self.rhs )
		# x = linalg.solve( self.system, self.rhs )
		# x = scipy.sparse.linalg.spsolve( self.system, self.rhs )
		### Return a nicely formatted chain of bezier curves.
		x = array( x[:self.total_dofs] ).reshape(-1,4).T

		result = []
		for i in range(num):
			P = x[:, i*dim:(i+1)*dim ]
			result.append( P )		
		
		return result	
	
	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1 ):
		mag0, mag1 = bundle0.magnitudes[1], bundle1.magnitudes[0]
		
		vec0 = (bundle0.control_points[2]-bundle0.control_points[3])[:2]
		vec1 = (bundle1.control_points[1]-bundle1.control_points[0])[:2]
		cos_theta = dot(vec0, vec1)/( mag(vec0)*mag(vec1) )
		sin_theta = (1.-cos_theta**2) ** 0.5
		
		dim = 2
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
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0)+i*4, i] = -1
				
				R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				
				## tell the angle from vec0 to vec1 is positive or negative.
				if cross(vec0, vec1) >= 0:
					R[sum(dofs0):sum(dofs0)+dim, dim:] = asarray([[-cos_theta, sin_theta], [cos_theta, -sin_theta]])
					R[-dim*2:-dim, dim:] = asarray([[-sin_theta, -cos_theta], [sin_theta, cos_theta]])
				else:
					R[sum(dofs0):sum(dofs0)+dim, dim:] = asarray([[-cos_theta, -sin_theta], [cos_theta, sin_theta]])
					R[-dim*2:-dim, dim:] = asarray([[sin_theta, -cos_theta], [-sin_theta, cos_theta]])
			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= mag1
			R[ sum(dofs0):, dim: ] *= mag0
			
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
		
		fixed = bundle0.control_points[-1][:2]
		is_fixed = bundle0.constraints[1][1]
		if is_fixed == True or is_fixed == 'True':

			fixed = asarray(fixed)
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - constraint_X' ) = 0
			lambda2 * ( P4y' - constraint_Y' ) = 0
			'''
			R2 = zeros( ( dofs, dim ) )
			for i in range( dim ):
				R2[i*4+3, i] = 1
		
			R = concatenate((R, R2), axis=1)
			rhs = concatenate((rhs, fixed))
	
		return R.T, rhs
		
		
	def system_for_curve( self, bundle ):
		'''
		## A is computed using Sage, integral of (tbar.T * tbar) with respect to t.
		#	A = asarray( [[	 1./7,	1./6,  1./5, 1./4], [ 1./6,	 1./5, 1./4,  1./3], 
		#		[ 1./5, 1./4,  1./3,  1./2], [1./4,	 1./3,	1./2,	1.]] )
		## MAM is computed using Sage. MAM = M * A * M
		'''
		MAM = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,	3./35, 9./140,	1./35], [ 1./35, 9./140,  3./35,  1./14], [1./140,	1./35,	1./14,	 1./7]] )
		dim = 2

		Left = zeros((8, 8))

		for i in range(dim):		
			Left[ i*4:(i+1)*4, i*4:(i+1)*4 ] = MAM[:,:]
		
		
		return Left
		
		
	def compute_dofs_per_curve( self, bundle ):
	
		constraints = asarray(bundle.constraints)
		dofs = zeros( 2, dtype = int )
		'''
		assume open end points can only emerge at the endpoints
		'''
		for i, smoothness in enumerate(constraints[:,0]):
			if smoothness == 'C0': dofs[i] += 4			## C0
			elif smoothness == 'A': dofs[i] += 4		## fixed angle
			elif smoothness == 'C1': dofs[i] += 4		## C1
			elif smoothness == 'G1': dofs[i] += 4		## G1
			elif smoothness == 'None': dofs[i] += 4		## Free of constraint
			
		return dofs
		
	
	def constraint_number_per_joint(self, constraint ):	   
	
		assert len(constraint) == 2
		smoothness = constraint[0]	
		fixed = constraint[1]	 
		
		num = 0
		if smoothness == 'C0': num = 2			## C0
		elif smoothness == 'A': num = 4		## fixed angle
		elif smoothness == 'C1': num = 4		## C1
		elif smoothness == 'G1': num = 4		## G1
		
		if fixed != 0:
			num += 2
			
		return num
		
		
	def rhs_for_curve( self, bundle, transforms ):
		'''
		The rhs is computed according to the formula:
			rhs = sum(Ti * P.T * M.T * W_i * M)
		'''
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
		
		Right = zeros( 8 )
		temp = zeros( (3, 4) )
		for i in range( len( transforms ) ):
		
			T_i = mat( asarray(transforms[i]).reshape(3, 3) )
			W_i = asarray(W_matrices[i])

			temp = temp + dot(asarray(T_i*(controls.T)*M), W_i)

		R = temp[:2,:]
		
		Right[:] = concatenate((R[0, :], R[1, :]))
			
		return Right