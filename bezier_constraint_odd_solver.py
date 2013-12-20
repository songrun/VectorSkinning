from generate_chain_system import *

class BezierConstraintSolverOdd( BezierConstraintSolver ):
    '''
    Free direction, magnitude fixed (for G1 or A).
    '''
    
    def update_system_with_result_of_previous_iteration( self, solution ):
        ### Iterate only over the parts of the matrix that will change,
        ### such as the lagrange multipliers across G1 or A edges and the right-hand-side.
        
        raise NotImplementedError
    
    def solve( self ):
		dim = 2
		num = len(self.bundles)
		
		x = linalg.solve( self.system, self.rhs )
		### Return a nicely formatted chain of bezier curves.
		x = array( x[:self.total_dofs] ).reshape(-1,4).T

		result = []
		for i in range(num):
			P = x[:, i*dim:(i+1)*dim ]
			result.append( P )		

		return result	
    
    def lagrange_equations_for_curve_constraints( self, bundle0, bundle1 ):

		weight0, weight1 = bundle0.weight, bundle1.weight
		dim = 2
		dofs0 = self.compute_dofs_per_curve(bundle0)
		dofs1 = self.compute_dofs_per_curve(bundle1)
		dofs = sum(dofs0) + sum(dofs1)
		
		R = None
		
		smoothness = bundle0.constraints[1,0]
		if smoothness == 1:         ## C0
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			'''
			R = zeros( ( dofs, dim ) )
			for i in range( dim ):
				R[i*4+3, i] = 1
				R[sum(dofs0) + i*4, i] = -1

		elif smoothness == 2:        ## fixed angle
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				R[sum(dofs0)+i*4, i] = R[sum(dofs0)+i*4+1, i+dim] = -1
				R[sum(dofs0)+i*4, i+dim] = 1

			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= weight1
			R[ sum(dofs0):, dim: ] *= weight0
			
		elif smoothness == 3:        ## C1
			'''
			Boundary Conditions are as follows:
			lambda1 * ( P4x' - Q1x' ) = 0
			lambda2 * ( P4y' - Q1y' ) = 0
			lambda4 * ( P4x' - P3x' + Q1x' - Q2x') = 0
			lambda5 * ( P4y' - P3y' + Q1y' - Q2y') = 0
			'''
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				R[sum(dofs0)+i*4, i] = R[sum(dofs0)+i*4+1, i+dim] = -1
				R[sum(dofs0)+i*4, i+dim] = 1

			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= weight1
			R[ sum(dofs0):, dim: ] *= weight0

		elif smoothness == 4:        ## G1
			R = zeros( ( dofs, 2*dim ) )
			for i in range( dim ):
				R[i*4+3, i] = R[i*4+3, i+dim] = 1
				R[i*4+2, i+dim] = -1
				R[sum(dofs0)+i*4, i] = R[sum(dofs0)+i*4+1, i+dim] = -1
				R[sum(dofs0)+i*4, i+dim] = 1

			## add weights to lambda	 
			R[ :sum(dofs0), dim: ] *= weight1
			R[ sum(dofs0):, dim: ] *= weight0
		
		rhs = zeros(R.shape[1])
		
		fixed = bundle0.control_points[-1][:2]
		
		if bundle0.constraints[1, 1] != 0:

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
		# 	A = asarray( [[  1./7,  1./6,  1./5, 1./4], [ 1./6,  1./5, 1./4,  1./3], 
		# 		[ 1./5, 1./4,  1./3,  1./2], [1./4,  1./3,  1./2,   1.]] )
		## MAM is computed using Sage. MAM = M * A * M
		'''
		MAM = asarray( [[  1./7,  1./14,  1./35, 1./140], [ 1./14,  3./35, 9./140,  1./35], 
			[ 1./35, 9./140,  3./35,  1./14], [1./140,  1./35,  1./14,   1./7]] )
		dim = 2
		weight = bundle.weight

		Left = zeros((8, 8))

		for i in range(dim):		
			Left[ i*4:(i+1)*4, i*4:(i+1)*4 ] = MAM[:,:]*weight

		return Left
		
        
    def compute_dofs_per_curve( self, bundle ):
    
        constraints = asarray(bundle.constraints)
        dofs = zeros(2)
        '''
        assume open end points can only emerge at the endpoints
        '''
        for i, smoothness in enumerate(constraints[:,0]):
            if smoothness == 1: dofs[i] += 4        ## C0
            elif smoothness == 2: dofs[i] += 4      ## fixed angle
            elif smoothness == 3: dofs[i] += 4      ## C1
            elif smoothness == 4: dofs[i] += 4      ## G1
        
        return dofs
        
    
    def constraint_number_per_joint(self, constraint ):    
    
        assert len(constraint) == 2
        smoothness = constraint[0]  
        fixed = constraint[1]    
        
        num = 0
        if smoothness == 1: num = 2         ## C0
        elif smoothness == 2: num = 4       ## fixed angle
        elif smoothness == 3: num = 4       ## C1
        elif smoothness == 4: num = 4       ## G1
        
        if fixed != 0:
            num += 2
            
        return num
        
    def rhs_for_curve( self, bundle, transforms ):
		'''
		The rhs is computed according to the formula:
			rhs = sum(Ti * P.T * M.T * W_i * M)
		'''
		dim = 3
		Ws = bundle.Ws
		controls = bundle.control_points
		weight = bundle.weight

		Right = zeros( (dim, 4) )
		for i in range( len( transforms ) ):

			T_i = mat( asarray(transforms[i]).reshape(dim,dim) )
			W_i = Ws[i,0]	

			Right = Right + T_i * (controls.T) * M * mat( W_i ) * M

		Right = asarray(Right).reshape(-1)*weight
		Right = Right[:8]	

		return Right