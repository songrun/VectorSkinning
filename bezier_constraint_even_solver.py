from generate_chain_system import *

class BezierConstraintSolverEven( BezierConstraintSolver ):
	'''
	Fixed direction, magnitude free (for G1 or A).
	'''
	
	def update_system_with_result_of_previous_iteration( self, solution ):
		### Iterate only over the parts of the matrix that will change,
		### such as the system matrix involving fixed directions,
		### lagrange multipliers across G1 or A edges, and the right-hand-side.
		
		raise NotImplementedError
	
	def solve( self ):
		### Return a nicely formatted chain of bezier curves, un-substituting the fixed
		### directions.
		raise NotImplementedError
	
	def lagrange_equations_for_curve_constraints( self, bundle0, bundle1 ):
		raise NotImplementedError
	def system_for_curve( self, bundle ):
		raise NotImplementedError
	def constraint_equations_per_bundle( self, constraint ):
		raise NotImplementedError
		
	def compute_dofs_per_curve( self, constraint ):
		constraints = asarray(bundle.constraints)
		dofs = zeros(2)
		'''
		assume open end points can only emerge at the endpoints
		'''
		for i, smoothness in enumerate(constraints[:,0]):
			if smoothness == 1: dofs[i] += 4		## C0
			elif smoothness == 2: dofs[i] += 3		## fixed angle
			elif smoothness == 3: dofs[i] += 4		## C1
			elif smoothness == 4: dofs[i] += 3		## G1
		
		return dofs
		return dofs
	
	def constraint_number_per_joint(self, constraint ):	   
		assert len(constraint) == 2
		smoothness = constraint[0]	
		is_fixed = constraint[1]	
		
		num = 0
		if smoothness == 1: num = 2			## C0
		elif smoothness == 2: num = 2		## fixed angle
		elif smoothness == 3: num = 4		## C1
		elif smoothness == 4: num = 2		## G1
		
		if is_fixed == 1:
			num += 2
			
		return num		
		  
	def rhs_for_curve( self, bundle, handles ):
		raise NotImplementedError