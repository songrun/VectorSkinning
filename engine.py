from chain_computer import *

class Handle:
	'''
	Handle is a class of the handle applied to a shape. It has transform and position.
	'''
	id = -1
	position = None
	transform = identity(3)
	
	def __init__( self, position, id=-1, transform=identity(3) ):
		position = asarray( position )
		assert position == (2,)
		
		self.position = position
		self.id = id
		self.transform = transform
		
	def compare_shape( compared_handle ):
		assert isinstance( compared_handle, Handle )
		if compared_handle.id == self.id and array_equal( compared_handle.position, self.position):
			return True
		else: 
			return False
		
class Precomputed_parameters:
	all_vertices = None
	all_weights = None
	all_indices = None
	all_pts = None
	all_dts = None
	W_matrices = None
	
	def __init__(W_matrices, all_weights, all_vertices, all_indices, all_pts, all_dts):
		self.W_matrices = W_matrices
		self.all_weights = all_weights
		self.all_vertices = all_weights
		self.all_indices = all_indices
		self.all_pts = all_pts
		self.all_dts = all_dts

class Engine:

	transforms = []
	controls = []	
	precomputed_parameter_table = []
	is_ready = False
	
	def __init__( controls, handles ):
		setup_configuration( controls, handles )
		self.is_ready = True
	
	
	def setup_configuration( controls, handles ):
		
			
		