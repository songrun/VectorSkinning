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
		



		
			
		