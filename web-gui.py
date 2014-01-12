#!/opt/local/bin/python

## For anything
from twisted.internet import reactor

## For WebSockets
from autobahn.websocket import WebSocketServerFactory, WebSocketServerProtocol, listenWS

## All payloads are JSON-formatted.
import json
from chain_computer import *
from tictoc import tic, toc, tictoc_dec

kVerbose = 2

class WebGUIServerProtocol( WebSocketServerProtocol ):
	def connectionMade( self ):
		WebSocketServerProtocol.connectionMade( self )
		self.engine = self.factory.engine
		
		print 'CONNECTED'
	
	@tictoc_dec
	def onMessage( self, msg, binary ):
		### BEGIN DEBUGGING
		if kVerbose >= 2:
			if not binary:
				from pprint import pprint
				space = msg.find( ' ' )
				print msg[ :space ]
				pprint( json.loads( msg[ space+1 : ] ) )
		elif kVerbose >= 1:
			if not binary:
				print( msg[:72] ), '...'
		### END DEBUGGING
		
		if binary:
			print 'Received unknown message: binary of length', len( msg )
		
		elif msg.startswith( 'paths-info ' ):	   
			paths_info = json.loads( msg[ len( 'paths-info ' ): ] )
			boundary_path = max(paths_info, key=lambda e : e[u'bbox_area'])
			boundary_index = paths_info.index( boundary_path )
			
			self.engine.set_control_positions( paths_info, boundary_index )
			all_constraints = self.engine.all_constraints

			for i, constraints in enumerate( all_constraints ):
				for j, constraint in enumerate( constraints ):
	
					continuity = constraint[0]
					fixed = constraint[1]
					
					payload = [ i, j, { 'fixed': fixed, 'continuity': continuity} ]
					self.sendMessage( 'control-point-constraint ' + json.dumps( payload ) )

		
		elif msg.startswith( 'handle-positions ' ):
			paths_info = json.loads( msg[ len( 'handle-positions ' ): ] )

			self.engine.set_handle_positions( paths_info )
			self.engine.precompute_configuration()
			
			all_paths = self.engine.solve()	

			all_positions = make_chain_from_control_groups( all_paths )
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Generate the triangulation and the BBW weights.
			# self.engine ...
		
		elif msg.startswith( 'handle-transforms ' ):
			handle_transforms = json.loads( msg[ len( 'handle-transforms ' ): ] )
			
			tic( 'transform_change' )
			for handle_index, handle_transform in handle_transforms:
				self.engine.transform_change( handle_index, handle_transform )
			toc()
			
			tic( 'engine.solve()' )
			all_paths = self.engine.solve_transform_change()	
# 			debugger()
			toc()

			tic( 'make_chain_from_control_groups' )
			all_positions = make_chain_from_control_groups( all_paths )
			print 'paths-positions ', all_positions
			toc()
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Solve for the new curve positions given the updated transform matrix.
			# new_positions = engine ...
			
			## Send the new positions to the GUI.
			# self.sendMessage( 'paths-positions ' + json.dumps( new_positions ) )
		
		elif msg.startswith( 'control-point-constraint ' ):
			paths_info = json.loads( msg[ len( 'control-point-constraint ' ): ] )
			
			constraint = [None]*2
			constraint[0] = str( paths_info[2][ u'continuity' ] )
			constraint[1] = paths_info[w][ u'fixed' ]
				
			self.engine.constraint_change( paths_info[0], paths_info[1], constraint )
			
			all_paths = self.engine.solve()	

			all_positions = make_chain_from_control_groups( all_paths )
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Solve for the new curve positions given the updated control point constraint.
			# new_positions = self.engine ...
			
			## Send the new positions to the GUI.
			# self.sendMessage( json.dumps( new_positions ) )
		
		else:
			print 'Received unknown message:', msg

def make_chain_from_control_groups( all_paths ):
	
	all_positions = [] 	
	for path in all_paths:
		if len( path ) > 1:
			new_positions = concatenate( asarray(path)[:-1, :-1] )
			new_positions = concatenate( ( new_positions, path[-1] ) )
		else:
			new_positions = path[0]
		new_positions = new_positions.tolist()
		all_positions.append( new_positions )
		
	return all_positions


def setupWebSocket( address, engine ):
	'''
	Listen for WebSocket connections at the given address.
	'''
	
	factory = WebSocketServerFactory( address )
	factory.engine = engine
	factory.protocol = WebGUIServerProtocol
	listenWS( factory )
	
	print "Listening for WebSocket connections at:", address

if __name__ == '__main__':
	import os, sys
	
	## Create engine:
	# engine = ...
	engine = Engine()
	
	setupWebSocket( "ws://localhost:9123", engine )
	
	## Maybe you find this convenient
	if sys.argv[1:2] == ['open']:
		os.system( 'open web-gui.html' )
	
	reactor.run()
