#!/opt/local/bin/python

## For anything
from twisted.internet import reactor

## For WebSockets
from autobahn.websocket import WebSocketServerFactory, WebSocketServerProtocol, listenWS

## All payloads are JSON-formatted.
import json
from chain_computer import *

engine = None

class WebGUIServerProtocol( WebSocketServerProtocol ):
	def connectionMade( self ):
		WebSocketServerProtocol.connectionMade( self )
		self.engine = self.factory.engine
		
		print 'CONNECTED'
	
	def onMessage( self, msg, binary ):
		### BEGIN DEBUGGING
# 		if not binary:
# 			from pprint import pprint
# 			pprint( json.loads( msg[ msg.find( ' ' )+1 : ] ) )
		### END DEBUGGING
		
		if binary:
			print 'Received unknown message: binary of length', len( msg )
		
		elif msg.startswith( 'paths-info ' ):	   
			paths_info = json.loads( msg[ len( 'paths-info ' ): ] )
			boundary_path = max(paths_info, key=lambda e : e[u'bbox_area'])
			boundary_index = paths_info.index( boundary_path )
			
# 			test_pushback = [ path[u'cubic_bezier_chain'] for path in paths_info ]
# 			self.sendMessage( 'paths-positions ' + json.dumps( test_pushback ) )
# 			debugger()
			self.engine.set_control_positions( paths_info, boundary_index )
			all_constraints = self.engine.all_constraints
# 			print 'constraints: ', all_constraints
			for i, constraints in enumerate( all_constraints ):
				for j, constraint in enumerate( constraints ):
					fixed = False
					if constraint[1] == 1:	fixed = True
					continuity = 'C0'
					if constraint[0] == 1: continuity = 'C0'
					elif constraint[0] == 2: continuity = 'A'
					elif constraint[0] == 3: continuity = 'C1'
					elif constraint[0] == 4: continuity = 'G1'
					
					payload = [ i, j, { 'fixed': fixed, 'continuity': continuity} ]
					self.sendMessage( 'control-point-constraint ' + json.dumps( payload ) )
			## Precompute something with the paths.
			# self.engine ...
			
			## Report back assumed constraints.
			# for constraint in constraints:
			# payload = [ path_index, segment_index, { 'fixed': True|False, 'continuity': 'C0|C1|G1|A' } ]
			# self.sendMessage( 'control-point-constraint ' + json.dumps( payload ) )
		
		elif msg.startswith( 'handle-positions ' ):
			paths_info = json.loads( msg[ len( 'handle-positions ' ): ] )

			self.engine.set_handle_positions( paths_info )
			engine.precompute_configuration()
			
			all_paths = self.engine.solve()	
# 			if type(all_paths) in types.StringTypes():
# 				print all_paths
# 				return

			all_positions = [] 					
			for path in all_paths:
				new_positions = concatenate( asarray(path)[:-1, :-1] )
				new_positions = concatenate( ( new_positions, path[-1] ) )
				new_positions = new_positions.tolist()
				all_positions.append( new_positions )
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Generate the triangulation and the BBW weights.
			# self.engine ...
		
		elif msg.startswith( 'handle-transform ' ):
			paths_info = json.loads( msg[ len( 'handle-transform ' ): ] )
			
			self.engine.transform_change( paths_info[0], paths_info[1] )
			
			all_paths = self.engine.solve()	
# 			if type(all_paths) in types.StringTypes():
# 				print all_paths
# 				return

			all_positions = [] 					
			for path in all_paths:
				new_positions = concatenate( asarray(path)[:-1, :-1] )
				new_positions = concatenate( ( new_positions, path[-1] ) )
				new_positions = new_positions.tolist()
				all_positions.append( new_positions )
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Solve for the new curve positions given the updated transform matrix.
			# new_positions = engine ...
			
			## Send the new positions to the GUI.
			# self.sendMessage( 'paths-positions ' + json.dumps( new_positions ) )
		
		elif msg.startswith( 'control-point-constraint ' ):
			paths_info = json.loads( msg[ len( 'control-point-constraint ' ): ] )
			constraint = zeros(2)
			if paths_info[2][ u'continuity' ] == u'C0':	constraint[0] = 1
			elif paths_info[2][ u'continuity' ] == u'A':	constraint[0] = 2
			elif paths_info[2][ u'continuity' ] == u'C1':	constraint[0] = 3
			elif paths_info[2][ u'continuity' ] == u'G1':	constraint[0] = 4
			
			if paths_info[w][ u'fixed' ] == True:	constraint[1] = 1
				
			self.engine.constraint_change( paths_info[0], paths_info[1], constraint )
			
			all_paths = self.engine.solve()	
# 			if type(all_paths) in types.StringTypes():
# 				print all_paths
# 				return
			
			all_positions = [] 	
			for path in all_paths:
				new_positions = concatenate( asarray(path)[:-1, :-1] )
				new_positions = concatenate( ( new_positions, path[-1] ) )
				new_positions = new_positions.tolist()
				all_positions.append( new_positions )
			self.sendMessage( 'paths-positions ' + json.dumps( all_positions ) )

			## Solve for the new curve positions given the updated control point constraint.
			# new_positions = self.engine ...
			
			## Send the new positions to the GUI.
			# self.sendMessage( json.dumps( new_positions ) )
		
		else:
			print 'Received unknown message:', msg

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
