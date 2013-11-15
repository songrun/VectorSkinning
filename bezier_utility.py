from math import *
from numpy import *
from myarray import *

kEps = 1e-7

try:
   from pydb import debugger

   ## Also add an exception hook.
   import pydb, sys
   sys.excepthook = pydb.exception_hook

except ImportError:
   import pdb
   def debugger():
       pdb.set_trace()

## The bezier curve's coefficient matrix
M = matrix('-1. 3. -3. 1.; 3. -6. 3. 0.; -3. 3. 0. 0.; 1. 0. 0. 0.') 

def sample_cubic_bezier_curve_chain( cps, num_samples = 100 ):
	'''
	Given bezier curve chain control points 'cps' in the format ...,
	and a positive integer representing the number of samples per curve in the chain,
	returns a list of pairs ( samples, ts ) as would be returned by sample_cubic_bezier_curve().
	'''
	result = []
	
	assert len( cps ) % 3 == 0 and len( cps ) >= 6
	dim = len( cps[0] )
	
	cps = cps.tolist()
	cps.append( cps[0] )
	cps = asarray( cps ) 

	for i in range( len( cps )/3 ):
		split = cps[i*3: i*3+4]
		samples, ts = sample_cubic_bezier_curve( split, num_samples )
		result.append( ( samples, ts ) )

	## result in the shape of n by num_samples by dim, n is the number of bezier curves, dim is dimensions
	#result = asarray( result )[:,:,:2]

	return result

def sample_cubic_bezier_curve( P, num_samples = 100 ):
	'''
		a 4-by-k numpy.array P containing the positions of the control points as the rows,
		return two lists: sample points of the bezier curve denoted in P, and corresponding t values
	'''
	result = []
	ts = []
	tbar = ones( 4 )
	for t in linspace( 0, 1, num_samples ):
		ts.append( t )
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
 		tbar = tbar.reshape( (4,1) )
 		
 		point = dot( P.T, dot( M.T, tbar ) )
 		result.append( asarray(point).squeeze() )
 		
	return asarray( result ), asarray( ts )
	
def length_of_cubic_bezier_curve( P, num_samples = 100 ):
	'''
		a 4-by-k numpy.array P containing the positions of the control points as the rows,
		return a list of sample points of the bezier curve denoted in P
	'''

	P = asarray( P )
	assert P.shape[0] == 4
	
	samples = []
	tbar = ones( 4 )
	for t in linspace( 0, 1, num_samples ):
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
 		tbar = tbar.reshape( (4,1) )
 		point = dot( P.T, dot( M.T, tbar ) )
 		samples.append( asarray(point).reshape(-1) )
 	
 	samples = asarray(samples)	
 	lengths = [mag(samples[i]-samples[i+1]) for i in range(len(samples)-1)]	
	
	return sum( lengths )
	
def make_control_points_chain( controls, close = True ):
	
	if len( controls ) %3 != 0 or len(controls) < 4:
		print 'bad number of control points.'
		return None
	
	Cset = []
	for i in range( len( controls )//3 -1 ):
 		Cset.append(controls[i*3:i*3+4].tolist())
	if (close):
		last = controls[-3:].tolist() + [controls[0].tolist()]
 		Cset.append(last)
	Cset = asarray( Cset )
	
	return Cset