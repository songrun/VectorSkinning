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

## The bezier curve's coefficient matrix, P(t) = tbar*M*P
M = matrix('-1. 3. -3. 1.; 3. -6. 3. 0.; -3. 3. 0. 0.; 1. 0. 0. 0.') 

def sample_cubic_bezier_curve_chain( Cset, num_samples = 100 ):
	'''
	Given a set of bezier curve chain control points 
	and a positive integer representing the number of samples per curve in the chain,
	returns a list of pairs ( samples, ts ) as would be returned by sample_cubic_bezier_curve().
	
	if curves are open, add samples of the straight line that connect the begin and the end.
	'''
	result = []
	all_dts = []
	
	for P in Cset:
		samples, ts, dts = sample_cubic_bezier_curve_with_dt( P, num_samples )
		result.append( ( samples, ts ) )
		all_dts.append( dts )
		
	if array_equal(Cset[0][0], Cset[-1][-1]) == False:
		samples, ts, dts = sample_cubic_bezier_curve_with_dt( Cset[-1][-1], Cset[0][0], num_samples )	
		result.append( ( samples, ts ) )
		all_dts.append( dts )

	## result in the shape of n by num_samples by dim, n is the number of bezier curves, dim is dimensions
	#result = asarray( result )[:,:,:2]

	return result, all_dts

def sample_cubic_bezier_curve_with_dt( P, num_samples = 100 ):
	'''
	a 4-by-k numpy.array P containing the positions of the control points as the rows,
	return two lists: sample points of the bezier curve denoted in P, and corresponding t values
	'''
	if num_samples is None:
		num_samples = max(int(length_of_cubic_bezier_curve(P) / 1), 2)
	
	result = []
	ts = []
	tbar = ones( (4,1) )
	for t in linspace( 0, 1, num_samples ):
		ts.append( t )
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
 		
 		point = dot( P.T, dot( M.T, tbar ) )
 		result.append( asarray(point).squeeze() )
 	
 	dts = ones( num_samples-1 ) * (1./(num_samples-1) )
 		
	return asarray( result ), asarray( ts ), asarray( dts )
	
	
def sample_cubic_bezier_curve_with_ds( P, num_samples = 100 ):	
	pass
	
def sample_straight_line( begin, end, num_samples = 100 ):

	begin = asarray(begin)
	end = asarray(end)
	result = []
	
	if num_samples is None:
		num_samples = max(int(mag(begin - end) / 1), 2)
	
	ts = []	
	for t in linspace( 0, 1, num_samples ):
		ts.append( t )
		 
		point = (end - begin)*t + begin
		result.append(asarray(point).squeeze())
		
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
	Cset = []
	if close == True:
		if len( controls ) %3 != 0 or len(controls) < 3:
			print 'bad number of control points.'
			return None
			
		for i in range( len( controls )//3 -1 ):
			Cset.append(controls[i*3:i*3+4].tolist())
		
		last = controls[-3:].tolist() + [controls[0].tolist()]
		Cset.append(last)
		Cset = asarray( Cset )
	else:
		if len( controls ) %3 != 1 or len(controls) < 4:
			print 'bad number of control points.'
			return None
			
		for i in range( len( controls )//3 ):
			Cset.append(controls[i*3:i*3+4].tolist())
		
		Cset = asarray( Cset )
	
	return Cset