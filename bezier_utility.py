from math import *
from numpy import *
from myarray import *

kEps = 1e-7
kG1andAconstraints = True

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
	all_pts = []
	all_ts = []
	all_dts = []
	
	for P in Cset:
		samples, ts, dts = sample_cubic_bezier_curve( P, num_samples )
		all_pts.append( samples )
		all_ts.append( ts )
		all_dts.append( dts )

	return all_pts, all_ts, all_dts

def sample_cubic_bezier_curve( P, num_samples = 100 ):
	'''
	a 4-by-k numpy.array P containing the positions of the control points as the rows,
	return two lists: sample points of the bezier curve denoted in P, and corresponding t values
	'''
	if num_samples is None:
		num_samples = max(int(length_of_cubic_bezier_curve(P) / 1), 2)
	
	P = asarray( P )
	result = []
	ts = []
	tbar = ones( (4,1) )
	for t in linspace( 0, 1, num_samples ):
		ts.append( t )
		
		tbar[0] = t**3
		tbar[1] = t**2
		tbar[2] = t
		
		point = dot( P.T, dot( M.T, tbar ) )
		result.append( asarray(point).squeeze().tolist() )
	
	dts = ones( num_samples-1 ) * (1./(num_samples-1) )
		
	return asarray( result ), asarray( ts ), asarray( dts )
	
	
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
	'''
	make every control points a group from a chain of control points. Discard the remainder control points if there is any extra.
	'''
	controls = asarray(controls)	
	remain_num = len( controls ) % 3
	
	if remain_num == 2:
		print 'wrong num for controls: ', len(controls), ': ', controls, 'closed:', close  
	elif remain_num == 1:
		if close and not array_equal( controls[0], controls[-1] ):
			print 'wrong num for closed controls: ', len(controls), ': ', controls
	else:
		if not close:
			debugger()
			print 'wrong num for open controls: ', len(controls), ': ', controls
	
	num_segment = len( controls )//3
	if close:
		controls = concatenate( (controls[:num_segment*3], controls[0].reshape(1,2) ), axis=0 )
	elif not close and remain_num == 0:
		num_segment -= 1

	controls = controls[ :3*num_segment+1 ]			
				
	Cset = []
	for i in range( num_segment ):
		Cset.append(controls[i*3:i*3+4].tolist())
	
	Cset = asarray( Cset )
	
	return Cset
	
def make_constraints_from_control_points( control_group, close=True ):
	'''
	Make default constraints based on the following assumptions:
	- if upon loading they appear to be G1, they should stay G1
	- C1 is G1
	- 90 degrees defaults to fixed angle
	- pinning should never be done; user can add later or handles placed near control points have that effect	
	'''
	control_group = asarray( control_group )
	num = len( control_group )
	constraints = [ ['C0', False] for i in range( num ) ]
	
	for i in range( num ):
		dir1 = dir_allow_zero( control_group[i,-1] - control_group[i,-2] )
		dir2 = dir_allow_zero( control_group[(i+1)%num,1] - control_group[(i+1)%num,0] )
		
		if allclose( dir1, dir2, atol=1e-03 ) and mag(dir1) != 0 and mag(dir2) != 0:
			if kG1andAconstraints:
				## G1
				constraints[ (i+1)%num ][0] = 'G1' 
			else:
				## C1
				constraints[ (i+1)%num ][0] = 'C1'
		elif allclose( dot( dir1, dir2 ), 0, atol=1e-03 ) and mag(dir1) != 0 and mag(dir2) != 0:
			if kG1andAconstraints:
				## fixed angle
				constraints[ (i+1)%num ][0] = 'A'	 
			else:
				## C0
				constraints[ (i+1)%num ][0] = 'C0'

		
	if not close:
		constraints[0][0] = 'None'
		constraints.append( ['None', False] )
		
	return constraints	
	
	
def split_cublic_beizer_curve( controls, partition ):
	'''
	controls are the four control points of an cubic bezier curve
	partition is an array has each splitted curve's portion
	e.g. [0.5, 0.3, 0.2] is an partition
	'''
	controls = asarray( controls )
	partition = asarray( partition )

	assert len(controls.shape) == 2
	assert controls.shape[0] == 4
	assert sum( partition ) == 1
	assert partition.any() > 0
	
	# num = len( partition )
	
	ref = 1.0
	for i in range( len(partition) ):
		partition[i] = partition[i] / ref
		ref = ( 1 - partition[i] ) * ref
	
	result = []
	for i, k in enumerate( partition[:-1] ):
		
		r1 = controls[:-1]*(1.-k) + controls[1:]*k
		r2 = r1[:-1]*(1.-k) + r1[1:]*k
		r3 = r2[:-1]*(1.-k) + r2[1:]*k
		
		result.append([controls[0].tolist(), r1[0].tolist(), r2[0].tolist(), r3[0].tolist()])

		controls = array( [r3[-1], r2[-1], r1[-1], controls[-1]] )
	
	result.append(controls)
	
	return asarray(result)
	