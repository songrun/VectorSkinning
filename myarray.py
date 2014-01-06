from numpy import *

# a small set of helper functions to
# call common array creation functions
# these are useful to ensure that
# all arrays are created as double-precision
# floats, no matter what data are provided
# as argument. For example array([1,3,4]) normally returns
# an array with data of type int, but arrayf([1,3,4])
# always creates an array of floats 

kFloatType = float64

def arrayf( arg ):
	return array( arg, kFloatType )
def asarrayf( arg ):
	return asarray( arg, kFloatType )
def zerosf( arg ):
	return zeros( arg, kFloatType )
def onesf( arg ):
	return ones( arg, kFloatType )
def identityf( arg ):
	return identity( arg, kFloatType )
def emptyf( arg ):
	return empty( arg, kFloatType )

def mag2( vec ):
	vec = asarray(vec)
	return dot(vec,vec)
def mag( vec ):
	return sqrt(mag2(vec))
def dir( vec ):
	assert mag(vec) != 0
	vec = asarray(vec)
	return vec * 1./mag(vec)
def dir_allow_zero( vec ):
	vec = asarray(vec)
	if mag(vec) == 0:
		return vec
	else:
		return vec * 1./mag(vec)
def activate( test ):
	if test >= 0:
		return 1.
	else:
		return 0.
def activates( test ):
	result = []
	for each in test:
		result.append( activate( each ) )
		
	result = asarray( result )		
	return result		