from math import *
from numpy import *

def ease_cos( t, a, b ):
	return -(b-a)/2. * (math.cos(pi*t) - 1.) + a

def ease_dict( t, a, b ):
    # assert set( a.keys() == set( b.keys() )
    r = {}
    for key in a.keys():
        r[key] = ease_cos( t, a[key], b[key] )
    return r

def make_anim():
    Ma = {u'a': 1, u'c': 0, u'b': 0, u'e': 0, u'd': 1, u'f': 0}
    #Mb = {u'a': 19.105973797508184, u'c': 0, u'b': 0, u'e': -4218.687694819407, u'd': 19.105973797508184, u'f': -2842.634986208785}
    Mb = {"a":26.3113393432659,"b":0,"c":0,"d":26.3113393432659,"e":-5897.542066980955,"f":-3973.880276892746}
    
    fps = 60
    wait = 5
    transition = 5
    
    anim = []
    ## wait at the beginning.
    for i in xrange( wait*fps ): anim.append( Ma )
    
    ## go up
    for t in linspace( 0, 1, transition*fps ):
        anim.append( ease_dict( t, Ma, Mb ) )
    
    ## wait
    for i in xrange( wait*fps ): anim.append( Mb )
    
    ## go down
    for t in linspace( 1, 0, transition*fps ):
        anim.append( ease_dict( t, Ma, Mb ) )
    
    return anim
