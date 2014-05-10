import time

starts = []

def tic( msg = None ):
#     if msg is not None:
#         print ( '+' * (len( starts )+1) ), msg
    
    starts.append( ( time.clock(), msg ) )

def toc():
    end = time.clock()
    duration = end - starts[-1][0]
    
#     if starts[-1][1] is None:
#         print ( '=' * len( starts ) ), 'tictoc():', duration
#     else:
#         print ( '-' * len( starts ) ), starts[-1][1], duration
    
    del starts[-1]

from contextlib import contextmanager
@contextmanager
def tictoc( msg = None ):
    tic( msg )
    yield
    toc()

def tictoc_dec( func ):
	
	def wrapped( *args, **kwargs ):
		with tictoc( func.func_name ):
			return func( *args, **kwargs )
	
	return wrapped
