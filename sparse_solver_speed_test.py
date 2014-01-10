import timeit, os, sys

setup = '''
print "===> Starting setup"
import sys
excepthook = sys.excepthook
import generate_chain_system_config
generate_chain_system_config.kSystemSolvePackage = '%s'
generate_chain_system_config.kBuildDense = '%s'
import chain_computer, bezier_utility, generate_chain_system
bezier_utility.kG1andAconstraints = %s

sys.excepthook = excepthook

try:
    paths_info, skeleton_handle_vertices, constraint = chain_computer.get_test_%s()
    engine = chain_computer.Engine()
    boundary_path = max(paths_info, key=lambda e : e[u'bbox_area']) 
    boundary_index = paths_info.index( boundary_path )
    engine.set_control_positions( paths_info, boundary_index )
    engine.set_handle_positions( skeleton_handle_vertices )
    engine.precompute_configuration()
except Exception as e:
    print 'Setup died:'
    print e
    
    import time
    class Engine( object ):
        def solve( self ):
            return None
    engine = Engine()
print "===> Finished setup."
'''

statement = '''
all_paths = engine.solve()
'''

#solves = [ 'numpy-inv', 'numpy-solve', 'scipy', 'cvxopt' ]
solves = [ 'numpy-solve', 'scipy', 'cvxopt' ]
denses = [ 'numpy', 'scipy', 'cvxopt' ]
#whichs = ( 'simple_closed', 'pebble', 'alligator' )
#whichs = ( 'simple_closed', )
whichs = sys.argv[1:]

from numpy import array

N = 100
R = 3
for which in whichs:
    for g1 in ( True, False ):
        for solve in solves:
            for dense in denses:
                
                ## Fork before running this, because timeit does not actually re-import modules,
                ## which we require in order to actually run the different solvers.
                pid = os.fork()
                if pid != 0:
                    os.wait()
                    continue
                
                print '=======> Running', solve, 'build-dense', dense, which, 'G1-and-A', g1, 'with repeat', R, 'and number', N
                durations = timeit.repeat( statement, setup = setup % ( solve, dense, g1, which ), repeat = 3, number = N )
                seconds_per_calls = list(sorted( (array(durations)/N).tolist() ) )
                print '=======>', round( min( seconds_per_calls ), 3 ), 'Duration of', solve, 'build-dense', dense, which, 'G1-and-A', g1, 'was', durations, 'aka', seconds_per_calls, 'seconds per call'
                sys.exit()
