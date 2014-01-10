import timeit, os, sys

setup = '''
import generate_chain_system_config
generate_chain_system_config.kSystemSolvePackage = '%s'
generate_chain_system_config.kBuildDense = %s
import chain_computer, bezier_utility, generate_chain_system
bezier_utility.kG1andAconstraints = %s

paths_info, skeleton_handle_vertices, constraint = chain_computer.get_test_%s()
engine = chain_computer.Engine()
boundary_path = max(paths_info, key=lambda e : e[u'bbox_area']) 
boundary_index = paths_info.index( boundary_path )
engine.set_control_positions( paths_info, boundary_index )
engine.set_handle_positions( skeleton_handle_vertices )
engine.precompute_configuration()
'''

statement = '''
all_paths = engine.solve()
'''

solves = [ 'numpy-inv', 'numpy-solve', 'scipy', 'cvxopt' ]
whichs = ( 'simple_closed', 'pebble', 'alligator' )

from numpy import array

N = 100
R = 3
for which in whichs:
    for solve in solves:
        for dense in ( True, False ):
            for g1 in ( True, False ):
                
                ## Fork before running this, because timeit does not actually re-import modules,
                ## which we require in order to actually run the different solvers.
                pid = os.fork()
                if pid != 0:
                    os.wait()
                    continue
                
                print '=======> Running', solve, dense, which, 'with repeat', R, 'and number', N
                durations = timeit.repeat( statement, setup = setup % ( solve, dense, g1, which ), repeat = 3, number = N )
                seconds_per_calls = list(sorted( (array(durations)/N).tolist() ) )
                print '=======> Duration of', solve, dense, which, 'was', durations, 'aka', seconds_per_calls, 'seconds per call'
                sys.exit()
