import numpy

kDefaultSystemSolvePackage = 'scipy'
kDefaultBuildType = 'numpy'

def get_system_and_factor_funcs_for_system_size( system_size, G1orA ):
    '''
    Returns the best building and solving system combination for
    the given system size and whether the system will be used
    in a setting with G1orA constraints.
    If G1orA is False, then the system will be factored once and then back-substituted
    repeatedly.
    Otherwise, it will be sovled anew every time.
    '''
    
    ## At least up to system size 500, solver numpy-inv/scipy and matrix numpy/cvxopt are tied.
    ## Scipy is slightly slower for the solve-every-time case with matrices less than 200.
    ## For the G1 case, it looks like numpy matrix type has a slight advantage over cvxopt,
    ## perhaps because the cvxopt type grows indefinitely since it's coo.
    ## For systems of size 1000x1000 or 2000x2000, scipy solver is faster,
    ## but numpy matrix still wins.
    
    if G1orA and system_size < 200:
        return get_system_and_factor_funcs( 'numpy-inv', 'numpy' )
    else:
        ## This is the best up to even our largest examples (despite being a dense matrix)
        ## and down to our smallest systems (provided that the same matrix factorization
        ## is used for many right-hand-sides).
        return get_system_and_factor_funcs( 'scipy', 'numpy' )

def get_system_and_factor_funcs( solver_type = None, build_system_type = None ):
    '''
    Returns a triplet of functions:
        The first acts like numpy.zeros, returning a matrix-like object,
        The second converts the matrix-like object to an object suitable for passing to the third,
        The third is a factorization function:
            Given a system matrix 'system',
            returns a function that can be used to solve for 'x' in
                system * x = b
            as follows:
                x = compute_symbolic_factorization( system )( system )( b )
    '''
    
    ## Defaults
    if solver_type is None: solver_type = kDefaultSystemSolvePackage
    if build_system_type is None: build_system_type = kDefaultBuildType
    
    if 'scipy' == solver_type or 'scipy' == build_system_type: import scipy.sparse.linalg
    if 'cvxopt' == solver_type or 'cvxopt' == build_system_type: import cvxopt
    
    if 'numpy' == build_system_type:
        zeros_system_build_t = numpy.zeros
    elif 'scipy' == build_system_type:
        zeros_system_build_t = scipy.sparse.lil_matrix
    elif 'cvxopt' == build_system_type:
        zeros_system_build_t = lambda shape: cvxopt.spmatrix( [], [], [], shape )
    
    def scipy2cvx( m ):
        m = m.tocoo()
        cvx = cvxopt.spmatrix( numpy.asarray( m.data ), numpy.asarray( m.row, dtype = int ), numpy.asarray( m.col, dtype = int ) )
        return cvx
    def cvx2scipy( m ):
        return scipy.sparse.coo_matrix( ( numpy.asarray( m.V ).squeeze(), ( numpy.asarray( m.I ).squeeze(), numpy.asarray( m.J ).squeeze() ) ), m.size )
    
    if 'numpy-inv' == solver_type:
        if 'numpy' == build_system_type:
            to_system_solve_t = lambda x: x
        elif 'scipy' == build_system_type:
            to_system_solve_t = lambda x: asarray( x.todense() )
        elif 'cvxopt' == build_system_type:
            import cvxopt
            to_system_solve_t = lambda x: numpy.asarray( cvxopt.matrix( x ) )
        else:
            raise RuntimeError( "Unknown build dense type: " + str( build_system_type ) )
        
        def compute_symbolic_factorization( system ):
            '''
            Given a scipy.sparse system matrix 'system',
            return a function that can be used to solve for 'x' in
                system * x = b
            as follows:
                x = compute_symbolic_factorization( system )( system )( b )
            '''
            def compute_numeric_factorization( system ):
                try:
                    inverse = numpy.linalg.inv( system )
                except numpy.linalg.linalg.LinAlgError as e:
                    print e
                    inverse = numpy.eye( system.shape[0] )
                def solve( rhs ):
                    return numpy.dot( inverse, rhs )
                return solve
            return compute_numeric_factorization
    
    elif 'numpy-solve' == solver_type:
        if 'numpy' == build_system_type:
            to_system_solve_t = lambda x: x
        elif 'scipy' == build_system_type:
            to_system_solve_t = lambda x: asarray( x.todense() )
        elif 'cvxopt' == build_system_type:
            import cvxopt
            to_system_solve_t = lambda x: numpy.asarray( cvxopt.matrix( x ) )
        else:
            raise RuntimeError( "Unknown build dense type: " + str( build_system_type ) )
        
        def compute_symbolic_factorization( system ):
            '''
            Given a scipy.sparse system matrix 'system',
            return a function that can be used to solve for 'x' in
                system * x = b
            as follows:
                x = compute_symbolic_factorization( system )( system )( b )
            '''
            def compute_numeric_factorization( system ):
                def solve( rhs ):
                    return numpy.linalg.solve( system, rhs )
                return solve
            return compute_numeric_factorization
    
    elif 'scipy' == solver_type:
        if 'numpy' == build_system_type:
            to_system_solve_t = scipy.sparse.csc_matrix
        elif 'scipy' == build_system_type:
            to_system_solve_t = scipy.sparse.csc_matrix
        elif 'cvxopt' == build_system_type:
            to_system_solve_t = lambda x: cvx2scipy( x ).tocsc()
        else:
            raise RuntimeError( "Unknown build dense type: " + str( build_system_type ) )
        
        def compute_symbolic_factorization( system ):
            '''
            Given a scipy.sparse system matrix 'system',
            return a function that can be used to solve for 'x' in
                system * x = b
            as follows:
                x = compute_symbolic_factorization( system )( system )( b )
            '''
            def compute_numeric_factorization( system ):
                factorized = scipy.sparse.linalg.factorized( system )
                def solve( rhs ):
                    ## spsolve() is slower than the factorized above.
                    # return scipy.sparse.linalg.spsolve( system, rhs )
                    return factorized( rhs )
                return solve
            return compute_numeric_factorization
    
    elif 'cvxopt' == solver_type:
        ## UPDATE: cholmod dies with our system for some reason.
        #import cvxopt.cholmod as cvxopt_solver
        import cvxopt.umfpack as cvxopt_solver
        
        if 'numpy' == build_system_type:
            to_system_solve_t = lambda x: cvxopt.sparse( cvxopt.matrix( x ) )
        elif 'scipy' == build_system_type:
            to_system_solve_t = scipy2cvx
        elif 'cvxopt' == build_system_type:
            to_system_solve_t = lambda x: x
        else:
            raise RuntimeError( "Unknown build dense type: " + str( build_system_type ) )
        
        def compute_symbolic_factorization( system ):
            '''
            Given a scipy.sparse system matrix 'system',
            return a function that can be used to solve for 'x' in
                system * x = b
            as follows:
                x = compute_symbolic_factorization( system )( system )( b )
            '''
            symbolic_factorization = cvxopt_solver.symbolic( system )
            def compute_numeric_factorization( system ):
                # assert abs(asarray( cvxopt.matrix( (system - system.T) ) )).max() < 1e-8
                full_factorization = cvxopt_solver.numeric( system, symbolic_factorization )
                def solve( rhs ):
                    x = cvxopt.matrix( rhs )
                    ## If cvxopt_solver is cvxopt.umfpack
                    cvxopt_solver.solve( system, full_factorization, x )
                    ## If cvxopt_solver is cvxopt.cholmod
                    ## UPDATE: Neither of these work.
                    # cvxopt_solver.solve( full_factorization, x, sys = 1 )
                    # cvxopt_solver.spsolve( system, x, sys = 0 )
                    return numpy.asarray( x ).reshape( rhs.shape )
                return solve
            return compute_numeric_factorization
    
    else:
        raise RuntimeError( "Unknown system solve package " + str( solver_type ) )
    
    return zeros_system_build_t, to_system_solve_t, compute_symbolic_factorization
