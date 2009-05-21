import scipy

from sfepy.base.base import *
from sfepy.solvers.solvers import LinearSolver

def try_imports(imports, fail_msg=None):
    for imp in imports:
        try:
            exec imp
            break
        except:
            pass
        else:
            if fail_msg is not None:
                raise ValueError(failmsg)
    return locals()

class ScipyDirect(LinearSolver):
    name = 'ls.scipy_direct'

    def process_conf(conf):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_1100 = {
            'name' : 'dls1100',
            'kind' : 'ls.scipy_direct',

            'method' : 'superlu',
            'presolve' : False,
            'warn' : True,
        }
        """
        get = conf.get_default_attr

        method = get('method', 'auto')
        presolve = get('presolve', False)
        warn = get('warn', True)

        common = LinearSolver.process_conf(conf)
        return Struct(**locals()) + common
    process_conf = staticmethod(process_conf)

    def __init__(self, conf, **kwargs):
        LinearSolver.__init__(self, conf, **kwargs)

        um = self.sls = None
        
        aux = try_imports(['import scipy.linsolve as sls',
                           'import scipy.splinalg.dsolve as sls',
                           'import scipy.sparse.linalg.dsolve as sls'],
                          'cannot import scipy sparse direct solvers!')
        self.sls = aux['sls']
        aux = try_imports(['import scipy.linsolve.umfpack as um',
                           'import scipy.splinalg.dsolve.umfpack as um',
                           'import scipy.sparse.linalg.dsolve.umfpack as um',
                           'import scikits.umfpack as um'])
        if um in aux:
            um = aux['um']

        if um is not None:
            is_umfpack = hasattr(um, 'UMFPACK_OK')
        else:
            is_umfpack = False

        method = self.conf.method
        if method == 'superlu':
            self.sls.use_solver(useUmfpack=False)
        elif method == 'umfpack':
            if not is_umfpack and self.conf.warn:
                output('umfpack not available, using superlu!')
        elif method != 'auto':
            raise ValueError('uknown solution method! (%s)' % method)

        if method != 'superlu' and is_umfpack:
            self.sls.use_solver(useUmfpack=True,
                                assumeSortedIndices=True)

        self.solve = None
        if self._presolve() and hasattr(self, 'mtx'):
            if self.mtx is not None:
                self.solve = self.sls.factorized(self.mtx)

    def __call__(self, rhs, conf=None, mtx=None, status=None):
        conf = get_default(conf, self.conf)
        mtx = get_default(mtx, self.mtx)
        status = get_default(status, self.status)

        if self.solve is not None:
            # Matrix is already prefactorized.
            return self.solve(rhs)
        else:
            return self.sls.spsolve(mtx, rhs)
                
    def _presolve(self):
        if hasattr(self, 'presolve'):
            return self.presolve
        else:
            try:
                return self.conf.presolve
            except:
                return False

class Umfpack(ScipyDirect):
    """This class stays for compatability with old input files. Use ScipyDirect
    isntead."""
    name = 'ls.umfpack'

    def __init__(self, conf, **kwargs):
        conf.method = 'umfpack'
        ScipyDirect.__init__(self, conf, **kwargs)

##
# c: 22.02.2008
class ScipyIterative( LinearSolver ):
    name = 'ls.scipy_iterative'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_110 = {
            'name' : 'ls110',
            'kind' : 'ls.scipy_iterative',

            'method' : 'cg',
            'i_max'   : 1000,
            'eps_a'   : 1e-12,
        }
        """
        get = conf.get_default_attr

        method = get( 'method', 'cg' )
        i_max = get( 'i_max', 100 )
        eps_a = get( 'eps_a', 1e-8 )

        common = LinearSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )
    
    ##
    # c: 22.02.2008, r: 23.06.2008
    def __init__( self, conf, **kwargs ):
        if scipy.version.version < '0.7.0.dev3861':
            import scipy.linalg as la
        else:
            if scipy.version.version < '0.7.0.dev3998':
                import scipy.splinalg.isolve as la
            else:
                import scipy.sparse.linalg.isolve as la

        LinearSolver.__init__( self, conf, **kwargs )

        try:
            solver = getattr( la, self.conf.method )
        except AttributeError:
            output( 'scipy solver %s does not exist!' % self.conf.method )
            output( 'using cg instead' )
            solver = la.cg
        self.solver = solver
        
    ##
    # c: 22.02.2008, r: 22.02.2008
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = get_default( conf, self.conf )
        mtx = get_default( mtx, self.mtx )
        status = get_default( status, self.status )

        sol, info = self.solver( mtx, rhs, tol = conf.eps_a,
                                 maxiter = conf.i_max )
        
        return sol

##
# c: 02.05.2008, r: 02.05.2008
class PyAMGSolver( LinearSolver ):
    name = 'ls.pyamg'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_102 = {
            'name' : 'ls102',
            'kind' : 'ls.pyamg',

            'method' : 'smoothed_aggregation_solver',
            'eps_a'   : 1e-12,
        }
        """
        get = conf.get_default_attr

        method = get( 'method', 'smoothed_aggregation_solver' )
        eps_a = get( 'eps_a', 1e-8 )

        common = LinearSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

    ##
    # c: 02.05.2008, r: 02.05.2008
    def __init__( self, conf, **kwargs ):
        try:
            import pyamg
        except ImportError:
            msg =  'cannot import pyamg!'
            raise ImportError( msg )

        LinearSolver.__init__( self, conf, mg = None, **kwargs )

        try:
            solver = getattr( pyamg, self.conf.method )
        except AttributeError:
            output( 'pyamg.%s does not exist!' % self.conf.method )
            output( 'using pyamg.smoothed_aggregation_solver instead' )
            solver = pyamg.smoothed_aggregation_solver
        self.solver = solver

        if hasattr( self, 'mtx' ):
            if self.mtx is not None:
                self.mg = self.solver( self.mtx )
        
    ##
    # c: 02.05.2008, r: 02.05.2008
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = get_default( conf, self.conf )
        mtx = get_default( mtx, self.mtx )
        status = get_default( status, self.status )

        if (self.mg is None) or (mtx is not self.mtx):
            self.mg = self.solver( mtx )
            self.mtx = mtx

        sol = self.mg.solve( rhs, tol = conf.eps_a )
        
        return sol

class PETScKrylovSolver( LinearSolver ):
    """PETSc Krylov subspace solver.

    The solver and preconditioner types are set upon the solver object
    creation. Tolerances can be overriden when called by passing a conf
    object.

    Convergence is reached when 'rnorm < max(eps_r * rnorm_0, eps_a)', where
    rnorm is the norm of rezidual, i.e. ||Ax - b||.
    """
    name = 'ls.petsc'

    def process_conf( conf ):
        """
        Missing items are set to default values.
        
        Example configuration, all items:
        
        solver_120 = {
            'name' : 'ls120',
            'kind' : 'ls.petsc',

            'method' : 'cg', # ksp_type
            'precond' : 'icc', # pc_type
            'eps_a' : 1e-12, # abstol
            'eps_r' : 1e-12, # rtol
            'i_max' : 1000, # maxits
        }
        """
        get = conf.get_default_attr

        method = get( 'method', 'cg' )
        precond = get( 'precond', 'icc' )
        eps_a = get( 'eps_a', 1e-8 )
        eps_r = get( 'eps_r', 1e-8 )
        i_max = get( 'i_max', 100 )

        common = LinearSolver.process_conf( conf )
        return Struct( **locals() ) + common
    process_conf = staticmethod( process_conf )

    def __init__( self, conf, **kwargs ):
        try:
            import petsc4py
            petsc4py.init([])
            from petsc4py import PETSc
        except ImportError:
            msg = 'cannot import petsc4py!'
            raise ImportError( msg )

        LinearSolver.__init__( self, conf, petsc = PETSc, pmtx = None, **kwargs )

        ksp = PETSc.KSP().create()

        ksp.setType( self.conf.method )
        ksp.getPC().setType( self.conf.precond )

        if hasattr( self, 'mtx' ):
            if self.mtx is not None:
                self.pmtx, self.sol, self.rhs = self.set_matrix( self.mtx )
                ksp.setOperators( self.pmtx ) # set the matrix
                ksp.setFromOptions()

        self.ksp = ksp


    def set_matrix( self, mtx ):
        pmtx = self.petsc.Mat().createAIJ( mtx.shape,
                                           csr = (mtx.indptr,
                                                  mtx.indices,
                                                  mtx.data) )
        sol, rhs = pmtx.getVecs()
        return pmtx, sol, rhs
                
    def __call__( self, rhs, conf = None, mtx = None, status = None ):
        conf = get_default( conf, self.conf )
        mtx = get_default( mtx, self.mtx )
        status = get_default( status, self.status )

        ksp = self.ksp
        if (self.pmtx is None) or (mtx is not self.mtx):
            self.pmtx, self.sol, self.rhs = self.set_matrix( mtx )
            self.ksp.setOperators( self.pmtx )
            self.ksp.setFromOptions() # PETSc.Options() not used yet...
            self.mtx = mtx

        ksp.setTolerances( atol = conf.eps_a, rtol = conf.eps_r,
                           max_it = conf.i_max )

        # Set PETSc rhs, solve, get solution from PETSc solution.
        self.rhs[...] = rhs
        ksp.solve( self.rhs, self.sol )
        sol = self.sol[...].copy()
        
        return sol
