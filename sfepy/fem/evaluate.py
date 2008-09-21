from sfepy.base.base import *
import extmods.fem as fem
from sfepy.terms import Term, DataCaches
from region import Region
from equations import Equation, build_args
from integrals import Integrals, quadratures

##
# 02.10.2007, c
class Evaluator( Struct ):
    pass

##
# 02.10.2007, c
class BasicEvaluator( Evaluator ):
    ##
    # 02.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        Evaluator.__init__( self, problem = problem, data = kwargs )
        if mtx is None:
            self.mtx = problem.mtx_a
        else:
            self.mtx = mtx

    ##
    # c: 11.04.2008, r: 11.04.2008
    def set_term_args( self, **kwargs ):
        self.data = kwargs

    ##
    # 02.12.2005, c
    # 25.07.2006
    # 02.10.2007
    def eval_residual( self, vec ):
        status = 0
        try:
            pb = self.problem
            vec_r = eval_residuals( vec, pb.equations, pb.conf.fe.chunk_size,
                                  **self.data )
        except StopIteration, exc:
            vec_r = None
            status = exc.args[0]
            print 'error %d in term "%s" of equation "%s"!'\
                  % (status, exc.args[1].name, exc.args[2].desc)

        return vec_r, status
            
    ##
    # 02.12.2005, c
    # 25.07.2006
    # 02.10.2007
    def eval_tangent_matrix( self, vec, mtx = None ):
        status = 0
        try:
            pb = self.problem
            if mtx is None:
                mtx = self.mtx
            mtx.data[:] = 0.0
            mtx = eval_tangent_matrices( vec, mtx,
                                       pb.equations, pb.conf.fe.chunk_size,
                                       **self.data )
        except StopIteration, exc:
            status = exc.args[0]
            print ('error %d in term "%s" of derivative of equation "%s"'
                   + ' with respect to variable "%s"!')\
                  % (status, exc.args[1].name, exc.args[2].desc, exc.args[3] )

        return mtx, status

    ##
    # 02.12.2005, c
    # 09.12.2005
    # 25.07.2006
    # 02.10.2007
    def update_vec( self, vec, delta ):
        self.problem.update_vec( vec, delta )

##
# 04.10.2007, c
class LCBCEvaluator( BasicEvaluator ):

    ##
    # 04.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        BasicEvaluator.__init__( self, problem, mtx, **kwargs )
        self.op_lcbc = problem.variables.get_lcbc_operator()

    ##
    # 04.10.2007, c
    def eval_residual( self, vec ):
        vec_r, status = BasicEvaluator.eval_residual( self, vec )
        vec_rr = self.op_lcbc.T * vec_r
        return vec_rr, status
            
    ##
    # 04.10.2007, c
    def eval_tangent_matrix( self, vec, mtx = None ):
        mtx, status = BasicEvaluator.eval_tangent_matrix( self, vec, mtx )
        mtx_r = self.op_lcbc.T * mtx * self.op_lcbc
        mtx_r = mtx_r.tocsr()
        mtx_r.sort_indices()
##         import pylab
##         from sfepy.base.plotutils import spy
##         spy( mtx_r )
##         pylab.show()
        print mtx_r.__repr__()
        return mtx_r, status

    ##
    # 04.10.2007, c
    def update_vec( self, vec, delta_r ):
        delta = self.op_lcbc * delta_r
        BasicEvaluator.update_vec( self, vec, delta )

def assemble_vector( vec, equation, variables, materials,
                    chunk_size = 1000, **kwargs ):
    get_a_dof_conn = variables.get_a_dof_conn

    for term in equation.terms:
##         print '>>>>>>', term.name, term.sign
        var_names = term.get_variable_names()
        args = build_args( term, variables, materials, **kwargs )
        vn = term.get_virtual_name( variables = variables )
        dc_type = term.get_dof_conn_type()
##         print args
##         print vn
##         print dc_type
        for ig in term.iter_groups():
            dc = get_a_dof_conn( vn, True, dc_type, ig )
##             print vn, dc.shape
#            pause()
            for vec_in_els, iels, status in term( chunk_size = chunk_size,
                                                **args ):
                if status != 0:
                    raise StopIteration( status, term, equation )

                check_finitness = False
                if check_finitness and (not nm.isfinite( vec_in_els ).all()):
                    print term.name, term.sign, ig
                    debug()

                if vec.dtype == nm.float64:
                    fem.assemble_vector( vec, vec_in_els, iels, term.sign, dc )

                else:
                    assert vec.dtype == nm.complex128
                    sign = nm.array( term.sign, dtype = nm.complex128 )
                    fem.assemble_vector_complex( vec.real, vec.imag,
                                                 vec_in_els.real,
                                                 vec_in_els.imag,
                                                 iels,
                                                 float( sign.real ),
                                                 float( sign.imag ), dc )

def assemble_matrix( mtx, equation, variables, materials,
                    chunk_size = 1000, group_can_fail = True, **kwargs ):
    if not sp.isspmatrix_csr( mtx ):
        raise TypeError, 'must be CSR matrix!'
    tmd = (mtx.data, mtx.indptr, mtx.indices)

    get_a_dof_conn = variables.get_a_dof_conn

    for term in equation.terms:
#        print '>>>>>>', term.name, term.sign
        var_names = term.get_variable_names()
        args = build_args( term, variables, materials, **kwargs )
        vn = term.get_virtual_name( variables = variables )
        sns = term.get_state_names( variables = variables )
        dc_type = term.get_dof_conn_type()

        for ig in term.iter_groups():
            rdc = get_a_dof_conn( vn, True, dc_type, ig )
#            print vn, rdc.shape
            for sn in sns:
                cdc = get_a_dof_conn( sn, False, dc_type, ig )
#                print sn, cdc.shape
#                pause()
                for mtx_in_els, iels, status in term( diff_var = sn,
                                                    chunk_size = chunk_size,
                                                    **args ):
                    if status != 0:
                        raise StopIteration( status, term, equation,
                                             var_name_col )

                    if mtx.dtype == nm.float64:
                        fem.assemble_matrix( tmd[0], tmd[1], tmd[2], mtx_in_els,
                                             iels, term.sign, rdc, cdc )

                    else:
                        assert mtx.dtype == nm.complex128
                        sign = nm.array( term.sign, dtype = nm.complex128 )
                        fem.assemble_matrix_complex( tmd[0].real, tmd[0].imag,
                                                     tmd[1], tmd[2],
                                                     mtx_in_els.real,
                                                     mtx_in_els.imag,
                                                     iels,
                                                     float( sign.real ),
                                                     float( sign.imag ),
                                                     rdc, cdc )

##
# 01.10.2007, c
def eval_term_op( state, term_desc, problem, **kwargs ):
    """Convenience wrapper of eval_term() in a context of ProblemDefinition
    instance."""
    return eval_term( state, term_desc, problem.conf,
                     problem.domain, problem.variables, problem.materials,
                     chunk_size = problem.domain.shape.n_el, **kwargs )

##
# c: 03.01.2006, r: 05.03.2008
def eval_term( state, term_desc, conf, domain, variables, materials,
              funmod = None, chunk_size = 1000, term_prefixes = None,
              caches = None, ret_caches = False,
              override = True, new_geometries = True,
              dw_mode = 'vector', tangent_matrix = None,
              **kwargs ):
    """Evaluate a term. May not succeed!"""
    if term_prefixes is None:
        term_prefixes = {}
    if caches is None:
        caches = DataCaches()

    equation = Equation.from_desc( 'tmp', term_desc, term_prefixes )
    equation.setup_terms( domain.regions, variables, materials, caches, kwargs )
    for cache in caches.itervalues():
        cache.set_mode( override = override )

    if new_geometries:
        i_names = equation.get_term_integral_names()
        integrals = Integrals.from_conf( conf.integrals, i_names )
        integrals.set_quadratures( quadratures )
        
        geometries = {}
        equation.describe_geometry( geometries, variables, integrals )

    variables.data_from_state( state )
    # itype according to the first term in term_desc!
    term = equation.terms[0]
    if term.itype == 'dw':

        variables.setup_dof_conns()
        if dw_mode == 'vector':
            residual = variables.create_stripped_state_vector()
            assemble_vector( residual, equation, variables, materials,
                            chunk_size, group_can_fail = False, **kwargs )
            if variables.has_lcbc:
                op_lcbc = variables.op_lcbc
                residual = op_lcbc.T * residual
            ret_val = residual

        elif dw_mode == 'matrix':
            if tangent_matrix is None:
                tangent_matrix = variables.create_matrix_graph()

            tangent_matrix.data[:] = 0.0
            assemble_matrix( tangent_matrix, equation, variables, materials,
                            chunk_size, group_can_fail = False, **kwargs )
            if variables.has_lcbc:
                op_lcbc = variables.op_lcbc
                tangent_matrix = op_lcbc.T * tangent_matrix * op_lcbc
                tangent_matrix = tangent_matrix.tocsr()
                tangent_matrix.sort_indices()
            ret_val = tangent_matrix

        else:
            print dw_mode
            raise ValueError

    elif term.itype == 'd':
        val = 0.0

        for term in equation.terms:
            args = build_args( term, variables, materials, **kwargs )
            for ig in term.iter_groups():
                for aux, iels, status in term( chunk_size = chunk_size,
                                               **args ):
                    val += term.sign * aux
            ret_val = val

    elif term.itype == 'di':
        val = None

        for term in equation.terms:
            args = build_args( term, variables, materials, **kwargs )
            for ig in term.iter_groups():
                for aux, iels, status in term( chunk_size = chunk_size,
                                               **args ):
                    if val is None:
                        val = term.sign * aux
                    else:
                        val += term.sign * aux
            ret_val = val

    elif (term.itype == 'de') or (term.itype == 'dq'):
        val = None

        for term in equation.terms:
            args = build_args( term, variables, materials, **kwargs )
            for ig in term.iter_groups():
                for aux, iels, status in term( chunk_size = chunk_size,
                                               **args ):
                    if val is None:
                        val = term.sign * aux
                    else:
                        val = nm.concatenate( (val, term.sign * aux), axis = 0 )
        ret_val = val

    else:
        raise NotImplementedError, 'unknown term int. type: %s' % term.itype

    if ret_caches:
        return ret_val, caches
    else:
        return ret_val

##
# 21.11.2005, c
# 27.11.2005
# 02.12.2005
# 09.12.2005
# 14.12.2005
# 21.03.2006
# 24.07.2006
# 25.07.2006
# 22.08.2006
# 11.10.2006
# 16.02.2007
# 03.09.2007
def eval_residuals( state, equations, chunk_size = 1000,
                   **kwargs ):

    variables = equations.variables
    materials = equations.materials

    variables.data_from_state( state )

    residual = variables.create_stripped_state_vector()
    equations.invalidate_term_caches()

    for equation in equations:
        assemble_vector( residual, equation, variables, materials,
                        chunk_size = chunk_size, **kwargs )

    return residual

##
# 22.11.2005, c
# 26.11.2005
# 02.12.2005
# 09.12.2005
# 14.12.2005
# 21.03.2006
# 24.07.2006
# 25.07.2006
# 22.08.2006
# 11.10.2006
# 16.02.2007
# 03.09.2007
def eval_tangent_matrices( state, tangent_matrix, equations, chunk_size = 1000,
                         **kwargs ):

    variables = equations.variables
    materials = equations.materials

    variables.data_from_state( state )

    for equation in equations:
        assemble_matrix( tangent_matrix, equation, variables, materials,
                        chunk_size = chunk_size, **kwargs )

    return tangent_matrix
