from sfepy.base.base import *
import extmods.fem as fem
from sfepy.terms import Term, DataCaches
from region import Region
from equations import Equation, Equations
from integrals import Integrals
from variables import Variables
from fields import setup_dof_conns, setup_extra_data

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

    def eval_residual( self, vec, is_full = False ):
        if not is_full:
            vec = self.make_full_vec( vec )
        try:
            pb = self.problem
            vec_r = eval_residuals(vec, pb.equations,
                                   pb.conf.fe.chunk_size, **self.data)

        except StopIteration, exc:
            status = exc.args[0]
            output( 'error %d in term "%s" of equation "%s"!'\
                    % (status, exc.args[1].name, exc.args[2].desc) )
            raise ValueError

        return vec_r
            
    def eval_tangent_matrix( self, vec, mtx = None, is_full = False ):
        if isinstance( vec, str ) and vec == 'linear':
            return get_default( mtx, self.mtx )
        
        if not is_full:
            vec = self.make_full_vec( vec )
        try:
            pb = self.problem
            if mtx is None:
                mtx = self.mtx
            mtx.data[:] = 0.0
            mtx = eval_tangent_matrices(vec, mtx, pb.equations,
                                        pb.conf.fe.chunk_size, **self.data)

        except StopIteration, exc:
            status = exc.args[0]
            output( ('error %d in term "%s" of derivative of equation "%s"'
                     + ' with respect to variable "%s"!')\
                    % (status,
                       exc.args[1].name, exc.args[2].desc, exc.args[3] ) )
            raise ValueError

        return mtx

    ##
    # 02.12.2005, c
    # 09.12.2005
    # 25.07.2006
    # 02.10.2007
    def update_vec( self, vec, delta ):
        self.problem.update_vec( vec, delta )

    def strip_state_vector( self, vec ):
        return self.problem.equations.strip_state_vector(vec,
                                                         follow_epbc=False)

    def make_full_vec( self, vec ):
        return self.problem.equations.make_full_vec(vec)

##
# 04.10.2007, c
class LCBCEvaluator( BasicEvaluator ):

    ##
    # 04.10.2007, c
    def __init__( self, problem, mtx = None, **kwargs ):
        BasicEvaluator.__init__( self, problem, mtx, **kwargs )
        self.op_lcbc = problem.equations.get_lcbc_operator()

    ##
    # 04.10.2007, c
    def eval_residual( self, vec, is_full = False ):
        if not is_full:
            vec = self.make_full_vec( vec )
        vec_r = BasicEvaluator.eval_residual( self, vec, is_full = True )
        vec_rr = self.op_lcbc.T * vec_r
        return vec_rr
            
    ##
    # 04.10.2007, c
    def eval_tangent_matrix( self, vec, mtx = None, is_full = False ):
        if isinstance( vec, str ) and vec == 'linear':
            return get_default( mtx, self.mtx )

        if not is_full:
            vec = self.make_full_vec( vec )
        mtx = BasicEvaluator.eval_tangent_matrix( self, vec, mtx = mtx,
                                                  is_full = True )
        mtx_r = self.op_lcbc.T * mtx * self.op_lcbc
        mtx_r = mtx_r.tocsr()
        mtx_r.sort_indices()
##         import pylab
##         from sfepy.base.plotutils import spy
##         spy( mtx_r )
##         pylab.show()
##         print mtx_r.__repr__()
        return mtx_r

    ##
    # 04.10.2007, c
    def update_vec( self, vec, delta_r ):
        delta = self.op_lcbc * delta_r
        BasicEvaluator.update_vec( self, vec, delta )

    def strip_state_vector( self, vec ):
        vec = BasicEvaluator.strip_state_vector( self, vec )
        return self.op_lcbc.T * vec
    
def assemble_vector(vec, equation, chunk_size=1000, **kwargs):

    for term in equation.terms:
        ## print '>>>>>>', term.name, term.sign
        vvar = term.get_virtual_variable()
        dc_type = term.get_dof_conn_type()
        ## print vvar
        ## print dc_type

        for ig in term.iter_groups():
            dc = vvar.get_dof_conn(dc_type, ig, active=True)
            ## print vvar.name, dc.shape

            for vec_in_els, iels, status in term( chunk_size = chunk_size,
                                                  **kwargs ):
                if status != 0:
                    raise StopIteration( status, term, equation )

                check_finitness = False
                if check_finitness and (not nm.isfinite( vec_in_els ).all()):
                    print term.name, term.sign, ig
                    debug()

                assert_( vec_in_els.shape[2] == dc.shape[1] )

                if vec.dtype == nm.float64:
                    fem.assemble_vector( vec, vec_in_els, iels, term.sign, dc )

                else:
                    assert_( vec.dtype == nm.complex128 )
                    sign = nm.array( term.sign, dtype = nm.complex128 )
                    fem.assemble_vector_complex( vec.real, vec.imag,
                                                 vec_in_els.real,
                                                 vec_in_els.imag,
                                                 iels,
                                                 float( sign.real ),
                                                 float( sign.imag ), dc )

def assemble_matrix(mtx, equation, chunk_size=1000,
                    group_can_fail=True, **kwargs):
    """Assemble tangent matrix. Supports backward time difference of state
    variables."""
    if not sp.isspmatrix_csr( mtx ):
        raise TypeError, 'must be CSR matrix!'
    tmd = (mtx.data, mtx.indptr, mtx.indices)

    for term in equation.terms:
        ## print '>>>>>>', term.name, term.sign
        vvar = term.get_virtual_variable()
        svars = term.get_state_variables(unknown_only=True)
        dc_type = term.get_dof_conn_type()

        for ig in term.iter_groups():
            rdc = vvar.get_dof_conn(dc_type, ig, active=True)
            ## print vvar.name, rdc.shape

            for svar in svars:
                is_trace = term.arg_traces[svar.name]
                ## print dc_type, ig, is_trace
                cdc = svar.get_dof_conn(dc_type, ig, active=True,
                                        is_trace=is_trace)
                ## print svar.name, cdc.shape
                ## pause()
                for mtx_in_els, iels, status in term( diff_var = svar.name,
                                                      chunk_size = chunk_size,
                                                      **kwargs ):
                    if status != 0:
                        raise StopIteration( status, term, equation,
                                             var_name_col )

                    assert_( mtx_in_els.shape[2:] == (rdc.shape[1],
                                                      cdc.shape[1]) )

                    sign = term.sign
                    if term.arg_derivatives[svar.name]:
                        sign *= 1.0 / term.dt

                    if mtx.dtype == nm.float64:
                        fem.assemble_matrix( tmd[0], tmd[1], tmd[2], mtx_in_els,
                                             iels, sign, rdc, cdc )

                    else:
                        assert_( mtx.dtype == nm.complex128 )
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
                      problem.domain, problem.fields, problem.materials,
                      problem.get_timestepper(), problem.ebcs,
                      problem.epbcs, problem.lcbcs, problem.functions,
                      chunk_size = problem.domain.shape.n_el, **kwargs )

##
# c: 03.01.2006, r: 05.03.2008
def eval_term( state, term_desc, conf, domain, fields, materials, ts,
               ebcs, epbcs, lcbcs, functions,
               funmod = None, chunk_size = 1000, term_prefixes = None,
               caches = None, ret_caches = False,
               override = True, new_geometries = True,
               dw_mode = 'vector', tangent_matrix = None,
               copy_materials = True, update_materials = True, **kwargs ):
    """Evaluate a term. May not succeed!"""
    if term_prefixes is None:
        term_prefixes = {}
    if caches is None:
        caches = DataCaches()

    if update_materials:
        if copy_materials:
            materials = materials.semideep_copy()

    variables = Variables.from_conf(conf.variables, fields)
    equation = Equation.from_desc('tmp', term_desc, variables, domain.regions,
                                  materials, caches=caches,
                                  user=kwargs, term_prefixes=term_prefixes)

    for cache in caches.itervalues():
        cache.set_mode( override = override )

    if 'call_mode' in kwargs:
        itype = kwargs['call_mode'].split( '_' )[0]
    else:
        # itype according to the first term in term_desc!
        itype = equation.terms[0].itype

    if new_geometries:
        equations = Equations([equation], caches=caches)
        equations.collect_conn_info()

        if itype == 'dw':
            setup_dof_conns(equations.conn_info, single_term=True)

        else:
            setup_extra_data(equations.conn_info)

        integrals = Integrals.from_conf(conf.integrals)
        equations.describe_geometry(integrals, verbose=False)
        variables = equations.variables

    else:
        variables.setup_dof_info()

    variables.set_data(state, ignore_unknown=True)

    if update_materials:
        materials.time_update(ts, domain,
                              [equation], verbose=False)

    if itype == 'dw':
        equations.time_update(ts, domain.regions,
                              ebcs, epbcs, lcbcs, functions)

        if dw_mode == 'vector':
            residual = variables.create_stripped_state_vector()
            assemble_vector(residual, equation,
                            chunk_size, group_can_fail=False, **kwargs)
            if variables.has_lcbc:
                op_lcbc = variables.op_lcbc
                residual = op_lcbc.T * residual
            ret_val = residual

        elif dw_mode == 'matrix':
            if tangent_matrix is None:
                tangent_matrix = equations.create_matrix_graph()

            tangent_matrix.data[:] = 0.0
            assemble_matrix(tangent_matrix, equation,
                            chunk_size, group_can_fail=False, **kwargs)
            if variables.has_lcbc:
                op_lcbc = variables.op_lcbc
                tangent_matrix = op_lcbc.T * tangent_matrix * op_lcbc
                tangent_matrix = tangent_matrix.tocsr()
                tangent_matrix.sort_indices()
            ret_val = tangent_matrix

        else:
            print dw_mode
            raise ValueError

    elif itype == 'd':
        kwargs.setdefault( 'call_mode', 'd_eval' )

        val = 0.0
        for term in equation.terms:
            for ig in term.iter_groups():
                for aux, iels, status in term(chunk_size=chunk_size, **kwargs):
                    val += term.sign * aux
            ret_val = val

    elif itype == 'di':

        val = None
        for term in equation.terms:
            for ig in term.iter_groups():
                for aux, iels, status in term(chunk_size=chunk_size, **kwargs):
                    if val is None:
                        val = term.sign * aux
                    else:
                        val += term.sign * aux
            ret_val = val

    elif (itype == 'de') or (itype == 'dq'):

        val = None
        for term in equation.terms:
            for ig in term.iter_groups():
                for aux, iels, status in term(chunk_size=chunk_size, **kwargs):
                    if val is None:
                        val = term.sign * aux
                    else:
                        val = nm.concatenate( (val, term.sign * aux), axis = 0 )
        ret_val = val

    else:
        raise NotImplementedError, 'unknown term integration type: %s' % itype

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
def eval_residuals(state, equations, chunk_size=1000, **kwargs):
    equations.invalidate_term_caches()

    equations.set_variables_from_state(state)
    residual = equations.create_stripped_state_vector()

    for equation in equations:
        assemble_vector(residual, equation,
                        chunk_size=chunk_size, **kwargs)

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
def eval_tangent_matrices(state, tangent_matrix, equations,
                          chunk_size=1000, **kwargs):
    equations.set_variables_from_state(state)

    for equation in equations:
        assemble_matrix(tangent_matrix, equation,
                        chunk_size=chunk_size, **kwargs)

    return tangent_matrix

def assemble_by_blocks(conf_equations, problem, ebcs=None, epbcs=None,
                       dw_mode='matrix'):
    """Instead of a global matrix, return its building blocks as defined in
    `conf_equations`. The name and row/column variables of each block have to
    be encoded in the equation's name, as in:

    ebcs, epbcs must be either lists of BC names, or BC configuration
    dictionaries 

    Caveat: problem.variables must already contain all the variables used in
    conf_equations!

    conf_equations = {
      'A,v,u' : "dw_lin_elastic_iso.i1.Y2( inclusion.lame, v, u )",
    }
    """
    if isinstance( ebcs, list ) and isinstance( epbcs, list ):
        bc_mode = 0
    elif isinstance( ebcs, dict ) and isinstance( epbcs, dict ):
        bc_mode = 1
    else:
        raise TypeError('bad BC!')
    
    dummy = problem.create_state_vector()

    indx = problem.variables.get_indx
    matrices = {}
    for key, mtx_term in conf_equations.iteritems():
        ks = key.split( ',' )
        mtx_name, var_names = ks[0], ks[1:]
        output( mtx_name, var_names )

        problem.set_equations( {'eq': mtx_term}, single_term=True )

        if bc_mode == 0:
            problem.select_bcs( ebc_names = ebcs, epbc_names = epbcs )
        else:
            problem.time_update( conf_ebc = ebcs, conf_epbc = epbcs )

        ir = indx( var_names[0], stripped = True, allow_dual = True )
        ic = indx( var_names[1], stripped = True, allow_dual = True )

        mtx = eval_term_op( dummy, mtx_term, problem, dw_mode = 'matrix' )
        matrices[mtx_name] = mtx[ir,ic]

    return matrices
