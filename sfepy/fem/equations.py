import time
from copy import copy

import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import output, assert_, get_default, iter_dict_of_lists
from sfepy.base.base import debug, OneTypeList, Container, Struct
from sfepy.fem import Materials, Variables, setup_dof_conns
from extmods.mesh import create_mesh_graph
from sfepy.terms import Terms, Term

"""
Note:
- create found materials, variables from configuration/input file data
  ... no - user should be able to create objects even if they are not
  used in equations
"""

def parse_definition(equation_def):
    """
    Parse equation definition string to create term description list.
    """
    from parseEq import create_bnf

    term_descs = []
    bnf = create_bnf(term_descs)
    try:
        bnf.parseString(equation_def)
    except:
        raise ValueError('cannot parse equation! (%s)' % equation_def)

    return term_descs

def get_expression_arg_names(expression, strip_dots=True):
    """
    Parse expression and return set of all argument names. For arguments
    with attribute-like syntax (e.g. materials), if `strip_dots` is
    True, only base argument names are returned.
    """
    args = ','.join(aux.args for aux in parse_definition(expression))
    args = [arg.strip() for arg in args.split(',')]

    if strip_dots:
        for ii, arg in enumerate(args[:]):
            aux = arg.split('.')
            if len(aux) == 2:
                args[ii] = aux[0]

    return set(args)

##
# 21.07.2006, c
class Equations( Container ):

    @staticmethod
    def from_conf(conf, variables, regions, materials, integrals,
                  setup=True, user=None,
                  make_virtual=False, verbose=True):

        objs = OneTypeList(Equation)

        conf = copy(conf)

        ii = 0
        for name, desc in conf.iteritems():
            if verbose:
                output('equation "%s":' %  name)
                output(desc)
            eq = Equation.from_desc(name, desc, variables, regions,
                                    materials, integrals, user=user)
            objs.append(eq)
            ii += 1

        obj = Equations(objs, setup=setup,
                        make_virtual=make_virtual, verbose=verbose)

        return obj

    def __init__(self, equations, setup=True,
                 make_virtual=False, verbose=True):
        Container.__init__(self, equations)

        self.variables = Variables(self.collect_variables())
        self.materials = Materials(self.collect_materials())

        self.domain = self.get_domain()

        self.active_bcs = set()

        if setup:
            self.setup(make_virtual=make_virtual, verbose=verbose)

    def get_domain(self):
        domain = None

        for eq in self:
            for term in eq.terms:
                if term.has_region:
                    domain = term.region.domain

        return domain

    def setup(self, make_virtual=False, verbose=True):
        self.collect_conn_info()

        # This uses the conn_info created above.
        self.dof_conns = {}
        setup_dof_conns(self.conn_info, dof_conns=self.dof_conns,
                        make_virtual=make_virtual, verbose=verbose)

    def collect_materials(self):
        """
        Collect materials present in the terms of all equations.
        """
        materials = []
        for eq in self:
            materials.extend(eq.collect_materials())

        # Make the list items unique.
        materials = list(set(materials))

        return materials

    def reset_materials(self):
        """
        Clear material data so that next materials.time_update() is
        performed even for stationary materials.
        """
        self.materials.reset()

    def collect_variables(self):
        """
        Collect variables present in the terms of all equations.
        """
        variables = []
        for eq in self:
            variables.extend(eq.collect_variables())

        # Make the list items unique.
        variables = list(set(variables))

        return variables

    def get_variable(self, name):
        var = self.variables.get(name,
                                 msg_if_none='unknown variable! (%s)' % name)
        return var

    def collect_conn_info(self):
        """
        Collect connectivity information as defined by the equations.
        """
        self.conn_info = {}

        for eq in self:
            eq.collect_conn_info(self.conn_info)

        ## print_structs(self.conn_info)
        ## pause()

        return self.conn_info

    def get_variable_names( self ):
        """Return the list of names of all variables used in equations."""
        vns = set()
        for eq in self:
            for term in eq.terms:
                vns.update( term.get_variable_names() )
        return list( vns )

    def invalidate_term_caches(self):
        """
        Invalidate evaluate caches of variables present in equations.
        """
        for var in self.variables:
            var.invalidate_evaluate_cache()

    def print_terms(self):
        """
        Print names of equations and their terms.
        """
        output('equations:')
        for eq in self:
            output('  %s:' % eq.name)
            for term in eq.terms:
                output('    %+.2e * %s.%d.%s(%s)'
                       % (term.sign, term.name, term.integral.order,
                          term.region.name, term.arg_str))

    def time_update(self, ts, ebcs=None, epbcs=None, lcbcs=None,
                    functions=None, problem=None, verbose=True):
        """
        Update the equations for current time step.

        The update involves creating the mapping of active DOFs from/to
        all DOFs for all state variables, the setup of linear
        combination boundary conditions operators and the setup of
        active DOF connectivities.

        Parameters
        ----------
        ts : TimeStepper instance
            The time stepper.
        ebcs : Conditions instance, optional
            The essential (Dirichlet) boundary conditions.
        epbcs : Conditions instance, optional
            The periodic boundary conditions.
        lcbcs : Conditions instance, optional
            The linear combination boundary conditions.
        functions : Functions instance, optional
            The user functions for boundary conditions, materials, etc.
        problem : ProblemDefinition instance, optional
            The problem that can be passed to user functions as a context.
        verbose : bool
            If False, reduce verbosity.

        Returns
        -------
        graph_changed : bool
            The flag set to True if the current time step set of active
            boundary conditions differs from the set of the previous
            time step.
        """
        self.variables.time_update(ts, functions, verbose=verbose)

        active_bcs = self.variables.equation_mapping(ebcs, epbcs, ts, functions,
                                                     problem=problem)
        graph_changed = active_bcs != self.active_bcs
        self.active_bcs = active_bcs

        self.variables.setup_lcbc_operators(lcbcs, ts, functions)
        self.variables.setup_adof_conns()

        for eq in self:
            for term in eq.terms:
                term.time_update(ts)

        return graph_changed

    def time_update_materials(self, ts, problem=None, verbose=True):
        """
        Update data materials for current time and possibly also state.
        """
        self.materials.time_update(ts, self, problem=problem, verbose=verbose)

    def setup_initial_conditions(self, ics, functions):
        self.variables.setup_initial_conditions(ics, functions)

    def get_graph_conns(self, any_dof_conn=False, rdcs=None, cdcs=None):
        """
        Get DOF connectivities needed for creating tangent matrix graph.

        Parameters
        ----------
        any_dof_conn : bool
            By default, only volume DOF connectivities are used, with
            the exception of trace surface DOF connectivities. If True,
            any kind of DOF connectivities is allowed.
        rdcs, cdcs : arrays, optional
            Additional row and column DOF connectivities, corresponding
            to the variables used in the equations.

        Returns
        -------
        rdcs, cdcs : arrays
            The row and column DOF connectivities defining the matrix
            graph blocks.
        """
        if rdcs is None:
            rdcs = []
            cdcs = []

        elif cdcs is None:
            cdcs = copy(rdcs)

        else:
            assert_(len(rdcs) == len(cdcs))
            if rdcs is cdcs: # Make sure the lists are not the same object.
                rdcs = copy(rdcs)

        adcs = self.variables.adof_conns

        # Only volume dof connectivities are used, with the exception of trace
        # surface dof connectivities.
        shared = set()
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
            rvar, cvar = info.virtual, info.state
            if (rvar is None) or (cvar is None):
                continue

            is_surface = rvar.is_surface or cvar.is_surface

            dct = info.dc_type.type
            if not (dct in ('volume', 'scalar') or is_surface
                    or info.is_trace or any_dof_conn):
                continue


            rreg_name = info.get_region_name(can_trace=False)
            creg_name = info.get_region_name()

            for rig, cig in info.iter_igs():
                rname = rvar.get_primary_name()
                rkey = (rname, rreg_name, dct, rig, False)
                ckey = (cvar.name, creg_name, dct, cig, info.is_trace)

                dc_key = (rkey, ckey)
                ## print dc_key

                if not dc_key in shared:
                    try:
                        rdcs.append(adcs[rkey])
                        cdcs.append(adcs[ckey])
                    except:
                        debug()
                    shared.add(dc_key)

        return rdcs, cdcs

    def create_matrix_graph(self, any_dof_conn=False, rdcs=None, cdcs=None,
                            shape=None):
        """
        Create tangent matrix graph, i.e. preallocate and initialize the
        sparse storage needed for the tangent matrix. Order of DOF
        connectivities is not important.

        Parameters
        ----------
        any_dof_conn : bool
            By default, only volume DOF connectivities are used, with
            the exception of trace surface DOF connectivities. If True,
            any kind of DOF connectivities is allowed.
        rdcs, cdcs : arrays, optional
            Additional row and column DOF connectivities, corresponding
            to the variables used in the equations.
        shape : tuple, optional
            The required shape, if it is different from the shape
            determined by the equations variables. This may be needed if
            additional row and column DOF connectivities are passed in.

        Returns
        -------
        matrix : csr_matrix
            The matrix graph in the form of a CSR matrix with
            preallocated structure and zero data.
        """
        if not self.variables.has_virtuals():
            output('no matrix (no test variables)!')
            return None

        shape = get_default(shape, self.variables.get_matrix_shape())

        output( 'matrix shape:', shape )
        if nm.prod( shape ) == 0:
            output( 'no matrix (zero size)!' )
            return None

        rdcs, cdcs = self.get_graph_conns(any_dof_conn=any_dof_conn,
                                          rdcs=rdcs, cdcs=cdcs)

        if not len(rdcs):
            output('no matrix (empty dof connectivities)!')
            return None

        output( 'assembling matrix graph...' )
        tt = time.clock()

        nnz, prow, icol = create_mesh_graph(shape[0], shape[1],
                                            len(rdcs), rdcs, cdcs)

        output( '...done in %.2f s' % (time.clock() - tt) )
        output( 'matrix structural nonzeros: %d (%.2e%% fill)' \
                % (nnz, float( nnz ) / nm.prod( shape ) ) )
        ## print ret, prow, icol, nnz

        data = nm.zeros( (nnz,), dtype = self.variables.dtype )
        matrix = sp.csr_matrix( (data, icol, prow), shape )
        ## matrix.save( 'matrix', format = '%d %d %e\n' )
        ## pause()

        return matrix

    ##
    # c: 02.04.2008, r: 02.04.2008
    def init_time( self, ts ):
        pass

    ##
    # 08.06.2007, c
    def advance( self, ts ):
        for eq in self:
            for term in eq.terms:
                term.advance(ts)

        self.variables.advance(ts)

    ##
    # Interface to self.variables.
    def create_state_vector(self):
        return self.variables.create_state_vector()

    def create_stripped_state_vector(self):
        return self.variables.create_stripped_state_vector()

    def strip_state_vector(self, vec, follow_epbc=False):
        """
        Strip a full vector by removing EBC dofs.

        Notes
        -----
        If 'follow_epbc' is True, values of EPBC master dofs are not simply
        thrown away, but added to the corresponding slave dofs, just like when
        assembling. For vectors with state (unknown) variables it should be set
        to False, for assembled vectors it should be set to True.
        """
        return self.variables.strip_state_vector(vec, follow_epbc=follow_epbc)

    def make_full_vec(self, svec, force_value=None):
        """
        Make a full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.
        """
        return self.variables.make_full_vec(svec, force_value)

    def set_variables_from_state(self, vec, step=0):
        """
        Set data (vectors of DOF values) of variables.

        Parameters
        ----------
        data : array
            The state vector.
        step : int
            The time history step, 0 (default) = current.
        """
        self.variables.set_data(vec, step=step)

    def get_state_parts(self, vec=None):
        """
        Return parts of a state vector corresponding to individual state
        variables.

        Parameters
        ----------
        vec : array, optional
            The state vector. If not given, then the data stored in the
            variables are returned instead.

        Returns
        -------
        out : dict
            The dictionary of the state parts.
        """
        return self.variables.get_state_parts(vec)

    def set_data(self, data, step=0, ignore_unknown=False):
        """
        Set data (vectors of DOF values) of variables.

        Parameters
        ----------
        data : array
            The dictionary of {variable_name : data vector}.
        step : int, optional
            The time history step, 0 (default) = current.
        ignore_unknown : bool, optional
            Ignore unknown variable names if `data` is a dict.
        """
        self.variables.set_data(data, step=step,
                                ignore_unknown=ignore_unknown)

    def apply_ebc(self, vec, force_values=None):
        """
        Apply essential (Dirichlet) boundary conditions to a state vector.
        """
        self.variables.apply_ebc(vec, force_values=force_values)

    def apply_ic(self, vec, force_values=None):
        """
        Apply initial conditions to a state vector.
        """
        self.variables.apply_ic(vec, force_values=force_values)

    def state_to_output(self, vec, fill_value=None, var_info=None,
                        extend=True):
        return self.variables.state_to_output(vec,
                                              fill_value=fill_value,
                                              var_info=var_info,
                                              extend=extend)

    def get_lcbc_operator(self):
        return self.variables.get_lcbc_operator()

    def evaluate(self, mode='eval', dw_mode='vector', term_mode=None,
                 asm_obj=None):
        """
        Parameters
        ----------
        mode : one of 'eval', 'el_avg', 'qp', 'weak'
            The evaluation mode.
        """
        if mode == 'weak':
            out = asm_obj

        else:
            out = {}

        for eq in self:
            eout = eq.evaluate(mode=mode, dw_mode=dw_mode, term_mode=term_mode,
                               asm_obj=asm_obj)
            if mode != 'weak':
                out[eq.name] = eout

        if (len(self) == 1) and (mode != 'weak'):
            out = out.popitem()[1]

        return out

    def eval_residuals(self, state, by_blocks=False, names=None):
        """
        Evaluate (assemble) residual vectors.

        Parameters
        ----------
        state : array
            The vector of DOF values. Note that it is needed only in
            nonlinear terms.
        by_blocks : bool
            If True, return the individual blocks composing the whole
            residual vector. Each equation should then correspond to one
            required block and should be named as `'block_name,
            test_variable_name, unknown_variable_name'`.
        names : list of str, optional
            Optionally, select only blocks with the given `names`, if
            `by_blocks` is True.

        Returns
        -------
        out : array or dict of array
            The assembled residual vector. If `by_blocks` is True, a
            dictionary is returned instead, with keys given by
            `block_name` part of the individual equation names.
        """
        self.set_variables_from_state(state)

        if by_blocks:
            names = get_default(names, self.names)

            out = {}

            get_indx = self.variables.get_indx
            for name in names:
                eq = self[name]
                key, rname, cname = [aux.strip()
                                     for aux in name.split(',')]

                ir = get_indx(rname, stripped=True, allow_dual=True)

                residual = self.create_stripped_state_vector()
                eq.evaluate(mode='weak', dw_mode='vector', asm_obj=residual)

                out[key] = residual[ir]

        else:
            out = self.create_stripped_state_vector()

            self.evaluate(mode='weak', dw_mode='vector', asm_obj=out)

        return out

    def eval_tangent_matrices(self, state, tangent_matrix,
                              by_blocks=False, names=None):
        """
        Evaluate (assemble) tangent matrices.

        Parameters
        ----------
        state : array
            The vector of DOF values. Note that it is needed only in
            nonlinear terms.
        tangent_matrix : csr_matrix
            The preallocated CSR matrix with zero data.
        by_blocks : bool
            If True, return the individual blocks composing the whole
            matrix. Each equation should then correspond to one
            required block and should be named as `'block_name,
            test_variable_name, unknown_variable_name'`.
        names : list of str, optional
            Optionally, select only blocks with the given `names`, if
            `by_blocks` is True.

        Returns
        -------
        out : csr_matrix or dict of csr_matrix
            The assembled matrix. If `by_blocks` is True, a dictionary
            is returned instead, with keys given by `block_name` part
            of the individual equation names.
        """
        self.set_variables_from_state(state)

        if by_blocks:
            names = get_default(names, self.names)

            out = {}

            get_indx = self.variables.get_indx
            for name in names:
                eq = self[name]
                key, rname, cname = [aux.strip()
                                     for aux in eq.name.split(',')]

                ir = get_indx(rname, stripped=True, allow_dual=True)
                ic = get_indx(cname, stripped=True, allow_dual=True)

                tangent_matrix.data[:] = 0.0
                eq.evaluate(mode='weak', dw_mode='matrix',
                            asm_obj=tangent_matrix)

                out[key] = tangent_matrix[ir, ic]

        else:
            tangent_matrix.data[:] = 0.0

            self.evaluate(mode='weak', dw_mode='matrix', asm_obj=tangent_matrix)

            out = tangent_matrix

        return out

##
# 21.07.2006, c
class Equation( Struct ):

    @staticmethod
    def from_desc(name, desc, variables, regions, materials, integrals,
                  user=None):
        term_descs = parse_definition(desc)
        terms = Terms.from_desc(term_descs, regions, integrals)

        terms.setup()
        terms.assign_args(variables, materials, user)

        obj = Equation(name, terms)

        return obj

    def __init__(self, name, terms):
        Struct.__init__(self, name = name)

        if isinstance(terms, Term): # single Term
            terms = Terms([terms])

        self.terms = terms

        self.terms.setup()

    def collect_materials(self):
        """
        Collect materials present in the terms of the equation.
        """
        materials = []
        for term in self.terms:
            materials.extend(term.get_materials(join=True))

        return materials

    def collect_variables(self):
        """
        Collect variables present in the terms of the equation.

        Ensures that corresponding primary variables of test/parameter
        variables are always in the list, even if they are not directly
        used in the terms.
        """
        variables = []
        for term in self.terms:
            var_names = term.get_variable_names()

            aux = term.get_args_by_name(var_names)
            for var in aux:
                variables.append(var)
                pvar = var.get_primary()
                if pvar is not None:
                    variables.append(pvar)

        return variables

    def collect_conn_info(self, conn_info):

        for term in self.terms:
            key = (self.name,) + term.get_conn_key()

            conn_info[key] = term.get_conn_info()

    def evaluate(self, mode='eval', dw_mode='vector', term_mode=None,
                 asm_obj=None):
        """
        Parameters
        ----------
        mode : one of 'eval', 'el_avg', 'qp', 'weak'
            The evaluation mode.
        """
        if mode == 'eval':
            val = 0.0
            for term in self.terms:
                aux, status = term.evaluate(mode=mode,
                                            term_mode=term_mode,
                                            standalone=False,
                                            ret_status=True)
                val += aux

            out = val

        elif mode in ('el_avg', 'el', 'qp'):

            vals = []
            for term in self.terms:
                val, iels, status = term.evaluate(mode=mode,
                                                  term_mode=term_mode,
                                                  standalone=False,
                                                  ret_status=True)
                vals.append(val)

            if len(vals) == 1:
                vals = vals[0]

            out = vals

        elif mode == 'weak':

            if dw_mode == 'vector':

                for term in self.terms:
                    val, iels, status = term.evaluate(mode=mode,
                                                      term_mode=term_mode,
                                                      standalone=False,
                                                      ret_status=True)
                    term.assemble_to(asm_obj, val, iels, mode=dw_mode)

            elif dw_mode == 'matrix':

                for term in self.terms:
                    svars = term.get_state_variables(unknown_only=True)

                    for svar in svars:
                        val, iels, status = term.evaluate(mode=mode,
                                                          term_mode=term_mode,
                                                          diff_var=svar.name,
                                                          standalone=False,
                                                          ret_status=True)
                        term.assemble_to(asm_obj, val, iels,
                                         mode=dw_mode, diff_var=svar)

            else:
                raise ValueError('unknown assembling mode! (%s)' % dw_mode)

            out = asm_obj

        else:
            raise ValueError('unknown evaluation mode! (%s)' % mode)

        return out
