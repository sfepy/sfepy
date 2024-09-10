"""
Classes of equations composed of terms.
"""
from copy import copy

import numpy as nm
import scipy.sparse as sp

from sfepy.base.base import output, assert_, get_default, iter_dict_of_lists
from sfepy.base.base import OneTypeList, Container, Struct
from sfepy.base.timing import Timer
from sfepy.discrete import Materials, Variables, create_adof_conns
from sfepy.discrete.common.extmods.cmesh import create_mesh_graph
from sfepy.terms import Terms, Term
from sfepy.terms.terms_multilinear import ETermBase

def parse_definition(equation_def):
    """
    Parse equation definition string to create term description list.
    """
    from .parse_equations import create_bnf

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

class Equations(Container):

    @staticmethod
    def from_conf(conf, variables, regions, materials, integrals,
                  user=None, eterm_options=None, allow_derivatives=False,
                  verbose=True):

        objs = OneTypeList(Equation)

        conf = copy(conf)

        ii = 0
        for name, desc in conf.items():
            if verbose:
                output('equation "%s":' %  name)
                output(desc)
            eq = Equation.from_desc(name, desc, variables, regions,
                                    materials, integrals, user=user,
                                    eterm_options=eterm_options,
                                    allow_derivatives=allow_derivatives)
            objs.append(eq)
            ii += 1

        obj = Equations(objs)

        return obj

    def __init__(self, equations):
        Container.__init__(self, equations)

        self.variables = Variables(self.collect_variables())
        self.materials = Materials(self.collect_materials())

        self.domain = self.get_domain()

        self.active_bcs = set()

        self.collect_conn_info()

    def add_equation(self, equation):
        """
        Add a new equation.

        Parameters
        ----------
        equation : Equation instance
                    The new equation.
        """
        self.append(equation)
        self.variables.extend(
             set(equation.collect_variables() ) - set(self.variables)
        )
        self.materials.extend(
             set(equation.collect_materials() ) - set(self.materials)
        )
        equation.collect_conn_info(self.conn_info)
        if not self.domain:
           self.domain = self.get_domain()

    def create_subequations(self, var_names, known_var_names=None):
        """
        Create sub-equations containing only terms with the given virtual
        variables.

        Parameters
        ----------
        var_names : list
            The list of names of virtual variables.
        known_var_names : list
            The list of  names of (already) known state variables.

        Returns
        -------
        subequations : Equations instance
            The sub-equations.
        """
        from sfepy.discrete import FieldVariable

        known_var_names = get_default(known_var_names, [])

        objs = []
        for iv, var_name in enumerate(var_names):
            terms = [term.copy(name=term.name)
                     for eq in self for term in eq.terms
                     if term.get_virtual_name() == var_name]

            # Make parameter variables from known state variables in terms
            # arguments.
            for known_name in known_var_names:
                for term in terms:
                    if known_name in term.arg_names:
                        ii = term.arg_names.index(known_name)
                        state = self.variables[known_name]
                        par = FieldVariable(known_name, 'parameter',
                                            state.field,
                                            primary_var_name='(set-to-None)')
                        term.args[ii] = par
                        term._kwargs[known_name] = par
                        par.set_data(state())

            new_terms = Terms(terms)
            objs.append(Equation('eq_%d' % iv, new_terms))

        subequations = Equations(objs)

        return subequations

    def get_domain(self):
        domain = None

        for eq in self:
            for term in eq.terms:
                if term.has_region:
                    domain = term.region.domain

        return domain

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

        return self.conn_info

    def get_variable_names(self):
        """
        Return the list of names of all variables used in equations.
        """
        vns = set()
        for eq in self:
            for term in eq.terms:
                vns.update(term.get_variable_names())
        return list(vns)

    def get_variable_dependencies(self):
        """
        For each virtual variable get names of state/parameter variables that
        are present in terms with that virtual variable.

        The virtual variables define the actual equations and their
        dependencies define the variables needed to evaluate the equations.

        Returns
        -------
        deps : dict
            The dependencies as a dictionary with virtual variable names as
            keys and sets of state/parameter variables as values.
        """
        deps = {}
        for eq in self:
            for term in eq.terms:
                dep_list = deps.setdefault(term.get_virtual_name(), set())
                dep_list.update(term.get_state_names())

        return deps

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
                output('    %s' % term.get_str())

    def time_update(self, ts, ebcs=None, epbcs=None, lcbcs=None,
                    functions=None, problem=None, active_only=True,
                    verbose=True):
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
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        active_only : bool
            If True, the active DOF connectivities and matrix graph have
            reduced size and are created with the reduced (active DOFs only)
            numbering.
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
                                                     problem=problem,
                                                     active_only=active_only)
        graph_changed = active_bcs != self.active_bcs
        self.active_bcs = active_bcs

        if graph_changed or not self.variables.adof_conns:
            adcs = create_adof_conns(self.conn_info, self.variables.adi.indx,
                                     active_only=active_only)
            self.variables.set_adof_conns(adcs)

        self.variables.setup_lcbc_operators(lcbcs, ts, functions)

        for eq in self:
            for term in eq.terms:
                term.time_update(ts)

        return graph_changed

    def time_update_materials(self, ts, mode='normal', problem=None,
                              verbose=True):
        """
        Update data materials for current time and possibly also state.

        Parameters
        ----------
        ts : TimeStepper instance
            The time stepper.
        mode : 'normal', 'update' or 'force'
            The update mode, see
            :func:`sfepy.discrete.materials.Material.time_update()`.
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        verbose : bool
            If False, reduce verbosity.
        """
        self.materials.time_update(ts, self, mode=mode, problem=problem,
                                   verbose=verbose)

    def setup_initial_conditions(self, ics, functions=None):
        self.variables.setup_initial_conditions(ics, functions)

    def get_graph_conns(self, any_dof_conn=False, rdcs=None, cdcs=None,
                        active_only=True):
        """
        Get DOF connectivities needed for creating tangent matrix graph.

        Parameters
        ----------
        any_dof_conn : bool
            By default, only cell DOF connectivities are used, with
            the exception of trace facet DOF connectivities. If True,
            any kind of DOF connectivities is allowed.
        rdcs, cdcs : arrays, optional
            Additional row and column DOF connectivities, corresponding
            to the variables used in the equations.
        active_only : bool
            If True, the active DOF connectivities have reduced size and are
            created with the reduced (active DOFs only) numbering.

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

        # Only cell dof connectivities are used, with the exception of trace
        # facet dof connectivities.
        shared = set()
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
            rvar, cvar = info.virtual, info.state
            if (rvar is None) or (cvar is None):
                continue

            is_surface = rvar.is_surface or cvar.is_surface

            rreg_name = info.get_region_name(can_trace=False)
            creg_name = info.get_region_name()
            mreg_name = None if creg_name == rreg_name else rreg_name

            rname = rvar.get_primary_name()
            cname = cvar.get_primary_name()
            rkey = (rname, rreg_name, info.dof_conn_types[rname], None)
            ckey = (cvar.name, creg_name, info.dof_conn_types[cname],
                    mreg_name)

            dc_key = (rkey, ckey)

            if dc_key not in shared:
                rdc = adcs[rkey]
                cdc = adcs[ckey]
                if not active_only:
                    ii = nm.where(rdc < 0)
                    rdc = rdc.copy()
                    rdc[ii] = -1 - rdc[ii]

                    ii = nm.where(cdc < 0)
                    cdc = cdc.copy()
                    cdc[ii] = -1 - cdc[ii]

                rdcs.append(rdc)
                cdcs.append(cdc)

                shared.add(dc_key)

        return rdcs, cdcs

    def create_matrix_graph(self, any_dof_conn=False, rdcs=None, cdcs=None,
                            shape=None, active_only=True, verbose=True):
        """
        Create tangent matrix graph, i.e. preallocate and initialize the
        sparse storage needed for the tangent matrix. Order of DOF
        connectivities is not important.

        Parameters
        ----------
        any_dof_conn : bool
            By default, only cell region DOF connectivities are used, with
            the exception of trace facet DOF connectivities. If True,
            any DOF connectivities are used.
        rdcs, cdcs : arrays, optional
            Additional row and column DOF connectivities, corresponding
            to the variables used in the equations.
        shape : tuple, optional
            The required shape, if it is different from the shape
            determined by the equations variables. This may be needed if
            additional row and column DOF connectivities are passed in.
        active_only : bool
            If True, the matrix graph has reduced size and is created with the
            reduced (active DOFs only) numbering.
        verbose : bool
            If False, reduce verbosity.

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

        output('matrix shape:', shape, verbose=verbose)
        size = nm.prod(shape, dtype=nm.int64)
        if size == 0:
            output('no matrix (zero size)!')
            return None

        rdcs, cdcs = self.get_graph_conns(any_dof_conn=any_dof_conn,
                                          rdcs=rdcs, cdcs=cdcs,
                                          active_only=active_only)

        if not len(rdcs):
            output('no matrix (empty dof connectivities)!')
            return None

        output('assembling matrix graph...', verbose=verbose)
        timer = Timer(start=True)

        nnz, prow, icol = create_mesh_graph(shape[0], shape[1],
                                            len(rdcs), rdcs, cdcs)

        output('...done in %.2f s' % timer.stop(), verbose=verbose)
        output('matrix structural nonzeros: %d (%.2e%% fill)' \
               % (nnz, 100.0 * float(nnz) / size), verbose=verbose)

        data = nm.zeros((nnz,), dtype=self.variables.dtype)
        matrix = sp.csr_matrix((data, icol, prow), shape)

        return matrix

    def init_time(self, ts):
        pass

    def advance(self, ts):
        for eq in self:
            for term in eq.terms:
                term.advance(ts)

        self.variables.advance(ts)

    ##
    # Interface to self.variables.
    def create_vec(self):
        return self.variables.create_vec()

    def create_reduced_vec(self):
        return self.variables.create_reduced_vec()

    def reduce_vec(self, vec, follow_epbc=False):
        """
        Get the reduced DOF vector, with EBC and PBC DOFs removed.

        Notes
        -----
        If 'follow_epbc' is True, values of EPBC master dofs are not simply
        thrown away, but added to the corresponding slave dofs, just like when
        assembling. For vectors with state (unknown) variables it should be set
        to False, for assembled vectors it should be set to True.
        """
        return self.variables.reduce_vec(vec, follow_epbc=follow_epbc)

    def make_full_vec(self, svec, force_value=None):
        """
        Make a full DOF vector satisfying E(P)BCs from a reduced DOF
        vector.
        """
        return self.variables.make_full_vec(svec, force_value)

    def set_data(self, data, step=0, ignore_unknown=False):
        """
        Set data (vectors of DOF values) of variables.

        Parameters
        ----------
        data : dict
            The dictionary of {variable_name : data vector}.
        step : int, optional
            The time history step, 0 (default) = current.
        ignore_unknown : bool, optional
            Ignore unknown variable names if `data` is a dict.
        """
        self.variables.set_data(data, step=step,
                                ignore_unknown=ignore_unknown)

    def init_state(self, vec=None):
        self.variables.init_state(vec=vec)

    def apply_ebc(self, vec=None, force_values=None):
        """
        Apply essential (Dirichlet) boundary conditions to equations' variables,
        or a given vector.
        """
        self.variables.apply_ebc(vec=vec, force_values=force_values)

    def apply_ic(self, vec=None, force_values=None):
        """
        Apply initial conditions to equations' variables, or a given vector.
        """
        self.variables.apply_ic(vec=vec, force_values=force_values)

    def set_state(self, vec, reduced=False, force=False, preserve_caches=False):
        self.variables.set_state(vec, reduced=reduced, force=force,
                                 preserve_caches=preserve_caches)

    def get_lcbc_operator(self):
        return self.variables.get_lcbc_operator()

    def evaluate(self, names=None, mode='eval', dw_mode='vector',
                 term_mode=None, diff_vars=None, asm_obj=None,
                 select_term=None):
        """
        Evaluate the equations.

        Parameters
        ----------
        names : str or sequence of str, optional
            Evaluate only equations of the given name(s).
        mode : one of 'eval', 'el_avg', 'qp', 'weak'
            The evaluation mode.
        dw_mode : one of 'vector', 'matrix', 'sensitivity'
            The particular evaluation mode if `mode` is ``'weak'``.
        term_mode : str
            The term evaluation mode, used mostly if `mode` is ``'eval'`` in
            some terms.
        diff_vars : list of str
            The names of parameters with respect to the equations are
            differentiated if `dw_mode` is ``'sensitivity'``.
        asm_obj : ndarray or spmatrix
            The object for storing the evaluation result in the ``'weak'`` mode.
        select_term : function(term)
            Optional boolean function returning True for terms that should be
            evaluated.

        Returns
        -------
        out : dict or result
            The evaluation result. In 'weak' mode it is the
            `asm_obj`. Otherwise, it is a dict of results with equation names
            as keys or a single result for a single equation.
        """
        if names is None:
            eqs = self
            single = (len(eqs) == 1)

        else:
            single = isinstance(names, str)
            if single:
                names = [names]

            eqs = [self[eq] for eq in names]

        if mode == 'weak':
            extras = []
            for eq in eqs:
                out = eq.evaluate(mode=mode, dw_mode=dw_mode,
                                  term_mode=term_mode, diff_vars=diff_vars,
                                  asm_obj=asm_obj,
                                  select_term=select_term)
                if isinstance(out, tuple): extras.extend(out[1])

            out = asm_obj
            for extra in extras:
                out = out + extra

        else:
            out = {}
            for eq in eqs:
                eout = eq.evaluate(mode=mode, dw_mode=dw_mode,
                                   term_mode=term_mode,
                                   select_term=select_term)
                out[eq.name] = eout

            if single:
                out = out.popitem()[1]

        return out

    def eval_residuals(self, state, by_blocks=False, names=None,
                       select_term=None):
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
        select_term : function(term)
            Optional boolean function returning True for terms that should be
            evaluated.

        Returns
        -------
        out : array or dict of array
            The assembled residual vector. If `by_blocks` is True, a
            dictionary is returned instead, with keys given by
            `block_name` part of the individual equation names.
        """
        self.set_state(state, force=True)

        if by_blocks:
            names = get_default(names, self.names)

            out = {}

            get_indx = self.variables.get_indx
            for name in names:
                eq = self[name]
                key, rname, cname = [aux.strip()
                                     for aux in name.split(',')]

                ir = get_indx(rname, reduced=True, allow_dual=True)

                residual = self.create_reduced_vec()
                eq.evaluate(mode='weak', dw_mode='vector', asm_obj=residual,
                            select_term=select_term)

                out[key] = residual[ir]

        else:
            out = self.create_reduced_vec()

            self.evaluate(mode='weak', dw_mode='vector', asm_obj=out,
                          select_term=select_term)

        return out

    def eval_tangent_matrices(self, state, tangent_matrix,
                              by_blocks=False, names=None, select_term=None):
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
        select_term : function(term)
            Optional boolean function returning True for terms that should be
            evaluated.

        Returns
        -------
        out : csr_matrix or dict of csr_matrix
            The assembled matrix. If `by_blocks` is True, a dictionary
            is returned instead, with keys given by `block_name` part
            of the individual equation names.
        """
        self.set_state(state, force=True)

        if by_blocks:
            names = get_default(names, self.names)

            out = {}

            get_indx = self.variables.get_indx
            for name in names:
                eq = self[name]
                key, rname, cname = [aux.strip()
                                     for aux in eq.name.split(',')]

                ir = get_indx(rname, reduced=True, allow_dual=True)
                ic = get_indx(cname, reduced=True, allow_dual=True)

                tangent_matrix.data[:] = 0.0
                aux = eq.evaluate(mode='weak', dw_mode='matrix',
                                  asm_obj=tangent_matrix,
                                  select_term=select_term)

                out[key] = aux[ir, ic]

        else:
            tangent_matrix.data[:] = 0.0

            out = self.evaluate(mode='weak', dw_mode='matrix',
                                asm_obj=tangent_matrix,
                                select_term=select_term)

        return out

class Equation(Struct):

    @staticmethod
    def from_desc(name, desc, variables, regions, materials, integrals,
                  user=None, eterm_options=None, allow_derivatives=False):
        term_descs = parse_definition(desc)
        terms = Terms.from_desc(term_descs, regions, integrals)

        terms.setup(allow_derivatives=allow_derivatives)
        terms.assign_args(variables, materials, user)

        if eterm_options is not None:
            for term in terms:
                if isinstance(term, ETermBase):
                    term.set_verbosity(eterm_options.get('verbosity', 0))
                    term.set_backend(**eterm_options.get('backend_args', {}))

        obj = Equation(name, terms, setup=False)

        return obj

    def __init__(self, name, terms, setup=True):
        Struct.__init__(self, name=name)

        if isinstance(terms, Term): # A single term.
            terms = Terms([terms])

        self.terms = terms

        if setup:
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
                 diff_vars=None, asm_obj=None, select_term=None):
        """
        Evaluate the equation.

        Parameters
        ----------
        mode : one of 'eval', 'el_eval', 'el_avg', 'qp', 'weak'
            The evaluation mode.
        dw_mode : one of 'vector', 'matrix', 'sensitivity'
            The particular evaluation mode if `mode` is ``'weak'``.
        term_mode : str
            The term evaluation mode, used mostly if `mode` is ``'eval'`` in
            some terms.
        diff_vars : list of str
            The names of parameters with respect to the equation is
            differentiated if `dw_mode` is ``'sensitivity'``.
        asm_obj : ndarray or spmatrix
            The object for storing the evaluation result in the ``'weak'`` mode.
        select_term : function(term)
            Optional boolean function returning True for terms that should be
            evaluated.

        Returns
        -------
        out : result
            The evaluation result. In 'weak' mode it is the
            `asm_obj`.
        """
        terms = self.terms
        if select_term is not None:
            terms = [term for term in self.terms if select_term(term)]

        if mode in ('eval', 'el_eval', 'el_avg', 'qp'):
            val = 0.0
            for term in terms:
                aux, status = term.evaluate(mode=mode,
                                            term_mode=term_mode,
                                            standalone=False,
                                            ret_status=True)
                val += aux

            out = val

        elif mode == 'weak':

            if dw_mode == 'vector':

                for term in terms:
                    val, iels, status = term.evaluate(mode=mode,
                                                      term_mode=term_mode,
                                                      standalone=False,
                                                      ret_status=True)
                    term.assemble_to(asm_obj, val, iels, mode=dw_mode)

                out = asm_obj

            elif dw_mode == 'matrix':

                extras = []
                for term in terms:
                    svars = term.get_state_variables(unknown_only=True)

                    for svar in svars:
                        val, iels, status = term.evaluate(mode=mode,
                                                          term_mode=term_mode,
                                                          diff_var=svar.name,
                                                          standalone=False,
                                                          ret_status=True)
                        extra = term.assemble_to(asm_obj, val, iels,
                                                 mode=dw_mode, diff_var=svar)
                        if extra is not None: extras.append(extra)

                out = (asm_obj, extras) if len(extras) else asm_obj

            elif dw_mode == 'sensitivity':
                # Differentiation w.r.t. material parameters.
                if diff_vars is None: diff_vars = ()

                for ic, diff_var in enumerate(diff_vars):
                    for term in terms:
                        if not (term.diff_info and
                                (diff_var in term.get_material_names(part=1))):
                            continue
                        val, iels, status = term.evaluate(mode=mode,
                                                          term_mode=term_mode,
                                                          diff_var=diff_var,
                                                          standalone=False,
                                                          ret_status=True)
                        term.assemble_to(asm_obj[:, ic], val, iels)

                out = asm_obj

            else:
                raise ValueError('unknown assembling mode! (%s)' % dw_mode)

        else:
            raise ValueError('unknown evaluation mode! (%s)' % mode)

        return out
