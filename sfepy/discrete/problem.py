import os
import os.path as op
from copy import copy
from functools import partial

import numpy as nm

from sfepy.base.base import (
    dict_from_keys_init, select_by_names, is_string, is_integer, is_sequence,
    output, get_default, Struct, IndexedStruct)
import sfepy.base.ioutils as io
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.conf import transform_variables, transform_materials
from sfepy.base.timing import Timer
from .functions import Functions
from sfepy.discrete.fem.mesh import Mesh
from sfepy.discrete.fem.meshio import check_format_suffix
from sfepy.discrete.fem.fields_base import set_mesh_coors
from sfepy.discrete.common.fields import fields_from_conf
from .variables import Variables, Variable
from .materials import Materials, Material
from .equations import Equations
from .integrals import Integrals
from sfepy.discrete.conditions import Conditions
from sfepy.discrete.evaluate import create_evaluable, eval_equations
from sfepy.solvers.ts import TimeStepper
from sfepy.discrete.evaluate import Evaluator
from sfepy.solvers import Solver, NonlinearSolver
from sfepy.solvers.solvers import use_first_available
from sfepy.solvers.ts_solvers import (StationarySolver, ElastodynamicsBaseTS,
                                      transform_equations_ed)

def make_is_save(options):
    """
    Given problem options, return a callable that determines whether to save
    results of a time step.
    """
    class IsSave(Struct):
        def __init__(self, save_times):
            if is_sequence(save_times):
                save_times = nm.asarray(save_times)

            self.save_times0 = save_times
            self.reset()

        def reset(self, ts=None):
            self.ilast = 0
            self.save_times = self.save_times0
            if ts is not None:
                if is_integer(self.save_times0):
                    self.save_times = nm.linspace(ts.t0, ts.t1,
                                                  self.save_times0)

        def __call__(self, ts):
            if is_string(self.save_times) and self.save_times == 'all':
                return True

            elif isinstance(self.save_times, nm.ndarray):
                if (self.ilast < len(self.save_times)
                    and (ts.time + (1e-14 * ts.dt)
                         >= self.save_times[self.ilast])):
                    self.ilast += 1
                    return True

            elif callable(self.save_times):
                return self.save_times(ts)

            return False

    save_times = options.get('save_times', 'all')
    is_save = IsSave(save_times)

    return is_save

def prepare_matrix(problem, state):
    """
    Pre-assemble tangent system matrix.
    """
    problem.update_materials()

    ev = problem.get_evaluator()
    try:
        mtx = ev.eval_tangent_matrix(state(), is_full=True)

    except ValueError:
        output('matrix evaluation failed, giving up...')
        raise

    return mtx

##
# 29.01.2006, c
class Problem(Struct):
    """
    Problem definition, the top-level class holding all data necessary to solve
    a problem.

    It can be constructed from a :class:`ProblemConf
    <sfepy.base.conf.ProblemConf>` instance using `Problem.from_conf()` or
    directly from a problem description file using `Problem.from_conf_file()`

    For interactive use, the constructor requires only the `equations`,
    `nls` and `ls` keyword arguments, see below.

    Parameters
    ----------
    name : str
        The problem name.
    conf : ProblemConf instance, optional
        The :class:`ProblemConf <sfepy.base.conf.ProblemConf>` describing the
        problem.
    functions : Functions instance, optional
        The user functions for boundary conditions, materials, etc.
    domain : Domain instance, optional
        The solution :class:`Domain <sfepy.discrete.common.domain.Domain>`.
    fields : dict, optional
        The dictionary of :class:`Field <sfepy.discrete.common.fields.Field>`
        instances.
    equations : Equations instance, optional
        The :class:`Equations <sfepy.discrete.equations.Equations>` to solve.
        This argument is required when `auto_conf` is True.
    auto_conf : bool
        If True, fields and domain are determined by `equations`.
    active_only : bool
        If True, the (tangent) matrices and residual vectors (right-hand sides)
        contain only active DOFs, see below.

    Notes
    -----
    The Problem is by default created with `active_only` set to True. Then the
    (tangent) matrices and residual vectors (right-hand sides) have reduced
    sizes and contain only the active DOFs, i.e., DOFs not constrained by EBCs
    or EPBCs.

    Setting `active_only` to False results in full-size vectors and
    matrices. Then the matrix size non-zeros structure does not depend on the
    actual E(P)BCs applied. It must be False when using parallel PETSc solvers.

    The active DOF connectivities contain all DOFs, with the E(P)BC-constrained
    ones stored as `-1 - <DOF number>`, so that the full connectivities can be
    reconstructed for the matrix graph creation. However, the negative entries
    mean that the assembled matrices/residuals have zero values at positions
    corresponding to constrained DOFs.

    The resulting linear system then provides a solution increment, that has to
    be added to the initial guess used to compute the residual, just like in
    the Newton iterations. The increment of the constrained DOFs is
    automatically zero.

    When solving with a direct solver, the diagonal entries of a matrix at
    positions corresponding to constrained DOFs has to be set to ones, so that
    the matrix is not singular, see
    :func:`sfepy.discrete.evaluate.apply_ebc_to_matrix()`, which is called
    automatically in
    :func:`sfepy.discrete.evaluate.Evaluator.eval_tangent_matrix()`. It is
    not called automatically in :func:`Problem.evaluate()`. Note that setting
    the diagonal entries to one might not be necessary with iterative solvers,
    as the zero matrix rows match the zero residual rows, i.e. if the reduced
    matrix would be regular, then the right-hand side (the residual) is
    orthogonal to the kernel of the matrix.
    """

    @staticmethod
    def from_conf_file(conf_filename, required=None, other=None,
                       init_fields=True, init_equations=True,
                       init_solvers=True):

        _required, _other = get_standard_keywords()
        if required is None:
            required = _required
        if other is None:
            other = _other

        conf = ProblemConf.from_file(conf_filename, required, other)

        obj = Problem.from_conf(conf, init_fields=init_fields,
                                init_equations=init_equations,
                                init_solvers=init_solvers)
        return obj

    @staticmethod
    def from_conf(conf, init_fields=True, init_equations=True,
                  init_solvers=True):
        if conf.options.get('absolute_mesh_path', False):
            conf_dir = None
        else:
            conf_dir = op.dirname(conf.funmod.__file__)

        functions = Functions.from_conf(conf.functions)

        if conf.get('filename_mesh') is not None:
            from sfepy.discrete.fem.domain import FEDomain

            mesh = Mesh.from_file(conf.filename_mesh, prefix_dir=conf_dir)
            domain = FEDomain(mesh.name, mesh)

            refine = conf.options.get('refinement_level', 0)
            if refine > 0:
                for ii in range(refine):
                    output('refine %d...' % ii)
                    domain = domain.refine()
                    output('... %d nodes %d elements'
                           % (domain.shape.n_nod, domain.shape.n_el))

            if conf.options.get('ulf', False):
                domain.mesh.coors_act = domain.mesh.coors.copy()

            if conf.options.get('mesh_eps') is not None:
                import sfepy.discrete.fem.mesh as msh
                import sfepy.discrete.fem.periodic as per
                msh.set_accuracy(conf.options.mesh_eps)
                per.set_accuracy(conf.options.mesh_eps)

        elif conf.get('filename_domain') is not None:
            from sfepy.discrete.iga.domain import IGDomain
            domain = IGDomain.from_file(conf.filename_domain)

        else:
            raise ValueError('missing filename_mesh or filename_domain!')

        active_only = conf.options.get('active_only', True)
        obj = Problem('problem_from_conf', conf=conf, functions=functions,
                      domain=domain, auto_conf=False,
                      active_only=active_only)

        allow_empty = conf.options.get('allow_empty_regions', False)
        obj.set_regions(conf.regions, obj.functions,
                        allow_empty=allow_empty)

        obj.clear_equations()

        if init_fields:
            obj.set_fields(conf.fields)

            if init_equations:
                obj.set_equations(conf.equations)

        if init_solvers:
            obj.set_conf_solvers(conf.solvers, conf.options)

        return obj

    def __init__(self, name, conf=None, functions=None,
                 domain=None, fields=None, equations=None, auto_conf=True,
                 active_only=True):
        if conf is None:
            self.conf = Struct(options={}, ics={},
                               ebcs={}, epbcs={}, lcbcs={}, materials={},
                               solvers={})

        else:
            self.conf = conf.copy()
            self.conf.set_default('options', {})

        self.active_only = active_only
        self.name = name
        self.functions = functions

        self.reset()
        self.ls_conf = self.nls_conf = self.tsc_conf = self.ts_conf = None
        self.conf_variables = self.conf_materials = None

        if auto_conf:
            if equations is None:
                raise ValueError('missing equations in auto_conf mode!')

            if fields is None:
                variables = equations.variables
                fields = {}
                for field in [var.get_field() for var in variables]:
                    fields[field.name] = field

            if domain is None:
                domain = list(fields.values())[0].domain

        self.equations = equations
        self.fields = fields
        self.domain = domain
        if auto_conf:
            self.set_ics(self.conf.ics)

        self.setup_default_output(conf=self.conf)
        self.not_active_only_modify_matrix = True
        self.not_active_only_modify_matrix = self.conf.options.get(
            'not_active_only_modify_matrix', True,
        )

    def reset(self):
        if hasattr(self.conf, 'options'):
            self.setup_hooks(self.conf.options)

        else:
            self.setup_hooks()

        self.mtx_a = None
        self.solver = None
        self.ts = self.get_default_ts()
        self.clear_equations()

        self._restart_filenames = []

    def setup_hooks(self, options=None):
        """
        Setup various hooks (user-defined functions), as given in `options`.

        Supported hooks:

          - `matrix_hook`

            - check/modify tangent matrix in each nonlinear solver
              iteration

          - `nls_iter_hook`

            - called prior to every iteration of nonlinear solver, if the
              solver supports that
            - takes the Problem instance (`self`) as the first
              argument
        """
        hook_names = ['nls_iter_hook', 'matrix_hook']
        for hook_name in hook_names:
            setattr(self, hook_name, None)
            if options is not None:
                hook = options.get(hook_name, None)
                if hook is not None:
                    hook = self.conf.get_function(hook)
                    setattr(self, hook_name, hook)

    def copy(self, name=None):
        """
        Make a copy of Problem.
        """
        if name is None:
            name = self.name + '_copy'
        obj = self.__class__(name, conf=self.conf, functions=self.functions,
                      domain=self.domain, fields=self.fields,
                      equations=self.equations, auto_conf=False,
                      active_only=self.active_only)

        obj.ebcs = self.ebcs
        obj.epbcs = self.epbcs
        obj.lcbcs = self.lcbcs
        obj.ics = self.ics

        obj.set_conf_solvers(self.conf.solvers, self.conf.options)

        obj.setup_output(output_filename_trunk=self.ofn_trunk,
                         output_dir=self.output_dir,
                         output_format=self.output_format,
                         file_format=self.file_format,
                         split_results_by=self.split_results_by,
                         linearization=self.linearization)

        return obj

    def create_subproblem(self, var_names, known_var_names):
        """
        Create a sub-problem with equations containing only terms with the
        given virtual variables.

        Parameters
        ----------
        var_names : list
            The list of names of virtual variables.
        known_var_names : list
            The list of  names of (already) known state variables.

        Returns
        -------
        subpb : Problem instance
            The sub-problem.
        """
        subpb = Problem(self.name + '_' + '_'.join(var_names), conf=self.conf,
                        functions=self.functions, domain=self.domain,
                        fields=self.fields, auto_conf=False,
                        active_only=self.active_only)
        subpb.set_conf_solvers(self.conf.solvers, self.conf.options)

        subeqs = self.equations.create_subequations(var_names,
                                                    known_var_names)
        subpb.set_equations_instance(subeqs, keep_solvers=True)
        subpb.set_ics(self.ics)

        return subpb

    def setup_default_output(self, conf=None, options=None):
        """
        Provide default values to `Problem.setup_output()`
        from `conf.options` and `options`.
        """
        conf = get_default(conf, self.conf)

        if options and getattr(options, 'output_filename_trunk', None):
            default_output_dir, of = op.split(options.output_filename_trunk)
            default_trunk = io.get_trunk(of)

        else:
            default_trunk = None
            default_output_dir = conf.options.get('output_dir', None)

        if options and getattr(options, 'output_format', None):
            default_output_format = options.output_format

        else:
            default_output_format = conf.options.get('output_format', None)

        default_file_format = conf.options.get('file_format', None)
        default_split_results_by = conf.options.get('split_results_by', None)
        default_float_format = conf.options.get('float_format', None)
        default_linearization = Struct(kind='strip')

        self.setup_output(output_filename_trunk=default_trunk,
                          output_dir=default_output_dir,
                          output_format=default_output_format,
                          file_format=default_file_format,
                          float_format=default_float_format,
                          split_results_by=default_split_results_by,
                          linearization=default_linearization)

    def setup_output(self, output_filename_trunk=None, output_dir=None,
                     output_format=None, file_format=None, float_format=None,
                     split_results_by=None, linearization=None):
        """
        Sets output options to given values, or uses the defaults for
        each argument that is None.
        """
        self.output_modes = {'vtk' : 'sequence',
                             'h5' : 'single', 'h5x' : 'single',
                             'msh' : 'sequence'}

        self.ofn_trunk = get_default(output_filename_trunk,
                                     op.basename(self.domain.name))

        self.set_output_dir(output_dir)

        self.output_format = get_default(output_format, 'vtk')
        self.file_format = file_format
        if self.file_format is not None:
            check_format_suffix(self.file_format, self.output_format)

        self.float_format = get_default(float_format, None)
        self.split_results_by = get_default(split_results_by, None)
        self.linearization = get_default(linearization, Struct(kind='strip'))

        if ((self.output_format == 'h5') and
            (self.linearization.kind == 'adaptive')):
            self.linearization.kind = None

    def set_output_dir(self, output_dir=None):
        """
        Set the directory for output files.

        The directory is created if it does not exist.
        """
        self.output_dir = get_default(output_dir, os.curdir)

        if self.output_dir and not op.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def set_regions(self, conf_regions=None,
                     conf_materials=None, functions=None, allow_empty=False):
        conf_regions = get_default(conf_regions, self.conf.regions)
        functions = get_default(functions, self.functions)

        self.domain.create_regions(conf_regions, functions,
                                   allow_empty=allow_empty)

    def set_materials(self, conf_materials=None):
        """
        Set definition of materials.
        """
        self.conf_materials = get_default(conf_materials, self.conf.materials)

    def select_materials(self, material_names, only_conf=False):
        if type(material_names) == dict:
            conf_materials = transform_materials(material_names)

        else:
            conf_materials = select_by_names(self.conf.materials, material_names)

        if not only_conf:
            self.set_materials(conf_materials)

        return conf_materials

    def set_fields(self, conf_fields=None):
        conf_fields = get_default(conf_fields, self.conf.fields)
        self.fields = fields_from_conf(conf_fields, self.domain.regions)

    def set_variables(self, conf_variables=None):
        """
        Set definition of variables.
        """
        self.conf_variables = get_default(conf_variables, self.conf.variables)
        self.reset()

    def select_variables(self, variable_names, only_conf=False):
        if type(variable_names) == dict:
            conf_variables = transform_variables(variable_names)

        else:
            conf_variables = select_by_names(self.conf.variables, variable_names)

        if not only_conf:
            self.set_variables(conf_variables)

        return conf_variables

    def clear_equations(self):
        self.integrals = None
        self.equations = None
        self.ebcs = None
        self.epbcs = None
        self.lcbcs = None
        self.ics = None

    def set_equations(self, conf_equations=None, user=None,
                      keep_solvers=False, make_virtual=False):
        """
        Set equations of the problem using the `equations` problem
        description entry.

        Fields and Regions have to be already set.
        """
        conf_equations = get_default(conf_equations,
                                     self.conf.get('equations', None))

        self.set_variables(self.conf_variables)
        variables = Variables.from_conf(self.conf_variables, self.fields)

        self.set_materials(self.conf_materials)
        materials = Materials.from_conf(self.conf_materials, self.functions)

        self.integrals = self.get_integrals()

        default_user = vars(self.conf)
        if user is not None:
            default_user.update(user)
        user = default_user
        eterm_options = self.conf.options.get('eterm', {})
        transform = self.conf.options.get('auto_transform_equations', False)
        equations = Equations.from_conf(conf_equations, variables,
                                        self.domain.regions,
                                        materials, self.integrals,
                                        user=user,
                                        eterm_options=eterm_options,
                                        allow_derivatives=transform)

        self.equations = equations
        self.set_ics(self.conf.ics)

        if not keep_solvers:
            self.solver = None

    def set_equations_instance(self, equations, keep_solvers=False):
        """
        Set equations of the problem to `equations`.
        """
        self.mtx_a = None
        self.clear_equations()
        self.equations = equations

        if not keep_solvers:
            self.solver = None

    def get_integrals(self, names=None):
        """
        Get integrals, initialized from problem configuration if available.

        Parameters
        ----------
        names : list, optional
            If given, only the named integrals are returned.

        Returns
        -------
        integrals : Integrals instance
            The requested integrals.
        """
        conf_integrals = self.conf.get('integrals', {})
        integrals = Integrals.from_conf(conf_integrals)

        if names is not None:
            integrals.update([integrals[ii] for ii in names
                              if ii in integrals.names])

        return integrals

    def update_materials(self, ts=None, mode='normal', verbose=True):
        """
        Update materials used in equations.

        Parameters
        ----------
        ts : TimeStepper instance
            The time stepper.
        mode : 'normal', 'update' or 'force'
            The update mode, see :func:`Material.time_update()
            <sfepy.discrete.materials.Material.time_update()>`.
        verbose : bool
            If False, reduce verbosity.
        """
        if self.equations is not None:
            self.update_time_stepper(ts)
            self.equations.time_update_materials(self.ts, mode=mode,
                                                 problem=self, verbose=verbose)


    def update_equations(self, ts=None, ebcs=None, epbcs=None,
                         lcbcs=None, functions=None, create_matrix=False,
                         is_matrix=True, any_dof_conn=None):
        """
        Update equations for current time step.

        The tangent matrix graph is automatically recomputed if the set
        of active essential or periodic boundary conditions changed
        w.r.t. the previous time step.

        Parameters
        ----------
        ts : TimeStepper instance, optional
            The time stepper. If not given, `self.ts` is used.
        ebcs : Conditions instance, optional
            The essential (Dirichlet) boundary conditions. If not given,
            `self.ebcs` are used.
        epbcs : Conditions instance, optional
            The periodic boundary conditions. If not given, `self.epbcs`
            are used.
        lcbcs : Conditions instance, optional
            The linear combination boundary conditions. If not given,
            `self.lcbcs` are used.
        functions : Functions instance, optional
            The user functions for boundary conditions, materials,
            etc. If not given, `self.functions` are used.
        create_matrix : bool
            If True, force the matrix graph computation.
        is_matrix : bool
            If False, the matrix is not created. Has precedence over
            `create_matrix`.
        any_dof_conn : bool or None, default False
            If True, all DOF connectivities are used to pre-allocate the matrix
            graph. If False, only cell region connectivities are used. If None,
            the value is, if available, taken from `conf.options`.
        """
        self.update_time_stepper(ts)
        functions = get_default(functions, self.functions)

        graph_changed = self.equations.time_update(
                                       self.ts,
                                       ebcs, epbcs, lcbcs,
                                       functions, self,
                                       active_only=self.active_only,
                                       verbose=self.conf.get('verbose', True))
        self.graph_changed = graph_changed

        if (is_matrix
            and ((self.active_only and graph_changed)
                 or (self.mtx_a is None) or create_matrix)):
            any_dof_conn = get_default(any_dof_conn,
                                       self.conf.options.get('any_dof_conn',
                                                             False))
            chunk_size = self.conf.options.get('graph_cell_chunk_size', 200000)
            self.mtx_a = self.equations.create_matrix_graph(
                any_dof_conn=any_dof_conn,
                active_only=self.active_only,
                chunk_size=chunk_size,
            )
            ## import sfepy.base.plotutils as plu
            ## plu.spy(self.mtx_a)
            ## plu.plt.show()

    def set_bcs(self, ebcs=None, epbcs=None, lcbcs=None):
        """
        Update boundary conditions.
        """
        if isinstance(ebcs, Conditions):
            self.ebcs = ebcs

        else:
            conf_ebc = get_default(ebcs, self.conf.ebcs)
            self.ebcs = Conditions.from_conf(conf_ebc, self.domain.regions)

            conf_dgebc = self.conf.get("dgebcs", {})
            self.ebcs.extend(Conditions.from_conf(conf_dgebc,
                                                  self.domain.regions))

        if isinstance(epbcs, Conditions):
            self.epbcs = epbcs

        else:
            conf_epbc = get_default(epbcs, self.conf.epbcs)
            self.epbcs = Conditions.from_conf(conf_epbc, self.domain.regions)

            conf_dgepbc = self.conf.get("dgepbcs", {})
            self.ebcs.extend(Conditions.from_conf(conf_dgepbc,
                                                  self.domain.regions))

        if isinstance(lcbcs, Conditions):
            self.lcbcs = lcbcs

        else:
            conf_lcbc = get_default(lcbcs, self.conf.lcbcs)
            self.lcbcs = Conditions.from_conf(conf_lcbc, self.domain.regions)

    def time_update(self, ts=None,
                    ebcs=None, epbcs=None, lcbcs=None,
                    functions=None, create_matrix=False, is_matrix=True,
                    any_dof_conn=None):
        self.set_bcs(get_default(ebcs, self.ebcs),
                     get_default(epbcs, self.epbcs),
                     get_default(lcbcs, self.lcbcs))
        self.update_equations(ts=ts,
                              ebcs=self.ebcs,
                              epbcs=self.epbcs,
                              lcbcs=self.lcbcs,
                              functions=functions,
                              create_matrix=create_matrix,
                              is_matrix=is_matrix,
                              any_dof_conn=any_dof_conn)

    def set_ics(self, ics=None):
        """
        Set the initial conditions to use.
        """
        if isinstance(ics, Conditions):
            self.ics = ics

        else:
            conf_ics = get_default(ics, self.conf.ics)
            self.ics = Conditions.from_conf(conf_ics, self.domain.regions)

    def select_bcs(self, ebc_names=None, epbc_names=None,
                   lcbc_names=None, create_matrix=False, any_dof_conn=None):

        if ebc_names is not None:
            conf_ebc = select_by_names(self.conf.ebcs, ebc_names)
        else:
            conf_ebc = None

        if epbc_names is not None:
            conf_epbc = select_by_names(self.conf.epbcs, epbc_names)
        else:
            conf_epbc = None

        if lcbc_names is not None:
            conf_lcbc = select_by_names(self.conf.lcbcs, lcbc_names)
        else:
            conf_lcbc = None

        self.set_bcs(conf_ebc, conf_epbc, conf_lcbc)
        self.update_equations(ts=self.ts,
                              ebcs=self.ebcs,
                              epbcs=self.epbcs,
                              lcbcs=self.lcbcs,
                              functions=self.functions,
                              create_matrix=create_matrix,
                              any_dof_conn=any_dof_conn)

    def create_state(self):
        return self.get_initial_state()

    def get_mesh_coors(self, actual=False):
        return self.domain.get_mesh_coors(actual=actual)

    def set_mesh_coors(self, coors, update_fields=False, actual=False,
                       clear_all=True, extra_dofs=False):
        """
        Set mesh coordinates.

        Parameters
        ----------
        coors : array
            The new coordinates.
        update_fields : bool
            If True, update also coordinates of fields.
        actual : bool
            If True, update the actual configuration coordinates,
            otherwise the undeformed configuration ones.
        """
        set_mesh_coors(self.domain, self.fields, coors,
                       update_fields=update_fields, actual=actual,
                       clear_all=clear_all, extra_dofs=extra_dofs)

    def refine_uniformly(self, level):
        """
        Refine the mesh uniformly `level`-times.

        Notes
        -----
        This operation resets almost everything (fields, equations, ...)
        - it is roughly equivalent to creating a new Problem
        instance with the refined mesh.
        """
        if level == 0: return

        domain = self.domain
        for ii in range(level):
            domain = domain.refine()

        self.domain = domain
        self.set_regions(self.conf.regions, self.functions)
        self.clear_equations()

        self.set_fields(self.conf.fields)
        self.set_equations(self.conf.equations, user={'ts' : self.ts})

    def get_dim(self, get_sym=False):
        """Returns mesh dimension, symmetric tensor dimension (if `get_sym` is
        True).
        """
        dim = self.domain.mesh.dim
        if get_sym:
            return dim, (dim + 1) * dim // 2
        else:
            return dim

    def init_time(self, ts):
        self.update_time_stepper(ts)
        self.equations.init_time(ts)

        self._restart_filenames = []

    def advance(self, ts=None):
        self.update_time_stepper(ts)
        self.equations.advance(self.ts)

    def get_meshes_from_region(self, reg_names):
        meshes = {}
        for rname in reg_names:
            if rname not in meshes:
                reg = self.domain.regions[rname]
                is_surface = reg.kind != 'cell'
                mesh = Mesh.from_region(reg, self.domain.mesh,
                                        localize=True, is_surface=is_surface)
                meshes[rname] = mesh

        return meshes

    def save_state(self, filename, state=None, out=None,
                   fill_value=None, post_process_hook=None,
                   linearization=None, split_results_by=None, **kwargs):
        """
        Parameters
        ----------
        split_results_by : None, 'region', 'variable'
            If 'region' or 'variable', data of each region/variable are
            stored in a separate file.
            If None, it is set to the application option value.
        linearization : Struct or None
            The linearization configuration for higher order
            approximations. If its kind is 'adaptive', `split_results_by` is
            assumed 'variable'.
        """
        linearization = get_default(linearization, self.linearization)
        if linearization.kind != 'adaptive':
            split_results_by = get_default(split_results_by,
                                           self.split_results_by)

        else:
            split_results_by = 'variable'

        if (out is None) and (state is not None):
            extend = split_results_by not in ['variable', 'region']
            out = state.create_output(fill_value=fill_value,
                                      extend=extend,
                                      linearization=linearization)

            if post_process_hook is not None:
                out = post_process_hook(out, self, state, extend=extend)

        if linearization.kind == 'adaptive':
            for key, val in out.items():
                mesh = val.get('mesh', self.domain.mesh)
                aux = io.edit_filename(filename, suffix='_' + val.var_name)
                mesh.write(aux, io='auto', out={key : val},
                           float_format=self.float_format, **kwargs)
                if hasattr(val, 'levels'):
                    output('max. refinement per group:', val.levels)

        elif split_results_by == 'variable':
            if self.equations is None:
                varnames = {}
                for key, val in out.items():
                    varnames[val.var_name] = 1
                varnames = list(varnames.keys())
                outvars = self.create_variables(varnames)
                itervars = outvars.__iter__
            else:
                itervars = self.equations.variables.iter_state

            rnames = [var.field.region.name for var in itervars()]
            meshes = self.get_meshes_from_region(rnames)

            for var in itervars():
                vout = {}
                for key, val in out.items():
                    try:
                        if val.var_name == var.name:
                            vout[key] = val

                    except AttributeError:
                        msg = 'missing var_name attribute in output!'
                        raise ValueError(msg)

                if len(vout) > 0:
                    mesh = meshes[var.field.region.name]
                    aux = io.edit_filename(filename, suffix='_' + var.name)
                    mesh.write(aux, io='auto', out=vout,
                               float_format=self.float_format, **kwargs)

        elif split_results_by == 'region':
            rnames = [v.region_name for v in out.values()]
            meshes = self.get_meshes_from_region(rnames)

            for rname, mesh in meshes.items():
                rout = {key: val for key, val in out.items()
                        if val.region_name == rname}

                aux = io.edit_filename(filename, suffix='_' + rname)
                mesh.write(aux, io='auto', out=rout,
                           float_format=self.float_format, **kwargs)

        else:
            mesh = out.pop('__mesh__', self.domain.mesh)
            mesh.write(filename, io='auto', out=out,
                       float_format=self.float_format, **kwargs)

    def save_ebc(self, filename, ebcs=None, epbcs=None,
                 force=True, default=0.0):
        """
        Save essential boundary conditions as state variables.

        Parameters
        ----------
        filename : str
            The output file name.
        ebcs : Conditions instance, optional
            The essential (Dirichlet) boundary conditions. If not given,
            `self.conf.ebcs` are used.
        epbcs : Conditions instance, optional
            The periodic boundary conditions. If not given, `self.conf.epbcs`
            are used.
        force : bool
            If True, sequential nonzero values are forced to individual `ebcs`
            so that the conditions are visible even when zero.
        default : float
            The default constant value of state vector.
        """
        output('saving ebc...')
        variables = self.get_variables(auto_create=True)

        if ebcs is None:
            ebcs = Conditions.from_conf(self.conf.ebcs, self.domain.regions)

        if epbcs is None:
            epbcs = Conditions.from_conf(self.conf.epbcs, self.domain.regions)

        try:
            variables.equation_mapping(ebcs, epbcs, self.ts, self.functions,
                                       problem=self)
        except:
            output('cannot make equation mapping!')
            raise

        variables.init_state()
        variables.fill_state(default)

        if force:
            vals = dict_from_keys_init(variables.state)
            for ii, key in enumerate(vals.keys()):
                vals[key] = ii + 1

            variables.apply_ebc(force_values=vals)

        else:
            variables.apply_ebc()

        out = variables.create_output(extend=True)
        self.save_state(filename, out=out, fill_value=default)
        output('...done')

    def save_regions(self, filename_trunk, region_names=None):
        """
        Save regions as meshes.

        Parameters
        ----------
        filename_trunk : str
            The output filename without suffix.
        region_names : list, optional
            If given, only the listed regions are saved.
        """
        filename = '%s.mesh' % filename_trunk
        self.domain.save_regions(filename, region_names=region_names)

    def save_regions_as_groups(self, filename_trunk, region_names=None):
        """
        Save regions in a single mesh but mark them by using different
        element/node group numbers.

        See :func:`Domain.save_regions_as_groups()
        <sfepy.discrete.fem.domain.Domain.save_regions_as_groups()>` for more
        details.

        Parameters
        ----------
        filename_trunk : str
            The output filename without suffix.
        region_names : list, optional
            If given, only the listed regions are saved.
        """
        filename = '%s.%s' % (filename_trunk, self.output_format)
        self.domain.save_regions_as_groups(filename,
                                           region_names=region_names)

    def get_evaluator(self, reuse=False):
        """
        Either create a new Evaluator instance (reuse == False),
        or return an existing instance, created in a preceding call to
        Problem.init_solvers().
        """
        if reuse:
            try:
                ev = self.evaluator
            except AttributeError:
                raise AttributeError('call Problem.init_solvers() or'\
                                     ' set reuse to False!')
        else:
            assemble = self.conf.options.get('assemble_equations_fun', None)
            if assemble is not None:
                assemble = partial(assemble, self)

            UserEvaluator = self.conf.options.get('user_evaluator', None)
            Eval = UserEvaluator if UserEvaluator is not None else Evaluator
            ev = self.evaluator = Eval(self, matrix_hook=self.matrix_hook,
                                       assemble=assemble)

        return ev

    def get_ebc_indices(self):
        """
        Get indices of E(P)BC-constrained DOFs in the full global state vector.
        """
        variables = self.get_variables()

        ebc_indx = []
        epbc_indx = []
        for ii, variable in enumerate(variables.iter_state(ordered=True)):
            eq_map = variable.eq_map
            offset = variables.di.indx[variable.name].start
            ebc_indx.append(eq_map.eq_ebc + offset)
            epbc_indx.append((eq_map.master + offset, eq_map.slave + offset))

        ebc_indx = nm.concatenate(ebc_indx)
        epbc_indx = nm.concatenate(epbc_indx, axis=1)
        return ebc_indx, epbc_indx

    def set_conf_solvers(self, conf_solvers=None, options=None):
        """
        Choose which solvers should be used. If solvers are not set in `options`,
        use the ones named `ls`, `nls`, `ts` and optionally `tsc`. If such
        solver names do not exist, use the first of each required solver kind
        listed in `conf_solvers`.
        """
        conf_solvers = get_default(conf_solvers, self.conf.solvers)
        self.solver_confs = {}

        for key, val in conf_solvers.items():
            self.solver_confs[val.name] = val

        def _find_suitable(prefix):
            cands = []
            for key, val in self.solver_confs.items():
                if val.kind.find(prefix) == 0:
                    if val.name == prefix[:-1]:
                        return val
                    else:
                        cands.append(val)
            if len(cands) > 0:
                return cands[0]
            else:
                return None

        def _get_solver_conf(kind):
            try:
                key = options[kind]
                if key is None:
                    conf = None
                else:
                    conf = self.solver_confs[key]
            except:
                conf = _find_suitable(kind + '.')
            return conf

        self.ts_conf = _get_solver_conf('ts')
        if self.ts_conf is None:
            self.ts_conf = Struct(name='no ts', kind='ts.stationary')
        self.tsc_conf = _get_solver_conf('tsc')
        self.nls_conf = _get_solver_conf('nls')
        self.ls_conf = _get_solver_conf('ls')

        info = 'using solvers:'
        if self.ts_conf:
            info += '\n                ts: %s' % self.ts_conf.name
        if self.tsc_conf:
            info += '\n               tsc: %s' % self.tsc_conf.name
        if self.nls_conf:
            info += '\n               nls: %s' % self.nls_conf.name
        if self.ls_conf:
            info += '\n                ls: %s' % self.ls_conf.name
        if info != 'using solvers:':
            output(info)

    def get_solver_conf(self, name):
        return self.solver_confs[name]

    def init_solvers(self, status=None, ls_conf=None, nls_conf=None,
                     tsc_conf=None, ts_conf=None, force=False):
        """
        Create and initialize solver instances.

        Parameters
        ----------
        status : dict-like, IndexedStruct, optional
            The user-supplied object to hold the time-stepping/nonlinear solver
            convergence statistics.
        ls_conf : Struct, optional
            The linear solver options.
        nls_conf : Struct, optional
            The nonlinear solver options.
        tsc_conf : Struct, optional
            The time step contoller options.
        ts_conf : Struct, optional
            The time stepping solver options.
        force : bool
            If True, re-create the solver instances even if they already exist
            in `self.nls` attribute.
        """
        if (self.solver is None) or force:
            ls_conf = get_default(ls_conf, self.ls_conf,
                                  'you must set linear solver!')
            nls_conf = get_default(nls_conf, self.nls_conf,
                                   'you must set nonlinear solver!')

            fb_list = []
            for ii in range(100):
                fb_list.append((ls_conf.kind, ls_conf))
                if hasattr(ls_conf, 'fallback'):
                    ls_conf = self.solver_confs[ls_conf.fallback]
                else:
                    break

            if len(fb_list) > 1:
                ls = use_first_available(fb_list, context=self)
            else:
                ls = Solver.any_from_conf(ls_conf, context=self)

            ev = self.get_evaluator()

            if self.conf.options.get('ulf', False):
                self.nls_iter_hook = ev.new_ulf_iteration

            if status is None:
                status = IndexedStruct()

            status.set_default('nls_status', IndexedStruct())

            nls = Solver.any_from_conf(nls_conf, fun=ev.eval_residual,
                                       fun_grad=ev.eval_tangent_matrix,
                                       lin_solver=ls,
                                       iter_hook=self.nls_iter_hook,
                                       status=status.nls_status, context=self)

            ts_conf = get_default(ts_conf, self.ts_conf)
            if ts_conf is None:
                self.set_solver(nls, status=status)

            else:
                tsc_conf = get_default(tsc_conf, self.tsc_conf)
                if tsc_conf is not None:
                    tsc = Solver.any_from_conf(tsc_conf)

                else:
                    tsc = None

                tss = Solver.any_from_conf(ts_conf, nls=nls, tsc=tsc,
                                           context=self, status=status)
                self.set_solver(tss)

    def get_default_ts(self, t0=None, t1=None, dt=None, n_step=None,
                       step=None):
        t0 = get_default(t0, 0.0)
        t1 = get_default(t1, 1.0)
        dt = get_default(dt, 1.0)
        n_step = get_default(n_step, 1)

        ts = TimeStepper(t0, t1, dt, n_step, step=step)

        return ts

    def update_time_stepper(self, ts):
        if ts is not None:
            self.ts = ts

    def get_timestepper(self):
        return self.ts

    def set_solver(self, solver, status=None):
        """
        Set a time-stepping or nonlinear solver to be used in
        :func:`Problem.solve()` call.

        Parameters
        ----------
        solver : NonlinearSolver or TimeSteppingSolver instance
            The nonlinear or time-stepping solver.
        status : dict-like, optional
            The user-supplied object to hold the solver convergence statistics.

        Notes
        -----
        A copy of the solver is used, and the nonlinear solver functions are
        set to those returned by :func:`Problem.get_nls_functions()`, if not
        set already. If a nonlinear solver is set, a default StationarySolver
        instance is created automatically as the time-stepping solver. Also
        sets `self.ts` attribute.

        If `self.conf.options.auto_transform_equations` is True (the default is
        False), the problem equations are automatically transformed to a form
        suitable for the given solver. Implemented for
        :class:`ElastodynamicsBaseTS
        <sfepy.solvers.ts_solvers.ElastodynamicsBaseTS>`-based solvers. If it
        is False, `solver.var_names` have to be defined.
        """
        transform = self.conf.options.get('auto_transform_equations', False)
        if isinstance(solver, ElastodynamicsBaseTS):
            self.solver = solver.copy()
            if transform:
                if self.get('orig_equations', None) is None:
                    equations, var_names = transform_equations_ed(
                        self.equations, self.get_materials(),
                    )
                    self.orig_equations = self.equations
                    self.ed_var_names = var_names
                    self.equations = equations

                    output('transformed equations of elastodynamics solvers:')
                    self.equations.print_terms()

                else:
                    var_names = self.ed_var_names

                self.solver.conf.var_names = var_names

            elif solver.conf.get('var_names', None) is None:
                raise ValueError('specify names of dynamic variables by'
                                 ' defining `solver.conf.var_names`!')

        elif isinstance(solver, NonlinearSolver):
            self.solver = StationarySolver({}, nls=solver.copy(),
                                           ts=self.get_default_ts(),
                                           status=status)

        else:
            self.solver = solver.copy()

        self.ts = self.solver.ts
        self.status = get_default(solver.status, IndexedStruct())

        # Assign the nonlinear solver functions.
        nls = self.get_nls()
        if nls.fun is None:
            fun, fun_grag, iter_hook = self.get_nls_functions()
            nls.fun = fun
            nls.fun_grad = fun_grag
            nls.iter_hook = iter_hook

    def try_presolve(self, mtx):
        ls = self.get_ls()

        timer = Timer(start=True)
        ls.presolve(mtx)
        tt = timer.stop()
        output('presolve: %.2f [s]' % tt)

    def get_solver(self):
        return self.get_tss()

    def get_tss(self):
        tss = get_default(None, self.solver, 'solver is not set!')
        return tss

    def get_tss_functions(self, update_bcs=True, update_materials=True,
                          save_results=True,
                          step_hook=None, post_process_hook=None):
        """
        Get the problem-dependent functions required by the time-stepping
        solver during the solution process.

        Parameters
        ----------
        update_bcs : bool, optional
            If True, update the boundary conditions in each `prestep_fun` call.
        update_materials : bool, optional
            If True, update the values of material parameters in each
            `prestep_fun` call.
        save_results : bool, optional
            If True, save the results in each `poststep_fun` call.
        step_hook : callable, optional
            The optional user-defined function that is called in each
            `poststep_fun` call before saving the results.
        post_process_hook : callable, optional
            The optional user-defined function that is passed in each
            `poststep_fun` to :func:`Problem.save_state()`.

        Returns
        -------
        init_fun : callable
            The initialization function called before the actual time-stepping.
        prestep_fun : callable
            The function called in each time (sub-)step prior to the nonlinear
            solver call.
        poststep_fun : callable
            The function called at the end of each time step.
        """
        is_save = make_is_save(self.conf.options)

        def init_fun(ts, vec0):
            if not ts.is_quasistatic:
                self.init_time(ts)

            is_save.reset(ts)

            restart_filename = self.conf.options.get('load_restart', None)
            if restart_filename is not None:
                variables = self.load_restart(restart_filename, ts=ts)
                self.advance(ts)
                ts.advance()
                vec0 = variables.get_state(self.active_only)

            return vec0

        def prestep_fun(ts, vec):
            if update_bcs:
                self.time_update(ts)
                variables = self.equations.variables
                # EBCs are applied in Evaluator.eval_residual().
                variables.set_state(vec, self.active_only, apply_ebc=False)

            if update_materials:
                self.update_materials(ts, verbose=self.conf.get('verbose', True))

            return vec

        def poststep_fun(ts, vec):
            variables = self.equations.variables
            # EBCs need to be applied here because of algebraically computed
            # variables in elastodynamics solvers.
            variables.set_state(vec, self.active_only, apply_ebc=True)
            if step_hook is not None:
                step_hook(self, ts, variables)

            restart_filename = self.get_restart_filename(ts=ts)
            if restart_filename is not None:
                self.save_restart(restart_filename, ts=ts)

            if save_results and is_save(ts):
                if not isinstance(self.get_solver(), StationarySolver):
                    suffix = ts.suffix % ts.step

                else:
                    suffix = None

                filename = self.get_output_name(suffix=suffix)
                self.save_state(filename, variables,
                                post_process_hook=post_process_hook,
                                split_results_by=None,
                                ts=ts,
                                file_format=self.file_format)

            self.advance(ts)
            return vec

        return init_fun, prestep_fun, poststep_fun

    def get_nls_functions(self):
        """
        Returns functions to be used by a nonlinear solver to evaluate the
        nonlinear function value (the residual) and its gradient (the tangent
        matrix) corresponding to the problem equations.

        Returns
        -------
        fun : function
            The function ``fun(x)`` for computing the residual.
        fun_grad : function
            The function ``fun_grad(x)`` for computing the tangent matrix.
        iter_hook : function
            The optional (user-defined) function to be called before each
            nonlinear solver iteration iteration.
        """
        ev = self.get_evaluator()
        return ev.eval_residual, ev.eval_tangent_matrix, self.nls_iter_hook

    def get_nls(self):
        tss = self.get_tss()
        return tss.nls

    def get_ls(self):
        nls = self.get_nls()
        return nls.lin_solver

    def is_linear(self):
        nls = self.get_nls()
        return nls.conf.get('is_linear', False)

    def set_linear(self, is_linear):
        nls = self.get_nls()
        nls.conf.is_linear = is_linear

    def set_default_state(self, vec=None):
        """
        Return variables with an initialized state.

        A convenience function that obtains the problem equations' variables,
        initializes the state ones with zeros (default) or using `vec` and then
        returns the variables.
        """
        variables = self.get_variables()
        variables.init_state(vec=vec)
        return variables

    def get_initial_state(self, vec=None):
        """
        Create a zero state and apply initial conditions.
        """
        self.equations.setup_initial_conditions(self.ics, self.functions)
        self.equations.init_state(vec=vec)
        self.equations.apply_ic()

        return self.equations.variables

    def solve(self, state0=None, status=None, force_values=None,
              var_data=None, update_bcs=True, update_materials=True,
              save_results=True,
              step_hook=None, post_process_hook=None,
              post_process_hook_final=None,
              report_nls_status=False, log_nls_status=False,
              verbose=True):
        """
        Solve the problem equations by calling the top-level solver.

        Before calling this function the top-level solver has to be set, see
        :func:`Problem.set_solver()`. Also, the boundary conditions and the
        initial conditions (for time-dependent problems) has to be set, see
        :func:`Problem.set_bcs()`, :func:`Problem.set_ics()`.

        Parameters
        ----------
        state0 : array, optional
            If given, the initial state - then the initial conditions stored in
            the Problem instance are ignored. By default, the initial state is
            created and the initial conditions are applied automatically.
        status : dict-like, optional
            The user-supplied object to hold the solver convergence statistics.
        force_values : dict of floats or float, optional
            If given, the supplied values override the values of the essential
            boundary conditions.
        var_data : dict, optional
            A dictionary of {variable_name : data vector} used to initialize
            parameter variables.
        update_bcs : bool, optional
            If True, update the boundary conditions in each `prestep_fun` call.
            See :func:`Problem.get_tss_functions()`.
        update_materials : bool, optional
            If True, update the values of material parameters in each
            `prestep_fun` call. See :func:`Problem.get_tss_functions()`.
        save_results : bool, optional
            If True, save the results in each `poststep_fun` call. See
            :func:`Problem.get_tss_functions()`.
        step_hook : callable, optional
            The optional user-defined function that is called in each
            `poststep_fun` call before saving the results. See
            :func:`Problem.get_tss_functions()`.
        post_process_hook : callable, optional
            The optional user-defined function that is passed in each
            `poststep_fun` to :func:`Problem.save_state()`. See
            :func:`Problem.get_tss_functions()`.
        post_process_hook_final : callable, optional
            The optional user-defined function that is called after the
            top-level solver returns.
        report_nls_status: bool, optional
            If True, print summary non-linear solver info.
        log_nls_status: bool, optional
            If True, log non-linear solver info into a CSV file.

        Returns
        -------
        variables : Variables
            The variables with the final time step state.
        """
        if status is None:
            status = IndexedStruct()

        if self.solver is None:
            self.init_solvers(status=status)

        tss = self.get_solver()

        self.equations.set_data(var_data, ignore_unknown=True)

        report_nls_status = getattr(
            self.conf.options, 'report_nls_status', report_nls_status)
        log_nls_status = getattr(
            self.conf.options, 'log_nls_status', log_nls_status)

        if self.conf.options.get('block_solve', False):
            variables = self.block_solve(state0, status=status,
                                         save_results=save_results,
                                         step_hook=step_hook,
                                         post_process_hook=post_process_hook,
                                         report_nls_status=report_nls_status,
                                         log_nls_status=log_nls_status,
                                         verbose=verbose)

        else:
            if state0 is not None:
                variables = self.set_default_state(vec=state0)

            else:
                variables = self.get_initial_state()

            self.time_update(tss.ts)

            variables.apply_ebc(force_values=force_values)

            if self.is_linear():
                mtx = prepare_matrix(self, variables) # Updates materials.
                self.try_presolve(mtx)

            init_fun, prestep_fun, poststep_fun = self.get_tss_functions(
                update_bcs=update_bcs, update_materials=update_materials,
                save_results=save_results,
                step_hook=step_hook, post_process_hook=post_process_hook)

            tss.set_dof_info(variables.adi)
            vec = tss(variables.get_state(self.active_only, force=True),
                      init_fun=init_fun,
                      prestep_fun=prestep_fun,
                      poststep_fun=poststep_fun,
                      status=status,
                      log_nls_status=log_nls_status)

            time_stats = status.get('time_stats')
            if time_stats is not None:
                output('====== time stats ======')
                for key in time_stats.keys():
                    output('%12s: %.8f [s]' % ('nls ' + key, time_stats[key]))

            if report_nls_status:
                step_stats = status.get('step_stats')
                if step_stats is not None:
                    output('====== step stats ======')
                    output.prefix, aux = '', output.prefix
                    # output('  step, time,      cond,  nit, ls_nit,      err0,'
                    output('  step,         time, cond,  nit, ls_nit,      err0,'
                           '       err, elapsed [s]')  # 78 chars long
                    if len(step_stats) > 1 and step_stats[0].get("step") > 0:
                        s0 = IndexedStruct(step=0, step_time=0.0, condition=0,
                                           n_iter=0, ls_n_iter=0,
                                           err=0, err0=0, time=0)
                        step_stats = [s0] + step_stats

                    for step in step_stats:
                        msg = f'{step.get("step") + 1:6}, '
                        msg += f'{step.get("step_time"):.6e}, '
                        msg += f'{step.get("condition"):4}, '
                        msg += f'{step.get("n_iter"):4}, '
                        msg += f'{step.get("ls_n_iter"):6}, '
                        msg += f'{step.get("err0"):.3e}, '
                        msg += f'{step.get("err"):.3e}, '
                        msg += f'{step.get("time"):.4f}'

                        output(msg)
                    output.prefix, aux = aux, output.prefix

            output('solved in %d steps in %.2f seconds'
                   % (status['n_step'], status['time']), verbose=verbose)

            variables.set_state(vec, self.active_only)

        if post_process_hook_final is not None: # User postprocessing.
            post_process_hook_final(self, variables)

        return variables

    def block_solve(self, state0=None, status=None, save_results=True,
                    step_hook=None, post_process_hook=None,
                    report_nls_status=False, log_nls_status=False,
                    verbose=True):
        """
        Call :func:`Problem.solve()` sequentially for the individual matrix
        blocks of a block-triangular matrix. It is called by
        :func:`Problem.solve()` if the `'block_solve'` option is set to True.
        """
        from sfepy.base.base import invert_dict, get_subdict
        from sfepy.base.resolve_deps import resolve

        if not isinstance(self.get_solver(), StationarySolver):
            msg = 'The block solve can be used only for stationary problems!'
            raise ValueError(msg)

        def replace_virtuals(deps, pairs):
            out = {}
            for key, val in deps.items():
                out[pairs[key]] = val

            return out

        if state0 is not None:
            variables = self.set_default_state(vec=state0)

        else:
            variables = self.get_initial_state()

        vtos = variables.get_dual_names()
        vdeps = self.equations.get_variable_dependencies()
        sdeps = replace_virtuals(vdeps, vtos)

        sorder = resolve(sdeps)

        stov = invert_dict(vtos)
        vorder = [[stov[ii] for ii in block] for block in sorder]

        parts = variables.get_state_parts()
        solved = []
        for ib, block in enumerate(vorder):
            output('solving for %s...' % sorder[ib], verbose=verbose)

            subpb = self.create_subproblem(block, solved)
            subpb.conf.options.block_solve = False

            subpb.equations.print_terms()

            subvars0 = subpb.get_initial_state()

            vals = get_subdict(parts, block)
            subvars0.set_state_parts(vals)

            subvars = subpb.solve(state0=subvars0(), status=status,
                                  save_results=False, step_hook=step_hook,
                                  post_process_hook=post_process_hook,
                                  report_nls_status=report_nls_status,
                                  log_nls_status=log_nls_status,
                                  verbose=verbose)

            variables.set_state_parts(subvars.get_state_parts())

            solved.extend(sorder[ib])
            output('...done', verbose=verbose)

        if step_hook is not None:
            step_hook(self, None, variables)

        if save_results:
            self.save_state(self.get_output_name(), variables,
                            post_process_hook=post_process_hook,
                            split_results_by=None)

        return variables

    def create_evaluable(self, expression, try_equations=True, auto_init=False,
                         preserve_caches=False, copy_materials=True,
                         integrals=None,
                         ebcs=None, epbcs=None, lcbcs=None,
                         ts=None, functions=None,
                         mode='eval', var_dict=None, strip_variables=True,
                         extra_args=None, active_only=True,
                         eterm_options=None, verbose=True, **kwargs):
        """
        Create evaluable object (equations and corresponding variables)
        from the `expression` string. Convenience function calling
        :func:`create_evaluable()
        <sfepy.discrete.evaluate.create_evaluable()>` with defaults provided
        by the Problem instance `self`.

        The evaluable can be repeatedly evaluated by calling
        :func:`eval_equations() <sfepy.discrete.evaluate.eval_equations()>`,
        e.g. for different values of variables.

        Parameters
        ----------
        expression : str
            The expression to evaluate.
        try_equations : bool
            Try to get variables from `self.equations`. If this fails,
            variables can either be provided in `var_dict`, as keyword
            arguments, or are created automatically according to the
            expression.
        auto_init : bool
            Set values of all variables to all zeros.
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        copy_materials : bool
            Work with a copy of `self.equations.materials` instead of
            reusing them. Safe but can be slow.
        integrals : Integrals instance, optional
            The integrals to be used. Automatically created as needed if
            not given.
        ebcs : Conditions instance, optional
            The essential (Dirichlet) boundary conditions for 'weak'
            mode. If not given, `self.ebcs` are used.
        epbcs : Conditions instance, optional
            The periodic boundary conditions for 'weak'
            mode. If not given, `self.epbcs` are used.
        lcbcs : Conditions instance, optional
            The linear combination boundary conditions for 'weak'
            mode. If not given, `self.lcbcs` are used.
        ts : TimeStepper instance, optional
            The time stepper. If not given, `self.ts` is used.
        functions : Functions instance, optional
            The user functions for boundary conditions, materials
            etc. If not given, `self.functions` are used.
        mode : one of 'eval', 'el_avg', 'qp', 'weak'
            The evaluation mode - 'weak' means the finite element
            assembling, 'qp' requests the values in quadrature points,
            'el_avg' element averages and 'eval' means integration over
            each term region.
        var_dict : dict, optional
            The variables (dictionary of (variable name) : (Variable instance))
            to be used in the expression. Use this if the name of a variable
            conflicts with one of the parameters of this method.
        strip_variables : bool
            If False, the variables in `var_dict` or `kwargs` not present in
            the expression are added to the actual variables as a context.
        extra_args : dict, optional
            Extra arguments to be passed to terms in the expression.
        active_only : bool
            If True, in 'weak' mode, the (tangent) matrices and residual
            vectors (right-hand sides) contain only active DOFs.
        eterm_options : dict, optional
            The einsum-based terms evaluation options.
        verbose : bool
            If False, reduce verbosity.
        **kwargs : keyword arguments
            Additional variables can be passed as keyword arguments, see
            `var_dict`.

        Returns
        -------
        equations : Equations instance
            The equations that can be evaluated.
        variables : Variables instance
            The corresponding variables. Set their values and use
            :func:`eval_equations() <sfepy.discrete.evaluate.eval_equations()>`.

        Examples
        --------
        `problem` is Problem instance.

        >>> out = problem.create_evaluable('ev_integrate.i1.Omega(u)')
        >>> equations, variables = out

        `vec` is a vector of coefficients compatible with the field
        of 'u' - let's use all ones.

        >>> vec = nm.ones((variables['u'].n_dof,), dtype=nm.float64)
        >>> variables['u'].set_data(vec)
        >>> vec_qp = eval_equations(equations, variables, mode='qp')

        Try another vector:

        >>> vec = 3 * nm.ones((variables['u'].n_dof,), dtype=nm.float64)
        >>> variables['u'].set_data(vec)
        >>> vec_qp = eval_equations(equations, variables, mode='qp')
        """
        from sfepy.discrete.equations import get_expression_arg_names

        variables = Variables(get_default(var_dict, {}).values())
        var_context = get_default(var_dict, {})

        if try_equations and self.equations is not None:
            # Make a copy, so that possible variable caches are preserved.
            for key, var in self.equations.variables.as_dict().items():
                if key in variables:
                    continue
                var = var.copy(name=key)
                if not preserve_caches:
                    var.clear_evaluate_cache()
                variables[key] = var

        elif var_dict is None:
            possible_var_names = get_expression_arg_names(expression)
            variables = self.create_variables(possible_var_names)

        materials = self.get_materials()
        if copy_materials or (materials is None):
            possible_mat_names = get_expression_arg_names(expression)
            materials = self.create_materials(possible_mat_names)

        else:
            materials = Materials(objs=materials._objs)

        _kwargs = copy(kwargs)
        for key, val in kwargs.items():
            if isinstance(val, Variable):
                if val.name != key:
                    msg = 'inconsistent variable name! (%s == %s)' \
                          % (val.name, key)
                    raise ValueError(msg)
                var_context[key] = variables[key] = val.copy(name=key)
                _kwargs.pop(key)

            elif isinstance(val, Material):
                if val.name != key:
                    msg = 'inconsistent material name! (%s == %s)' \
                          % (val.name, key)
                    raise ValueError(msg)
                materials[val.name] = val
                _kwargs.pop(key)

        kwargs = _kwargs

        ebcs = get_default(ebcs, self.ebcs)
        epbcs = get_default(epbcs, self.epbcs)
        lcbcs = get_default(lcbcs, self.lcbcs)
        ts = get_default(ts, self.get_timestepper())
        functions = get_default(functions, self.functions)
        integrals = get_default(integrals, self.get_integrals())

        out = create_evaluable(expression, self.fields, materials,
                               variables, integrals,
                               ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs,
                               ts=ts, functions=functions,
                               auto_init=auto_init,
                               mode=mode, extra_args=extra_args,
                               active_only=active_only,
                               eterm_options=eterm_options,
                               verbose=verbose,
                               kwargs=kwargs)

        if not strip_variables:
            variables = out[1]
            variables.extend([var for var in var_context.values()
                              if var not in variables])

        equations = out[0]
        mode = 'update' if not copy_materials else 'normal'
        equations.time_update_materials(self.ts, mode=mode, problem=self,
                                        verbose=verbose)

        return out

    def evaluate(self, expression, try_equations=True, auto_init=False,
                 preserve_caches=False, copy_materials=True, integrals=None,
                 ebcs=None, epbcs=None, lcbcs=None, ts=None, functions=None,
                 mode='eval', dw_mode='vector', term_mode=None,
                 var_dict=None, strip_variables=True, ret_variables=False,
                 any_dof_conn=False,
                 active_only=True, eterm_options=None, verbose=True,
                 extra_args=None, **kwargs):
        """
        Evaluate an expression, convenience wrapper of
        :func:`Problem.create_evaluable` and
        :func:`eval_equations() <sfepy.discrete.evaluate.eval_equations>`.

        Parameters
        ----------
        dw_mode : 'vector' or 'matrix'
            The assembling mode for 'weak' evaluation mode.
        term_mode : str
            The term call mode - some terms support different call modes
            and depending on the call mode different values are
            returned.
        ret_variables : bool
            If True, return the variables that were created to evaluate
            the expression.
        other : arguments
            See docstrings of :func:`Problem.create_evaluable()`.

        Returns
        -------
        out : array
            The result of the evaluation.
        variables : Variables instance
            The variables that were created to evaluate
            the expression. Only provided if `ret_variables` is True.
        """
        aux = self.create_evaluable(expression,
                                    try_equations=try_equations,
                                    auto_init=auto_init,
                                    preserve_caches=preserve_caches,
                                    copy_materials=copy_materials,
                                    integrals=integrals,
                                    ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs,
                                    ts=ts, functions=functions,
                                    mode=mode, var_dict=var_dict,
                                    strip_variables=strip_variables,
                                    extra_args=extra_args,
                                    active_only=active_only,
                                    eterm_options=eterm_options,
                                    verbose=verbose,
                                    **kwargs)
        equations, variables = aux

        out = eval_equations(equations, variables,
                             preserve_caches=preserve_caches,
                             mode=mode, dw_mode=dw_mode, term_mode=term_mode,
                             active_only=active_only,
                             any_dof_conn=any_dof_conn,
                             verbose=verbose)

        if ret_variables:
            out = (out, variables)

        return out

    def eval_equations(self, names=None, preserve_caches=False,
                   mode='eval', dw_mode='vector', term_mode=None,
                   active_only=True, any_dof_conn=False, verbose=True):
        """
        Evaluate (some of) the problem's equations, convenience wrapper of
        :func:`eval_equations() <sfepy.discrete.evaluate.eval_equations>`.

        Parameters
        ----------
        names : str or sequence of str, optional
            Evaluate only equations of the given name(s).
        preserve_caches : bool
            If True, do not invalidate evaluate caches of variables.
        mode : one of 'eval', 'el_avg', 'qp', 'weak'
            The evaluation mode - 'weak' means the finite element
            assembling, 'qp' requests the values in quadrature points,
            'el_avg' element averages and 'eval' means integration over
            each term region.
        dw_mode : 'vector' or 'matrix'
            The assembling mode for 'weak' evaluation mode.
        term_mode : str
            The term call mode - some terms support different call modes
            and depending on the call mode different values are
            returned.
        verbose : bool
            If False, reduce verbosity.

        Returns
        -------
        out : dict or result
            The evaluation result. In 'weak' mode it is the vector or sparse
            matrix, depending on `dw_mode`. Otherwise, it is a dict of results
            with equation names as keys or a single result for a single
            equation.
        """
        return eval_equations(self.equations, self.equations.variables,
                              names=names, preserve_caches=preserve_caches,
                              mode=mode, dw_mode=dw_mode, term_mode=term_mode,
                              active_only=active_only,
                              any_dof_conn=any_dof_conn, verbose=verbose)

    def get_materials(self):
        if self.equations is not None:
            materials = self.equations.materials

        else:
            materials = None

        return materials

    def create_materials(self, mat_names=None):
        """
        Create materials with names in `mat_names`. Their definitions
        have to be present in `self.conf.materials`.

        Notes
        -----
        This method does not change `self.equations`, so it should not
        have any side effects.
        """
        if mat_names is not None:
            conf_materials = self.select_materials(mat_names, only_conf=True)

        else:
            conf_materials = self.conf.materials

        materials = Materials.from_conf(conf_materials, self.functions)

        return materials

    def get_variables(self, auto_create=False):
        if self.equations is not None:
            variables = self.equations.variables

        elif auto_create:
            variables = self.create_variables()

        else:
            variables = None

        return variables

    def create_variables(self, var_names=None):
        """
        Create variables with names in `var_names`. Their definitions
        have to be present in `self.conf.variables`.

        Notes
        -----
        This method does not change `self.equations`, so it should not
        have any side effects.
        """
        if var_names is not None:
            conf_variables = self.select_variables(var_names, only_conf=True)

        else:
            conf_variables = self.conf.variables

        variables = Variables.from_conf(conf_variables, self.fields)

        return variables

    def get_output_name(self, suffix=None, extra=None, mode=None):
        """
        Return default output file name, based on the output directory,
        output format, step suffix and mode. If present, the extra
        string is put just before the output format suffix.
        """
        out = op.join(self.output_dir, self.ofn_trunk)

        if suffix is not None:
            if mode is None:
                mode = self.output_modes[self.output_format]

            if mode == 'sequence':
                out = '.'.join((out, suffix))

        if extra is not None:
            out = '.'.join((out, extra, self.output_format))
        else:
            out = '.'.join((out, self.output_format))

        return out

    def remove_bcs(self):
        """
        Convenience function to remove boundary conditions.
        """
        self.time_update(ebcs={}, epbcs={}, lcbcs={})

    def get_restart_filename(self, ts=None):
        """
        If restarts are allowed in problem definition options, return the
        restart file name, based on the output directory and time step.
        """
        if self.conf.options.get('save_restart', None) is None:
            return

        suffix = 'restart'
        if ts is not None:
            suffix += '-' + ts.suffix % ts.step

        aux = self.get_output_name(extra=suffix)
        iext = len(aux) - len('.' + self.output_format)
        restart_filename = aux[:iext] + '.h5'

        return restart_filename

    def save_restart(self, filename, ts=None):
        """
        Save the current state and time step to a restart file.

        Parameters
        ----------
        filename : str
            The restart file name.
        ts : TimeStepper instance, optional
            The time stepper. If not given, a default one is created.

        Notes
        -----
        Does not support terms with internal state.
        """
        import tables as pt

        if ts is None:
            ts = self.get_default_ts()

        fd = pt.open_file(filename, mode='w', title='SfePy restart file')

        tgroup = fd.create_group('/', 'ts', 'ts')
        for key, val in ts.get_state().items():
            fd.create_array(tgroup, key, val, key)

        variables = self.get_variables()
        if variables.r_vec is not None:
            fd.create_array('/', 'r_vec', variables.r_vec,
                            'reduced state vector')

        for var in variables.iter_state():
            vgroup = fd.create_group('/', var.name, var.name)

            history_length = len(var.data)
            fd.create_array(vgroup, 'history_length', history_length,
                            'history length')
            for ii in range(history_length):
                data = var(step=-ii)
                fd.create_array(vgroup, 'data_%d' % ii, data, 'data')

        fd.close()

        mode = self.conf.options.get('save_restart', None)

        if (mode == -1) and len(self._restart_filenames):
            last_filename = self._restart_filenames.pop()
            if last_filename != filename:
                try:
                    os.remove(last_filename)

                except OSError:
                    pass

        self._restart_filenames.append(filename)

    def load_restart(self, filename, ts=None):
        """
        Load the current state and time step from a restart file.

        Alternatively, a regular output file in the HDF5 format can be used in
        place of the restart file. In that case the restart is only
        approximate, because higher order field DOFs (if any) were stripped
        out. Files with the adaptive linearization are not supported. Use with
        caution!

        Parameters
        ----------
        filename : str
            The restart file name.
        ts : TimeStepper instance, optional
            The time stepper. If not given, a default one is created.
            Otherwise, it is modified in place.

        Returns
        -------
        variables : Variables instance
            The loaded variables.
        """
        import tables as pt

        if ts is None:
            ts = self.get_default_ts()

        variables = self.get_variables()

        output('loading restart file "%s"...' % filename)

        fd = pt.open_file(filename, mode='r')

        if fd.title == 'SfePy restart file':
            ts_state = {}
            for val in fd.root.ts._f_walknodes():
                ts_state[val.name] = val.read()

            ts.set_state(**ts_state)

            vec = variables.create_vec()
            for var in variables.iter_state():
                vgroup = fd.root._f_get_child(var.name)

                history_length = vgroup.history_length.read()
                for ii in range(0, history_length):
                    data = vgroup._f_get_child('data_%d' % ii).read()
                    var.locked = False
                    var.set_data(data, step=-ii)
                    var.locked = True
                    if ii == 0:
                        variables.set_vec_part(vec, var.name, data)

            variables.init_state(vec=vec)

            if '/r_vec' in fd:
                r_vec = fd.root.r_vec.read()
                variables.r_vec = r_vec

            fd.close()

        elif fd.title == 'SfePy output file':
            from sfepy.discrete.fem.meshio import MeshIO

            output('WARNING: using a SfePy output file in place of a restart'
                   ' file discards higher order DOFs! Use with caution!')

            fd.close()
            io = MeshIO.any_from_filename(filename)

            out = io.read_data(step=ts.step)

            vec = variables.create_vec()
            for var in variables.iter_state():
                val = out[var.name]
                var.set_from_mesh_vertices(val.data)
                variables.set_vec_part(vec, var.name, var())

            variables.init_state(vec=vec)

        else:
            raise IOError('unknown file type! ("%s" in ("%s", "%s"))'
                          % (fd.title,
                             'SfePy restart file', 'SfePy output file'))

        output('...done')

        return variables
