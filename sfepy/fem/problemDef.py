import os
import os.path as op
import time
from copy import copy

import numpy as nm

from sfepy.base.base import dict_from_keys_init, select_by_names
from sfepy.base.base import output, get_default, get_default_attr, Struct
import sfepy.base.ioutils as io
from sfepy.base.conf import ProblemConf, get_standard_keywords
from sfepy.base.conf import transform_variables, transform_materials
from functions import Functions
from mesh import Mesh
from domain import Domain
from fields_base import fields_from_conf
from variables import Variables, Variable
from materials import Materials, Material
from equations import Equations
from integrals import Integrals
from sfepy.fem.state import State
from sfepy.fem.conditions import Conditions
from sfepy.fem.evaluate import create_evaluable, eval_equations
import fea as fea
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.evaluate import BasicEvaluator, LCBCEvaluator
from sfepy.solvers import Solver
from sfepy.solvers.ls import ScipyDirect
from sfepy.solvers.nls import Newton

##
# 29.01.2006, c
class ProblemDefinition( Struct ):
    """
    Problem definition, the top-level class holding all data necessary to solve
    a problem.

    It can be constructed from a :class:`ProblemConf` instance using
    `ProblemDefinition.from_conf()` or directly from a problem
    description file using `ProblemDefinition.from_conf_file()`

    For interactive use, the constructor requires only the `equations`,
    `nls` and `ls` keyword arguments.
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

        obj = ProblemDefinition.from_conf(conf,
                                          init_fields=init_fields,
                                          init_equations=init_equations,
                                          init_solvers=init_solvers)
        return obj

    @staticmethod
    def from_conf(conf, init_fields=True, init_equations=True,
                  init_solvers=True):
        if conf.options.get_default_attr('absolute_mesh_path', False):
            conf_dir = None
        else:
            conf_dir = op.dirname(conf.funmod.__file__)

        functions = Functions.from_conf(conf.functions)

        mesh = Mesh.from_file(conf.filename_mesh, prefix_dir=conf_dir)

        trans_mtx = conf.options.get_default_attr('mesh_coors_transform', None)

        if trans_mtx is not None:
            mesh.transform_coors(trans_mtx)

        domain = Domain(mesh.name, mesh)
        if get_default_attr(conf.options, 'ulf', False):
            domain.mesh.coors_act = domain.mesh.coors.copy()

        obj = ProblemDefinition('problem_from_conf', conf=conf,
                                functions=functions, domain=domain,
                                auto_conf=False, auto_solvers=False)

        obj.set_regions(conf.regions, obj.functions)

        obj.clear_equations()

        if init_fields:
            obj.set_fields( conf.fields )

            if init_equations:
                obj.set_equations(conf.equations, user={'ts' : obj.ts})

        if init_solvers:
            obj.set_solvers( conf.solvers, conf.options )

        return obj

    def __init__(self, name, conf=None, functions=None,
                 domain=None, fields=None, equations=None, auto_conf=True,
                 nls=None, ls=None, ts=None, auto_solvers=True):
        self.name = name
        self.conf = conf
        self.functions = functions

        self.reset()

        self.ts = get_default(ts, self.get_default_ts())

        if auto_conf:
            if equations is None:
                raise ValueError('missing equations in auto_conf mode!')

            self.equations = equations

            if fields is None:
                variables = self.equations.variables
                fields = {}
                for field in [var.get_field() for var in variables]:
                    fields[field.name] = field

            self.fields = fields

            if domain is None:
                domain = self.fields.values()[0].domain

            self.domain = domain

            if conf is None:
                self.conf = Struct(ebcs={}, epbcs={}, lcbcs={})

        else:
            self.domain = domain
            self.fields = fields
            self.equations = equations

        if auto_solvers:
            if ls is None:
                ls = ScipyDirect({})

            if nls is None:
                nls = Newton({}, lin_solver=ls)

            ev = self.get_evaluator()
            nls.fun = ev.eval_residual
            nls.fun_grad = ev.eval_tangent_matrix

            self.solvers = Struct(name='solvers', ls=ls, nls=nls)

        self.setup_output()

    def reset(self):
        if hasattr(self.conf, 'options'):
            self.setup_hooks(self.conf.options)

        else:
            self.setup_hooks()

        self.mtx_a = None
        self.solvers = None
        self.clear_equations()

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
            - takes the ProblemDefinition instance (`self`) as the first
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

        iter_hook = self.nls_iter_hook
        if iter_hook is not None:
            self.nls_iter_hook = lambda *args, **kwargs: \
                                 iter_hook(self, *args, **kwargs)


    def copy(self, name=None):
        """
        Make a copy of ProblemDefinition.
        """
        if name is None:
            name = self.name + '_copy'
        obj = ProblemDefinition(name, conf=self.conf,
                                functions=self.functions,
                                domain=self.domain, fields=self.fields,
                                equations=self.equations,
                                auto_conf=False, auto_solvers=False)

        obj.ebcs = self.ebcs
        obj.epbcs = self.epbcs
        obj.lcbcs = self.lcbcs

        obj.set_solvers(self.conf.solvers, self.conf.options)

        obj.setup_output(output_filename_trunk=self.ofn_trunk,
                         output_dir=self.output_dir,
                         output_format=self.output_format,
                         file_per_var=self.file_per_var,
                         linearization=self.linearization)

        return obj

    def setup_default_output(self, conf=None, options=None):
        """
        Provide default values to `ProblemDefinition.setup_output()`
        from `conf.options` and `options`.
        """
        conf = get_default(conf, self.conf)

        if options and options.output_filename_trunk:
            default_output_dir, of = op.split(options.output_filename_trunk)
            default_trunk = io.get_trunk(of)

        else:
            default_trunk = None
            default_output_dir = get_default_attr(conf.options,
                                                  'output_dir', None)

        if options and options.output_format:
            default_output_format = options.output_format

        else:
            default_output_format = get_default_attr(conf.options,
                                                    'output_format', None)

        default_file_per_var = get_default_attr(conf.options,
                                                'file_per_var', None)
        default_float_format = get_default_attr(conf.options,
                                                'float_format', None)
        default_linearization = Struct(kind='strip')

        self.setup_output(output_filename_trunk=default_trunk,
                          output_dir=default_output_dir,
                          file_per_var=default_file_per_var,
                          output_format=default_output_format,
                          float_format=default_float_format,
                          linearization=default_linearization)

    def setup_output(self, output_filename_trunk=None, output_dir=None,
                     output_format=None, float_format=None,
                     file_per_var=None, linearization=None):
        """
        Sets output options to given values, or uses the defaults for
        each argument that is None.
        """
        self.output_modes = {'vtk' : 'sequence', 'h5' : 'single'}

        self.ofn_trunk = get_default(output_filename_trunk,
                                     io.get_trunk(self.domain.name))

        self.set_output_dir(output_dir)

        self.output_format = get_default(output_format, 'vtk')
        self.float_format = get_default(float_format, None)
        self.file_per_var = get_default(file_per_var, False)
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

    def set_regions( self, conf_regions=None,
                     conf_materials=None, functions=None):
        conf_regions = get_default(conf_regions, self.conf.regions)
        functions = get_default(functions, self.functions)

        self.domain.create_regions(conf_regions, functions)

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
            self.set_variables( conf_variables )

        return conf_variables

    def clear_equations( self ):
        self.integrals = None
        self.equations = None
        self.ebcs = None
        self.epbcs = None
        self.lcbcs = None

    def set_equations(self, conf_equations=None, user=None,
                      keep_solvers=False, make_virtual=False):
        """
        Set equations of the problem using the `equations` problem
        description entry.

        Fields and Regions have to be already set.
        """
        conf_equations = get_default(conf_equations,
                                     self.conf.get_default_attr('equations',
                                                                None))

        self.set_variables()
        variables = Variables.from_conf(self.conf_variables, self.fields)

        self.set_materials()
        materials = Materials.from_conf(self.conf_materials, self.functions)

        self.integrals = self.get_integrals()
        equations = Equations.from_conf(conf_equations, variables,
                                        self.domain.regions,
                                        materials, self.integrals,
                                        user=user,
                                        make_virtual=make_virtual)

        self.equations = equations

        if not keep_solvers:
            self.solvers = None

    def set_equations_instance(self, equations, keep_solvers=False):
        """
        Set equations of the problem to `equations`.
        """
        self.reset()
        self.equations = equations

        if not keep_solvers:
            self.solvers = None

    def set_solvers(self, conf_solvers=None, options=None):
        """
        Choose which solvers should be used. If solvers are not set in
        `options`, use first suitable in `conf_solvers`.
        """
        conf_solvers = get_default( conf_solvers, self.conf.solvers )
        self.solver_confs = {}
        for key, val in conf_solvers.iteritems():
            self.solver_confs[val.name] = val
        
        def _find_suitable( prefix ):
            for key, val in self.solver_confs.iteritems():
                if val.kind.find( prefix ) == 0:
                    return val
            return None

        def _get_solver_conf( kind ):
            try:
                key = options[kind]
                conf = self.solver_confs[key]
            except:
                conf = _find_suitable( kind + '.' )
            return conf

        self.ts_conf = _get_solver_conf( 'ts' )
        self.nls_conf = _get_solver_conf( 'nls' )
        self.ls_conf = _get_solver_conf( 'ls' )
        info =  'using solvers:'
        if self.ts_conf:
            info += '\n                ts: %s' % self.ts_conf.name
        if self.nls_conf:
            info += '\n               nls: %s' % self.nls_conf.name
        if self.ls_conf:
            info += '\n                ls: %s' % self.ls_conf.name
        if info != 'using solvers:':
            output( info )

    def set_solvers_instances(self, ls=None, nls=None):
        """
        Set the instances of linear and nonlinear solvers that will be
        used in `ProblemDefinition.solve()` call.
        """
        if (ls is not None) and (nls is not None):
            if not (nls.lin_solver is ls):
                raise ValueError('linear solver not used in nonlinear!')

        self.solvers = Struct(name='solvers', ls=ls, nls=nls)

    ##
    # Utility functions below.
    ##

    ##
    # 17.10.2007, c
    def get_solver_conf( self, name ):
        return self.solver_confs[name]

    ##
    # c: 13.06.2008, r: 13.06.2008
    def get_default_ts( self, t0 = None, t1 = None, dt = None, n_step = None,
                      step = None ):
        t0 = get_default( t0, 0.0 )
        t1 = get_default( t1, 1.0 )
        dt = get_default( dt, 1.0 )
        n_step = get_default( n_step, 1 )

        ts = TimeStepper(t0, t1, dt, n_step, step=step)

        return ts

    def get_integrals(self, names=None, kind=None):
        """
        Get integrals, initialized from problem configuration if available.

        Parameters
        ----------
        names : list, optional
            If given, only the named integrals are returned.
        kind : 'v' or 's', optional
            If given, only integrals of the given kind are returned.

        Returns
        -------
        integrals : Integrals instance
            The requested integrals.
        """
        conf_integrals = self.conf.get_default_attr('integrals', {})
        integrals = Integrals.from_conf(conf_integrals)

        if names is not None:
            integrals.update([integrals[ii] for ii in names
                              if ii in integrals.names])

        if kind is not None:
            integrals.update([ii for ii in integrals
                              if ii.kind.startswith(kind)])

        return integrals

    def update_time_stepper(self, ts):
        if ts is not None:
            self.ts.set_from_ts(ts)

    def update_materials(self, ts=None):
        """
        Update materials used in equations.
        """
        if self.equations is not None:
            self.update_time_stepper(ts)
            self.equations.time_update_materials(self.ts, self)


    def update_equations(self, ts=None, ebcs=None, epbcs=None,
                         lcbcs=None, functions=None, create_matrix=False):
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
        """
        self.update_time_stepper(ts)
        functions = get_default(functions, self.functions)

        graph_changed = self.equations.time_update(self.ts,
                                                   ebcs, epbcs, lcbcs,
                                                   functions, self)
        self.graph_changed = graph_changed

        if graph_changed or (self.mtx_a is None) or create_matrix:
            self.mtx_a = self.equations.create_matrix_graph()
            ## import sfepy.base.plotutils as plu
            ## plu.spy( self.mtx_a )
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

        if isinstance(epbcs, Conditions):
            self.epbcs = epbcs

        else:
            conf_epbc = get_default(epbcs, self.conf.epbcs)
            self.epbcs = Conditions.from_conf(conf_epbc, self.domain.regions)

        if isinstance(lcbcs, Conditions):
            self.lcbcs = lcbcs

        else:
            conf_lcbc = get_default(lcbcs, self.conf.lcbcs)
            self.lcbcs = Conditions.from_conf(conf_lcbc, self.domain.regions)

    def time_update(self, ts=None,
                    ebcs=None, epbcs=None, lcbcs=None,
                    functions=None, create_matrix=False):
        self.set_bcs(ebcs, epbcs, lcbcs)
        self.update_equations(ts, self.ebcs, self.epbcs, self.lcbcs,
                              functions, create_matrix)

    def setup_ic( self, conf_ics = None, functions = None ):
        conf_ics = get_default(conf_ics, self.conf.ics)
        ics = Conditions.from_conf(conf_ics, self.domain.regions)

        functions = get_default(functions, self.functions)

        self.equations.setup_initial_conditions(ics, functions)

    def select_bcs(self, ebc_names=None, epbc_names=None,
                   lcbc_names=None, create_matrix=False):

        if ebc_names is not None:
            conf_ebc = select_by_names( self.conf.ebcs, ebc_names )
        else:
            conf_ebc = None

        if epbc_names is not None:
            conf_epbc = select_by_names( self.conf.epbcs, epbc_names )
        else:
            conf_epbc = None

        if lcbc_names is not None:
            conf_lcbc = select_by_names( self.conf.lcbcs, lcbc_names )
        else:
            conf_lcbc = None

        self.set_bcs(conf_ebc, conf_epbc, conf_lcbc)
        self.update_equations(self.ts, self.ebcs, self.epbcs, self.lcbcs,
                              self.functions, create_matrix)

    def get_timestepper( self ):
        return self.ts

    def create_state(self):
        return State(self.equations.variables)

    ##
    # 26.07.2006, c
    # 22.08.2006
    def get_mesh_coors( self ):
        return self.domain.get_mesh_coors()

    def set_mesh_coors(self, coors, update_fields=False, actual=False,
                       clear_all=True):
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
        fea.set_mesh_coors(self.domain, self.fields, coors,
                           update_fields=update_fields, actual=actual,
                           clear_all=clear_all)

    def refine_uniformly(self, level):
        """
        Refine the mesh uniformly `level`-times.

        Notes
        -----
        This operation resets almost everything (fields, equations, ...)
        - it is roughly equivalent to creating a new ProblemDefinition
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

    def get_dim( self, get_sym = False ):
        """Returns mesh dimension, symmetric tensor dimension (if `get_sym` is
        True).
        """
        dim = self.domain.mesh.dim
        if get_sym:
            return dim, (dim + 1) * dim / 2
        else:
            return dim

    ##
    # c: 02.04.2008, r: 02.04.2008
    def init_time( self, ts ):
        self.update_time_stepper(ts)
        self.equations.init_time( ts )
        self.update_materials()

    def advance(self, ts=None):
        self.update_time_stepper(ts)
        self.equations.advance(self.ts)

    def save_state(self, filename, state=None, out=None,
                   fill_value=None, post_process_hook=None,
                   linearization=None, file_per_var=False, **kwargs):
        """
        Parameters
        ----------
        file_per_var : bool or None
            If True, data of each variable are stored in a separate
            file. If None, it is set to the application option value.
        linearization : Struct or None
            The linearization configuration for higher order
            approximations. If its kind is 'adaptive', `file_per_var` is
            assumed True.
        """
        linearization = get_default(linearization, self.linearization)
        if linearization.kind != 'adaptive':
            file_per_var = get_default(file_per_var, self.file_per_var)

        else:
            file_per_var = True

        extend = not file_per_var
        if (out is None) and (state is not None):
            out = state.create_output_dict(fill_value=fill_value,
                                           extend=extend,
                                           linearization=linearization)

            if post_process_hook is not None:
                out = post_process_hook(out, self, state, extend=extend)

        if linearization.kind == 'adaptive':
            for key, val in out.iteritems():
                mesh = val.get('mesh', self.domain.mesh)
                aux = io.edit_filename(filename, suffix='_' + val.var_name)
                mesh.write(aux, io='auto', out={key : val},
                           float_format=self.float_format, **kwargs)
                if hasattr(val, 'levels'):
                    output('max. refinement per group:', val.levels)

        elif file_per_var:
            meshes = {}

            if self.equations is None:
                varnames = {}
                for key, val in out.iteritems():
                    varnames[val.var_name] = 1
                varnames = varnames.keys()
                outvars = self.create_variables(varnames)
                itervars = outvars.__iter__
            else:
                itervars = self.equations.variables.iter_state

            for var in itervars():
                rname = var.field.region.name
                if meshes.has_key( rname ):
                    mesh = meshes[rname]
                else:
                    mesh = Mesh.from_region(var.field.region, self.domain.mesh,
                                            localize=True,
                                            is_surface=var.is_surface)
                    meshes[rname] = mesh

                vout = {}
                for key, val in out.iteritems():
                    try:
                        if val.var_name == var.name:
                            vout[key] = val

                    except AttributeError:
                        msg = 'missing var_name attribute in output!'
                        raise ValueError(msg)

                aux = io.edit_filename(filename, suffix='_' + var.name)
                mesh.write(aux, io='auto', out=vout,
                           float_format=self.float_format, **kwargs)
        else:
            self.domain.mesh.write(filename, io='auto', out=out,
                                   float_format=self.float_format, **kwargs)

    def save_ebc(self, filename, force=True, default=0.0):
        """
        Save essential boundary conditions as state variables.
        """
        output('saving ebc...')
        variables = self.get_variables(auto_create=True)

        ebcs = Conditions.from_conf(self.conf.ebcs, self.domain.regions)
        epbcs = Conditions.from_conf(self.conf.epbcs, self.domain.regions)

        try:
            variables.equation_mapping(ebcs, epbcs, self.ts, self.functions,
                                       problem=self)
        except:
            output('cannot make equation mapping!')
            raise

        state = self.create_state()
        state.fill(default)

        if force:
            vals = dict_from_keys_init(variables.state)
            for ii, key in enumerate(vals.iterkeys()):
                vals[key] = ii + 1

            state.apply_ebc(force_values=vals)

        else:
            state.apply_ebc()

        out = state.create_output_dict(extend=True)
        self.save_state(filename, out=out, fill_value=default)
        output('...done')

    def save_regions( self, filename_trunk, region_names = None ):
        """Save regions as meshes."""

        if region_names is None:
            region_names = self.domain.regions.get_names()

        output( 'saving regions...' )
        for name in region_names:
            region = self.domain.regions[name]
            output( name )
            aux = Mesh.from_region(region, self.domain.mesh)
            aux.write( '%s_%s.mesh' % (filename_trunk, region.name),
                       io = 'auto' )
        output( '...done' )

    def save_regions_as_groups(self, filename_trunk, region_names=None):
        """Save regions in a single mesh but mark them by using different
        element/node group numbers.

        If regions overlap, the result is undetermined, with exception of the
        whole domain region, which is marked by group id 0.

        Region masks are also saved as scalar point data for output formats
        that support this.

        Parameters
        ----------
        filename_trunk : str
            The output filename without suffix.
        region_names : list, optional
            If given, only the listed regions are saved.
        """

        output( 'saving regions as groups...' )
        aux = self.domain.mesh.copy()
        n_ig = c_ig = 0

        n_nod = self.domain.shape.n_nod

        # The whole domain region should go first.
        names = get_default(region_names, self.domain.regions.get_names())
        for name in names:
            region = self.domain.regions[name]
            if region.all_vertices.shape[0] == n_nod:
                names.remove(region.name)
                names = [region.name] + names
                break

        out = {}
        for name in names:
            region = self.domain.regions[name]
            output(region.name)

            aux.ngroups[region.all_vertices] = n_ig
            n_ig += 1

            mask = nm.zeros((n_nod, 1), dtype=nm.float64)
            mask[region.all_vertices] = 1.0
            out[name] = Struct(name = 'region',
                               mode = 'vertex', data = mask,
                               var_name = name, dofs = None)

            if region.has_cells():
                for ig in region.igs:
                    ii = region.get_cells(ig)
                    aux.mat_ids[ig][ii] = c_ig
                    c_ig += 1

        try:
            aux.write( '%s.%s' % (filename_trunk, self.output_format), io='auto',
                       out=out)
        except NotImplementedError:
            # Not all formats support output.
            pass

        output( '...done' )

    ##
    # c: 03.07.2007, r: 27.02.2008
    def save_field_meshes( self, filename_trunk ):

        output( 'saving field meshes...' )
        for field in self.fields:
            output( field.name )
            field.write_mesh( filename_trunk + '_%s' )
        output( '...done' )

    def get_evaluator(self, reuse=False):
        """
        Either create a new Evaluator instance (reuse == False),
        or return an existing instance, created in a preceding call to
        ProblemDefinition.init_solvers().
        """
        if reuse:
            try:
                ev = self.evaluator
            except AttributeError:
                raise AttributeError('call ProblemDefinition.init_solvers() or'\
                      ' set reuse to False!')
        else:
            if self.equations.variables.has_lcbc:
                ev = LCBCEvaluator(self, matrix_hook=self.matrix_hook)
            else:
                ev = BasicEvaluator(self, matrix_hook=self.matrix_hook)

        self.evaluator = ev

        return ev

    def init_solvers(self, nls_status=None, ls_conf=None, nls_conf=None,
                     mtx=None, presolve=False):
        """Create and initialize solvers."""
        ls_conf = get_default( ls_conf, self.ls_conf,
                               'you must set linear solver!' )

        nls_conf = get_default( nls_conf, self.nls_conf,
                              'you must set nonlinear solver!' )

        if presolve:
            tt = time.clock()
        if get_default_attr(ls_conf, 'needs_problem_instance', False):
            extra_args = {'problem' : self}
        else:
            extra_args = {}

        ls = Solver.any_from_conf(ls_conf, mtx=mtx, presolve=presolve,
                                  **extra_args)
        if presolve:
            tt = time.clock() - tt
            output('presolve: %.2f [s]' % tt)

        if get_default_attr(nls_conf, 'needs_problem_instance', False):
            extra_args = {'problem' : self}
        else:
            extra_args = {}
        ev = self.get_evaluator()

        if get_default_attr(self.conf.options, 'ulf', False):
            self.nls_iter_hook = ev.new_ulf_iteration

        nls = Solver.any_from_conf(nls_conf, fun=ev.eval_residual,
                                   fun_grad=ev.eval_tangent_matrix,
                                   lin_solver=ls, iter_hook=self.nls_iter_hook,
                                   status=nls_status, **extra_args)

        self.solvers = Struct( name = 'solvers', ls = ls, nls = nls )

    ##
    # c: 04.04.2008, r: 04.04.2008
    def get_solvers( self ):
        return getattr( self, 'solvers', None )

    ##
    # c: 04.04.2008, r: 04.04.2008
    def is_linear( self ):
        nls_conf = get_default(None, self.nls_conf,
                               'you must set nonlinear solver!')
        aux = Solver.any_from_conf(nls_conf)
        if aux.conf.problem == 'linear':
            return True
        else:
            return False

    ##
    # c: 13.06.2008, r: 13.06.2008
    def set_linear( self, is_linear ):
        nls_conf = get_default( None, self.nls_conf,
                              'you must set nonlinear solver!' )
        if is_linear:
            nls_conf.problem = 'linear'
        else:
            nls_conf.problem = 'nonlinear'

    def solve(self, state0=None, nls_status=None,
              ls_conf=None, nls_conf=None, force_values=None,
              var_data=None):
        """Solve self.equations in current time step.

        Parameters
        ----------
        var_data : dict
            A dictionary of {variable_name : data vector} used to initialize
            parameter variables.
        """
        solvers = self.get_solvers()
        if solvers is None:
            self.init_solvers(nls_status, ls_conf, nls_conf)
            solvers = self.get_solvers()

        if state0 is None:
            state0 = State(self.equations.variables)

        else:
            if isinstance(state0, nm.ndarray):
                state0 = State(self.equations.variables, vec=state0)

        self.equations.set_data(var_data, ignore_unknown=True)

        self.update_materials()
        state0.apply_ebc(force_values=force_values)

        vec0 = state0.get_reduced()

        vec = solvers.nls(vec0)

        state = state0.copy(preserve_caches=True)
        state.set_reduced(vec, preserve_caches=True)

        return state

    def create_evaluable(self, expression, try_equations=True, auto_init=False,
                         preserve_caches=False, copy_materials=True,
                         integrals=None,
                         ebcs=None, epbcs=None, lcbcs=None,
                         ts=None, functions=None,
                         mode='eval', var_dict=None, extra_args=None,
                         verbose=True, **kwargs):
        """
        Create evaluable object (equations and corresponding variables)
        from the `expression` string. Convenience function calling
        :func:`create_evaluable()
        <sfepy.fem.evaluate.create_evaluable()>` with defaults provided
        by the ProblemDefinition instance `self`.

        The evaluable can be repeatedly evaluated by calling
        :func:`eval_equations() <sfepy.fem.evaluate.eval_equations()>`,
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
            The variables (dictionary of (variable name) : (Variable
            instance)) to be used in the expression. Use this if the
            name of a variable conflicts with one of the parameters of
            this method.
        extra_args : dict, optional
            Extra arguments to be passed to terms in the expression.
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
            :func:`eval_equations() <sfepy.fem.evaluate.eval_equations()>`.

        Examples
        --------
        `problem` is ProblemDefinition instance.

        >>> out = problem.create_evaluable('dq_state_in_volume_qp.i1.Omega(u)')
        >>> equations, variables = out

        `vec` is a vector of coefficients compatible with the field
        of 'u' - let's use all ones.

        >>> vec = nm.ones((variables['u'].n_dof,), dtype=nm.float64)
        >>> variables['u'].data_from_any(vec)
        >>> vec_qp = eval_equations(equations, variables, mode='qp')

        Try another vector:

        >>> vec = 3 * nm.ones((variables['u'].n_dof,), dtype=nm.float64)
        >>> variables['u'].data_from_any(vec)
        >>> vec_qp = eval_equations(equations, variables, mode='qp')
        """
        from sfepy.fem.equations import get_expression_arg_names

        variables = get_default(var_dict, {})

        if try_equations and self.equations is not None:
            # Make a copy, so that possible variable caches are preserved.
            for key, var in self.equations.variables.as_dict().iteritems():
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
        if materials is not None:
            if copy_materials:
                materials = materials.semideep_copy()

            else:
                materials = Materials(objs=materials._objs)

        else:
            possible_mat_names = get_expression_arg_names(expression)
            materials = self.create_materials(possible_mat_names)

        _kwargs = copy(kwargs)
        for key, val in kwargs.iteritems():
            if isinstance(val, Variable):
                if val.name != key:
                    msg = 'inconsistent variable name! (%s == %s)' \
                          % (val.name, key)
                    raise ValueError(msg)
                variables[val.name] = val
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
                               variables.itervalues(), integrals,
                               ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs,
                               ts=ts, functions=functions,
                               auto_init=auto_init,
                               mode=mode, extra_args=extra_args, verbose=verbose,
                               kwargs=kwargs)

        if copy_materials:
            equations = out[0]
            equations.time_update_materials(self.ts, self, verbose=verbose)

        return out

    def evaluate(self, expression, try_equations=True, auto_init=False,
                 preserve_caches=False, copy_materials=True, integrals=None,
                 ebcs=None, epbcs=None, lcbcs=None,
                 ts=None, functions=None,
                 mode='eval', dw_mode='vector', term_mode=None,
                 var_dict=None, ret_variables=False, verbose=True,
                 extra_args=None, **kwargs):
        """
        Evaluate an expression, convenience wrapper of
        `self.create_evaluable()` and
        :func:`eval_equations() <sfepy.fem.evaluate.eval_equations()>`.

        Parameters
        ----------
        ... : arguments
            See docstrings of `self.create_evaluable()`.
        dw_mode : 'vector' or 'matrix'
            The assembling mode for 'weak' evaluation mode.
        term_mode : str
            The term call mode - some terms support different call modes
            and depending on the call mode different values are
            returned.
        ret_variables : bool
            If True, return the variables that were created to evaluate
            the expression.

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
                                    extra_args=extra_args,
                                    verbose=verbose, **kwargs)
        equations, variables = aux

        out = eval_equations(equations, variables,
                             preserve_caches=preserve_caches,
                             mode=mode, dw_mode=dw_mode, term_mode=term_mode)

        if ret_variables:
            out = (out, variables)

        return out

    def get_time_solver(self, ts_conf=None, **kwargs):
        """
        Create and return a TimeSteppingSolver instance.

        Notes
        -----
        Also sets `self.ts` attribute.
        """
        ts_conf = get_default(ts_conf, self.ts_conf,
                              'you must set time-stepping solver!')
        ts_solver = Solver.any_from_conf(ts_conf, **kwargs)
        self.ts = ts_solver.ts

        return ts_solver

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

    def init_variables( self, state ):
        """Initialize variables with history."""
        self.equations.variables.init_state( state )

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
