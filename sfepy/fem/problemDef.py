import os.path as op

import sfepy
from sfepy.base.base import *

import sfepy.base.ioutils as io
from sfepy.base.conf import ProblemConf, get_standard_keywords, transform_variables
from functions import Functions
from mesh import Mesh
from domain import Domain
from fields import fields_from_conf
from variables import Variables, Variable
from materials import Materials
from equations import Equations
from integrals import Integrals
from sfepy.fem.conditions import Conditions
from sfepy.fem.evaluate import create_evaluable, eval_equations
import fea as fea
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.evaluate import BasicEvaluator, LCBCEvaluator, evaluate
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

        domain = Domain(mesh.name, mesh)

        obj = ProblemDefinition('prblem_from_conf', conf=conf,
                                functions=functions, domain=domain,
                                auto_conf=False, auto_solvers=False)

        obj.set_regions(conf.regions, conf.materials, obj.functions)

        obj.clear_equations()

        if init_fields:
            obj.set_fields( conf.fields )

	    if init_equations:
		obj.set_equations( conf.equations )

        if init_solvers:
            obj.set_solvers( conf.solvers, conf.options )


        obj.ts = None

        return obj

    def __init__(self, name, conf=None, functions=None,
                 domain=None, fields=None, materials=None,
                 equations=None, auto_conf=True,
                 nls=None, ls=None, ts=None, auto_solvers=True):
        self.name = name
        self.conf = conf
        self.functions = functions

        self.reset()

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

            if materials is None:
                materials = Materials(self.equations.collect_materials())

            self.domain = domain
            self.materials = materials

            if conf is None:
                self.conf = Struct(ebcs={}, epbcs={}, lcbcs={})

        else:
            self.domain = domain
            self.fields = fields
            self.materials = materials
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
        self.mtx_a = None
        self.solvers = None
        self.clear_equations()

    def copy(self, name=None):
        """
        Make a copy of ProblemDefinition.
        """
        if name is None:
            name = self.name + '_copy'
        obj = ProblemDefinition(name, conf=self.conf,
                                functions=self.functions,
                                domain=self.domain, fields=self.fields,
                                materials=self.materials,
                                equations=self.equations,
                                auto_conf=False, auto_solvers=False)

        obj.ebcs = self.ebcs
        obj.epbcs = self.epbcs
        obj.lcbcs = self.lcbcs

        obj.set_solvers(self.conf.solvers, self.conf.options)

        return obj

    def setup_output(self, output_filename_trunk=None, output_dir=None,
                     output_format=None, float_format=None):
        """
        Sets output options to given values, or uses the defaults for
        each argument that is None.
        """
        self.output_modes = {'vtk' : 'sequence', 'h5' : 'single'}

	self.ofn_trunk = get_default(output_filename_trunk,
                                     io.get_trunk(self.domain.name))

        self.output_dir = get_default(output_dir, '.')

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.output_format = get_default(output_format, 'vtk')

        try:
            _float_format = get_default_attr(self.conf.options,
                                             'float_format', None)

        except AttributeError:
            _float_format = None

        self.float_format = get_default(float_format, _float_format)


    def set_regions( self, conf_regions=None,
                     conf_materials=None, functions=None):
        conf_regions = get_default(conf_regions, self.conf.regions)
        conf_materials = get_default(conf_materials, self.conf.materials)
        functions = get_default(functions, self.functions)

        self.domain.create_regions(conf_regions, functions)

        materials = Materials.from_conf(conf_materials, functions)
        self.materials = materials

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
                      cache_override=None,
                      keep_solvers=False, make_virtual=False):
        """
        Set equations of the problem. Regions, Variables and Materials
        have to be already set.
        """
        conf_equations = get_default(conf_equations,
                                     self.conf.get_default_attr('equations',
                                                                None))
        if cache_override is None:
            cache_override = get_default_attr( self.conf.fe,
                                               'cache_override', True )

	self.set_variables()
        variables = Variables.from_conf(self.conf_variables, self.fields)

        self.integrals = Integrals.from_conf(self.conf.integrals)
        equations = Equations.from_conf(conf_equations, variables,
                                        self.domain.regions,
                                        self.materials, self.integrals,
                                        user=user,
                                        cache_override=cache_override,
                                        make_virtual=make_virtual)

        self.equations = equations

        if not keep_solvers:
            self.solvers = None

    ##
    # c: 16.10.2007, r: 20.02.2008
    def set_solvers( self, conf_solvers = None, options = None ):
        """If solvers are not set in options, use first suitable in
        conf_solvers."""
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

    ##
    # Utility functions below.
    ##

    ##
    # 17.10.2007, c
    def get_solver_conf( self, name ):
        return self.solver_confs[name]
    
    ##
    # 29.01.2006, c
    # 25.07.2006
    def create_state_vector( self ):
        return self.equations.create_state_vector()

    ##
    # c: 13.06.2008, r: 13.06.2008
    def get_default_ts( self, t0 = None, t1 = None, dt = None, n_step = None,
                      step = None ):
        t0 = get_default( t0, 0.0 )
        t1 = get_default( t1, 1.0 )
        dt = get_default( dt, 1.0 )
        n_step = get_default( n_step, 1 )
        ts = TimeStepper( t0, t1, dt, n_step )
        ts.set_step( step )
        return ts

    def reset_materials(self):
        """Clear material data so that next materials.time_update() is
        performed even for stationary materials."""
        self.materials.reset()

    def update_materials(self, ts=None):
        if ts is None:
            ts = self.get_default_ts(step=0)

        self.materials.time_update(ts, self.domain, self.equations)

    def update_equations(self, ts=None, ebcs=None, epbcs=None,
                         lcbcs=None, functions=None, create_matrix=False):
        """
        Assumes same EBC/EPBC/LCBC nodes for all time steps. Otherwise set
        create_matrix to True.
        """
        if ts is None:
            ts = self.get_default_ts(step=0)
        functions = get_default(functions, self.functions)

        self.equations.time_update(ts, ebcs, epbcs, lcbcs, functions)

        if (self.mtx_a is None) or create_matrix:
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
        if ts is None:
            ts = self.get_default_ts( step = 0 )

        self.ts = ts
        self.update_materials(ts)
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

    ##
    # 29.01.2006, c
    # 25.07.2006
    # 19.09.2006
    def apply_ebc( self, vec, force_values = None ):
        """Apply essential (Dirichlet) boundary conditions."""
        self.equations.apply_ebc( vec, force_values )

    def apply_ic( self, vec, force_values = None ):
        """Apply initial conditions."""
        self.equations.apply_ic( vec, force_values )

    ##
    # 25.07.2006, c
    def update_vec( self, vec, delta ):
        self.equations.update_vec( vec, delta )

    ##
    # c: 18.04.2006, r: 07.05.2008
    def state_to_output( self, vec, fill_value = None, var_info = None,
                       extend = True ):
        """
        Transforms state vector 'vec' to an output dictionary, that can be
        passed as 'out' kwarg to Mesh.write(). 'vec' must have full size,
        i.e. all fixed or periodic values must be included.

        Example:
        >>> out = problem.state_to_output( state )
        >>> problem.save_state( 'file.vtk', out = out )

        Then the  dictionary entries a formed by components of the state vector
        corresponding to the unknown variables, each transformed to shape
        (n_mesh_nod, n_dof per node) - all values in extra nodes are removed.
        """
        return self.equations.state_to_output(vec, fill_value,
                                              var_info, extend)

    ##
    # 26.07.2006, c
    # 22.08.2006
    def get_mesh_coors( self ):
        return self.domain.get_mesh_coors()

    def set_mesh_coors(self, coors, update_state=False):
        if self.equations is None:
            geom = {}
        else:
            geom = self.equations.geometries
        fea.set_mesh_coors(self.domain, self.fields, geom,
                           coors, update_state)

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
        self.equations.init_time( ts )

    ##
    # 08.06.2007, c
    def advance( self, ts ):
        self.equations.advance( ts )

    ##
    # c: 01.03.2007, r: 23.06.2008
    def save_state( self, filename, state = None, out = None,
                   fill_value = None, post_process_hook = None,
                   file_per_var = False, **kwargs ):
        extend = not file_per_var
        if (out is None) and (state is not None):
            out = self.state_to_output( state,
                                      fill_value = fill_value, extend = extend )
            if post_process_hook is not None:
                out = post_process_hook( out, self, state, extend = extend )

        if file_per_var:
            import os.path as op

            meshes = {}
            for var in self.equations.variables.iter_state():
                rname = var.field.region.name
                if meshes.has_key( rname ):
                    mesh = meshes[rname]
                else:
                    mesh = Mesh.from_region( var.field.region, self.domain.mesh,
                                            localize = True )
                    meshes[rname] = mesh
                vout = {}
                for key, val in out.iteritems():
                    if val.var_name == var.name:
                        vout[key] = val
                base, suffix = op.splitext( filename )
                mesh.write(base + '_' + var.name + suffix,
                           io='auto', out=vout,
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

        ebcs = Conditions.from_conf(self.conf.ebcs)
        epbcs = Conditions.from_conf(self.conf.epbcs)

        try:
            ts = TimeStepper.from_conf(self.conf.ts)
            ts.set_step(0)

        except:
            ts = None

        try:
            variables.equation_mapping(ebcs, epbcs,
                                       self.domain.regions, ts, self.functions)
        except Exception, e:
            output( 'cannot make equation mapping!' )
            output( 'reason: %s' % e )

        state = variables.create_state_vector()
        state.fill(default)

        if force:
            vals = dict_from_keys_init(variables.state)
            for ii, key in enumerate(vals.iterkeys()):
                vals[key] = ii + 1

            variables.apply_ebc(state, force_values=vals)

        else:
            variables.apply_ebc(state)

        out = variables.state_to_output(state, extend=True)
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
            aux = Mesh.from_region( region, self.domain.mesh, self.domain.ed,
                                   self.domain.fa )
            aux.write( '%s_%s.mesh' % (filename_trunk, region.name),
                       io = 'auto' )
        output( '...done' )

    def save_regions_as_groups(self, filename_trunk):
	"""Save regions in a single mesh but mark them by using different
        element/node group numbers.

        If regions overlap, the result is undetermined, with exception of the
        whole domain region, which is marked by group id 0.

        Region masks are also saved as scalar point data for output formats
        that support this.
        """

        output( 'saving regions as groups...' )
        aux = self.domain.mesh.copy()
        n_ig = c_ig = 0

        n_nod = self.domain.shape.n_nod

        # The whole domain region should go first.
        names = self.domain.regions.get_names()
        for region in self.domain.regions:
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
    # created:       02.01.2008
    # last revision: 27.02.2008
    def save_region_field_meshes( self, filename_trunk ):

        output( 'saving regions of fields...' )
        for field in self.fields:
            fregion = self.domain.regions[field.region_name]
            output( 'field %s: saving regions...' % field.name )

            for region in self.domain.regions:
                if not fregion.contains( region ): continue
                output( region.name )
                aux = Mesh.from_region_and_field( region, field )
                aux.write( '%s_%s_%s.mesh' % (filename_trunk,
                                              region.name, field.name),
                           io = 'auto' )
            output( '...done' )
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
                ev = LCBCEvaluator(self)
            else:
                ev = BasicEvaluator(self)

        self.evaluator = ev
        
        return ev

    def init_solvers(self, nls_status=None, ls_conf=None, nls_conf=None,
                     mtx=None, presolve=False):
        """Create and initialize solvers."""
        ls_conf = get_default( ls_conf, self.ls_conf,
                               'you must set linear solver!' )

        nls_conf = get_default( nls_conf, self.nls_conf,
                              'you must set nonlinear solver!' )
        
        ls = Solver.any_from_conf(ls_conf, mtx=mtx, presolve=presolve)

        if get_default_attr(nls_conf, 'needs_problem_instance', False):
            extra_args = {'problem' : self}
        else:
            extra_args = {}
        ev = self.get_evaluator()
        nls = Solver.any_from_conf( nls_conf, fun = ev.eval_residual,
                                    fun_grad = ev.eval_tangent_matrix,
                                    lin_solver = ls, status = nls_status,
                                    **extra_args )

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
            state = self.create_state_vector()
        else:
            state = state0.copy()

        self.equations.set_data(var_data, ignore_unknown=True)

        self.apply_ebc( state, force_values = force_values )

        ev = self.evaluator

        vec0 = ev.strip_state_vector( state )
        vec = solvers.nls( vec0 )
        state = ev.make_full_vec( vec )
        
        return state

    def create_evaluable(self, expression, try_equations=True, auto_init=False,
                         copy_materials=True, integrals=None,
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
        copy_materials : bool
            Work with a copy of `self.materials`. Safe but can be slow.
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
        if try_equations and self.equations is not None:
            variables = self.equations.variables.as_dict()

        else:
            if var_dict is None:
                from sfepy.fem.equations import parse_definition

                args = ','.join(aux.args
                                for aux in parse_definition(expression))
                aux = [arg.strip() for arg in args.split(',')]
                possible_var_names = set(aux)

                variables = self.create_variables(possible_var_names)

            else:
                variables = var_dict

        _kwargs = copy(kwargs)
        for key, val in kwargs.iteritems():
            if isinstance(val, Variable):
                variables[val.name] = val
                _kwargs.pop(key)
        kwargs = _kwargs

        if copy_materials:
            materials = self.materials.semideep_copy()

        else:
            materials = self.materials

        ebcs = get_default(ebcs, self.ebcs)
        epbcs = get_default(epbcs, self.epbcs)
        lcbcs = get_default(lcbcs, self.lcbcs)
        ts = get_default(ts, self.get_timestepper())
        functions = get_default(functions, self.functions)
        integrals = get_default(integrals,
                                Integrals.from_conf(self.conf.integrals))

        out = create_evaluable(expression, self.fields, materials,
                               variables.itervalues(), integrals,
                               update_materials=copy_materials,
                               ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs,
                               ts=ts, functions=functions,
                               auto_init=auto_init,
                               mode=mode, extra_args=extra_args, verbose=verbose,
                               kwargs=kwargs)

        return out

    def evaluate(self, expression, try_equations=True, auto_init=False,
                 copy_materials=True, integrals=None,
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
                                    copy_materials=copy_materials,
                                    integrals=integrals,
                                    ebcs=ebcs, epbcs=epbcs, lcbcs=lcbcs,
                                    ts=ts, functions=functions,
                                    mode=mode, var_dict=var_dict,
                                    extra_args=None,
                                    verbose=verbose, **kwargs)
        equations, variables = aux

        out = eval_equations(equations, variables,
                             mode=mode, dw_mode=dw_mode, term_mode=term_mode)

        if ret_variables:
            out = (out, variables)

        return out

    ##
    # c: 06.02.2008, r: 04.04.2008
    def get_time_solver( self, ts_conf = None, **kwargs ):
        ts_conf = get_default( ts_conf, self.ts_conf,
                             'you must set time-stepping solver!' )
        
        return Solver.any_from_conf( ts_conf, **kwargs )


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
	variables.setup_dof_info()

	return variables

    def get_output_name(self, suffix=None, extra=None, mode=None):
        """Return default output file name, based on the output format,
        step suffix and mode. If present, the extra string is put just before
        the output format suffix.
        """
        out = self.ofn_trunk
        if suffix is not None:
            if mode is None:
                mode = self.output_modes[self.output_format]

            if mode == 'sequence':
                out = '.'.join((self.ofn_trunk, suffix))

        if extra is not None:
            out = '.'.join((out, extra, self.output_format))
        else:
            out = '.'.join((out, self.output_format))

        return out
