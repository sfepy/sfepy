import os.path as op

import sfepy
from sfepy.base.base import *

import sfepy.base.ioutils as io
from sfepy.base.conf import ProblemConf, get_standard_keywords
from functions import Functions
from mesh import Mesh
from domain import Domain
from fields import Fields
from variables import Variables
from materials import Materials
from equations import Equations
from integrals import Integrals
import fea as fea
from sfepy.solvers.ts import TimeStepper
from sfepy.fem.evaluate import BasicEvaluator, LCBCEvaluator, eval_term_op
from sfepy.solvers import Solver

##
# 29.01.2006, c
class ProblemDefinition( Struct ):
    """
    Problem definition, the top-level class holding all data necessary to solve
    a problem.

    Contains: mesh, domain, materials, fields, variables, equations, solvers
    """

    def from_conf_file(conf_filename,
                       required=None, other=None,
                       init_fields = True,
                       init_variables = True,
                       init_equations = True,
                       init_solvers = True):

        _required, _other = get_standard_keywords()
        if required is None:
            required = _required
        if other is None:
            other = _other
            
        conf = ProblemConf.from_file(conf_filename, required, other)

        obj = ProblemDefinition.from_conf(conf,
                                          init_fields=init_fields,
                                          init_variables=init_variables,
                                          init_equations=init_equations,
                                          init_solvers=init_solvers)
        return obj
    from_conf_file = staticmethod(from_conf_file)
    
    def from_conf( conf,
                   init_fields = True,
                   init_variables = True,
                   init_equations = True,
                   init_solvers = True ):
        if conf.options.get_default_attr('absolute_mesh_path', False):
            conf_dir = None
        else:
            conf_dir = op.dirname(conf.funmod.__file__)

        functions = Functions.from_conf(conf.functions)
            
        mesh = Mesh.from_file(conf.filename_mesh, prefix_dir=conf_dir)

        eldesc_dir = op.join( sfepy.base_dir, 'eldesc' )
        domain = Domain.from_mesh( mesh, eldesc_dir )
        domain.setup_groups()
        domain.fix_element_orientation()
        domain.setup_neighbour_lists()

        obj = ProblemDefinition( conf = conf,
                                 functions = functions,
                                 domain = domain,
                                 eldesc_dir = eldesc_dir )

	# Default output file trunk and format.
	obj.ofn_trunk = io.get_trunk( conf.filename_mesh )
        obj.output_format = 'vtk'

        obj.set_regions(conf.regions, conf.materials, obj.functions)

        if init_fields:
            obj.set_fields( conf.fields )

            if init_variables:
                obj.set_variables( conf.variables )

                if init_equations:
                    obj.set_equations( conf.equations )

        if init_solvers:
            obj.set_solvers( conf.solvers, conf.options )


        obj.ts = None
        
        return obj
    from_conf = staticmethod( from_conf )

    ##
    # 18.04.2006, c
    def copy( self, **kwargs ):
        if 'share' in kwargs:
            share = kwargs['share']
            
        obj = ProblemDefinition()
        for key, val in self.__dict__.iteritems():
##             print key
            if key in share:
                obj.__dict__[key] = val
            else:
                obj.__dict__[key] = copy( val )
        return obj

    def set_regions( self, conf_regions=None,
                     conf_materials=None, functions=None):
        conf_regions = get_default(conf_regions, self.conf.regions)
        conf_materials = get_default(conf_materials, self.conf.materials)
        functions = get_default(functions, self.functions)

        self.domain.create_regions(conf_regions, functions)

        materials = Materials.from_conf(conf_materials, functions)
        materials.setup_regions(self.domain.regions)
        self.materials = materials

    ##
    # c: 23.04.2007, r: 09.07.2008
    def set_fields( self, conf_fields = None ):
        conf_fields = get_default( conf_fields, self.conf.fields )
        fields = Fields.from_conf( conf_fields )
        fields.read_interpolants( self.eldesc_dir )
        fields.setup_approximations( self.domain )
##         print fields
##         print fields[0].aps
##         pause()
        fields.setup_global_base()
        fields.setup_coors()
        self.fields = fields

##         self.save_field_meshes( '.' )
##         pause()

    ##
    # c: 26.07.2006, r: 14.04.2008
    def set_variables( self, conf_variables = None ):
        conf_variables = get_default( conf_variables, self.conf.variables )
        variables = Variables.from_conf( conf_variables, self.fields )
        variables.setup_dof_info() # Call after fields.setup_global_base().
        self.variables = variables
        self.mtx_a = None
        self.solvers = None
        self.clear_equations()
##         print variables.di
##         pause()
        

    def select_variables( self, variable_names ):
        conf_variables = select_by_names( self.conf.variables, variable_names )
        self.set_variables( conf_variables )

    def clear_equations( self ):
        self.integrals = None
        self.geometries = None
        self.equations = None
    
    ##
    # c: 18.04.2006, r: 13.06.2008
    def set_equations( self, conf_equations = None, user = None,
                       cache_override = None,
                       keep_solvers = False, make_virtual = False,
                       single_term = False ):
        conf_equations = get_default( conf_equations,
                                      self.conf.get_default_attr('equations',
                                                                 None) )
        equations = Equations.from_conf( conf_equations )
        equations.setup_terms( self.domain.regions, self.variables,
                               self.materials, user )

        # This uses the actual conn_info created in equations.setup_terms().
        self.variables.setup_dof_conns(make_virtual=make_virtual,
                                       single_term=single_term)

        i_names = equations.get_term_integral_names()
        self.integrals = Integrals.from_conf( self.conf.integrals, i_names )
        self.integrals.set_quadratures( fea.collect_quadratures() )

        self.geometries = {}
        equations.describe_geometry( self.geometries, self.variables,
                                     self.integrals )

##         print self.geometries
##         pause()

        if cache_override is None:
            cache_override = get_default_attr( self.conf.fe,
                                               'cache_override', True )
        equations.set_cache_mode( cache_override )

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
        return self.variables.create_state_vector()

    def update_bc( self, ts, conf_ebc, conf_epbc, conf_lcbc, functions,
                   create_matrix = False ):
        """Assumes same EBC/EPBC/LCBC nodes for all time steps. Otherwise set
        create_matrix to True."""
        self.variables.equation_mapping(conf_ebc, conf_epbc,
                                        self.domain.regions, ts, functions)
        self.variables.setup_lcbc_operators(conf_lcbc, self.domain.regions)
                
        self.variables.setup_a_dof_conns()
        if (self.mtx_a is None) or create_matrix:
            self.mtx_a = self.variables.create_matrix_graph()
##             import sfepy.base.plotutils as plu
##             plu.spy( self.mtx_a )
##             plu.pylab.show()

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

    def update_materials(self, ts=None):
        if ts is None:
            ts = self.get_default_ts(step=0)

        self.materials.time_update(ts, self.domain,
                                   self.equations, self.variables)

    def update_equations( self, ts = None ):
        if ts is None:
            ts = self.get_default_ts( step = 0 )
        self.equations.time_update( ts )
        self.variables.time_update( ts )

    def time_update( self, ts = None,
                     conf_ebc = None, conf_epbc = None, conf_lcbc = None,
                     functions = None, create_matrix = False ):
        if ts is None:
            ts = self.get_default_ts( step = 0 )

        self.ts = ts
        conf_ebc = get_default( conf_ebc, self.conf.ebcs )
        conf_epbc = get_default( conf_epbc, self.conf.epbcs )
        conf_lcbc = get_default( conf_lcbc, self.conf.lcbcs )
        functions = get_default(functions, self.functions)
        self.update_bc(ts, conf_ebc, conf_epbc, conf_lcbc, functions,
                       create_matrix)
        self.update_materials(ts)
        self.update_equations(ts)

    def setup_ic( self, conf_ics = None, functions = None ):
        conf_ics = get_default(conf_ics, self.conf.ics)
        functions = get_default(functions, self.functions)
        self.variables.setup_initial_conditions(conf_ics,
                                                self.domain.regions, functions)

    def select_bcs( self, ts = None, ebc_names = None, epbc_names = None,
                    lcbc_names = None ):
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

        self.time_update( ts, conf_ebc, conf_epbc, conf_lcbc )

    def get_timestepper( self ):
        return self.ts

    ##
    # 29.01.2006, c
    # 25.07.2006
    # 19.09.2006
    def apply_ebc( self, vec, force_values = None ):
        """Apply essential (Dirichlet) boundary conditions."""
        self.variables.apply_ebc( vec, force_values )

    def apply_ic( self, vec, force_values = None ):
        """Apply initial conditions."""
        self.variables.apply_ic( vec, force_values )

    ##
    # 25.07.2006, c
    def update_vec( self, vec, delta ):
        self.variables.update_vec( vec, delta )
        
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
        return self.variables.state_to_output( vec, fill_value,
                                               var_info, extend )

    ##
    # 26.07.2006, c
    # 22.08.2006
    def get_mesh_coors( self ):
        return self.domain.get_mesh_coors()

    ##
    # created: 26.07.2006
    # last revision: 21.12.2007
    def set_mesh_coors( self, coors, update_state = False ):
        fea.set_mesh_coors( self.domain, self.fields, self.geometries,
                          coors, update_state )

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
        self.variables.advance( ts )

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

        float_format = get_default_attr( self.conf.options,
                                         'float_format', None )

        if file_per_var:
            import os.path as op

            meshes = {}
            for var in self.variables.iter_state():
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
                mesh.write( base + '_' + var.name + suffix,
                            io = 'auto', out = vout,
                            float_format = float_format, **kwargs )
        else:
            self.domain.mesh.write( filename, io = 'auto', out = out,
                                    float_format = float_format, **kwargs )

    ##
    # c: 19.09.2006, r: 27.02.2008
    def save_ebc( self, filename, force = True, default = 0.0 ):
        output( 'saving ebc...' )
        state = self.create_state_vector()
        state.fill( default )
        if force:
            vals = dict_from_keys_init( [self.variables.names[ii]
                                      for ii in self.variables.state] )
            for ii, key in enumerate( vals.iterkeys() ):
                vals[key] = ii + 1
            self.apply_ebc( state, force_values = vals )
        else:
            self.apply_ebc( state )
        self.save_state( filename, state, fill_value = default )
        output( '...done' )

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

    def get_evaluator( self, reuse = False, **kwargs ):
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
            if self.variables.has_lcbc:
                ev = LCBCEvaluator( self, **kwargs )
            else:
                ev = BasicEvaluator( self, **kwargs )

        self.evaluator = ev
        
        return ev

    def init_solvers( self, nls_status = None, ls_conf = None, nls_conf = None,
                      mtx = None, **kwargs ):
        """Create and initialize solvers."""
        ls_conf = get_default( ls_conf, self.ls_conf,
                               'you must set linear solver!' )

        nls_conf = get_default( nls_conf, self.nls_conf,
                              'you must set nonlinear solver!' )
        
        ls = Solver.any_from_conf( ls_conf, mtx = mtx )

        if get_default_attr(nls_conf, 'needs_problem_instance', False):
            extra_args = {'problem' : self}
        else:
            extra_args = {}
        ev = self.get_evaluator( **kwargs )
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

    def solve( self, state0 = None, nls_status = None,
               ls_conf = None, nls_conf = None, force_values = None,
               **kwargs ):
        """Solve self.equations in current time step."""
        solvers = self.get_solvers()
        if solvers is None:
            self.init_solvers( nls_status, ls_conf, nls_conf, **kwargs )
            solvers = self.get_solvers()
        else:
            if kwargs:
                ev = self.get_evaluator( reuse = True )
                ev.set_term_args( **kwargs )
            
        if state0 is None:
            state = self.create_state_vector()
        else:
            state = state0.copy()

        self.apply_ebc( state, force_values = force_values )

        ev = self.evaluator

        vec0 = ev.strip_state_vector( state )
        vec = solvers.nls( vec0 )
        state = ev.make_full_vec( vec )
        
        return state

    def evaluate(self, expression, state=None, **kwargs):
        """
        Evaluate an expression, convenience wrapper of eval_term_op().

        Parameters
        ----------
        state .. state vector
        kwargs .. dictionary of {variable_name : data vector}, and other
        arguments supported by eval_term()


        Examples
        --------
        >>> problem.evaluate("di_volume_integrate.i1.Omega(Psi)",
        ... Psi=data['n'].data)
        array([ 5.68437535])
        """
        if state is None:
            vargs = {}
            for key, val in kwargs.iteritems():
                if self.variables.has_key(key):
                    vargs[key] = val
            out = eval_term_op(vargs, expression, self, **kwargs)

        else:
            out = eval_term_op(state, expression, self, **kwargs)
            
        return out

    ##
    # c: 06.02.2008, r: 04.04.2008
    def get_time_solver( self, ts_conf = None, **kwargs ):
        ts_conf = get_default( ts_conf, self.ts_conf,
                             'you must set time-stepping solver!' )
        
        return Solver.any_from_conf( ts_conf, **kwargs )


    def init_variables( self, state ):
        """Initialize variables with history."""
        self.variables.init_state( state )

    def get_output_name( self, suffix = None ):
        """Return default output file name."""
        if suffix is None:
            return '.'.join( (self.ofn_trunk, self.output_format) )
        else:
            return '.'.join( (self.ofn_trunk, suffix, self.output_format) )
