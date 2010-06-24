from sfepy.base.base import *
from sfepy.fem import Materials, Variables
from extmods.fem import raw_graph
from sfepy.terms import Terms, Term, DataCaches

"""
Note:
- create found materials, variables from configuration/input file data
  ... no - user should be able to create objects even if they are not
  used in equations
"""

def parse_definition(equation_def, itps):
    """
    Parse equation definition string to create term description list.
    """
    from parseEq import create_bnf

    term_descs = []
    bnf = create_bnf(term_descs, itps)
    try:
        bnf.parseString(equation_def)
    except:
        raise ValueError('cannot parse equation! (%s)' % equation_def)

    return term_descs

##
# 21.07.2006, c
class Equations( Container ):

    @staticmethod
    def from_conf(conf, variables, regions, materials, caches=None,
                  user=None, term_prefixes=None):

        objs = OneTypeList(Equation)

        conf = copy(conf)
        tps = conf.pop('namespaces', {})

        if caches is None:
            caches = DataCaches()

        ii = 0
        for name, desc in conf.iteritems():
            output('equation "%s":' %  name)
            output(desc)
            eq = Equation.from_desc(name, desc, variables, regions,
                                    materials, caches=caches,
                                    user=user, term_prefixes=tps)
            objs.append(eq)
            ii += 1

        obj = Equations(objs, caches=caches)

        return obj

    def __init__(self, equations, caches=None):
        Container.__init__(self, equations)

        self.variables = Variables(self.collect_variables())
        self.variables.setup_dof_info() # Call after fields.setup_global_base().

        self.caches = get_default(caches, DataCaches())

        self.clear_geometries()

    def clear_geometries(self):
        self.geometries = {}

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

    def describe_geometry(self, integrals, verbose=True):
        output('describing geometries...', verbose=verbose)
        tt = time.clock()
        for eq in self:
            eq.describe_geometry(self.geometries, integrals)
        output('...done in %.2f s' % (time.clock() - tt), verbose=verbose)
        
    def get_variable_names( self ):
        """Return the list of names of all variables used in equations."""
        vns = set()
        for eq in self:
            for term in eq.terms:
                vns.update( term.get_variable_names() )
        return list( vns )

    ##
    # 27.02.2007, c
    def invalidate_term_caches( self ):
        for cache in self.caches.itervalues():
            cache.clear()

    ##
    # c: 07.05.2008, r: 07.05.2008
    def reset_term_caches( self ):
        for cache in self.caches.itervalues():
            cache.reset()

    ##
    # 02.03.2007, c
    def set_cache_mode( self, cache_override ):
        for cache in self.caches.itervalues():
            cache.set_mode( cache_override )

    def time_update(self, ts, regions, ebcs=None, epbcs=None, lcbcs=None,
                    functions=None):
        self.variables.time_update(ts, functions)

        self.variables.equation_mapping(ebcs, epbcs, regions, ts, functions)
        self.variables.setup_lcbc_operators(lcbcs, regions)
        self.variables.setup_adof_conns()

        for eq in self:
            for term in eq.terms:
                term.time_update(ts)

    def setup_initial_conditions(self, ics, regions, functions):
        self.variables.setup_initial_conditions(ics, regions, functions)

    def create_matrix_graph( self, var_names = None, vvar_names = None ):
        """
        Create tangent matrix graph. Order of dof connectivities is not
        important here...
        """
        if not self.variables.has_virtuals():
            output('no matrix (no test variables)!')
            return None

        shape = self.variables.get_matrix_shape()
        output( 'matrix shape:', shape )
        if nm.prod( shape ) == 0:
            output( 'no matrix (zero size)!' )
            return None

        adcs = self.variables.adof_conns

        # Only volume dof connectivities are used, with the exception of trace
        # surface dof connectivities.
        shared = set()
        rdcs = []
        cdcs = []
        for key, ii, info in iter_dict_of_lists(self.conn_info,
                                                return_keys=True):
            dct = info.dc_type.type
            if not (dct in ('volume', 'scalar') or info.is_trace):
                continue

            rvar, cvar = info.virtual, info.state
            if (rvar is None) or (cvar is None):
                continue

            rreg_name = info.get_region_name(can_trace=False)
            creg_name = info.get_region_name()

            for rig, cig in info.iter_igs():
                rname = rvar.get_primary_name()
                rkey = (rname, rreg_name, dct, rig)
                ckey = (cvar.name, creg_name, dct, cig)

                dc_key = (rkey, ckey)
                ## print dc_key

                if not dc_key in shared:
                    try:
                        rdcs.append(adcs[rkey])
                        cdcs.append(adcs[ckey])
                    except:
                        debug()
                    shared.add(dc_key)

        ## print shared
        for ii in range(len(rdcs)):
            if (rdcs[ii].ndim == 1) and (cdcs[ii].ndim == 2):
                rdcs[ii] = _fix_scalar_dc(rdcs[ii], cdcs[ii])

            elif (cdcs[ii].ndim == 1) and (rdcs[ii].ndim == 2):
                cdcs[ii] = _fix_scalar_dc(cdcs[ii], rdcs[ii])

            elif (cdcs[ii].ndim == 1) and (rdcs[ii].ndim == 1):
                rdcs[ii] = nm.array(rdcs[ii], ndmin=2)
                cdcs[ii] = nm.array(cdcs[ii], ndmin=2)

        ##     print rdcs[ii], cdcs[ii]
        ## pause()

        if not shared:
            # No virtual, state variable -> no matrix.
            output( 'no matrix (empty dof connectivities)!' )
            return None

        output( 'assembling matrix graph...' )
        tt = time.clock()

        ret, prow, icol = raw_graph( int( shape[0] ), int( shape[1] ),
                                    len( rdcs ), rdcs, cdcs )
        output( '...done in %.2f s' % (time.clock() - tt) )
        nnz = prow[-1]
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
        for cache in self.caches.itervalues():
            cache.init_time( ts )

    ##
    # 08.06.2007, c
    def advance( self, ts ):
        for cache in self.caches.itervalues():
            cache.advance( ts.step + 1 )

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

    def strip_state_vector(self, vec, follow_epbc=True):
        """
        Strip a full vector by removing EBC dofs. If 'follow_epbc' is True,
        values of EPBC master dofs are not simply thrown away, but added to the
        corresponding slave dofs, just like when assembling.
        """
        return self.variables.strip_state_vector(vec, follow_epbc=follow_epbc)

    def make_full_vec(self, svec, var_name=None, force_value=None):
        """
        Make a full vector satisfying E(P)BC
        from a stripped vector. For a selected variable if var_name is set.
        """
        return self.variables.make_full_vec(svec, var_name=var_name,
                                            force_value=force_value)

    def set_variables_from_state(self, vec, step=0):
        """
        Set data (vectors of DOF values) of variables.

        Paramters
        ---------
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

##
# 21.07.2006, c
class Equation( Struct ):

    @staticmethod
    def from_desc(name, desc, variables, regions, materials,
                  caches=None, user=None, term_prefixes=None):
        if term_prefixes is None: term_prefixes = {}

        itps = invert_dict(term_prefixes, True)

        term_descs = parse_definition(desc, itps)
        terms = Terms.from_desc(term_descs, regions)

        terms.setup()
        terms.assign_args(variables, materials, user)

        obj = Equation(name, terms, caches=caches)

        return obj

    def __init__(self, name, terms, caches=None):
        obj = Struct.__init__(self, name = name, terms = terms)

        self.terms.setup()

        caches = get_default(caches, DataCaches)
        self.terms.assign_caches(caches)

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

    def describe_geometry(self, geometries, integrals):
        for term in self.terms:
            term.describe_geometry(geometries, integrals)
