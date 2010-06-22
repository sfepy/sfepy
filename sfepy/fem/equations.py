from sfepy.base.base import *
from sfepy.fem import Materials, Variables
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

class ConnInfo(Struct):
    mirror_map = {}

    def get_region(self, can_trace=True):
        if self.is_trace and can_trace:
            return self.mirror_region
        else:
            return self.region

    def get_region_name(self, can_trace=True):
        if self.is_trace and can_trace:
            reg = self.mirror_region
        else:
            reg = self.region

        if reg is not None:            
            return reg.name
        else:
            return None

    def iter_igs(self):
        if self.region is not None:
            for ig in self.region.igs:
                if self.virtual_igs is not None:
                    ir = self.virtual_igs.index(ig)
                    rig = self.virtual_igs[ir]
                else:
                    rig = None

                if not self.is_trace:
                    ii = ig
                else:
                    ii = self.ig_map_i[ig]

                if self.state_igs is not None:
                    ic = self.state_igs.index(ii)
                    cig = self.state_igs[ic]
                else:
                    cig = None
                    
                yield rig, cig

        else:
            yield None, None

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
            eq.collect_conn_info(self.conn_info, self.variables)

        ## print_structs(self.conn_info)
        ## pause()

        return self.conn_info

    def describe_geometry(self, integrals):
        output('describing geometries...')
        tt = time.clock()
        for eq in self:
            eq.describe_geometry(self.geometries, self.variables, integrals)
        output('...done in %.2f s' % (time.clock() - tt))
        
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

    def time_update( self, ts ):
        for eq in self:
            for term in eq.terms:
                term.time_update( ts )

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
                variables.append(var.get_primary())

        return variables

    def collect_conn_info(self, conn_info, variables):

        for term in self.terms:
            key = (self.name,) + term.get_conn_key()

            vn = term.get_virtual_name()
            sns = term.get_state_names()
            pns = term.get_parameter_names()

            term_var_names = term.get_variable_names()
            pn_map = variables.get_primary_names(term_var_names)

##             print key
##             print vn, sns, pns
##             pause()
                
            dc_type = term.get_dof_conn_type()
            tgs = term.get_geometry_types()

            v_igs = v_tg = None
            if vn is not None:
                field = variables[vn].get_field()
                if field is not None:
                    v_igs = field.igs()
                    v_tg = tgs[vn]

            region = term.get_region()
            if region is not None:
                is_any_trace = reduce(lambda x, y: x or y,
                                      term.arg_traces.values())
                if is_any_trace:
                    for reg in region.domain.regions:
                        if (reg is not region) and \
                               (len(reg.igs) == len(region.igs)) and \
                               nm.all(region.all_vertices == reg.all_vertices):
                            mirror_region = reg
                            break
                    else:
                        raise ValueError('trace: cannot find mirror region! (%s)' \
                                         % region)

                    ig_map = {}
                    ig_map_i = {}
                    for igr in region.igs:
                        for igc in mirror_region.igs:
                            if nm.all(region.vertices[igr] ==
                                      mirror_region.vertices[igc]):
                                ig_map[igc] = igr
                                ig_map_i[igr] = igc
                                break
                        else:
                            raise ValueError('trace: cannot find group! (%s)' \
                                             % geom_request)
                else:
                    mirror_region = ig_map = ig_map_i = None
                
            vals = []
            aux_pns = []
            for sn in sns:
                # Allow only true state variables.
                if not variables[sn].is_state():
                    aux_pns.append(sn)
                    continue
                
                
                field = variables[sn].get_field()
                if field is not None:
                    s_igs = field.igs()
                else:
                    s_igs = None
                is_trace = term.arg_traces[sn]

                if sn in tgs:
                    ps_tg = tgs[sn]
                else:
                    ps_tg = v_tg

                val = ConnInfo(virtual = vn, virtual_igs = v_igs,
                               state = sn, state_igs = s_igs,
                               primary = sn, primary_igs = s_igs,
                               has_virtual = True,
                               has_state = True,
                               is_trace = is_trace,
                               dc_type = dc_type,
                               v_tg = v_tg,
                               ps_tg = ps_tg,
                               region = region,
                               mirror_region = mirror_region,
                               all_vars = term_var_names,
                               ig_map = ig_map, ig_map_i = ig_map_i)
                if region is not None:
                    ConnInfo.mirror_map[region.name] = (mirror_region,
                                                        ig_map, ig_map_i)
                vals.append(val)

            pns += aux_pns
            for pn in pns:
                field = variables[pn].get_field()
                if field is not None:
                    p_igs = field.igs()
                else:
                    p_igs = None
                is_trace = term.arg_traces[pn]

                if pn in tgs:
                    ps_tg = tgs[pn]
                else:
                    ps_tg = v_tg

                val = ConnInfo(virtual = vn, virtual_igs = v_igs,
                               state = None, state_igs = [],
                               primary = pn_map[pn], primary_igs = p_igs,
                               has_virtual = vn is not None,
                               has_state = False,
                               is_trace = is_trace,
                               dc_type = dc_type,
                               v_tg = v_tg,
                               ps_tg = ps_tg,
                               region = region,
                               mirror_region = mirror_region,
                               all_vars = term_var_names,
                               ig_map = ig_map, ig_map_i = ig_map_i)
                if region is not None:
                    ConnInfo.mirror_map[region.name] = (mirror_region,
                                                        ig_map, ig_map_i)
                vals.append(val)

            if vn and (len(vals) == 0):
                # No state, parameter variables, just the virtual one.
                val = ConnInfo(virtual = vn, virtual_igs = v_igs,
                               state = pn_map[vn], state_igs = v_igs,
                               primary = pn_map[vn], primary_igs = v_igs,
                               has_virtual = True,
                               has_state = False,
                               is_trace = False,
                               dc_type = dc_type,
                               v_tg = v_tg,
                               ps_tg = v_tg,
                               region = region,
                               mirror_region = None,
                               all_vars = term_var_names,
                               ig_map = None, ig_map_i = None)
                vals.append(val)
            
            conn_info[key] = vals

    def describe_geometry(self, geometries, variables, integrals):
        for term in self.terms:
            term.describe_geometry(geometries, variables, integrals)
