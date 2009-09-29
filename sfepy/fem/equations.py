from sfepy.base.base import *
from parseEq import create_bnf
from materials import Materials
from sfepy.terms import Terms, Term, term_table, DataCaches, cache_table

"""
Note:
- create found materials, variables from configuration/input file data
  ... no - user should be able to create objects even if they are not
  used in equations
"""

def parse_terms( regions, desc, itps ):
    """
    Parse equation given by 'desc' into terms. Assign each term its region.
    """
    # Parse.
    term_descs = []
    bnf = create_bnf( term_descs, itps )
    try:
        bnf.parseString( desc )
    except:
        print 'cannot parse:\n', desc
        raise
    
    # Construct terms.
    terms = OneTypeList( Term )
    for td in term_descs:
##         print td
##         pause()
        try:
            constructor = term_table[td.name]
        except:
            msg = "term '%s' is not in %s" % (td.name,
                                              sorted( term_table.keys() ))
            raise ValueError( msg )
        region = regions[td.region]
        arg_names = []
        arg_steps = {}
        arg_derivatives = {}
        arg_traces = {}
        for arg in td.args:
            trace = False
            derivative = None

            if isinstance(arg[1], int):
                name, step = arg

            else:
                kind = arg[0]
                name, step = arg[1]
                if kind == 'd':
                    derivative = arg[2]
                elif kind == 'tr':
                    trace = True

            arg_names.append( name )
            arg_steps[name] = step
            arg_derivatives[name] = derivative
            arg_traces[name] = trace

        term = constructor( region, td.name, td.sign )
        term.arg_names = arg_names
        term.arg_steps = arg_steps
        term.arg_derivatives = arg_derivatives
        term.arg_traces = arg_traces
        term.integral_name = td.integral

        terms.append( term )

    return terms

def setup_term_args( terms, variables, materials, user = None ):
    """terms ... can be both Terms or Term class
       - checks term argument existence in variables, materials, user
       - checks equality of field and term subdomain lists (igs)"""
    terms.classify_args( variables )
    for term in terms:
        igs = term.char_fun.igs
        vns = term.get_variable_names()
        for name in vns:
            if name not in variables.names:
                msg = 'variable "%s" not found' % name
                raise IndexError(msg)

            field = variables[name].field

            if term.arg_traces[name]:
                if not nm.all(nm.setmember1d(term.region.all_vertices,
                                             field.region.all_vertices)):
                    msg = ('%s: incompatible regions: (term, trace of field %s)'
                           + '(%s in %s)') %\
                           (term.name, field.name,
                            term.region.all_vertices, field.region.all_vertices)
                    raise ValueError(msg)
            else:
                if not set( igs ).issubset( set( field.aps.igs ) ):
                    msg = ('%s: incompatible regions: (term, field)'
                           + ' (%s(%s) in %s(%s)') %\
                             (term.name, igs, name, field.igs(), field.name)
                    raise ValueError(msg)

        mns = term.get_material_names()
        for name in mns:
            if name not in materials.names:
                msg = 'material "%s" not found' % name
                raise IndexError(msg)

            mat = materials[name]

            if not set( igs ).issubset( set( mat.igs ) ):
                msg= ('%s: incompatible regions: (term, material)'
                      + ' (%s(%s) in %s(%s)') %\
                      (term.name, igs, name, mat.igs, mat.name)
                raise ValueError(msg)

    if user is None:
        return

    uns = terms.get_user_names()
    uks = user.keys()
    for name in uns:
        if name not in uks:
            output( 'user data "%s" not found' % name )
            raise IndexError
#        print '********* ok'

##
# 24.07.2006, c
def build_args( term, variables, materials, **kwargs ):
    args = kwargs
    vns = term.get_variable_names()
    for vn in vns:
        args[vn] = variables[vn]
    mns = term.get_material_names()
    for mn in mns:
        args[mn] = materials[mn]
    return args

class ConnInfo(Struct):
    mirror_map = {}

    def get_region(self, can_trace=True):
        if self.is_trace and can_trace:
            return self.mirror_region
        else:
            return self.region

    def iter_igs(self):
        for ig in self.region.igs:
            ir = self.virtual_igs.index(ig)
            if not self.is_trace:
                ic = self.state_igs.index(ig)
            else:
                ic = self.state_igs.index(self.ig_map_i[ig])
            yield self.virtual_igs[ir], self.state_igs[ic]

##
# 21.07.2006, c
class Equations( Container ):

    ##
    # c: 18.04.2006, r: 20.02.2008
    def from_conf( conf ):
        objs = OneTypeList( Equation )

        conf = copy( conf )
        tps = conf.pop( 'namespaces', {} )
        itps = invert_dict( tps, True )

        ii = 0
        for name, desc in conf.iteritems():
            output( 'equation "%s":' %  name )
            output( desc )
            eq = Equation( name = name,
                           desc = desc,
                           itps = itps )
            objs.append( eq )
            ii += 1

        obj = Equations( objs, itps = itps )
        
        return obj
    from_conf = staticmethod( from_conf )

    def setup_terms( self, regions, variables, materials, caches = None,
                     user = None ):
        """Parse equations and create term instances.

        Grabs references to materials and variables."""
        if caches is None:
            self.caches = DataCaches()
        else:
            self.caches = caches

        self.materials = materials
        self.variables = variables

        conn_info = {}

        for eq in self:
            eq.setup_terms( regions, variables, materials, self.caches, user )
            eq.collect_conn_info(conn_info, variables)

##         print_structs(conn_info)
##         pause()
        
        variables.conn_info = conn_info
        
    ##
    # c: ??, r: 26.02.2008
    def describe_geometry( self, geometries, variables, integrals ):
        output( 'describing geometries...' )
        tt = time.clock()
        for eq in self:
            eq.describe_geometry( geometries, variables, integrals )
        output( '...done in %.2f s' % (time.clock() - tt) )
        

    ##
    # 24.08.2006, c
    # 24.04.2007
    def get_term_geometries( self ):
        tgs = set()
        for eq in self:
            tgs.update( eq.get_term_geometries() )
        return tgs

    ##
    # 16.11.2007, c
    def get_term_integral_names( self ):
        i_names = set()
        for eq in self:
            i_names.update( eq.get_term_integral_names() )
        return i_names

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

    ##
    # 25.07.2006, c
    # 28.08.2006
    # 12.02.2007
    def from_desc( name, desc, term_prefixes = None ):
        if term_prefixes is None: term_prefixes = {}

        obj = Equation( name = name, desc = desc,
                        itps = invert_dict( term_prefixes, True ) )
        return obj
    from_desc = staticmethod( from_desc )

    ##
    # 21.07.2006, c
    # 25.07.2006
    # 01.08.2006
    # 11.08.2006
    # 12.02.2007
    # 27.02.2007
    def parse_terms( self, regions ):
        terms = parse_terms( regions, self.desc, self.itps )
        self.terms = Terms( terms )


    ##
    # 21.07.2006, c
    # 24.07.2006
    # 22.08.2006
    # 25.08.2006
    # 27.11.2006
    # 20.02.2007
    def setup_term_args( self, variables, materials, user = None ):
        """- checks term argument existence in variables, materials, user
           - checks compatability of field and term subdomain lists (igs)"""
        setup_term_args( self.terms, variables, materials, user )

    ##
    # 29.11.2006, c
    # 27.02.2007
    # 08.06.2007
    # 11.06.2007
    def assign_term_caches( self, caches ):
        """
        History sizes for a particular cache instance are taken as maximum
        of history_sizes requirements of all terms using the instance.
        """
        for term in self.terms:
            if not hasattr( term, 'use_caches' ): continue

##             print term.name
            for name, arg_lists in term.use_caches.iteritems():
##                 print term.arg_names
##                print name, arg_lists
                for args in arg_lists:
                    # Order should be handled in terms...
                    args = copy( args )
                    if len(args) and (type( args[-1] ) == dict):
                        history_sizes = args.pop()
                    else:
                        history_sizes = None
                    ans = [term.get_arg_name( arg, full = True ) for arg in args]
                    cname = '_'.join( [name] + ans )
##                     print term.name, name, arg_lists, args, self.name, cname
##                     print history_sizes
##                     debug()
                    if caches.has_key( cname ):
                        caches[cname].merge_history_sizes( history_sizes )

                    else:
##                    print 'new'
                        try:
                            constructor = cache_table[name]
                        except:
                            raise RuntimeError, 'cache not found! %s in %s'\
                                  % (name, sorted( cache_table.keys() ))
                        cache = constructor( cname, ans, history_sizes )
                        caches.insert_cache( cache )
                caches.insert_term( cname, term.name, ans )
            term.caches = caches

    def setup_terms( self, regions, variables, materials, caches,
                     user = None ):
        """Parse equation and create term instances."""
        self.parse_terms( regions )
        self.setup_term_args( variables, materials, user )
        self.assign_term_caches( caches )

    def collect_conn_info(self, conn_info, variables):

        for term in self.terms:
            key = (self.name, term.name, term.integral_name, term.region.name)
            key += tuple(term.arg_names)

            vn = term.get_virtual_name()
            sns = term.get_state_names()
            pns = term.get_parameter_names()

            term_var_names = term.get_variable_names()
            pn_map = variables.get_primary_names(term_var_names)

##             print key
##             print vn, sns, pns
##             pause()
                
            dc_type = term.get_dof_conn_type()
            tgs = term.get_geometry()
            
            if vn is not None:
                v_igs = variables[vn].field.igs()
                v_tg = tgs[vn]
            else:
                v_igs = v_tg = None

            region = term.region

            is_any_trace = reduce(lambda x, y: x or y, term.arg_traces.values())
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
                
                s_igs = variables[sn].field.igs()
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
                ConnInfo.mirror_map[region.name] = (mirror_region,
                                                    ig_map, ig_map_i)
                vals.append(val)

            pns += aux_pns
            for pn in pns:
                p_igs = variables[pn].field.igs()
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
    ##
    #
    def describe_geometry( self, geometries, variables, integrals ):
        for term in self.terms:
            term.describe_geometry( geometries, variables, integrals )

    ##
    # 24.04.2007, c
    def get_term_geometries( self ):
        tgs = set()
        for term in self.terms:
            for tg in term.get_geometry():
                tgs.add( tg )
        return tgs

    ##
    # 16.11.2007, c
    def get_term_integral_names( self ):
        i_names = set()
        for term in self.terms:
            i_names.add( term.integral_name )
        return i_names
