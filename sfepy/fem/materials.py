from sfepy.base.base import *
from functions import ConstantFunction


##
# 21.07.2006, c
class Materials( Container ):

    def from_conf(conf, functions, wanted=None):
        """Construct Materials instance from configuration."""
        if wanted is None:
            wanted = conf.keys()

        objs = OneTypeList(Material)
        for key, mc in conf.iteritems():
            if key not in wanted: continue

            fun = get_default_attr(mc, 'function', None)
            vals = get_default_attr(mc, 'values', None)
            if (fun is not None) and (vals is not None):
                msg = 'material can have function or values but not both! (%s)' \
                      % mc
                raise ValueError(msg)
            elif vals is not None: # => fun is None
                fun = ConstantFunction(vals, functions = functions)
            else: # => vals is None
                fun = functions[fun]

            if (fun is None):
                msg = 'material has no values! (%s)' % mc
                raise ValueError(msg)

            kind = get_default_attr(mc, 'kind', 'time-dependent')
            flags = get_default_attr(mc, 'flags', {})
            mat =  Material(mc.name, mc.region, kind, fun, flags)
            objs.append(mat)

        obj = Materials( objs )
        return obj
    from_conf = staticmethod( from_conf )

    def semideep_copy(self):
        """Copy materials, while external data (e.g. region) remain shared."""
        others = copy(self)
        others.update(OneTypeList(Material))
        for mat in self:
            other = mat.copy(name=mat.name)
            other.reset()
            others.append(other)
        return others

    ##
    # 22.08.2006, c
    def setup_regions( self, regions ):
        for mat in self:
            mat.setup_regions( regions )

    def reset(self):
        """Clear material data so that next materials.time_update() is
        performed even for stationary materials."""
        for mat in self:
            mat.reset()

    def time_update(self, ts, domain, equations, variables, verbose=True):
        """Update material parameters for given time, domain, and equations."""
        if verbose: output('updating materials...')
        tt = time.clock()
        for mat in self:
            if verbose: output(' ', mat.name)
            mat.time_update(ts, domain, equations, variables)
        if verbose: output('...done in %.2f s' % (time.clock() - tt))

##
# 21.07.2006, c
class Material( Struct ):
    """
    A class holding constitutive and other material parameters.

    Example input:

    material_2 = {
       'name' : 'm',
       'region' : 'Omega',
       'values' : {'E' : 1.0},
    }

    Material parameters are passed to terms using the dot notation,
    i.e. 'm.E' in our example case.
    
    """
    def __init__(self, name, region_name, kind, function, flags):
        Struct.__init__(self, name = name,
                        region_name = region_name,
                        kind = kind,
                        function = function,
                        flags = flags,
                        region = None)
        self.reset()
    ##
    # 22.08.2006, c
    def setup_regions( self, regions ):
        region = regions[self.region_name]
        self.igs = region.igs
        self.region = region 


    def iter_qps(self, equations, variables, only_new=True):
        """
        Iterate groups of quadrature points, where the material data should be
        evaluated.
        """
        # Quadrature point function values.
        for equation in equations:
            for term in equation.terms:
                names = [ii.split('.')[0] for ii in term.names.material]
                if self.name not in names: continue

                key = (term.region.name, term.integral.name)
                if only_new and (key in self.datas): continue

                # Any term has at least one variable, all variables used
                # in a term share the same integral.
                var_name = term.names.variable[0]
                var = variables[var_name]

                aps = var.field.aps

                qps = aps.get_physical_qps(term.region, term.integral)

                yield key, term.igs(), qps

    def set_data(self, key, ig, qps, data):
        """
        Set the material data in quadrature points.

        Parameters
        ----------
        key : tuple
            The (region_name, integral_name) data key.
        ig : int
            The element group id.
        qps : Struct
            Information about the quadrature points.
        data : dict.
            The material data. Changes the shape of data!
        """
        datas = self.datas[key]

        # Restore shape to (n_el, n_qp, ...) until the C
        # core is rewritten to work with a bunch of physical
        # point values only.
        if qps.is_uniform:
            if data is not None:
                n_qp = qps.el_indx[ig][1] - qps.el_indx[ig][0]
                for val in data.itervalues():
                    val.shape = (val.shape[0] / n_qp, n_qp,
                                 val.shape[1], val.shape[2])
        else:
            raise NotImplementedError

        datas[ig] = data

    def time_update(self, ts, domain, equations, variables):
        """All material parameters are evaluated in all physical QPs."""
        self.data = None
        if self.datas and (self.kind == 'stationary'): return

        self.datas = {}
        # Quadrature point function values.
        for key, igs, qps in self.iter_qps(equations, variables):
            for ig in self.igs:
                if ig not in igs: continue

                self.datas.setdefault(key, {})
                if (qps.n_qp[ig] == 0):
                    self.set_data(key, ig, qps, None)
                    continue

                data = self.function(ts, qps.values[ig], mode='qp',
                                     region=self.region, ig=ig,
                                     **self.extra_args)

                self.set_data(key, ig, qps, data)

        # Special function values (e.g. flags).
        datas = self.function(ts, domain.get_mesh_coors(), mode='special',
                              region=self.region)
        if datas is not None:
            self.datas['special'] = datas
            self.special_names.update(datas.keys())

        # Special constant values.
        if self.flags.get('special_constant'):
            datas = self.function(None, None)
            self.datas['special_constant'] = datas
            self.constant_names.update(datas.keys())

    def set_all_data( self, datas ):
        self.mode = 'user'
        self.datas = datas
        self.data = None

    def set_function(self, function):
        self.function = function
        self.reset()

    def reset(self):
        """Clear all data created by a call to time_update()."""
        self.datas = None
        self.special_names = set()
        self.constant_names = set()
        self.data = None
        self.extra_args = {}

    ##
    # 01.08.2007, c
    def set_extra_args(self, **extra_args):
        """Extra arguments passed tu the material function."""
        self.extra_args = extra_args
        
    def get_data( self, key, ig, name ):
        """`name` can be a dict - then a Struct instance with data as
        attributes named as the dict keys is returned."""
##         print 'getting', self.name, name

        if isinstance( name, str ):
            return self._get_data( key, ig, name )
        else:
            out = Struct()
            for key, item in name.iteritems():
                setattr( out, key, self._get_data( key, ig, item ) )
            return out
                       
    def _get_data( self, key, ig, name ):
        if name is None:
            msg = 'material arguments must use the dot notation!\n'\
                  '(material: %s, region: %s)' % (self.name, region_name)
            raise ValueError( msg )

        if self.datas is None:
            raise ValueError( 'material data not set! (call time_update())' )

        if name in self.special_names:
            # key, ig ignored.
            return self.datas['special'][name]

        else:
            datas = self.datas[key]

            if isinstance( datas[ig], Struct ):
                return getattr( datas[ig], name )
            elif datas[ig] is not None:
                return datas[ig][name]

    def get_constant_data(self, name):
        """Get constant data by name."""
        if name in self.constant_names:
            # no key, ig.
            return self.datas['special_constant'][name]
        else:
            raise ValueError('material %s has no constant %s!'
                             % (self.name, name))

    ##
    # 01.08.2007, c
    def reduce_on_datas( self, reduce_fun, init = 0.0 ):
        """For non-special values only!"""
        out = {}.fromkeys(self.datas[self.datas.keys()[0]][0].keys(), init)
        for datas in self.datas.itervalues():
            for data in datas.itervalues():
                for key, val in data.iteritems():
                    out[key] = reduce_fun(out[key], val)

        return out
