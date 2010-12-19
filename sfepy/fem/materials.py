import time
from copy import copy

import numpy as nm

from sfepy.base.base import Struct, Container, OneTypeList
from sfepy.base.base import output, get_default_attr, get_default
from functions import ConstantFunction


##
# 21.07.2006, c
class Materials( Container ):

    @staticmethod
    def from_conf(conf, functions, wanted=None):
        """
        Construct Materials instance from configuration.
        """
        if wanted is None:
            wanted = conf.keys()

        objs = OneTypeList(Material)
        for key, mc in conf.iteritems():
            if key not in wanted: continue

            mat = Material.from_conf(mc, functions)
            objs.append(mat)

        obj = Materials( objs )
        return obj

    def semideep_copy(self):
        """Copy materials, while external data (e.g. region) remain shared."""
        others = copy(self)
        others.update(OneTypeList(Material))
        for mat in self:
            other = mat.copy(name=mat.name)
            other.reset()
            others.append(other)
        return others

    def reset(self):
        """Clear material data so that next materials.time_update() is
        performed even for stationary materials."""
        for mat in self:
            mat.reset()

    def time_update(self, ts, equations, problem=None,
                    force=False, verbose=True):
        """
        Update material parameters for given time, problem, and equations.
        """
        if verbose: output('updating materials...')
        tt = time.clock()
        for mat in self:
            if verbose: output(' ', mat.name)
            mat.time_update(ts, equations, problem, force=force)
        if verbose: output('...done in %.2f s' % (time.clock() - tt))

##
# 21.07.2006, c
class Material( Struct ):
    """
    A class holding constitutive and other material parameters.

    Example input:

    material_2 = {
       'name' : 'm',
       'values' : {'E' : 1.0},
    }

    Material parameters are passed to terms using the dot notation,
    i.e. 'm.E' in our example case.
    """
    @staticmethod
    def from_conf(conf, functions):
        """
        Construct Material instance from configuration.
        """
        kind = get_default_attr(conf, 'kind', 'time-dependent')
        flags = get_default_attr(conf, 'flags', {})

        function = get_default_attr(conf, 'function', None)
        values = get_default_attr(conf, 'values', None)

        if isinstance(function, str):
            function = functions[function]

        obj =  Material(conf.name, kind, function, values, flags)

        return obj
    
    def __init__(self, name, kind='time-dependent',
                 function=None, values=None, flags=None, **kwargs):
        """
        Parameters
        ----------
        name : str
            The name of the material.
        kind : 'time-dependent' or 'stationary'
            The kind of the material.
        function : function
            The function for setting up the material values.
        values : dict
            Constant material values.
        flags : dict, optional
            Special flags.
        **kwargs : keyword arguments, optional
            Constant material values passed by their names.
        """
        Struct.__init__(self, name=name, kind=kind)

        if (function is not None) and ((values is not None) or len(kwargs)):
            msg = 'material can have function or values but not both! (%s)' \
                  % self.name
            raise ValueError(msg)

        self.flags = get_default(flags, {})

        if hasattr(function, '__call__'): 
            self.function = function

        elif (values is not None) or len(kwargs): # => function is None
                all_values = {}
                if values is not None:
                    all_values.update(values)
                all_values.update(kwargs)

                self.function = ConstantFunction(all_values)

        else: # => both values and function are None
            msg = 'material %s: neither function nor values given! (%s)' \
                  % self.name
            raise ValueError(msg)

        self.reset()

    def iter_terms(self, equations, only_new=True):
        """
        Iterate terms for which the material data should be evaluated.
        """
        if equations is None: raise StopIteration

        for equation in equations:
            for term in equation.terms:
                names = [ii[0] for ii in term.names.material]
                if self.name not in names: continue

                key = term.get_qp_key()
                if only_new and (key in self.datas): continue

                self.datas.setdefault(key, {})

                yield key, term

    def set_data(self, key, ig, qps, data, indx):
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
        data : dict
            The material data. Changes the shape of data!
        indx : array
            The indices of elements in the group `ig`.
        """
        datas = self.datas[key]

        # Restore shape to (n_el, n_qp, ...) until the C
        # core is rewritten to work with a bunch of physical
        # point values only.
        group_data = {}
        if qps.is_uniform:
            if data is not None:
                n_qp = qps.el_indx[ig][1] - qps.el_indx[ig][0]
                for key, val in data.iteritems():
                    aux = val[indx]
                    aux.shape = (aux.shape[0] / n_qp, n_qp,
                                 aux.shape[1], aux.shape[2])
                    group_data[key] = aux
        else:
            raise NotImplementedError

        datas[ig] = group_data

    def set_data_from_variable(self, var, name, equations):
        for key, term in self.iter_terms(equations):
            qps = term.get_physical_qps()
            for ig in term.igs():
                data = var.evaluate_at(qps.values[ig])
                data.shape = data.shape + (1,)

                self.set_data(key, ig, qps, {name : data})

    def update_data(self, key, ts, equations, term, problem=None):
        """
        Update the material parameters in quadrature points.

        Parameters
        ----------
        key : tuple
            The (region_name, integral_name) data key.
        ts : TimeStepper
            The time stepper.
        equations : Equations
            The equations for which the update occurs.
        term : Term
            The term for which the update occurs.
        problem : ProblemDefinition, optional
            The problem definition for which the update occurs.
        """
        self.datas.setdefault(key, {})

        qps = term.get_physical_qps()
        coors = nm.concatenate(qps.values.values(), axis=0)

        data = self.function(ts, coors, mode='qp',
                             equations=equations, term=term, problem=problem,
                             group_indx=qps.group_indx,
                             **self.extra_args)

        for ig, indx in qps.group_indx.iteritems():
            if (qps.n_qp[ig] == 0):
                self.set_data(key, ig, qps, None, None)

            else:
                self.set_data(key, ig, qps, data, indx)

    def update_special_data(self, ts, equations, problem=None):
        """
        Update the special material parameters.

        Parameters
        ----------
        ts : TimeStepper
            The time stepper.
        equations : Equations
            The equations for which the update occurs.
        problem : ProblemDefinition, optional
            The problem definition for which the update occurs.
        """
        # Special function values (e.g. flags).
        datas = self.function(ts, problem.get_mesh_coors(), mode='special',
                              problem=problem, equations=equations,
                              **self.extra_args)
        if datas is not None:
            self.datas['special'] = datas
            self.special_names.update(datas.keys())

    def update_special_constant_data(self, equations=None, problem=None):
        """
        Update the special constant material parameters.

        Parameters
        ----------
        equations : Equations
            The equations for which the update occurs.
        problem : ProblemDefinition, optional
            The problem definition for which the update occurs.
        """
        # Special constant values.
        if self.flags.get('special_constant'):
            datas = self.function(None, None, mode='special_constant',
                                  problem=problem, equations=equations)
            self.datas['special_constant'] = datas
            self.constant_names.update(datas.keys())

    def time_update(self, ts, equations, problem=None, force=False):
        """
        Evaluate material parameters in physical quadrature points.

        Setting `force` to True forces the update to be done. If `force`
        is False, do nothing, if ``self.mode == 'user'`` or ``self.kind
        == 'stationary'`` and the parameters are already set.
        """
        if not force and ((self.mode == 'user')
                          or self.datas and (self.kind == 'stationary')):
            return

        self.datas = {}
        for key, term in self.iter_terms(equations):
            self.update_data(key, ts, equations, term, problem=problem)

        self.update_special_data(ts, equations, problem=problem)
        self.update_special_constant_data(equations, problem=problem)

    def get_keys(self, region_name=None):
        """
        Get all data keys.

        Parameters
        ----------
        region_name : str
            If not None, only keys with this region are returned.
        """
        if not self.datas:
            keys = None

        elif region_name is None:
            keys = self.datas.keys()

        else:
            keys = [key for key in self.datas.keys()
                    if (isinstance(key, tuple) and key[0] == region_name)]

        return keys

    def set_all_data( self, datas ):
        """
        Use the provided data, set mode to 'user'.
        """
        self.mode = 'user'
        self.datas = datas

    def set_function(self, function):
        self.function = function
        self.reset()

    def reset(self):
        """
        Clear all data created by a call to ``time_update()``, set ``self.mode``
        to ``None``.
        """
        self.mode = None
        self.datas = {}
        self.special_names = set()
        self.constant_names = set()
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
                  '(material: %s, key: %s)' % (self.name, key)
            raise ValueError( msg )

        if not self.datas:
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
