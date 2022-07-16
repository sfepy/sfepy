from __future__ import absolute_import

from sfepy.base.base import (Struct, Container, OneTypeList, assert_,
                             output, get_default, basestr)
from sfepy.base.timing import Timer
from .functions import ConstantFunction, ConstantFunctionByRegion
import six
import numpy as nm

class Materials(Container):

    @staticmethod
    def from_conf(conf, functions, wanted=None):
        """
        Construct Materials instance from configuration.
        """
        if wanted is None:
            wanted = list(conf.keys())

        objs = OneTypeList(Material)
        for key, mc in six.iteritems(conf):
            if key not in wanted: continue

            mat = Material.from_conf(mc, functions)
            objs.append(mat)

        obj = Materials(objs)
        return obj

    def reset(self):
        """Clear material data so that next materials.time_update() is
        performed even for stationary materials."""
        for mat in self:
            mat.reset()

    def time_update(self, ts, equations, mode='normal', problem=None,
                    verbose=True):
        """
        Update material parameters for given time, problem, and equations.

        Parameters
        ----------
        ts : TimeStepper instance
            The time stepper.
        equations : Equations instance
            The equations using the materials.
        mode : 'normal', 'update' or 'force'
            The update mode, see :func:`Material.time_update()`.
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        verbose : bool
            If False, reduce verbosity.
        """
        if verbose: output('updating materials...')
        timer = Timer(start=True)
        for mat in self:
            if verbose: output(' ', mat.name)
            mat.time_update(ts, equations, mode=mode, problem=problem)
        if verbose: output('...done in %.2f s' % timer.stop())

class Material(Struct):
    """
    A class holding constitutive and other material parameters.

    Example input::

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
        kind = conf.get('kind', 'time-dependent')
        flags = conf.get('flags', {})

        function = conf.get('function', None)
        values = conf.get('values', None)

        if isinstance(function, basestr):
            function = functions[function]

        obj = Material(conf.name, kind, function, values, flags)

        return obj

    def __init__(self, name, kind='time-dependent',
                 function=None, values=None, flags=None, **kwargs):
        """
        A material is defined either by a function, or by a set of constant
        values, potentially distinct per region. Therefore, either `function`
        must be specified, or a combination of `values` and `**kwargs`.

        For constant materials, `**kwargs` are simply combined with `values`
        into a dictionary mapping material parameter names to parameter values.
        The parameter values may either be specified as a constant value, or as
        another dictionary mapping region names to constant values (see
        :py:class:`sfepy.discrete.functions.ConstantFunctionByRegion`).

        Special material parameters, that are not evaluated in quadrature
        points - for example flags or geometry independent data - are denoted
        by parameter names starting with '.' - in this case the `values`
        argument need to be used, or a function that returns the parameters
        when ``mode == 'special'``.

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
        Struct.__init__(self, name=name, kind=kind, is_constant=False)

        if kwargs:
            if values is None:
                values = kwargs
            else:
                values = dict(values)
                values.update(kwargs)

        if (function is not None) and (values is not None):
            raise ValueError(
                f'material {self.name}: use either "function" or "values"'
                ' arguments but not both!'
            )

        if (function is None) and (values is None):
            raise ValueError(
                f'material {self.name}: neither "function" nor "values"'
                ' arguments (or keyword arguments) given!'
            )

        self.flags = get_default(flags, {})

        if function is not None:
            if not hasattr(function, '__call__'):
                raise TypeError(
                    f'material {self.name}: "function" needs to be callable!'
                )
            self.function = function

        else: # => function is None
            assert_(all(isinstance(k, str) for k in values.keys()))
            isbyregion = list(
                (not k.startswith('.')) and isinstance(v, dict)
                for k,v in values.items()
            )

            if all(isbyregion):
                self.function = ConstantFunctionByRegion(values)
                self.is_constant = True

            elif not any(isbyregion):
                self.function = ConstantFunction(values, no_tile=True)
                self.is_constant = True

            else:
                raise ValueError(
                    f'material {self.name}: Either all parameter values need to '
                    'be specified by region, or none at all.'
                )

        self.reset()

    def iter_terms(self, equations, only_new=True):
        """
        Iterate terms for which the material data should be evaluated.
        """
        if equations is None: return

        for equation in equations:
            for term in equation.terms:
                names = [ii[0] for ii in term.names.material]
                if self.name not in names: continue

                key = term.get_qp_key()
                if only_new and (key in self.datas): continue

                self.datas.setdefault(key, {})

                yield key, term

    def set_data(self, key, qps, data):
        """
        Set the material data in quadrature points.

        Parameters
        ----------
        key : tuple
            The (region_name, integral_name) data key.
        qps : Struct
            Information about the quadrature points.
        data : dict
            The material data.
        """
        # Restore shape to (n_el, n_qp, ...) until the C
        # core is rewritten to work with a bunch of physical
        # point values only.
        new_data = {}
        if data is not None:
            for dkey, val in six.iteritems(data):
                if val.ndim != 3:
                    raise ValueError('material parameter array must have'
                                     " three dimensions! ('%s' has %d)"
                                     % (dkey, val.ndim))
                qps_shape = qps.get_shape(val.shape)
                if qps_shape[0] == 0:
                    new_data[dkey] = nm.tile(val, (1, qps_shape[1], 1, 1))
                else:
                    new_data[dkey] = val.reshape(qps_shape)

        self.datas[key] = new_data

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
        problem : Problem, optional
            The problem definition for which the update occurs.
        """
        self.datas.setdefault(key, {})

        qps = term.get_physical_qps()
        coors = qps.values
        data = self.function(ts, coors, mode='qp',
                             equations=equations, term=term, problem=problem,
                             **self.extra_args)

        self.set_data(key, qps, data)

    def update_special_data(self, ts, equations, problem=None):
        """
        Update the special material parameters.

        Parameters
        ----------
        ts : TimeStepper
            The time stepper.
        equations : Equations
            The equations for which the update occurs.
        problem : Problem, optional
            The problem definition for which the update occurs.
        """
        if 'special' in self.datas: return

        # Special function values (e.g. flags).
        datas = self.function(ts, None, mode='special',
                              problem=problem, equations=equations,
                              **self.extra_args)
        if datas is not None:
            self.datas['special'] = datas
            self.special_names.update(list(datas.keys()))

    def update_special_constant_data(self, equations=None, problem=None):
        """
        Update the special constant material parameters.

        Parameters
        ----------
        equations : Equations
            The equations for which the update occurs.
        problem : Problem, optional
            The problem definition for which the update occurs.
        """
        if 'special_constant' in self.datas: return
        if not self.flags.get('special_constant'): return

        # Special constant values.
        datas = self.function(None, None, mode='special_constant',
                              problem=problem, equations=equations)
        self.datas['special_constant'] = datas
        self.constant_names.update(list(datas.keys()))

    def time_update(self, ts, equations, mode='normal', problem=None):
        """
        Evaluate material parameters in physical quadrature points.

        Parameters
        ----------
        ts : TimeStepper instance
            The time stepper.
        equations : Equations instance
            The equations using the materials.
        mode : 'normal', 'update' or 'force'
            The update mode. In 'force' mode, ``self.datas`` is cleared and all
            updates are redone. In 'update' mode, existing data are preserved
            and new can be added. The 'normal' mode depends on other
            attributes: for stationary (``self.kind == 'stationary'``)
            materials and materials in 'user' mode, nothing is done if
            ``self.datas`` is not empty. For time-dependent materials
            (``self.kind == 'time-dependent'``, the default) that are not
            constant, i.e., are given by a user function, 'normal' mode behaves
            like 'force' mode. For constant materials it behaves like 'update'
            mode - existing data are reused.
        problem : Problem instance, optional
            The problem that can be passed to user functions as a context.
        """
        if mode == 'force':
            self.datas = {}

        elif self.datas:
            if mode == 'normal':
                if (self.mode == 'user') or (self.kind == 'stationary'):
                    return

                elif not self.is_constant:
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
            keys = list(self.datas.keys())

        else:
            keys = [key for key in self.datas.keys()
                    if (isinstance(key, tuple) and key[0] == region_name)]

        return keys

    def set_all_data(self, datas):
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
        Clear all data created by a call to ``time_update()``, set
        ``self.mode`` to ``None``.
        """
        self.mode = None
        self.datas = {}
        self.special_names = set()
        self.constant_names = set()
        self.extra_args = {}

    def set_extra_args(self, **extra_args):
        """Extra arguments passed tu the material function."""
        self.extra_args = extra_args

    def get_data(self, key, name):
        """`name` can be a dict - then a Struct instance with data as
        attributes named as the dict keys is returned."""

        if isinstance(name, basestr):
            return self._get_data(key, name)
        else:
            out = Struct()
            for key, item in six.iteritems(name):
                setattr(out, key, self._get_data(key, item))
            return out

    def _get_data(self, key, name):
        if name is None:
            msg = 'material arguments must use the dot notation!\n'\
                  '(material: %s, key: %s)' % (self.name, key)
            raise ValueError(msg)

        if not self.datas:
            raise ValueError('material data not set! (call time_update())')

        if name in self.special_names:
            # key ignored.
            return self.datas['special'][name]

        else:
            datas = self.datas[key]

            if isinstance(datas, Struct):
                return getattr(datas, name)

            elif datas:
                return datas[name]

    def get_constant_data(self, name):
        """Get constant data by name."""
        if name in self.constant_names:
            # no key.
            return self.datas['special_constant'][name]
        else:
            raise ValueError('material %s has no constant %s!'
                             % (self.name, name))

    def reduce_on_datas(self, reduce_fun, init=0.0):
        """For non-special values only!"""
        out = {}.fromkeys(list(self.datas[list(self.datas.keys())[0]].keys()), init)
        for data in six.itervalues(self.datas):
            for key, val in six.iteritems(data):
                out[key] = reduce_fun(out[key], val)

        return out
