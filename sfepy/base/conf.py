"""
Problem description file handling.

Notes
-----

Short syntax: key is suffixed with '__<number>' to prevent collisions with long
syntax keys -> both cases can be used in a single input.
"""
from __future__ import absolute_import
import re
import numpy as nm

from sfepy.base.base import (Struct, IndexedStruct, dict_to_struct,
                             output, copy, update_dict_recursively,
                             import_file, assert_, get_default, basestr)
from sfepy.base.parse_conf import create_bnf
import six

# TODO this requires EBC to be specified but in some examples only
#   EPBCs specified is valid
_required = ['filename_mesh|filename_domain', 'field_[0-9]+|fields',
             '(ebc_[0-9]+|ebcs|dgebc_[0-9]+|dgebcs)', 'equations',
             'region_[0-9]+|regions', 'variable_[0-9]+|variables',
             'material_[0-9]+|materials',
             'solver_[0-9]+|solvers']

_other = ['epbc_[0-9]+|epbcs', 'dgepbc_[0-9]+|dgepbcs',
          'lcbc_[0-9]+|lcbcs', 'nbc_[0-9]+|nbcs',
          'ic_[0-9]+|ics', 'function_[0-9]+|functions', 'options',
          'integral_[0-9]+|integrals']

def get_standard_keywords():
    return copy(_required), copy(_other)

def tuple_to_conf(name, vals, order):
    """
    Convert a configuration tuple `vals` into a Struct named `name`, with
    attribute names given in and ordered by `order`.

    Items in `order` at indices outside the length of `vals` are ignored.
    """
    conf = Struct(name=name)
    for ii, key in enumerate(order[:len(vals)]):
        setattr(conf, key, vals[ii])
    return conf

def transform_variables(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf, ['kind', 'field'])
            if len(conf) >= 3:
                kind = c2.kind.split()[0]
                if kind == 'unknown':
                    c2.order = conf[2]
                elif kind == 'test':
                    c2.dual = conf[2]
                elif kind == 'parameter':
                    if isinstance(conf[2], basestr) or (conf[2] is None):
                        c2.like = conf[2]
                    else:
                        c2.like = None
                        c2.special = conf[2]
                if len(conf) == 4:
                    c2.history = conf[3]
            d2['variable_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            d2['variable_'+c2.name] = c2
    return d2

def transform_conditions(adict, prefix):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            if len(conf) == 2:
                c2 = tuple_to_conf(key, conf, ['region', 'dofs'])

            else:
                c2 = tuple_to_conf(key, conf, ['region', 'times', 'dofs'])

            d2['%s_%s__%d' % (prefix, c2.name, ii)] = c2

        else:
            c2 = transform_to_struct_1(conf)
            d2['%s_%s' % (prefix, c2.name)] = c2

    return d2

def transform_dgebcs(adict):
    return transform_conditions(adict, "dgebc")

def transform_ebcs(adict):
    return transform_conditions(adict, 'ebc')

def transform_ics(adict):
    return transform_conditions(adict, 'ic')

def transform_lcbcs(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            if len(conf) >= 4:
                if isinstance(conf[1], dict):
                    c2 = tuple_to_conf(key, conf, ['region', 'dofs',
                                                   'dof_map_fun', 'kind'])
                    c2.arguments = conf[4:]

                else:
                    c2 = tuple_to_conf(key, conf, ['region', 'times', 'dofs',
                                                   'dof_map_fun', 'kind'])
                    c2.arguments = conf[5:]

            else:
                msg = 'LCBC syntax has to be: region[, times], dofs,' \
                      ' dof_map_fun, kind[, other arguments]'
                raise SyntaxError(msg)

            d2['lcbc_%s__%d' % (c2.name, ii)] = c2

        else:
            c2 = transform_to_struct_1(conf)
            c2.set_default('dof_map_fun', None)
            c2.set_default('arguments', ())
            d2['lcbc_%s' % (c2.name)] = c2

    return d2

def transform_epbcs(adict, prefix="epbc"):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            if len(conf) == 3:
                c2 = tuple_to_conf(key, conf, ['region', 'dofs', 'match'])

            else:
                c2 = tuple_to_conf(key, conf,
                                   ['region', 'times', 'dofs', 'match'])

            d2['%s_%s__%d' % (prefix, c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            d2['%s_%s' % (prefix, c2.name)] = c2
    return d2

def transform_dgepbcs(adict):
    return transform_epbcs(adict, "dgepbc")

def transform_regions(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, basestr):
            c2 = Struct(name=key, select=conf)
            d2['region_%s__%d' % (c2.name, ii)] = c2
        elif isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf, ['select', 'kind'])
            if len(conf) == 3:
                c2.parent = conf[2]
            if len(conf) == 4:
                c2.parent = conf[2]
                c2.extra_options = conf[3]
            d2['region_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            d2['region_'+c2.name] = c2
    return d2

def transform_integrals(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, int):
            c2 = Struct(name=key, order=conf)
            d2['integral_%s__%d' % (c2.name, ii)] = c2

        elif isinstance(conf, tuple):
            if len(conf) == 2: # Old tuple version with now-ignored 'kind'.
                conf = conf[1]
                c2 = Struct(name=key, order=conf)

            elif len(conf) == 3:
                c2 = tuple_to_conf(key, conf, ['order', 'vals', 'weights'])

            d2['integral_%s__%d' % (c2.name, ii)] = c2

        else:
            c2 = transform_to_struct_1(conf)
            d2['integral_'+c2.name] = c2

    return d2

def transform_fields(adict):
    dtypes = {'real' : nm.float64, 'complex' : nm.complex128}
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf,
                               ['dtype', 'shape', 'region', 'approx_order',
                                'space', 'poly_space_base'])
            if c2.dtype in dtypes:
                c2.dtype = dtypes[c2.dtype]
            d2['field_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            c2.set_default('dtype', nm.float64)
            if c2.dtype in dtypes:
                c2.dtype = dtypes[c2.dtype]
            d2['field_'+c2.name] = c2
    return d2

def transform_materials(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, basestr):
            c2 = Struct(name=key, function=conf)
            d2['material_%s__%d' % (c2.name, ii)] = c2

        elif isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf,
                               ['values', 'function', 'kind'])
            if len(conf) == 4:
                c2.flags = conf[3]
            d2['material_%s__%d' % (c2.name, ii)] = c2

        else:
            c2 = transform_to_struct_1(conf)
            d2['material_'+conf['name']] = c2

    return d2

def transform_solvers(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf, ['kind','params'])
            for param, val in six.iteritems(c2.params):
                setattr(c2, param, val)
            delattr(c2, 'params')
            d2['solvers_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            d2['solvers_'+c2.name] = c2
    return d2

def transform_functions(adict):
    d2 = {}
    for ii, (key, conf) in enumerate(six.iteritems(adict)):
        if isinstance(conf, tuple):
            c2 = tuple_to_conf(key, conf, ['function'])
            d2['function_%s__%d' % (c2.name, ii)] = c2
        else:
            c2 = transform_to_struct_1(conf)
            d2['function_'+c2.name] = c2
    return d2

def transform_to_struct_1(adict):
    return dict_to_struct(adict, flag=(1,))
def transform_to_i_struct_1(adict):
    return dict_to_struct(adict, flag=(1,), constructor=IndexedStruct)
def transform_to_struct_01(adict):
    return dict_to_struct(adict, flag=(0,1))
def transform_to_struct_10(adict):
    return dict_to_struct(adict, flag=(1,0))

transforms = {
    'options'   : transform_to_i_struct_1,
    'solvers'   : transform_solvers,
    'integrals' : transform_integrals,
    'regions'   : transform_regions,
    'fields'    : transform_fields,
    'variables' : transform_variables,
    'ebcs'      : transform_ebcs,
    'epbcs'     : transform_epbcs,
    'dgebcs'    : transform_dgebcs,
    'dgepbcs'   : transform_dgepbcs,
    'nbcs'      : transform_to_struct_01,
    'lcbcs'     : transform_lcbcs,
    'ics'       : transform_ics,
    'materials' : transform_materials,
    'functions' : transform_functions,
}

def dict_from_string(string, allow_tuple=False, free_word=False):
    """
    Parse `string` and return a dictionary that can be used to
    construct/override a ProblemConf instance.
    """
    if string is None:
        return {}

    if isinstance(string, dict):
        return string

    parser = create_bnf(allow_tuple=allow_tuple, free_word=free_word)

    out = {}
    for r in parser.parseString(string, parseAll=True):
        out.update(r)

    return out

def dict_from_options(options):
    """
    Return a dictionary that can be used to construct/override a ProblemConf
    instance based on `options`.

    See ``--conf`` and ``--options`` options of the ``simple.py`` script.
    """
    override = dict_from_string(options.conf)
    if options.app_options:
        if not 'options' in override:
            override['options'] = {}

        override_options = dict_from_string(options.app_options)
        override['options'].update(override_options)

    return override

##
# 27.10.2005, c
class ProblemConf(Struct):
    """
    Problem configuration, corresponding to an input (problem description
    file). It validates the input using lists of required and other keywords
    that have to/can appear in the input. Default keyword lists can be obtained
    by sfepy.base.conf.get_standard_keywords().

    ProblemConf instance is used to construct a Problem instance via
    Problem.from_conf(conf).
    """

    @staticmethod
    def from_file(filename, required=None, other=None, verbose=True,
                  define_args=None, override=None, setup=True):
        """
        Loads the problem definition from a file.

        The filename can either contain plain definitions, or it can contain
        the define() function, in which case it will be called to return the
        input definitions.

        The job of the define() function is to return a dictionary of
        parameters. How the dictionary is constructed is not our business, but
        the usual way is to simply have a function define() along these lines
        in the input file::

            def define():
                options = {
                    'save_eig_vectors' : None,
                    'eigen_solver' : 'eigen1',
                }
                region_2 = {
                    'name' : 'Surface',
                    'select' : 'nodes of surface',
                }
                return locals()

        Optionally, the define() function can accept additional arguments
        that should be defined using the `define_args` tuple or dictionary.
        """
        funmod = import_file(filename, package_name=False)

        if "define" in funmod.__dict__:
            if define_args is None:
                define_dict = funmod.__dict__["define"]()

            else:
                if isinstance(define_args, str):
                    define_args = dict_from_string(define_args)

                if isinstance(define_args, dict):
                    define_dict = funmod.__dict__["define"](**define_args)

                else:
                    define_dict = funmod.__dict__["define"](*define_args)

        else:
            define_dict = funmod.__dict__

        obj = ProblemConf(define_dict, funmod=funmod, filename=filename,
                          required=required, other=other, verbose=verbose,
                          override=override, setup=setup)

        return obj

    @staticmethod
    def from_file_and_options(filename, options, required=None, other=None,
                              verbose=True, define_args=None, setup=True):
        """
        Utility function, a wrapper around ProblemConf.from_file() with
        possible override taken from `options`.
        """
        override = dict_from_options(options)
        obj = ProblemConf.from_file(filename, required=required, other=other,
                                    verbose=verbose, define_args=define_args,
                                    override=override, setup=setup)
        return obj

    @staticmethod
    def from_module(module, required=None, other=None, verbose=True,
                    override=None, setup=True):
        obj = ProblemConf(module.__dict__, module, module.__name__,
                          required, other, verbose, override, setup=setup)

        return obj

    @staticmethod
    def from_dict(dict_, funmod, required=None, other=None, verbose=True,
                  override=None, setup=True):
        obj = ProblemConf(dict_, funmod, None, required, other, verbose,
                          override, setup=setup)

        return obj

    def __init__(self, define_dict, funmod=None, filename=None,
                 required=None, other=None, verbose=True, override=None,
                 setup=True):
        if override:
            if isinstance(override, Struct):
                override = override.__dict__
            define_dict = update_dict_recursively(define_dict, override, True)

        self.__dict__.update(define_dict)
        self.verbose = verbose

        if setup:
            self.setup(funmod=funmod, filename=filename,
                       required=required, other=other)


    def setup(self, define_dict=None, funmod=None, filename=None,
              required=None, other=None):

        define_dict = get_default(define_dict, self.__dict__)

        self._filename = filename

        self.validate(required=required, other=other)

        self.transform_input_trivial()
        self._raw = {}
        for key, val in six.iteritems(define_dict):
            if isinstance(val, dict):
                self._raw[key] = copy(val)

        self.transform_input()
        self.funmod = funmod

    def _validate_helper(self, items, but_nots):
        keys = list(self.__dict__.keys())
        left_over = keys[:]
        if but_nots is not None:
            for item in but_nots:
                match = re.compile('^' + item + '$').match
                for key in keys:
                    if match(key):
                        left_over.remove(key)

        missing = []
        if items is not None:
            for item in items:
                found = False
                match = re.compile('^' + item + '$').match
                for key in keys:
                    if match(key):
                        found = True
                        left_over.remove(key)
                if not found:
                    missing.append(item)
        return left_over, missing

    def validate(self, required=None, other=None):
        required_left_over, required_missing \
                            = self._validate_helper(required, other)
        other_left_over, other_missing \
                         = self._validate_helper(other, required)

        assert_(required_left_over == other_left_over)

        if other_left_over and self.verbose:
            output('left over:', other_left_over)

        if required_missing:
            raise ValueError('required missing: %s' % required_missing)

        return other_missing

    def transform_input_trivial(self):
        """Trivial input transformations."""

        ##
        # Unordered inputs.
        tr_list = ['([a-zA-Z0-9]+)_[0-9]+']
        # Keywords not in 'required', but needed even empty (e.g. for
        # running tests).
        for key in transforms.keys():
            if key not in self.__dict__:
                self.__dict__[key] = {}

        keys = list(self.__dict__.keys())
        for item in tr_list:
            match = re.compile(item).match
            for key in keys:
                obj = match(key)
                if obj:
                    new = obj.group(1) + 's'
                    result = {key : self.__dict__[key]}
                    try:
                        self.__dict__[new].update(result)
                    except:
                        self.__dict__[new] = result

                    del self.__dict__[key]

    def transform_input(self):
        keys = list(self.__dict__.keys())
        for key, transform in six.iteritems(transforms):
            if not key in keys: continue
            self.__dict__[key] = transform(self.__dict__[key])

    def get_raw(self, key=None):
        if key is None:
            return self._raw
        else:
            return self._raw[key]

    def get_item_by_name(self, key, item_name):
        """
        Return item with name `item_name` in configuration group given
        by `key`.
        """
        val = getattr(self, key)
        for item in six.itervalues(val):
            if item.name == item_name:
                return item

    def get_function(self, name):
        """
        Get a function object given its name.

        It can be either in `ProblemConf.funmod`, or a `ProblemConf`
        attribute directly.

        Parameters
        ----------
        name : str or function or None
            The function name or directly the function.

        Returns
        -------
        fun : function or None
            The required function, or None if `name` was `None`.
        """
        if name is None:
            fun = None

        elif callable(name):
            import inspect
            if not (inspect.isfunction(name) or inspect.ismethod(name)):
                msg = '`name` has to have `str` or `function` type! (got %s)'
                raise ValueError(msg % type(name))
            fun = name

        else:
            try:
                fun = getattr(self.funmod, name)

            except AttributeError:
                try:
                    fun = getattr(self, name)

                except AttributeError:
                    raise ValueError('function %s cannot be found!' % name)

        return fun

    def edit(self, key, newval):
        self.__dict__[key] = transforms[key](newval)

    def update_conf(self, conf):
        """
        Update configuration by values in another problem configuration.

        Values that are dictionaries are updated in-place by ``dict.update()``.

        Parameters
        ----------
        conf : ProblemConf instance
            The other configuration.
        """
        for x in conf.__dict__:
            his = conf.__dict__[x]
            my = getattr(self, x, None)
            if isinstance(my, dict) and isinstance(his, dict):
                my.update(his)
            else:
                setattr(self, x, his)

    def add_missing(self, conf):
        """
        Add missing values from another problem configuration.

        Missing keys/values are added also to values that are dictionaries.

        Parameters
        ----------
        conf : ProblemConf instance
            The other configuration.
        """
        for x in conf.__dict__:
            his = conf.__dict__[x]
            my = getattr(self, x, None)
            if isinstance(my, dict) and isinstance(his, dict):
                for key in his:
                    if key not in my:
                        my[key]=his[key]
            elif my is None:
                setattr(self, x, his)
