"""
Some utilities for work with units of physical quantities.
"""
try:
    import sympy as sm
except ImportError:
    sm = None

import numpy as nm

from sfepy.base.base import get_default, invert_dict, is_sequence, Struct

default_units_of_basic_quantities = {
    'length' : 'm',
    'time' : 's',
    'mass' : 'kg',
    'temperature' : 'C',
}

# Cannot use N for Newton as it is a sympy function...
derived_units = {
    'Newton' : '(kg * m) / s**2',
    'Pa' : 'kg / (m * s**2)', # N / m**2
    'J' : '(kg * m**2) / s**2', # N m
}

units_of_quantities = {
    'density' : 'kg / m**3',
    'force' : 'Newton',
    'stress' : 'Pa',
    'energy' : 'J',
    'thermal_expandability' : 'Pa / C',
}

prefixes = {
    'p'  : 1e-12,
    'n'  : 1e-9,
    'mu' : 1e-6,
    'm'  : 1e-3,
    'c'  : 1e-2,
    ''   : 1e0,
    'd'  : 1e1,
    'k'  : 1e3,
    'M'  : 1e6,
    'G'  : 1e9,
    'T'  : 1e12,
}

inv_prefixes = invert_dict(prefixes)

num_prefixes = [str(ii) for ii in range(-20, 21)]
prefixes.update(dict(((key, 10**int(key)) for key in num_prefixes)))

class Unit(Struct):
    """
    A unit of a physical quantity. The prefix and coefficient of the unit
    are determined from to its name.

    Examples
    --------
    Construct some units:

    >>> from sfepy.mechanics.units import Unit
    >>> unit = Unit('mm')
    >>> print unit
    Unit:mm
      coef:
        0.001
      name:
        mm
      prefix:
        m
      prefix_length:
        1
      unit:
        m
    >>> unit = Unit('kg')
    >>> print unit
    Unit:kg
      coef:
        1000.0
      name:
        kg
      prefix:
        k
      prefix_length:
        1
      unit:
        g

    Get prefixes for a coefficient:

    >>> Unit.get_prefix(100.0)
    ('d', 10.0)
    >>> Unit.get_prefix(100.0, omit=('d',))
    ('k', 0.10000000000000001)
    """

    @staticmethod
    def get_prefix(coef, bias=0.1, omit=None):
        """
        Get the prefix and numerical multiplier corresponding to a numerical
        coefficient, omitting prefixes in omit tuple.
        """
        if omit is None:
            omit = num_prefixes

        values = [val for key, val in prefixes.items() if key not in omit]
        coefs = nm.array(values, dtype=nm.float64)
        coefs.sort()
        ii = nm.searchsorted(coefs, bias*coef, side='left')

        if ii == len(coefs):
            ii = ii - 1

        cc = coefs[ii]
        prefix = inv_prefixes[cc]
        mul = coef / cc

        return prefix, mul

    def __init__(self, name):
        self.name = name

        aux = sorted(list(prefixes.keys()), reverse=True)
        for prefix in aux:
            lp = len(prefix)
            if (prefix == name[:lp]) and (lp < len(name)): break

        self.prefix = prefix
        self.prefix_length = len(prefix)
        self.unit = name[self.prefix_length:]

        self.coef = prefixes[self.prefix]

class Quantity(Struct):
    """
    A physical quantity in a given set of basic units.

    Examples
    --------

    Construct the stress quantity:

    >>> from sfepy.mechanics.units import Unit, Quantity
    >>> units = ['m', 's', 'kg', 'C']
    >>> unit_set = [Unit(key) for key in units]
    >>> q1 = Quantity('stress', unit_set)
    >>> q1()
    '1.0 Pa'

    Show its unit using various prefixes:

    >>> q1('m')
    '1000.0 mPa'
    >>> q1('')
    '1.0 Pa'
    >>> q1('k')
    '0.001 kPa'
    >>> q1('M')
    '1e-06 MPa'

    Construct the stress quantity in another unit set:

    >>> units = ['mm', 's', 'kg', 'C']
    >>> unit_set = [Unit(key) for key in units]
    >>> q2 = Quantity('stress', unit_set)
    >>> q2()
    '1.0 kPa'

    Show its unit using various prefixes:

    >>> q2('m')
    '1000000.0 mPa'
    >>> q2('')
    '1000.0 Pa'
    >>> q2('k')
    '1.0 kPa'
    >>> q2('M')
    '0.001 MPa'
    """

    def __init__(self, name, unit_set):
        """
        Create a quantity in the given unit set. The name must be listed in the
        units_of_quantities dictionary.
        """
        self.name = name
        self.unit_set = unit_set

        self.unit_name = units_of_quantities[self.name]

        unit_expr = sm.sympify(self.unit_name)
        unit_expr = unit_expr.subs(derived_units)

        self.symbolic_value = sm.sympify(str(unit_expr))
        atoms = self.symbolic_value.atoms(sm.Symbol)
        self.def_units = [Unit(atom.name) for atom in atoms]

        self.def_names, self.def_units = self._get_dicts(self.def_units)
        self.names, self.units = self._get_dicts(self.unit_set)

        self.def_coef = float(self.symbolic_value.subs(self.def_names))

        coef_dict = {}
        for key, val in self.def_units.items():
            coef_dict[val.name] = self.units[key].coef
        self.coef_dict = coef_dict

        self.raw_coef = float(self.symbolic_value.subs(self.coef_dict))
        self.coef = self.raw_coef / self.def_coef

    def _get_dicts(self, unit_set):
        """Get auxiliary dictionaries for a unit set."""
        name_dict = {}
        unit_dict = {}
        for unit in unit_set:
            name_dict[unit.name] = unit.coef
            unit_dict[unit.unit] = unit

        return name_dict, unit_dict

    def __call__(self, prefix=None, omit=None):
        """Get the quantity units."""
        if prefix is None:
            if omit is None:
                omit = ['c', 'd'] + num_prefixes
            prefix, mul = Unit.get_prefix(self.coef, omit=omit)

        else:
            coef = prefixes[prefix]
            mul = self.coef / coef

        return '%s %s%s' % (mul, prefix, self.unit_name)

def get_consistent_unit_set(length=None, time=None, mass=None,
                            temperature=None):
    """
    Given a set of basic units, return a consistent set of derived units for
    quantities listed in the units_of_quantities dictionary.
    """
    defaults = default_units_of_basic_quantities
    length = get_default(length, defaults['length'])
    time = get_default(time, defaults['time'])
    mass = get_default(mass, defaults['mass'])
    temperature = get_default(temperature, defaults['temperature'])

    unit_set = [Unit(length), Unit(time), Unit(mass), Unit(temperature)]

    derived_units = {}

    for quantity_name in units_of_quantities.keys():
        quantity = Quantity(quantity_name, unit_set)

        derived_units[quantity_name] = quantity()

    return derived_units

def apply_unit_multipliers(values, unit_kinds, unit_multipliers):
    """
    Apply time, length and mass unit multipliers to given values with units
    corresponding to unit kinds.

    Returns
    -------
    new_values : list
        The new values with applied unit multipliers
    """
    a, b, c = unit_multipliers

    coefs = {
        'one' : 1.0 / 1.0,
        'time' : 1.0 / a,
        'length' : 1.0 / b,
        'mass' : 1.0 / c,
        'frequency' : 1.0 / (1.0 / a),
        'wave_number' : 1.0 / (1.0 / b),
        'velocity' : 1.0 / (b / a),
        'acceleration' : 1.0 / (b / a**2),
        'density' : 1.0 / (c / b**3),
        'force' : 1.0 / ((c * b) / a**2),
        'stress' : 1.0 / (c / (b * a**2)),
        'dyn_viscosity' : 1.0 / (c / (b * a)),
        'compressibility' : 1.0 / ((b * a**2) / c),
        'energy' : 1.0 / (c * b**2 / a**2),
    }
    avalues = [nm.asanyarray(val) if is_sequence(val) else val
               for val in values]
    new_values = [coefs[unit_kind] * avalues[ii]
                  for ii, unit_kind in enumerate(unit_kinds)]
    return new_values

def apply_units_to_pars(pars, pars_kinds, unit_multipliers):
    """
    Apply units in `unit_multipliers` to `pars` according to their kinds.

    Parameters
    ----------
    pars : dict
        The input parameters given as `name : value` items.
    pars_kinds : dict
        The kinds of the parameters given as `name : kind` items, with kinds
        defined in :func:`apply_unit_multipliers()`.
    unit_multipliers : tuple
        The time, length and mass unit multipliers.

    Returns
    -------
    new_pars : dict
        The output parameters.
    """
    keys = pars.keys()
    values = list(pars.values())
    kinds = [pars_kinds[key] for key in keys]

    new_values = apply_unit_multipliers(values, kinds, unit_multipliers)
    new_pars = dict(zip(keys, new_values))

    return new_pars
