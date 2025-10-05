"""
Various global options/parameters.

Notes
-----
Inspired by rcParams of matplotlib.
"""

def validate_bool(val):
    """
    Convert b to a boolean or raise a ValueError.
    """
    if type(val) is str:
        val = val.lower()

    if val in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
        return True

    elif val in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
        return False

    else:
        raise ValueError('Could not convert "%s" to boolean!' % val)

default_goptions = {
    'verbose' : [True, validate_bool],
    'check_term_finiteness' : [False, validate_bool],
}

class ValidatedDict(dict):
    """
    A dictionary object including validation.
    """
    validate = dict([(key, validator) for key, (default, validator) in \
                     default_goptions.items()])

    def __setitem__(self, key, val):
        try:
            cval = self.validate[key](val)
            dict.__setitem__(self, key, cval)
        except KeyError:
            raise KeyError('%s is not a valid option.'
                           ' See goptions.keys() for a list of valid'
                           ' options.' % (key,))

    def keys(self):
        """
        Return sorted list of keys.
        """
        ks = dict.keys(self)
        ks.sort()
        return ks

    def values(self):
        """
        Return values in order of sorted keys.
        """
        return [self[key] for key in self.keys()]

goptions = ValidatedDict([(key, val[0])
                          for key, val in default_goptions.items()])
