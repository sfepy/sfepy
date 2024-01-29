from __future__ import absolute_import
from __future__ import print_function
import sys
import os
from copy import copy, deepcopy
from types import MethodType
from .getch import getch

import numpy as nm
import scipy.sparse as sp
import six

real_types = [nm.float64]
complex_types = [nm.complex128]

nm.set_printoptions(threshold=100)

from sfepy.base.goptions import goptions

sfepy_config_dir = os.path.expanduser('~/.sfepy')
if not os.path.exists(sfepy_config_dir):
    os.makedirs(sfepy_config_dir)

if sys.version_info[0] < 3:
    PY3 = False
    basestr = basestring
else:
    PY3 = True
    basestr = str

def get_debug():
    """
    Utility function providing ``debug()`` function.
    """
    try:
        import IPython

    except ImportError:
        debug = None

    else:
        old_excepthook = sys.excepthook

        def debug(frame=None, frames_back=1):
            if IPython.__version__ >= '0.11':
                from IPython.core.debugger import Pdb

                try:
                    ip = get_ipython()

                except NameError:
                    from IPython.frontend.terminal.embed \
                         import InteractiveShellEmbed
                    ip = InteractiveShellEmbed()

                colors = ip.colors

            else:
                from IPython.Debugger import Pdb
                from IPython.Shell import IPShell
                from IPython import ipapi

                ip = ipapi.get()
                if ip is None:
                    IPShell(argv=[''])
                    ip = ipapi.get()

                colors = ip.options.colors

            sys.excepthook = old_excepthook

            if frame is None:
                frame = sys._getframe(frames_back)

            Pdb(colors).set_trace(frame)

    if debug is None:
        import pdb
        debug = lambda frame=None, frames_back=1: pdb.set_trace()

    debug.__doc__ = """
    Start debugger on line where it is called, roughly equivalent to::

        import pdb; pdb.set_trace()

    First, this function tries to start an `IPython`-enabled
    debugger using the `IPython` API.

    When this fails, the plain old `pdb` is used instead.

    With IPython, one can say in what frame the debugger can stop.
    """

    return debug

debug = get_debug()

def debug_on_error():
    """
    Start debugger at the line where an exception was raised.
    """
    try:
        from IPython.core import ultratb

        except_hook = ultratb.FormattedTB(mode='Verbose',
                                          color_scheme='Linux', call_pdb=1)

    except ImportError:
        def except_hook(etype, value, tb):
            if hasattr(sys, 'ps1') or not sys.stderr.isatty():
                # We are in interactive mode or we don't have a tty-like
                # device, so we call the default hook.
                sys.__excepthook__(etype, value, tb)

            else:
                import traceback, pdb
                # We are NOT in interactive mode, print the exception...
                traceback.print_exception(etype, value, tb)
                print()
                # ...then start the debugger in post-mortem mode.
                pdb.post_mortem(tb)

    sys.excepthook = except_hook

def import_file(filename, package_name=None, can_reload=True):
    """
    Import a file as a module. The module is explicitly reloaded to
    prevent undesirable interactions.
    """
    from sfepy import base_dir

    top_dir = os.path.normpath(os.path.join(base_dir, '..'))
    path = os.path.dirname(os.path.normpath(os.path.realpath(filename)))

    if (package_name is None) and (top_dir == path[:len(top_dir)]):
        package_name = path[len(top_dir) + 1:].replace(os.sep, '.')
        path = top_dir

    if not path in sys.path:
        sys.path.append(path)
        remove_path = True

    else:
        remove_path = False

    name = os.path.splitext(os.path.basename(filename))[0]

    if package_name:
        mod = __import__('.'.join((package_name, name)), fromlist=[name])

    else:
        mod = __import__(name)

    if (name in sys.modules) and can_reload:
        if PY3:
            import importlib
            importlib.reload(mod)

        else:
            reload(mod)

    if remove_path:
        sys.path.remove(path)

    return mod

def try_imports(imports, fail_msg=None):
    """
    Try import statements until one succeeds.

    Parameters
    ----------
    imports : list
        The list of import statements.
    fail_msg : str
        If not None and no statement succeeds, a `ValueError` is raised with
        the given message, appended to all failed messages.

    Returns
    -------
    locals : dict
        The dictionary of imported modules.
    """
    msgs = []
    for imp in imports:
        try:
            exec(imp)
            break

        except Exception as inst:
            msgs.append(str(inst))

    else:
        if fail_msg is not None:
            msgs.append(fail_msg)
            raise ValueError('\n'.join(msgs))

    return locals()

def python_shell(frame=0):
    import code
    frame = sys._getframe(frame+1)
    code.interact(local=frame.f_locals)

def ipython_shell(frame=0):
    from IPython.terminal.embed import InteractiveShellEmbed
    ipshell = InteractiveShellEmbed()

    ipshell(stack_depth=frame+1)

def shell(frame=0):
    """
    Embed an IPython (if available) or regular Python shell in the given frame.
    """
    try:
        ipython_shell(frame=frame+2)

    except ImportError:
        python_shell(frame=frame+1)

def assert_(condition, msg='assertion failed!'):
    if not condition:
        raise ValueError(msg)

##
# c: 06.04.2005, r: 05.05.2008
def pause(msg=None):
    """
    Prints the line number and waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging.
    """
    f = sys._getframe(1)
    ff = f.f_code
    if (msg):
        print('%s, %d: %s(), %d: %s' % (ff.co_filename, ff.co_firstlineno,
                                        ff.co_name, f.f_lineno, msg))
    else:
        print('%s, %d: %s(), %d' % (ff.co_filename, ff.co_firstlineno,
                                    ff.co_name, f.f_lineno))
    spause()

##
# Silent pause.
# 18.02.2005, c
# 12.02.2007
def spause(msg=None):
    """
    Waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging. This function is called from pause().
    """
    if (msg):
        print(msg)
    sys.stdout.flush()
    ch = getch()
    if ch == 'q':
        sys.exit()

##
# 02.01.2005
class Struct(object):
    # 03.10.2005, c
    # 26.10.2005
    def __init__(self, **kwargs):
        if kwargs:
            self.__dict__.update(kwargs)

    def _format_sequence(self, seq, threshold):
        threshold_half = threshold // 2

        if len(seq) > threshold:
            out = ', '.join(str(ii) for ii in seq[:threshold_half]) \
                  + ', ..., ' \
                  + ', '.join(str(ii) for ii in seq[-threshold_half:])

        else:
            out = str(seq)

        return out

    # 08.03.2005
    def __str__(self):
        """Print instance class, name and items in alphabetical order.

        If the class instance has '_str_attrs' attribute, only the attributes
        listed there are taken into account. Other attributes are provided only
        as a list of attribute names (no values).

        For attributes that are Struct instances, if
        the listed attribute name ends with '.', the attribute is printed fully
        by calling str(). Otherwise only its class name/name are printed.

        Attributes that are NumPy arrays or SciPy sparse matrices are
        printed in a brief form.

        Only keys of dict attributes are printed. For the dict keys as
        well as list or tuple attributes only several edge items are
        printed if their length is greater than the threshold value 20.
        """
        return self._str()

    def _str(self, keys=None, threshold=20):
        ss = '%s' % self.__class__.__name__
        if hasattr(self, 'name'):
            ss += ':%s' % self.name
        ss += '\n'

        if keys is None:
            keys = list(self.__dict__.keys())

        str_attrs = sorted(Struct.get(self, '_str_attrs', keys))
        printed_keys = []
        for key in str_attrs:
            if key[-1] == '.':
                key = key[:-1]
                full_print = True
            else:
                full_print = False

            printed_keys.append(key)

            try:
                val = getattr(self, key)

            except AttributeError:
                continue

            if isinstance(val, Struct):
                if not full_print:
                    ss += '  %s:\n    %s' % (key, val.__class__.__name__)
                    if hasattr(val, 'name'):
                        ss += ':%s' % val.name
                    ss += '\n'

                else:
                    aux = '\n' + str(val)
                    aux = aux.replace('\n', '\n    ')
                    ss += '  %s:\n%s\n' % (key, aux[1:])

            elif isinstance(val, dict):
                sval = self._format_sequence(list(val.keys()), threshold)
                sval = sval.replace('\n', '\n    ')
                ss += '  %s:\n    dict with keys: %s\n' % (key, sval)

            elif isinstance(val, list):
                sval = self._format_sequence(val, threshold)
                sval = sval.replace('\n', '\n    ')
                ss += '  %s:\n    list: %s\n' % (key, sval)

            elif isinstance(val, tuple):
                sval = self._format_sequence(val, threshold)
                sval = sval.replace('\n', '\n    ')
                ss += '  %s:\n    tuple: %s\n' % (key, sval)

            elif isinstance(val, nm.ndarray):
                ss += '  %s:\n    %s array of %s\n' \
                      % (key, val.shape, val.dtype)

            elif isinstance(val, sp.spmatrix):
                ss += '  %s:\n    %s spmatrix of %s, %d nonzeros\n' \
                      % (key, val.shape, val.dtype, val.nnz)

            else:
                aux = '\n' + str(val)
                aux = aux.replace('\n', '\n    ')
                ss += '  %s:\n%s\n' % (key, aux[1:])

        other_keys = sorted(set(keys).difference(set(printed_keys)))
        if len(other_keys):
            ss += '  other attributes:\n    %s\n' \
                  % '\n    '.join(key for key in other_keys)

        return ss.rstrip()

    def __repr__(self):
        ss = "%s" % self.__class__.__name__
        if hasattr(self, 'name'):
            ss += ":%s" % self.name
        return ss

    ##
    # 28.08.2007, c
    def __add__(self, other):
        """Merge Structs. Attributes of new are those of self unless an
        attribute and its counterpart in other are both Structs - these are
        merged then."""
        new = copy(self)
        for key, val in six.iteritems(other.__dict__):
            if hasattr(new, key):
                sval = getattr(self, key)
                if issubclass(sval.__class__, Struct) and \
                        issubclass(val.__class__, Struct):
                    setattr(new, key, sval + val)
                else:
                    setattr(new, key, sval)
            else:
                setattr(new, key, val)
        return new

    ##
    # 28.08.2007, c
    def __iadd__(self, other):
        """Merge Structs in place. Attributes of self are left unchanged
        unless an attribute and its counterpart in other are both Structs -
        these are merged then."""
        for key, val in six.iteritems(other.__dict__):
            if hasattr(self, key):
                sval = getattr(self, key)
                if issubclass(sval.__class__, Struct) and \
                       issubclass(val.__class__, Struct):
                    setattr(self, key, sval + val)
            else:
                setattr(self, key, val)
        return self

    def str_class(self):
        """
        As __str__(), but for class attributes.
        """
        return self._str(list(self.__class__.__dict__.keys()))

    # 08.03.2005, c
    def str_all(self):
        ss = "%s\n" % self.__class__
        for key, val in six.iteritems(self.__dict__):
            if issubclass(self.__dict__[key].__class__, Struct):
                ss += "  %s:\n" % key
                aux = "\n" + self.__dict__[key].str_all()
                aux = aux.replace("\n", "\n    ")
                ss += aux[1:] + "\n"
            else:
                aux = "\n" + str(val)
                aux = aux.replace("\n", "\n    ")
                ss += "  %s:\n%s\n" % (key, aux[1:])
        return(ss.rstrip())

    ##
    # 09.07.2007, c
    def to_dict(self):
        return copy(self.__dict__)

    def get(self, key, default=None, msg_if_none=None):
        """
        A dict-like get() for Struct attributes.
        """
        out = getattr(self, key, default)

        if (out is None) and (msg_if_none is not None):
            raise ValueError(msg_if_none)

        return out

    def update(self, other, **kwargs):
        """
        A dict-like update for Struct attributes.
        """
        if other is None: return

        if not isinstance(other, dict):
            other = other.to_dict()
        self.__dict__.update(other, **kwargs)

    def set_default(self, key, default=None):
        """
        Behaves like dict.setdefault().
        """
        return self.__dict__.setdefault(key, default)

    def copy(self, deep=False, name=None):
        """Make a (deep) copy of self.

        Parameters:

        deep : bool
            Make a deep copy.
        name : str
            Name of the copy, with default self.name + '_copy'.
        """
        if deep:
            other = deepcopy(self)
        else:
            other = copy(self)

        if hasattr(self, 'name'):
            other.name = get_default(name, self.name + '_copy')

        return other
#
# 12.07.2007, c
class IndexedStruct(Struct):

    ##
    # 12.07.2007, c
    def __getitem__(self, key):
        return getattr(self, key)

    ##
    # 12.07.2007, c
    def __setitem__(self, key, val):
        setattr(self, key, val)

    def setdefault(self, key, default=None):
        return self.set_default(key, default=default)

##
# 14.07.2006, c
class Container(Struct):

    def __init__(self, objs=None, **kwargs):
        Struct.__init__(self, **kwargs)

        if objs is not None:
            self._objs = objs
            self.update()
        else:
            self._objs = []
            self.names = []

    def update(self, objs=None):
        if objs is not None:
            self._objs = objs

        self.names = [obj.name for obj in self._objs]

    def __setitem__(self, ii, obj):
        try:
            if isinstance(ii, basestr):
                if ii in self.names:
                    ii = self.names.index(ii)
                else:
                    ii = len(self.names)

            elif not isinstance(ii, int):
                    raise ValueError('bad index type! (%s)' % type(ii))

            if ii >= len(self.names):
                self._objs.append(obj)
                self.names.append(obj.name)

            else:
                self._objs[ii] = obj
                self.names[ii] = obj.name

        except (IndexError, ValueError) as msg:
            raise IndexError(msg)

    def __getitem__(self, ii):
        try:
            if isinstance(ii, basestr):
                ii = self.names.index(ii)
            elif not isinstance(ii, int):
                raise ValueError('bad index type! (%s)' % type(ii))

            return  self._objs[ii]

        except (IndexError, ValueError) as msg:
            raise IndexError(msg)

    def __iter__(self):
        return self._objs.__iter__()

    def __add__(self, other):
        """
        Add items of `other` to `self`.
        """
        new = Container()
        objs = self._objs + other._objs
        new.update(objs)

        return new

    def __iadd__(self, other):
        """
        Add items of `other` to `self` in place.
        """
        self.extend(copy(other._objs))
        self.update()

        return self

    ##
    # 18.07.2006, c
    def __len__(self):
        return len(self._objs)

    def insert(self, ii, obj):
        self._objs.insert(ii, obj)
        self.names.insert(ii, obj.name)

    def append(self, obj):
        self[len(self.names)] = obj

    def extend(self, objs):
        """
        Extend the container items by the sequence `objs`.
        """
        for obj in objs:
            self.append(obj)

    def get(self, ii, default=None, msg_if_none=None):
        """
        Get an item from Container - a wrapper around
        Container.__getitem__() with defaults and custom error message.

        Parameters
        ----------
        ii : int or str
            The index or name of the item.
        default : any, optional
            The default value returned in case the item `ii` does not exist.
        msg_if_none : str, optional
            If not None, and if `default` is None and the item `ii` does
            not exist, raise ValueError with this message.
        """
        try:
            out = self[ii]

        except (IndexError, ValueError):
            if default is not None:
                out = default

            else:
                if msg_if_none is not None:
                    raise ValueError(msg_if_none)

                else:
                    raise

        return out

    def remove_name(self, name):
        ii = self.names.index[name]
        del self.names[ii]
        del self._objs[ii]

    ##
    # dict-like methods.
    def itervalues(self):
        return self._objs.__iter__()

    def iterkeys(self):
        return self.get_names().__iter__()

    def iteritems(self):
        for obj in self._objs:
            yield obj.name, obj

    ##
    # 20.09.2006, c
    def has_key(self, ii):
        if isinstance(ii, int):
            if (ii < len(self)) and (ii >= (-len(self))):
                return True
            else:
                return False
        elif isinstance(ii, basestr):
            try:
                self.names.index(ii)
                return True
            except:
                return False
        else:
            raise IndexError('unsupported index type: %s' % ii)

    ##
    # 12.06.2007, c
    def print_names(self):
        print([obj.name for obj in self._objs])

    def get_names(self):
        return [obj.name for obj in self._objs]

    def as_dict(self):
        """
        Return stored objects in a dictionary with object names as keys.
        """
        out = {}
        for key, val in self.iteritems():
            out[key] = val

        return out

##
# 30.11.2004, c
# 01.12.2004
# 01.12.2004
class OneTypeList(list):

    def __init__(self, item_class, seq=None):
        self.item_class = item_class

        if seq is not None:
            for obj in seq:
                self.append(obj)

    def __setitem__(self, key, value):
        if (type(value) in (list, tuple)):
            for ii, val in enumerate(value):
                if not isinstance(val, self.item_class):
                    raise TypeError
        else:
            if not isinstance(value, self.item_class):
                raise TypeError
        list.__setitem__(self, key, value)

    ##
    # 21.11.2005, c
    def __getitem__(self, ii):
        if isinstance(ii, int):
            return list.__getitem__(self, ii)
        elif isinstance(ii, basestr):
            ir = self.find(ii, ret_indx=True)
            if ir:
                return list.__getitem__(self, ir[0])
            else:
                raise IndexError(ii)
        else:
            raise IndexError(ii)

    def __str__(self):
        ss = "[\n"
        for ii in self:
            aux = "\n" + ii.__str__()
            aux = aux.replace("\n", "\n  ")
            ss += aux[1:] + "\n"
        ss += "]"
        return(ss)

    def find(self, name, ret_indx=False):
        for ii, item in enumerate(self):
            if item.name == name:
                if ret_indx:
                    return ii, item
                else:
                    return item
        return None

    ##
    # 12.06.2007, c
    def print_names(self):
        print([ii.name for ii in self])

    def get_names(self):
        return [ii.name for ii in self]

class Output(Struct):
    """
    Factory class providing output (print) functions. All SfePy
    printing should be accomplished by this class.

    Examples
    --------
    >>> from sfepy.base.base import Output
    >>> output = Output('sfepy:')
    >>> output(1, 2, 3, 'hello')
    sfepy: 1 2 3 hello
    >>> output.prefix = 'my_cool_app:'
    >>> output(1, 2, 3, 'hello')
    my_cool_app: 1 2 3 hello
    """

    def __init__(self, prefix, filename=None, quiet=False, combined=False,
                 append=False, **kwargs):
        Struct.__init__(self, **kwargs)

        self.prefix = prefix

        self.set_output(filename=filename, quiet=quiet,
                        combined=combined, append=append)

    def __call__(self, *argc, **argv):
        """Call self.output_function.

        Parameters
        ----------
        argc : positional arguments
            The values to print.
        argv : keyword arguments
            The arguments to control the output behaviour. Supported keywords
            are listed below.
        verbose : bool (in **argv)
            No output if False.
        """
        verbose = argv.get('verbose', goptions['verbose'])
        if verbose:
            self.output_function(*argc, **argv)

    def set_output(self, filename=None, quiet=False, combined=False,
                   append=False):
        """
        Set the output mode.

        If `quiet` is `True`, no messages are printed to screen. If
        simultaneously `filename` is not `None`, the messages are logged
        into the specified file.

        If `quiet` is `False`, more combinations are possible. If
        `filename` is `None`, output is to screen only, otherwise it is
        to the specified file. Moreover, if `combined` is `True`, both
        the ways are used.

        Parameters
        ----------
        filename : str or file object
            Print messages into the specified file.
        quiet : bool
            Do not print anything to screen.
        combined : bool
            Print both on screen and into the specified file.
        append : bool
            Append to an existing file instead of overwriting it. Use with
            `filename`.
        """
        if not isinstance(filename, basestr):
            # filename is a file descriptor.
            append = True

        self.level = 0

        def output_none(*argc, **argv):
            pass

        def output_screen(*argc, **argv):
            format = '%s' + ' %s' * (len(argc) - 1)
            msg = format % argc

            if msg.startswith('...'):
                self.level -= 1

            print(self._prefix + ('  ' * self.level) + msg)

            if msg.endswith('...'):
                self.level += 1

        def print_to_file(filename, msg):
            if isinstance(filename, basestr):
                fd = open(filename, 'a')

            else:
                fd = filename

            print(self._prefix + ('  ' * self.level) + msg, file=fd)

            if isinstance(filename, basestr):
                fd.close()

            else:
                fd.flush()

        def output_file(*argc, **argv):
            format = '%s' + ' %s' * (len(argc) - 1)
            msg = format % argc

            if msg.startswith('...'):
                self.level -= 1

            print_to_file(filename, msg)

            if msg.endswith('...'):
                self.level += 1

        def output_combined(*argc, **argv):
            format = '%s' + ' %s' * (len(argc) - 1)
            msg = format % argc

            if msg.startswith('...'):
                self.level -= 1

            print(self._prefix + ('  ' * self.level) + msg)

            print_to_file(filename, msg)

            if msg.endswith('...'):
                self.level += 1

        def reset_file(filename):
            if isinstance(filename, basestr):
                output_dir = os.path.dirname(filename)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                fd = open(filename, 'w')
                fd.close()

            else:
                raise ValueError('cannot reset a file object!')

        if quiet is True:
            if filename is not None:
                if not append:
                    reset_file(filename)

                self.output_function = output_file

            else:
                self.output_function = output_none

        else:
            if filename is None:
                self.output_function = output_screen

            else:
                if not append:
                    reset_file(filename)

                if combined:
                    self.output_function = output_combined

                else:
                    self.output_function = output_file

    def get_output_function(self):
        return self.output_function

    def set_output_prefix(self, prefix):
        assert_(isinstance(prefix, basestr))
        if len(prefix) > 0:
            prefix += ' '
        self._prefix = prefix

    def get_output_prefix(self):
        return self._prefix[:-1]
    prefix = property(get_output_prefix, set_output_prefix)

output = Output('sfepy:')

def configure_output(options):
    """
    Configure the standard :func:`output()` function using
    `output_log_name` and `output_screen` attributes of `options`.

    Parameters
    ----------
    options : Struct or dict
        The options with `output_screen` and `output_log_name` items. Defaults
        are provided if missing.
    """
    output_screen = options.get('output_screen', True)
    output_log_name = options.get('output_log_name', None)

    output.set_output(filename=output_log_name, quiet=not output_screen,
                      combined=output_screen and (output_log_name is not None))

def print_structs(objs):
    """Print Struct instances in a container, works recursively. Debugging
    utility function."""
    if isinstance(objs, dict):
        for key, vals in six.iteritems(objs):
            print(key)
            print_structs(vals)
    elif isinstance(objs, list):
        for vals in objs:
            print_structs(vals)
    else:
        print(objs)

def iter_dict_of_lists(dol, return_keys=False):
    for key, vals in six.iteritems(dol):
        for ii, val in enumerate(vals):
            if return_keys:
                yield key, ii, val
            else:
                yield val

##
# 19.07.2005, c
# 26.05.2006
# 17.10.2007
def dict_to_struct(*args, **kwargs):
    """Convert a dict instance to a Struct instance."""
    try:
        level = kwargs['level']
    except:
        level = 0

    try:
        flag = kwargs['flag']
    except:
        flag = (1,)

    # For level 0 only...
    try:
        constructor = kwargs['constructor']
    except:
        constructor = Struct

    out = []
    for arg in args:
        if type(arg) == dict:
            if flag[level]:
                aux = constructor()
            else:
                aux = {}

            for key, val in six.iteritems(arg):
                if type(val) == dict:
                    try:
                        flag[level + 1]
                    except:
                        flag = flag + (0,)
                    val2 = dict_to_struct(val, level=level + 1, flag=flag)
                    if flag[level]:
                        aux.__dict__[key] = val2
                    else:
                        aux[key] = val2
                else:
                    if flag[level]:
                        aux.__dict__[key] = val
                    else:
                        aux[key] = val
            out.append(aux)
        else:
            out.append(arg)

    if len(out) == 1:
        out = out[0]

    return out

def structify(obj):
    """
    Convert a (nested) dict `obj` into a (nested) Struct.
    """
    out = Struct(**obj)
    for key, val in out.to_dict().items():
        if isinstance(val, dict):
            out.__dict__[key] = structify(val)
    return out

def is_string(var):
    return isinstance(var, basestr)

def is_integer(var):
    if PY3:
        return isinstance(var, int)

    else:
        return isinstance(var, (int, long))

##
# 23.01.2006, c
def is_sequence(var):
    try:
        from collections.abc import Sequence
    except ImportError:
        from collections import Sequence
    if isinstance(var, basestr):
        return False
    return isinstance(var, Sequence)

##
# 17.10.2007, c
def is_derived_class(cls, parent):
    return issubclass(cls, parent) and (cls is not parent)

##
# 23.10.2007, c
def insert_static_method(cls, function):
    setattr(cls, function.__name__, staticmethod(function))

##
# 23.10.2007, c
def insert_method(instance, function):
    if PY3:
        meth = MethodType(function, instance)
    else:
        meth = MethodType(function, instance, type(instance))
    setattr(instance, function.__name__, meth)

def use_method_with_name(instance, method, new_name):
    setattr(instance, new_name, method)

def insert_as_static_method(cls, name, function):
    setattr(cls, name, staticmethod(function))

def find_subclasses(context, classes, omit_unnamed=False, name_attr='name'):
    """Find subclasses of the given classes in the given context.

    Examples
    --------

    >>> solver_table = find_subclasses(vars().items(),
                                       [LinearSolver, NonlinearSolver,
                                        TimeSteppingSolver, EigenvalueSolver,
                                        OptimizationSolver])
    """
    var_dict = list(context.items())
    table = {}

    for key, var in var_dict:
        try:
            for cls in classes:
                if is_derived_class(var, cls):
                    if hasattr(var, name_attr):
                        key = getattr(var, name_attr)
                        if omit_unnamed and not key:
                            continue

                    elif omit_unnamed:
                        continue

                    else:
                        key = var.__class__.__name__

                    table[key] = var
                    break

        except TypeError:
            pass
    return table

def load_classes(filenames, classes, package_name=None, ignore_errors=False,
                 name_attr='name'):
    """
    For each filename in filenames, load all subclasses of classes listed.
    """
    table = {}
    for filename in filenames:
        if not ignore_errors:
            mod = import_file(filename, package_name=package_name,
                              can_reload=False)

        else:
            try:
                mod = import_file(filename, package_name=package_name,
                                  can_reload=False)

            except:
                output('WARNING: module %s cannot be imported!' % filename)
                output('reason:\n', sys.exc_info()[1])
                continue

        table.update(find_subclasses(vars(mod), classes, omit_unnamed=True,
                                     name_attr=name_attr))

    return table

def update_dict_recursively(dst, src, tuples_too=False,
                            overwrite_by_none=True):
    """
    Update `dst` dictionary recursively using items in `src` dictionary.

    Parameters
    ----------
    dst : dict
        The destination dictionary.
    src : dict
        The source dictionary.
    tuples_too : bool
        If True, recurse also into dictionaries that are members of tuples.
    overwrite_by_none : bool
        If False, do not overwrite destination dictionary values by None.

    Returns
    -------
    dst : dict
        The destination dictionary.
    """
    def tuplezip(a):
        if isinstance(a[0], dict) and isinstance(a[1], dict):
            return update_dict_recursively(a[0], a[1], True)
        return a[1]

    for key in src:
        if key in dst:
            if isinstance(src[key], dict) and isinstance(dst[key], dict):
                dst[key] = update_dict_recursively(dst[key],
                                                   src[key], tuples_too)
                continue

            if tuples_too and isinstance(dst[key], tuple) \
                   and isinstance(src[key], tuple):
                out = map(tuplezip, zip(src[key], dst[key]))
                out = tuple(out)
                dst[key] = out[:len(dst[key])]
                continue

        if overwrite_by_none or not src[key] is None:
            dst[key] = src[key]

    return dst

def edit_tuple_strings(str_tuple, old, new, recur=False):
    """
    Replace substrings `old` with `new` in items of tuple
    `str_tuple`. Non-string items are just copied to the new tuple.

    Parameters
    ----------
    str_tuple : tuple
        The tuple with string values.
    old : str
        The old substring.
    new : str
        The new substring.
    recur : bool
        If True, edit items that are tuples recursively.

    Returns
    -------
    new_tuple : tuple
        The tuple with edited strings.
    """
    new_tuple = []
    for item in str_tuple:
        if isinstance(item, basestr):
            item = item.replace(old, new)

        elif recur and isinstance(item, tuple):
            item = edit_tuple_strings(item, old, new, recur=True)

        new_tuple.append(item)

    return tuple(new_tuple)

def edit_dict_strings(str_dict, old, new, recur=False):
    """
    Replace substrings `old` with `new` in string values of dictionary
    `str_dict`. Both `old` and `new` can be lists of the same length - items
    in `old` are replaced by items in `new` with the same index.

    Parameters
    ----------
    str_dict : dict
        The dictionary with string values or tuples containing strings.
    old : str or list of str
        The old substring or list of substrings.
    new : str or list of str
        The new substring or list of substrings.
    recur : bool
        If True, edit tuple values recursively.

    Returns
    -------
    new_dict : dict
        The dictionary with edited strings.
    """
    if isinstance(old, basestr):
        new_dict = {}
        for key, val in six.iteritems(str_dict):
            if isinstance(val, basestr):
                new_dict[key] = val.replace(old, new)

            elif isinstance(val, tuple):
                new_dict[key] = edit_tuple_strings(val, old, new, recur=recur)

            else:
                raise ValueError('unsupported value! (%s)' % type(val))

    else:
        assert_(len(old) == len(new))

        new_dict = dict(str_dict)
        for ii, _old in enumerate(old):
            new_dict.update(edit_dict_strings(new_dict, _old, new[ii],
                                              recur=recur))

    return new_dict

def invert_dict(d, is_val_tuple=False, unique=True):
    """
    Invert a dictionary by making its values keys and vice versa.

    Parameters
    ----------
    d : dict
        The input dictionary.
    is_val_tuple : bool
        If True, the `d` values are tuples and new keys are the tuple items.
    unique : bool
        If True, the `d` values are unique and so the mapping is
        one to one. If False, the `d` values (possibly) repeat, so the inverted
        dictionary will have as items lists of corresponding keys.

    Returns
    -------
    di : dict
        The inverted dictionary.
    """
    di = {}

    for key, val in six.iteritems(d):
        if unique:
            if is_val_tuple:
                for v in val:
                    di[v] = key
            else:
                di[val] = key

        else:
            if is_val_tuple:
                for v in val:
                    item = di.setdefault(v, [])
                    item.append(key)

            else:
                item = di.setdefault(val, [])
                item.append(key)

    return di

def remap_dict(d, map):
    """
    Utility function to remap state dict keys according to var_map.
    """
    out = {}
    for new_key, key in six.iteritems(map):
        out[new_key] = d[key]

    return out

##
# 24.08.2006, c
# 05.09.2006
def dict_from_keys_init(keys, seq_class=None):

    if seq_class is None:
        return {}.fromkeys(keys)

    out = {}
    for key in keys:
        out[key] = seq_class()
    return out

##
# 16.10.2006, c
def dict_extend(d1, d2):
    for key, val in six.iteritems(d1):
        val.extend(d2[key])

def get_subdict(adict, keys):
    """
    Get a sub-dictionary of `adict` with given `keys`.
    """
    return dict((key, adict[key]) for key in keys if key in adict)

def set_defaults(dict_, defaults):
    for key, val in six.iteritems(defaults):
        dict_.setdefault(key, val)

##
# c: 12.03.2007, r: 04.04.2008
def get_default(arg, default, msg_if_none=None):
    if arg is None:
        out = default
    else:
        out = arg

    if (out is None) and (msg_if_none is not None):
        raise ValueError(msg_if_none)

    return out

##
# c: 28.04.2008, r: 28.04.2008
def get_default_attr(obj, attr, default, msg_if_none=None):
    if hasattr(obj, attr):
        out = getattr(obj, attr)
    else:
        out = default

    if (out is None) and (msg_if_none is not None):
        raise ValueError(msg_if_none)

    return out

def get_arguments(omit=None):
    """Get a calling function's arguments.

    Returns:

    args : dict
        The calling function's  arguments.
    """
    from inspect import getargvalues, stack
    if omit is None:
        omit = []

    _args, _, _, _vars = getargvalues(stack()[1][0])

    args = {}
    for name in _args:
        if name in omit: continue
        args[name] = _vars[name]

    return args

def check_names(names1, names2, msg):
    """Check if all names in names1 are in names2, otherwise raise IndexError
    with the provided message msg.
    """
    names = set(names1)
    both = names.intersection(names2)
    if both != names:
        missing = ', '.join(ii for ii in names.difference(both))
        raise IndexError(msg % missing)

##
# c: 27.02.2008, r: 27.02.2008
def select_by_names(objs_all, names, replace=None, simple=True):
    objs = {}
    for key, val in six.iteritems(objs_all):
        if val.name in names:
            if replace is None:
                objs[key] = val
            else:
                new_val = copy(val)
                old_attr = getattr(val, replace[0])
                if simple:
                    new_attr = old_attr % replace[1]
                    setattr(new_val, replace[0], new_attr)
                else:
                    new_attr = replace[1].get(val.name, old_attr)
                    setattr(new_val, replace[0], new_attr)
                objs[key] = new_val
    return objs

def ordered_iteritems(adict):
    keys = list(adict.keys())
    order = nm.argsort(keys)
    for ii in order:
        key = keys[ii]
        yield key, adict[key]

def dict_to_array(adict):
    """
    Convert a dictionary of nD arrays of the same shapes with
    non-negative integer keys to a single (n+1)D array.
    """
    keys = list(adict.keys())
    ik = nm.array(keys, dtype=nm.int32)
    assert_((ik >= 0).all())

    if ik.shape[0] == 0:
        return nm.zeros((0,), dtype=nm.int32)

    aux = nm.asarray(adict[ik[0]])
    out = nm.empty((ik.max() + 1,) + aux.shape, dtype=aux.dtype)
    out.fill(-1)
    for key, val in six.iteritems(adict):
        out[key] = val

    return out

def as_float_or_complex(val):
    """
    Try to cast val to Python float, and if this fails, to Python
    complex type.
    """
    success = False
    try:
        out = float(val)
    except:
        pass
    else:
        success = True

    if not success:
        try:
            out = complex(val)
        except:
            pass
        else:
            success = True

    if not success:
        raise ValueError('cannot cast %s to float or complex!' % val)

    return out
