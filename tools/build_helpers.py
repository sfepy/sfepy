"""
Build helpers for setup.py.

Notes
-----
The original version of this file was adapted from NiPy project [1].

[1] http://nipy.sourceforge.net/
"""

# Standard library imports
import os
import shutil
import glob
import fnmatch
import logging
from setuptools import Command
from skbuild.command.clean import clean
from packaging import version

from install_data import install_data

import sfepy.version as INFO

CYTHON_MIN_VERSION = INFO.CYTHON_MIN_VERSION

log = logging.getLogger(__name__)

class NoOptionsDocs(Command):
    user_options = [('None', None, 'this command has no options'),]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

def get_sphinx_make_command():
    if os.name in ['posix']:
        return 'make'

    elif os.name in ['nt']:
        return 'make.bat'

    else:
        raise ValueError('unsupported system! (%s)' % os.name)

class SphinxHTMLDocs(NoOptionsDocs):
    description = """generate html docs by Sphinx"""

    def run(self):
        os.chdir('doc')

        try:
            cmd = get_sphinx_make_command()
            os.system(cmd + ' html')

        finally:
            os.chdir('..')

class SphinxPDFDocs(NoOptionsDocs):
    description = """generate pdf docs by Sphinx"""

    def run(self):
        cwd = os.getcwd()

        os.chdir('doc')

        try:
            cmd = get_sphinx_make_command()
            os.system(cmd + ' latex')
            os.chdir('_build/latex')
            os.system(cmd + ' all-pdf')
            os.chdir(cwd)
            try:
                os.remove('doc/sfepy_manual.pdf')

            except:
                pass
            os.rename('doc/_build/latex/SfePy.pdf', 'doc/sfepy_manual.pdf')

        finally:
            os.chdir(cwd)

class DoxygenDocs(NoOptionsDocs):
    description = """generate docs by Doxygen"""

    def run(self):
        try:
            shutil.rmtree('doc/html')

        except OSError:
            pass

        fd_in = open('doc/doxygen.config', 'r')
        fd_out = open('doc/doxygenrc', 'w')
        for line in fd_in:
            aux = line.split('=')
            if len(aux) and (aux[0] == 'PROJECT_NUMBER'):
                line = '='.join([aux[0], INFO.__version__])

            fd_out.write(line)

        fd_out.close()
        fd_in.close()

        os.system('doxygen doc/doxygenrc')

def recursive_glob(top_dir, pattern):
    """
    Utility function working like `glob.glob()`, but working recursively
    and returning generator.

    Parameters
    ----------
    topdir : str
        The top-level directory.
    pattern : str or list of str
        The pattern or list of patterns to match.
    """
    if isinstance(pattern, list):
        for pat in pattern:
            for fn in recursive_glob(top_dir, pat):
                yield fn

    else:
        for dirpath, dirnames, filenames in os.walk(top_dir):
            for fn in [fn for fn in filenames
                       if fnmatch.fnmatchcase(fn, pattern)]:
                yield os.path.join(dirpath, fn)

class Clean(clean):
    """
    Command class to clean, enhanced to clean also files generated during
    `python setup.py build_ext --inplace`.
    """

    def run(self):
        clean.run(self)

        print('extra clean:')
        suffixes = ['*.pyc', '*.o', '*.so', '*.pyd',
                    '*_wrap.c', '*.bak', '*~', '*%']
        for filename in recursive_glob('sfepy', suffixes):
            print(filename)
            os.remove(filename)

        for filename in glob.glob('*.pyc'):
            print(filename)
            os.remove(filename)

        for _filename in recursive_glob('sfepy', ['*.pyx']):
            filename = _filename.replace('.pyx', '.c')
            print(filename)
            try:
                os.remove(filename)
            except OSError:
                pass

            filename = _filename.replace('.pyx', '.html')
            print(filename)
            try:
                os.remove(filename)
            except OSError:
                pass

# The command classes used by setup.py.
cmdclass = {
    'htmldocs' : SphinxHTMLDocs,
    'pdfdocs' : SphinxPDFDocs,
    'doxygendocs' : DoxygenDocs,
    'clean': Clean,
    'install_data': install_data
}

def have_good_cython():
    try:
        from Cython.Compiler.Version import version as cversion
    except ImportError:
        return False
    return version.parse(cversion) >= version.parse(CYTHON_MIN_VERSION)

def package_check(pkg_name, version=None, optional=False, checker=version.parse,
                  version_getter=None, messages=None, show_only=False):
    """
    Check if package `pkg_name` is present, and in correct version.

    Parameters
    ----------
    pkg_name : str or sequence of str
       The name of the package as imported into python. Alternative names
       (e.g. for different versions) may be given in a list.
    version : str, optional
       The minimum version of the package that is required. If not given, the
       version is not checked.
    optional : bool, optional
       If False, raise error for absent package or wrong version;
       otherwise warn
    checker : callable, optional
       If given, the callable with which to return a comparable thing from a
       version string. The default is ``packaging.version.parse()``.
    version_getter : callable, optional:
       If given, the callable that takes `pkg_name` as argument, and returns
       the package version string - as in::

          ``version = version_getter(pkg_name)``

       The default is equivalent to::

          mod = __import__(pkg_name); version = mod.__version__``
    messages : dict, optional
       If given, the dictionary providing (some of) output messages.
    show_only : bool
       If True, do not raise exceptions, only show the package name and version
       information.
    """
    if version_getter is None:
        def version_getter(pkg_name):
            mod = __import__(pkg_name)
            return mod.__version__
    if messages is None:
        messages = {}

    msgs = {
         'available' : '%s is available',
         'missing' : '%s is missing',
         'opt suffix' : '; you may get run-time errors',
         'version' : '%s is available in version %s',
         'version old' : '%s is available in version %s, but >= %s is needed',
         'no version' : '%s is available, cannot determine version',
    }
    msgs.update(messages)

    if isinstance(pkg_name, str):
        names = [pkg_name]

    else:
        names = pkg_name

    import_ok = False
    for pkg_name in names:
        try:
            __import__(pkg_name)
        except ImportError:
            pass
        else:
            import_ok = True

    pkg_info = pkg_name + (' (optional)' if optional else '')

    if not import_ok:
        if not (optional or show_only):
            raise RuntimeError(msgs['missing'] % pkg_name)
        log.warn(msgs['missing'] % pkg_info + msgs['opt suffix'])
        return

    if not version:
        if show_only:
            log.info(msgs['available'] % pkg_info)

        return

    try:
        have_version = version_getter(pkg_name)

    except AttributeError:
        raise RuntimeError(msgs['no version'] % pkg_info)

    if not have_version:
        if optional or show_only:
            log.warn(msgs['no version'] % pkg_info)

        else:
            raise RuntimeError(msgs['no version'] % pkg_info)

    elif checker(have_version) < checker(version):
        if optional or show_only:
            log.warn(msgs['version old'] % (pkg_info, have_version, version)
                     + msgs['opt suffix'])

        else:
            raise RuntimeError(msgs['version old'] % (pkg_info, have_version,
                                                      version))

    elif show_only:
        log.info(msgs['version'] % (pkg_info, have_version))
