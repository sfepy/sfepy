"""
Build helpers for setup.py.

Includes package dependency checks and monkey-patch to numpy.distutils
to work with Cython.

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
from distutils.cmd import Command
from os.path import join as pjoin, dirname
from distutils.command.clean import clean
from distutils.version import LooseVersion
from distutils.dep_util import newer_group
from distutils.errors import DistutilsError

from numpy.distutils.misc_util import appendpath
from numpy.distutils import log

import sfepy.version as INFO

CYTHON_MIN_VERSION = INFO.CYTHON_MIN_VERSION

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
    Distutils Command class to clean, enhanced to clean also files
    generated during `python setup.py build_ext --inplace`.
    """

    def run(self):
        clean.run(self)

        print 'extra clean:'
        suffixes = ['*.pyc', '*.o', '*.so', '*_wrap.c', '*.bak', '*~', '*%']
        for filename in recursive_glob('sfepy', suffixes):
            print filename
            os.remove(filename)

        for filename in glob.glob('*.pyc'):
            print filename
            os.remove(filename)

# The command classes for distutils, used by setup.py.
cmdclass = {
    'htmldocs' : SphinxHTMLDocs,
    'pdfdocs' : SphinxPDFDocs,
    'doxygendocs' : DoxygenDocs,
    'clean': Clean,
}

def have_good_cython():
    try:
        from Cython.Compiler.Version import version
    except ImportError:
        return False
    return LooseVersion(version) >= LooseVersion(CYTHON_MIN_VERSION)

def generate_a_pyrex_source(self, base, ext_name, source, extension):
    '''
    Monkey patch for numpy build_src.build_src method

    Uses Cython instead of Pyrex.
    '''
    good_cython = have_good_cython()
    if self.inplace or not good_cython:
        target_dir = dirname(base)
    else:
        target_dir = appendpath(self.build_src, dirname(base))
    target_file = pjoin(target_dir, ext_name + '.c')
    depends = [source] + extension.depends
    sources_changed = newer_group(depends, target_file, 'newer')
    if self.force or sources_changed:
        if good_cython:
            # add distribution (package-wide) include directories, in order
            # to pick up needed .pxd files for cython compilation
            incl_dirs = extension.include_dirs[:]
            dist_incl_dirs = self.distribution.include_dirs
            if not dist_incl_dirs is None:
                incl_dirs += dist_incl_dirs
            import Cython.Compiler.Main
            log.info("cythonc:> %s" % (target_file))
            self.mkpath(target_dir)
            options = Cython.Compiler.Main.CompilationOptions(
                defaults=Cython.Compiler.Main.default_options,
                include_path=incl_dirs,
                output_file=target_file)
            cython_result = Cython.Compiler.Main.compile(source,
                                                       options=options)
            if cython_result.num_errors != 0:
                raise DistutilsError("%d errors while compiling "
                                     "%r with Cython"
                                     % (cython_result.num_errors, source))
        elif sources_changed and os.path.isfile(target_file):
            raise DistutilsError("Cython >=%s required for compiling %r"
                                 " because sources (%s) have changed" %
                                 (CYTHON_MIN_VERSION, source,
                                  ','.join(depends)))
        else:
            raise DistutilsError("Cython >=%s required for compiling %r"
                                 " but not available" %
                                 (CYTHON_MIN_VERSION, source))
    return target_file

def package_check(pkg_name, version=None,
                  optional=False,
                  checker=LooseVersion,
                  version_getter=None,
                  messages=None
                  ):
    """
    Check if package `pkg_name` is present, and correct version

    Parameters
    ----------
    pkg_name : str
       name of package as imported into python
    version : {None, str}, optional
       minimum version of the package that we require. If None, we don't
       check the version.  Default is None
    optional : {False, True}, optional
       If False, raise error for absent package or wrong version;
       otherwise warn
    checker : callable, optional
       callable with which to return comparable thing from version
       string.  Default is ``distutils.version.LooseVersion``
    version_getter : {None, callable}:
       Callable that takes `pkg_name` as argument, and returns the
       package version string - as in::

          ``version = version_getter(pkg_name)``

       If None, equivalent to::

          mod = __import__(pkg_name); version = mod.__version__``
    messages : None or dict, optional
       dictionary giving output messages
    """
    if version_getter is None:
        def version_getter(pkg_name):
            mod = __import__(pkg_name)
            return mod.__version__
    if messages is None:
        messages = {}
    msgs = {
         'missing': 'Cannot import package "%s" - is it installed?',
         'missing opt': 'Missing optional package "%s"',
         'opt suffix' : '; you may get run-time errors',
         'version too old': 'You have version %s of package "%s"'
                            ' but we need version >= %s', }
    msgs.update(messages)
    try:
        __import__(pkg_name)
    except ImportError:
        if not optional:
            raise RuntimeError(msgs['missing'] % pkg_name)
        log.warn(msgs['missing opt'] % pkg_name +
                 msgs['opt suffix'])
        return
    if not version:
        return
    try:
        have_version = version_getter(pkg_name)
    except AttributeError:
        raise RuntimeError('Cannot find version for %s' % pkg_name)
    if checker(have_version) < checker(version):
        if optional:
            log.warn(msgs['version too old'] % (have_version,
                                                pkg_name,
                                                version)
                     + msgs['opt suffix'])
        else:
            raise RuntimeError(msgs['version too old'] % (have_version,
                                                          pkg_name,
                                                          version))
