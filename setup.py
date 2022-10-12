#!/usr/bin/env python
"""SfePy: Simple finite elements in Python

SfePy (simple finite elements in Python) is a software, distributed
under the BSD license, for solving systems of coupled partial
differential equations by the finite element method. The code is based
on NumPy and SciPy packages.
"""

DOCLINES = __doc__.split("\n")

import os
import sys
sys.path.append('./tools')

from build_helpers import (generate_a_pyrex_source, package_check, log,
                           cmdclass, INFO)
# monkey-patch numpy distutils to use Cython instead of Pyrex
from numpy.distutils.command.build_src import build_src

build_src.generate_a_pyrex_source = generate_a_pyrex_source

VERSION = INFO.__version__

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
"""

DOWNLOAD_URL = "http://sfepy.org/doc-devel/downloads.html"

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('sfepy')

    config.add_data_files(('sfepy', ('VERSION', 'LICENSE')))
    config.add_data_dir(('sfepy/meshes', 'meshes'))
    config.add_data_dir(('sfepy/examples', 'sfepy/examples'))
    config.add_data_dir(('sfepy/tests', 'sfepy/tests'))

    config.get_version('sfepy/version.py')  # sets config.version

    return config


def _cython_version(pkg_name):
    from Cython.Compiler.Version import version

    return version

def _igakit_version(pkg_name):
    return '0.1'

def _pymetis_version(pkg_name):
    import pymetis

    return pymetis.version

def _scikit_umfpack_version(pkg_name):
    try:
        import scikits.umfpack; scikits.umfpack
        try:
            return scikits.umfpack.__version__

        except AttributeError:
            return '<0.3.1'

    except:
        return None

def check_versions(show_only=False):
    # Cython is a build dependency.
    package_check('cython', INFO.CYTHON_MIN_VERSION,
                  version_getter=_cython_version,
                  show_only=show_only)

    # Check hard and soft dependencies.
    package_check('numpy', INFO.NUMPY_MIN_VERSION,
                  show_only=show_only)
    package_check('scipy', INFO.SCIPY_MIN_VERSION,
                  show_only=show_only)
    package_check('matplotlib', INFO.MATPLOTLIB_MIN_VERSION,
                  show_only=show_only)
    package_check('pyparsing', INFO.PYPARSING_MIN_VERSION,
                  show_only=show_only)
    package_check('tables', INFO.PYTABLES_MIN_VERSION,
                  show_only=show_only)
    package_check('sympy', INFO.SYMPY_MIN_VERSION, optional=True,
                  messages={'opt suffix' : '; some tests are going to fail!'},
                  show_only=show_only)
    package_check('igakit', INFO.IGAKIT_MIN_VERSION, optional=True,
                  version_getter=_igakit_version,
                  show_only=show_only)
    package_check('petsc4py', INFO.PETSC4PY_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('mpi4py', INFO.MPI4PY_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('slepc4py', INFO.SLEPC4PY_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('pymetis', INFO.PYMETIS_MIN_VERSION, optional=True,
                  version_getter=_pymetis_version,
                  show_only=show_only)
    package_check('scikits.umfpack', INFO.SCIKIT_UMFPACK_MIN_VERSION,
                  optional=True,
                  version_getter=_scikit_umfpack_version,
                  show_only=show_only)
    package_check('meshio', INFO.MESHIO_MIN_VERSION,
                  show_only=show_only)
    package_check('psutil', INFO.PSUTIL_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('pyvista', INFO.PYVISTA_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('opt_einsum', INFO.OPT_EINSUM_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('jax', INFO.JAX_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('dask', INFO.DASK_MIN_VERSION, optional=True,
                  show_only=show_only)

def setup_package():
    if not 'sdist' in sys.argv[1:]:
        # Import setuptools to find a C compiler on windows.
        import setuptools; setuptools

    from numpy.distutils.core import setup

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, 'sfepy'))  # to retrive version

    # Write the version file.
    fd = open('VERSION', 'w')
    fd.write(VERSION)
    fd.close()

    # Create version.h file.
    filename_in = 'sfepy/discrete/common/extmods/version.h.in'
    filename_out = 'sfepy/discrete/common/extmods/version.h'
    fdi = open(filename_in, 'r')
    fdo = open(filename_out, 'w')
    for line in fdi:
        if line.find('VERSION "0.0.0"') >= 0:
            aux = line.split()
            aux[2] = VERSION
            line = ' '.join(aux) + '\n'
        fdo.write(line)
    fdi.close()
    fdo.close()

    install_requires = [
        'matplotlib',
        'meshio',
        'numpy',
        'pyparsing',
        'pyvista',
        'scipy',
        'sympy',
        'tables',
    ]

    try:
        setup(name='sfepy',
              maintainer="Robert Cimrman",
              maintainer_email="cimrman3@ntc.zcu.cz",
              description=DOCLINES[0],
              long_description="\n".join(DOCLINES[2:]),
              url="http://sfepy.org",
              download_url=DOWNLOAD_URL,
              license='BSD',
              classifiers=list(filter(None, CLASSIFIERS.split('\n'))),
              platforms=["Linux", "Mac OS-X", 'Windows'],
              entry_points={
                  'console_scripts': [
                      'sfepy-convert=sfepy.scripts.convert_mesh:main',
                      'sfepy-mesh=sfepy.scripts.gen_mesh:main',
                      'sfepy-probe=sfepy.scripts.probe:main',
                      'sfepy-run=sfepy.scripts.simple:main',
                      'sfepy-test=sfepy.scripts.run_tests:main',
                      'sfepy-view=sfepy.scripts.resview:main',
                  ],
              },
              install_requires=install_requires,
              cmdclass=cmdclass,
              configuration=configuration)

    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    check_versions()
    setup_package()

    from sfepy import Config
    site_config = Config()
    log.info('\nUsing Python {}.'.format(site_config.python_version()))

    log.info('\nRequired and optional packages found:\n')
    check_versions(show_only=True)
