#!/usr/bin/env python
"""SfePy: Simple finite elements in Python

SfePy (simple finite elements in Python) is a software, distributed
under the BSD license, for solving systems of coupled partial
differential equations by the finite element method. The code is based
on NumPy and SciPy packages.
"""
import glob
import os

from skbuild import setup
from setuptools import find_packages

import sys
sys.path.append('./tools')
from build_helpers import INFO, cmdclass, logging, log, package_check

from sfepy import site_config, version

DOCLINES = __doc__.split("\n")

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
            return '0.0.0'

    except:
        return None


def _primme_version(pkg_name):
    import primme

    return '3.0.0'


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
    package_check('primme', INFO.PRIMME_MIN_VERSION, optional=True,
                  version_getter=_primme_version, show_only=show_only)
    package_check('oct2py', INFO.OCT2PY_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('mumpspy', INFO.MUMPSPY_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('mumps', INFO.PYTHON_MUMPS_MIN_VERSION, optional=True,
                  show_only=show_only)
    package_check('ipctk', INFO.IPCTK_MIN_VERSION, optional=True,
                  show_only=show_only)


def data_dir_walk(dir_name: str, prefix: str) -> list:
    """
    Generate instructions for setup() to add all files in a tree rooted at `dirname`
    as data_files.
    """
    data_files = []
    for root, dirs, files in os.walk(dir_name):
        full_paths = [os.path.join(root, fname) for fname in files]
        data_files.append((os.path.join(prefix, root), full_paths))

    return data_files


def compose_data_files() -> list:
    data_files = [
        ('sfepy', ['LICENSE', 'VERSION', 'site_cfg_template.py']),
    ]
    test_files = [('sfepy/tests', glob.glob('sfepy/tests/*.py'))]
    mesh_data_files = data_dir_walk('meshes', 'sfepy')
    example_files = data_dir_walk('examples', 'sfepy')

    return data_files + test_files + mesh_data_files + example_files


def cmake_bool(py_bool: bool) -> str:
    return "ON" if py_bool else "OFF"


def compose_cmake_args() -> list:
    cmake_args = [f'-DCMAKE_C_FLAGS={" ".join(site_config.compile_flags())}']

    # Debug flags are always added explicitly, so they won't be taken from cmake cache.
    debug_flags = set(site_config.debug_flags())
    cmake_args.append(f"-DDEBUG_FMF={cmake_bool('DEBUG_FMF' in debug_flags)}")
    cmake_args.append(f"-DDEBUG_MESH={cmake_bool('DEBUG_MESH' in debug_flags)}")

    # On Azure CI images, Ninja isn't found automaticvally. Since this is harmless
    # on other platforms, we specify it here.
    # See https://conda-forge.org/docs/maintainer/maintainer_faq.html#mfaq-windows-cmake
    cmake_args.extend(["-G", "Ninja"])

    return cmake_args


def setup_package():
    # Write the version file.
    fd = open('VERSION', 'w')
    fd.write(VERSION)
    fd.close()

    # Create version.h file.
    # There is probably a way to do it with CMake but we'll get to it later.
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
        'numpy<2',
        'pyparsing',
        'pyvista',
        'scipy',
        'sympy',
        'tables',
        'packaging',
    ]

    setup(
        name='sfepy',
        version=version.__version__,
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
        packages=find_packages(),
        data_files=compose_data_files(),
        setup_requires=['cython'],
        cmake_args=compose_cmake_args(),
        cmake_languages=('C')
    )


if __name__ == '__main__':
    check_versions()
    setup_package()

    logging.basicConfig(level=logging.INFO, force=True)

    log.info('Using Python {}.'.format(site_config.python_version()))

    log.info('Required and optional packages found:')
    check_versions(show_only=True)
