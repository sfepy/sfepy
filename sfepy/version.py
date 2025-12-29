# SfePy version
__version__ = '2025.4'

# "Minimal" supported versions.
NUMPY_MIN_VERSION = '1.3'
SCIPY_MIN_VERSION = '0.7'
MATPLOTLIB_MIN_VERSION = '0.99.0'
PYPARSING_MIN_VERSION = '1.5.0'
PYTABLES_MIN_VERSION = '2.1.2'
MAYAVI_MIN_VERSION = '3.3.0'
SYMPY_MIN_VERSION = '0.7.3'
IGAKIT_MIN_VERSION = '0.1'
PETSC4PY_MIN_VERSION = '3.4'
SLEPC4PY_MIN_VERSION = '3.4'
MPI4PY_MIN_VERSION = '1.3.1'
PYMETIS_MIN_VERSION = '2014.1'
SCIKIT_UMFPACK_MIN_VERSION = '0.1'
MESHIO_MIN_VERSION='4.0.0'
PSUTIL_MIN_VERSION='1.0.0'
PYVISTA_MIN_VERSION='0.23.0'
OPT_EINSUM_MIN_VERSION='3.0.0'
JAX_MIN_VERSION='0.2.0'
DASK_MIN_VERSION='2.0.0'
PRIMME_MIN_VERSION='3.0.0'
OCT2PY_MIN_VERSION='5.6.0'
PYTHON_MUMPS_MIN_VERSION='0.0.3'
MUMPSPY_MIN_VERSION='1.2.0'
IPCTK_MIN_VERSION='1.3.0'
PYPARDISO_MIN_VERSION='0.4.6'

CYTHON_MIN_VERSION = '0.14.1'

def get_version(version=__version__):
    """
    Get sfepy version.

    Appends current git commit hash to `version` if `is_release` is set to True
    in the site configuration.
    """
    import os.path as op
    from sfepy.config import top_dir, site_config

    if not site_config.is_release():
        # Append current git commit hash to __version__.
        master = op.join(top_dir, '.git/refs/heads/master')
        if op.isfile(master):
            fd = open(master, 'r')
            version += '+git.%s' % fd.readline().strip()
            fd.close()

    return version

__version__ = get_version(__version__)
