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
from setuptools import find_packages, Extension
from Cython.Build import cythonize
import numpy as np

import sys
sys.path.append('./tools')
from build_helpers import INFO, cmdclass

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
"""
print(find_packages())
test_ext = Extension(
    'sfepy.discrete.common.extmods.mappings',
    sources=['sfepy/discrete/common/extmods/mappings*.pyx'],
    include_dirs=[np.get_include()]
)
print(test_ext)
"""

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


def data_dir_walk(dir_name: str, prefix: str) -> list[tuple[str, list[str]]]:
    """
    Generate instructions for setup() to add all files in a tree rooted at `dirname`
    as data_files.
    """
    data_files = []
    for root, dirs, files in os.walk(dir_name):
        full_paths = [os.path.join(root, fname) for fname in files]
        data_files.append((os.path.join(prefix, root), full_paths))

    return data_files


mesh_data_files = data_dir_walk('meshes', 'sfepy')

setup(
    name='sfepy',
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
    data_files=[
        ('sfepy', ['LICENSE', 'VERSION']),
        ('sfepy/tests/', glob.glob('sfepy/tests/*.py'))
    ] + mesh_data_files,
    setup_requires=['cython'],
)