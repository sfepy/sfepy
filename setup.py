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

VERSION = '2010.3'

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

DOWNLOAD_URL = "http://code.google.com/p/sfepy/wiki/Downloads?tm=2"

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('sfepy')
    aux_scripts = [
        'blockgen.py',
        'config.py',
        'convert.py',
        'convert_mesh.py',
        'cylindergen.py',
        'edit_identifiers.py',
        'edit_neu.py',
        'evalForms.py',
        'eval_tl_forms.py',
        'gen_term_table.py',
        'hfm3_mesh.py',
        'mesh_to_vtk.py',
        'neu_mesh.py',
        'spymatrix.py',
    ]
    aux_scripts = [os.path.join('script', ii) for ii in aux_scripts]

    config.add_data_files(('sfepy', ('VERSION', 'INSTALL', 'README', 'LICENSE',
                                     'RELEASE_NOTES.txt', 'AUTHORS',
                                     'site_cfg_template.py')))
    config.add_data_files(('../../../share/sfepy/script', aux_scripts))
    config.add_data_dir(('../../../share/sfepy/meshes', 'meshes'))
    config.add_data_dir(('../../../share/sfepy/examples', 'examples'))
    config.add_data_dir(('../../../share/sfepy/doc', 'doc'))
    config.add_data_dir(('../../../share/sfepy/tests', 'tests'))

    config.get_version('sfepy/version.py') # sets config.version
    ## print config

    return config

def setup_package():

    from numpy.distutils.core import setup
    from numpy.distutils.misc_util import Configuration

    old_path = os.getcwd()
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)
    sys.path.insert(0, os.path.join(local_path, 'sfepy')) # to retrive version

    # Write the version file.
    fd = open('VERSION', 'w')
    fd.write(VERSION)
    fd.close()

    # Create version.h file.
    filename_in = 'sfepy/fem/extmods/version.h.in'
    filename_out = 'sfepy/fem/extmods/version.h'
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

    main_scripts = [
        'eigen.py',
        'extractor.py',
        'findSurf.py',
        'genPerMesh.py',
        'homogen.py',
        'isfepy',
        'postproc.py',
        'probe.py',
        'runTests.py',
        'schroedinger.py',
        'sfepy_gui.py',
        'simple.py',
        'site_cfg_template.py',
    ]

    try:
        setup(name = 'sfepy',
              maintainer = "Robert Cimrman",
              maintainer_email = "cimrman3@ntc.zcu.cz",
              description = DOCLINES[0],
              long_description = "\n".join(DOCLINES[2:]),
              url = "http://sfepy.org",
              download_url = DOWNLOAD_URL,
              license = 'BSD',
              classifiers = filter(None, CLASSIFIERS.split('\n')),
              platforms = ["Linux", "Mac OS-X", 'Windows'],
              scripts = main_scripts,
#              cmdclass = {'install_scripts' : install_scripts},
              configuration = configuration)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
