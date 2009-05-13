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

VERSION = '2009.2-release'

CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: MacOS :: MacOS X
"""

DOWNLOAD_URL = "http://code.google.com/p/sfepy/wiki/Downloads?tm=2"

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

os.environ['NO_SFEPY_IMPORT']='SfePy/setup.py'

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
        'convert_doc.py',
        'convert_mesh.py',
        'edit_identifiers.py',
        'edit_neu.py',
        'evalForms.py',
        'eval_tl_forms.py',
        'genDocsXML.py',
        'hfm3_mesh.py',
        'mesh_to_vtk.py',
        'neu_mesh.py',
        'sfepyconverter.py',
        'spymatrix.py',
    ]
    aux_scripts = [os.path.join('script', ii) for ii in aux_scripts]

    config.add_data_files(('sfepy', ('VERSION', 'INSTALL', 'README', 'LICENSE',
                                     'RELEASE_NOTES.txt',
                                     'site_cfg_template.py')))
    config.add_data_files(('../../../share/sfepy', aux_scripts))

    config.get_version('sfepy/version.py') # sets config.version
    print config
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

    main_scripts = [
        'eigen.py',
        'extractor.py',
        'findSurf.py',
        'gen',
        'genDocs.py',
        'genPerMesh.py',
        'isfepy',
        'postproc.py',
        'probe.py',
        'runTests.py',
        'schroedinger.py',
        'simple.py'
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
              platforms = ["Linux", "Mac OS-X"],
              scripts = main_scripts,
#              cmdclass = {'install_scripts' : install_scripts},
              configuration = configuration)
    finally:
        del sys.path[0]
        os.chdir(old_path)

    return

if __name__ == '__main__':
    setup_package()
