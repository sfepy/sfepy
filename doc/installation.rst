.. include:: links.inc

.. _introduction_installation:

Installation
============

.. only:: html

   .. contents:: Table of Contents
      :local:
      :backlinks: top

Supported Platforms
-------------------

*SfePy* is known to work on various flavors of recent Linux, Intel-based MacOS
and Windows. The release 2019.4 was the last with Python 2.7 support. *SfePy*
should work with any recent Python 3.x that is supported by `NumPy`_ and
`SciPy`_.

Note: Depending on Python installation and OS used, replacing ``python`` by
``python3`` might be required in all the commands below
(e.g. in :ref:`compilation`) in order to use Python 3.

.. _Python_distribution:

Notes on selecting Python Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is only matter of taste to use either native OS Python installation, `pip`,
or any other suitable distribution. On all supported platforms we could
recommend multi-platform scientific Python distributions `Anaconda`_ as
easy-to-use, stable and up-to-date Python distribution with all the required
dependencies (including pre-built `sfepy` package). See also `Notes on
Multi-platform Python Distributions`_ for further details.

On different platforms the following options can be recommended:

- **Linux**: `Anaconda`_, OS native installation, if available, or `pip`_.

- **macOS**: `Anaconda`_.

- **Windows**: free versions of commercial scientific Python distributions
  `Anaconda`_ or `Enthought Deployment Manager`_ . In addition a completely free
  open-source portable distribution `WinPython`_ can be used.

.. _installing_sfepy:

Installing SfePy
----------------

The released versions of SfePy can be installed as follows.

- Using `pip`_::

    pip install sfepy

- Using `Anaconda`_: install `sfepy` from `conda-forge`_ channel::

    conda install -c conda-forge sfepy

  See `Notes on Multi-platform Python Distributions`_ for additional notes.

If the installation succeeded, proceed with `Testing Installation`_.

.. _running_sfepy_docker_images:

Using SfePy Docker Images
---------------------------

Besides the classical installation we also provide official Docker images
with ready-to-run Anaconda and *SfePy* installation.

Before you start using *SfePy* images, you need to first install and configure
Docker on your computer. To do this follow official
`Docker documentation <https://docs.docker.com/get-docker/>`__.

Currently available all-in-one image is:

- `sfepy/sfepy-desktop
  <https://hub.docker.com/r/sfepy/sfepy-desktop>`__ -
  an Ubuntu based container containing a full desktop environment in
  officially supported flavors accessible via any modern web browser.

For available runtime options and further information see
`sfepy-docker <https://github.com/sfepy/sfepy-docker>`__ project on Github.

.. _installing_from_sources:

Installing SfePy from Sources
-----------------------------

The latest stable release can be obtained from the `download`_ page. Otherwise,
download the development version of the code from `SfePy git repository`_::

    git clone git://github.com/sfepy/sfepy.git

In case you wish to use a specific release instead of the latest master
version, use::

    git tag -l

to see the available releases - the release tags have form
`release_<year>.<int>`.

See the `download`_ page for additional download options.

Requirements
^^^^^^^^^^^^

Installation prerequisites, required to build *SfePy*:

- a C compiler suite,
- `Python`_ 3.x,
- `NumPy`_,
- `Cython`_.
- `Cmake`_
- `scikit-build`_
- `ninja`_

Python packages required for using *SfePy*:

- `Pyparsing`_,
- `SciPy`_,
- `meshio`_ for reading and writing mesh files,
- `scikit-umfpack`_ for enabling `UMFPACK`_ solver for SciPy >= 0.14.0,
- `Matplotlib`_ for various plots,
- `PyTables`_ for storing results in HDF5 files,
- `SymPy`_ for some tests and functions,
- `igakit`_ for generating IGA domains,
- `petsc4py`_ and `mpi4py`_ for running parallel examples and using parallel
  solvers from `PETSc`_,
- `slepc4py`_ for eigenvalue problem solvers from `SLEPc`_,
- `pymetis`_ for mesh partitioning using `Metis`_,
- `Read the Docs`_ `Sphinx`_ theme for building documentation,
- `psutil`_ for memory requirements checking,
- `PyVista`_ for post-processing.

Make sure the dependencies of those packages are also installed (e.g `igakit`_
reguires FORTRAN compiler, `scikit-umfpack`_ does not work without UMFPACK,
`petsc4py`_ without PETSc etc.). All dependencies of `meshio`_ need to be
installed for full mesh file format support (when using pip: ``pip install
meshio[all]``).

*SfePy* should work both with bleeding edge (Git) and last released versions of
`NumPy` and `SciPy`. Please, submit an issue at `Issues`_ page in case this
does not hold.

Other dependencies/suggestions:

- To be able to (re)generate the documentation `Sphinx`_, `numpydoc`_ and
  `LaTeX`_ are needed (see :ref:`how_to_regenerate_documentation`).
- If `doxygen` is installed, the documentation of data structures and functions
  can be automatically generated by running::

   python setup.py doxygendocs

- Mesh generation tools use `pexpect` and `gmsh` or `tetgen`.
- `IPython`_ is recommended over the regular Python shell to fluently follow
  some parts of primer/tutorial (see :ref:`using-ipython`).
- `MUMPS`_ library for using MUMPS linear direct solver
  (real and complex arithmetic, parallel factorization)

.. _compilation:

Compilation of C Extension Modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the *SfePy* top-level directory:

#. Look at ``site_cfg_template.py`` and follow the instructions
   therein. Usually no changes are necessary.

#. For in-place use, compile the extension modules::

     python setup.py build_ext --inplace

   After a successful compilation, *SfePy* can be used in-place. However, the
   the ``sfepy-*`` commands, such as ``sfepy-run`` are only available after
   installing the package. Their functionality can be accessed by invoking
   directly the corresponding scripts in ``sfepy/scripts/``.

Installation
^^^^^^^^^^^^

*SfePy* can be used without any installation by running its main scripts
and examples from the top-level directory of the distribution or can be
installed *locally* or *system-wide*:

- system-wide (may require root privileges)::

    pip install .

- local::

    pip install --user .

- development (editable install)::

    pip install -e .

  The editable install allows working in-place and at the same time the
  ``sfepy-*`` commands are available.

If all went well, proceed with `Testing Installation`_.

.. _testing_installation:

Testing Installation
--------------------

After building and/or installing *SfePy* you should check if all the
functions are working properly by running the automated tests.

Running Automated Test Suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The test suite is based on `pytest`_. Install it by::

  pip install pytest

or::

  conda install pytest

when working in `Anaconda`_. If *SfePy* was installed, it can be tested using
the command::

  sfepy-test

that accepts all of the `pytest`_ options, for example::

  sfepy-test -vv --durations=0 -m 'not slow' -k test_assembling.py

The tests output directory can also be specified::

  sfepy-test --output-dir=output-tests

In general. the tests can be run using::

  python -c "import sfepy; sfepy.test()"

in the *SfePy* top-level directory in case of the in-place build and anywhere
else when testing the installed package. Additional `pytest`_ options can be
passed as arguments to ``sfepy.test()``, for example::

  python -c "import sfepy; sfepy.test('-vv', '--durations=0', '-m not slow', '-k test_assembling.py')"

Analogously to ``sfepy-test``, the tests output directory can be specified
using::

  python -c "import sfepy; sfepy.test('--output-dir=output-tests')"

See `pytest usage instructions
<https://docs.pytest.org/en/latest/how-to/usage.html>`_ for other options and
usage patterns.

To test an in-place build (e.g. in a cloned git repository), the following
simpler command can be used in the sources top-level directory::

  python -m pytest sfepy/tests
  python -m pytest -v sfepy/tests/test_assembling.py

which will also add the current directory to ``sys.path``. If the top-level
directory is already in ``sys.path`` (e.g. using ``export PYTHONPATH=.``), the
simplest way of invoking pytest is::

  pytest sfepy/tests
  pytest -v sfepy/tests/test_assembling.py

Debugging
---------

If something goes wrong,  edit the ``site_cfg.py`` config file and set
``debug_flags = '-DDEBUG_FMF'``
to turn on bound checks in the low level C functions, and recompile the code::

    python setup.py clean
    python setup.py build_ext --inplace

Then re-run your code and report the output to the `SfePy mailing list`_.

.. _using-ipython:

Using IPython
-------------

We generally recommend to use (a customized) `IPython`_ interactive shell over
the regular Python interpreter when following :doc:`tutorial` or
:doc:`primer` (or even for any regular interactive work with *SfePy*).

Install `IPython`_ (as a generic part of your selected distribution) and then
customize it to your choice.

Depending on your IPython usage, you can customize your `default` profile or
create a *SfePy* specific new one as follows:

#. Create a new *SfePy* profile::

     ipython profile create sfepy

#. Open the ``~/.ipython/profile_sfepy/ipython_config.py`` file in a text
   editor and add/edit after the ``c = get_config()`` line:

   .. sourcecode:: python

      exec_lines = [
          'import numpy as nm',
          'import matplotlib as mpl',
          'mpl.use("WXAgg")',
      #
      # Add your preferred SfePy customization here...
      #
      ]

      c.InteractiveShellApp.exec_lines = exec_lines
      c.TerminalIPythonApp.gui = 'wx'
      c.TerminalInteractiveShell.colors = 'Linux' # NoColor, Linux, or LightBG

   Please note, that generally it is not recommended to use `star` (*)
   imports here.

#. Run the customized IPython shell::

     ipython --profile=sfepy


.. _multi_platform_distributions_notes:

Notes on Multi-platform Python Distributions
--------------------------------------------

Anaconda
^^^^^^^^

We highly recommend this scientific-oriented Python distribution.

(Currently regularly tested by developers on *SfePy* releases
with Python 3.6 64-bit on Ubuntu 16.04 LTS, Windows 8.1+ and macOS 10.12+.)

Download appropriate `Anaconda`_ Python 3.x installer package and follow
install instructions. We recommend to choose *user-level* install option (no
admin privileges required).

Anaconda can be used for:

#. installing the latest release of *SfePy*  directly from the `conda-forge`_
   channel (see `sfepy-feedstock`_). In this case, follow the instructions
   in `Installing SfePy`_.

   Installing/upgrading *SfePy*  from the conda-forge channel can also be
   achieved by adding `conda-forge`_ to your channels with::

     conda config --add channels conda-forge

   Once the `conda-forge`_ channel has been enabled, *SfePy* can be installed
   with::

     conda install sfepy

   It is possible to list all of the versions of *SfePy*  available on your
   platform with::

     conda search sfepy --channel conda-forge

#. installing the *SfePy* dependencies only - then proceed with the
   :ref:`installing_from_sources` instructions.

   In this case, install the missing/required packages using built-in `conda`
   package manager::

     conda install wxpython

   See `conda help` for further information.

Occasionally, you should check for distribution and/or installed packages
updates (there is no built-in automatic update mechanism available)::

  conda update conda
  conda update anaconda
  conda update <package>

or try::

  conda update --all


Compilation of C Extension Modules on Windows
"""""""""""""""""""""""""""""""""""""""""""""

To build *SfePy* extension modules, included `mingw-w32/64`_ compiler tools
should work fine. If you encounter any problems, we recommend to install and
use Microsoft `Visual C++ Build Tools`_ instead (see `Anaconda FAQ`_).
