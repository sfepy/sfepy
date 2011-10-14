Introduction
============

.. contents:: Table of Contents
   :local:
   :backlinks: top

*SfePy* is a finite element analysis software written almost entirely in
`Python <http://python.org>`_, with exception of the most time demanding
routines - those are written in C and wrapped by `Cython
<http://cython.org/>`_ or written directly in Cython.

*SfePy* is a free software released under the `New BSD License
<http://www.opensource.org/licenses/bsd-license.php>`_.  It relies on
`NumPy/SciPy <http://scipy.org>`_ (an excellent collection of tools for scientific
computations in Python).

*SfePy* was originally developed as a flexible framework to quickly implement
and test the mathematical models developed during our various
research projects. It has evolved, however, to a
rather full-featured (yet small) finite element code with many weak forms
to build equations so there is a chance it might serve you as well.

New users should start by going through the :doc:`tutorial` and then the
more focused :doc:`primer`.

To find more information regarding the code itself, go to http://sfepy.org
where you can find:

* releases
* mailing lists
* issue tracking
* git repository (bleeding edge code)
* further documentation, examples, and some of the research projects
  the code has been developed for.

To discuss in real time, join our IRC channel #sfepy at freenode.

.. _introduction_installation:

Installation
------------

Requirements
^^^^^^^^^^^^

Installation prerequisites:

* recent numpy, scipy (with umfpack wrapper, or umfpack scikit), cython

Dependencies:

* matplotlib, pyparsing, umfpack, pytables
* some tests and functions use sympy
* schroedinger.py requires pysparse, pexpect, gmsh (2D), tetgen (3D)
* log.py (live plotting) requires multiprocessing, matplotlib with GTKAgg
* isfepy requires ipython, matplotlib with WXAgg
* postproc.py requires mayavi2 
* to be able to (re)generate the documentation: sphinx, numpydoc, LaTeX, see
  :ref:`how_to_regenerate_documentation`

*SfePy* is known to work on various flavours of Linux, on Intel Macs and Windows.

On Linux, consult the package manager of your favourite distribution. For
example in Gentoo::

    $ emerge -va pytables pyparsing numpy scipy matplotlib ipython 

in Debian::

    $ apt-get install python-tables python-pyparsing python-matplotlib python-scipy 

On Windows, all the required packages are part of the `Enthought Python
Distribution (EPD) <http://www.enthought.com/products/epd.php>`_, which is free
for academic purposes. A completely free `Python(x,y)
<http://www.pythonxy.com>`_ can be used too, but pyparsing has to
be installed manually. Instructions for installing Python(x,y) can be found in
`Running on Windows using Python(x,y)`_.

*SfePy* can be used without any installation by running the scripts from the
top-level directory of the distribution (TOPDIR), or can be installed locally or
system-wide.

*SfePy* should work both with bleeding edge (SVN) and last released versions of
its dependencies; see INSTALL file in the tarball. Submit an `issue
<http://code.google.com/p/sfepy/issues/entry>`_ in
case this does not hold.

Generic Installation Instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the latest source release or the development version from our git
repository http://github.com/sfepy/sfepy::

    git clone git://github.com/sfepy/sfepy.git

See the *Downloads* tab at http://sfepy.org for additional download options.

In-place compilation of C extension modules
"""""""""""""""""""""""""""""""""""""""""""

1. Look at site_cfg_template.py and follow the instructions
   therein. Usually no changes are necessary.

2. Run::

    python setup.py build_ext --inplace

Installation
""""""""""""

(As mentioned above, this step is not required to use *SfePy*.)

* System-Wide (may require root privileges)::

    python setup.py install

* Local (requires write access to ``<installation prefix>``)::

    python setup.py install --root=<installation prefix>

See also INSTALL and RELEASE_NOTES files in the tarball.

If all went well, proceed with `Checking the SfePy installation`_.

Checking the SfePy installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After installing *SfePy* you can check if all the functionalities are working by
running the automated tests. From the source directory type::

    ./runTests.py

If a particular test fails, please run it in debug mode::

    ./runTests.py --debug tests/failing_test_name.py

and report the output to the sfepy-devel mailing list.

Platform-specific notes
^^^^^^^^^^^^^^^^^^^^^^^

Fedora 8
""""""""
Notes on using umfpack (contributed by David Huard).

entry in numpy site.cfg::

    [umfpack]
    library_dirs=/usr/lib64
    include_dirs = /usr/include/suitesparse

Comment by david.huard, Mar 26, 2008::

> Of course, suitesparse and suitesparse-devel must be installed. 

Intel Mac
"""""""""

(thanks to Dominique Orban for his advice)

To build *SfePy* on an Intel Mac the following options need to be set in
site_cfg.py::

    opt_flags = '-g -O2 -fPIC -DPIC -fno-strict-aliasing -fno-common -dynamic' 
    link_flags = '-dynamiclib -undefined dynamic_lookup -fPIC -DPIC' 

Installation on Ubuntu
""""""""""""""""""""""

(tested on Jaunty Jackalope 9.04 and Lucid Lynx 10.04)

Prerequisites
+++++++++++++

First, you have to install the dependencies packages::

    sudo aptitude install python-scipy python-matplotlib python-tables
    python-pyparsing libsuitesparse-dev python-setuptools

Then download and install the umfpack scikits in some local dir. In the
following example it will be installed in $HOME/local::

    svn checkout http://svn.scipy.org/svn/scikits/trunk/umfpack
    cd umfpack
    mkdir -p ${HOME}/local/lib/python2.6/site-packages
    python setup.py install --prefix=${HOME}/local

Add to your .bashrc the line::

    export PYTHONPATH="${HOME}/local"

then re-open a terminal and if the scikits was installed correctly importing
scikits.umfpack in python should give no error::

    $ python
    >>> import scikits.umfpack
    >>>

Next Download sympy 6.7 or later. Extract the contents.

cd sympy-0.6.7

python setup.py install --prefix=${HOME}/local

Installing SfePy
++++++++++++++++

Now download the latest *SfePy* tarball release (or the latest development
version).

Uncompress the archive and enter the *SfePy* dir, then type::

    python setup.py build_ext --inplace

after a few minutes the compilation finishes.

Finally you can test *SfePy* with::

    ./runTests.py

If some test fails see `Checking the SfePy installation`_ section for further
details.


Running on Windows using Python(x,y)
""""""""""""""""""""""""""""""""""""

(tested on Windows 7)

Here we provide instructions for using *SfePy* on Windows through
`Python(x,y)`_. We will also use
`msysgit <http://code.google.com/p/msysgit>`_ to install the umfpack scikit to
speed performance.

This procedure was tested on a Windows 7 machine. It should work in
theory for any Windows machine supported by Python(x,y) and msysgit, but your
milage may vary.

There several steps, but hopefully it is straightforward to follow this
procedure. If you have any questions or difficulties please feel free to ask on
the sfepy-devel mailing list (see http://sfepy.org). Also, if you have any
suggestions for improving or streamlining this process, it would be very
beneficial as well!

We assume the installation to be done in C:/ - substitute your path
where appropriate.

Steps to get a working *SfePy* on Windows using Python(x,y)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#. Minimum 4 Gigabytes of free disk space is required, Due to the
   installed size of python(x,y) and msysgit.

#. Download the latest python(x,y) (http://www.pythonxy.com/) windows
   installer, and make a *Full installation* in the default installation
   directory.

#. Download the latest pyparsing
   (http://sourceforge.net/projects/pyparsing) windows installer (Python
   version 2.6) and install it in the default installation directory.

#. Download the latest mysysgit (http://code.google.com/p/msysgit/)
   windows installer (get the file that begins with
   "msysGit-fullinstall") and install it in the default installation
   directory.

#. Download the latest umfpackpy (http://code.google.com/p/umfpackpy/)
   zip archive and follow the instructions below:

   a) Extract the *umfpackpy_0.1.zip* to your convenient location in
      Hard disk, Lets assume it's extracted in C:/. Now there will be
      two files on the extracted folder, *ez_setup.py* and
      *scikits.umfpack-5.1.0-py2.6-win32.egg*.

   b) Start msys and write the following code to go to the extracted folder::

      $ cd /c/umfpackpy_0.1/

   c) Write the following code on msys to install the UMFPACK library
      for python::

      $ ez_setup.py scikits.umfpack-5.1.0-py2.6-win32.egg

#. Either download the latest sfepy (http://code.google.com/p/sfepy/)
   tarball and extract it to your convenient location in Hard disk,
   Lets assume it's extracted in C:/.

   Or, If you want to use the latest features and contribute to the
   development of *SfePy*, clone the git development repository

      * In msys, type::

          cd /c/
          git clone git://github.com/sfepy/sfepy.git

   Then follow the instructions below::

   a) Start msys and write the following code to go to the extracted folder::

      $ cd /c/sfepy_folder_name/

   b) Write the following code on msys to Compile SfePy C extensions::

      $ python setup.py build_ext --inplace --compiler=mingw32

#. You should now have a working copy of SfePy on Windows, Please help
   aid SfePy development by running the built-in tests. Run the
   *runTests.py* in python IDLE or Write the following code on msys::

     $ runTests.py --filter-less

   * Report any failures to the sfepy-devel mailing list
   * See `Checking the SfePy installation`_ for further details.
