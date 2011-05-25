Introduction
============

*SfePy* is a finite element analysis software written almost entirely in
`Python <http://python.org>`_, with exception of the most time demanding
routines - those are written in C and wrapped by `SWIG <http://www.swig.org>`_.

*SfePy* is a free software released under the `New BSD License
<http://www.opensource.org/licenses/bsd-license.php>`_.  It relies on
`NumPy/SciPy <http://scipy.org>`_ (an excellent collection of tools for scientific
computations in Python).

*SfePy* was originally developed as a flexible framework to quickly implement
and test the mathematical models developed within the `ZČU Projects
<http://sfepy.kme.zcu.cz/index.cgi/Projects>`_. It has evolved, however, to a
rather full-featured (yet small) finite element code with many weak forms
to build equations so there is a chance it might serve you as well.

New users should start by going through the :doc:`tutorial`.

To find more information regarding the code itself, go to http://sfepy.org
where you can find:

* releases
* mailing lists
* issue tracking
* git repository (bleeding edge code)
* further documentation and examples

To discuss in real time, join our IRC channel #sfepy at freenode.

.. _introduction_installation:

Installation
------------

Requirements
^^^^^^^^^^^^

Installation prerequisites:

* recent numpy, scipy (with umfpack wrapper, or umfpack scikit), swig 

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

(tested on Windows XP Pro 32-bit)

Here we provide instructions for using *SfePy* on Windows through
`Python(x,y)`_. We will also use
`msysgit <http://code.google.com/p/msysgit>`_ to install the umfpack scikit to
speed performance. 

You will probably need a few gigabytes of free disk space due to the installed
size of Python(x,y) and msysgit. 

This procedure was tested on a Windows XP 32-bit machine. It should work in
theory for any Windows machine supported by Python(x,y) and msysgit, but your
milage may vary.

There are many steps, but hopefully it is straightforward to follow this
procedure. If you have any questions or difficulties please feel free to ask on
the sfepy-devel mailing list (see http://sfepy.org). Also, if you have any
suggestions for improving or streamlining this process, it would be very
beneficial as well!

Steps to get a working *SfePy* on Windows using Python(x,y)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#. Download and install current release of Python(x,y) from the *Downloads* tab
   at http://www.pythonxy.com/

    * Version at time of writing is 2.6.5.1
    * Detailed Python(x,y) installation notes

        #. Double click on the Python(x,y) installation file
        #. Click *Okay* at warning if not installing with an administrator
           account
        #. Click *I Agree* at license notification window
        #. In the *Python* subtree of the *Choose components* window,
           additionally select the following packages:

            #. ETS
            #. SymPy

        #. You can optionally select the following packages:

            #. Sphinx - needed to rebuild the documentation
            #. winpdb - a graphical python debugger - useful for solving
               problems with *SfePy* runs

        #. You can optionally choose whether or not to install Eclipse - it is
           not needed by *SfePy*
        #. In the *Other* subtree of the *Choose components* window,
           additionally select SWIG
        #. Choose an installation directory

            * Probably ``C:\pythonxy`` is good unless you have a reason to change it

        #. It will take a few minutes to install Python(x,y)

#. Download and install msysgit from *Downloads* tab at
   http://code.google.com/p/msysgit/

    * Version at time of writing is 1.7.0.2
    * Make sure to get the file that begins with ``msysGit-fullinstall``
    * Detailed msysgit installation notes

        #. Double click on the file beginning with
           ``msysGit-fullinstall`` to start the installation process
    	#. Choose an installation directory

            * Probably the default, ``C:\msysgit``, is best unless you
              have a reason to change it

    #. After clicking ok, the msysgit files will be extracted and then a
       terminal window will open showing git being compiled
     
#. Open an msys terminal (or use the one that opens after installing msysgit)

    * The msys terminal is opened with ``C:\msysgit\msys.bat``

#. Install pyparsing using easy_install

    * In the msys terminal, type the following command::

          easy_install http://pypi.python.org/packages/source/p/pyparsing/pyparsing-1.5.2.tar.gz

        * Note: ``easy_install pyparsing`` should also work, but appears to have
          a problem so it may be better to enter the full URL as above

#. Decide where to put *SfePy* and UMFPACK scikit sources

    * Probably ``C:\src`` is good
    * In msys, this path would be ``/c/src``

#. Create the source directory and change to it

    * In the msys window, type the following commands::

        mkdir /c/src
        cd /c/src

#. Download the UMFPACK scikit source code

    * For this task, we will use the svn support built in to git
    * In the msys window, type the following commands::

        git svn clone http://svn.scipy.org/svn/scikits/trunk/umfpack umfpack-scikit

#. Create the source directory inside ``umfpack-scikit`` to store the source code
   for UMFPACK and AMD

    * In msys, type:: 
    
        mkdir umfpack-scikit/src
        cd umfpack-scikit/src

#. Download UMFPACK, AMD, and UFconfig source code

    * In msys, type the following commands::

        curl -O http://www.cise.ufl.edu/research/sparse/umfpack/current/UMFPACK.tar.gz
        curl -O http://www.cise.ufl.edu/research/sparse/amd/current/AMD.tar.gz
        curl -O http://www.cise.ufl.edu/research/sparse/UFconfig/current/UFconfig.tar.gz

#. Extract the UMFPACK, AMD, and UFconfig sources

    * In msys, type the following::

        tar zxf UMFPACK.tar.gz
        tar zxf AMD.tar.gz
        tar zxf UFconfig.tar.gz

#. Edit ``UFconfig.mk``

    * We need to set some configuration options in ``UFconfig/UFconfig.mk``
    * Use your favorite text editor to edit this file
    * Find the line that reads ``UMFPACK_CONFIG =``
    * Modify this line to the following:

        * ``UMFPACK_CONFIG = -DNCHOLMOD -DNBLAS``
        * Note: we are disabling BLAS and CHOLMOD to make it easier to compile
          UMFPACK. This may have some performance penalty associated with it. If
          you have experience compiling BLAS/LAPACK/ATLAS on Windows, please
          send us a message on the sfepy-devel mailing list!

#. Now change to the UMFPACK directory and make the library:

    * In msys, type::

        cd UMFPACK
        make library

#. Copy ``UFconfig.h`` to ``UMFPACK/Include``

    * In msys, type:: 
    
        cp ../UFconfig/UFconfig.h Include/

#. Now we need to make a ``site.cfg`` in the umfpack-scikit directory
   corresponding to our current setup

    * In msys, type::

        cd /c/src/umfpack-scikit
        cp site.cfg.example site.cfg

   * Using your favorite text editor, change the all the paths to point to the
     UMFPACK and AMD directories (non-msys paths)

     * E.g., ``include_dirs = /Users/stefan/src/UMFPACK/Include`` ->
       ``include_dirs = c:/src/umfpack-scikit/src/UMFPACK/Include``

#. Now it's time to install the UMFPACK scikit!

    * In msys, type::

        python setup.py install

    * Congratulations, you should now have a working UMFPACK scikit on Windows!

#. Decide which version of *SfePy* to use

    * If you want to use the stable released version, grab the tarball from the
      *Downloads* tab at and extract it in ``C:/src``

        * In msys, type:: 

            cd /c/src
            curl -O http://sfepy.googlecode.com/files/sfepy-release-2010.2.tgz
            tar zxf sfepy-release-2010.2.tgz

    * If you want to use the latest features and contribute to the development
      of *SfePy*, clone the git development repository

        * In msys, type::

            cd /c/src
            git clone git://github.com/sfepy/sfepy.git

#. Compile *SfePy* C extensions

   * In msys, change to the *SfePy* directory you downloaded in the preceding
     step with the ``cd`` command
   * Type:: 

       python setup.py build_ext --inplace --compiler=mingw32

 #. Run *SfePy* tests

    * Congratulations! You should (hopefully) now have a working copy of *SfePy*
      on Windows 
    * Please help aid *SfePy* development by running the builtin tests

        * In msys, in the *SfePy* source directory, type::

            ./runTests.py --filter-less

        * Report any failures to the sfepy-devel mailing list
        * See `Checking the SfePy installation`_ for further details

