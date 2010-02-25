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
* schroedinger.py requires pysparse, pexpect, gmsh (2D), tetgen (3D)
* log.py (live plotting) requires multiprocessing, matplotlib with GTKAgg
* isfepy requires ipython, matplotlib with WXAgg
* postproc.py requires mayavi2 

*SfePy* is known to work on various flavours of Linux, on Intel Macs and Windows.

On Linux, consult the package manager of your favourite distribution. For
example in Gentoo::

    $ emerge -va pytables pyparsing numpy scipy matplotlib ipython 

in Debian::

    $ apt-get install python-tables python-pyparsing python-matplotlib python-scipy 

On Windows, all the required packages are part of the `Enthought Python
Distribution (EPD) <http://www.enthought.com/products/epd.php>`_, which is free
for academic purposes. A completely free `Python(x,y)
<http://www.pythonxy.com/foreword.php>`_ can be used too, but pyparsing has to
be installed manually.

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
repository http://git.sympy.org/?p=sfepy.git::

    git clone git://git.sympy.org/sfepy.git

See the *Downloads* tab at http://sfepy.org for additional download options.

In-place compilation of C extension modules
"""""""""""""""""""""""""""""""""""""""""""

There are two methods of compiling the C extension modules.

1. Using the Makefile (Linux, Mac OS X)

    1. Look at site_cfg_template.py and follow the instructions therein.
    2. (Optionally) edit the Makefile:

        * compiler, swig command, ... 

    3. ::
    
        make 

2. Using distutils (all platforms)::

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

    ./runTests.py --debug test/failing_test_name.py

and report the output to the sfepy-devel mailing list.

Platform-specific notes
^^^^^^^^^^^^^^^^^^^^^^^

Using umfpack on fedora 8
"""""""""""""""""""""""""

(contributed by David Huard)

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

(revision http://hg.sympy.org/sfepy/rev/609196c918be is needed) 

Installation on Ubuntu (tested on Jaunty Jackalope 9.04)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Prerequisites
+++++++++++++

First, you have to install the dependencies packages::

    sudo aptitude install python-scipy python-matplotlib python-tables
    python-pyparsing libsuitesparse-dev 

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

Installing SfePy
++++++++++++++++

Now download the latest *SfePy* tarball release (or the latest development
version).

Uncompress the archive and enter the *SfePy* dir, then type::

    make

after a few minutes the compilation finishes.

Finally you can test *SfePy* with::

    ./runTests.py

If some test fails see `Checking the SfePy installation`_ section for further
details.

