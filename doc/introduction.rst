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
`NumPy/SciPy <http://scipy.org>`_ (an excellent collection of tools for
scientific computations in Python).

*SfePy* was originally developed as a flexible framework to quickly implement
and test the mathematical models developed during our various
research projects. It has evolved, however, to a
rather full-featured (yet small) finite element code with many weak forms
to build equations so there is a chance it might serve you as well.

New users should start by going through the :doc:`tutorial` and then the
more focused :doc:`primer`.

Features:

* solution of linear, nonlinear problems
* multi-platform (Linux, Mac OS X, Windows)
* collection of modules (a library):
    * FE engine, problem description facilities,
      interfaces to various solvers, postprocessing utilities
    * usable to build custom applications
* "black box" PDE solver:
    * no real programming involved
    * just prepare a problem description file (in Python!) and solve it
    * highly customizable behaviour (with a bit of coding)

To find more information regarding the code itself, go to http://sfepy.org
where you can find:

* releases
* mailing lists
* issue tracking
* git repository (bleeding edge code)
* further documentation, examples, and some of the research projects
  the code has been developed for.

To discuss in real time, join our IRC channel #sfepy at freenode.
