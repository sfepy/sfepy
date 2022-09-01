.. image:: https://travis-ci.org/sfepy/sfepy.svg?branch=master
    :target: https://travis-ci.org/sfepy/sfepy
    :alt: CI

========================================
SfePy (Simple finite elements in Python)
========================================

SfePy [1]_ is a software for solving systems of coupled partial differential
equations (PDEs) by the finite element method in 1D, 2D and 3D. It can be
viewed both as black-box PDE solver, and as a Python package which can be used
for building custom applications. The word "simple" means that complex FEM
problems can be coded very easily and rapidly. Its source code can be found on
github [2]_.

SfePy is based primarily on NumPy [3]_ and SciPy [4]_. See the INSTALL file for
installation instructions or refer to [1]_.

- License: New BSD License, see the LICENSE file.

- Authors: Robert Cimrman and Contributors, see the AUTHORS file.

How to start
------------

SfePy uses so called "problem definition files" (also referred to as "input
files" or "problem description files"") that describe the partial differential
equations (PDEs), boundary conditions, function spaces and other ingredients of
the finite element (FE) formulation of a PDE-related problem, see [1]_. The
PDEs are given in weak formulation as usual in the FE context, see [5]_, where
each equation is composed of one or more terms. To see which terms are
available consult ``doc/sfepy_manual.pdf``, or see them online [6]_.

In order to solve a problem, a problem description file has to be created.
There is also an interactive solution support for advanced users [7]_. When
starting to solve a new problem, it is best to have a look at example problem
definition files in the ``sfepy/examples/`` directory - copy the one that is
similar to the problem at hand, and modify it.

While a problem definition file describes a mathematical problem, it does not
contain a discretized solution domain (a FE mesh). The FE mesh must be provided
in another file in one of the supported formats, notably the legacy VTK format
[8]_. SfePy does not provide meshing tools, but it can use a number of standard
formats. The results are almost exclusively stored in legacy VTK files, or
custom HDF5 files. Many standard open-source tools can be used to display the
VTK files, namely paraview [9]_ or pyvista [10]_. The latter is supported
directly within SfePy.

Once an input file and a corresponding mesh file are prepared, the solution of
the problem can be sought, see the documentation.

References
----------

.. [1] https://sfepy.org - main SfePy development site (releases, mailing lists,
       wiki, issue tracking, downloads, documentation, examples)
.. [2] https://github.com/sfepy - master git repository
.. [3] http://numpy.org
.. [4] http://scipy.org
.. [5] http://en.wikipedia.org/wiki/Weak_formulation
.. [6] https://docs.sfepy.org/doc/terms_overview.html
.. [7] https://docs.sfepy.org/doc/tutorial.html#interactive-example-linear-elasticity
.. [8] http://www.vtk.org/VTK/img/file-formats.pdf
.. [9] http://paraview.org/
.. [10] https://docs.pyvista.org/
