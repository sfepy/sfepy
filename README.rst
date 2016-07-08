========================================
SfePy (Simple finite elements in Python)
========================================

A finite element analysis software based primarily on NumPy and SciPy.

See the INSTALL file for installation instructions or refer to SfePy
documentation site [1].

How to start
------------

SfePy uses so called "problem definition files" (also referred to as
"input files" or "problem description files"") that describe the partial
differential equations (PDEs), boundary conditions, function spaces and
other ingredients of the finite element (FE) formulation of a
PDE-related problem, see [1]. The PDEs are given in weak formulation as
usual in the FE context, see [2], where each equation is composed of one
or more terms. To see which terms are available consult
'doc/sfepy_manual.pdf', or see them online [3].

In order to solve a problem, a problem description file has to be
created. There is also an interactive solution support for advanced
users [4].

When starting to solve a new problem, it is best to have a look at
example problem definition files in the 'examples/' directory - copy the
one that is similar to the problem at hand, and modify it. Two of the
examples are commented: 'examples/diffusion/poisson.py' and
'examples/linear_elasticity/linear_elastic.py'.

While a problem definition file describes a mathematical problem, it does not
contain a discretized solution domain (a FE mesh). The FE mesh must be provided
in another file in one of the supported formats, notably the legacy VTK format
[5]. SfePy does not provide meshing tools, but it can use a number of standard
formats. The results are almost exclusively stored in legacy VTK files, or
custom HDF5 files. Many standard open-source tools can be used to display the
VTK files, namely paraview [6], or mayavi [7]. The latter is supported directly
within SfePy, via the postproc.py script.

Once an input file and a corresponding mesh file are prepared, the
solution of the problem can be attempted by the 'simple.py' script, see
Examples section in the INSTALL file.

[1] http://docs.sfepy.org/
[2] http://en.wikipedia.org/wiki/Weak_formulation
[3] http://docs.sfepy.org/doc/terms_overview.html
[4] http://docs.sfepy.org/doc/tutorial.html#interactive-example-linear-elasticity
[5] http://www.vtk.org/pdf/file-formats.pdf
[6] http://paraview.org/
[7] http://code.enthought.com/projects/mayavi/

Links
-----

 - http://sfepy.org ... main SfePy development site (releases, mailing lists,
                        wiki, issue tracking, downloads, examples
 - http://docs.sfepy.org ... documentation
 - http://github.com/sfepy/ ... master git repository
 - http://sfepy.kme.zcu.cz ... projects solved within SfePy

License: New BSD License, see the LICENSE file.

--
Robert Cimrman and Contributors
