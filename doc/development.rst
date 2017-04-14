.. include:: links.inc

Development
===========

The SfePy development takes place in the `sfepy/sfepy`_ repository on Github.
The users and developers can also communicate using the `mailing list`_.

General Information
-------------------

We are interested in any contribution. There are many ways how you can
contribute:

- You can report bugs using `our mailing list`_. You can also add a bug
  report (or a comment) into the `issues`_.
- You can contribute interesting examples/tutorials.
- You can blog about how you use SfePy (let us know!).
- You can help with improving our documentation and these pages.
- ...

To get acquainted with SfePy, you can start by reading the :ref:`sec-tutorial`
and :ref:`sec-primer` sections of the documentation and trying out the examples
that come with the sources. Your first contribution could be pinpointing
anything that is not clear in the docs.

We also recommend reading the :ref:`how_to_contribute` section of our
:ref:`developer_guide`.

Topics
------

Several specific topics that we wish to address in the future are listed below.
If you are a domain expert in some of those, and would like to contribute
code/advice to our project, do not hesitate to contact us (either directly:
cimrman3(at)ntc.zcu.cz, or on `our mailing list`_)

- finish/improve IGA implementation (see :ref:`isogeometric_analysis`):

  - support multiple patches
  - efficient quadrature formulas
  - local refinement?

- discretization methods:

  - implement vector elements (Nedelec, Raviart-Thomas, ...)
  - implement the discontinuous Galerkin method

- improve parallelization (see :ref:`solving_problems_in_parallel`):

  - cluster installation with fast BLAS
  - parallel code speed-up
  - remove (some of) the serial parts
  - preconditioning for multi-physics problems

- better defaults/recommendations for iterative solvers (PETSc) with respect
  to large problems

- material models: plasticity, viscoplasticity, damage, ...

- automatic differentiation:

  - for tangent matrices
  - for identification of (material) parameters

- core data structures & programming:

  - using octree-based(?) mesh representation for local refinement
  - continue with/improve the current hanging nodes implementation
  - exploit lazy evaluation

- visualization of large data

See also the `enhacement issues <https://github.com/sfepy/sfepy/issues?q=is%3Aissue+is%3Aopen+label%3Aenhancement>`_.
