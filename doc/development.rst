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
If you would like to contribute code/advice to our project with respect to
these topics, do not hesitate to contact us (either directly:
cimrman3(at)ntc.zcu.cz, or on `our mailing list`_)

- finish/improve IGA implementation (see :ref:`isogeometric_analysis`):

  - support multiple patches
  - efficient quadrature formulas
  - local refinement?

- discretization methods:

  - implement vector elements (Nedelec, Raviart-Thomas, ...)
  - implement the discontinuous Galerkin method

- material models: plasticity, viscoplasticity, damage, ...
- improve parallelization (see :ref:`solving_problems_in_parallel`):

  - cluster installation with fast BLAS
  - parallel code speed-up
  - remove (some of) the serial parts
  - preconditioning for multi-physics problems

- solvers:

  - better defaults/recommendations for iterative solvers (`PETSc`_) with
    respect to large problems
  - dynamics/time-stepping solvers, interface PETSc time-steppers
  - interface more sparse linear solvers (or enable via PETSc), for example
    `BDDCML`_
  - interface more eigenvalue problem solvers

- visualization of large data
- automatic differentiation:

  - for tangent matrices
  - for identification of (material) parameters

- core data structures & programming:

  - using octree-based(?) mesh representation for local refinement
  - continue with/improve the current hanging nodes implementation
  - exploit lazy evaluation

See also the `enhacement issues <https://github.com/sfepy/sfepy/issues?q=is%3Aissue+is%3Aopen+label%3Aenhancement>`_.

Developer Guide
---------------

.. toctree::
   :maxdepth: 2

   developer_guide
