.. include:: links.inc

Possible Topics
===============

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
