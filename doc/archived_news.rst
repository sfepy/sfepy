.. _archived_news:

Archived News
=============

* **29.09.2023** Version 2023.3 released (example demonstrating nonlinear
  diffusion and volume forces terms, new sensitivity term, proof-of-concept
  JAX-based auto-differentiable terms (elasticity, mass) + elastodynamics
  parameter identification example, etc.), see :ref:`release notes
  <2023.2-2023.3>`.

* **29.06.2023** Version 2023.2 released (replace numpy.distutils by
  scikit-build (cmake), interface for the PRIMME eigenvalue problem solver, new
  declarative modal analysis example, new terms: nonlinear diffusion and
  nonlinear volume forces, etc.), see :ref:`release notes <2023.1-2023.2>`.

* **29.03.2023** Version 2023.1 released (support meshes with cells of
  different topological dimensions, unified treatment of volume and surface
  integrals, updated elastodynamics solvers, reciprocal mass matrix algorithm,
  seismic load and piezo-elastodynamics examples, use GitHub Actions for CI,
  etc.), see :ref:`release notes <2022.4-2023.1>`.

* **29.12.2022** Version 2022.4 released (new examples: transient heat equation
  with several material domains and scalar Helmholtz equation, adaptive time
  step control for elastodynamics solvers, central difference elastodynamics
  solver, etc.), see :ref:`release notes <2022.3-2022.4>`.

* **29.09.2022** Version 2022.3 released (reorganized scripts, PyPI support,
  speed-up of linear elastodynamics problems, improved homogenization-based
  recovery of micro-scale solutions, etc.), see :ref:`release notes
  <2022.2-2022.3>`.

* **29.06.2022** Version 2022.2 released (custom testing code replaced by
  pytest, improved pyvista-based visualization script, etc.), see :ref:`release
  notes <2022.1-2022.2>`.

* **30.03.2022** Version 2022.1 released (new handling of state variables data
  and State class removal, many new sensitivity analysis terms based on
  multi-linear term implementation, etc.), see :ref:`release notes
  <2021.4-2022.1>`.

* **29.12.2021** Version 2021.4 released (improved pyvista-based visualization
  script resview.py, gallery images generated using resview.py, homogenization
  tools: new parallel recovery of multiple microstructures, new "dry water"
  flow example, etc.), see :ref:`release notes <2021.3-2021.4>`.

* **29.09.2021** Version 2021.3 released (unified volume and surface
  integration terms, improved pyvista-based visualization script, etc.), see
  :ref:`release notes <2021.2-2021.3>`.

* **29.06.2021** Version 2021.2 released (new sensitivity analysis terms,
  positive FE basis based on Bernstein polynomials, smaller memory footprint of
  terms with constant material parameters, etc.), see :ref:`release notes
  <2021.1-2021.2>`.

* **29.03.2021** Version 2021.1 released (non-square homogenized coefficient
  matrices, new implementation of multi-linear terms, improved handling of
  Dirichlet and periodic boundary conditions in common nodes, terms in the term
  table document linked to examples, etc.), see :ref:`release notes
  <2020.4-2021.1>`.

* **27.12.2020** Version 2020.4 released (Ogden hyperelastic term, serendipity
  finite element basis of orders 1-3, etc.), see :ref:`release notes
  <2020.3-2020.4>`.

* **29.09.2020** Version 2020.3 released (new script for visualizations based
  on pyvista, generalized Yeoh hyperelastic term + example, etc.), see
  :ref:`release notes <2020.2-2020.3>`.

* **30.06.2020** Version 2020.2 released (discontinuous Galerkin method
  implementation and examples, new website look, memory usage
  improvements, etc.), see :ref:`release notes <2020.1-2020.2>`.

* **31.03.2020** Version 2020.1 released (reading/writing of additional mesh
  formats by using `meshio <https://github.com/nschloe/meshio>`_), Python 3
  only from now on, etc.), see :ref:`release notes <2019.4-2020.1>`.

* **23.12.2019** Version 2019.4 released (support surface terms in 1D problems,
  improved Gmsh mesh format support (write .msh files), new updating procedure
  in nonlinear homogenization, improved/faster log plotter, etc.), see
  :ref:`release notes <2019.3-2019.4>`.

* **24.09.2019** Version 2019.3 released (interface to eigenvalue problem
  solvers in SLEPc, new Python 3 enabled Timer class and other Python 3
  compatibility fixes, etc.), see :ref:`release notes <2019.2-2019.3>`.

* **02.07.2019** Version 2019.2 released (improved support for time-dependent
  homogenization problems, Python 3.7 compatibility, etc.), see :ref:`release
  notes <2019.1-2019.2>`.

* **26.03.2019** Version 2019.1 released (automatic fallback for linear
  solvers, quadratic eigenvalue problem solver, etc.), see :ref:`release notes
  <2018.4-2019.1>`.

* **27.12.2018** Version 2018.4 released (better support for eigenvalue
  problems, improved MUMPS solver interface, support for logging and plotting
  of complex values, etc.), see :ref:`release notes <2018.3-2018.4>`.

* **17.09.2018** Version 2018.3 released (easier setting of values of
  variables, new script for outline edge extraction, new example:
  homogenization of a piezoelectric heterogeneous structure, etc.), see
  :ref:`release notes <2018.2-2018.3>`.

* **19.06.2018** Version 2018.2 released (generalized-alpha and velocity Verlet
  elastodynamics solvers, terms for dispersion in fluids, caching of reference
  coordinates for faster repeated use of probes, new wrapper of MUMPS linear
  solver for parallel runs, etc.), see :ref:`release notes <2018.1-2018.2>`.

* **06.03.2018** Version 2018.1 released (major update of time-stepping solvers
  and solver handling, Newmark and Bathe elastodynamics solvers, interface to
  MUMPS linear solver, new examples: iron plate impact problem
  (elastodynamics), incompressible Mooney-Rivlin material model
  (hyperelasticity) as a script, etc.), see :ref:`release notes
  <2017.4-2018.1>`.

* **29.12.2017** Version 2017.4 released (basic support for penalty-based
  contacts, support for user-defined contexts in all solvers and
  preconditioners, new example: dispersion analysis of heterogeneous periodic
  materials, etc.), see :ref:`release notes <2017.3-2017.4>`.

* **19.09.2017** Version 2017.3 released (support preconditioning in SciPy and
  PyAMG based linear solvers, user-defined preconditioners for PETSc linear
  solvers, parallel multiscale (macro-micro) homogenization-based computations,
  improved tutorial and installation instructions, etc.), see :ref:`release
  notes <2017.2-2017.3>`.

* **19.05.2017** Version 2017.2 released (simplified and unified implementation
  of some homogenized coefficients, support for saving custom structured data to
  HDF5 files, new tutorial on preparing meshes using FreeCAD/OpenSCAD and Gmsh,
  etc.), see :ref:`release notes <2017.1-2017.2>`.

* **28.02.2017** Version 2017.1 released (spline-box parametrization of an
  arbitrary field, conda-forge recipe, fixes for Python 3.6,
  etc.), see :ref:`release notes <2016.4-2017.1>`.

* **07.12.2016** Version 2016.4 released (support tensor product element meshes
  with one-level hanging nodes, improve homogenization support for large
  deformations, parallel calculation of homogenized coefficients and related
  sub-problems, evaluation of second derivatives of Lagrange basis functions,
  etc.), see :ref:`release notes <2016.3-2016.4>`.

* **30.09.2016** Version 2016.3 released (Python 3 support, testing with Travis
  CI, new classes for homogenized coefficients, using argparse instead of
  optparse, etc.), see :ref:`release notes <2016.2-2016.3>`.

* **12.05.2016** Version 2016.2 released (partial shell10x element
  implementation, parallel computation of homogenized coefficients, clean up of
  elastic terms, read support for msh mesh file format of gmsh, etc.), see
  :ref:`release notes <2016.1-2016.2>`.

* **24.02.2016** Version 2016.1 released (major simplification of finite
  element field code, automatic checking of shapes of term arguments, improved
  mesh parametrization code and documentation, support for fieldsplit
  preconditioners of PETSc, etc.), see :ref:`release notes <2015.4-2016.1>`.

* **01.12.2015** Version 2015.4 released (basic support for restart files,
  new type of linear combination boundary conditions, balloon inflation
  example, etc.), see :ref:`release notes <2015.3-2015.4>`.

* **23.09.2015** Version 2015.3 released (preliminary support for parallel
  computing, unified evaluation of basis functions (= isogeometric analysis
  fields can be evaluated in arbitrary points), (mostly) fixed finding of
  reference element coordinates of physical points, several new or improved
  examples, etc.), see :ref:`release notes <2015.2-2015.3>`.

* **29.05.2015** Version 2015.2 released (major code simplification (removed
  element groups), time stepping solvers updated for interactive use, improved
  finding of reference element coordinates of physical points, reorganized
  examples, reorganized installation on POSIX systems (sfepy-run script),
  etc.), see :ref:`release notes <2015.1-2015.2>`.

* **26.02.2015** Version 2015.1 released (support for multiple fields in
  isogeometric analysis, redesigned handling of solver parameters, new modal
  analysis example, etc.), see :ref:`release notes <2014.4-2015.1>`.

* **28.11.2014** Version 2014.4 released (preliminary support for 1D problems,
  data probes using pyVTK library, etc.), see :ref:`release notes
  <2014.3-2014.4>`.

* **25.09.2014** Version 2014.3 released (isogeometric analysis (IGA) speed-up
  by C implementation of NURBS basis evaluation, generalized linear combination
  boundary conditions that work between different fields/variables and support
  non-homogeneous periodic conditions, non-constant essential boundary
  conditions given by a function in IGA, reorganized and improved
  documentation, etc.), see :ref:`release notes <2014.2-2014.3>`.

* **23.05.2014** Version 2014.2 released (preliminary support for isogeometric
  analysis, improved postprocessing and visualization script for time-dependent
  problems with adaptive time steps, three new terms, etc.), see :ref:`release
  notes <2014.1-2014.2>`.

* **25.02.2014** Version 2014.1 released (sfepy.fem was split to separate
  FEM-specific and general modules, lower memory usage by creating active DOF
  connectivities directly from field connectivities, new handling of field and
  variable shapes, clean up: many obsolete modules were removed, all module
  names follow naming conventions, etc.), see :ref:`release notes
  <2013.4-2014.1>`.

* **22.11.2013** Version 2013.4 released (simplified quadrature definition,
  equation sequence solver, initial support for 'plate'
  integration/connectivity type, script for visualization of quadrature points
  and weights, etc.), see :ref:`release notes <2013.3-2013.4>`.

* **18.09.2013** Version 2013.3 released (implementation of Mesh topology data
  structures in C, implementation of regions based on C Mesh, MultiProblem
  solver for conjugate solution of subproblems, new advanced examples
  (vibro-acoustics, Stokes flow with slip conditions), etc.), see :ref:`release
  notes <2013.2-2013.3>`.

* **22.05.2013** Version 2013.2 released (automatic testing of term calls (many
  terms fixed w.r.t. corner cases), new elastic contact plane term + example,
  translated low level base functions from Cython to C for reusability,
  improved gallery http://docs.sfepy.org/gallery/gallery, etc.), see
  :ref:`release notes <2013.1-2013.2>`.

* **27.02.2013** Version 2013.1 released (unified use of stationary and
  evolutionary solvers, new implicit adaptive time stepping solver, elements of
  set and nodes of set region selectors, simplified setting of variables data,
  etc.), see :ref:`release notes <2012.4-2013.1>`.

* **21.11.2012** Version 2012.4 released (initial support for hierarchical
  basis on quadrilateral and brick elements, unified C/Cython structures for
  reference mappings, new linear combination boundary condition: edge
  direction, new examples showing some advanced features, etc.), see
  :ref:`release notes <2012.3-2012.4>`.

* **12.09.2012** Version 2012.3 released (several new terms, material
  parameters can be defined per region using region names, base function values
  can be defined per element, support for global options, etc.), see
  :ref:`release notes <2012.2-2012.3>`.

* **29.05.2012** Version 2012.2 released (reimplement acoustic band gaps code
  using the homogenization engine, high order quadrature rules, unify dot
  product and mass terms, lots of other term updates/fixes, update the PDE
  solver application, etc.), see :ref:`release notes <2012.1-2012.2>`.

* **27.02.2012** Version 2012.1 released (initial version of linearizer of
  higher order solutions, rewrite variable and evaluate cache history handling,
  lots of term updates/fixes/simplifications, move web front page to sphinx
  docs, etc.), see :ref:`release notes <2011.4-2012.1>`.

* **05.12.2011** Version 2011.4 released (cython used instead of swig to
  interface C code, many terms unified thanks to new optional material
  term argument type, updated Lagrangian formulation for large
  deformations, automatic generation of gallery of examples, etc.), see
  :ref:`release notes <2011.3-2011.4>`.

* **10.08.2011** Version 2011.3 released (major update of terms aiming at
  easier usage and definition while retaining original C functions,
  overriding problem description items on command line, improved
  developer guide, Primer tutorial - a step-by-step walk-through of the
  process to solve a simple mechanics problem, etc.), see
  :ref:`release notes <2011.2-2011.3>`.

* **31.05.2011** Version 2011.2 released (experimental implementation of
  terms aiming at easier usage and definition of new terms,
  Mooney-Rivlin membrane term, update build system to use exclusively
  setup.py, allow switching boundary conditions on/off depending on
  time, support for variable time step solvers, etc.), see
  :ref:`release notes <2011.1-2011.2>`.

* **24.03.2011** Version 2011.1 released (discontinuous approximations,
  user-defined material nonlinearities, improved surface approximations,
  speed-up mesh reading, extensive clean-up - less code, many bugfixes
  and many more updates), see
  :ref:`release notes <2010.4-2011.1>`.

* **06.12.2010** Version 2010.4 released (higher order elements,
  refactoring of geometries (reference mappings), transparent DOF vector
  synchronization with variables, interface variables defined on a
  surface region, many bugfixes and many more updates), see
  :ref:`release notes <2010.3-2010.4>`.

* **06.08.2010** Version 2010.3 released (significantly rewritten code for
  better interactive use, cleaner and simpler high level interface, new
  examples, tests, simplified but more powerful homogenization engine,
  many bugfixes), see :ref:`release notes <2010.2-2010.3>`.

* **10.05.2010** Version 2010.2 released (significantly updated
  documentation, new mesh readers, conversion formulas for elastic
  constants, basic tensor transformations, stress tensor conversion, new
  examples, tests, many new terms and bugfixes), see :ref:`release notes
  <2010.1-2010.2>`.

* **01.03.2010** Version 2010.1 released (new sphinx-based documentation,
  refactoring of base functions (polynomial spaces) and element geometry
  description, interpolation between different meshes, terms for
  describing perfusion and active fibres in the total Lagrangian
  formulation (applicable, for example, to active muscle tissue models)
  new tests, many new terms and bugfixes), see :ref:`release notes
  <2009.4-2010.1>`.

* **24.11.2009** Version 2009.4 released (greatly improved postprocessing
  and visualization capabilities, unified handling of user-defined
  functions, new tests, terms, many bugfixes), see :ref:`release notes
  <2009.3-2009.4>`.

* **21.07.2009** Version 2009.3 released (_Windows installation_, updated
  postproc.py - visualization, new tests, terms, solvers, bugfixes), see
  :ref:`release notes <2009.2-2009.3>`.

* **12.05.2009** Version 2009.2 released (new top level scripts
  (_isfepy_ - customized IPython shell, _postproc.py_ - mayavi2 based
  result viewer, _probe.py_), automatic html documentation generation
  via doxygen, new solvers, new mesh readers, extended syntax of
  equations for boundary traces of variables, short syntax for almost
  all input elements and other improvements), see :ref:`release notes
  <2009.1-2009.2>`.

* **02.03.2009** Version 2009.1 released (new solvers, new mesh readers,
  unified homogenization framework, dispersion analysis, phase velocity
  computation for phononic materials and other improvements), see
  :ref:`release notes <2008.4-2009.1>`.

* **04.12.2008** Version 2008.4 released (framework for running
  parametric studies, greatly improved support for time-dependent
  problems, live plotting using multiprocessing module, type of term
  arguments determined fully at run-time, new terms and other
  improvements), see :ref:`release notes <00.50.00-2008.4>`.

* **02.09.2008** Version 00.50.00 released (finite strain elasticity:
  total Lagrangian (TL) formulation, solving problems in complex
  numbers, generalized equations to allow linear combination of terms,
  run-time type of state term arguments, refactoring to follow Python
  coding style guidelines and other improvements), see :ref:`release
  notes <00.46.02-00.50.00>`.

* **01.07.2008** Version 00.46.02 released (alternative short syntax for
  specifying essential boundary conditions, variables and regions,
  manufactured solutions tests using !SymPy and other improvements),
  see :ref:`release notes <00.41.03-00.46.02>`.

* **26.03.2008** Version 00.41.03 released (works on 64 bits, support for
  various mesh formats, new solvers and other improvements), see
  :ref:`release notes <00.35.01-00.41.03>`.
