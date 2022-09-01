sfepy-run
=========

synopsis
--------

sfepy-run [-h] [--version] [-a {bvp,homogen,bvp-mM,evp,phonon}]
          [--debug] [--debug-mpi] [-c "key : value, ..."]
          [-O "key : value, ..."] [-d "key : value, ..."] [-o filename]
          [--format format] [--save-restart mode]
          [--load-restart filename] [--log file] [-q] [--save-ebc]
          [--save-ebc-nodes] [--save-regions]
          [--save-regions-as-groups] [--save-field-meshes]
          [--solve-not] [--phonon-band-gaps] [--phonon-dispersion]
          [--phonon-plot] [--phonon-phase-velocity] [--list what]
          [filename_in]

description
-----------

Solve partial differential equations given in a SfePy problem definition file.

Example problem definition files can be found in ``sfepy/examples/`` directory
of the SfePy top-level directory.

In the examples below it is supposed that sfepy is installed. When using the
in-place build, replace ``sfepy-run`` by ``python3 sfepy/scripts/simple.py``.

The supported application kinds (--app option) are:

- bvp - boundary value problem. Example::

    sfepy-run sfepy/examples/diffusion/poisson.py

- homogen - calculation of local microscopic problems (correctors) and
  homogenized coefficients. Example::

    sfepy-run sfepy/examples/homogenization/perfusion_micro.py

- bvp-mM - micro-macro boundary value problem. Solve a coupled two-scale
  problem in parallel using MPI. One computational node is solving a
  macroscopic equation while the others are solving local microscopic problems
  and homogenized coefficients. The --app option is required in this case.
  Example::

    mpiexec -n 4 sfepy-run --app=bvp-mM --debug-mpi sfepy/examples/homogenization/nonlinear_hyperelastic_mM.py

- evp - eigenvalue problem. Example::

    sfepy-run sfepy/examples/quantum/well.py

- phonon - phononic band gaps. Example::

    sfepy-run sfepy/examples/phononic/band_gaps.py --phonon-plot

Both normal and parametric study runs are supported. A parametric study allows
repeated runs for varying some of the simulation parameters - see
``sfepy/examples/diffusion/poisson_parametric_study.py`` file.

positional arguments
--------------------

| filename_in         SfePy problem description file

optional arguments
------------------

-h, --help            show this help message and exit
--version             show program's version number and exit
-a {bvp,homogen,bvp-mM,evp,phonon}, --app {bvp,homogen,bvp-mM,evp,phonon}
                      override application kind, normally determined
                      automatically. The supported kinds are: bvp (boundary
                      value problem), homogen (correctors, homogenized
                      coefficients), bvp-mM (micro-macro boundary value
                      problem, homogenized coefficients computed in parallel
                      using MPI), evp (eigenvalue problem), phonon (phononic
                      band gaps)
--debug               automatically start debugger when an exception is
                      raised
--debug-mpi           log MPI communication (mM mode only)
-c "key : value, ...", --conf "key : value, ..."
                      override problem description file items, written as
                      python dictionary without surrounding braces
-O "key : value, ...", --options "key : value, ..."
                      override options item of problem description, written
                      as python dictionary without surrounding braces
-d "key : value, ...", --define "key : value, ..."
                      pass given arguments written as python dictionary
                      without surrounding braces to define() function of
                      problem description file
-o filename           basename of output file(s) [default: <basename of
                      input file>]
--format format       output file format, one of: {vtk, h5} [default: vtk]
--save-restart mode   if given, save restart files according to the given
                      mode.
--load-restart filename
                      if given, load the given restart file
--log file            log all messages to specified file (existing file will
                      be overwritten!)
-q, --quiet           do not print any messages to screen
--save-ebc            save a zero solution with applied EBCs (Dirichlet
                      boundary conditions)
--save-ebc-nodes      save a zero solution with added non-zeros in EBC
                      (Dirichlet boundary conditions) nodes - scalar
                      variables are shown using colors, vector variables
                      using arrows with non-zero components corresponding to
                      constrained components
--save-regions        save problem regions as meshes
--save-regions-as-groups
                      save problem regions in a single mesh but mark them by
                      using different element/node group numbers
--save-field-meshes   save meshes of problem fields (with extra DOF nodes)
--solve-not           do not solve (use in connection with --save-*)
--phonon-band-gaps    detect frequency band gaps
--phonon-dispersion   analyze dispersion properties (low frequency domain)
--phonon-plot         plot frequency band gaps, assumes -b
--phonon-phase-velocity
                      compute phase velocity (frequency-independent mass
                      only)
--list what           list data, what can be one of: {terms, solvers}
