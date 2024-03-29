# created: 20.07.2007 (-1)

.. _2023.4-2024.1:

from 2023.4 to 2024.1
=====================

- merge pull request #1056 from rc/l2-constant-field

  - fix docstring of FEField.get_dofs_in_region()
  - fix FEField.get_econn() error message
  - replace FEField._setup_geometry() by new _find_geometry() - update
    .__init__()
  - new sfepy/discrete/fem/fields_l2.py

    - new L2ConstantVolumeField, .__init__(), ._create_interpolant(),
      .setup_extra_data(), .get_coor(), .get_econn(), .get_data_shape(),
      .get_dofs_in_region(), .create_mapping(), .create_output()
    - new L2ConstantSurfaceField

  - update Field.from_conf()
  - update ExpressionArg.get_bf()
  - update L2ConstantVolumeField.__init__(), .get_econn()
  - update balloon.py example to optionally use L2 constant field - update
    plot_radius(), define()
  - fix ExpressionArg.get_dofs() for argument time derivatives
  - new L2ConstantVolumeField.get_base()
  - fix SurfaceFluxOperatorTerm.get_fargs()
  - fix L2ConstantVolumeField.create_mapping() for volume field + surface
    integral
  - fix ExpressionArg.get_bf() for single cell meshes
  - new has_time_derivatives solver option of ElastodynamicsBaseTS - update
    ._common_parameters, .__call__()
  - simplify ProblemConf.get_function()
  - update piezo_elastodynamic.py example to use L2 constant field

    - new plot_voltage()
    - update post_process(), define()

  - gen_gallery.py: update custom view for piezo_elastodynamic.py example
  - fix FEField.get_qp() for boundary QP requests
  - re-fix SurfaceFluxOperatorTerm.get_fargs() for different state and virtual
    fields
  - docs: sync module index of developer guide with current sources

- merge pull request #1057 from rc/update-fields-docs

  - new tools/gen_field_table.py - new typeset_field_table(), typeset(),
    gen_field_table(), setup(), main()
  - docs: document all field configuration options in users guide
  - add gen_field_table.py to sphinx extensions in doc/conf.py
  - docs: add field table to users guide
  - clean up whitespace in sfepy/discrete/dg/fields.py
  - add/update docstrings of fields for field table
  - docs: explain field table columns in users guide
  - gen_field_table.py: update typeset_field_table()

- merge pull request #1059 from fix-read-probe-results

  - fix read_results()

- merge pull request #1060 from rc/update-mesh-info

  - update create_bnf() to allow '-' in vertex set names
  - show_mesh_info.py: print vertex and cell groups in show_mesh_info()
  - show_mesh_info.py: print data stats in show_mesh_info()

- merge pull request #1061 from rc/fix-hyperelasticity-in-subdomains

  - fix HyperElasticFamilyData.__call__() for subdomain integrals

- merge pull request #1062 from rc/docs-update-support

  - docs: update support

- merge pull request #1063 from rc/passive-exp-tl-fibres-term

  - new FibresExponentialTLTerm (dw_tl_fib_e) - new .get_fargs(),
    .stress_function(), .tan_mod_function(), .get_eval_shape()

- merge pull request #1064 from rc/update-mechanics-tensors

  - new StressTransform.get_1pk_from_2pk()
  - update docstring of StressTransform.get_cauchy_from_2pk()
  - use numpy.linalg.det() directly in  StressTransform.__init__()

- merge pull request #1065 from rc/hyperelasticity-on-facet-regions

  - fix dw_he_rtm() for facet_extra integration
  - new Term.arg_geometry_types attribute to allow per-variable integration

    - update Term.setup_geometry_types() - geometry type can be overridden per
      variable/mode for each integration

  - set .arg_geometry_types in TL/UL terms with scalar variables - update
    BulkPressureTLTerm, VolumeTLTerm, BulkPressureULTerm, VolumeULTerm,
    CompressibilityULTerm
  - add facet_extra integration to HyperElasticBase
  - add facet_extra integration to DeformationGradientTerm, update .get_fargs()

- merge pull request #1067 from vlukes/fix_elastic_terms

  - fix order of arguments

- merge pull request #1067 from vlukes/update_esdlinearelasticterm

  - update ESDLinearElasticTerm for non-symmetric material tensor

- merge pull request #1068 from rc/fix-deprecation-warnings

  - fix deprecation warnings in ScalarDotGradIScalarTerm, GenYeohTLTerm
  - fix deprecation warning in ETermBase.get_eval_shape()
  - use recommended jax.config import
  - fix deprecation warnings in test_compose_sparse()

- merge pull request #1069 from rc/report-solver-status

  - new report_status parameter of Solver
  - report status in Newton.__call__()
  - report status in ScipyBroyden.__call__()
  - report status in PETScNonlinearSolver.__call__()
  - new IndexedStruct.setdefault() for dict compatibility
  - new standard_nls_call() decorator - apply it to .__call__() of Newton,
    ScipyBroyden, PETScNonlinearSolver
  - patch poststep_fun() in standard_ts_call() to collect nls time stats
  - report nls time stats in Problem.solve()
  - new Timer.add()
  - new Timers - new .__init__(), .create(), .start(), .stop(), .get_dts(),
    .get_totals()
  - tweak time stats reports of ScipyBroyden, PETScNonlinearSolver
  - fix time stats totals report in Newton.__call__() by using Timers
  - fix AdaptiveTimeSteppingSolver.solve_step() to use nls.status
  - update define() of balloon.py example to report nls status
  - replace patched poststep_fun() by _TimingNLS in standard_ts_call()

    - nls time stats collected after every nls call
    - fixes total nls time stats for adaptive time-stepping

  - special-case balloon.py example in test_examples()
  - test_install.py: update main() for changed output

- merge pull request #1071 from rc/fix-multilinear-terms-cell-basis

  - fix cell dependent basis detection in ExpressionArg.get_bf()

- merge pull request #1073 from vlukes/fix-mixed_mesh-example

  - move meshes from '2d' to '3d', fix format

- merge pull request #1050 from vlukes/wedge_elements

  - new wedge finite element
  - new PolySpace for wedge element
  - remove forgotten comment
  - remove trailing whitespace
  - fix wedge quadrature
  - fix quadrature points for 3_6 element
  - clean-up in LagrangeWedgePolySpace._eval_base()
  - update facet mapping for wedge elements (tri and quad facets)
  - fix FEMapping.get_physical_qps() for wedge elements
  - fix FEField.get_qp()
  - new linear elastic example: test wedge elements
  - add wedge elements example to tests
  - update wedge quadrature rules
  - fix wedge surface integration for higher quadratures
  - update example: set integration order to '2'
  - remove 'six' dependency
  - add forgotten comma
  - update example docstring
  - fix quadrature rule and quadrature tests
  - example: update docstring and remove obsolete import
  - fix docstring

- merge pull request #1070 from vlukes/report-step-stats

  - report status of nonlinear steps
  - fix nls stats in ossen solver
  - print step info only if `report_status` is True
  - log time_stats and step_stats into a file
  - move log_file handle to status structure
  - fix deprecation warning: np.sum(generator)
  - use fixed formatting of step stats in Problem.solve()
  - move report_status and log_status to global options
  - update hyperela. example: switch on nls report and log
  - fix previous wrong deletion
  - change log file suffix: .txt -> .csv
  - move 'verbose' argument to the last position

- merge pull request #1075 from rc/skip-lin-precision-check

  - skip linear system precision check in Newton.__call__() if lin_red is None
  - update solver settings in balloon.py example
  - update test_examples()

- merge pull request #1076 from rc/update-test-install

  - test_install.py: update main() for changed output
    - remove checks of examples already tested by test_declarative_examples.py

- merge pull request #1074 from rc/visualize-probes

  - new parse_vector(), parse_scalar()
  - parse probe parameters in read_header() - new .parse_report() of probes -
    update PointsProbe, LineProbe, RayProbe, CircleProbe
  - new ensure3d(), read_probes_as_annotations()
  - resview.py: new --probes option, update pv_plot(), main() to show probes
  - support probe labels in read_probes_as_annotations(), pv_plot()
  - resview.py: new --no-probe-labels option in main()

- merge pull request #1078 from rc/resview-show-time

  - update read_mesh() to return time
  - update pv_plot() to show step, time and filename (if multiple) as text
  - cache actual step in read_mesh(), read .h5 steps on demand, update
    pv_plot()
  - resview.py: new --no-step-time option in main(), update pv_plot()

- merge pull request #1082 from vlukes/update-evp-app

  - update EVPSolverApp.save_results(): force set_state() to apply LCBCs

- merge pull request #1083 from vlukes/matcoefs-orthotropic

  - new function for stiffness of an orthotropic material
  - fix spelling
  - rename mat. function to 'stiffness_from_yps_ortho3'

- merge pull request #1084 from rc/fix-project-link

  - web: fix project link

- merge pull request #1085 from vlukes/multi-node-combination-lcbc

  - new LCBC operator: NodalCombinationXOperator
  - fix error message
  - place NodalCombinationXOperator class in front of LCBCOperators
  - fix make_global_operator() and finalize() for shared dofs
  - remove six
  - remove __future__ import
  - update: constraints can be generated by dof_map_fun function
  - rename oprerator and add docstring
  - fix operator matrix values
  - new example for 'multi_node_combination' lcbc
  - add multi-node-combination BC to tests
  - fix docstring of MultiNodeLCOperator()
  - fix missing "=" in the math expression of MultiNodeLCOperator()

- merge pull request #1089 from vlukes/fix-dw_vm_dot_s-term

  - fix VectorDotScalarTerm term for fmode=2

- merge pull request #1092 from rc/update-gallery-custom-views

  - gen_gallery.py: fix missing pv_plot() option in generate_images()
  - gen_gallery.py: add custom view for wedge_mesh.py example
  - gen_gallery.py: add custom view for multi_node_lcbcs.py example

.. _2023.3-2023.4:

from 2023.3 to 2023.4
=====================

- merge pull request #1029 from rc/fix-modal-analysis-declarative-docstring

  - fix docstring of modal_analysis_declarative.py example

- merge pull request #1030 from rc/clean-up-poisson-nonlinear-material

  - remove superfluous declarations from poisson_nonlinear_material.py example

- merge pull request #1031 from rc/gen-gallery-pattern-option

  - gen_gallery.py: new --pattern option - update generate_images(),
    generate_rst_files(), main()

- merge pull request #1027 from vlukes/new_structural_elements

  - simplify ConnInfo, remove dof_conn_type structure
  - merge FEField.surface_data and FEField.point_data into FEField.extra_data
  - update arguments of FEField.setup_extra_data()
  - allow 1_2 cells in 2D and 3D
  - new LinearSpringTerm, LinearTrussTerm and LinearTrussInternalForceTerm
  - new region extra opt. 'finalize': to avoid region.finalize() call
  - add missing tdim argument in domain.get_conn() calls
  - new 2D example: truss bridge
  - new 3D example with solid and structural elements
  - new meshes: 2D and 3D bridges
  - update term tests: do not test truss and spring terms
  - add custom views to the gallery generator

- merge pull request #1032 from rc/sem

  - new _get_table(), update PolySpace.any_from_args()
  - new register_poly_space()
  - new get_lgl_nodes(), eval_lagrange1d_basis(), SEMTensorProductPolySpace
  - new H1SEMVolumeField, H1SEMSurfaceField
  - test SEM basis in test_poly_spaces.py - update _gen_common_data(),
    test_partition_of_unity(), test_continuity()
  - new sfepy/examples/miscellaneous/refine_evp.py

- merge pull request #1034 from rc/resview-fix-plot-positions

  - resview.py: fix different positions of single slot plots in pv_plot()

- merge pull request #1033 from vlukes/multi_tdim_problem

  - update field.get_econn() argument and function calls
  - remove Term.dof_conn_info, update ConnInfo and Term.geometry_types
  - fix 'volume' integration to 'cell'
  - simplify ConnInfo.get_region_name()
  - remove FEField.get_connectivity()
  - update ConnInfo struct: remove unused arttibutes v_tg, ps_tg
  - update ConnInfo: dof_conn_type replaced by dof_conn_types list
  - update regions: create a new cmesh with a lower dim. form the region
  - rename region extra option: 'tdim' -> 'mesh_dim'
  - update vector fields: num. of components according to region.cmesh.tdim
  - update mappings: new coor. transformation to lower dimensions
  - update vibro-acoustic example: 'ls.cm_pb' solver is no more needed
  - rename variable in Equations.get_graph_conns(): cpname -> cname
  - explicit precedense of operators in the conditions
  - fix docstring of FEField.get_econn()
  - update example: rearrange the code into define(), 'ls.scipy_direct' ->
    'ls.auto_direct'
  - fix IGField.setup_extra_data(): wrong info.dof_conn_types handling
  - fix Variable.evaluate(): swop region/mirror region
  - update FESurface class: make meconn, mleconn dicts for various mirror
    regions
  - update Field.get_econn(): get connectivity for FEPhantomSurface
  - fix mirror region connectivity

- merge pull request #1037 from vlukes/fix_broken_examples, closes #1035

  - fix nonlinear homogenization example: follow changes in PR #1033
  - fix for issue #1035: new attribute field_dim in Region class
  - nonlinear_homogenization example: fix the previous fix

- merge pull request #1038 from vlukes/fix_field_base_create_mesh

  - fix FEField.create_mesh()

- merge pull request #1039 from vlukes/fix_create_basis_context

  - fix H1NodalMixin.create_basis_context()

- merge pull request #1036 from rc/gallery-add-interactive-examples - closes
  #371, #905, #942

  - gen_gallery.py: support script/interactive examples

    - commands can be defined in custom views dict
    - new _get_image_names(), _apply_commands()
    - update _get_fig_filenames(), generate_images(), generate_rst_files()

  - gen_gallery.py: report failed interactive examples in generate_images()

    - update _apply_commands() to use subprocess.run(), raise exception on
      error

  - gen_gallery.py: add custom view commands for interactive examples - update
    omit_images
  - fix LogPlotter.apply_commands() for matplotlib >= 3.6.0
  - fix probe_results() in time_poisson_interactive.py example
  - gen_yeoh_tl_up_interactive.py: new --no-show option, update docstring
  - hyperelastic_tl_up_interactive.py: new --no-show option, update docstring
  - add basic docstring to linear_elastic_interactive.py example
  - clean up post_process() in linear_elastic_probes.py example
  - update docstring of shell10x_cantilever.py example
  - shell10x_cantilever_interactive.py: new --no-show option, update docstring
  - live_plot.py: new --output-dir, --plot-log options, add docstring
  - use clip transform in band_gaps_rigid.py example for nicer plots
  - gen_gallery.py: fix acoustics/vibro_acoustic3d.py custom view
  - gen_gallery.py: fix generate_gallery() for no images
  - gen_gallery.py: add custom view commands for two homogenization examples

    - support linear_homogenization.py, perfusion_micro.py
    - update omit_images

  - update docstring of linear_homogenization.py example
  - update docstring of perfusion_micro.py example
  - add basic docstring to linear_elastic_mM.py example
  - gen_gallery.py: new --no-thumbnails option, update main()
  - gen_gallery.py: fix generate_images() for any separator
  - fix docs build without petsc4py, mpi4py
  - gen_gallery.py: allow plot failures in generate_images()
  - gen_gallery.py: allow import failures in generate_rst_files()
  - gen_gallery.py: update custom views of 1D results
  - gen_gallery.py: clean up custom views
  - resview.py: ensure nonzero plot shifts in pv_plot() for all axes
  - gen_gallery.py: fix result name for its2D_4.py example custom view
  - gen_gallery.py: fix commands of dispersion_analysis.py example custom view
  - remove debug code in print_camera_position()

- merge pull request #1042 from vlukes/fix_integer_division

  - fix ccontres.pyx: replace floating-point divisin by integer division

- merge pull request #1044 from rc/print-ccore-exceptions

  - print exception message in errclear()

- merge pull request #1045 from vlukes/multi_cell_mesh

  - fix element orientation function
  - update Domain.create_region(): borrow vertices from another region
  - update FEDomain.get_conn()
  - remove DofInfo.ptr, ptr[-1] replaced by DofInfo.n_dof_total
  - update eq. mapping to be able to share dofs between variables
  - fix surface connectivity
  - fix FEField.create_mapping() for surface regions
  - new example: a beam consisting of hexa and tetra elements
  - fix FEField.create_mapping(): connectivity for 'surface_extra' integration
  - use fem.utils.prepare_remap()
  - fix example docstring

- merge pull request #1047 from rc/fix-probes-write-results

  - fix write_results()
  - fix Field.evaluate_at() docstring

- merge pull request #1049 from flothesof/patch-1

  - better eigenvalue solver for modal_analysis_declarative.py
  - Update docstring numerical values

- merge pull request #1051 from rc/fix-docs

  - fix comment in define() of refine_evp.py example
  - docs: fix command in primer
  - gen_gallery.py: prefix custom images to prevent name clashes

    - new _make_fig_name()
    - update _get_fig_filenames(), _apply_commands(), generate_images()

  - gen_gallery.py: new run_resview_plot(), update generate_images()

- merge pull request #1052 from rc/jax-he-tl-terms-proof-of-concept

  - new NeoHookeanTLADTerm (dw_tl_he_neohook_ad)
  - new OgdenTLADTerm (dw_tl_he_ogden_ad) (WIP)
  - do not test dw_tl_he_ogden_ad in test_term_call_modes() - singular matrix
    for zero displacements

.. _2023.2-2023.3:

from 2023.2 to 2023.3
=====================

- merge pull request #1006 from yosefm/deprecation-warnings

  - Fix some deprecations on NumPy 1.25 and newer setuptools
  - Explicit dependency, just for safety

- merge pull request #1009 from yosefm/cython_version

  - Restrict Cython

- merge pull request #1008 from bubulk/sfepy-docker

  - Updated sfepy-docker doc.

- merge pull request from #1012 from peppe988/my-feature

  - example for the new non linear terms
  - Update Non_linear_general_poisson_equation.py
  - Update and rename Non_linear_general_poisson_equation.py to
    poisson_nonlinear_material.py

- merge pull request #1013 from rc/test-poisson-nonlinear-material

  - test poisson_nonlinear_material.py in test_examples()
  - tweak formatting of poisson_nonlinear_material.py example

- merge pull request #1015 from rc/fix-solver-init

  - do not modify conf argument in Solver.__init__()

- merge pull request #1016 from rc/fix-active-tl-fibres-qp-mode

  - fix FibresActiveTLTerm.get_eval_shape() for qp mode

- merge pull request #1017 from vlukes/piezo_flow_example

  - add a new example application

- merge pull request #1020 from vlukes/new_v_dot_grad_s_sensitivity_term

  - fix docstring
  - new sensitivity term: de_sd_v_dot_grad_s
  - fix format

- merge pull request #1022 from ostueker/ostueker-patch-1

  - include spaces when joining cflags, configure_cppflags, configure_cflags

- merge pull request #1024 from vlukes/mumps_multiple_rhs, closes #1023

  - implement multiple RHS for mumps solver
  - multiple RHS for SchurMumps() solver
  - always make a copy of rhs
  - fix for using ls_mumps standalone

- merge pull request #1025 from rc/differentiable-terms

  - new Term.diff_info attribute, update .evaluate()
  - check 'eval' mode availability in Term.evaluate()
  - support sensitivity dw_mode in Equations.evaluate(), Equation.evaluate() -
    new diff_vars argument
  - optionally clear constant_matrices in
    ElastodynamicsBaseTS.clear_lin_solver()

    - update .__call__()
    - update BatheTS.clear_lin_solver()

  - clear linear solver in Problem.solve()
  - allow passing ElastodynamicsBaseTS instance to init_fun(), poststep_fun() -
    update .__call__()
  - store pack(), unpack() in ElastodynamicsBaseTS instance in
    .get_initial_vec()
  - fix typo in comment
  - new sfepy/terms/terms_jax.py - JAX-based terms with automatic
    differentiation

    - new LinearElasticLADTerm (dw_lin_elastic_l_ad)
    - new LinearElasticYPADTerm (dw_lin_elastic_yp_ad)
    - new MassADTerm (dw_mass_ad)
    - new get_strain(), get_stress(), ceval_elasticity_l(),
      ceval_elasticity_yp(),   ceval_mass()

  - fix process_terms(), typeset_term_syntax() in tools/gen_term_table.py
  - new sfepy/examples/linear_elasticity/elastodynamic_identification.py
    example

    - new define(), apply_sensor(), update_pars(), eval_fun(), eval_jac_fd(),
    eval_jac(), parse_args(), main()
    - new NewmarkSATS, .create_nlst()

  - docs: sync module index of developer guide with current sources
  - clear linear solver explicitly in eval_fun(), revert Problem.solve() -
    update elastodynamic_identification.py example
  - test_install.py: test elastodynamic_identification.py example
  - use argument type names in *ADTerm.diff_info - update LinearElasticLADTerm,
    LinearElasticYPADTerm, MassADTerm
  - new part argument in Term.get_material_names()
  - update Equation.evaluate(), Term.evaluate() for actual Term.diff_info
  - elastodynamic_identification.py: simplify eval_jac_fd()
  - elastodynamic_identification.py: clear linear solver in eval_jac()
  - elastodynamic_identification.py: print() -> output()
  - elastodynamic_identification.py: implement proportional damping in define()
  - elastodynamic_identification.py: new --alpha, --beta options, identify
    damping - update efun_grad_par(), update_pars(), eval_fun(), eval_jac(),
    parse_args(), main()

- merge pull request #1026 from rc/fix-primme-n-eigs-none

  - fix PrimmeEigenvalueSolver.__call__() for no n_eigs provided

.. _2023.1-2023.2:

from 2023.1 to 2023.2
=====================

- merge pull request #977 from rc/docs-update-support

  - docs: update support

- merge pull request #978 from vlukes/update_users_guide

  - update User's guide: new 'match_dofs' lcbcs
  - fix punctation

- merge pull request #979 from rc/fix-piezo-ed-example-doctring

  - fix docstring of piezo_elastodynamic.py example

- merge pull request #980 from rc/fix-presolve-in-homogenization

  - do not presolve in MiniAppBase.init_solvers()
  - replace ls.scipy_direct by ls.auto_direct with presolve in hom. examples
  - replace ls.auto_iterative by ls.auto_direct with presolve in example
  - use direct solver for rs correctors in example
  - set is_linear to True in hom. examples where appropriate
  - set empty solvers dict to default conf in Problem.__init__()
  - set automatically use_presolve to True in MiniAppBase.init_solvers()

- merge pull request #960 from yosefm/numpy-distutils-byebye, closes #944

  - initial make files
  - File removed in master
  - Add another depenency of terms to CMakeLists
  - Don't ignore cmake
  - Build file depending on another directory
  - Change installation target dir
  - Make install testable
  - Change order of commands to what's apparently expected. See
    https://scikit-build.readthedocsa.io/en/latest/cmake-modules/PythonExtensions.html
    'Amend the configuration of the library target (created using add_library)
    ... '
  - sfepy_common should be static Otherwise it's not found on runtime. Also, no
    need to install it then.
  - More builds Should get us done with sfepy.discrete, but tests fail to
    import (but not `from ... import`)
  - Get tems.extmods fully built So that sfepy.discrete can import
  - Install tests as data files so that tests actually run. some of them fail
    on not finding their meshes, but that's for next commit.
  - Add meshes directory to install Maybe would be preferable to follow the
    setuptools current advice and move it into the package, but this branch is
    already shaking the boat pretty hard, so let's leave it for a later PR.
  - Fix some missing links
  - Add mechanics/extmods first time all tests pass
  - Move common definitions to main CMake file
  - clean up a bit
  - Clean up early experiments
  - Compose compile flags as in numpy distutils
  - distutils only adds cflags for unices
  - Handle debug flags in site_cfg
  - Responsibility for link flags moves to CMake fully Until anyone demands it
    back. It would be better to do this through CMakeLists anyway.
  - Install examples too Leave wider pyproject.toml reform for later :)
  - Feature completion
  - Use setup_test.py as the new setup.py
  - Add skbuild to pyproject.toml
  - Remove 3.9 annotations
  - Steal distutils hack to install data with package.
  - Editable install not well-supported by skbuild
  - Fix: doesn't write the version file.
  - Run coverage using installed package
  - Fix relative paths
  - Windows peculiariuty?
  - Removals
  - One more removal
  - Observe line length in an skbuild-acceptable way
  - disable C++ in CMake As it will not be found in the conda recipe and is not
    actually needed, since the conda recipe always specified onlt a C compiler.
  - Languages also inb setup.py
  - Use Ninja generator explicitly
  - Add new build dependencies also to docs
  - Keep support for space-separated flag strings
  - Manually include compiled files in manifest
  - replace global-include by recursive-include, exclude _skbuild
  - Catch manifest up with removed file in master

- merge pull request #987 from rc/fix-missing-mumps_ls-attribute

  - fix no mumps_ls attribute in MUMPSSolver.clear() (via .__del__())

- merge pull request #988 from rc/fix-resview-for-pyvista-0.39.0

  - fix make_grid_from_mesh() for pyvista>=0.39.0

- merge pull request #989 from rc/remove-save-field-meshes

  - remove unused Problem.save_field_meshes() - Field.write_mesh() was removed
    already in 2011.1
  - remove save_field_meshes from default options
  - simple.py: remove --save-field-meshes option

- merge pull request #990 from rc/document-file-format-option

  - docs: document file_format option in users guide, clean up

- merge pull request #991 from rc/update-lobpcg-init

  - initialize lobpcg() in LOBPCGEigenvalueSolver.__call__() by random data -
    seed fixed for repeatability

- merge pull request #992 from rc/resview-factor-by-warp

  - compute factor in pv_plot() from warp field if given

- merge pull request #994 from rc/resview-window-title

  - resview.py: set window title to filename(s), new make_title() - update
    pv_plot(), main()

- merge pull request #995 from vlukes/primme_eigensolver

  - new eigenvalue solver: PRIMME
  - set default value of 'tol' to 0
  - conf.which() -> conf.which.upper()
  - return eigenvectors?
  - add eig.primme to test_eigenvalue_solvers()

- merge pull request #982 from vlukes/modal_analysis_declarative

  - new example: simple modal analysis in declarative description
  - update docstring: running the simulation
  - add the example to tests
  - increase solver tolerance
  - new mesh (wheelset.vtk) and eigenvalue solver (primme)
  - new mesh file
  - update tests: modal_analysis_declarative.py uses the optional PRIMME solver
  - update mesh
  - update docstring of modal_analysis_declarative.py
  - gen_gallery.py: add custom view for modal_analysis_declarative.py example

- merge pull request #993 from rc/remove-ordering-variables-attributes

  - remove Variable.reset(), ._count, ._orders, ._all_var_names

    - update .__init__(), ._set_kind()
    - update Variables.from_conf()

  - update Variables.setup_ordering()
  - do not call .setup_ordering() in Variables.from_conf(), .__init__() - it is
    called in .setup_ordering()
  - update transform_equations_ed() for new Variables.setup_ordering()
  - fix typo in Variables.setup_ordering()

- merge pull request #996 from vlukes/update_modal_analysis_example update
  example: add comparison to ansys results

  - update example: add comparison to Ansys results

- merge pull request #997 from rc/remove-distutils-remains, closes #984

  - update site_cfg_template.py
  - replace distutils in tools/build_helpers.py, remove
    generate_a_pyrex_source(). update Clean
  - replace LooseVersion by parse_version()
  - replace find_executable() with which() in test_gen_mesh_from_geom()
  - fix warning and wrong type in CMesh.create_new()
  - fix error messages in GmshIO
  - fix matplotlib warning in time_heat_equation_multi_material.py exemple,
    update post_process_hook()
  - update warning filters in pytest.ini

- merge pull request #1001 from peppe988/my-feature

  - new non linear term
  - Update terms_diffusion.py
  - Update terms_volume.py
  - updated
  - solved a few whitespace errors
  - allow callable in place of shape in make_term_args()
  - update arg_shapes of NonlinearDiffusionTerm, NonlinearVolumeForceTerm

    - add missing import
    - test_term_call_modes.py::test_term_call_modes now passes

  - fixed latex mode and empty lines
  - fix indent

.. _2022.4-2023.1:

from 2022.4 to 2023.1
=====================

- merge pull request #917 from rc/fix-typos

  - fix typos

- merge pull request #918 from rc/convert-mesh-force-3d

  - convert_mesh.py: new --3d option in main()

- merge pull request #919 from rc/resview-camera-position-handling

  - resview.py: unify camera position handling in animation and normal modes

    - set default value of --view option to None
    - --camera-position and --view have precedence over --2d-view
    - new _get_cpos(), update main()

- merge pull request #920 from vlukes/on_screen_animations

  - fix the off-screen animation issue mentioned in #919

- merge pull request #921 from rc/mass-term-amm-rmm

  - new sfepy/terms/terms_mass.py: new MassTerm (de_mass)
  - update MassTerm docstring
  - docs: sync module index of developer guide with current sources
  - new ReciprocalMassMatrixSolver (ls.rmm) (WIP) - new .__init__(),
    .init_rmm(), .__call__()
  - update MassTerm.get_function() for residual mode
  - rename ReciprocalMassMatrixSolver -> RMMSolver, update docstring
  - obey active_only in RMMSolver.init_rmm(), update .__call__()
  - use dedicated NoNLS for RMM in CentralDifferenceTS - update
    ._create_nlst_a(), .create_nlst()
  - explain a0 in RMMSolver.__call__()
  - new output_dir define() argument in elastodynamic.py example
  - update elastodynamic.py example to show RMMSolver use, support mass
    lumping - new mass_beta, mass_lumping, fast_rmm define() arguments
  - set tscedl as default in define() of elastodynamic.py example
  - new dt, edt_safety define() arguments, update docstring of elastodynamic.py
  - update docstring of MassTerm
  - update docstring of RMMSolver
  - update docstrings of VelocityVerletTS, CentralDifferenceTS
  - optimize beta == 0 branch in fun() in MassTerm.get_function()
  - new test_rmm_solver(), update define() in test_ed_solvers.py
  - improve code for ths storing in test_ed_solvers()
  - update make_term_args(), _test_single_term() for MassTerm - fixes
    test_term_call_modes()

- merge pull request #922 from rc/variable-has-ebc

  - new FieldVariable.has_ebc(), update Variables.has_ebc(), add verbose
    argument

- merge pull request #923 from rc/resview-animation-cpos

  - resview.py: fix and simplify setting camera position in animation mode

- merge pull request #924 from vlukes/match_dofs

  - new MatchDOFsOperator: tie DOFs of two fields

- merge pull request #928 from rc/sfepy-view-png-anim

  - resview.py: support saving animation as PNG sequence in main()

- merge pull request #925 from vlukes/multi_topo

  - new copy_coors argument in CMesh.from_data()
  - read mesh contaning cells of various topological dimensions
  - cmesh linked to regions according to the topological dimension of cells
  - update FEDomain.fix_element_orientation(): check "real" cells only
  - cleanup: remove six
  - link cmeshs to fields by topological dimensions
  - update mappings
  - update DG field
  - update IGA domain
  - update recover_micro_hook()
  - fix domain.region_leaf(): 'KW_COG' -> 'E_COG'
  - fix assert_equal(): add nm.number to assert_base_types
  - move back to cmesh (cmesh_highest -> cmesh)
  - fix _check_region(): region.cmesh can not be None
  - update docstrings

- merge pull request #929 from vlukes/new_example

  - add link to 'reinforced shell beam' example

- merge pull request #930 from vlukes/fix_multi_topo

  - new tdim argument of function Region.light_copy()

- merge pull request #935 from vlukes/update_mirror_regions

  - fix region.copy(): call light_copy with tdim=self.tdim argument
  - fig setup_mirror_region(): return mirror_name if ret_name is True
  - update mirror region: mirror with tdim = dim - 1
  - replace 'is_trace' by 'trace_region'
  - fields_base.py clean-up

- merge pull request #934 from yosefm/numpy-dep-py39

  - Fix build dependency for Python 3.9

- merge pull request #936 from rc/multilinear-terms-surface-extra

  - support surface_extra integration in ExpressionBuilder - update
    .add_virtual_arg(), .add_state_arg(), ExpressionArg.get_bf()
  - new SurfaceFluxOperatorTerm (de_surface_flux)

- merge pull request #939 from yosefm/fix-pyvista-deprecation

  - Pyvista deprecation: must provide inplace parameter

- merge pull request #910 from flothesof/add_github_ci

  - add CI script
  - reformat ci file
  - try other branch name
  - disable build on push to master and add documentation
  - add link to CI file in description
  - add changes proposed by @rc

- merge pull request #940 from yosefm/fix-elastodynamic-solver-return

  - Return last vec from elastodynamic solvers

- merge pull request #943 from yosefm/remove-setuptools-limitation

  - Remove setuptools limitation

- merge pull request #946 from rc/update-variable-ts-print-info

  - use supplied n_step for print info in VariableTimeStepper.set_from_data()

- merge pull request #945 from yosefm/explain-missing-odes

  - Explain the relations between variables in example Elastodynamic problem
    defines 3 "independent" variables that are dependent by solver magic.
  - CR tweaks

- merge pull request #937 from vlukes/mappings_refactoring

  - unification of VolumeMapping and SurfaceMapping
  - element diameters calculated in domain, not in CMapping
  - remove refmaps.c, replaced by DMapping class and eval_mapping_data_in_qp()
  - employ DMapping in FE and IGA mappings
  - update FE fields for new mapping
  - use CMapping instead of DMapping when callig C functions
  - translete DMapping to CMapping in fargs for C functions
  - CMapping.shape removed, use DMapping.n_el instead
  - add docstring to eval_mapping_data_in_qp()
  - call einsum with `optimize=True`
  - rename DMapping -> PyCMapping
  - rena coorIn to coor in eval_mapping_data_in_qp()
  - improve docstrings
  - clean-up in CMapping
  - update docstrings
  - rename (coor, weight) to (coors, weights)
  - move raise_if_too_large() back
  - fix PyCMapping.bfg computation
  - improve calculation efficiency - use dets_fast and invs_fast
  - use dets_fast() in invs_fast()
  - fix 2x2 inversion
  - test_install.py: tweak tolerance to accommodate different rounding errors

- merge pull request #947 from rc/fix-surface-region-check

  - fix region check in Term.setup_integration() for 1D, 2D

- merge pull request #950 from rc/fix-region-copy

  - fix Region.copy() to preserve region kind

- merge pull request #949 from vlukes/update_de_grad_term

  - de_grad term: grad can be multiplied by vector or matrix
  - more readable code

- merge pull request #951 from rc/pass-variable-to-setter

  - update FieldVariable._get_setter() to pass variable to setter functions
  - update examples for variable argument of setter functions
  - docs: update users guide for variable argument of setter functions

- merge pull request #952 from rc/fix-typo

  - fix typo in debug code in divgrad_build_gtg()

- merge pull request #953 from rc/seismic-load-example

  - new apply_ebc argument of Variables.set_state()
  - apply EBCs in prestep_fun(), poststep_fun() - update
    Problem.get_tss_functions()

    - prestep_fun(), poststep_fun() return vec
    - fixes seismic loading, hyperelastic.py with active_only set to False etc.

  - update prestep_fun(). poststep_fun() calls to return vec
  - new sfepy/examples/linear_elasticity/seismic_load.py
  - update docstring and comments in elastodynamic.py example
  - gen_gallery.py: add custom view for seismic_load.py example
  - clean up get_ebcs() in seismic_load.py
  - test seismic_load.py in test_examples()

- merge pull request #954 from vlukes/mappings_refactoring_ii

  - VolumeMapping and SurfaceMapping -> FEMapping
  - fix docstring of Field.get_mapping()
  - fields_base.py: remove six
  - FEMapping.get_mapping(): simplify handling of bf
  - fix naming: base function -> basis function
  - eval_mapping_data_in_qp(): fix docstring

- merge pull request #956 from rc/fix-mat-by-region-cell-facet, closes #955

  - fix get_constants() in ConstantFunctionByRegion for mixed cell/facet
    regions
  - make Region.get_cell_indices() consistent with its definition, simplify

- merge pull request #957 from vlukes/rename_integration

  - rename integration: volume->cell, surface->facet,
    surface_extra->facet_extra
  - update test_term_call_modes
  - rename integration II.
  - update fields and variables
  - remove unused argument 'integration' from get_econn()
  - update region.set_kind(): edge, face -> facet
  - fix test: volume -> cell
  - get_data_shape(): fix docstrings

- merge pull request #961 from rc/get-blocks-stats

  - new get_blocks_stats()
  - new test_get_blocks_stats()

- merge pull request #962 from vlukes/fix_link_to_releases

  - fix download link: -> https://github.com/sfepy/sfepy/tags

- merge pull request #958 from vlukes/unification_of_fefields

  - unification of VolumeField and SurfaceField, only FEField remains
  - rearange function ordering in FEField
  - fix format
  - replace reference to VolumeField and SurfaceField by FEField
  - FEField._setup_geoemtry(): fix region.kind check
  - enable einsum optimization
  - FEField: setup bubble dofs for surface field
  - docs: replace volume -> cell, surface -> facet
  - docs: fix typo

- merge pull request #964 from rc/make-any_dof_conn-configurable

  - update docstrings of Equations.get_graph_conns(), .create_matrix_graph() -
    volume -> cell, surface -> facet
  - add any_dof_conn argument to relevant functions, support it in conf.options

    - update eval_equations()
    - update Problem.update_equations(), .time_update(), .select_bcs(),
    .evaluate(), .eval_equations()

  - docs: update users guide

- merge pull request #966 from rc/resview-plot-on-step-change-only

  - return early from pv_plot() when step stays the same

- merge pull request #967 from rc/plot-times-tight-layout

  - plot_times.py: use tight layout

- merge pull request #969 from lokik/scalar_bar_limits

  - Allow one side of scalar_bar_limits to be computed in resview

- merge pull request #970 from vlukes/update_homogenization

  - update CoefExprPar class: variable setting is not required
  - fix recovery for 2D plates

- merge pull request #971 from rc/resview-fix-for-no-mat-id

  - fix pv_plot() for meshes without mat_id array

- merge pull request #968 from rc/update-elastodynamics, closes #972

  - remove Python 2 compatibility code
  - new TimeSteppingSolver.set_dof_info(), use in Problem.solve()
  - update gen_multi_vec_packing() for extra variables (multiphysics)
  - support extra variables in ElastodynamicsBaseTS

    - update .get_matrices(), .get_a0(), .get_initial_vec()
    - original unpacking supported in .get_initial_vec() for compatibility

  - update NewmarkTS for extra variables - remove ._create_nlst_a(), update
    .create_nlst(), .step()
  - simplify gen_multi_vec_packing() - update
    ElastodynamicsBaseTS.get_initial_vec(), NewmarkTS.step()
  - fix NewmarkTS.create_nlst()
  - update create_arg_parser() for d+name derivatives, new allow_derivatives
    arg - update create_bnf() (no change in behaviour)
  - add allow_derivatives arg to Term.setup(), .setup_formal_args()
  - update Term.assemble_to() to ignore d+name derivatives
  - do not setup terms twice when calling Equation.from_desc() - new setup
    argument of Equation.__init__()
  - new transform_equations_ed()
  - support transforming elastodynamics equations in Problem.set_solver()
  - update ElastodynamicsBasicTSC for extra variables

    - extra variables are not used for time step control
    - update .get_scaled_errors, .get_initial_dt(), .__call__()

  - new sfepy/examples/multi_physics/piezo_elastodynamic.py
  - test piezo_elastodynamic.py in test_examples()
  - gen_gallery.py: add custom view for piezo_elastodynamic.py example
  - docs: update latex_elements in conf.py
  - remove debug print in transform_equations_ed()
  - new auto_transform_equations option in Problem.set_solver(), transform only
    once
  - docs: update users guide
  - simplify seismic_load.py example by using equations transformation
  - fix Problem.set_solver() for setting nls
  - set default value of auto_transform_equations to False - update
    seismic_load.py, piezo_elastodynamic.py examples, docs
  - update docstring of apply_ebc_to_matrix()
  - apply EBCs to vec in Evaluator when active_only is False

    - update .eval_residual(), .eval_tangent_matrix()
    - the tangent matrix is usually computed just after the residual, so in
      that case applying EBCs there is not needed, but the overhead is minimal

  - update Problem.get_tss_functions() - do not apply EBCs in prestep_fun()
  - change default output directory in piezo_elastodynamic.py example
  - new active_only define() argument in elastodynamics examples - update
    elastodynamic.py, seismic_load.py, piezo_elastodynamic.py
  - set default value of allow_derivatives to False - update
    create_arg_parser(), Term.setup(), .setup_formal_args()
  - add allow_derivatives argument to Equations.from_conf(),
    Equation.from_desc() - update Problem.set_equations(), Terms.setup()
  - update Problem.set_solver() to check for solver.var_names
  - new ElastodynamicsBaseTS._common_parameters, use in subclasses
  - new var_names solver option in ElastodynamicsBaseTS._common_parameters -
    update .__init__(), .get_initial_vec()
  - add var_names to elastodynamics solver options in elastodynamic.py example
  - update docstring of elastodynamic.py example
  - make var_names optional in ElastodynamicsBaseTS._common_parameters - update
    .__init__(), .get_initial_vec()
  - add var_names to elastodynamics solver options in test_ed_solvers.py
  - new test_active_only() in test_ed_solvers.py - update get_ic() in define()
    for 2D

- merge pull request #975 from vlukes/fix_nonlinhomog_example

  - replace term.integration by term.act_integration
  - fix MPI logging

.. _2022.3-2022.4:

from 2022.3 to 2022.4
=====================

- merge pull request #880 from vlukes/update_hdf5meshio_write

  - update HDF5MeshIO.write(): dump region_name

- merge pull request #878 from rc/fix-warnings

  - fix warning in test_base_functions_delta()
  - fix file not closed in Mesh3DMeshIO.read()
  - store connectivity in XYZMeshIO.write() as integers, update xyz meshes
  - fix sparse matrix comparison warning in test_hdf5_meshio()
  - fix MumpsSolver.__init__()
  - fix MUMPSSolver.__init__() for MUMPS not installed
  - silence SparseEfficiencyWarning in apply_ebc_to_matrix()
  - update gmres callback setup in ScipyIterative.__call__()
  - silence SparseEfficiencyWarning in project_to_facets()
  - silence SparseEfficiencyWarning in ScipyDirect.presolve()
  - new pytest.ini
  - clean up coefficients.py
  - fix remaining invalid escape sequence warnings
  - update pytest.ini
  - fix ScipyIterative.__call__() for scipy < 1.4.0
  - update pytest.ini: do not raise error on warnings

- merge pull request #881 from rc/sfepy-test

  - new sfepy/scripts/run_tests.py
  - new sfepy-test entry point in setup.py
  - docs: update installation/testing instructions

- merge pull request #888 from rc/allow-mesh-in-place-filename

  - allow Mesh instance in place of filename in
    PDESolverApp.setup_output_info()

- merge pull request #892 from rc/fix-eig-arguments

  - rename method argument of eig() to solver_kind, fixes shadowing method
    parameter
  - update eig() calls

- merge pull request #887 from rc/no-save-in-corrector-calls

  - do not save results in problem.solve() calls in corrector classes - update
    CorrNN, CorrN, CorrOne, CorrEqPar

- merge pull request #891 from heczis/fix_simple_phonon_args, closes #883, #884

  - Fix typo in docstring
  - Do not plot without band gaps or dispersion
  - Make band-gaps exclusive to dispersion and phase-velocity
  - Fix AcousticBandGapsApp.setup_options
  - Use parser.error instead of raise
  - Make all band-gaps related options available simultaneously
  - Facelift plotting functions
  - Remove show and new_axes from plotting functions
  - Separate log files for band_gaps and dispersion

- merge pull request #895 from rc/fix-simplify-log-plot

  - use plt.tight_layout() in every LogPlotter.apply_commands() call
  - unregister Log.terminate() in atexit on successful termination
  - do not use threading in LogPlotter.__call__(), update .apply_commands()
  - use spawn start method instead of fork in Log.__call__() to start
    LogPlotter - update .__init__()

- merge pull request #897 from rc/fix-log-plotter-labels

  - plot labels and vlines in every LogPlotter.apply_commands() call - fixes
    not plotting labels and vlines when show_legends was False

- merge pull request #906 from rc/fix-sfepy-test-explicit-paths

  - fix test() for explicit paths passed to run_tests.py (sfepy-test)

- merge pull request #904 from flothesof/transient_poisson_multimaterial_example

  - first version of example
  - finish docstring of example
  - take into account comments from PR
  - adding post processing
  - fix line lengths
  - fix run examples sphinx formatting
  - make usage path relative to sfepy package dir like in other examples
  - adding example to list of examples to be tested
  - removing plt.show call to make example testable
  - add custom gallery view for time_heat_equation_multi_material.py example
  - change ordering in test

- merge pull request #908 from fix-pyvista-mesh-style

  - fix pv_plot() for pyvista 0.37.0
  - show edges with glyphs in pv_plot() only when requested explicitly

- merge pull request #902 from rc/elastodynamics-time-step-control

  - new TimeStepController
  - add tsc argument to TimeSteppingSolver.__init__(), .__call__()
  - new FixedTCS, ElastodynamicsBasicTCS in new
    sfepy/solvers/ts_controllers.py - update solver_table
  - create VariableTimeStepper in ElastodynamicsBaseTS.__init__() when
    adapting - for non-fixed TimeStepController subclasses
  - support time step control in GeneralizedAlphaTS, new .step(), split
    .__call__()
  - remove .dts attribute of VariableTimeStepper
  - de-struct arguments of
    {TimeStepController,ElastodynamicsBasicTCS}.__call__() - update
    GeneralizedAlphaTS.__call__()
  - new TimeStepController.get_initial_dt()
  - fix adaptive time stepping in GeneralizedAlphaTS.__call__() - do not update
    dt in-place when the step is accepted
  - add comment to ElastodynamicsBaseTS.__init__()
  - new TimesSequenceTCS, new .__init__(), get_initial_dt(), .__call__()
  - move ElastodynamicsBaseTS._create_nlst_*() to particular classes

    - move ._create_nlst_a() to NewmarkTS
    - move ._create_nlst_u() to BatheTS

  - split/move GeneralizedAlphaTS.__call__() into new
    ElastodynamicsBaseTS.__call__() - new GeneralizedAlphaTS.__init__(), update
    .step()
  - new define(), test_ed_solvers() in new
    sfepy/tests/test_elastodynamic_solvers.py
  - pass prestep_fun to .step() in ElastodynamicsBaseTS.__call__() - update
    GeneralizedAlphaTS.step()
  - new {VelocityVerletTS, NewmarkTS, BatheTS}.step() replacing .__call__() -
    use ElastodynamicsBaseTS.__call__() for all elastodynamics solvers
  - new ElastodynamicsBasicTCS.get_scaled_errors(), update .__call__()
  - fix typo (TCS -> TSC)
  - new ElastodynamicsPIDTSC, new .__init__(), .__call__() - update docstrings
  - rename test_elastodynamic_solvers.py -> test_ed_solvers.py
  - new ElastodynamicsBaseTS.clear_lin_solver(), update .__call__()
  - new BatheTS.clear_lin_solver() - fix for adaptive time step control
  - test time step controllers in test_ed_solvers(), update define()
  - mark test_ed_solvers() as slow
  - halve t1, clean up test_ed_solvers()
  - gen_solver_table.py: update for TimeStepController
  - docs: sync module index of developer guide with current sources
  - docs: add time step controllers table to users guide
  - update Problem.set_conf_solvers(), .init_solvers() for time step
    controllers

    - update .__init__()
    - support time step controllers in problem definition files

  - simplify test_ed_solvers() by using Problem.init_solvers() - rename
    _list_elastodynamic_solvers() -> _list_solvers()
  - illustrate time step control in elastodynamic.py example, new define()
  - gen_gallery.py: add custom view for elastodynamic.py example
  - use error_order parameter in ElastodynamicsBasicTSC.__call__() - update
    ElastodynamicsPIDTSC._parameters
  - new ElastodynamicsBasicTSC.get_initial_dt(), guess_dt0 parameter - new
    eval_scaled_norm(), update .get_scaled_errors()
  - tweak ElastodynamicsBasicTSC parameter description
  - new ElastodynamicsLinearTSC, new .__init__(), .__call__()
  - update test_ed_solvers(), define()
  - show ElastodynamicsLinearTSC use in elastodynamic.py example, update
    define()

- merge pull request #912 from rc/gallery-camera-position

  - gen_gallery.py: support camera_position option in custom settings - update
    resview_plot(), generate_images()

- merge pull request #909 from flothesof/helmholtz_example, closes #907

  - first version helmholtz
  - improve example following review
  - mesh improvement around source region
  - add to test_declarative_examples
  - format long lines
  - fix LaTeX in docstring
  - gen_gallery.py: add custom view for helmholtz_apartment.py example
  - merge materials into a single materials dictionary to reduce number of
    terms in equations
  - remove unused import

- merge pull request #913 from rc/misc-fixes

  - fix test_eigenvalue_solvers() to obey can_fail when solving
  - docs: sync module index of developer guide with current sources

- merge pull request #914 from rc/central-difference-tss

  - new CentralDifferenceTS - new ._create_nlst_a(), .create_nlst(), .step()
  - test CentralDifferenceTS in test_ed_solvers(), update define()
  - store u as function in CentralDifferenceTS like in other ED solvers
  - add ts.central_difference to elastodynamic.py example, fix line lengths
  - fix line length, update comments in define() of test_ed_solvers.py
  - update CentralDifferenceTS docstring

- merge pull request #915 from rc/fix-multilinear-terms-output-axes-order

  - fix undefined output axes order in ExpressionBuilder.build()

.. _2022.2-2022.3:

from 2022.2 to 2022.3
=====================

- merge pull request #839 from rc/fix-active-bcs

  - fix active BCs in EquationMap.map_equations() for no key set

- merge pull request #840 from rc/fix-test-docstring

  - fix docstring of test()

- merge pull request #841 from rc/improve-material-constructor-doc-pr737 -
  finishes and closes #737 by @burnpanck

  - tried to improve documentation of Material constructor. Made values= and
    kwargs completely equivalent instead of almost equivalent
  - fix Material.__init__() formatting
  - fix link in Material.__init__() docstring
  - document special material parameters in Material.__init__() docstring

- merge pull request #842 from rc/remove-simple-homog-mpi

  - simple.py: new bvp-mM mode, new --debug-mpi option
  - document application kinds in simple.py docstring
  - remove simple_homog_mpi.py
  - fix micro_filename path in nonlinear_hyperelastic_mM.py example
  - update materials in nonlinear_homogenization.py for new Material
    constructor

    - see #841 (#737)
    - fix def_mat()

  - fix get_slaves() for Python 3
  - docs: sync module index of developer guide with current sources
  - add basic description to nonlinear_hyperelastic_mM.py example

- merge pull request #843 rc/resview-any-sfepy-format

  - new make_grid_from_mesh(), use it in read_mesh()
  - resview.py: update read_mesh() to support any SfePy mesh format
  - distinguish sfepy and pyvista mesh variables in read_data()
  - resview.py: update print_camera_position() to work with subprocess.call()

- merge pull request #845 from rc/reorganize-scripts-entry-points

  - remove old scripts edit_identifiers.py, eval_ns_forms.py, eval_tl_forms.py
  - move script/convert_mesh.py -> sfepy/scripts/convert_mesh.py
  - new sfepy-convert entry point in setup.py
  - update docstring of convert_mesh.py
  - move simple.py -> sfepy/scripts/simple.py
  - new sfepy-run entry point in setup.py, remove sfepy-run script
  - update docstring of simple.py
  - move probe.py -> sfepy/scripts/probe.py
  - new sfepy-probe entry point in setup.py
  - update docstring of probe.py
  - move resview.py -> sfepy/scripts/resview.py
  - new sfepy-view entry point in setup.py
  - update docstring of resview.py
  - move extractor.py -> sfepy/scripts/extractor.py
  - update docstring of extractor.py
  - move mesh generator scripts to sfepy/scripts/

    - move script/blockgen.py -> sfepy/scripts/blockgen.py,
      script/cylindergen.py -> sfepy/scripts/cylindergen.py,
      script/gen_iga_patch.py -> sfepy/scripts/gen_iga_patch.py,

  - blockgen.py: split main() into new gen_block(), add_args()
  - cylindergen.py: split main() into new gen_cylinder(), add_args()
  - gen_iga_patch.py: split main() into new gen_iga_patch(), add_args()
  - new sfepy/scripts/gen_mesh.py
  - new sfepy-mesh entry point in setup.py
  - convert_mesh.py: new --extract-edges option in main()

    - move merge_lines(), extract_edges() into mesh_tools.py
    - remove extract_edges.py

  - convert_mesh.py: new --extract-surface, --print-surface options in main()

    - move _get_facets(), get_surface_faces(), surface_graph(),
      surface_components() into mesh_tools.py
    - remove extract_surface.py

  - convert_mesh.py: new --tile option in main()

    - remove tile_periodic_mesh.py

  - move script/show_mesh_info.py -> sfepy/scripts/show_mesh_info.py
  - show_mesh_info.py: split main() into new show_mesh_info(), add_args()
  - gen_mesh.py: add show_mesh_info.py as subcommand in main()
  - move script/dg_plot_1D.py -> sfepy/examples/dg/dg_plot_1D.py
  - update DG examples for moved dg_plot_1D.py
  - move script/save_basis.py -> sfepy/scripts/save_basis.py
  - move script/plot_* scripts into sfepy/scripts/
  - move script/gen_mesh_prev.py -> sfepy/scripts/gen_mesh_prev.py

    - to be integrated with resview.py after fixing

  - rename script/ -> tools/ (not-to-be-installed/development scripts)
  - move build_helpers.py -> tools/build_helpers.py
  - move test_install.py -> sfepy/scripts/test_install.py
  - remove homogen.py compatibility script
  - update docstrings for script/ -> tools/
  - new sfepy/scripts/__init__.py, update sfepy/setup.py
  - update setup.py for reorganized scripts, update get_basic_info()

    - do not install unnecessary files

  - test_install.py: update main() for reorganized scripts
  - update MANIFEST.in
  - fix get_basic_info() to check for LICENSE instead of generated VERSION
  - update installed/sdist contents in setup.py, MANIFEST.in
  - update sphinx extension paths in doc/conf.py
  - sync_module_docs.py: update for reorganized scripts
  - docs: sync module index of developer guide with current sources
  - simple.py: fix phonon options handling in main()
  - docs: update INSTALL, README.rst for reorganized scripts
  - docs: update developer guide for reorganized scripts
  - docs: update FAQ for reorganized scripts
  - docs: update installation doc for reorganized scripts
  - docs: update preprocessing doc for reorganized scripts
  - docs: update primer for reorganized scripts
  - docs: update tutorial for reorganized scripts
  - docs: update release tasks for reorganized scripts
  - docs: update users guide for reorganized scripts
  - docs: update requirements list
  - update docstrings of examples for reorganized scripts
  - gen_gallery.py: update for reorganized scripts
  - docs: make using salome doc obsolete, update for reorganized scripts
  - docs: update man pages (sfepy-run, sfepy-view)

- merge pull request #849 from rc/update-authors

  - update AUTHORS

- merge pull request #850 from rc/materials-fixes

  - pass ts to .update_material() in Problem.get_tss_functions()'s
    prestep_fun()
  - docs: fix FAQ example
  - docs: auto-document special methods to have Material.__init__() in docs
  - docs: fix typo in FAQ

- merge pull request #851 from rc/fix-multilinear-terms-with-computed-arguments

  - fix ETermBase.make_function() for repeated call of de_non_penetration_p
    term - treat correctly ExpressionArg arguments

- merge pull request #852 from rc/fix-examples

  - fix default output_dir in define() of stokes_slip_bc.py example
  - improve plots in main() of linear_viscoelastic.py example
  - fix solve_problem() in shell10x_cantilever_interactive.py example - use
    AutoDirect solver
  - clean up several examples

- merge pull request #856 from vlukes/fix_splineregion2d, closes #854

  - fix SplineRegion2D.find_ts(): depreciated 'xtol' parameter

- merge pull request #855 from rc/fix-test-install-parsing

  - test_install.py: fix parsing in report_tests()

- merge pull request #853 from rc/support-pip-install, closes #320

  - new pyproject.toml, update MANIFEST.in
  - specify install_requires in setup_package() of setup.py
  - docs: update release tasks
  - docs: update installation instructions
  - docs: fix venv testing instructions in release tasks
  - docs: simplify installation instructions
  - docs: fix/update TestPyPI upload instructions in release tasks

- merge pull request #858 from rc/add-missing-tutorial-image

  - docs: add missing tutorial image

- merge pull request #859 from rc/fix-for-docs-generation

  - use absolute import in test_dg_terms_calls.py
  - prepend data_dir to mesh name in recovery_micro() of
    piezo_elasticity_micro.py
  - new define() in linear_viscoelastic.py example, silence output on import -
    update main()

- merge pull request #861 from vlukes/update_mumps_interface, closes #860

  - update interface for mumps library version 5.4.x
  - fix for library version 5.3

- merge pull request #862 from rc/speedup-ls-reuse-unify-presolve

  - new LinearSolver.clear()
  - new use_mtx_digest parameter in ScipyDirect, new .clear()

    - update .__init__(), .__call__(), .presolve()
    - new use_mtx_digest parameter in ScipySuperLU, ScipyUmfpack
    - update presolve behaviour

  - new use_mtx_digest parameter in MUMPSSolver, new .clear()

    - update .__init__(), .__call__(), .presolve(), .__del__()
    - new use_mtx_digest parameter in ScipySuperLU, ScipyUmfpack
    - update presolve behaviour

  - clear linear solver in ElastodynamicsBaseTS.get_a0() - allows factorization
    reuse without matrix digest (conf.use_mtx_digest = False)
  - use separate linear solver instances in BatheTS

    - update .__init__(), .create_nlst1(), .create_nlst2()
    - allows factorization reuse without matrix digest (conf.use_mtx_digest =
      False)

  - use and describe use_mtx_digest in elastodynamic.py example

- merge pull request #863 from rc/fix-generalized-alpha-tss

  - fix GeneralizedAlphaTS: missing M coefficient, wrong state update

    - new ._create_nlst_a()
    - update .create_nlst(), .__call__()

- merge pull request #864 from rc/decrease-sdist-size

  - gen_gallery.py: set screenshot size to (800, 600) in resview_plot()
  - exclude PDF docs (sfepy_manual.pdf) from source tarball

- merge pull request #866 from rc/docs-update-enthought-link, closes #688

  - docs: replace Enthought Canopy link by Enthought Deployment Manager

- merge pull request #867 from rc/fix-primer-interactive, closes #629

  - docs: fix/update interactive mode part of primer

- merge pull request #868 from vlukes/update_homog_recovery

  - update CorrSolution: new get_output() for exporting states
  - unify recover_micro_hook() and recover_micro_hook_eps()
  - new 'div' and 'cauchy_strain' modes of fields.evaluate_at()
  - update piezo-elasticity example
  - fix recover_micro_hook(): slice -> nm.arange
  - update linear elasticity example
  - fix micro cell center
  - fix format
  - add docstring to 'recover_micro_hook()'
  - new mechanics.tensors.get_cauchy_strain()
  - region_mode and eval_mode as define() arguments
  - update docstring of recover_micro_hook()
  - fix typo
  - add docstring to get_causchy_strain()

- merge pull request #869 from vlukes/save_files_per_region

  - update save_state(): split file by region, file_per_var -> file_split_by
  - update examples
  - fix save_state()
  - change `file_split_by` to `split_results_by`

- merge pull request #870 from rc/gen-gallery-all-sources

  - update gen_gallery.py to include all examples sources in docs

    - new omit_images
    - update _omit(), generate_images(), generate_rst_files()

- merge pull request #871 from vlukes/fix_save_files_by_region

  - fix 'extend' value in save_state()
  - do not split h5 files

- merge pull request #872 from rc/fix-custom-views-files-by-region

  - fix custom views in gen_gallery.py for files by region

- merge pull request #873 from vlukes/update_evaluate_at

  - update evalauate_at(): return three-dimensional array

- merge pull request #874 from vlukes/fix_recovery_micro_hook

  - use region_name if available

.. _2022.1-2022.2:

from 2022.1 to 2022.2
=====================

- merge pull request #806 from vlukes/fix_meshio_vtkint64

  - fix meshio: int64 required by the newer paraview

- merge pull request #805 from vlukes/update_ev_div_grad_terms

  - new EGradTerm
  - update DivTerm, GradTerm: new optional material parameter

- merge pull request #809 from rc/de-lin-convect-term new elinearconvectterm
  term (de_lin_convect)

  - new ELinearConvectTerm term (de_lin_convect)

- merge pull request #811 from vlukes/suppres_meshio_warnings

  - suppress meshio warnings

- merge pull request #812 from vlukes/suppres_meshio_warnings suppress meshio
  warnings - update

  - suppress meshio warnings - update

- merge pull request #803 from rc/use-pytest, closes #363

  - move utility TestCommon methods into module scope - move report(),
    eval_coor_expression(), compare_vectors(), assert_equal()
  - update test_assembling.py for pytest
  - update test_base.py for pytest, clean up
  - update test_cmesh.py for pytest
  - new conftest.py: new output_dir fixture
  - import sfepy.base.testing as tst instead of st
  - update test_conditions.py for pytest
  - remove run_tests_py related code in test_dg_field.py, clean up
  - remove run_tests_py related code in DGTermTestEnvornment - prepare
    DGTermTestEnvornment using fixture
  - remove TestCommon, fix compare_vectors(), assert_equal()
  - new run_declaratice_example(), remove tests/tests_basic.py

    - move NLSStatus from tests/tests_basic.py to sfepy/base/testing.py
    - make TestInput.check_conditions() standalone function
    - remove TestDummy, TestInput, TestInputEvolutionary
    - remove TestLCBC (WIP)

  - new test_declarative_examples.py: new test_examples(), test_examples_dg()
  - remove tests/test_input*, tests/test_dg_input* - replaced by
    test_declarative_examples.py
  - fix error checking in updated tests - fix test_assembling.py, test_base.py,
    test_cmesh.py, test_conditions.py
  - remove return statements in test_dg_terms_calls.py
  - update test_domain.py for pytest
  - update test_eigenvalue_solvers.py for pytest
  - update test_elasticity_small_strain.py for pytest
  - update test_fem.py for pytest
  - update test_functions.py for pytest
  - update test_high_level.py for pytest
  - update test_homogenization_engine.py for pytest
  - update test_homogenization_perfusion.py for pytest
  - update test_hyperelastic_tlul.py for pytest
  - update test_io.py for pytest, clean up
  - update test_laplace_unit_disk.py for pytest, clean up
  - update test_laplace_unit_square.py for pytest, clean up
  - fix LCBCOperators.add_from_bc() for no arguments
  - support output_dir in define() of stokes_slip_bc.py example
  - support output_dir in run() of laplace_shifted_periodic.py example
  - update test_lcbcs.py for pytest, use output_dir for all output
  - new test_elasticity_rigid() in test_lcbcs.py
  - remove tests/test_lcbc_2d.py, tests/test_lcbc_3d.py - replaced by
    test_elasticity_rigid()
  - update MultiProblem solver to use master problem output_dir - update
    .init_subproblems(), .__call__()
  - do not save mesh in gen_two_bodies() and log in two_bodies_contact.py
    example
  - do not save state twice in test_ebc_functions(), update
    test_region_functions()
  - new save_lcbc_vecs argument of define() in stokes_slip_bc.py example
  - update test_stokes_slip_bc() for save_lcbc_vecs argument
  - do not save state twice in test_solving()
  - update test_linalg.py for pytest
  - update test_linear_solvers.py for pytest
  - update test_linearization.py for pytest
  - update test_log.py for pytest
  - update test_matcoefs.py for pytest
  - update test_mesh_expand.py for pytest
  - update test_mesh_generators.py for pytest
  - update test_mesh_interp.py for pytest
  - update test_mesh_smoothing.py for pytest
  - update test_meshio.py for pytest
  - update test_msm_laplace.py for pytest, clean up
  - update test_msm_symbolic.py for pytest
  - update test_normals.py for pytest
  - update test_parsing.py for pytest
  - check EPBC application in test_epbcs(), fix error report
  - remove tests/test_periodic_bc_2d.py, tests/test_periodic_bc_3d.py - covered
    by test_epbcs()
  - update test_poly_spaces.py for pytest
  - do not save state twice in _test_msm_symbolic()
  - fix return in test_hdf5_meshio()
  - update test_projections.py for pytest
  - update test_quadratures.py for pytest
  - update test_ref_coors.py for pytest
  - update test_refine_hanging.py for pytest
  - update test_regions.py for pytest
  - fix ls_red parameter default in SemismoothNewton parameters
  - update test_semismooth_newton.py for pytest
  - update test_sparse.py for pytest
  - update test_splinebox.py for pytest
  - update test_tensors.py for pytest
  - update test_term_call_modes.py for pytest
  - update test_term_consistency.py for pytest, fix test_surface_evaluate() -
    clean up
  - update test_term_sensitivity.py for pytest
  - update test_units.py for pytest
  - update test_volume.py for pytest
  - warnings clean up: do not use numpy.int, numpy.bool
  - warnings clean up: do not use numpy.object, numpy.float
  - warnings clean up: update imports in ScipyDirect, ScipyIterative
  - warnings clean up: update import in LOBPCGEigenvalueSolver
  - warnings clean up: remove unused arpack_eigs()
  - warnings clean up: update and simplify dot_sequences() - do not use
    numpy.core.umath_tests.matrix_multiply()
  - do not print warp violation messages - only set error in _v_describe(),
    dq_finite_strain(), dq_tl_finite_strain_surface()
  - update test_install.py for pytest, clean up - update report_tests(), main()
  - new pytest_addoption(): new --output-dir option for output_dir fixture
  - new pytest_configure(): new slow marker
  - mark slow tests
  - new test() in sfepy/__init__.py for installed sfepy testing
  - remove run_tests.py, update setup.py, sfepy-run
  - move tests/ -> sfepy/tests/
  - obey output_dir in test_solution() of test_homogenization_perfusion.py
  - add output_dir argument to homogenization-based material functions - update
    get_homog_coefs_linear(), get_homog_coefs_nonlinear()
  - add output_dir argument to homogenization micro recovery functions - update
    recover_micro_hook_init(), recover_micro_hook(), recover_micro_hook_eps()
  - obey output_dir in micro-problem of linear_elastic_mM.py example - update
    post_process(), get_homog()
  - verify that output_dir is directory in output_dir() fixture
  - docs: update testing info
  - new remove_prefix argument of run_declaratice_example()
  - update test_examples(), test_examples_dg() to remove common name prefix

    - new inedir(), examples_dir
    - remove common prefix in examples, examples_dg lists

  - update docstrings of postproc.py, resview.py for renamed result files and
    pytest
  - docs: update users guide for pytest-related changes
  - docs: update developer guide for moved tests/
  - move examples/ -> sfepy/examples/
  - update examples for moved examples/
  - assume example files relative to sfepy.base_dir in
    run_declaratice_example() - add comment to examples_dir in
    test_declarative_examples.py
  - use sfepy.base_dir in multi-file examples
  - update tests for moved examples/
  - update test_install.py for moved tests/ and examples/
  - update gen_gallery.py, gen_term_table.py for moved examples/
  - update setup.py for moved tests/ and examples/
  - docs: update FAQ for moved examples/
  - update docstrings of examples for moved examples/
  - use sfepy.base_dir in more multi-file examples
  - docs: update for moved examples/
  - update docstrings of simple.py, simple_homog_mpi.py for moved examples/
  - sync_module_docs.py: omit sfepy/examples/ directory
  - docs: sync module index of developer guide with current sources, include
    tests
  - update MANIFEST.in for moved examples/
  - CI: update appveyor config for pytest, use Python 3.9
  - test_install.py: update some output directories

- merge pull request #814 from rc/update-testing-docs

  - docs: update testing instructions

- merge pull request #815 from rc/elastic-wave-speeds

  - new youngpoisson_from_wave_speeds(), wave_speeds_from_youngpoisson()
  - new test_wave_speeds()

- merge pull request #818 from heczis/master

  - Fix legend placement

- merge pull request #820 from rc/remove-phonon

  - run AcousticBandGapsApp from simple.py, new --app option

    - add phonon.py-specific options with --phonon- prefix
    - determine app_mode in main()

  - remove phonon.py
  - update scripts for no phonon.py
  - docs: update for no phonon.py

- merge pull request #821 from rc/fix-multilinear-terms-for-complex-args

  - fix multilinear terms for complex arguments - new ETermBase._eval(),
    .eval_real(), .eval_complex()

- merge pull request #823 from vlukes/issue_760, closes #760

  - update posprocessing instructions (#760)
  - examples: update posprocessing instructions (#760)
  - linear_viscoelastic example: plot deformed mesh

- merge pull request #824 from vlukes/resview_isosurface, closes #759

  - resview isosurfaces (#759)

- merge pull request #826 from vlukes/update_resview

  - resview.py: new 'camera_position' and 'window_size' options
  - fix option names: replace '_' by '-'
  - update print of camera position

- merge pull request #828 from vlukes/resview_grid

  - resview.py: new 'camera_position' and 'window_size' options
  - fix option names: replace '_' by '-'
  - update print of camera position
  - resview.py: plot figures in grid
  - fix spelling
  - fix settings of position vectors
  - rename '--position-vector' to '--grid-vector'
  - update users' guide

- merge pull request #829 from vlukes/remove_mayavi_dependency, closes #108,
  #687

  - remove mayavi dependent files
  - update documentation
  - update setup
  - update 'sfepy-run'
  - update INSTALL instructions
  - update tutorials
  - resview.py: report saving figure to file
  - test resview.py script
  - remove doc files
  - update examples

- merge pull request #830 from rc/fix-eval-complex

  - fix weak mode in eval_complex() for two complex arguments

- merge pull request #831 from vlukes/term_docs_cleanup, closes #773

  - update term docs - simplify descriptions
  - update terms notation
  - update docs: define evaluation modes
  - update 'ev_integrate_mat' docs: remove 'el_avg', 'qp' modes
  - update term docs: replace '\bar{p}' with 'p'

- merge pull request #832 from heczis/fix_omit_facets

  - Fix omit_facets

- merge pull request #833 from vlukes/evaluable_EScalarDotMGradScalarTerm

  - Update EScalarDotMGradScalarTerm to be evaluable

- merge pull request #834 from vlukes/fix_homog_phononic

  - Fix compute_eigenmomenta()
  - fix: init state on compute_eigenmomenta()

- merge pull request #835 from vlukes/fix_multiproblem_output

  - Fix output setting in subproblems

- merge pull request #837 from rc/fix-test-examples-iga

  - fix igakit-dependent examples list - add navier_stokes2d_iga.py
  - return pytest result from test()

- merge pull request #836 from rc/fix-gen-gallery-resview

  - resview.py: fix pv_plot() for 1D meshes
  - script/gen_gallery: fix generate_images() for current resview.py - update
    custom plots

.. _2021.4-2022.1:

from 2021.4 to 2022.1
=====================

- merge pull request #762 from rc/newton-step-reduction

  - new step_red parameter of Newton solver, update .__call__()

- merge pull request #765 from rmcgibbo/patch-1

  - Fix for python3.10 collections.Sequence -> collections.abc.Sequence See for
    example: https://github.com/wireservice/agate/pull/737/files

- merge pull request #767 from rc/fix-auto-solver-dict-conf

  - fix AutoFallbackSolver.__new__() for dict conf

- merge pull request #768 from rc/fix-its2d_5-example

  - fix stress_strain() in its2D_5.py example

- merge pull request #764 from rc/move-state-into-variables, closes #378

  - rename stripped -> reduced, strip -> reduce in DOF vector context

    - rename Equations.create_stripped_state_vector() ->
      .create_reduced_state_vector()
    - rename Equations.strip_state_vector() -> .reduce_vec()
    - rename Variables.create_stripped_state_vector() ->
      .create_reduced_state_vector()
    - rename Variables.strip_state_vector() -> .reduce_vec()
    - update Variables.get_indx(), .check_vector_size(),
      .get_state_part_view(), .set_state_part() argument names
    - update affected code

  - new Variables.vec, .rvec attributes to hold state vector, new .init_state()

    - update .apply_ebc(), .apply_ic(), .has_ebc() - vec argument optional, use
      self.vec by default

  - new Equations.init_state(), update .apply_ebc(), .apply_ic()
  - remove Variables.set_from_state(), update Oseen.__call__() (its only use)
  - remove unused Equations.state_to_output()
  - rename state_vector -> vec in DOF vector context

    - rename Equations.create_state_vector() -> .create_vec(),
      .create_reduced_state_vector() -> .create_reduced_vec()
    - rename Variables.create_state_vector() -> .create_vec(),
      .create_reduced_state_vector() -> .create_reduced_vec()
    - update affected code

  - rename state_part{_view} -> vec_part in DOF vector context

    - rename Variables.get_state_part_view() -> .get_vec_part(),
      .set_state_part() -> .set_vec_part()
    - update affected code

  - rename Variables.state_to_output() -> .create_output() - update affected code
  - remove unused variable in FieldVariable.save_as_mesh()
  - fix Variable.advance() for multi-variable problems - replace deque .data
    attribute by list
  - add step argument to Variable.set_constant(), use in .init_data()
  - new .locked attribute of Variable, update .set_data()
  - implement state handling in Variables, regroup/rename functions

    - lock state variables in .init_state()
    - new .fill_state(), .invalidate_evaluate_caches()
    - invalidate caches in .apply_ebc(), .apply_ic()
    - new .set_reduced_state(), .get_reduced_state(), .set_full_state(),
      .set_state(), .get_state(), .set_state_parts()
    - update .get_state_parts(), .create_output()
    - rename .check_vector_size() -> .check_vec_size()

  - update Equations for Variables.set_state()

    - remove Equations.set_variables_from_state(), .get_state_parts()
    - use new .set_state() in .eval_residuals(), .eval_tangent_matrices()
    - update .init_state()

  - update Problem for Variables replacing State

    - initialize ics in .__init__(), .create_subproblem(), .set_equations()
    - remove .setup_ics()
    - update .create_state() to call .get_initial_state()
    - update .save_state(), .save_ebc(), .get_tss_functions(),
      .get_initial_state(), .solve(), .block_solve(), .save_restart(),
      .load_restart()

  - update affected code for Variables replacing State

    - update EVPSolverApp.save_results()
    - update Evaluator.new_ulf_iteration()
    - update create_mass_matrix()
    - update .__call__() of CorrNN, CorrN, CorrSetBCS, CorrEqPar,
      PressureRHSVector
    - update MultiProblem.init_subproblems(), .__call__()

  - update tests for Variables replacing State
  - update examples for Variables replacing State
  - do not import State in sfepy/discrete/__init__.py
  - docs: update users' guide
  - update SimpleEVP.save() for Variables replacing State
  - remove sfepy/discrete/state.py
  - docs: sync module index of developer guide with current sources
  - update yet more examples for Variables replacing State
  - docs: update FAQ for Variables replacing State
  - docs: update primer for Variables replacing State, Python 3 and current
    code
  - simplify and clean up verify_incompressibility()
  - replace stripped -> reduced in SchurEVP.post_process()
  - docs: update tutorial, sync it with linear_elastic_interactive.py example

- merge pull request #772 from rc/plot-mesh-kwargs

  - add color, **plot_kwargs arguments to plot_mesh()

- merge pull request #775 from rc/pass-eterm-options

  - add eterm_options argument to Problem.evaluate() - update
    .create_evaluable(), create_evaluable()
  - update ETermBase.make_function() for recent jax (0.2.21)

- merge pull request #777 from rc/check-errors-in-mesh-graph

  - check errors in mesh_graph(), exit early

- merge pull request #774 from vlukes/new_sensitivity_terms

  - new multilinear term: ENonSymElasticTerm - nonsymmetric gradient
  - new multilinear sensitivity term: ESDLinearElasticTerm
  - new multilinear sensitivity term: ESDPiezoCouplingTerm
  - new multilinear terms: EDiffusionTerm, ESDDiffusionTerm
  - new multilinear sensitivity term: ESDStokesTerm
  - new multilinear sensitivity term: ESDDivGradTerm
  - new multilinear sensitivity term: ESDDotTerm
  - update docstrings of multilinear terms, see #773
  - update docstrings of sensitivity terms, see #773
  - update term notation table
  - update term table generator - compact virtual/state/parameter arguments
  - rename 'parameter_mesh_velocity' to 'parameter_mv'
  - new multilinear terms: ELinearTractionTerm, ESDLinearTractionTerm
  - fix term docstring: EIntegrateOperatorTerm
  - move multilinear sensitivity terms to 'terms_senstitivity.py'
  - update test_term_sensitivity: use smaller mesh and gc.collect()
  - multilinear terms: small fixes
  - multilinear terms: use `v()` for nonsymmetric vector storage
  - use ':' only for symmetric gradient

- merge pull request #779 from rc/fix-polyspace-for-0d

  - fix shape check of quadrature coordinates in PolySpace.eval_base() for 0D -
    fixes surface integrals in 1D
  - fix BernsteinSimplexPolySpace._eval_base() for 0D

- merge pull request #780 from rc/fix-problem-solve-initial-state-arg

  - fix initial state argument handling in Problem.solve(), .block_solve() -
    new .set_default_state()
  - use Problem.set_default_state() where suitable
  - add docstring to Problem.set_default_state()

- merge pull request #781 from rc/update-faq-structural

  - docs: update FAQ with structural elements entry
  - docs: fix indentation in FAQ

- merge pull request #783 from rc/multilinear-terms-use-actual-ndarray-args -
  closes #782

  - use actual ndarray arguments in ETermBase.make_function()
  - clean up whitespace

- merge pull request #784 from rc/update-faq-element-matrix

  - docs: update FAQ with element matrix entry

- merge pull request #786 from rc/fix-balloon-example

  - fix extract_time_history() for Python 3
  - update plot_radius() in balloon.py example for recent matplotlib

- merge pull request #787 from rc/array-values-apply-unit-multipliers

  - support array-like values in apply_unit_multipliers()
  - remove unused elastic_constants_relations

- merge pull request #788 from rc/fix-resview-inf-factor

  - fix pv_plot() for data with zero norm

- merge pull request #789 from vlukes/update_div_grad

  - DivGradTerm: allow different aprroximations

- merge pull request #791 from vlukes/fix_div_grad_term fix divgrad_build_gtg()

  - fix divgrad_build_gtg()

- merge pull request #790 from rc/multilinear-terms-use-actual-mapping-data

  - use actual reference mapping data in get_einsum_ops()

    - do not store the mapping data in ExpressionArg
    - update ExpressionArg.from_term_arg(), .get_dofs()
    - new ExpressionArg.get_bf()
    - update ETermBase.make_function()

- merge pull request #793 from rc/remove-homogen

  - run HomogenizationApp from simple.py - update main()
  - remove homogen.py
  - update scripts for no homogen.py
  - docs: update for no homogen.py
  - new homogen.py proxy for compatibility, calls simple.py - re-add it to
    setup.py

- merge pull request #797 from vlukes/update_to_meshio5

  - update fem/meshio.py to meshio5
  - remove meshio version

.. _2021.3-2021.4:

from 2021.3 to 2021.4
=====================

- merge pull request #729 from brylie/patch-1

  - Include example Docker compose file for simplicity
  - Use `code/` subdirectory inside of container

- merge pull request #735 from burnpanck/patch-1

  - Clarify orientation of geometry elements
  - Added a note clarifying that the figure showing geometry elements applies
    to right-handed coordinate systems.

- merge pull request #738 from rc/fix-plot-condition-numbers

  - script/plot_condition_numbers.py: fix division in main()
  - script/plot_condition_numbers.py: resolve numpy warning in main()

- merge pull request #741 from rc/plot-mass-condition-numbers

  - script/plot_condition_numbers.py: add smass, vmass types to --matrix option

    - use choices in --matrix
    - remove useless dw_laplace coefficient in main()

  - script/plot_condition_numbers.py: remove order_fix in main()
  - script/plot_condition_numbers.py: new --output-dir, --no-show options
  - script/plot_condition_numbers.py: tweak plot parameters in main()

- merge pull request #734 from rc/resview-sfepy-h5-no-steps

  - support sfepy HDF5 meshes with no time steps in read_mesh() - new
    add_mat_id_to_grid()
  - update XDMF support in read_mesh() for current meshio, fix 2D meshes
  - change add_mat_id_to_grid() argument to cell_groups, update read_mesh()
  - determine default --position-vector from mesh bounding box in pv_plot()

    - set --position-vector default to None in main()
    - fixes placement of plots for non-3D meshes

- merge pull request #742 from rc/create-faq, closes #413, #706

  - new doc/faq.rst
  - docs: link faq.rst in index.rst
  - docs: add faq.rst to documentation.rst contents

- merge pull request #744 from vlukes/update_faq

  - new advice: use msh22 instead of msh4

- merge pull request #747 from rc/fix-plot-mesh-3d

  - script/plot_mesh.py: do not scale axis in 3D - fixes main() for current
    matplotlib (3.4.3)

- merge pull request #750 from vlukes/multiproc_recovery

  - update multiproc_proc module: import Pool
  - new parallel recovery of multiple microstructures

- merge pull request #748 from rc/update-faq

  - docs: add examples of regions by functions (imperative API) into FAQ
  - docs: add visualization of various FEM-related information into FAQ
  - docs: fix cell region definition in FAQ
  - docs: add PYTHONPATH setting for inplace builds into FAQ

- merge pull request #751 from flothesof/master

  - "dry water" flow example

- merge pull request #752 from vlukes/reinit_homogenization_engine

  - allow to reinitialize the homogenization engine

- merge pull request #753 from rc/pass-kwargs-conf-to-auto-solver

  - respect keyword arguments in 'auto' solvers - pass kwargs from
    AutoFallbackSolver.__new__() to use_first_available()
  - fix AutoFallbackSolver.__new__() to use the original solver configuration

- merge pull request #754 from vlukes/resview_streamlines

  - new streamlines plotting
  - new option for the XY plane view
  - update visualization parameters in docstring
  - relabel option

- merge pull request #749 from zitkat/dg-multiple-limiters

  - Allow multiple limiters for multiple fields.
  - Fix sfepy imports and remove hardcoded verbose flag.

- merge pull request #755 from rc/resview-fix-streamlines

  - resview.py: fix streamline vectors name in pv_plot()
  - resview.py: use sphere streamlines source in 3D in pv_plot()

- merge pull request #756 from vlukes/resview_gallery

  - update scalar bars positioning and range
  - update files caching
  - default glyphs for vector fields, new empty black field
  - use resview for gallery generation
  - adjustment of some gallery plots

.. _2021.2-2021.3:

from 2021.2 to 2021.3
=====================

- merge pull request #709 from vlukes/volume_surface_terms

  - update term integration: guess integration from region.kind
  - update terms: unify `ev_/dw_` volume/surface terms
  - update examples, tests, etc.: new term names
  - update terms: unify `dw_volume/surface_dot` terms
  - update examples, tests, etc.: `dw_dot`
  - new compatibility file: allow for old term names
  - rename term: `d_region` --> `ev_volume`
  - update term docstrings: `\cal{D}` stands for a volume or surface domain
  - unify multilinear terms
  - update Develeoper Guide
  - fix User's Guide

- merge pull request #711 from rc/volume-surface-terms-follow-up

  - docs: use only unified volume/surface terms in users' guide
  - fix docstring of poisson_functions.py example
  - use dw_dot and dw_integrate in interactive examples

- merge pull request #710 from vlukes/rename_d_terms

  - rename evaluation terms starting with `d_` to `ev_`
  - update user's guide: Term Overview
  - rename term: `ev_sd_volume_dot` to `ev_sd_dot`

- merge pull request #714 from rc/evp-solvers-fixes

  - fix QuadraticEVPSolver.__init__() for no solver in conf (use default)
  - fix MatlabEigenvalueSolver.__call__() for None value of eigenvectors arg

- merge pull request #716 from livoire13/example-interactive-iga

  - new iga-interactive example

- merge pull request #715 from rc/web-update-citing

  - web: update citing information
  - web: update related projects
  - web: improve plain text citation

- merge pull request #718 from rc/update-test-install

  - fix imperative_burgers_1D.py for unified surface/volume terms - use
    Term.new() in main()
  - imperative_burgers_1D.py: new --output-dir, --plot options - new
    parse_args(), update main()
  - laplace_iga_interactive.py: rename --output_dir option to --output-dir
  - test_install.py: test laplace_iga_interactive.py, imperative_burgers_1D.py

- merge pull request #719 from vlukes/fix_evp_corrector

  - fix TCorrectorsViaPressureEVP.setup_equations()

- merge pull request #721 from rc/fix-fill-in-message

  - fix fill-in message in Equations.create_matrix_graph(),
    Mesh.create_conn_graph()

- merge pull request #724 from rc/include-nls-tolerances

  - include tolerances in error checking in conv_test() - fixes 'or' eps mode
    with eps_a == 0.0 and zero initial residual

- merge pull request #725 from vlukes/fix_resview fix resview

  - fix format
  - fix colorbar position

- merge pull request #726 from rc/resview-auto-factor

  - resview.py: allow automatic percentage-based scaling factor in pv_plot() -
    update parse_options() for '%' character
  - resview.py: update help message and docstring
  - resview.py: use sentences in docstring, facilitate copy/paste of examples
  - resview.py: fix factor calculation in pv_plot()
  - resview.py: report scaling factor in pv_plot()

- merge pull request #727 from rc/resview-sfepy-h5

  - resview.py: new make_cells_from_conn(), used in read_mesh()
  - resview.py: support custom sfepy .h5, .h5x formats in read_mesh()
  - support h5x format in Problem.setup_output()

.. _2021.1-2021.2:

from 2021.1 to 2021.2
=====================

- merge pull request #680 from vlukes/update_doc

  - new example application
  - add supporting project

- merge pull request #681 from vlukes/new_example_app

  - new example application

- merge pull request #678 from rc/remove-med-mesh-files

  - remove unused meshes/various_formats/{med_2d_tri_quad.med,
    med_3d_tet_hex.med}
  - update filename_meshes in tests/test_meshio.py

- merge pull request #682 from vlukes/new_sd_term

  - d_sd_div_grad: only one optional material paremeter
  - update from_term_arg(): allow for (ndarray, arg_name) arguments
  - new SDPiezoCouplingTerm(ETermBase) sensitivity term
  - new tests of sensitivity terms

- merge pull request #683 from vlukes/make_de_convect_evaluable

  - update multilinear terms: make de_convect evaluable
  - update sensitivity tests: enable test of convect term
  - update docstring of de_convect term

- merge pull request #684 from vlukes/new_d_sd_surface_ltr

  - fix geme_mulAVSB3py() for calculating over all cells
  - new sensitivity term: SDLinearTractionTerm
  - update tests to check SDLinearTractionTerm
  - fix formating
  - update docstring for SDLinearTractionTerm() and rename: nel, nqp -> n_el,
    n_qp

- merge pull request #690 from rc/no-mat-update-in-init-time

  - do not update materials in Problem.init_time() - allows calling
    Problem.evaluate() in material update functions in time zero with
    'quasistatic' set to False
  - update comments in Problem.solve()

- merge pull request #692 from rc/fix-saving-animations, closes #266

  - update Viewer.encode_animation() for current ffmpeg (4.2.4)
  - postproc.py: update --ffmpeg-options default, update docstring
  - docs: update postproc.py help message in users guide

- merge pull request #694 from vlukes/new_example_application

  - new example application

- merge pull request #695 from rc/link-example-applications

  - docs: link example application in index.rst, remove applications section

- merge pull request #696 from rc/fix-numpy-indexing-warning

  - fix numpy future warning in TransformToPlane.tensor_plane_stress() - using
    a non-tuple sequence for multidimensional indexing

- merge pull request #697 from rc/revisit-facet-orientations

  - describe quad face orientation groups in facets.py docstring
  - add missing permutations in GeometryElement.get_conn_permutations()
  - permute common quad face in _gen_common_data(), new _permute_quad_face()

    - all keys from _quad_ori_groups in sfepy.discrete.fem.facets are tested
    - update Test.test_continuity(), .test_gradients()

  - stop after all possible orientations were tested in _gen_common_data()

    - new _get_possible_oris()

- merge pull request #698 from rc/bernstein-basis-new-set-dofs

  - new BernsteinTensorProductPolySpace

    - new .__init__(), ._define_nodes(), _eval_base()

  - check shape of quadrature coordinates in PolySpace.eval_base()
  - new local argument of VolumeField.get_econn(), create missing surface data
  - move global basis functions from H1NodalMixin to new GlobalNodalLikeBasis

    - move ._setup_facet_orientations(), ._setup_edge_dofs(),
      ._setup_face_dofs(), ._setup_facet_dofs(), ._setup_bubble_dofs()

  - new H1BernsteinVolumeField, H1BernsteinSurfaceField (WIP)

    - new H1BernsteinVolumeField.create_basis_context() - same "hack" as in
      H1HierarchicVolumeField

  - new local argument of IGField.get_econn()
  - split/fix IGField.set_dofs(), new Field.set_dofs()

    - split IGField.set_dofs() into new IGField.get_surface_basis() (IGField-
      specific) and project_to_facets() (general)
    - fixes wrong shape of projected function return value
    - Field.set_dofs() calls project_to_facets()

  - new GlobalNodalLikeBasis.get_surface_basis()

    - Field.set_dofs() works with H1BernsteinVolumeField

  - update H1NodalMixin.set_dofs(), H1HierarchicVolumeField.set_dofs()

    - according to Field.set_dofs()

  - update tests/test_poly_spaces.py for Bernstein basis

    - update _gen_common_data(), test_partition_of_unity(), test_continuity()

  - new define() in sinbc.py example, update docstring

- merge pull request #699 from vlukes/material_shape_update

  - update materials: allow for (1, N, M) material shape
  - update terms: replace FMF_SetCell() by FMF_SetCellX1() for material arrays
  - update piezo term
  - modify examples to use the newly allowed material shape
  - update users guide

- merge pull request #700 from vlukes/const_material

  - update material functions
  - update CMappings
  - new FMF_PtrCellX1() macro
  - new make_full_mat_array(): (1, nqp, n, m) -> (nel, nqp, n, m)
  - update terms
  - rename and move make_full_mat_array() -> Term.tile_mat()

- merge pull request #701 from rc/bernstein-basis-simplex-naive

  - new BernsteinSimplexPolySpace (naive implementation)

    - new .__init__(), ._define_nodes(), ._get_barycentric(), _eval_base()

  - clean up BernsteinTensorProductPolySpace._eval_base()
  - update tests/test_poly_spaces.py for simplex Bernstein basis

    - update _gen_common_data(), test_partition_of_unity()

- merge pull request #702 from rc/constant-materials-small-fixes

  - use single cell dummy material in DotProductVolumeTerm.get_fargs()
  - fix float comparison in PhysicalQPs.get_shape()

- merge pull request #703 from vlukes/fix_d_sd_div_grad

  - fix opt_material in SDDivGradTerm

- merge pull request #704 from vlukes/update_coefs_base

  - update CoefMN.set_variables_default()

- merge pull request #705 from rc/fix-more-c-funs-for-const-mats

  - fix d_of_ns_min_grad() for constant materials
  - fix dw_st_grad_div() for constant materials
  - fix actBfT(), sym2nonsym() for constant materials

.. _2020.4-2021.1:

from 2020.4 to 2021.1
=====================

- merge pull request #664 from vlukes/update_coefs_base

  - update CoefNN: allow different dimensions for row and column
  - rename CoefNN to CoefMN

- merge pull request #665 from vlukes/update_coefs_phono

  - fix volumes in DensityVolumeInfo class

- merge pull request #668 from rc/multilinear-terms, closes #666

- new sfepy/terms/terms_multilinear.py with ETermBase and multilinear terms
  (see git log)

  - implemented terms: de_volume_integrate, de_laplace, de_volume_dot,
    de_surface_dot, de_s_dot_mgrad_s, de_non_penetration_p, de_div_grad,
    de_convect, de_div, de_stokes, de_lin_elastic, de_cauchy_stress
  - docs: add multi-linear terms section to developer guide (WIP)
  - remove tests/test_input_stokes_slip_bc_penalty.py
  - script/gen_gallery.py: remove stokes_slip_bc_penalty.py custom
    visualization
  - script/gen_term_table.py: make new section for multi-linear terms - update
    typeset_term_tables()
  - docs: add multi-linear terms section to terms overview, explain de prefix
  - sync test_stokes_slip_bc() with stokes_slip_bc.py example

- merge pull request #675 from rc/fix-term-table

  - script/gen_term_table.py: fix typeset_term_tables() (de_* terms twice)

- merge pull request #671 from rc/propagate-ebcs-via-epbcs

  - propagate ebcs via epbcs in EquationMap.map_equations()
  - simplify some expressions in EquationMap.map_equations()

- merge pull request #676 from antonykamp/374-doc-link-terms-example - closes
  #374

  - get term-example dict, generate links
  - Added omits, edited reference style
  - shorten long term argument names to allow line breaks in term table

.. _2020.3-2020.4:

from 2020.3 to 2020.4
=====================

- merge pull request #648 from rc/fix-mat-init-msg

  - fix error message in Material.__init__()

- merge pull request #650 from vlukes/update_homog_app

  - add 'define_args' argument to get_homog_coefs_nonlinear()
  - optional multiplier in updating_coors structure
  - new 'id' identifier in micro_states

- merge pull request #652 from rc/fix-inline-defs, closes #651

  - mark private inline functions in mesh.c with static keyword - mark
    _det3x3(), _tri_area(), _aux_hex()

- merge pull request #654 from heczis/add_ogden_term

  - Add Ogden term

- merge pull request #656 from rc/resview-warp-scalar

  - resview.py: support warp by scalar in pv_plot()

- merge pull request #658 from rc/plot-cmesh-label-local-1d

  - improve 1D mesh entities label placement in label_local_entities()

- merge pull request #660 from vlukes/fix_splinebox

  - fix splinebox - Python3 related issue

- merge pull request #661 from vlukes/fix_homog_app

  - fix HomogenizationApp: Python3 related issue

- merge pull request #662 from rc/serendipity-basis

  - new SerendipityTensorProductPolySpace
  - new H1SNodalVolumeField, H1SNodalSurfaceField

    - new H1SNodalVolumeField.create_basis_context() - same "hack" as in
      H1HierarchicVolumeField

  - update test_continuity() for serendipity basis
  - new test_partition_of_unity() in tests/test_poly_spaces.py - tests
    lagrange, serendipity bases
  - new script/gen_serendipity_basis.py
  - new sfepy/discrete/fem/_serendipity.py, generated by
    script/gen_serendipity_basis.py
  - clean up sfepy/discrete/fem/poly_spaces.py
  - limit 3_8 serendipity basis order to 2 in _gen_common_data()

    - update docstring of tests/test_poly_spaces.py, remove obsolete
      information

.. _2020.2-2020.3:

from 2020.2 to 2020.3
=====================

- merge pull request #625 from vlukes/mesh_write_binary

  - update meshio: pass options in kwargs to meshio library
  - update mesh conversion script: write mesh in ascii format

- merge pull request #626 from rc/fix-elastodynamics-docstring, closes #598

  - fix equation in docstring of elastodynamic.py example

- merge pull request #582 from vlukes/pyvista

  - new script for visualisations using pyvista
  - update requirements: pyvista
  - update users_guide
  - add pyvista to instalation requirements
  - update names in installation docs

- merge pull request #634 from vlukes/fix_web

  - various fixes in web doc

- merge pull request #635 from vlukes/mumps5.2

  - update mumps_struc_c for v5.2
  - absolute import

- merge pull request #637 from vlukes/example_bd2b

  - new example application

- merge pull request #638 from vlukes/update_web

  - links to sfepy_examples: http --> https

- merge pull request #639 from heczis/add_gen_yeoh_term

  - Add generalized Yeoh hyperelastic term (total Lagrangian formulation)
  - Add example for generalized Yeoh term

- merge pull request #640 from vlukes/meshio_file_format_variants

  - fix meshio interface to accept ascii/binary file format variants
  - check meshio version, 'binary' argument in write() function from v4.0.3

- merge pull request #641 from vlukes/fix_meshio

  - fix #640

- merge pull request #643 from rc/get-virtual-dof-conn

  - update fieldvariable.get_dof_conn() for no primary variable

.. _2020.1-2020.2:

from 2020.1 to 2020.2
=====================

- merge pull request #589 from heczis/fix_meshio_import

  - Fix import from meshio

- merge pull request #590 from rc/fix-meshio-optional-formats

  - fix update_supported_formats() to skip optional meshio formats
  - fix test_read_meshes() to initialize self.meshes and catch read failures
  - docs: mention installation of meshio dependencies

- merge pull request #588 from vlukes/new_web_look

  - update doc web pages: use sphinx-rtd-theme
  - drop `download.php`, link to github releases
  - use new sfepy logo
  - web: update gen_gallery.py, use os.path.splitext() to remove suffix

- merge branch 'fix-xyz-mesh-io'

  - fix xyzmeshio.read() for single entity meshes

- merge pull request #591 from vlukes/update_web

  - web: fix css file
  - new example application

- merge pull request #592 from bubulk/docker-doc

  - Add simple mention of docker images to installation doc.

- merge pull request #594 from zitkat/dg-sfepy-integration

  - Allow term arguments to be callables.
  - Add refine for 1D meshes to FEDomain
  - Ensure propagation of verbose parameter to: * materials update * equation
    time_update
  - Override C implementation of mesh_get_facet_normals for 1D mesh.
  - Rework propagating verbose flag in problem.py: * use self.conf.get instead
    of get_default_attr * use default values of called methods
  - Add cell groups treatment to refine for 1D mesh.

- merge pull request #599 from rc/docs-problem-description-functions

  - docs: add problem description functions examples to users guide

- merge pull request #596 from zitkat/dg-sfepy-bc-integration

  - Add mechanism for treating DG BC conditions:

    - new syntax dgebc and dgepbc for problem description file
    - new classes and mechanism for parsing BCs including  gradient
    - specific evaluation in map_equations method

- merge pull request #600 from vlukes/fix_vtk_1d

  - fix saving 1D meshes in vtk and vtu formats

- merge pull request #595 from rc/gc-hom-engine

  - force garbage collection (free memory) in HomogenizationWorker.__call__()

- merge pull request #602 from vlukes/fix_test_region

  - fix test_regions: test only existing cell and vertex groups

- merge pull request #603 from vlukes/fix_meshio_read_data, closes #601

  - fix meshio read_data()

- merge pull request #604 from rc/gc-run-tests

  - force garbage collection (free memory) after each test file in run_tests()

- merge pull request #611 from zitkat/i610, closes #610

  - Force nm.int64 in nm.prod in Equations.create_matrix_graph to prevent
    overflow on win64 platform.

- merge pull request #606 from rc/simplify-gitignore

  - simplify .gitignore - add new non-python files with 'git add -f'

- merge pull request #609 from zitkat/mumps-libname, closes #608

  - Add fallback to lib*.dll name for loading mumps on win32.

- merge pull request #613 from vlukes/fix_coefs

  - fix saving coefficients

- merge pull request #614 from vlukes/fix_meshio_cell_types, closes #607

  - fix meshio cell_types
  - update documentation of Preprocessing

- merge pull request #593 from zitkat/gmshio-thesecond

  - Add custom GmshIO class to enable reading and writing of ElementNodeData: *
    includes reading multiple time steps * omits "gmsh:ref" data to make
    working with gmsh comfortable
  - new file_format application option
  - add gmsh-dg variant description to _supported_formats
  - use standard vertex and cell group names in MeshioLibIO.read_data()

- merge pull request #616 from rc/fix-mesh-generator-scripts

  - move suffix check from Problem.setup_output() to new check_format_suffix()
  - script/blockgen.py: check format and suffix, update, fix for current MeshIO
  - script/cylindergen.py: check format and suffix, update, fix for current
    MeshIO

- merge pull request #618 from vlukes/update_ceofs_phonoic

  - update coefs_phononic: allow for more general usage

- merge pull request #620 from rc/misc-fixes

  - fix plot_edges() for 1D
  - move time stamping in Timer.start(), .stop() closer to measured code

- merge pull request #621 from rc/guard-memory-fem-mapping

  - new raise_if_too_large()
  - guard memory use in VolumeMapping.get_mapping()
  - CI: add psutil to test environments
  - new Config.refmap_memory_factor(), update site_cfg_template.py
  - new PSUTIL_MIN_VERSION in sfepy/version.py, add psutil as optional to
    setup.py
  - docs: document optional psutil requirement

- merge pull request #619 from rc/speed-up-set-dofs

  - major speed up by not calling hstack() on Field.get_dofs_in_region()
    results - update EquationMap.map_equations(), H1NodalMixin.set_dofs(),
    IGField.set_dofs(), FieldVariable._get_setter()

- merge pull request #612 from zitkat/dg-method-main

  - discontinuous Galerkin method implementation and examples (see git log)

.. _2019.4-2020.1:

from 2019.4 to 2020.1
=====================

- merge branch 'dispersion-misc-updates'

  - dispersion_analysis.py: allow custom save_eigenvectors()
  - update read_log() for logs prior to int-log-labels branch merge
  - dispersion_analysis.py: increase log precision to .12e

- merge pull request #571 from vlukes/fix_trace_normals

  - fix `dw_surface_dot` term: nvec orient. - outwards to the parent region of
    the virtual variable

- merge pull request #575 from rc/fix-splinebox-python-3.7

  - fix argument type in SplineRegion2D.create_spb() - fixes Travis CI failure
    for Python 3.7, numpy-1.18.1

- merge pull request #574 from rc/docs-update-space-def, closes #572

  - docs: update function space definition

- merge branch 'fix-comparison'

  - fix comparison operator in NurbsPatch.elevate()

- merge pull request #577 from bubulk/multiprocessing-spawn

  - Update dirty fix for multiprocessing with default 'spawn' method.

- merge branch 'fix-read-log'

  - fix read_log() for empty log group

- merge pull request #579 from rc/drop-python-2.7-support

  - CI: remove Python 2.7
  - docs: update installation instructions (drop Python 2.7 support)
  - remove PysparseEigenvalueSolver (Pytnon 2.7 only), update docs, tests,
    setup.py

- merge branch 'fix-slepc-eigenvectors'

  - fix eigenvectors returned by eig.slepc - fix
    SLEPcEigenvalueSolver.__call__()

- merge branch 'dispersion-dict-pars'

  - add allow_tuple, free_word arguments to dict_from_string()
  - new apply_units_to_pars()
  - dispersion_analysis.py: replace apply_units() by apply_units_to_pars()
    - new pars_kinds dict
    - update define() for parameters in Struct
  - dispersion_analysis.py: reorder definitions

- merge pull request #581 from rc/homogenized-coefs-dtype

  - force dtype in {CoefNN, CoefN, CoefOne}.set_variables_default()
    - fixes numpy casting TypeError when adding inplace a complex128 array to
    a float64 one

- merge pull request #580 from vlukes/meshio, closes #460

  - use meshio to read/write mesh files
  - meshio: update examples
  - meshio: update tests
  - fix numpy.savetxt() issue, see https://github.com/numpy/numpy/issues/10018
  - fix meshes to pass meshio tests
  - fix comsol writer: flatten mat_ids arrays
  - update requirements: meshio
  - update meshio: new handling with file formats
  - put ansys_cdb back
  - update script for mesh conversion and add checks to any_from_filename()
  - new xdmf3 extension to h5 format
  - CI: install h5py
  - CI: install netCDF4
  - CI: install meshio
  - update installation instructions: add meshio to requirements
  - fix vtk_cell_types in postprocessing code
  - update create_file_source()
  - new MeshioLibIO.read_data(), ._get_dimension()
  - change format to double in all medit meshes, fixes precision issues
  - skip meshio write-failed formats in test_write_read_meshes()
  - hdf5 + xdmf -> hdf5-xdmf format with .h5x and .xdmf suffixes

- merge branch 'misc-fixes'

  - fix classifiers in setup.py for Python 3
  - fix count type in get_log_freqs()

- merge pull request #583 from rc/fix-ones-dim-saving

  - fix saving correctors in OnesDim

- merge pull request #586 from vlukes/example_poropiezo

  - new poropiezo example

- merge pull request #587 from rc/web-update-front-pages

  - web: add PUCGen link
  - web: update support
  - web: add link to sfepy docker image

.. _2019.3-2019.4:

from 2019.3 to 2019.4
=====================

- merge pull request #562 from vlukes/fix_mumps, closes #561

  - fix MUMPSSolver.__init__(): raise if no mpi4py found
  - check material arguments in Term.classify_args()
  - use conf.funmod for default user parameters of terms in
    Problem.set_equations()
  - use directly conf in Problem.set_equations() to support define()

- merge branch '1d-surface-mapping'

  - add 0_1 geometry to geometry_data, new _get_grid_0_1()
  - skip 0_1 geometry in CMesh.set_local_entities()
  - create surface facet also for 1D meshes in FEDomain.__init__()
  - fix n_facet for 1D cells in Region.update_shape()
  - update Lagrange polynomial spaces for 0_1 geometry
  - add 0_1 geometry to quadrature_tables, skip it in test_quadratures()
  - update SurfaceMapping.get_base(), .get_mapping(), _s_describe() for 0_1
    geometry

- merge pull request #559 from zitkat/misc-updates

  - Add Functionize decorator Functionize decorator converts python function to
    sfepy.discrete.functions.Function obejct
  - Refactor - more descriptive name for Functionize - make_sfepy_function
  - Add documentation.

- merge pull request #560 from zitkat/gmsh-write

  - Add write mesh method
  - method for writing element node data.
  - Add correct interpolation scheme output.
  - Add msh output to problem ouput_modes.
  - Update gmsh data reading
  - Better dimensions treatment when loading meshes from gmsh format
  - Format and extend documentation.

- merge branch 'xyz-meshio'

  - new XYZMeshIO
  - new test meshes in xyz format (2_4, 3_4 cells)
  - test XYZMeshIO, skip in test_write_read_meshes()

- merge pull request #567 from vlukes/update_homog

  - update homogenization: avoid reloading coefficients and correctors
  - update CorrEval class: save corrector values
  - update CoefOne: allow for "Corr1 + Corr2" as in CoefNN
  - update get_mesh_coors(): add `actual` parameter which goes into
    `domain.get_mesh_coors()`
  - new updating procedure in nonlinear homogenization
  - update nonlinear homog. example
  - fix CorrEval, CoefSum, CoefEval for multiprocessing evaluation
  - update CoefOne: make `set_variables` optional
  - update 'nls_iter_hook': simplify hook management
  - uppdate get_homog_coefs_nonlinear(): store the actual time step to
    `macro_data`
  - update recover_micro_hook_eps(): indicate the total number of recovered
    microstructures
  - fix CoefOne and CorrMiniApp classes
  - update DeformationGradientTerm: allow to choose actual or undeformed
    reference configuration

- merge branch 'matlab-evp-solver'

  - new MatlabEigenvalueSolver
  - new matlab_eig() in sfepy/solvers/matlab_eig.m
  - update tests/test_eigenvalue_solvers.py for eig.matlab

- merge branch 'matlab-evp-solver-2'

  - add method parameter to MatlabEigenvalueSolver
  - update matlab_eig() for method parameter, simplify logic

- merge branch 'misc-fixes'

  - fix Term.check_shapes() for spaces in shape specifications
  - fix geme_invert3x3() for very small cells - fixes zero basis gradient
    (singular matrix) problem in reference mappings for meshes with very
    small cells (about 1e-6 edge size in 3D and smaller)

- merge branch 'misc-fixes-2'

  - fix geme_invert4x4() for very small cells, do not throw error
  - fix dispersion_analysis.py for single requested eigenvalue
  - fix collect_term(), create_bnf() for leading minus complex term
    coefficients - example: -1-1j was parsed as -(1-1j) => -1+1j
  - fix FEField.get_base() for subdomains and basis transform (iels not in key)
  - do not save results twice in laplace_refine_interactive.py example
  - add '.' to sys.path in script/show_mesh_info.py

- merge branch 'misc-fixes-3'

  - fix FEField.create_mapping() for basis transform and subdomains - WIP -
    raise exception for surface integration
  - docs: describe active_only option in users guide

- merge pull request #569 from rc/verify-tractions

  - new verify_tractions() post-process hook in linear_elastic_tractions.py
    example

- merge branch 'speed-up-log-plotter'

  - update Log, LogPlotter to send/receive last values only
  - apply plotting commands in LogPlotter just before redrawing canvas
  - clear axes in LogPlotter.apply_commands()
  - pass plot_kwargs to LogPlotter.__call__() and in add_axis command
  - update log parameters in live_plot.py example

- merge branch 'int-log-labels'

  - use int labels in Log, simplify LogPlotter by passing keys with commands
  - test multiple groups and lines in test_log_create(), fix test_log_rw()

- merge branch 'fix-qeigen'

  - fix Timer import

- merge branch 'misc-fixes-4'

  - fix barycentric array shape in eval_basis_lagrange() for 0_1 geometry
  - include matlab_eig.m in source tarball
  - script/gen_release_notes.py: fix for Python 3, improve formatting
  - update Msh2MeshIO docstrings, fixes PDF documentation build
  - fix offset in plot_log() for given groups

.. _2019.2-2019.3:

from 2019.2 to 2019.3
=====================

- merge branch 'dispersion-brillouin-stepper'

  - dispersion_analysis.py: update get_stepper() - new BrillouinStepper
  - dispersion_analysis.py: allow passing wdir as argument to
    assemble_matrices()
  - dispersion_analysis.py: new --stepper option - update
    process_evp_results(), main()

- merge branch 'python3-fixes'

  - fix enum() for Python 3
  - fix load_library() for Python 3 - copy dec() from ioutils.py to avoid sfepy
    modules dependence
  - remove btrace_python
  - get python version from site_cfg.py in Makefile

- merge branch 'test-install-fixes'

  - add '.' to sys.path in interactive examples
  - test_install.py: update for changing numpy output - add match_numbers
    argument to report() - clean up
  - fix Python 3 string encoding problem in save_raw_bg_logs()
  - fix band_gaps_rigid.py example for Python 3

- merge pull request #545 from vlukes/fix_schur_mumps, closes #544

  - fix ls.schurs_mumps: 'bloc' -> 'block', use active dof info

- merge branch 'parallel-timing-stats'

  - new sfepy/base/timing.py, new Timer
  - store current dt in Timer
  - use time.perf_counter() in Timer in Python 3
  - use Timer in poisson_parallel_interactive.py example to gather timing stats
  - new call_in_rank_order(), update view_petsc_local()
  - update stats items in poisson_parallel_interactive.py example
  - poisson_parallel_interactive.py: new --new-stats option - use
    call_in_rank_order() to save stats - new save_stats()
  - poisson_parallel_interactive.py: new --stats option
  - poisson_parallel_interactive.py: fix save_stats() for Python 3
  - fix partition_mesh() for Python 3
  - biot_parallel_interactive.py: new --stats, --new-stats options - use Timer,
    update analogously to poisson_parallel_interactive.py

- merge branch 'fix-for-mayavi-4.7.1'

  - fix get_data_ranges()

- merge branch 'dispersion-define-kwargs'

  - dispersion_analysis.py: new --define-kwargs option, update
    assemble_matrices()
  - dispersion_analysis.py: new mesh_eps argument in define(), update docstring

- merge branch 'slepc-evp-solver'

  - update EigenvalueSolver.__init__() to pass on additional arguments
  - new SLEPcEigenvalueSolver, init_slepc_args()
  - update tests/test_eigenvalue_solvers.py for eig.slepc
  - update installation and development docs
  - check slepc4py version in setup.py, update sfepy/version.py
  - allow failing of eig.slepc in tests/test_eigenvalue_solvers.py (optional
    dep.)

- merge branch 'fix-for-numpydoc-0.9.1'

  - fix See Also sections for numpydoc-0.9.1

- merge branch 'python3-metaclass-update'

  - remove Python 2 metaclass attributes from Solver subclasses
  - add SolverMeta metaclass to Solver in Python 3 compatible way - add to base
    class only

- merge pull request #549 from rc/complex-mat-pars, closes #547

  - update ConstantFunction, ConstantFunctionByRegion for complex dtype
  - use previously set variables and materials in Problem.set_equations() -
    initialize default .conf_variables, .conf_materials in .__init__()
  - fix Region.get_cell_indices() for no common cells
  - update tests/test_functions.py to test complex material parameters - update
    test_material_functions()

- merge pull request #554 from vlukes/update_convert_mesh, closes #553

  - update convert_mesh.py: new '3d' option - write only cells of dimension 3

- merge pull request #555 from rc/ci-fix-igakit-download

  - CI: fix igakit download

- merge pull request #550 from rc/fix-site-cfg-template-python-version

  - fix default python_version in site_cfg_template.py to work with Makefile -
    update Config.python_version()
  - update Makefile to inform about site_cfg.py Python version setting
  - report Python version in setup.py

- merge pull request #556 from vlukes/update_convert_mesh

  - update convert_mesh.py: new '--cell-dim' option, write only cells of a
    given dim

- merge pull request #557 from rc/use-timer, closes #548

  - allow starting Timer on creation in .__init__(), improve .stop() message
  - use Timer instead of time.clock(), clean up
  - remove unused mark_time()

- merge pull request #558 from rc/fix-set-bcs-corr, closes #551

  - do not set problem variables in CorrSetBCS

.. _2019.1-2019.2:

from 2019.1 to 2019.2
=====================

- merge branch 'allow-zero-eigs'

  - new EigenvalueSolver._ret_zero(), use it in standard_call()
  - update test_eigenvalue_solvers() to test zero eigs

- merge branch 'plot-log-nbins'

  - allow specifying numbers of bins in x, y axes in plot_log()
  - script/plot_logs.py: new --nbins option

- merge branch 'misc-updates'

  - fix output_array_stats() for empty arrays
  - add velocity, acceleration, dyn_viscosity to apply_unit_multipliers()
  - allow corrector saving in ShapeDim, OnesDim calls
  - dispersion_analysis.py: attach options to Problem, clean up
  - dispersion_analysis.py: allow custom steppers, new get_stepper()
  - fix/update OnesDim, VolumeFractions for complex variables

- merge pull request #504 from bubulk/refine-hanging-segfault

  - fix mei_next()

- merge pull request #506 from vlukes/homog_ts

  - fix ls.py: avoid conflict between presolve option and presolve() method
  - update coefs: remove obsolote FM coefs, e.g. CoefFMSymSym can be repalced
    by TSCoef(CoefSymSym)
  - update set_conf_solvers(): allow to switch off a solver (e.g. ts in
    homogenization)
  - fix TSTime coefficient

- merge pull request #510 from vlukes/fix_mat_opt, closes #509

  - fix material_opt.py: Python 2/3 compatibility issue - zip()

- merge pull request #511 from vlukes/fix_mat_opt

  - fix material_opt.py: Python 2/3 compatibility issue - part II.

- merge pull request #512 from rc/fallback-solvers-fixes

  - fix use_first_available() to catch and report all relevant exceptions
  - clean up AutoFallbackSolver.__new__()
  - fix use_first_available() for Python 3

- merge pull request #513 from vlukes/update_homog

  - update homogenization: improved saving/dumping of time variable correctors
  - remove verify_correctors() in TCorrectorsViaPressureEVP class
  - fix recover_micro_hook_eps()
  - update saving homog. correctors, remove 'dump_*' options
  - update homog. examples
  - reimplement matching periodic planes - allow skew geometries
  - fix match_plane_by_dir() for 2D meshes

- merge pull request #515 from vlukes/fix_mumps

  - fix mumps solver for "small and dense" matrices

- merge pull request #516 from vlukes/nonsym_prestress_term

  - update LinearPrestressTerm term: allow for non-symmetric form

- merge pull request #517 from vlukes/fix_meshio

  - fix VTKMeshIO.write(): incorrect tensor data shape

- merge pull request #519 from rc/docs-update-citing-support

  - docs: update sfepy citing section
  - docs: update support section, move past support to a separate page
  - docs: add article full text link

- merge pull request #520 from vlukes/fix_mumps_solver

  - fix SchurMumps and MUMPSParallelSolver: add `memory_relaxation` parameter

- merge pull request #521 from vlukes/fix_splinebox

  - fix SplineBox: remove forgotten prints

- merge pull request #523 from heczis/fix_scipy_misc

  - Remove use of scipy.misc

- merge pull request #525 from vlukes/update_splinebox

  - update SplineBox.write_control_net(): use control point values as
    coordinates

- merge pull request #526 from vlukes/find_map

  - reimplement find_map(): use the scipy.spatial.cKDtree implementation

- merge pull request #527 from vlukes/remote_mirror

  - new remote mirror regions: allow mirror regions with no common vertices

- merge pull request #528 from vlukes/mirror_region_misc

  - update User's Guide: new `mesh_eps` option
  - fix get_conn_key() for mirror integration
  - update User's Guide: new `mirror region` option in the regions definition

- merge pull request #529 from rc/trace-regions

  - update create_bnf(), create_arg_parser() for trace regions
  - update test_parse_equations()
  - update Term.setup_formal_args() for trace regions

- merge pull request #530 from vlukes/mirror_region_misc

  - update region.setup_mirror_region() for `tr(reg, var)` syntax

- merge pull request #531 from vlukes/multi_traces

  - allow for multiple traces

- merge pull request #532 from vlukes/fix_traces

  - fix traces for evaluation mode

- merge branch 'misc-fixes'

  - fix compute_nodal_normals() for higher order nodes
  - fix PDESolverApp.call() to obey save_results option

- merge pull request #533 from bubulk/ci-update

  - Updated CI .yml config files.
  - Temporary disable conda update for Windows.

- merge pull request #534 from rc/fix-for-python-3.7

  - fix Material.iter_terms(), get_data_name() for Python 3.7 (PEP 479)

.. _2018.4-2019.1:

from 2018.4 to 2019.1
=====================

- merge branch 'change-complex-output-names'

  - update convert_complex_output(): convert name into real.name, imag.name

- merge branch 'evp-gallery'

  - update EVPSolverApp.call() to return (problem, evp) tuple
  - script/gen_gallery.py: make is_scalar_bar configurable in custom views
  - script/gen_gallery.py: update generate_images() for eigenvalue problems
  - script/gen_gallery.py: add custom views for quantum examples
  - add single line descriptions to quantum examples

- merge branch 'plot-log-swap-axes'

  - add swap_axes argument to plot_log()
  - script/plot_logs.py: new --swap-axes option
  - new draw_data(), update LogPlotter.process_command()
  - update read_log(), plot_log() for complex data - use draw_data()
  - improve log header saving with empty axes labels in Log.__init__()

- merge branch 'dispersion-kappa-mode'

  - dispersion_analysis.py: fix kappa mode - mask non-physical eigenmodes - fix
    saving eigenvectors
  - dispersion_analysis.py: support logging standard waves in kappa mode - log
    plot: use viridis colormap, plot points only
  - dispersion_analysis.py: update logging in omega mode according to kappa
    mode
  - dispersion_analysis.py: allow empty regions
  - dispersion_analysis.py: output eigenvalue indices
  - dispersion_analysis.py: update docstring, change default range
  - new square_1m.mesh, square_2m.mesh for kappa mode

- merge branch 'update-import-file'

  - explicitly name path to remove in import_file()

- merge pull request #493 from vlukes/vtk_point_probe

  - new VTK point probe

- merge pull request #495 from bubulk/tests-fix

  - Misc ASV data changes.
  - Removed ASV data.
  - Add tests for missing igalib module to tests suite.
  - "Fixed" NamedTemporyryFile behaviour on windows platform.
  - Fix py27/win64/numpy numpy.array.shape[] type inconsistency (int-long).
  - Add fix for broken multipocessing implementation on windows platform.
  - Updated appveryor.yml
  - Fixed appveryor.yml #2
  - Updated appveryor.yml (#3): trying to fix IGAkit build.

- merge pull request #482 from {vlukes,rc}/update_linear_solvers

  - update MUMPS solver: consider matrix symmetry, verbose flag, ...
  - update linear direct solvers: new scipy_superlu, scipy_umfpack
  - new 'fallback' option to linear solvers, new 'auto_direct' solver
  - update examples: use fallback option and AutoDirect/AutoIterative solvers
  - update doc: new "virtual" solver and fallback option
  - fix ls.schur_mumps: convert matrix to COO format
  - rename ls_fallback() -> use_first_available(), improve failure reporting

- merge pull request #499 from rc/block-compat

  - add block() from NumPy 1.14.1 to compat.py
  - use block() in _build_cauchy_strain_op()

- merge pull request #500 from rc/surface-integrate-mat-term

  - new IntegrateSurfaceMatTerm (ev_surface_integrate_mat)
  - fix docstring of IntegrateMatTerm
  - new alias of ev_integrate_mat: IntegrateVolumeMatTerm
    (ev_volume_integrate_mat)
  - rename IntegrateMatTerm -> IntegrateVolumeMatTerm (remove the alias) -
    ev_integrate_mat -> ev_volume_integrate_mat
  - update for ev_volume_integrate_mat
  - docs: update users guide for ev_volume_integrate_mat

- merge pull request #501 from vlukes/update_micmac

  - update get_correctors_from_file() - allow corrector file names to be
    different from corrector names

- merge pull request #502 from bubulk/fix-c-defs

  - Removed "obsolete" (and trouble-making) __SDIR__ C pre-processor cmd line
    defs.

- merge branch 'update-show-authors'

  - fix script/show_authors.py for unicode, update for several names

- merge branch 'write-log'

  - new write_log()
  - move header writing from write_log() to new _write_header() - update
    Log.__init__()
  - new tests/test_log.py - new test_log_create(), test_log_rw()

- merge branch 'fix-log-plot-kwargs'

  - fix saving plot kwargs in Log.__init__()
  - update test_log_create()

- merge branch 'quadratic-evp'

  - new sfepy/solvers/qeigen.py: QuadraticEVPSolver, LQuadraticEVPSolver - new
    standard_call()
  - update solver_table
  - update script/gen_solver_table.py
  - docs: update users guide
  - docs: sync module index of developer guide with current sources
  - dispersion_analysis.py: update for LQuadraticEVPSolver
  - dispersion_analysis.py: move _max_diff_csr() to sfepy/linalg/utils.py
  - add debug parameter to LQuadraticEVPSolver
  - move QuadraticEVPSolver into sfepy/solvers/solvers.py

- merge branch 'dispersion-split'

  - dispersion_analysis.py: move assembling to new assemble_matrices()
  - dispersion_analysis.py: put coefficients other than omega, kappa into
    equations
  - dispersion_analysis.py: remove _le suffix from problem-dependent functions
  - dispersion_analysis.py: update symmetry checks in assemble_matrices()

- merge branch 'quadratic-evp-status'

  - update LQuadraticEVPSolver to store matrix info in status in debug mode
  - update LQuadraticEVPSolver to store solution errors in status in debug mode

- merge branch 'dispersion-generalize'

  - dispersion_analysis.py: use new build_evp_matrices() - assemble_matrices()
    returns dict of blocks
  - dispersion_analysis.py: use new process_eigs()
  - dispersion_analysis.py: set output_dir to problem in assemble_matrices()
  - dispersion_analysis.py: set wave vector direction using material function -
    new get_wdir() - update define(), set_wave_dir(), assemble_matrices()
  - dispersion_analysis.py: use new setup_n_eigs()
  - dispersion_analysis.py: update, rename process_eigs() ->
    process_evp_results()
  - fix Term.get_str() for complex coefficients
  - do not check argument count when saving figures in Log.__call__()

- merge branch 'lin-convect2-term'

  - new LinearConvect2Term (dw_lin_convect2)

- merge pull request #503 from vlukes/update_homog

  - update homog. app - call setup_options() properly
  - update CorrEval

.. _2018.3-2018.4:

from 2018.3 to 2018.4
=====================

- merge branch 'fix-viewer-for-mayavi-4.6.2', closes #466

  - fix get_data_ranges() for Mayavi 4.6.2
  - fix hanging of ClosingHandler.object_button_quit_changed()

- merge pull request #468 from vlukes/mumps_update

  - fix mumps solver and update for mumps library version 5.1.2

- merge pull request #467 from vlukes/update_homog_example

  - update linear_homogenization.py: set option 'is_linear' to True for
    speed-up

- merge pull request #470 from bubulk/ci-appveyor, closes #469

  - Changed 'version' statement.

- merge pull request #471 from vlukes/coefs_update

  - update coeffficients.py: saving of complex valued coefficients to .txt

- merge pull request #472 from vlukes/fix_splinebox

  - fix scipy FutureWarning warning in scipy.linalg.lstsq()

- merge pull request #473 from vlukes/fix_homperf_example

  - fix the url of the cited paper in
    examples/homogenization/perfusion_micro.py

- merge pull request #474 from vlukes/new_schur

  - update mumps solver - rename functions: set_b -> set_rhs, set_A_* ->
    set_mtx_*
  - new schur_mumps linear solver
  - remove old Schur linear solver
  - update mumps wrapper - fix various library versions

- merge pull request #475 from vlukes/update_mumps_sym

  - update mumps: use a symmetric solver if the matrix is symmetric

- merge pull request #480 from vlukes/fix_gallery_link

  - fix links to the gallery

- merge pull request #481 from vlukes/fix_download_page

  - fix downloads counter after HW upgrade

- merge branch 'fix-refmaps-small-cells'

  - fix reference mappings in C for cells with very small volumes - update
    _v_describe(), dq_finite_strain(), dq_tl_finite_strain_surface()

- merge branch 'dispersion-post-processing'

  - dispersion_analysis.py: support eigenvector post-processing - update
    save_eigenvectors()
  - dispersion_analysis.py: new --post-process option - new
    compute_von_mises(), update define_le()
  - support complex arrays in output_array_stats()
  - dispersion_analysis.py: output nonzero stats of matrices
  - dispersion_analysis.py: new --save-regions option

- merge branch 'example-boundary-fluxes'

  - new meshes/2d/cross-51-0.34.mesh
  - update poisson_neumann.py example to compute boundary fluxes, change mesh

- merge branch 'example-material-opt-clean-up'

  - remove matplotlib backend selection (Log does not need it anymore)
  - add a rudimentary docstring to material_opt.py example

- merge pull request #484 from vlukes/fix_mumps

  - fix ls.mumps for symmetric complex valued matrices

- merge branch 'log-complex-values'

  - support str.format() style string formatting in Log - enables logging
    complex values - update .__init__(), .add_group(), .__call__()
  - fix default plot_kwargs in Log.add_group()
  - fix Log.add_group() use in live_plot.py example
  - use str.format() style format in live_plot.py example
  - log and plot complex values in live_plot.py example
  - automatically plot and label real and imaginary parts in LogPlotter - the
    imaginary part line has alpha reduced to one half of the real part - update
    Log.plot_data(), LogPlotter.process_command()

- merge pull request #485 from lokik/parse_conf

  - Allow empty dict in dict_from_string
  - parse_conf grammar improvments
  - Code cleaning

- merge branch 'evp-solver-app', closes #479

  - new EVPSolverApp, based on SchroedingerApp
  - simple.py: update for EVPSolverApp
  - tweak EVPSolverApp output
  - update examples/quantum/ for EVPSolverApp - new get_exact() in examples -
    new report_eigs() in quantum_common.py
  - add usage examples to quantum_common.py example
  - remove schroedinger.py - quantum examples can be run using simple.py
  - remove sfepy/physics/ - SchroedingerApp obsoleted by EVPSolverApp
  - update setup.py files
  - update sfepy-run
  - update test_install.py
  - docs: update for no schroedinger.py and sfepy/physics/
  - preserve dtype in EVPSolverApp.make_full()
  - docs: sync module index of developer guide with current sources
  - add rudimentary docstrings to quantum examples
  - use short syntax in quantum_common.py, replace eig.pysparse with eig.scipy
  - docs: update Pysparse description
  - update solve_pde() for EVPSolverApp, support status arg in
    EVPSolverApp.call()
  - new tests/test_input_{boron,hydrogen,oscillator,well}.py

- merge branch 'fix-extend-cell-data'

  - fix extend_cell_data() for complex values

- merge pull request #488 from rc/mumps-parallel-check-main-scope

  - do not execute code in ls_mumps_parallel.py on import - add main scope
    check

- merge branch 'update-plot-mesh'

  - update text alignment in label_global_entities(), label_local_entities()
  - script/plot_mesh.py: new --{vertex,edge,face,cell,wireframe}-opts options
  - script/plot_mesh.py: new --no-axes, --no-show, figname options

- merge pull request #489 from rc/http-to-https

  - docs: http -> https in sfepy-related links

- merge pull request #490 from heczis/update_docstring

  - Update example in the docstring of Problem.create_evaluable

- merge pull request #492 from vlukes/fix_mumps_presolve

  - fix presolve() in  mumps solver

- merge branch 'allow-n-eigs-none'

  - allow unspecified n_eigs in EVPSolverApp.solve_eigen_problem()
  - update saving complex eigs in EVPSolverApp.save_results() for easier
    loading
  - split dense and sparse solvers in ScipyEigenvalueSolver - 'method' can be
    {'eig', 'eigh', 'eigs', 'eigsh'}, default is 'eigs' - basic support of
    n_eigs argument for dense problems
  - update tests and examples for current ScipyEigenvalueSolver

.. _2018.2-2018.3:

from 2018.2 to 2018.3
=====================

- merge branch 'update-contacts'

  - update ContactTerm.get_fargs() to use field approximation for basis
    gradient
  - new ContactInfo.__init__(), .update(), ContactTerm.get_contact_info(),
    update ContactTerm.get_fargs() to use ContactInfo
  - update ContactTerm.get_fargs() to return gap function in evaluation mode

    - update .call_function(), .eval_real()
    - rename .function() -> .function_weak()
    - new .integrate(), .function(), .get_eval_shape()

  - fix surface region case in extend_cell_data()
  - two_bodies_contact.py: save gap function in new post_process()
  - fix actual order of tensor product rules in QuadraturePoints
  - plot_quadrature.py: print info in _get_bqp()

- merge branch 'misc-updates'

  - new structify()
  - new triangulate()
  - script/convert_mesh.py: use triangulate()

- merge branch 'update-term-implementation-docs'

  - docs: remove NewTerm references from developer guide
  - docs: update text
  - docs: describe dw_s_dot_mgrad_s term implementation in developer guide

- merge branch 'fix-vtk-write'

  - fix tensor reshaping in VTKMeshIO.write() for single cell meshes

- merge branch 'set-variable-from-function'

  - new FieldVariable.set_from_function()
  - move FieldVariable.set_from_mesh_vertices()
  - rename FieldVariable.set_data_from_qp() -> .set_from_qp()
  - rename Variables.set_data_from_state() -> .set_from_state()
  - update for renames
  - test FieldVariable.set_from_qp(), .set_from_function() in test_variables()
  - shorten FieldVariable.set_from_function()

- merge pull request #463 from vlukes/fix_extract_surface

  - fix extract_surface.py: remove useless semicolon

- merge pull request #462 from vlukes/fix_ofn_trunk

  - fix Problem.setup_output(): allow dots in domain name

- merge pull request #464 from vlukes/new_extract_edges

  - new_extract_edges.py: extract outline edges of a given mesh
  - update setup.py: add extract_edges.py to aux_scripts list

- merge pull request #461 from vlukes/piezo_homog_example

  - new example - homogenization of a piezoelectric heterogeneous structure
  - update recover_micro_hook*(): suppress outputs when recovering
    microstructures

- merge pull request #465 from vlukes/update_recovery

  - update recovery_micro_hook*()

- merge branch 'update-term-for-no-cells'

  - update Term.evaluate() for no cells in region - do not call term functions

- merge branch 'gen-release-notes'

  - new script/gen_release_notes.py
  - docs: update release tasks
  - add gen_release_notes.py into scripts to install in setup.py

.. _2018.1-2018.2:

from 2018.1 to 2018.2
=====================

- merge branch 'new-time-stepping-solvers'

  - new GeneralizedAlphaTS
  - fix cache name in ElastodynamicsBaseTS._create_nlst_a()
  - new VelocityVerletTS
  - new is_string(), fix IsSave.__call__() for future NumPy
  - add new solvers to elastodynamic.py example, use save_times option

- merge branch 'update-ntc-links'

  - update links to NTC web pages

- merge branch 'stokes-wave-terms'

  - new StokesWaveTerm (dw_stokes_wave), StokesWaveDivTerm (dw_stokes_wave_div)
  - new expand_basis()
  - use expand_basis() in terms

- merge branch 'remove-rcm'

  - remove CloseNodesIterator
  - remove sfepy/linalg/extmods/
  - update setup.py, __init__.py in sfepy/linalg/
  - update FieldVariable.set_from_other() for no CloseNodesIterator
  - remove tests/test_permutations.py
  - docs: sync module index of developer guide with current sources

- merge branch 'misc-updates'

  - update apply_unit_multipliers() for compressibility
  - fix ScipyIterative.__call__() to use eps_a from configuration options -
    update for future SciPy (1.1.x, atol argument)
  - change default eps_a of ScipyIterative
  - set eps_a of iterative scipy solvers in test_linear_solvers.py

- merge branch 'fix-dw_lin_elastic_iso'

  - simplify stiffness_from_lame() to consistently add two dimensions
  - update LinearElasticIsotropicTerm.get_fargs() for current
    stiffness_from_lame() - fixes wrong material shape in single quadrature
    point per cell case
  - update get_pars() in material_nonlinearity.py example

- merge pull request #451 from BubuLK/gen-solvers-table, closes #447

  - add gen_solver_table extension.
  - add preliminary version of gen_solver_table.py
  - updated 'Solvers' section in users guide. Changed "SfePy" typography
    according to (updated) install guide...
  - add solver_table.rst (gen_solver_table.py output) to .gitignore.
  - added generated solver_table (per-partes) to solvers section.
  - update/fix docstrings of time-stepping solvers for solver table
  - update docstrings of nonlinear solvers for solver table
  - update docstrings of eigenvalue problem solvers for solver table
  - docs: update introductory paragraphs of solver descriptions in users guide
  - add link to the corresponding solver class documentation.

- merge branch 'remove-google-analytics'

  - remove google analytics code

- merge branch 'fix-orient-small-cells'

  - fix orient_elements() for cells with very small volumes

- merge branch 'fix-allow-empty-regions'

  - fix Region.finalize() to set .is_empty also for non-cell regions
  - update Region.finalize() to check for emptiness in any case
  - update test_operators()

- merge branch 'cache-refcoors-in-probes'

  - reuse reference coordinates in probes if possible - new
    Probe.get_actual_cache(), update .__init__(), .probe()
  - share geometry in probes in time_poisson_interactive.py - update gen_lines()
  - add docstring to Probe.get_actual_cache(), clean up
  - clean up time_poisson_interactive.py

    - rename gen_lines() -> gen_probes()
    - show the last figure with --show option

- merge pull request #452 from vlukes/new_mumps

  - MUMPS linear solver: new wrapper for C library using ctypes
  - update installation instructions

- merge branch 'dispersion-analysis3'

  - dispersion_analysis.py: new --no-show option
  - dispersion_analysis.py: quote command line in saved options
  - dispersion_analysis.py: free unused memory after each step - allows running
    several examples in parallel without hogging memory
  - dispersion_analysis.py: support custom standard waves logs - new
    save_materials_le(), get_std_wave_fun_le()
  - dispersion_analysis.py: update docstring, update defaults for --range option
  - add show_legends argument to plot_log(), Log.__init__(), update .plot_data()
  - script/plot_logs.py: new --no-legends option
  - dispersion_analysis.py: new --no-legends option
  - dispersion_analysis.py: fix eigensolver calls for all eigenvalues case

- merge pull request #457 from vlukes/upadte_micro_recovery

  - update recover_micro_hook_eps()

- merge pull request #458 from lokik/mmaster

  - fix in script/gen_iga_patch.py

.. _2017.4-2018.1:

from 2017.4 to 2018.1
=====================

- merge branch 'fix-vc++-9-compilation'

  - fix for Visual C++ for Python 9.0

- merge branch 'fix-mtx-comparison'

  - fix matrix comparisons for in-place changes - use sha1 hash instead of only
    id

    - initialize .mtx_digest in LinearSolver.__init__()
    - new _get_cs_matrix_hash(), _is_new_matrix()
    - fix ScipyDirect, PyAMGSolver, PETScKrylovSolver
    - update PyAMGKrylovSolver

  - fix _is_new_matrix() for non-CSR matrices (e.g. PETSc.Mat)

- merge branch 'petsc-ls-options'

  - allow setting additional PETSc options in PETScKrylovSolver
  - update biot_short_syntax.py example

- merge branch 'force-reuse-ls-option'

  - new force_reuse option of PyAMGSolver, PETScKrylovSolver

    - update _is_new_matrix()
    - allow reusing solver objects without checking matrix digests

  - new test_ls_reuse() in tests/test_linear_solvers.py
  - fix linear solver call in poisson_parallel_interactive.py example

- merge pull request #441 from vlukes/fix_recovery

  - fix recover_micro_hook_eps()

- merge branch 'reimplement-advect-div-free-term'

  - new ScalarDotMGradScalarTerm (dw_s_dot_mgrad_s)
  - subclass AdvectDivFreeTerm from ScalarDotMGradScalarTerm, update docstring

- merge branch 'dispersion-analysis2'

  - dispersion_analysis.py: new --conf option for alternative problem
    descriptions

    - new apply_units_le(), set_wave_dir_le()
    - rename define() -> define_le()

  - dispersion_analysis.py: eliminate zeros in matrices

- merge branch 'fix-self-contacts'

  - fix evaluateContactConstraints() for self-contacts
  - turn on global search in ContactTerm.get_fargs()

- merge branch 'numpy-1.14.0-fixes'

  - fix einsum() call in add_eas_dofs()
  - fix transform_basis(), workaround for NumPy 1.14.0
  - use _cmp() in test_consistent_sets() to fix float comparison

- merge pull request #444 from vlukes/mumps_solver

  - new interface to MUMPS linear solver
  - update User's Guide: MUMPS linear solver
  - fix Solver.process_conf()

- merge pull request #443 from heczis/add_hyper_example

  - Add new interactive hyperelastic example

- merge pull request #446 from rc/rewrite-time-stepping-add-dynamics

  - set initial conditions of parameter variables using setter functions

    - update Variables.setup_initial_conditions()
    - new FieldVariable._get_setter()

  - update FieldVariable.time_update() for ._get_setter()
  - new ZeroTerm (dw_zero)
  - update TimeSteppingSolver arguments
  - new examples/linear_elasticity/elastodynamic.py + test
  - new NewmarkTS
  - set quasistatic option of default ts_conf in Problem.set_conf_solvers()
  - update Problem.get_time_solver(), NewmarkTS.__init__() for context argument
  - new output_array_stats()
  - update TimeSteppingSolver to have a generic interface

    - update TimeSteppingSolver.__init__(), .__call__() signatures
    - remove .init_time()

  - update PDESolverApp.call() for generic interface (WIP)

    - new init_fun(), prestep_fun(), poststep_fun()
    - remove init_hook(), step_hook()

  - update StationarySolver for generic interface
  - update PDESolverApp.call() for active_only, new get_vec(), set_state() (WIP)
  - update PDESolverApp.call(): presolve, update init_fun() to return vector
    (WIP)
  - update SimpleTimeSteppingSolver for generic interface

    - new .solve_step0()
    - remove .init_time(), .solve_step()

  - update SimpleTimeSteppingSolver to be base for AdaptiveTimeSteppingSolver

    - new .solve_step(), .output_step_info()

  - update AdaptiveTimeSteppingSolver for generic interface

    - adapt_fun in options has to be a function, not a function name string
    - update .solve_step()
    - remove .__call__()
    - new .output_step_info()
    - update adapt_time_step()

  - update VariableTimeStepper for current AdaptiveTimeSteppingSolver

    - allow setting current step in .set_step()
    - update .advance(), .iter_from_current()
    - new .iter_from()

  - fix Equations.time_update() for no active_only and changing EBC DOFs -
    update Problem.update_equations()

  - fix error reporting in test_term_call_modes()
  - update extract_time_history() to use actual step times
  - update linear_elastic_damping.py example for current TS solvers - use HDF5
    for output
  - update balloon.py example to use ts.adaptive
  - update adapt_time_step() docstring
  - update laplace_time_ebcs.py example for current TS solvers
  - new standard_ts_call(), decorate .__call__() of TS solvers
  - fix FieldVariable.__init__() to initialize .history attribute
  - use quasistatic TimeStepper in StationarySolver.__init__()
  - add status argument to TimeSteppingSolver.__init__() for consistence
  - move prepare_save_data(), prepare_matrix() into problem.py
  - move get_initial_state() -> Problem.get_initial_state()
  - move time-stepping solver management from PDESolverApp to Problem

    - update PDESolverApp.call()
    - replace Problem.nls with Problem.solver - a TS solver instance
    - update .solve()
    - new .get_tss(), .get_tss_functions(), .get_nls(), .get_ls()
    - update Problem.__init__(), .reset() .set_equations(),
      .set_equations_instance(), .set_conf_solvers(), .set_solver(),
      .try_presolve(), .get_solver(), .is_linear(), .set_linear()
    - update .init_solvers() - add ts_conf argument
    - remove .get_time_solver()
    - new State.get_vec(), .set_vec()

  - unite BasicEvaluator with LCBCEvaluator in Evaluator, clean up
  - update Problem.get_evaluator(), PETScParallelEvaluator for Evaluator
  - update tests to use nls.fun() instead of (Basic)Evaluator
  - update solve_pde(), PDESolverApp.call() to have status argument
  - remove nls, ls, ts, auto_solvers arguments of Problem.__init__() - update
    .from_conf(), .copy(), .create_subproblem()
  - initialize solver configuration attributes in Problem.__init__()
  - update Problem.set_solver() to set nonlinear solver functions - new
    .get_nls_functions()
  - initialize .ts in Problem.reset(), .set_solver()
  - update Evaluator to check for LCBCs at run-time (remove .has_lcbc,
    .mtx_lcbc)
  - update tests for new solver status handling
  - remove make_implicit_step()
  - update SchurGeneralized for no active DOF indices at solver creation time -
    update SchurComplement parameters initialization
  - new Problem.block_solve() replacing EquationSequenceSolver - update
    Problem.solve()
  - update make_l2_projection_data(), make_h1_projection_data() for current
    Problem
  - update tests for current Problem
  - update (interactive) examples for current Problem
  - update test_install.py for changed output
  - update MultiProblem for no active DOF indices at solver creation time
  - update for not passing time stepper in user argument in
    Problem.set_equations() - update Term.assign_args(), .time_update(),
    create_evaluable()
  - replace prepare_save_data() with new make_is_save() - new IsSave

    - update is_sequence()
    - update Problem.get_tss_functions()

  - update examples for save_times option and time-stepper default verbosity
    change
  - docs: update users guide for save_times option
  - new TimeStepper.set_substep_time(), .restore_step_time()

  - new BatheTS
  - simplify NewmarkTS, BatheTS by subclassing new ElastodynamicsBaseTS
  - docs: update tutorial
  - docs: update solvers sections in users guide
  - clean up: move solver related-functions together in Problem
  - docs: tweak for users guide changes in master

- merge branch 'fix-test-ls-reuse'

  - update test_ls_reuse() for #446

- merge branch 'remove-petsc-worker'

  - remove unused sfepy/solvers/petsc_worker.py - see
    125d59dd82c0f2e4c88031c7c58e2dfa255c8cf8
  - update sfepy/solvers/__init__.py
  - docs: sync module index of developer guide with current sources

- merge pull request #448 from vlukes/update_nonlin_homog_example

  - update Problem.get_evaluator(): allow a user evaluator Class specified in
    problem options
  - update nonlinear homog. example: adapt to the altered solvers

- merge branch 'fix-hdf5-saving-only-some-steps', closes #445

  - allow saving without step 0 in HDF5MeshIO.write()
  - update HDF5MeshIO for not saving all steps

    - new HDF5MeshIO._get_step_group_names()
    - update .read_times()
    - update ._get_step_group() .read_data_header() for no step 0
    - update .read_time_history() for missing steps

- merge pull request #449 from heczis/update_interactive_example

  - update hyperelastic_tl_up_interactive.py example for current Problem

- merge branch 'update-web-docs'

  - docs: update support section
  - docs: move (old) featured applications under examples
  - update script/gen_gallery.py for current Problem

.. _2017.3-2017.4:

from 2017.3 to 2017.4
=====================

- merge pull request #418 from lokik/python3-compatibility

  - Python 3 compatibility: region.py
  - Python 3 compatibility: update_dict_recursively

- merge pull request #420 from lokik/master

  - sfepy.base.parse_conf: fix and test

- merge branch 'contacts'

  - implements a penalty based contact term
  - new examples/linear_elasticity/two_bodies_contact2d.py
  - new sfepy/terms/terms_contact.py - new ContactTerm
  - allow strings as special material arguments in Term.check_shapes()
  - new sfepy/mechanics/extmods/__init__.py
  - new sfepy/mechanics/extmods/contres.{c, h}
  - new sfepy/mechanics/extmods/ccontres.pyx
  - new sfepy/mechanics/extmods/setup.py
  - new get_longest_edge_and_gps() cython function
  - new get_AABB() cython function
  - new init_global_search() cython function
  - new evaluate_contact_constraints() cython function
  - new assemble_contact_residual_and_stiffness() cython function
  - update Term.evaluate(), .assemble_to() for sparse vector data
  - clean up sfepy/base/plotutils.py
  - fix print_matrix_diff() for matrices in CSC format
  - update spy(): make dots visible, fix axes limits, do not shift points -
    with matplotlib 1.5.1
  - return extra matrix from Term.assemble_to(), do not assemble it there
  - update Equations, Equation for extra matrices

    - update Equations.evaluate(), .eval_tangent_matrices()
    - update Equation.evaluate()

  - rename two_bodies_contact2d.py -> two_bodies_contact.py
  - update two_bodies_contact.py for 3D, generate bodies by new gen_two_bodies()
  - update/fix contres.{c, h} for 3D
  - update ContactTerm for 3D
  - add docstring to Term.assemble_to(), clean up
  - new test for two_bodies_contact.py example
  - move active DOF treatment from ContactTerm.get_fargs() to
    Term.assemble_to()
  - update _test_single_term() for dynamic connectivity terms -
    test_term_call_modes() tests pass for dw_contact
  - docs: sync module index of developer guide with current sources

- merge branch 'plot-boundary-quadratures'

  - script/plot_quadratures.py: new --boundary option, update _get_bqp(),
    plot_quadrature()
  - script/plot_quadratures.py: new ---show-labels, --print-qp options, new
    label_points(), update plot_quadrature()

- merge pull request #422 from vlukes/update_doc

  - update users_guide: add links to PyAMG and PETSc documentation
  - update users_guide: remove 'ls.petsc_parallel' section

- merge branch 'solver-context'

  - update Solver.__init__() and subclasses: add context argument, update
    LinearSolver.__call__() signature
  - update Problem.init_solvers() to pass self as context
  - update linear solvers for context argument

    - update standard_call(), petsc_call() decorators
    - update ScipyIterative, PyAMGKrylovSolver, PETScKrylovSolver,
      SchurGeneralized, MultiProblem

  - update Oseen for context argument
  - update setup_precond() in biot_short_syntax.py example for context argument

- merge branch 'embed-shell'

  - add frame argument to python_shell()
  - new shell(), ipython_shell()

- merge pull request #424 from BubuLK/c-cython-warnings, closes #406

  - Updated array/pointer cdef to follow new Cython parser conventions (removed
    Cython warnings).
  - Removed unused variables definitions.
  - Replace abs()->fabs().
  - Updated cmesh.{c,h} explicit casting.
  - Updated explicit casting (discrete/common).
  - Updated explicit casting (iga), fixed typo.
  - Added cython explicit type casting (to remove warnings).
  - Add "unreachable code" explicit marks () to disable warnings.
  - Removed duplicated compiler options defines (sfepy_common).
  - Fixed comparsion bug (from unreachable code).
  - Cleanup mesh.c explicit type casting.
  - Updated OS detection defs.
  - Updated defs/casting according to issue comments.
  - remove unused variables in C code - closes pull request #423 from
    vlukes/c_clean_up
  - resolve remaining warnings

- merge branch 'scikit-umfpack-version'

  - update _scikit_umfpack_version()

- merge pull request #425 from BubuLK/c-compiler-flags

  - Fixed typo in sfepy_common library macros.

- merge branch 'allow-empty-regions-option'

  - new allow_empty_regions Problem configuration option, update
    Problem.from_conf(), .set_regions()
  - docs: update users guide

- merge pull request #428 from BubuLK/Appveyor-IGA

  - Add IGAkit install/build to Appveyor config.

- merge pull request #426 from rc/small-fixes

  - fix real definition in parse_conf.py, new cmplx definition
  - initialize time_stats in Newton.__call__()
  - fix spelling: rezidual -> residual
  - print time_stats in fixed order in Newton.__call__()
  - fix LogPlotter.__call__() docstring

- merge branch 'small-fixes-2'

  - add show_mesh_info.py into scripts to install in setup.py
  - script/show_mesh_info.py: fix misleading description, show real centre
  - clean up sfepy/mechanics/units.py
  - fix density definition in units_of_quantities, more num_prefixes
  - improve value shape checking, error reporting in H1NodalMixin.set_dofs()

- merge branch 'active-only-option'

  - new active_only Problem configuration option - update Problem.from_conf()
  - obey active_only in Problem.copy(), .create_subproblem()
  - add active_only argument to  Problem.evaluate() and related functions

    - update Problem.create_evaluable(), .eval_equations()
    - update create_evaluable(), eval_equations(), eval_in_els_and_qp(),
      assemble_by_blocks()

  - new Problem.get_ebc_indices()
  - update PETScParallelEvaluator() for Problem.get_ebc_indices()
  - fix matrix diagonal in BasicEvaluator.eval_tangent_matrix() for no
    active_only
  - fix Problem.solve() for no active_only
  - fix making full vector in BasicEvaluator for no active_only, fix
    .eval_residual(), .eval_tangent_matrix()
  - update PETScParallelEvaluator for apply_ebc_to_matrix() call in
    BasicEvaluator
  - fix create_adof_conns(), apply_ebc_to_matrix() for EPBCs

    - update Problem.get_ebc_indices()
    - update BasicEvaluator.eval_tangent_matrix()

  - fix equations in poisson_periodic_boundary_condition.py example
  - set active_only to False in poisson_periodic_boundary_condition.py example
  - move apply_ebc_to_matrix() into sfepy/discrete/evaluate.py
  - update Problem docstring

- merge pull request #430 from rc/update-convert-mesh

  - script/convert_mesh.py: new --2d option
  - docs: document --2d option in preprocessing section
  - docs: add data files
  - docs: mention legacy VTK reader 2D detection feature
  - fix printing writable mesh formats in for_format()

- merge branch 'fix-poisson_parallel_interactive'

  - poisson_parallel_interactive.py: fix for moved apply_ebc_to_matrix() - see
    active-only-option branch

- merge pull request #431 from vlukes/new_surface_grad_term

  - new ev_surface_grad and ev_surface_div terms

- merge pull request #432 from vlukes/update_complex_eval_at

  - update fields.evaluate_at() for evaluating complex fields
  - new test for fields.evaluate_at()

- merge branch 'small-fixes-3'

  - output last step KSP stats in PETScNonlinearSolver
  - fix ANSYSCDBMeshIO for meshes with both tetrahedra and hexahedra

- merge pull request #433 from vlukes/save_mesh_per_matid

  - update mesh conversion: extract cells by material id
  - update preprocessing tutorial - new "save-per-mat" arugment to conversion
    script

- merge pull request #434 from heczis/abaqus_ax_elm

  - Add reading of axisymmetric elements to AbaqusMeshIO

- merge branch 'show-mesh-info-euler-ncomp'

  - script/show_mesh_info.py: show Euler characteristic
  - script/show_mesh_info.py: show medians of volumes, update formatting
  - script/show_mesh_info.py: show Euler characteristics of mesh volume and
    surface
  - script/show_mesh_info.py: show numbers of volume/surface components
  - script/show_mesh_info.py: update output formatting

- merge pull request #436 from vlukes/new_vtk_probe

  - new ProbeFromFile class: init VTK probe using a given file

- merge pull request #437 from vlukes/recovery_eps

  - new recover_micro_hook_eps(): recover a real sized microstructure

- merge pull request #435 from rc/tetgen-remesh-option

  - script/convert_mesh.py: new --remesh option

- merge pull request #438 from vlukes/update_recovery

  - update calls of microproblems: pass arguments to define() at the microlevel
  - fix recover_micro_hook_eps(): fix for incorrect microstructure size
  - let the code be friedly to pep8

- merge branch 'small-fixes-4'

  - check that all facets are on surface in SurfaceField._check_region(),
    improve messages
  - allow quoting command line items in save_options() - new quote_command_line
    argument
  - fix Mesh.create_conn_graph() to obey verbose argument

- merge branch 'dispersion-analysis'

  - new ElasticWaveTerm (dw_elastic_wave), _build_wave_strain_op()
  - new ElasticWaveCauchyTerm (dw_elastic_wave_cauchy),
    _build_cauchy_strain_op()
  - new which option instead of hard-coded value in ScipyEigenvalueSolver
  - new examples/linear_elasticity/dispersion_analysis.py
  - new apply_unit_multipliers()
  - new lame_from_stiffness(), youngpoisson_from_stiffness()
  - update test_stiffness_tensors()
  - update Log, LogPlotter for plots with varying line properties

    - update Log.__init__(), .add_group(), .plot_data()
    - update LogPlotter.process_command()

  - save plot properties header in Log.__init__(), update read_log(), plot_log()
  - add raw_log_save_name option to BandGaps
  - new save_raw_bg_logs(), update AcousticBandGapsApp.call() for saving raw
    logs
  - update plot_log() for plotting to given list of axes
  - script/plot_logs.py: update --rc option, update ParseRc to use eval()

    - values with commas work
    - update for current plot_log()

  - improve BandGaps.save_log()
  - script/plot_logs.py: new --groups option, update plot_log()

  - fix dense eigh() call in ScipyEigenvalueSolver.__call__()
  - set accuracy for periodic vertex matching in band_gaps_conf.py
  - use basic SI units in band_gaps.py example, small updates
  - use basic SI units in band_gaps_rigid.py example, small updates
  - update test_install.py for updated units in phononic examples

.. _2017.2-2017.3:

from 2017.2 to 2017.3
=====================

- merge branch 'regions-update'

  - allow '-' in region names - update parsing code + test
  - fix reading of vertex sets (nodal bcs) in
    HDF5MeshIO.read_mesh_from_hdf5() - group argument is no longer overwritten

- merge pull request #395 from lokik/master

  - svec (output buffer) argument for variables.strip_state_vector
  - code lint in discrete/variables

- merge branch 'improve-parallel'

  - improve information outputs
  - speed-up assemble_rhs_to_petsc(), assemble_mtx_to_petsc() by removing loops
  - measure and report global domain/fields setup time in parallel examples

- merge pull request #403 from vlukes/homog_mpi

  - update homog. engine - rearrange functions, define new class
    HomogenizationWorker
  - update engine.py to comply pep8
  - update homog. engine: define numdeps as dict instead of list
  - new sfepy/base/multiproc_mpi.py - classes and functions for MPI
    parallelization
  - update sfepy/base/multiproc.py - unify multiproc. modules
  - update homog. engine and application for MPI computation
  - update sfepy/base/multiproc.py and sfepy/homogenization/homogen_app.py to
    comply pep8
  - fix test_homogenization_engine.py - new structure of homog. engine

- merge branch 'iterative-ls-precond'

  - make preconditioners for ScipyIterative solver actually usable - change
    option precond (a matrix-like) to setup_precond (a callable)
  - add iteration callback to ScipyIterative, calls user callback if provided
  - fix PETScKrylovSolver.__call__() for no initial guess
  - simple.py: allow additional options (to use with PETSc options)
  - new init_petsc_args(), used in PETScKrylovSolver.__init__()
  - new Solver.set_field_split()
  - allow slices in PETScKrylovSolver.set_field_split()
  - set field split data in Problem.solve()
  - obey verbose option in ScipyIterative
  - new examples/multi_physics/biot_short_syntax.py + test
  - update PETScNonlinearSolver to return same status information as Newton

    - set manually the solution from the update in case the KSP did not
      converge

- merge branch 'expand-nodes-node-by-node', closes #404

  - arrayize function values in {H1HierarchicVolumeField,
    H1NodalMixin}.set_dofs()
  - fix shape for NumPy 1.13.0 in EquationMap.map_equations()
  - fix shape for NumPy 1.13.0 in FieldVariable.setup_initial_conditions()
  - update expand_nodes_to_equations() to use node-by-node ordering of DOFs
  - unify shape of values returned by Field.set_dofs() implementations

    - update H1HierarchicVolumeField, H1NodalMixin, IGField
    - update for node-by-node ordering: use (n_nodes, n_components)

  - update MRLCBCOperator.setup() for node-by-node ordering
  - update EBC/LCBC functions in examples for node-by-node ordering
  - update test_ebc_functions() for vector variables
  - docs: update users guide for node-by-node ordering

- merge pull request #407 from BubuLK/doc-tutorial, closes #379, #401

  - updated tutorial.rst accordint to issue #379.
  - updated primer.rst
  - updated linear_elasticity_interactive.py according to tutorial.rst
  - added sfepy-wrapper label to user_guide.rst.
  - updated tutorial/installation.rst.
  - updated tutorial/basic-usage.rst
  - removed "$" from cli examples.

- merge pull request #409 from BubuLK/Sphinx-conf

  - replaced custom 'ipython_console_highlighting.py' with standard one.
    Corrected IPython console outputs in tutorial.rst.
  - replaced deprecated pngmath extension with imgmath.
  - updated conf.py to new LaTeX customization scheme
    (latex_preamble->latex_elements).

- merge branch 'docs-fix-term-table', closes #399

  - force longtable in script/gen_term_table.py
  - add LaTeX page breaks around tables in script/gen_term_table.py

- merge pull request #410 from vlukes/tri_tetra_elements

  - new option to convert_mesh.py script: '-t' convert quad/hexa elements to
    tri/tetra

- merge branch 'improve-ls'

  - add solver name to messages in ScipyIterative.__call__()
  - reuse KSP instance in PETScKrylovSolver for multiple solves with one matrix

    - speed-up, especially for direct solver preconditioning
    - update .__init__(), .__call__()

  - allow additional options in ScipyIterative
  - prepare ScipyIterative for future scipy support of both rtol and atol
  - tweak verbosity levels in ScipyIterative, print number of iterations
  - fix ScipyDirect to obey presolve option
  - do not store matrix in ScipyDirect
  - update PyAMGSolver to use id of matrix for solver reuse check - do not
    store matrix
  - remove bit-rotten/obsolete PETScParallelKrylovSolver
  - update standard_call(), petsc_call() to return number of iterations in
    status

    - update ScipyIterative.__call__()
    - supported where possible, closes #216

  - return total number of linear solver iterations in status - update Newton,
    PETScNonlinearSolver
  - update test_solvers() to report numbers of linear solver iterations
  - update test_install.py for updated nls status

- merge branch 'fix-for-sympy-1.1'

  - fix Quantity.__init__() for sympy 1.1

- merge pull request #412 from vlukes/homog_mpi

  - update parallel MPI homogenization - add features for solution of
    multiscale problems
  - update get_homog_coefs_nonlinear() for MPI parallel computation
  - update MPI homogenization - improve efficiency of MPI communication
  - update multiproc_mpi.py: update classes RemoteDict and RemoteDictMaster,
    fix typos in logs
  - rearrange multiprocessing modules
  - update homog. functions - reaaranged multiprocessing modules
  - update save_mappings() and get_mappings() - rearranged multiprocessing
    modules
  - new `simple_homog_mpi` solver - allows to run parallel micro-macro coupled
    simulation

- merge branch 'misc-updates'

  - improve error message in Term.check_shapes() to include actual shapes
  - fix PETScNonlinearSolver.__call__() for no SNES.getFunctionNorm() - see
    https://bitbucket.org/petsc/petsc4py/commits/1ffe3970457cf66c4354ca2d4601852ea06999b5

- merge pull request #414 from vlukes/homog_mpi_fix

  - fix get_homog_coefs_nonlinear() in homogenization/micmac.py - mpi switch

- merge branch 'mesh-entity-volumes'

  - new mesh_get_volumes() C function - new _det3x3(), _tri_area(), _aux_hex()
  - new CMesh.get_volumes()
  - script/show_mesh_info.py: show only names of nodal BCs
  - script/show_mesh_info.py: new --detailed option, shows entity volumes
  - update mesh_get_volumes() for approximate bilinear face area computation
  - new test_entity_volumes() in tests/test_cmesh.py

- merge pull request #415 from vlukes/homog_mpi_fix

  - update homog. engine - replace 'chunk_size' option by 'chunks_per_worker'
  - update nonlin. homog. example - pep8 code style
  - update multiproc. code - 'thread' in names is obsolete, replaced by 'proc'

- merge branch 'pyamg-krylov'

  - allow additional method/solve options in PyAMGSolver, support callbacks -
    add iteration callback to PyAMGSolver, calls user callback if provided
  - new PyAMGKrylovSolver - interface to PyAMG Krylov solvers
  - update tests/test_linear_solvers.py

- merge branch 'petsc-user-precond'

  - support user-defined preconditioners in PETScKrylovSolver - new
    setup_precond option
  - update tests/test_linear_solvers.py

- merge pull request #417 from BubuLK/deployment-CI, closes #350

  - add updated configs for Travis/AppVeyor testing.
  - changed run_test.py call (ps->cmd).
  - removed (obsolete) x86 arch to speedup test.

- merge branch 'remove-shaper'

  - remove obsolete shaper.py
  - remove obsolete sfepy/optimize/ - remove __init__.py, free_form_def.py,
    setup.py, shape_optim.py
  - update setup.py, sfepy-run for no shaper.py
  - update problem description file transforms
  - docs: update users guide for no shaper.py
  - update sfepy/setup.py for no sfepy/optimize/

- merge branch 'update-log-live-plot', closes #131

  - use threading in LogPlotter.__call__() to call .poll_draw() -
    update.poll_draw() to sleep between canvas updates (replaces gobject
    timeout)
  - allow plt.tight_layout() failure in LogPlotter.process_command()
  - update Log.__init__() for no fixed matplotlib backend dependence
  - update LogPlotter for new sleep argument
  - update Log for new sleep argument
  - update live_plot.py example to use aggregate, sleep options, clean up

- miscellaneous updates:

  - docs: sync module index of developer guide with current sources
  - update version string in get_basic_info() to conform with PEP 440

.. _2017.1-2017.2:

from 2017.1 to 2017.2
=====================

- merge pull request #369 from rc/fix-variable-history-advance

  - initialize history of variables in get_initial_state()

    - update make_implicit_step()
    - update SimpleTimeSteppingSolver, AdaptiveTimeSteppingSolver

  - fix Variable.advance() to initialize current step data

- merge pull request #370 from heczis/master

  - remove the unused method Problem.init_variables

- merge pull request #368 from vlukes/update_homogen

  - update: simplified and unified implementation of some homogenized
    coefficients
  - fix the homogenization example: perfusion_micro.py

- merge pull request #372 from vlukes/fix_material_shape_change

  - fix changing of the material shape

- merge pull request #373 from vlukes/piezo_strain

  - new PiezoStrainTerm
  - update piezo-elasticity example

- merge pull request #376 from vlukes/fix_truncation

  - fix: avoid number truncation in region definitions

- merge branch 'docs-main-page'

  - update support section
  - remove link to obsolete wiki pages
  - add link to anaconda installation instructions to main page

- merge branch 'fix-vtk-source-mayavi-4.4', closes #292

  - update VTKMeshIO.read_data() to read cell data, small tweaks
  - fix GenericFileSource._reshape() for single-axis data
  - fix GenericSequenceFileSource

    - update .read_common()
    - remove .create_source()
    - new .file_changed()
    - initialize .io in GenericFileSource.__init__()

  - update create_file_source() to work around a Mayavi 4.4.x issue

- merge branch 'fix-coefs-to-latex'

  - fix Coefficients._save_dict_latex() for scalars and general data
  - clean up: raise exception with message

- merge pull request #377 from rc/fix-variable-state-data-sharing

  - fix data copying in Variable.advance() - bad interaction with State

- merge pull request #380 from heczis/fix_doc_python3

  - fix things to be compatible with both python 2 and 3

- merge branch 'term-report-missing-virtual'

  - new Term.get_str()
  - use Term.get_str() in Term, Equations
  - update Term.evaluate() to report missing virtual variable in 'weak' mode

- merge pull request #383 from vlukes/update_homog_doc

  - update homogenization examples, add references

- merge pull request #385 from vlukes/change_shape_ev_grad

  - change the shape of the gradient array provided by 'ev_grad', now: (n_el,
    n_qp, dim, n_c)

- merge pull request #386 from vlukes/replace_copydata_corr

  - replace corrector CopyData by the more general one with name CorrEval

- merge pull request #388 from vlukes/fix_meshio_msh

  - fix Msh2MeshIO.read() to discard '2_2' elements

- merge pull request #390 from rc/docs-sfepy-at-python-org

  - docs: update for sfepy(at)python.org

- merge pull request #391 from rc/fix-project-by-component

  - fix project_by_component() for general tensor shape
  - new test_project_tensors()

- merge branch 'fix-ansys-cdb'

  - fix ANSYSCDBMeshIO.read(), convert tetras as degenerate hexas to tetras
  - new look_ahead_line()
  - update ANSYSCDBMeshIO.read() to determine true number of fields

    - fixes reading files with wrong nblock/eblock information
    - update make_format()

  - fix remapping of nodal bcs in ANSYSCDBMeshIO.read() for qtetras, qhexas

- merge branch 'misc-updates'

  - copying subclasses of problem
  - numpy compatibility
  - problem.make_full_vec accept vec argument.
  - new Container.__add__(), .__iadd__(), test_container_add()
  - update Viewer.build_mlab_pipeline() to add mat_id to source if not filtered

- merge pull request #375 from {lokik,rc}/save-custom-data-to-hdf5

  - saving custom structured data to h5 file in problem.save_state
  - ioutils.enc and ioutils.dec utf strings compatibility
  - HDF5ContextManager
  - IGDomain reading and writing from HDF5 file
  - asserting equality of complex structures
  - HDF5 reading and writing
  - faster assert_equals.
  - default mesh argument for HDF5MeshIo.read()
  - storing data to hd5 using softlinks.
  - clean up iga/domain.py, iga/io.py
  - clean up and reorganize HDF5MeshIO
  - clean up ioutils.py, reorganize new functions/classes
  - move assert_equals() into TestCommon.assert_equal(), update
  - clean up and update test_hdf5_meshio()

- merge branch 'docs-devel-page'

  - update development tab, new topics section
  - update copyright info
  - add more topics to development tab

- merge pull request #393 from heczis/fix_generators_next

  - fix generators' next method calls in script/gen_gallery.py
  - fix generators' next method calls in sfepy/base/log.py
  - fix generators' next method calls in sfepy/application/application.py

- merge branch 'band-gaps-ranges'

  - new get_gap_ranges(), use in BandGaps.__call__()
  - simplify plot_gap() by using gap_ranges, move text output to plot_gaps() -
    update AcousticBandGapsApp.plot_band_gaps(), .plot_dispersion(), use tight
    layout
  - fix plot resources for matplotlib >= 1.5.1

- merge pull request #398 from vlukes/tutorial_preproc

  - new tutorial: preparing meshes using FreeCAD/OpenSCAD and Gmsh
  - new "merge" option in `convert_mesh.py` - remove duplicate vertices

- merge pull request #400 from vlukes/update_homog

  - update homogenization to allow saving "pi" correctors

- merge pull request #397 from BubuLK/doc-install, closes #366, #382

  - updated Install doc: - issue #382 - issue #366 (?) - other misc doc cleanup
  - doc cleanup - bugfixes - updates according to PR comments - updated
    sections structure
  - updated and re-structured install doc
  - updated Anaconda instructions
  - add link to conda-forge on downloads page
  - add link to install doc
  - add direct link conda-forge SfePy packages

- miscellaneous updates:

  - update mailing lists addresses in release tasks

.. _2016.4-2017.1:

from 2016.4 to 2017.1
=====================

- merge pull request #355 from heczis/dont_redefine_help

  - fix redefining help

- merge pull request #354 from lokik/master

  - new numpy version compatibility

- merge pull request #359 from heczis/python_cmd_in_test_install

  - python -> python2 in test_install.py

- merge pull request #360 from heczis/logging_in_test_install

  - use the logging module for output in test_install.py

- merge branch 'fix-data-from-qp-shape'

  - fix caching in Integrals.get()
  - fix vertex data reshaping in GenericFileSource.add_data_to_dataset() - new
    GenericFileSource._reshape()
  - fix data shape in FieldVariable.set_data_from_qp()
  - update nodal_stress() in its2D_3.py example

- merge pull request #361 from vlukes/fix_doc_splinebox

  - fix splinebox example - "Mesh parametrization"

- merge pull request #362 from vlukes/update_splinebox

  - update splinebox - parameterization of an arbitrary field
  - update splinebox.py to pass the pep8 check
  - new splinebox test - check field parametrization

- merge pull request #356 from rc/docs-conda-forge-install

  - docs: update installation instructions for conda-forge releases

- merge branch 'python-3.6-fixes'

  - update .travis.yml to test with Python 3.5, 3.6
  - fix integer division errors in shapes/indices
  - fix errclear() for Python 3.6, remove useless line from errput()
  - catch ValueError in Term.call_get_fargs(), .call_function()
  - fix more integer division errors

- merge branch 'fix-petsc-sub-precond-type'

  - fix default sub_precond value in PETScKrylovSolver
  - fix sub_precond argument in parallel examples
  - report number of iterations in PETScKrylovSolver, PETScParallelKrylovSolver

- merge branch 'fix-integer-divisions'

  - fix integer division errors

- merge branch 'problem-docstring'

  - improve active_only description in Problem docstring
  - describe arguments of Problem.__init__() in class docstring

-  merge pull request #364 from vc12345679/master

   - fix bug: "Python.h" Include Path Error - use
     `sysconfig.get_config_var('INCLUDEPY')`, instead of
     `sys.prefix+'include'+'python'+version`, to obtain include path of
     'Python.h'

- merge pull request #365 from vlukes/fix_save_regions

  - fix saving surface regions

- miscellaneous updates:

  - docs: update release tasks

.. _2016.3-2016.4:

from 2016.3 to 2016.4
=====================

- merge branch 'fix-lcbc-several-fields'

  - fix _s_describe() for zero area facets
  - fix/improve geme_invert3x3(), geme_invert4x4() for singular matrices
  - fix LCBCOperators.finalize() to keep correct ordering of variables - update
    _dict_to_di()
  - new test_stokes_slip_bc() in tests/test_lcbcs.py

- merge pull request #343 from 'rc/debug-on-error'

  - new debug_on_error()
  - run_tests.py: rename --debug option to --raise
  - run_tests.py: new --debug option - run debugger on exception
  - docs: update for updated options of run_tests.py
  - update .travis.yml for --raise
  - new --debug option in top level scripts - update extractor.py, homogen.py,
    phonon.py, postproc.py, probe.py, schroedinger.py, shaper.py, simple.py

- merge branch 'empty-fe-surface'

  - update FESurface.__init__() for empty region
  - update CMapping.describe() for empty region
  - update FieldVariable.evaluate() for empty region

- merge branch 'non-penetration-penalty-term'

  - new NonPenetrationPenaltyTerm (dw_non_penetration_p)
  - new examples/navier_stokes/stokes_slip_bc_penalty.py + test
  - small tweaks in stokes_slip_bc.py example, reference
    stokes_slip_bc_penalty.py

- parallel support:

  - update Domain.create_regions() for empty regions - new allow_empty argument
  - update Mesh.from_region() to preserve nodal BCs
  - update Region.setup_from_highest() to always succeed when allowed empty
  - update create_task_dof_maps() for easier debugging of partitioning problems

    - new save_inter_regions, output_dir arguments
    - update distribute_fields_dofs()

- merge branch 'parallel-examples-update'

  - new remove_files_patterns()
  - new save_options()
  - poisson_parallel_interactive.py: new --save-inter-regions options, save
    options
  - biot_parallel_interactive.py: new --save-inter-regions options, save
    options

- merge pull request #340 from 'lokik/master'

  - equations.add_equation method

- merge pull request #348 from vlukes/fix_pt_open, closes #342

  - fix pytables compatibility issue: openFile -> open_file, createGroup ->
    create_group, ...

- merge pull request #346 from vlukes/fix_set_coors

  - fix setting field coordinates for higher order elements

- merge pull request #337 from vlukes/update_tests

  - display the test file numbers and test numbers to get a better view in a
    debug mode

- merge pull request #347 from vlukes/nonsym_biot

  - new nonsymmetric mode of BiotTerm

- merge pull request #349 from rc/fix-biot-ccode

  - fix op_nonsym_biot()
  - fix dw_biot_grad() for compiling on windows

- merge pull request #351 from rc/he-clean-up

  - remove obsolete CorrectorsPermeability
  - fix insert_sub_reqs() for arbitrary order of leaf requirements
  - simplify insert_sub_reqs() - remove too strict circular dependency check
  - new tests/test_homogenization_engine.py: new test_dependencies()

- merge branch 'hanging-nodes'

  - conflicts: sfepy/discrete/fem/fields_base.py
  - support basis transforms in FEField, VolumeMapping, PolySpace

    - new FEField.basis_transform attribute, FEField.set_basis_transform()
    - update FEField.get_base(), .create_mapping()
    - update VolumeMapping.get_mapping()
    - new transform_basis()
    - update PolySpace.eval_base()

  - new sfepy/discrete/fem/refine_hanging.py - initial 2D version, WIP

    - new find_level_interface(), refine_region(), find_facet_substitutions(),
      refine(), do_connectivity_substitutions(), eval_basis_transform()

  - manage connectivity substitutions and unused DOFs in FEField

    - move do_connectivity_substitutions() into new FEField.substitute_dofs()
    - new FEField.econn0, .unused_dofs attributes
    - new FEField.restore_dofs(), .restore_substituted()

  - update EquationMap.map_equations() to omit unused field DOFs from active
    DOFs
  - update FieldVariable.get_full() to restore unused field DOFs
  - fix Region.cells setter for empty cell regions
  - fix CMesh.get_incident() for no incident entities
  - update FEField._setup_esurface() to setup .eedges in 3D
  - fix refine_edges_3_8
  - update Variable._set_kind() to always initialize .dof_name

    - use the variable name as the DOF name for parameter variables without a
      primary
    - update ._setup_dofs()

  - update EquationMap.map_equations() to obey unused DOFs in no EBC case

    - update ._init_empty()
    - new ._mark_unused()

  - fix PointsProbe.__init__() to force C-contiguous order
  - new tests/test_refine_hanging.py: new test_continuity() test

    - new eval_fun(), _gen_lines_2_4(), _gen_grid_3_8(), _build_filenames()

  - move body of FEField.substitute_dofs() into new
    H1NodalMixin._substitute_dofs()
  - move eval_basis_transform() -> H1NodalMixin._eval_basis_transform()
  - update FEField.substitute_dofs(), .restore_dofs() for storing substitutions

    - new .stored_subs
    - evaluate and set basis transform in FEField.substitute_dofs()

  - preserve indices of non-refined cells

    - update refine_region() - new _interleave_refined()
    - update find_level_interface(), find_facet_substitutions(), refine()

  - update refine_region() to preserve vertex groups of non-refined cells
  - new test_preserve_coarse_entities() in tests/test_refine_hanging.py
  - new examples/diffusion/laplace_refine_interactive.py
  - update test_install.py to test laplace_refine_interactive.py example

- merge branch 'hessian-lagrange-basis'

  - new LagrangeSimplexPolySpace._eval_hessian()
  - new LagrangeTensorProductPolySpace._eval_hessian(), update .__init__()
  - update PolySpace.eval_base(), LagrangePolySpace._eval_base() for 2.
    derivatives
  - new test_hessians() in tests/test_poly_spaces.py

- merge pull request #352 from vlukes/ulf_homog

  - fix set_mesh_coors() - initiate coors_act array
  - update Problem.solve() to allow disabling materials update in a given time
    step
  - update periodic.match_() for caching matching coordinates
  - clean-up: periodic.py
  - new multiproc module - global multiprocessing management
  - update saving of field mappings, allow sharing data among processes
  - update homog. engine - compute coefficients for multiple micro
    configurations at once
  - clean up in homogenization modules
  - update homogenization engine: volumes are calculated as the coefficients
  - fix the test of homogenization_perfusion.py
  - new non-linear homogenization example
  - update homogenization engine test: check splitting/merging chunks
  - fix band_gaps_app.py and the related test

- merge pull request #353 from vlukes/update_homog_example

  - update linear homogenization examples

- miscellaneous updates:

  - docs: update latest snapshot link, closes #344
  - add custom view for stokes_slip_bc_penalty.py example to
    script/gen_gallery.py
  - fix streamline position in plot_velocity() - regression by 63171ad
  - sfepy-run: fix --version option for Python 3
  - fix SurfaceMomentTerm - add .arg_shapes, make shift special material
    parameter
  - fix AcousticBandGapsApp.__init__() for non-file problem configuration
  - fix typo in phononic examples
  - docs: sync module index of developer guide with current sources
  - do not omit linear_elastic_mM.py in script/gen_gallery.py
  - docs: update release tasks

.. _2016.2-2016.3:

from 2016.2 to 2016.3
=====================

- merge pull request #330 from 'vlukes/fixdoc'

  - fix docstrings in DotProductVolumeTerm, VectorDotGradScalarTerm

- merge pull request #331 from 'takluyver/py3' and 'rc/py3', closes #164

  - Python 2.7 and 3.4 support with the same code
  - manually fix syntax in some support files
  - run python-modernize relative import fixer
  - run python-modernize print syntax fixer
  - run 2to3 exec syntax fixer
  - run python-modernize raise and except syntax fixers
  - fix dynamic creation of new methods
  - run python-modernize dict iteration fixers
  - fix Python version comparison
  - switch from deprecated os.path.walk to os.walk
  - import reload() on Python 3
  - run python-modernize dict.has_key fixer
  - fix StringIO import for test on Python 3
  - run 2to3 tuple parameters fixer
  - run python-modernize xrange fixer
  - update use of string functions
  - reverse type check of filename - 'file' is not a reliable type
  - run python-modernize reduce fixer
  - use six.integer_types to check for integers
  - misc fixes in sfepy.discrete
  - misc fixes in sfepy.discrete.fem.meshio
  - run python-modernize print syntax fixer on scripts and examples
  - run python-modernize relative import fixer on scripts, examples and tests
  - run python-modernize dict iteration fixers on scripts and examples
  - run python-modernize xrange fixer on scripts and examples
  - do not use six.iteritems() with Container subclasses
  - do not use relative imports in examples - fixes SystemError: Parent module
    '' not loaded, cannot perform relative import
  - change cmp= to key= when sorting lists
  - fix VariableTimeStepper.set_step()
  - fix exception instance not defined outside except block
  - implement point.__truediv__()
  - new enc(), dec() encoding utility functions
  - fix string IO in HDF5MeshIO
  - fix string IO in HDF5 (pytables) related functions
  - update test_install.py for Python 3
  - fix integer division
  - fix comparison of strings containing floats in test_units()
  - fix reporting of failed tests in test_install.py with Python 3 - new
    report_tests()
  - docs: update installation instructions for Python 3

- merge pull request #332 from 'rc/travis-ci', closes #321

  - automatic testing on Python 2.7 and 3.4 using Travis CI
  - new .travis.yml
  - update run_tests.py to return status
  - do not require DISPLAY in linear_elastic_probes.py example
  - allow failing of evp0 in tests/test_eigenvalue_solvers.py
    - update linear_elastic_probes.py example to run without vtk probes
  - update linear_elastic_mM.py example to regenerate coefficients - fixes race
    condition with several travis runs

- merge branch 'readme-rst'

  - show travis build status on github
  - rename README -> README.rst
  - update setup.py, sfepy/version.py for README.rst
  - README.rst: show travis status, update and fix text and links

- merge branch 'plot-cmesh'

  - fix plot_wireframe()
  - new plot_cmesh(), support **kwargs in plotting functions
  - script/plot_mesh.py: use plot_cmesh()

- merge pull request #333 from 'vlukes/hyperelast'

  - new term: NonsymElasticTerm - non-symmetric gradient
  - fix Term.get_approximation()
  - update hyperelastic terms - new "family data" implementation
  - new sym2nonsym() function in terms_op.c
  - new classes for homogenized coefficients: CoefNonSymNonSym and CoefNonSym
  - update "ev_integrate_mat" term - allow arbitrary shaped material

- merge branch 'remove-get-approximation'

  - remove FieldVariable.get_approximation(), Term.get_approximation()
  - update terms for no Term.get_approximation()

- merge pull request #334 from 'vlukes/fixdocs'

  - fix docstring of SDDotVolumeTerm, SDDivTerm, SDDivGradTerm, SDConvectTerm,
    SDGradDivStabilizationTerm, SDDiffusionTerm

- merge pull request #336 from 'heczis/issue_281_argparse', closes #281

  - use argparse instead of optparse in:

    - examples/linear_elasticity/its2D_interactive.py
    - examples/linear_elasticity/linear_elastic_interactive.py
    - examples/linear_elasticity/linear_viscoelastic.py
    - examples/linear_elasticity/modal_analysis.py
    - examples/homogenization/rs_correctors.py
    - examples/diffusion/laplace_shifted_periodic.py
    - examples/diffusion/time_poisson_interactive.py
    - examples/large_deformation/compare_elastic_materials.py
    - test_install.py
    - script/blockgen.py
    - script/gen_term_table.py
    - script/gen_mesh_prev.py
    - script/plot_mesh.py
    - script/plot_quadratures.py
    - script/plot_times.py
    - script/tile_periodic_mesh.py
    - script/save_basis.py
    - script/convert_mesh.py
    - script/cylindergen.py
    - script/extract_surface.py
    - script/plot_logs.py
    - script/gen_iga_patch.py
    - script/show_terms_use.py
    - script/gen_lobatto1d_c.py
    - script/sync_module_docs.py
    - script/gen_gallery.py
    - script/plot_condition_numbers.py
    - run_tests.py
    - extractor.py
    - probe.py
    - simple.py
    - schroedinger.py
    - homogen.py
    - postproc.py
    - phonon.py
    - shaper.py

  - fix import path
  - fix usage in examples/linear_elasticity/shell10x_cantilever_interactive.py

- miscellaneous updates:

  - support numeric prefixes in Quantity, Unit.get_prefix(), update prefixes
  - update script/show_authors.py for Python 3
  - fix docstring of get_local_ids(), add comments in mesh_build()
  - fix regression in FieldVariable.get_element_diameters()
  - add support for os.walk() keyword arguments to locate_files(),
    remove_files()
  - fix _gen_common_data() in tests/test_poly_spaces.py to permute connectivity

.. _2016.1-2016.2:

from 2016.1 to 2016.2
=====================

- merge pull request #309 from 'vlukes/splines'

  - bsplines: update to_ndarray(), fix draw() - evaluate curve/surface
    coordinates, if needed

- merge branch 'fix-gradient-items-ordering'

  - fix ordering of gradient items in evaluate_in_rc() - ordering corresponds
    to (n_coor, n_components, dim) as described in Field.evaluate_at()
    docstring
  - fix docstring of FieldVariable.evaluate_at()
  - change data shape of GradTerm values to correspond to Field.evaluate_at() -
    update GradTerm.get_fargs(), .get_eval_shape()
  - update test_field_gradient() to test fields with more than one component

    - test proportions of component gradients
    - update prepare_variable() so that components are multiples of the first
      one

- merge pull request #310 from 'vlukes/meshio_msh'

  - update meshio.py: support for msh file format (gmsh) - reading
  - MSH mesh format: update tests, add test meshes

- merge pull request #315 from 'vlukes/terms_cleanup'

  - rename 'd_diffusion_sa' term to 'd_sd_diffusion', remove unused functions
  - remove terms_acoustic from doc
  - make dw_lin_elastic_iso as the shortcut to dw_lin_elastic +
    stiffness_from_lame()
  - fix docstring of assemble_by_blocks()
  - update examples: replace dw_lin_elastic_iso
  - update tests and scripts: linear elastic terms
  - rename di_surface_moment to d_surface_moment
  - term table divided into: basic, sensitivity, large deformation, special
    terms

- merge pull request #316 from 'vlukes/parallel_homog'

  - new: parallel computation of homogenized coefficients
  - examples/phononic/band_gaps.py: no multiprocessing
  - update test_install.py - new test to check presence of lines in the output
  - new flush() method in OutputFilter class - needed for multiprocessing

- merge branch 'shell10x'

  - partial shell10x element implementation
  - new sfepy/discrete/structural/ for structural elements
  - fix docstring of IGField.get_data_shape()
  - replace 'plate' integration with 'custom'
  - update Field.from_conf() to scan sfepy/discrete/structural/
  - new sfepy/mechanics/shell10x.py - shell10x element implementation
    functions

    - new create_elastic_tensor(), create_transformation_matrix(),
      transform_asm_matrices(), create_local_bases(), create_rotation_ops(),
      create_strain_transform(), get_mapping_data(), get_dsg_strain(),
      create_strain_matrix(), add_eas_dofs(), rotate_elastic_tensor(),
      create_drl_transform(), lock_drilling_rotations()

  - new sfepy/discrete/structural/mappings.py - new Shell10XMapping
  - new sfepy/discrete/structural/fields.py - new Shell10XField
  - new sfepy/terms/terms_shells.py - new Shell10XTerm (dw_shell10x)
  - update ConcentratedPointLoadTerm.arg_shapes for general number of
    components
  - fix PhysicalQPs.get_shape() for no quadrature points (point integration
    terms)
  - fix Term.check_shapes() for scalar 'N' values
  - new Shell10XTerm.poly_space_base class attribute
  - update tests/test_term_call_modes.py to report shapes of all values
  - update tests/test_term_call_modes.py to test Shell10XTerm

    - update make_term_args(), Test
    - support custom integration and dim != tdim geometry
    - obey optional Term.poly_space_base attribute
    - use identity for 2D material matrix

  - new examples/linear_elasticity/shell10x_cantilever_interactive.py
  - update Term.from_desc() to pass integral to term constructor
  - new examples/linear_elasticity/shell10x_cantilever.py + test
  - add custom view for shell10x_cantilever.py example to script/gen_gallery.py

- merge branch 'mayavi-dataset-manager'

  - remove sfepy/postprocess/dataset_manager.py
  - use dataset_manager.py from mayavi

- merge branch 'no-symlinks'

  - remove scripts-common/
  - update setup.py to install main scripts into sfepy/script/
  - sfepy-run: update for no scripts-common/ symlinks, use explicit script
    names
  - docs: update for no posix only sfepy-run

- merge branch 'fix-windows-build' - closes #317, #318, #325

  - setup.py: import setuptools in setup_package() to find a C compiler on
    windows
  - fix __SDIR__ definition, new inline definition for windows in setup.py
    files
  - allow long shape in parse_shape()
  - remove inline directive for ravel_multi_index() on windows - fixes linker
    error

- materials:

  - merge pull request #311 from 'vlukes/sd_elastic'

    - remove unused variables in terms_[elastic, diffusion, basic].c
    - update d_sd_lin_elastic term: new, much faster implementation

  - fix .arg_shapes class attribute of DivGradTerm for no material
  - use special material for index in ScalarDotGradIScalarTerm - update
    .arg_shapes, .dw_fun()
  - fix arg_shapes in SDDotVolumeTerm

- docs:

  - remove no longer used terms_acoustic.rst
  - update support section
  - sync module index of developer guide with current sources

- miscellaneous updates:

  - fix mesh_get_centroids() for cells of lower topological dimension
  - fix VTKMeshIO.read() for cells of lower topological dimension than space -
    simplify vtk_inverse_cell_types
  - obey linearization kind in FEField.create_output() - allows adaptive
    linearization also for, e.g., Q1 fields
  - add common sources to dependencies of igac extension module

.. _2015.4-2016.1:

from 2015.4 to 2016.1
=====================

- merge pull request #307 from 'vlukes/mesh_generators'

  - fix gen_mesh_from_voxels()
  - new tests of gen_mesh_from_geom(), gen_tiled_mesh(), gen_mesh_from_voxels()

- merge branch 'auto-check-material-shapes'

  - implement general Term.check_shapes() - check term argument shapes at
    run-time
  - update terms to use generic variable size in .arg_shapes where appropriate

    - update IntegrateVolumeTerm, IntegrateSurfaceTerm, VolumeTerm,
      SurfaceTerm, VolumeSurfaceTerm, IntegrateMatTerm, SumNodalValuesTerm,
      GradTerm

  - remove .check_shapes() from all terms having it

    - remove it from DotProductVolumeTerm, BCNewtonTerm,
      VectorDotGradScalarTerm, VectorDotScalarTerm, LinearElasticIsotropicTerm,
      LinearPrestressTerm, LinearStrainFiberTerm, SurfaceTractionTLTerm,
      VolumeSurfaceTLTerm, ConcentratedPointLoadTerm
    - new .arg_shapes class attribute in DotProductSurfaceTerm,
      ConcentratedPointLoadTerm

  - update LinearPointSpringTerm for new .arg_shapes class attribute

    - change special material argument to a single float
    - update tests/test_elasticity_small_strain.py

  - update get_arg_kinds() to distinguish 'ts' argument, new _match_ts
  - update tests/test_term_call_modes.py for TimeStepper ('ts') term argument -
    update make_term_args()
  - new/fill-in .arg_shapes class attributes in time history terms

    - update BiotTHTerm, BiotETHTerm, DotSProductVolumeOperatorWTHTerm,
      DotSProductVolumeOperatorWETHTerm, LinearElasticTHTerm,
      LinearElasticETHTerm, CauchyStressTHTerm, CauchyStressETHTerm

  - change 'N' value to 1 in _parse_scalar_shape() in make_term_args() - fix
    for time history terms

- merge pull request #306 from 'vlukes/fix-gen_mesh_from_geom'

  - fix to_poly_file() in geom_tools.py
  - fix gen_mesh_from_geom(), remove gen_mesh_from_poly()

- merge branch 'remove-ts-explicit'

  - remove make_explicit_step(), ExplicitTimeSteppingSolver
  - remove MassOperator and sfepy/discrete/mass_operator.py
  - update time_poisson_explicit.py to use ts.simple

- merge pull request #304 from 'vlukes/splines'

  - new documentation to SplineBox and SplineRegion2D classes
  - update bspline.py, new SplineRegion2D in splinebox.py

- merge branch 'no-fea'

  - replace Interpolant by PolySpace in GeometryElement

    - GeometryElement.interp -> .poly_space
    - update FEDomain.__init__() and affected code

  - move set_mesh_coors() into sfepy/discrete/fem/fields_base.py - update
    Problem.set_mesh_coors()
  - move Approximation into FEField and subclasses, part 1

    - update volume fields
    - prepare for volume-only PolySpace in fields
    - remove imports of fea.py
    - sfepy/discrete/fem/fea.py:

      - move eval_nodal_coors(), _interp_to_faces() into
        sfepy/discrete/fem/fields_base.py
      - remove Interpolant
      - remove Approximation.eval_extra_coor(), .get_connectivity(),
        .get_poly_space()
      - move into FEField:

        - Approximation.clear_qp_base(), .get_qp(), .get_base()
        - Approximation._create_bqp(), .create_bqp() into new
          FEField.create_bqp()
        - Approximation.describe_geometry() into FEField.create_mapping()

    - new FEField attributes:

      - from Approximation .surface_data, .point_data, .ori, .efaces, .econn
      - .poly_space

    - update FEField.__init__(), ._setup_esurface(), .setup_coors(),
      .get_data_shape(), .linearize(), .interp_to_qp()
    - replace Interpolant by PolySpace in VolumeField
    - update VolumeField._create_interpolant(), ._init_econn(),
      ._setup_vertex_dofs(), .setup_extra_data(), .get_econn()
    - remove VolumeField._setup_approximations()
    - rename VolumeField._setup_surface_data(), ._setup_point_data() to
      .setup_surface_data(), .setup_point_data(), merge with
      Approximation.setup_surface_data(), .setup_point_data()
    - update H1NodalMixin._setup_facet_orientations(), ._setup_facet_dofs(),
      ._setup_bubble_dofs(), .create_basis_context()
    - update H1NodalVolumeField.interp_v_vals_to_n_vals()
    - update H1HierarchicVolumeField._init_econn(),
      ._setup_facet_orientations(), ._setup_facet_dofs(),
      ._setup_bubble_dofs(), .create_basis_context()
    - update eval_in_els_and_qp(), create_expression_output()

  - move Approximation into FEField and subclasses, part 2

    - update surface integration/mappings for volume-only PolySpace in fields
    - update FEField.get_data_shape(), .create_bqp(), .create_mapping()
    - update FEMapping.__init__() - new .indices attribute
    - new SurfaceMapping.set_basis_indices(), .get_base()

  - move Approximation into FEField and subclasses, part 3

    - update FEField.create_mesh(), VolumeField.average_qp_to_vertices()
    - update compute_nodal_normals()
    - update FieldVariable.get_element_diameters()
    - update describe_geometry() in membranes.py

  - move Approximation into FEField and subclasses, part 4

    - remove H1DiscontinuousField._setup_approximations()
    - update H1DiscontinuousField._setup_global_base()

  - move Approximation into FEField and subclasses, part 5

    - update surface fields
    - fix .surface_data, .point_data initialization in FEField.__init__()
    - update FEField.get_qp(), .create_mapping()
    - update SurfaceField._create_interpolant(), .setup_extra_data(),
      ._init_econn(), ._setup_vertex_dofs(), .get_econn(),
      .average_qp_to_vertices()
    - remove SurfaceField._setup_approximations()

  - remove sfepy/discrete/fem/fea.py
  - update FieldVariable.get_approximation() to return Field
  - new FEField.get_connectivity() convenience alias
  - update sfepy/parallel/ code for no Approximation

    - update create_task_dof_maps(), distribute_field_dofs(),
      distribute_fields_dofs(), get_local_ordering(), plot_partitioning(),
      plot_local_dofs()

  - script/save_basis.py: update save_basis_on_mesh() for no Approximation
  - update tests/test_poly_spaces.py for no Approximation

- merge branch 'active-fibres-update'

  - update HyperElasticBase to pass kwargs to stress and tangent modulus
    functions - update HyperElasticBase.compute_stress(), .compute_tan_mod()
  - new create_omega(), compute_fibre_strain()
  - update FibresActiveTLTerm - move fibre_function() into the class, cache data

    - remove fibre_function()
    - new _setdefault_fibre_data()
    - update FibresActiveTLTerm.get_fargs(), .stress_function(),
      tan_mod_function()

- merge branch 'parallel-pc-fieldsplit'

  - support 'fieldsplit' preconditioner in PETScKrylovSolver

    - new .set_field_split()
    - update .__init__(), .create_ksp()

  - setup 'fieldsplit' preconditioner in biot_parallel_interactive.py example
  - update docstring of biot_parallel_interactive.py example

- docs:

  - sync IGA section with current state
  - document refinement_level configuration option in users guide
  - stop omitting time_poisson_explicit.py in script/gen_gallery.py
  - sync module index of developer guide with current sources

- scripts:

  - new script/show_mesh_info.py
  - script/convert_mesh.py: update --list option to list also readable formats

    - rename & update output_writable_meshes() -> output_mesh_formats()

- examples and tests:

  - parallel examples: fix race condition when output directory does not exist
  - adjust final time in time_poisson_explicit.py

- miscellaneous updates:

  - merge branch 'fix-mat-by-region-for-surfaces2'

    - update Region.get_cells() to obey parent region
    - fix Region.get_cell_indices() for non-disjoint, non-subset cells

  - remove useless FEField._create_interpolant()
  - remove obsolete SurfaceLaplaceLayerTerm, SurfaceCoupleLayerTerm
  - update Problem.from_conf() to support refinement_level configuration option
  - fix argument name in docstring of FieldVariable.evaluate()
  - fix NEUMeshIO.read() for empty lines and 3_8 cells
  - fix Region.get_facet_indices() for safe numpy casting on windows10

.. _2015.3-2015.4:

from 2015.3 to 2015.4
=====================

- basic support for restart files

  - merge branch restart-files
  - simple.py: new --save-restart, --load-restart options
  - new Problem.get_restart_filename(), .save_restart(), .load_restart() -
    update .init_time() to initialize new ._restart_filenames attribute
  - update Problem.__init__(): update default conf to have options attribute
  - fix Variables.set_data() to use step argument
  - update SimpleTimeSteppingSolver.__call__() to support restart files
  - update Problem for restarting stationary problems - update .reset() to
    initialize ._restart_filenames attribute
  - update StationarySolver.__call__() to support restart files
  - new TimeStepper.get_state(), .set_state()
  - new VariableTimeStepper.get_state(), .set_state()
  - update VariableTimeStepper to have current state stored in .times, .dts

    - new .advance(), .iter_from_current()
    - update .set_step(), .__iter__()

  - update AdaptiveTimeSteppingSolver.__call__() to support restart files
  - update SimpleTimeSteppingSolver.__call__() for variables with history
  - update TimeStepper.iter_from(), new .advance()
  - check step in Variable.set_data()
  - docs: update user's guide - introduce restart files

- linear combination boundary conditions:

  - improve docstrings of MRLCBCOperator, ShiftedPeriodicOperator
  - update ShiftedPeriodicOperator.__init__() for pyflakes
  - support general linear combination of DOFs in a node

    - merge branch nodal-lcbcs
    - update LCBCOperators.make_global_operator() for rhs without column
      variable
    - new NodalLCOperator - general linear combination of DOFs in a node
    - new examples/linear_elasticity/nodal_lcbcs.py + test

  - fix IntegralMeanValueOperator for several DOFs per node

- examples and tests:

  - merge branch example-balloon

    - new meshes/3d/unit_ball.mesh
    - new examples/large_deformation/balloon.py + test

- miscellaneous updates:

  - do not check facet-only meshes in FEDomain.fix_element_orientation() -
    fixes segfault
  - remove unused doc/images/sfepy_gui.png
  - fix Mesh.copy() to provide default name
  - fix bubble DOFs setup in H1DiscontinuousField._setup_global_base() -
    initialize .bubble_remap attribute
  - clean up sfepy/postprocess/time_history.py
  - fix extract_time_history() for no element groups
  - docs: fix Green strain definitions
  - allow None as problem argument in MiniAppBase.__init__()
  - lib.lapack import bug fix (merge pull request #301 from rexfuzzle/master)
  - fix output order of scipy.linalg.lapack functions in ScipySGEigenvalueSolver
  - fix version comparison in dets_fast()
  - set SYMPY_MIN_VERSION to 0.7.3 in sfepy/version.py (merge branch
    sympy-lcbc-compat)
  - fix describe_geometry() to initialize base functions in membrane_geo
  - merge branch fix-mat-by-region-for-surfaces

    - allow non-cell regions in ConstantFunctionByRegion.get_constants()
    - update Region.get_cell_indices() to allow cells to be sutperset of region
      cells

  - fix access to mat_id in GenericFileSource.create_source(), .get_mat_id()

.. _2015.2-2015.3:

from 2015.2 to 2015.3
=====================

- preliminary support for parallel computing

  - merge branch parallel
  - allow constructing empty regions

    - update Domain.create_region(), Region.setup_from_highest(), .finalize()
    - new allow_empty argument

  - new sfepy/parallel/__init__.py
  - new sfepy/parallel/setup.py, update sfepy/setup.py
  - new sfepy/parallel/parallel.py - start PETSc-based parallelization

    - new get_inter_facets(), create_task_dof_maps(), distribute_field_dofs(),
      get_local_ordering(), get_sizes(), expand_dofs(), create_petsc_matrix(),
      apply_ebc_to_matrix(), assemble_to_petsc()

  - use cmesh.tdim in get_inter_facets()
  - new create_prealloc_data()
  - new partition_mesh()
  - new petsc_call() linear solver decorator
  - update PETScKrylovSolver for parallel use

    - allow passing in PETSc matrices and vectors
    - new sub_precond option
    - update .__init__() - new comm argument, setup converged_reasons there
    - update .create_ksp()
    - remove .set_matrix(), new .create_petsc_matrix()
    - update .__call__() - new comm argument, use petsc_call() decorator,
      return PETSc solution vector for PETSc right-hand side vector
    - update docstring

  - new view_petsc_local()
  - new create_local_petsc_vector()
  - new create_gather_scatter(), create_gather_to_zero()
  - new verify_task_dof_maps()
  - new is_matrix argument in Problem.time_update(), .update_equations()
  - update PETScKrylovSolver.__call__() to output actual solver and
    preconditioner
  - new distribute_fields_dofs() - support multiple fields
  - new get_composite_sizes()
  - split assemble_to_petsc() - new assemble_rhs_to_petsc(),
    assemble_mtx_to_petsc()
  - remove debug() call in Equations.get_graph_conns()
  - support non-reduced (full size) system assembling in Equations and
    Variables

    - update Equations.time_update(), .get_graph_conns(),
      .create_matrix_graph() - new active_only argument
    - update create_adof_conns(), Variables.equation_mapping() - new
      active_only argument

  - update Problem for non-reduced (full size) system assembling

    - new .active_only attribute
    - update .__init__(), .update_equations() - new active_only argument

  - new PETScNonlinearSolver
  - new sfepy/parallel/evaluate.py - new PETScParallelEvaluator
  - new setup_composite_dofs()
  - new create_petsc_system()
  - update setup.py - new petsc4py and pymetis version checks - update
    sfepy/version.py
  - filter-out -h, --help from sys.argv options passed to petsc4py

    - -h, --help is avaliable for user options, -help can be used to show PETSc
      options

  - new sfepy/parallel/plot_parallel_dofs.py - new mark_subdomains(),
    label_dofs(), plot_partitioning(), plot_local_dofs()

  - new examples/diffusion/poisson_parallel_interactive.py
  - new examples/multi_physics/biot_parallel_interactive.py
  - update setup.py - new mpi4py version check - update sfepy/version.py
  - docs:

    - sync module index of developer guide with current sources
    - parallel examples: add mpi4py requirement
    - update installation requirements
    - update user's guide - add basic parallel problem solving description -
      update developer guide, doc/index.rst

  - update test_install.py to test parallel examples

- (mostly) fix finding of reference coordinates

  - merge branch fix-find-ref-coors - (almost) solves #285

    - TODO: make get_xi_tensor() robust w.r.t. multiple solutions, prefer those
      inside a cell

  - rename get_ref_coors() -> get_ref_coors_convex(), refc_find_ref_coors() ->
    refc_find_ref_coors_convex(), find_ref_coors() -> find_ref_coors_convex()
  - new refc_find_ref_coors() C function, new find_ref_coors()
  - new get_potential_cells(), get_ref_coors_general()
  - remove strategy argument of get_ref_coors_convex(), use extrapolate flag
  - new get_ref_coors() - support 'general' and 'convex' strategies
  - update .evaluate_at() and its calls for new get_ref_coors()

    - update H1HierarchicVolumeField.evaluate_at()
    - update H1NodalMixin.evaluate_at()
    - update FieldVariable.evaluate_at(), .set_from_other()
    - update Probe.probe()

  - fix comment in get_xi_tensor()
  - new tests/test_ref_coors.py + cross3d.mesh

- allow field gradient evaluation in arbitrary points

  - merge branch evaluate-gradient
  - update evaluate_in_rc() for gradients
  - update .evaluate_at() and its calls - new mode argument, support gradients

    - update H1HierarchicVolumeField.evaluate_at()
    - update H1NodalMixin.evaluate_at()
    - update FieldVariable.evaluate_at(), .set_from_other()
    - update Probe.probe()

  - update test_invariance_qp() in tests/test_mesh_interp.py - test field
    gradients, increase field approximation order to 2
  - fix eval_lagrange_simplex() - initialize properly gradient output to zeros
  - make vector field gradient in evaluate_in_rc() compatible with ev_grad
    term - define ib type
  - update test_invariance_qp() to test several meshes with different cell
    types
  - update test_invariance_qp() for new prepare_variable(), update reporting
  - new test_field_gradient() in tests/test_mesh_interp.py

- unify evaluation of basis functions using basis-specific contexts, part 1

  - merge branch basis-context
  - hide basis-specific arguments in lagrange.c into new LagrangeContext struct

    - update get_barycentric_coors(), get_xi_simplex(), get_xi_tensor(),
      eval_lagrange_simplex(), eval_lagrange_tensor_product() C functions
    - update get_barycentric_coors(), eval_lagrange_simplex(),
      eval_lagrange_tensor_product(), evaluate_in_rc() cython functions

  - update evaluate_in_rc() for merged 'evaluate-gradient' branch
  - new print_context_lagrange() C function
  - update find_ref_coors_convex, find_ref_coors() (cython) for
    LagrangeContext - update refc_find_ref_coors_convex(),
    refc_find_ref_coors() C functions
  - new get_xi_dist() C function, update LagrangeContext - add tdim attribute
  - new CLagrangeContext cython class, wraps LagrangeContext C struct

    - update LagrangeContext cython definition
    - new CLagrangeContext.__cinit__(), .__dealloc__(), .__str__(), .cprint()

  - update refc_find_ref_coors_convex(), refc_find_ref_coors() - new qp_eps
    argument

    - new abstract BasisContext C struct
    - remove _get_xi_dist()
    - no more dependence on lagrange.h

  - new abstract CBasisContext cython class
  - update find_ref_coors_convex, find_ref_coors() (cython) for CBasisContext -
    and for updated refc_find_ref_coors_convex(), refc_find_ref_coors()
  - new H1NodalMixin.create_basis_context()
  - update get_ref_coors_convex(), get_ref_coors_general() - use
    Field.create_basis_context()
  - new H1HierarchicVolumeField.create_basis_context()

- move common extension modules to sfepy/discrete/common/extmods

  - merge branch 'common-extmods'
  - new sfepy/discrete/common/extmods/ - move common extension modules there

    - move sfepy/discrete/fem/extmods/{_fmfield.pxd, _fmfield.pyx,
      _geommech.pxd, _geommech.pyx, assemble.pyx, cmesh.pxd, cmesh.pyx,
      common.h, common_python.c, crefcoors.pyx, fmfield.c}

  - new sfepy/discrete/common/extmods/__init__.py
  - update setup.py files, new sfepy/discrete/common/extmods/setup.py
  - update imports

- unify evaluation of basis functions using basis-specific contexts, part 2

  - merge branch basis-context2
  - update LagrangeContext

    - add .eval_basis, .order, .is_bubble attributes
    - new .iel, .geo_ctx, .mesh_coors, .mesh_conn, .n_cell, .n_cp attributes
    - change .bc attribute type
    - new .mbfg attribute
    - update print_context_lagrange()

  - update CLagrangeContext

    - update for current LagrangeContext
    - update .__cinit__(), new .is_bubble property
    - new .evaluate()
    - update .__cinit__()
    - new ._geo_ctx, .mesh_coors, .mesh_conn, .base1d attributes
    - new .iel, .geo_ctx properties
    - new .mbfg attribute

  - new LagrangePolySpace

    - new .create_context()
    - new ._eval_base(), update derived classes

      - update LagrangeSimplexPolySpace:

        - update .__init__() - new init_context argument
        - remove ._eval_base()

      - update LagrangeSimplexBPolySpace:

        - update .__init__() - new init_context argument
        - update .create_context() - support **kwargs
        - remove ._eval_base()

      - update LagrangeTensorProductPolySpace:

        - update .__init__() - new init_context argument
        - remove ._eval_base()

  - update LagrangeSimplexPolySpace, LagrangeTensorProductPolySpace

    - inherit from LagrangePolySpace
    - new LagrangeSimplexBPolySpace.create_context()

  - new fmf_set_qp() C function
  - new eval_basis_lagrange() C function, update eval_lagrange_simplex() -
    support only a single point
  - update get_xi_dist() for current LagrangeContext - update get_xi_simplex(),
    get_xi_tensor()
  - update H1NodalMixin.create_basis_context()
  - update BasisContext - new .eval_basis, .iel attributes
  - move evaluate_in_rc() from bases.pyx to crefcoors.pyx, update for
    BasisContext - great simplification by using abstract
    BasisContext.eval_basis()

  - remove cython functions replaced by CLagrangeContext.evaluate()

    - remove get_barycentric_coors(), eval_lagrange_simplex(),
      eval_lagrange_tensor_product()

  - remove unused definitions in bases.pyx
  - move global_interp.py into sfepy/discrete/common/
  - move H1NodalMixin.evaluate_at() -> Field.evaluate_at(), remove
    H1HierarchicVolumeField.evaluate_at()
  - update imports
  - improve error message in fmf_copy()
  - allow higher-order cell-vertex connectivities in CMesh.from_data()
  - move buffer for cell coordinates into context

    - update BasisContext, LagrangeContext - new e_coors_max attribute
    - update CLagrangeContext - new e_coors_max attribute, update .__cinit__()

  - split meaning of .iel into .iel and new .is_dx attribute of BasisContext

    - update BasisContext, LagrangeContext - new .is_dx attribute
    - update CLagrangeContext.__cinit__()
    - update evaluate_in_rc()
    - update refc_find_ref_coors_convex(), refc_find_ref_coors() C functions -
      set ctx->is_dx and ctx->iel
    - update eval_basis_lagrange()

  - update evaluate_in_rc() - add basic shape check of output buffer
  - improve handling of get_xi_dist() failure in refc_find_ref_coors() -
    initialize imin, imin and d_min are in sync
  - update eval_basis_lagrange() for bubble functions. update
    eval_lagrange_simplex()
  - new NURBSContext C struct

    - new print_context_nurbs(), get_xi_dist(), eval_basis_nurbs()
    - update eval_bspline_basis_tp(), eval_nurbs_basis_tp() - new is_dx
      argument

  - new CNURBSContext cython class, wraps NURBSContext C struct

    - new CNURBSContext.__cinit__(), .__dealloc__(), .__str__(), .cprint(),
      .evaluate()

  - update eval_mapping_data_in_qp(), eval_variable_in_qp(), eval_in_tp_coors()
    for is_dx argument
  - update IGDomain - new .eval_mesh attribute, update .__init__()
  - new IGField.create_mesh(), .create_eval_mesh(), .create_basis_context()
  - new Field.create_eval_mesh()
  - fix Field.evaluate_at() - do not use FEM-specific attributes
  - update get_ref_coors_general() for field.create_eval_mesh()
  - update sfepy/discrete/common/__init__.py - import Field
  - rename & update test_ref_coors() -> test_ref_coors_fem()

    - remove field creation from Test.from_conf()
    - use basis context for basis evaluation
    - fix typo

  - new test_ref_coors_iga() in tests/test_ref_coors.py
  - update IGField to have approx_order attribute (make projections work)
  - new test_projection_iga_fem() in tests/test_projections.py

- merge pull request #298 from vlukes/el_eval

  - new evaluation mode 'el_eval', 'el' mode removed
  - update equations.evaluate()

- merge pull request #300 from vlukes/mat_optim

  - new example: material identification - multiscale analysis
  - new example: material identification - mesh
  - new tutorial: material identification

- merge branch test-conditions

  - move tests/test_ebcs.py -> tests/test_conditions.py
  - tests/test_conditions.py: remove configuration
  - new test_ebcs() in tests/test_conditions.py, update Test.from_conf()
  - new test_epbcs() in tests/test_conditions.py
  - update test_save_ebc() in tests/test_conditions.py
  - move checking of applied conditions from test_ebcs() into new check_vec()
  - new test_ics() in tests/test_conditions.py

- merge branch how-to-contribute

  - update generic installation instructions, add two section labels
  - update how to contribute instructions in developer guide
  - remove gitwash section from developer guide
  - remove doc/dev/gitwash/
  - update script/sync_module_docs.py
  - update developer guide with Vladimir's suggestions

- postprocessing and visualization:

  - fix VTK version issue related to SetSource()
  - fix VTK version issue related to Update()
  - fix plot_dofs option in save_basis_on_mesh() for no element groups
  - remove show argument of plot functions, use single code for 2D and 3D

    - support 1D where applicable
    - update plot_control_mesh()
    - update plot_wireframe(), plot_entities(), label_global_entities(),
      label_local_entities()
    - update plot_points(), plot_mesh(), plot_global_dofs(), plot_local_dofs(),
      plot_nodes()
    - update plot_geometry(), plot_edges(), plot_faces()
    - update plot_weighted_points(), plot_quadrature()
    - update import in script/plot_mesh.py

  - fix tetrahedralize_vtk_mesh()
  - use tight bounding box in figure.savefig() calls in plot_parallel_dofs.py
  - postproc.py: new --colormap option

    - new set_colormap()
    - update Viewer.build_mlab_pipeline(), .call_mlab()

  - allow user-defined file sources in place of filename in
    create_file_source()

- probes:

  - update Probe.__call__() to support keyword arguments
  - update Probe.probe() - new ret_points argument, update docstring
  - update VTK Probe.__call__() - new ret_points argument, update docstring,
    new .dim attribute

- materials:

  - remove unused Material.set_data_from_variable()
  - use expression materials when copy_materials=True in
    Problem.create_evaluable()

    - fixes evaluation of terms with materials not present in Problem.equations
    - initialize Problem.conf.materials in Problem.__init__()

  - remove unused Materials.semideep_copy()

- regions:

  - fix Region.get_facet_indices() for inner cell "corners"
  - fix "Cannot cast array data ..." (Windows 8.1, Numpy 1.8.2)

- terms:

  - new PiezoStressTerm (ev_piezo_stress)
  - new AdvectDivFreeTerm (dw_advect_div_free)

- setup:

  - remove unused IPYTHON_MIN_VERSION
  - update package_check() - new show_only argument, update messages and
    docstring
  - update setup.py to show names and versions of dependencies at the end - new
    check_versions()
  - update Clean.run() build helper to clean *.pyd files (dynamic libs on
    windows)
  - increase SYMPY_MIN_VERSION to 0.7.2
  - new scikits.umfpack and pysparse version checks - update sfepy/version.py

- docs:

  - update installation requirements - add scikit-umfpack
  - fix label of theoretical background section
  - fix indentation, include links.inc in user's guide

- examples and tests:

  - fix time_poisson_interactive.py for new copy_materials semantics
  - time_poisson_interactive.py: update circle probe parameters in gen_lines()
  - fix --show option in poisson_parallel_interactive.py example - update for
    c985650d356bcc1eb3608f4144fc84299e67c458
  - clean up piezo_elasticity.py example, use short syntax
  - update piezo_elasticity.py example to compute stresses - new post_process()
  - fix its2D_interactive.py for new copy_materials semantics, update
    stress_strain(), nodal_stress()
  - modal_analysis.py:

    - new --axis option
    - support optional positional argument for user-provided mesh
    - fix --show option with user-provided mesh
    - fix --n-eigs help message
    - fix paths in help/docstring
    - update --bc-kind option - new fixed mode, rename clamped -> cantilever
    - new --ignore option
    - improve reporting results
    - update test_install.py for current modal_analysis.py
    - fix line lengths

  - new examples/diffusion/time_advection_diffusion.py + test

- miscellaneous updates:

  - remove whitespaces in linalg/utils.py
  - replace dets_fast based on lapack_lite by numpy.linalg.det for numpy
    version >= 1.8
  - fix expand2d()
  - remove duplicate line in get_ref_coors()
  - fix docstring of DiffusionCoupling term
  - fix sympy.zeros() call for sympy 0.7.6
  - fix FieldVariable.setup_initial_conditions() for multiple conditions
  - fix test_eigenvalue_solvers() to obey Test.can_fail
  - clean up tests/test_msm_symbolic.py
  - move tests/sympy_operators.py -> sfepy/linalg/sympy_operators.py
  - new HDF5MeshIO.read_dimension()
  - fix docstring of configure_output()
  - fix exception message in FEDomain.fix_element_orientation()
  - update PysparseEigenvalueSolver to support extra options
  - script/show_authors.py: merge by names, show commit counts

.. _2015.1-2015.2:

from 2015.1 to 2015.2
=====================

- update time stepping solvers for interactive use

  - merge branch ts-interactive
  - rename Problem.set_solvers() -> .set_conf_solvers()
  - rename Problem.set_solvers_instances() -> .set_solver()

    - work only with nonlinear solver
    - remove .solvers, new .nls (the nonlinear solver instance)
    - rename .get_solvers() -> .get_solver()

  - update Problem.is_linear(), .set_linear() to use nonlinear solver instance
  - update ScipyDirect to presolve in new .presolve()

    - new LinearSolver.presolve()
    - remove {ScipyDirect, SchurGeneralized, MultiProblem}._presolve()

  - update Problem.init_solvers() to only create instances, new .try_presolve()
  - update Problem.time_update() to reuse existing conditions by default
  - make functions argument optional in Equations.setup_initial_conditions()
  - remove Problem.setup_ic(), new .set_ics(), .setup_ics()

    - new Problem.ics attribute - update .copy(), .clear_equations()
    - update get_initial_state()

  - make generator from .__call__() of time stepping solvers - yield step,
    time, state
  - move initial step initialization of time stepping solvers to new
    .init_time()

    - update StationarySolver, EquationSequenceSolver, SimpleTimeSteppingSolver,
      new .init_time() - call Problem.init_solvers() there
    - update make_implicit_step() - use Problem.try_presolve()
    - update make_explicit_step()
    - new TimeSteppingSolver.init_time()

  - updates for time stepping solvers as generators and .init_time()
  - move project_by_component() from example into sfepy code
  - new examples/diffusion/time_poisson_interactive.py

- solvers:

  - fix PETSc Krylov solvers for reusing with different matrices

    - new PETScKrylovSolver.create_ksp(), update .set_matrix(), .__call__()
    - update PETScParallelKrylovSolver.__call__()

- reorganize examples:

  - merge branch reorganize-examples
  - move Biot, piezo-, thermo-elasticity examples into examples/multi_physics/
  - move compare_elastic_materials.py into examples/large_deformation/, remove
    el3.mesh
  - clean up compare_elastic_materials.py, generate mesh, use short syntax,
    make it an executable script
  - move rs_correctors.py into examples/homogenization/, move osteonT1_11.mesh
    into meshes/2d/special/
  - clean up rs_correctors.py, use short syntax, add docstring, make it an
    executable script, change output to current directory
  - move laplace_shifted_periodic.py into examples/diffusion/
  - move linear_elasticity.py into examples/linear_elasticity/ - rename to
    linear_elastic_interactive.py
  - move modal_analysis.py into examples/linear_elasticity/
  - remove examples/standalone/interactive/__init__.py
  - move live_plot.py into examples/miscellaneous/
  - move thermal_electric.py into examples/multi_physics/
  - remove examples/standalone/__init__.py
  - update tests for moved examples
  - update test_install.py for moved examples
  - script/gen_gallery.py: update for moved examples
  - docs: update tutorial for moved examples

- remove element groups - major code simplification

  - merge branch no-ig
  - update CMesh and Mesh:

    - new CMesh.from_data(), new .vertex_groups attribute
    - start using CMesh internally in Mesh to store mesh data, no element
      groups - update Mesh._set_data(), ._set_shape_info()
    - remove CMesh.from_mesh(), .get_from_cell_group(), .get_igs()
    - new CMesh.create_new()
    - remove unused Mesh.setup_done attribute
    - remove unused make_point_cells()
    - remove unused write_bb()
    - rename Mesh._set_data() -> ._set_io_data(), update MeditMeshIO.read()
    - remove Mesh.mat_ids, .ngroups, make .coors @property alias to .cmesh.coors

      - update ._set_io_data()
      - new .coors() property

    - new Mesh.get_conn()
    - new Mesh._get_io_data()
    - update Mesh.write() - remove coors, igs arguments
    - update Mesh.from_data() for no element groups
    - remove unused Mesh.get_element_coors()
    - update Mesh.__init__() for CMesh

      - new cmesh argument, remove filename, prefix_dir, kwargs arguments
      - new Mesh._collect_descs()
      - update Mesh docstring

    - update Mesh.from_region(), .copy() for CMesh.create_new()

      - do not add facet entities in Mesh.from_region()

    - remove unused Mesh._append_region_faces(), .localize()
    - remove unused Mesh.from_surface()
    - update merge_mesh(), Mesh.__add__(), .create_conn_graph() for no element
      groups
    - remove unused get_min_edge_size()
    - remove unused Mesh.explode_groups()
    - remove unused make_inverse_connectivity()

  - mesh generators:

    - update tiled_mesh1d(), gen_tiled_mesh() for no element groups
    - update gen_extended_block_mesh() for no element groups

  - update MeshIO and subclasses:

    - do not split connectivities by mat_id in MeditMeshIO.read()
    - update MeditMeshIO.write(), VTKMeshIO.write() for Mesh._get_io_data()
    - do not split connectivities by mat_id in VTKMeshIO.read()
    - clean up tests/test_meshio.py
    - new split_conns_mat_ids(), update MeditMeshIO.read()
    - update MeshIO and remaining subclasses for current Mesh

      - do not split connectivities by mat_id in *.read()
      - use Mesh._get_io_data() in *.write()
      - update MeshIO docstring
      - update TetgenMeshIO.read()
      - update ComsolMeshIO.read(), .write()
      - update HDF5MeshIO.read(), .write()
      - update MEDMeshIO.read()
      - update Mesh3DMeshIO.read()
      - update mesh_from_groups()
      - update BDFMeshIO.read(), .write()
      - update NEUMeshIO.read()

  - remove unused meshio functions related to element groups

    - remove sort_by_mat_id(), sort_by_mat_id2(), split_by_mat_id(),
      join_conn_groups()

  - update Domain, FEDomain:

    - update FEDomain__init__() for CMesh in Mesh and no element groups

      - remove .setup_groups(), .iter_groups(), .get_cell_offsets()
      - remove .groups, .mat_ids_to_i_gs, .cell_offsets attributes

    - update orient_elements() C function for C Mesh, update cmesh.pyx
    - update FEDomain.fix_element_orientation() for CMesh
    - update FEDomain__init__() - set self.shape.n_el
    - update region_leaf() for no element groups, use CMesh

      - update create_bnf(): remove E_CI2, rename E_CI1 -> E_CI
      - update tests/test_regions.py

    - fix element orientation in FEDomain.__init__() before setting local
      entities
    - new FEDomain.get_conn(), update .__init__() to disallow multiple cell
      kinds
    - update remaining FEDomain methods for no element groups - update
      .get_element_diameters(), .create_surface_group(), .refine()
    - update refine_*() for no element groups

      - update refine_2_3(), refine_2_4(), refine_3_4(), refine_3_8()
      - update FEDomain.refine() for current refine_*()

    - update Domain.save_regions_as_groups() for no element groups
    - update FEDomain.get_evaluate_cache() for CMesh

  - update Region:

    - update Region for no element groups

      - update .__init__(), .set_kind(), .update_shape(), .get_entities(),
        .get_cells(), .get_facet_indices(), .setup_mirror_region(),
        .get_mirror_region(), .get_n_cells()
      - remove .igs property, ._igs, .ig_map, .ig_map_i attributes
      - remove .get_vertices(), .get_edges(), .get_faces(), .get_facets(),
        .iter_cells(), .get_cell_offsets()

    - remove unused _try_delete()
    - fix Region.get_facet_indices() to setup required connectivity
    - new Region.get_cell_indices()
    - remove unused Region.get_vertices_of_cells()

  - fields:

    - update VolumeField._check_region(), ._setup_geometry() for no element
      groups
    - remove unused FEField._setup_approximations() (re-implemented in
      subclasses)
    - update FEField and subclasses for no element groups

      - setup of global basis of volume nodal fields works
      - remove FEField.igs, .aps attributes
      - new FEField.ap attribute
      - update FEField._setup_global_base(), .setup_coors()
      - update VolumeField._setup_approximations(), ._init_econn(),
        ._setup_vertex_dofs()
      - update H1NodalMixin._setup_facet_orientations(), ._setup_facet_dofs(),
        ._setup_bubble_dofs()
      - update Field.get_mapping()
      - update FEField._setup_esurface(), ._get_facet_dofs(), .get_data_shape(),
        .interp_to_qp(), .create_mapping()
      - rename & update FEField.get_dofs_in_region_group() ->
        .get_dofs_in_region()
      - update VolumeField._setup_surface_data(), ._setup_point_data(),
        .get_econn(), .average_qp_to_vertices()
      - update H1HierarchicVolumeField._init_econn(), ._setup_facet_dofs(),
        ._setup_bubble_dofs(), .set_dofs()
      - update H1NodalMixin.set_dofs()

    - remove Field.get_dofs_in_region()
    - remove facet_desc argument of H1HierarchicVolumeField._setup_facet_dofs()
    - update FEField linearization functions for no element groups

      - update eval_in_els_and_qp()
      - update get_eval_expression(), create_expression_output()
      - update FEField.linearize()

    - update SurfaceField for no element groups

      - update ._check_region(), ._setup_geometry(), ._setup_approximations(),
        .setup_extra_data(), ._init_econn(), ._setup_vertex_dofs(),
        ._setup_bubble_dofs(), .get_econn(), .average_qp_to_vertices()

    - update H1DiscontinuousField for no element groups - update
      ._setup_approximations(), ._setup_global_base()
    - update FEField.create_mesh() for no element groups
    - update H1NodalVolumeField.interp_v_vals_to_n_vals() for no element groups
    - update get_facet_dof_permutations() for no element groups
    - update find_ref_coors(), evaluate_in_rc() for no element groups
    - update H1NodalMixin.evaluate_at() for no element groups - update
      get_ref_coors()

  - IGA:

    - update IGDomain.__init__() for no element groups
    - update IGField for no element groups

      - update .__init__(), .get_econn(), .get_data_shape(), .set_dofs(),
        .create_mapping()
      - rename and update .get_dofs_in_region_group() -> .get_dofs_in_region()

  - Term, Terms:

    - update ConnInfo for no element groups - remove ConnInfo.iter_igs()
    - remove vector_chunk_generator(), CharacteristicFunction()
    - remove unused Term.__call__(), ._call()
    - remove Terms.set_current_group()
    - update Term for no element groups

      - remove .char_fun attribute
      - remove .get_current_group(), .set_current_group(), .igs(),
        .iter_groups()
      - update .setup(), .get_conn_info(), .get_args(), .get_approximation(),
        .get_physical_qps(), .get_mapping(), .get_data_shape(), .get(),
        .evaluate(), .assemble_to()
      - update docstring of .assign_args(), .check_args()

    - fix Term.eval_complex() - force complex output also for real arguments

  - remove sfepy/terms/terms_new.py - NewTerm and subclasses

    - update FieldVariable for no NewTerm

      - remove .clear_bases(), .setup_bases(), .clear_current_group(),
        .set_current_group(), .val(), .val_qp(), .grad(), .grad_qp(),
        .iter_dofs(), .get_element_zeros(), .get_component_indices()
      - update .__init__()

    - remove compare_scalar_terms.py, compare_vector_terms.py examples

  - update PhysicalQPs for no element groups

    - update .__init__(), .get_shape()
    - remove .get_merged_values()

  - update Material for no element groups

    - update .set_data(), .set_data_from_variable(), .update_data(),
      .get_data(), ._get_data(), .get_constant_data(), .reduce_on_datas()

  - update eval_nodal_coors() for no element groups
  - update Approximation and subclasses for no element groups

    - update Approximation.__init__(), .eval_extra_coor(),
      .setup_surface_data(), .setup_point_data(), .get_connectivity(),
      .describe_geometry()
    - remove Approximation.ig attribute
    - update DiscontinuousApproximation.eval_extra_coor()
    - update SurfaceApproximation.__init__()
    - update FESurface for no element groups - update .__init__(),
      .setup_mirror_connectivity()

  - update Equations.get_graph_conns() for no element groups
  - update EquationMap.map_equations() for no element groups
  - update create_adof_conns() for no element groups
  - update Mapping.from_args() for no element groups - update
    get_physical_qps(), get_mapping_data()

  - update FieldVariable for no element groups

    - update .get_mapping(), .get_dof_conn(), .time_update(),
      .setup_initial_conditions(), .get_approximation(), .get_data_shape(),
      .evaluate(), .get_state_in_region(), .get_element_diameters()

  - miscellaneous updates:

    - clean up meshutils.[ch]
    - clean up sfepy/discrete/fem/fea.py
    - update set_mesh_coors() for current Mesh
    - update refine_reference() for no element groups
    - update compute_nodal_normals(), extend_cell_data() for no element groups
    - update StabilizationFunction.__call__() for no element groups
    - update ConstantFunctionByRegion.__init__() for no element groups

  - scripts:

    - script/extract_surface.py: update for no element groups

      - use Mesh.from_region()
      - allow surface mesh saving also for 2D meshes

    - script/save_basis.py: update save_basis_on_mesh() for no element groups

  - update homogenization recovery functions for no element groups

    - update get_output_suffix(), .recover_bones(), .recover_paraflow(),
      .recover_micro_hook()

  - terms:

    - update membrane describe_geometry() for no element groups
    - update TLMembraneTerm for no element groups - update .__init__(),
      .get_fargs()
    - update hyperelastic terms for no element groups

      - update HyperElasticBase.__init__(), .get_fargs()
      - update BulkPressureTLTerm.get_fargs()
      - update BulkPressureULTerm.get_fargs()

    - update contact surface terms for no element groups

      - update ContactPlaneTerm.get_fargs()
      - update ContactSphereTerm.get_fargs()

    - remove commented-out code in hyperelastic terms

  - examples and tests:

    - fix mesh_hook() in examples/diffusion/laplace_1d.py
    - update examples/linear_elasticity/its2D_3.py for no element groups
    - update examples/homogenization/rs_correctors.py for no element groups, use
      cells of channels regions directly
    - update tests for no element groups and current Mesh

  - add docstring to FieldVariable.get_state_in_region()
  - docs: update for current sources

- improve finding reference element coordinates of physical points

  - merge branch improve-find-ref-coors - partially fixes #285
  - update gtr_normalize_v3() for 2D - add dim argument
  - new mesh_get_facet_normals() C function
  - new CMesh.get_facet_normals()
  - update gtr_dot_v3() for 2D - add dim argument
  - update orient_elements() for current gtr_dot_v3()
  - fix tensor-product mode in inverse_element_mapping()
  - split cmesh.pyx -> cmesh.pyx + cmesh.pxd
  - add verbose argument to gtr_normalize_v3()
  - new refc_find_ref_coors() C function in new refcoors.c

    - new sfepy/discrete/fem/extmods/refcoors.[ch]
    - new _mul_c_add_v3(), _intersect_line_plane(), _intersect_line_triangle(),
      _get_cell_coors(), _get_tri_coors(), _get_xi_dist()

  - new find_ref_coors() in new crefcoors.pyx extension module

    - new sfepy/discrete/fem/extmods/crefcoors.pyx
    - update sfepy/discrete/fem/extmods/setup.py

  - update get_ref_coors() for new find_ref_coors() C function
  - remove old find_ref_coors() Cython function
  - move FEDomain.get_evaluate_cache() -> FEField.get_evaluate_cache() - update
    for new find_ref_coors()
  - set failed values to nan in H1NodalMixin.evaluate_at() if not returning
    status, fix and update docstring
  - update Probe.probe() for current code

- reorganize main scripts

  - merge branch reorganize-scripts
  - add main sfepy-run wrapper and scripts-common dir
  - add sfepy.in_source_tree dependency
  - update doc build scripts for man builder, add sfepy manpage
  - add scripts-common links and modify setup call parameters
  - update user's guide for sfepy-run command

- scripts:

  - postproc.py: new --parallel-projection option - update Viewer.call_mlab(),
    .render_scene()

- miscellaneous updates:

  - change default verbose value to False in FieldVariable.evaluate_at() +
    related
  - update FEDomain.get_evaluate_cache() for new verbose argument
  - force solution in make_l2_projection_data() by decreasing absolute
    tolerance
  - add boilerplate in eval_ns_forms.py
  - fix probe constructors to obey share_geometry - fix PointsProbe, LineProbe,
    RayProbe, CircleProbe
  - new Probe.get_evaluate_cache(), update .probe()
  - add nonlinear solver options argument to projection functions -
    project_by_component(), update make_l2_projection(),
    make_l2_projection_data()
  - fix memory leak - new CMesh.__dealloc__(), new mesh_free() C function
  - fix graph_components() (according to cs_graph_components() in scipy.sparse)
  - add igakit as optional dependency in setup.py, update sfepy/version.py
  - new mesh function: expand2d - converts tri/quad mesh to tetra/hexa
  - allow matrix_hook to handle residual vector
  - update splinebox.py: number of segments as an input argument to __init__()
  - remove unused variable in evaluate_in_rc()
  - update plot_wireframe() for non-string color specifications
  - fix DOWNLOAD_URL
  - fix has_several_times, has_several_steps in Viewer.call_mlab()
  - speed up GenericFileSource.create_dataset()
  - check number of dimensions of material array in Material.set_data()

- examples and tests:

  - clean up tests/test_linear_solvers.py
  - new example: Darcy flow in multiple compartments
  - new tests/test_input_darcy_flow_multicomp.py
  - update test_install.py to test postproc.py with HDF5 output files

- terms:

  - update LinearTractionTerm docstring - explain load shapes
  - new term: dw_vm_dot_s

    - replace mulATB_integrate() by mulAB_integrate() with multiple modes; new:
      actBfT()

- docs:

  - sync module index of developer guide with current sources
  - update release tasks for sfepy-run

.. _2014.4-2015.1:

from 2014.4 to 2015.1
=====================

- support multiple fields in IGA

  - merge branch iga-multifield
  - move NurbsPatch-related code in IGDomain.__init__() to
    NurbsPatch.__init__()
  - update IGMapping for nurbs patch other than domain.nurbs

    - update IGMapping.__init__() - new nurbs argument
    - update docstring

  - allow 'iga*' approximation order in parse_approx_order()
  - new NurbsPatch.elevate()
  - update IGField.__init__() for degree elevation of NURBS basis

    - new approx_order argument in IGField.__init__()
    - new parse_approx_order()

  - use self.nurbs instead of self.domain.nurbs in IGField

    - update .get_econn(), .set_dofs(), .create_mapping(), .create_output()
    - new ._get_facets()
    - remove IGDomain.facets

  - new meshes/iga/block2d.iga
  - new examples/navier_stokes/navier_stokes2d_iga.py + test

- IGA:

  - new eval_bspline_basis_tp() - use if all weights are one - new is_nurbs(),
    update igac.pyx
  - pre-compute Bernstein basis in eval_mapping_data_in_qp(),
    eval_variable_in_qp() - update eval_bspline_basis_tp(),
    eval_nurbs_basis_tp()
  - fix get_facet_axes() in 2D
  - fix get_surface_degrees()
  - fix connectivity types in create_connectivity_1d()
  - update get_bezier_element_entities() to return also vertices
  - new create_from_igakit(), update gen_patch_block_domain()
  - new NurbsPatch._to_igakit(), ._from_igakit(), update .elevate()
  - update IGField.__init__() for str or None approx_order values

- redesign handling of solver parameters

  - merge branch solver-parameters
  - new SolverMeta

    - new format_next(), typeset_to_indent(), make_option_docstring()
    - new par_template

  - update Solver to use SolverMeta, new generic Solver.process_conf()
  - update LinearSolver.get_tolerance() for no tolerances in conf, clean up
  - update basic linear solvers for SolverMeta, new ._parameters class
    attribute

    - update ScipyDirect, ScipyIterative, PyAMGSolver, PETScKrylovSolver,
      PETScParallelKrylovSolver - remove .process_conf()

  - update nonlinear solvers for SolverMeta, new ._parameters class attribute

    - update Newton, ScipyBroyden, Oseen - remove .process_conf()

  - support extra parameters in Solver.process_conf() - marked by '*' name in
    _parameters solver attribute

  - update optimization solvers for SolverMeta, new ._parameters class
    attribute

    - update FMinSteepestDescent, ScipyFMinSolver - remove .process_conf()

  - update remaining linear solvers for SolverMeta, new ._parameters class
    attribute

    - update SchurGeneralized, SchurComplement, MultiProblem - remove
      .process_conf()

  - remove 'needs_problem_instance' solver parameter
  - remove 'needs_problem_instance' from examples
  - rename Newton parameter 'problem' -> 'is_linear'

    - update Newton.parameters, .__call__()
    - update Problem.is_linear(), .set_linear()

  - docs: update for 'is_linear'
  - update examples for 'is_linear' (remove default setting, reformat if
    needed)
  - update SemismoothNewton for SolverMeta, new ._parameters class attribute -
    remove .process_conf()
  - update time-stepping solvers for SolverMeta, new ._parameters class
    attribute

    - update StationarySolver, EquationSequenceSolver,
      SimpleTimeSteppingSolver, ExplicitTimeSteppingSolver,
      AdaptiveTimeSteppingSolver - remove .process_conf()

  - update eigenvalue solvers for SolverMeta, new ._parameters class attribute

    - update ScipyEigenvalueSolver, ScipySGEigenvalueSolver,
      LOBPCGEigenvalueSolver, PysparseEigenvalueSolver - remove .process_conf()

  - use SchurComplement solver in linear_elastic_up.py example

- solvers:

  - remove ls.umfpack solver (class Umfpack), replace by ls.scipy_direct
  - clean up sfepy/solvers/ls.py
  - remove unneeded code in make_implicit_step()
  - fix default i_max of LOBPCGEigenvalueSolver, use verbose option
  - obey verbose option and eigenvectors argument in PysparseEigenvalueSolver

    - fix default i_max, eps_a
    - move imports to .__init__() for early failure

  - new Solver.build_solver_kwargs()
  - refactor ScipyEigenvalueSolver, ScipySGEigenvalueSolver

    - ScipyEigenvalueSolver uses dense or sparse scipy linalg functions, new
      method option, support extra options
    - ScipySGEigenvalueSolver uses dense lapack functions, remove force_n_eigs
      option
    - imports moved to .__init__() for early failure

  - use Solver.build_solver_kwargs() in ScipyFMinSolver
  - clean up sfepy/solvers/optimize.py
  - clean up sfepy/solvers/oseen.py
  - clean up sfepy/solvers/semismooth_newton.py

- scripts:

  - script/gen_iga_patch.py: new --cp-mode option, update
    gen_patch_block_domain()
  - script/gen_gallery.py: update views for centered scene

- miscellaneous updates:

  - new Region.from_cells()
  - clean up tests/sympy_operators.py
  - fix get_mem_usage() (workaround of pyparsing object with iteritems of str
    type)
  - update units_of_quantities (add density), prefixes (add 'p', 'T')
  - fix GenericFileSource.get_bounding_box()
  - fix mappings, QP and base function key collisions

    - use integral order instead of name in mapping keys, in QP and base
      function keys
    - make region name the first item in mapping keys
    - update Approximation, SurfaceApproximation, Field, FieldVariable, Term,
      NewTerm

  - update QuadraturePoints docstring
  - support all QuadraturePoints.__init__() arguments in Integral.__init__() -
    update Integral.get_qp()
  - update actor positions in Viewer.build_mlab_pipeline() to center scene
  - fix sympy.zeros() calls for sympy 0.7.6
  - new bspline functions (curve and surface)
  - update SplineBox: use bspline.py functions, new test
  - update BSpline.basis_function_dg() for degree 0
  - add linear solver argument to make_l2_projection(),
    make_l2_projection_data()
  - fix SchroedingerApp.solve_eigen_problem() to request eigenvectors - fix for
    updated PysparseEigenvalueSolver

- tests and examples:

  - new examples/standalone/interactive/modal_analysis.py
  - update test_install.py to test modal_analysis.py example
  - update examples/navier_stokes/navier_stokes2d.py - generated domain, new BC
  - remove unused rectangle_fine_quad.mesh, rectangle_fine_tri.mesh
  - regenerate IGA meshes with greville control points mode
  - new tests/test_input_linear_elastic_iga.py, tests/test_input_poisson_iga.py
  - new tests/test_eigenvalue_solvers.py

- docs:

  - sync module index of developer guide with current sources
  - update citing section
  - update list of applications
  - update list of IGA examples
  - update description of IGA field definition
  - update installation notes:

    - add Anaconda distribution info
    - update tested versions
    - bux-fixes
    - minor reorganization and clean up

.. _2014.3-2014.4:

from 2014.3 to 2014.4
=====================

- preliminary support for 1D problems

  - merge branch 1d
  - update CMesh for 1D (preliminary support)
  - update Region, Domain and FEDomain for 1D, update orient_elements() -
    update geometry_data
  - update diffusion and Navier-Stokes terms for 1D
  - update Term.geometries
  - update tests/test_term_call_modes.py for 1D
  - update C implementation of hyperelasticity for 1D
  - update VTKMeshIO.write() for 1D
  - add 1D reference elements meshes (1_2_1.mesh, 1_2_2.mesh)
  - new examples/diffusion/laplace_1d.py + test
  - update SDLinearElasticTerm for 1D
  - support custom material parameter values in tests/test_term_call_modes.py

    - update make_term_args()
    - add custom value for 'dw_s_dot_grad_i_s' in Test._test_single_term()

  - add custom view for laplace_1d.py example to script/gen_gallery.py
  - update GenericFileSource.get_bounding_box(), .create_dataset() for 1D

- merge pull request #276 (VTK probes)

  - simplify probes, use pyVTK library
  - Primer - use VTK probes
  - new postrocessing function based on VTK, example
  - update User's guide - postrocessing filters

- interactive example for Primer

  - merge branch its-interactive-example
  - update nodal_stress() to accept user-provided integrals
  - new examples/linear_elasticity/its2D_interactive.py
  - docs: add interactive example section to primer

- scripts:

  - script/gen_iga_patch.py: report number of DOFs per axis
  - fix path in script/gen_lobatto1d_c.py - bug introduced in
    ec186b8c69b0719f483e11a864092a21a608a203
  - clean up script/edit_identifiers.py
  - fix sys.path in auxiliary scripts and standalone examples - for cases where
    the current directory (i.e. top level sfepy directory) is not added
    automatically
  - script/gen_gallery.py:

    - omit examples that use VTK probes - bad interaction, probably with
      Mayavi, leads to segfault
    - simplify section titles, generate better labels

- miscellaneous updates:

  - fix number of returned values for no refinement in refine_reference()
  - fix make_l2_projection_data() for parameter variables
  - IGA: update plot_iso_lines() for curves in 2D

- tests and examples:

  - update test_install.py for current linear_elastic_probes.py
  - prevent modifying linear_elastic.py namespace in derived examples
  - new tests/test_input_linear_elastic_probes.py

- docs:

  - fix testing instructions for installed build
  - update user's guide for interactive probing example in primer
  - fix IGA section for current sources
  - clean up script/gen_term_table.py, use common help and docstring
  - update script/gen_term_table.py: simpler table, correct links to term
    classes
  - update doc/index.rst for 1D and new gallery labels
  - sync doc/introduction.rst with doc/index.rst

.. _2014.2-2014.3:

from 2014.2 to 2014.3
=====================

- speed-up IGA by C implementation of NURBS basis evaluation

  - merge branches iga-c and iga-c-2
  - new eval_bernstein_basis() in C in sfepy/discrete/iga/extmods/

    - new sfepy/discrete/iga/extmods/nurbs.[ch]
    - new sfepy/discrete/iga/extmods/igac.pyx
    - new setup.py files

  - support DEBUG_FMF define in sfepy/discrete/iga/extmods/setup.py
  - update IGDomain.__init__() to have nurbs.cs in FMField compatible shape
  - new ravel_multi_index(), unravel_index() in C
  - new eval_nurbs_basis_tp() in C
  - new array2puint2(), array2puint1() helper functions
  - change return value of compute_bezier_extraction_1d() to 3D array
  - change nurbs.degrees dtype to int32 in IGDomain.__init__()
  - new eval_mapping_data_in_qp(), eval_variable_in_qp() in C - use in
    sfepy/discrete/iga/mappings.py
  - new eval_in_tp_coors() in C
  - update NurbsPatch.__call__() to use eval_in_tp_coors()

- generalize linear combination boundary conditions to work between different
  fields/variables and to support non-homogeneous periodic conditions (non-zero
  right-hand side or shift)

  - merge branch 'lcbc' - closes #179, #267
  - update transform_lcbcs() for new LCBC syntax

    - the syntax: region[, times], dofs, dof_map_fun, kind[, other arguments]

  - update Conditions.from_conf(), LinearCombinationBC.__init__() for new
    syntax

    - new LinearCombinationBC.canonize_dof_names()

  - new are_disjoint()
  - new LinearCombinationBC.get_var_names()
  - check that regions are disjoint in Variables.setup_lcbc_operators()
  - update Variables.setup_lcbc_operators() for new operators (to be done)
  - remove FieldVariable.create_lcbc_operators()
  - new kind attribute of LCBC operators
  - update LCBCOperator for two variables in two regions

    - rename & update LCBCOperator -> MRLCBCOperator (model reduction)
    - new general LCBCOperator for tying two variables in two regions

  - reimplement LCBCOperators, use common LCBCOperator constructor signature

    - update .__init__(), .add_from_bc(), .append(), .finalize()
    - move make_global_lcbc_operator() -> new
      LCBCOperators.make_global_operator()

  - new ShiftedPeriodicOperator
  - provide some defaults for long syntax in transform_lcbcs()
  - store field in MRLCBCOperator
  - update LCBC operators new base classes - update IntegralMeanValueOperator,
    RigidOperator, NoPenetrationOperator, NormalDirectionOperator
  - fix MRLCBCOperator active equation map for non-scalar variables
  - update MRLCBCOperator.treat_pbcs()
  - update Variables.setup_lcbc_operators(), .__init__(),
    FieldVariable.__init__()

    - new Variables.has_lcbc_rhs
    - rename Variables.op_lcbc, .rhs_lcbc -> .mtx_lcbc, .vec_lcbc
    - remove FieldVariable.has_lcbc

  - set vec_lc to None in LCBCOperators.make_global_operator() when no rhs
  - update LCBCEvaluator, eval_equations() for renamed operators
  - add comment to State.get_reduced()
  - update meshes/2d/square_quad.mesh to have two cell groups
  - simplify transform_conditions() - LCBCs treated in transform_lcbcs()
  - docs: update essential boundary conditions description, add EPBCs, LCBCs
  - tests and examples:

    - update examples and tests for new LCBC syntax
    - new examples/diffusion/laplace_coupling_lcbcs.py + test
    - new examples/linear_elasticity/elastic_shifted_periodic.py + test
    - new examples/standalone/interactive/laplace_shifted_periodic.py
    - new tests/test_lcbcs.py, test_laplace_shifted_periodic()

  - add custom views for laplace_coupling_lcbcs.py, elastic_shifted_periodic.py
    to script/gen_gallery.py

- support non-constant essential boundary conditions given by a function in IGA

  - merge branch 'iga-ebc'
  - update IGField.get_econn() for surface connectivities
  - simplify IGField.get_dofs_in_region_group() by using .get_econn()
  - new get_facet_axes()
  - new solve() convenience wrapper
  - new create_boundary_qp()
  - new get_surface_degrees()
  - update IGField.set_dofs() for DOFs given by function, update docstring -
    initial implementation
  - examples:

    - update examples/linear_elasticity/linear_elastic_iga.py for EBCs by
      function

  - update custom view for linear_elastic_iga.py example in
    script/gen_gallery.py

- IGA:

  - fix IGDomain.__init__() for standalone use - reset regions
  - fix IGField.create_output() for no key argument
  - new sfepy/discrete/iga/domain_generators.py - new gen_patch_block_domain()
  - fix plot_nurbs_basis_1d to use ax argument

- regions:

  - new Region.kind_tdim attribute, new .set_kind_tdim(), update .set_kind()
  - make ig argument optional in Region.get_entities(), signature change
  - update Region.contains() to use other region kind entities
  - make sure entities are initialized in Region.get_entities()
  - update Region.setup_from_highest() - new allow_lower argument

    - works also when the highest dimension of entities is lower than requested
    - fixes selectors like 'vertices of surface *f vertices in (x < 0.0)'

  - new Region.finalize(), update Domain.create_region()
  - fix region definitions in tests and examples as reported by
    Region.finalize()

- postprocessing and visualization:

  - remove import lines from sfepy/postprocess/__init__.py, update affected
    imports
  - update Viewer.call_mlab(): do not show GUI if show option is False

- applications:

  - update SchroedingerApp.save_results() for custom mesh in out['__mesh__']
  - relax too strict assertion in DensityVolumeInfo.__call__()
  - new: micro-recovery hook allows to return some values to the "master"
    macro-problem

- scripts:

  - script/gen_iga_patch.py: use gen_patch_block_domain()
  - script/gen_iga_patch.py: plot 1D basis along patch mid-line

- miscellaneous updates:

  - remove prints and config.py file generation from setup.py files
  - support user-specified ebcs, epbcs in Problem.save_ebc(), update docstring
  - fix Newton.__init__(), Oseen.__init__() to obey log plot option
  - normalize paths in top_dir, data_dir and base_dir
  - update import_file() to use full package paths for sfepy modules - prevents
    name clashes for different fields.py files
  - add can_reload argument to import_file()
  - do not force module reload in load_classes()
  - fix issue #263, ConstantFunctionByRegion
  - clean up sfepy/discrete/conditions.py
  - use os.sep instead of '/' for the pathname separator in import_file() - fix
    portability issue
  - fix extend_cell_data() for several groups, new average_surface argument
  - change keys in transform_epbcs()
  - report unknown LCBC kind in LCBCOperators.add_from_bc()
  - clean up run_tests.py, sfepy/base/testing.py, sfepy/discrete/materials.py,
    tests/tests_basic.py, sfepy/mechanics/matcoefs.py
  - fix set_mesh_coors() - update also cmesh coordinates
  - update get_debug() to allow frame specification in debug()
  - update MED file support, fix this issue:
    https://groups.google.com/forum/#!msg/sfepy-devel/z1e83Xgl1_U/884jDfKzKPsJ
  - fix setup.py - add missing plot_logs.py script
  - remove unused sfepy/physics/energy.py

- examples:

  - do not force examples to be in a package - update import_file(),
    ProblemConf.from_file()
  - remove examples/diffusion/octahedron.py + test - no added value
    w.r.t. other diffusion examples
  - update examples/linear_elasticity/linear_elastic.py: comments -> docstring

- docs:

  - module index of developer guide with current sources
  - fix docs generation without igakit
  - fix typeset_term_table() for import_file() using full package paths
  - add missing figure
  - replace alphas8bit with alphas - word_free with alphas8bit as default
    argument value in get_standard_type_defs() and list_dict() breaks sphinx
    docs generation on platforms with non utf-8 locale
  - new section on exploring the code in developer guide
  - fix download page
  - migrate tutorial "Using Salome with SfePy" from google site, see #171
  - merge branch 'theory-docs'

    - new doc/ebcs_implementation.rst - explanation of EBCs implementation
    - move IGA section from tutorial to user's guide
    - new doc/theory.rst

      - move notes on solving PDEs from doc/tutorial.rst to new
        doc/solving_pdes_by_fem.rst
      - add doc/ebcs_implementation.rst to doc/theory.rst
      - update doc/index.rst

    - remove doc/notes.rst

  - simplify paragraph in doc/index.rst
  - simplify doc/introduction.rst
  - put Primer under Examples
  - merge pull request #275 (short syntax)

    - remove long syntax from the user's guide
    - new doc page with long syntax
    - tutorial with short syntax of keywords, closes #272

  - add introductory paragraph to examples

.. _2014.1-2014.2:

from 2014.1 to 2014.2
=====================

- preliminary support for isogeometric analysis

  - merge branch 'iga'
  - new iga.py: compute_bezier_extraction_1d(), _get_knots_tuple(),
    eval_bernstein_basis(), compute_bezier_extraction(), get_raveled_index(),
    tensor_product(), combine_bezier_extraction(), create_connectivity_1d(),
    create_connectivity(), compute_bezier_control(), get_unraveled_indices(),
    eval_nurbs_basis_tp(), eval_variable_in_qp(), get_bezier_topology(),
    get_patch_box_regions(), get_bezier_element_entities(),
    eval_mapping_data_in_qp()
  - new plot_nurbs.py: plot_parametric_mesh(), plot_control_mesh(),
    plot_iso_lines(), plot_nurbs_basis_1d(), plot_bezier_mesh(), _get_edges(),
    plot_bezier_nurbs_basis_1d()
  - new utils.py: create_linear_fe_mesh(), create_mesh_and_output(),
    save_basis()
  - new io.py: write_iga_data(), read_iga_data()
  - add filename_domain as alternative to required keywords
  - update PDESolverApp.setup_output_info() for filename_domain keyword
  - simplify Domain.has_faces()
  - remove unused Domain.get_conns()
  - initialize shape in Domain.__init__()
  - move common Domain code into new sfepy/discrete/common/domain.py

    - move & update region_leaf()
    - move region_op(), Domain.get_centroids(), .has_faces(), .reset_regions(),
      .create_region(), .create_regions(), .save_regions(),
      .save_regions_as_groups()
    - rename original Domain -> FEDomain, update FEDomain.__init__()
    - new Domain.__init__()

  - new sfepy/discrete/iga/domain.py - new IGDomain

    - new .from_file(), .__init__()
    - update Problem.from_conf() for FEDomain, IGDomain
    - new NurbsPatch - ._get_ref_coors_1d(), .__call__(), .evaluate()

  - move common Field code into new sfepy/discrete/common/fields.py

    - move Field.from_args(), .from_conf(), ._setup_kind()
    - move parse_approx_order(), setup_extra_data(), fields_from_conf()
    - rename original Field -> FEField, update FEField.__init__()
    - update VolumeField, SurfaceField parent class
    - new parse_shape(), update FEField.__init__()

  - new sfepy/discrete/iga/fields.py - new IGField

    - new .__init__(), .get_true_order(), .is_higher_order()
    - new .get_econn(), .setup_extra_data() - volume only
    - new IGField.get_dofs_in_region_group()
    - new IGField.set_dofs() - only for a constant value over an entire patch
      side
    - new IGField.create_mapping()
    - new IGField.create_output() - initial version with hard-coded values

      - output mesh corresponding to uniform parametric mesh is returned as
        out['__mesh__']

    - update Problem.save_state() for custom mesh in out['__mesh__']
    - update parse_approx_order(), Field.from_conf() for IGField

  - get dimension from domain in FieldVariable._set_field()
  - move FEField.get_dofs_in_region() into Field
  - move common Mappings code into new sfepy/discrete/common/mappings.py

    - move PhysicalQPs, get_mapping_data(), get_jacobian(), get_normals()
    - replace create_mapping() by new Mapping.from_args()
    - rename original Mapping -> FEMapping - move .get_physical_qps()
    - update VolumeMapping, SurfaceMapping parent class

  - move mapping handling code from FEField to Field

    - move .clear_mappings(), .save_mappings(), .get_mapping()
    - update Field.get_mapping() to use .create_mapping(), fix docstring
    - update FEField.create_mapping() to return both CMapping and Mapping

  - do not use domain shape in Field.get_mapping() - remove sharing of full
    group mappings
  - do not use mesh coordinates in Material.update_special_data()
  - new sfepy/discrete/iga/mappings.py - new IGMapping

    - new .__init__(), .get_geometry(), .get_physical_qps(), .get_mapping()
    - update Mapping.from_args() for IGMapping

  - use Field.get_econn() in FieldVariable.evaluate()

    - do not use Approximation.get_connectivity()
    - update VolumeField.get_econn() - add integration argument
    - update signature of SurfaceField.get_econn(), IGField.get_econn()

  - move FieldVariable.get_data_shape() functionality to Field subclasses

    - new FEField.get_data_shape(), IGField.get_data_shape()

  - remove Approximation.get_v_data_shape(), .get_s_data_shape()
  - update H1DiscontinuousField.extend_dofs(), .remove_extra_dofs() for FEField
  - new is_surface attribute of fields

    - update VolumeField._setup_geometry(), SurfaceField._setup_geometry()
    - update FieldVariable._set_field() - does not import SurfaceField from
      sfepy.discrete

  - update sfepy/discrete/__init__.py for common Domain, Field
  - update sfepy/discrete/fem/__init__.py for FEDomain
  - update domain creating code for FEDomain
  - update 'surface_extra' assembling cells in Term.get_assembling_cells()

    - for Field.get_econn() in FieldVariable.evaluate() - surface regions have
      no more a complete group DOF connectivity

  - update 'volume' Term.get_assembling_cells() for subdomains

    - for Field.get_econn() in FieldVariable.evaluate() - subdomain volume
      regions have no more a complete group DOF connectivity
    - trivial now, could be removed

  - update FieldVariable.setup_initial_conditions() to use Field.set_dofs()

  - tests and examples:

    - fix get_pars() in tests/test_functions.py for None in coors argument
    - update import in _integrate() in tests/test_term_consistency.py
    - new meshes/iga/patch2d.iga
    - new examples/diffusion/poisson_iga.py
    - new meshes/iga/block3d.iga
    - new examples/linear_elasticity/linear_elastic_iga.py
    - update test_install.py to test poisson_iga.py example

  - scripts:

    - new script/gen_iga_patch.py
    - add custom view for poisson_iga.py example to script/gen_gallery.py,
      sort custom dict by keys
    - add custom view for linear_elastic_iga.py example to script/gen_gallery.py

  - docs:

    - update developer guide for IGA, split sfepy.discrete contents list
    - new IGA section in tutorial, update links and index
    - update optional installation dependencies

- postprocessing and visualization:

  - postproc.py: update --layout option for row#n, col#n layouts - update
    get_position_counts(), Viewer.get_size_hint()
  - always use cell-to-point filter for domain-specific plots

    - point scalars work as color_name
    - update Viewer.build_mlab_pipeline()

  - update file sources for providing both steps and times

    - update FileSource.reset(), .get_step_range()
    - remove FileSource.set_step(), VTKFileSource.get_step_range(),
      VTKSequenceFileSource.get_step_range(),
      GenericFileSource.get_step_range()
      GenericSequenceFileSource.get_step_range()
    - update GenericFileSource.__init__(), .read_common(), .file_changed()
    - new FileSource.get_step_time(), .get_ts_info()
    - new VTKSequenceFileSource.__init__()
    - new GenericSequenceFileSource.read_common()

  - update postproc.py, Viewer for updated file sources, add slider for times

    - update Viewer.call_mlab()
    - update SetStep - new ._time_editor, .time, .time_low, .time_high,
      .is_adjust ._time_changed()
    - update .__source_changed(), ._step_changed(), ._file_changed_changed()

  - update Viewer.save_animation() for selected steps or times

    - update Viewer.get_animation_info(), .encode_animation()
    - remove FileSource.get_step_range()

  - remove animate, anim_format arguments from Viewer.__init__()

    - update view_file(), postproc.py
    - update Viewer.call_mlab(), make_animation()

  - update Viewer, ViewerGUI for saving snapshots in given ranges of steps/times

    - update Viewer.call_mlab(), make_animation()
    - update SetStep:

      - new .seq_* traits
      - new .init_seq_selection(), ._seq_n_step_changed(), ._seq_dt_changed()

    - new ClosingHandler
    - update ViewerGUI:

      - custom quit button
      - remove button_make_animation, button_make_snapshots buttons
      - new button_make_animation_steps, button_make_snapshots_steps
        button_make_animation_times, button_make_snapshots_times buttons + their
        handlers

  - postproc.py: new --time option

    - change default of --step to None, --step has precedence over --time
    - update view_file(), Viewer.call_mlab()

  - postproc.py: improve --list-ranges

    - for HDF5 result files, output union of ranges for all steps if --step or
      --time are not given
    - output in form suitable for --ranges

  - fix suffix returned by Viewer.get_animation_info()
  - remove auto_screenshot argument of Viewer.__init__() - update .call_mlab()
  - postproc.py: update --output option behaviour, update docstring

- input-output:

  - update defaults for missing output data attributes in HDF5MeshIO.write() -
    only mode and data are required
  - provide better default values in MeshIO.read_times()

- scripts:

  - new script/plot_logs.py - plot text file logs of variables saved by Log
  - update schroedinger.py: remove mesh generation options, simplify

    - remove --create-mesh, --2d, --mesh, --mesh-dir options
    - new --n-eigs, --tau options for updated examples/quantum/
    - update docstring

  - script/extract_surface.py: fix default value of --mesh option
  - update script/gen_gallery.py for current Viewer - save image explicitly
  - fix docstring formatting of script/convert_mesh.py

- terms:

  - make material argument of dw_tl_surface_traction optional - update
    SurfaceTractionTLTerm, dw_tl_surface_traction()
  - new VolumeSurfaceTLTerm (d_tl_volume_surface) - new d_tl_volume_surface()
  - move SurfaceTractionTLTerm.compute_family_data() to new common class

    - new HyperElasticSurfaceTLBase
    - update SurfaceTractionTLTerm, VolumeSurfaceTLTerm

  - new SurfaceFluxTLTerm (d_tl_surface_flux) - new d_tl_surface_flux()
  - fix setting cells of normals in dw_surface_v_dot_n_s(),
    dw_surface_s_v_dot_n()
  - fix arg_shapes of LaplaceTerm
  - new SurfaceFluxOperatorTerm (dw_surface_flux) - new dw_surface_flux()
  - fix mat_le_stress() for heterogeneous parameters in 2D

- solvers:

  - update cm_pb solver: add "master problem" to kwargs
  - improve linear solution precision warning message in Newton.__call__()
  - new eps_mode parameter of Newton solver, update conv_test()
  - clean up sfepy/solvers/nls.py, add module and conv_test() docstrings
  - update PETScParallelKrylovSolver - new log_dir option - closes #250
  - remove obsolete is_plot option from Newton and Oseen solvers

- miscellaneous updates:

  - fix create_adof_conns(): problem with trace
  - fix docstring of Region.get_facet_indices()
  - new get_simplex_volumes()
  - update read_log() for general (float) x values
  - fix read_log() for vlines before actual data
  - fix solve_pde() to preserve conf.options class
  - use tight layout in LogPlotter.process_command() after plotting labels
  - fix surface mode in extend_cell_data()
  - fix Problem.create_evaluable() to use copies of problem variables
  - update test_install.py for meshes/quantum/ and schroedinger.py changes
  - use tempfile in gen_misc_mesh()
  - remove mesh coordinate transformation code in Problem.from_conf()
  - remove tmp/.ignore
  - update make_l2_projection_data() for array data, new order argument, fix
    overwriting of target._variables
  - new get_condition_value(), update EquationMap.map_equations()
  - update EquationMap.map_equations() to use get_condition_value() for EPBCs

- tests and examples:

  - increase required precision in examples/acoustics/acoustics3d.py
  - allow trying mesh refinement in tests/test_linear_solvers.py
  - clean up tests/test_volume.py
  - update tests/test_volume.py to test TL volume terms, new test_volume_tl()
  - update examples/quantum/ - use meshes/quantum/square.mesh

    - update quantum_common.py to use square.mesh, add docstring, increase
      approximation order to 2, update common() arguments
    - new arguments to define() in examples: n_eigs, tau
    - remove unused meshes in meshes/quantum/
    - update meshes/quantum/square.geo
    - new meshes/quantum/square.mesh (corresponds to meshes/quantum/square.geo)

- docs:

  - update developer guide for sfepy/discrete/
  - update support section
  - update years in LICENSE
  - replace link by contents of LICENSE, new doc/license.rst
  - update year in doc/conf.py
  - remove obsolete Intel Mac installation instructions
  - update platform-specific notes
  - simplify INSTALL, update installation and links
  - new citing section in index

.. _2013.4-2014.1:

from 2013.4 to 2014.1
=====================

- new handling of field and variable shapes:

  - use region space dimension instead of field shape in vector LCBC operators
  - use DOF per node count instead of field shape in IntegralMeanValueOperator
  - update Field.from_args(), __init__() docstrings (shape parameter)
  - update Field: new n_components, val_shape attributes, new H1Mixin
  - inherit H1NodalMixin from H1Mixin, update .set_dofs()
  - inherit H1HierarchicVolumeField from H1Mixin, update .set_dofs()
  - remove n_components argument in Variable, FieldVariable constructors

    - update Variable.__init__(), FieldVariable.__init__()
    - update Variable._setup_dofs() - new n_components, val_shape arguments
    - update FieldVariable._set_field() - set n_components, val_shape from field

  - update sfepy.fem modules for no n_components argument
  - update scripts, examples and tests for no n_components argument
  - update docs for no n_components argument
  - remove Field.setup_dof_conns() (duplicated in VolumeField.setup_dof_conns())
  - remove dpn argument of {VolumeField, SurfaceField}.setup_dof_conns()

- create active DOF connectivities directly from field connectivities:

  - fields have no dof_conns dict
  - remove create_dof_conn(), setup_dof_conns(), Field.clear_dof_conns()
    VolumeField.setup_dof_conns(), SurfaceField.setup_dof_conns()
  - new VolumeField.get_econn(), SurfaceField.get_econn()
  - update Field.__init__()
  - remove Variables.setup_adof_conns(), FieldVariable.setup_adof_conns()
  - new create_adof_conns(), Variables.set_adof_conns()
  - update create_adof_conn(), Variable.__init__()
  - update sfepy.fem for create_adof_conns()

    - remove Equations.setup()
    - remove setup, make_virtual, verbose arguments of Equations.__init__()
    - remove setup, make_virtual arguments of Equations.from_conf()
    - update Equations.time_update(), create_evaluable(),
      ProblemDefinition.set_equations()

  - update Term.standalone_setup() for create_adof_conns()

- variables, equations:

  - update active DOF connectivities only when needed in Equations.time_update()

    - initialize adof_conns in Variables.__init__()

  - remove active argument of FieldVariable.get_dof_conn()
  - update Variable.get_primary(), .get_dual() for incomplete variables
  - selective equations evaluation - new ProblemDefinition.eval_equations()
  - update Equations.evaluate() - split weak branch, fix/update docstring
  - add verbosity argument to FieldVariable.evaluate_at()

- terms:

  - update LinearTractionTerm: new eval mode

- regions:

  - fix passing region kind argument in Region.from_facets()
  - add parent argument to Region.from_facets()
  - report no entities error in Region.setup_from_highest()

- solvers:

  - update ProblemDefinition to provide and store default nls_status

    - update .__init__(), .init_solvers() to use .set_solvers_instances()
    - update .set_solvers_instances() to set .nls_status
    - update .solve() to use nls_status and set .nls_status

  - fix make_implicit_step() to use nls_status argument in all time steps
  - fix time step decrease in adapt_time_step() - update current time
  - update VariableTimeStepper.set_time_step() - new update_time argument

- mechanics:

  - new get_t4_from_t2s() helper function
  - clean up, update tests/test_tensors.py

- homogenization:

  - clean up sfepy/homogenization/engine.py
  - clean up homogen.py

- large deformation:

  - prevent false memory errors after warp violation in dq_finite_strain()
  - fix false memory errors after warp violation in
    dq_tl_finite_strain_surface()
  - fix error reporting of finite strain family functions

- mesh generation:

  - vectorize gen_block_mesh(), new get_tensor_product_cells()

- scripts:

  - new script/plot_times.py - taken from plant cell mechanics project

- miscellaneous updates:

  - fix transform_asm_vectors()
  - re-add '2_2', '3_2' to vtk_cell_types
  - remove unused filename argument of IntegralMeanValueOperator.__init__()
  - save LCBC data only in the initial time step of LCBC application
  - fix ProblemDefinition.save_ebc() for no equations (--solve-not option)
  - move sfepy/physics/radial_mesh.py into specialized repository
  - new sfepy/base/mem_usage.py - get_mem_usage(), print_mem_usage()
  - new ProblemConf.update_conf(), .add_missing()
  - fix link in docstring of ProblemDefinition.update_materials()
  - update parsing lexical elements of problem configuration - new list_dict()
  - fix approx_order description in Field.__init__() docstring
  - update QuadraturePoints.from_table() docstring
  - fix IndexError in Mesh.from_region() / detected on Win8
  - update plot_dofs.py for 1D
  - update ConstantFunctionByRegion.get_constants() to use Region.get_cells()
  - update SplineBox.__init__(), new coors_shape attribute

- tests and examples:

  - fix typo in examples/piezo_elasticity/piezo.py
  - fix docstring of examples/navier_stokes/stokes_slip_bc.py
  - update test_install.py to test mesh generation/conversion scripts
  - update test_install.py to test probe.py, extractor.py

- general clean up:

  - remove kill_* shell scripts
  - remove unused Variable.get_full_state()
  - remove .hgignore
  - remove sfepy/mesh/femlab.py, sfepy/mesh/meshutils.py
  - remove friction_slip.py, sfepy/mechanics/friction.py
  - remove isfepy, sfepy/interactive/ - update Viewer docstring, docs
  - remove plotPerfusionCoefs.py
  - remove script/config.py
  - remove sfepy/base/progressbar.py - update sfepy/ for no progress bar

- split sfepy.fem - separate FEM-specific and general modules:

  - merge branch 'split_fem'
  - rename sfepy/fem/ -> sfepy/discrete/fem/
  - move general files from sfepy/discrete/fem/ to sfepy/discrete/
  - update setup.py files for new paths
  - split fem/dof_info.py -> common/dof_info.py, fem/lcbc_operators.py
  - update sfepy/discrete/__init__.py, sfepy/discrete/fem/__init__.py
  - update sfepy/ for new paths
  - update top level scripts for new paths
  - update examples for new paths
  - update tests for new paths
  - update scripts for new paths
  - update docs for new paths
  - update docstrings for new paths
  - update MANIFEST.in for new paths
  - require integral argument in Approximation.describe_geometry()
  - make Region.create_mapping() a standalone function create_mapping()

    - update get_physical_qps()

  - fix Region.iter_cells() for subdomains
  - add cell_offsets attribute to Domain, use it in Region

    - simplify Domain.get_cell_offsets()

  - move fem/region.py -> common/region.py
  - new Domain.get_evaluate_cache()
  - remove mesh argument from probe constructors, remove all mesh references

    - update Probe.probe() to use Domain.get_evaluate_cache()
    - simplify Probe.__init__()

  - move discrete/fem/probes.py -> discrete/probes.py
  - update examples for moved and updated probes

- rename ProblemDefinition -> Problem
- rename many modules (fix according to naming conventions):

  - genPerMesh.py -> script/tile_periodic_mesh.py
  - findSurf.py -> script/extract_surface.py
  - script/evalForms.py -> script/eval_ns_forms.py
  - runTests.py -> run_tests.py
  - sfepy/optimize/shapeOptim.py -> sfepy/optimize/shape_optim.py
  - sfepy/optimize/freeFormDef.py -> sfepy/optimize/free_form_def.py
  - termsAcoustic.{py, c, h} -> terms_acoustic.{py, c, h}
  - termsAdjointNavierStokes.{py, c, h} -> terms_adj_navier_stokes.{py, c, h}
  - termsBasic.{py, c, h} -> terms_basic.{py, c, h}
  - termsBiot.{py, c, h} -> terms_biot.{py, c, h}
  - termsElectric.{py, c, h} -> terms_electric.{py, c, h}
  - termsLaplace.{py, c, h} -> terms_diffusion.{py, c, h}
  - termsLinElasticity.{py, c, h} -> terms_elastic.{py, c, h}
  - termsNavierStokes.{py, c, h} -> terms_navier_stokes.{py, c, h}
  - termsPiezo.{py, c, h} -> terms_piezo.{py, c, h}
  - termsPoint.py -> terms_point.py
  - termsSurface.{py, c, h} -> terms_surface.{py, c, h}
  - termsVolume.{py, c, h} -> terms_volume.{py, c, h}
  - termsHyperElasticity.{c, h} -> terms_hyperelastic.{c, h}
  - formSDCC.{c, h} -> form_sdcc.{c, h}
  - tests/testsBasic.py -> tests/tests_basic.py

- docs:

  - fix gallery page title
  - new download page with the downloads counter
  - add debugging section to installation.rst
  - update development tab
  - fix and update description of defining material parameters by functions
  - delete findSurf.py related text
  - sync module index of developer guide with current sources

.. _2013.3-2013.4:

from 2013.3 to 2013.4
=====================

- fix nodal polynomial spaces for orders >= 3:

  - merge branch 'fix_nodal_basis' - closes #205
  - fix get_facet_dof_permutations(), H1NodalMixin._setup_facet_orientations()

    - permutations should come from facet nodes in reference order

  - increase approximation orders of 3_4, 3_8 in tests/test_poly_spaces.py
  - prevent building ori_map for higher order elements in FESurface.__init__()

- fields:

  - reimplement get_facet_dof_permutations() using simple iterators
  - update Field.from_conf() to set field kind from region kind

    - update surface field syntax in examples

- simplify integrals - remove integral kind:

  - update transform_integrals() - remove/ignore integral kind
  - clean up sfepy/fem/integrals.py, add module docstring
  - remove kind attribute of Integral
  - remove Integral.get_key()
  - update Term for no integral kind, remove .get_integral_info()

    - remove PointTermBase

  - update sfepy.fem modules for no integral kind
  - allow 'i' as integral name in equations
  - update examples for no integral kind
  - update tests for no integral kind
  - docs: update users guide and tutorial for no integral kind

- mesh, domain, regions:

  - add tdim argument to mesh_set_coors()
  - update CMesh to use max. topological dimension - new tdim attribute
  - update Region to use max. topological dimension - new tdim attribute
  - add tdim attribute to Domain.shape
  - update fields for max. topological dimension

- start support of 'plate' integration/connectivity type

  - merge branch 'plate_mapping'
  - update CMapping - new mtx_t attribute
  - update Approximation.describe_geometry() for plate integration
  - update get_physical_qps() - new map_kind argument
  - update Term.get_physical_qps() for plate integration
  - update get_shape_kind(), Term.iter_groups() for plate integration
  - update FieldVariable.get_data_shape(), .evaluate() for plate integration
  - update Approximation.get_connectivity() for plate integration
  - new transform_asm_vectors(), transform_asm_matrices()
  - update TLMembraneTerm.weak_function() - use transform functions

- equation sequence solver:

  - new Variables.get_dual_names()
  - clean up sfepy/fem/equations.py, add module docstring
  - new Equations.get_variable_dependencies()
  - new sfepy/base/resolve_deps.py - functions for resolving dependencies

    - new get_nums(), solvable(), remove_known(), try_block(), resolve()

  - update tests/test_base.py: new test_resolve_deps()
  - new Equations.create_subequations()
  - new ProblemDefinition.create_subproblem()
  - new EquationSequenceSolver for stationary problems with an equation sequence
  - new examples/thermo_elasticity/thermo_elasticity_ess.py + test

- scripts:

  - script/save_basis.py: fix mesh mode, make permutations directly from mesh
  - script/plot_mesh.py: fix help message
  - simple.py: --list option can list solvers, new print_solvers()
  - new script/plot_quadratures.py - closes #178

    - update QuadraturePoints.from_table() to store true order
    - new sfepy/postprocess/plot_quadrature.py

  - new script/show_terms_use.py
  - clean up script/eval_tl_forms.py, update for current sympy
  - clean up script/evalForms.py, update for current sympy

- terms:

  - new ContactSphereTerm (dw_contact_sphere)

    - rename sfepy/mechanics/contact_planes.py
      -> sfepy/mechanics/contact_bodies.py
    - new ContactSphere - elastic contact sphere

  - new ContactPlaneTerm._get_force_pars(), update .get_fargs()
  - fix empty array condition in ContactPlaneTerm.smooth_f()
  - fix docstring of VolumeSurfaceTerm
  - remove Term.has_integral, .itype, .raw_itype attributes (unused)
  - fix form_tlcc_strainGreen_VS() - Green strain computation for
    post-processing - out-of diagonal entries are no more doubled

- homogenization:

  - update AcousticMassTensor, AppliedLoadTensor zero division check, fix eye
    call
  - update MiniAppBase.init_solvers(): allow to use different solvers for the
    correctors

- miscellaneous updates:

  - new sfepy/linalg/check_derivatives.py - check_fx(), check_vfvx()
  - update HypermeshAsciiMeshIO reader: support for quad4 and tri3 cells
  - update mesh reading/writing: print infromation about supported cell types
  - generalize dict_to_array() for nD arrays
  - fix type errors under Windows: in nm.take() cells
  - fix dict_to_array() for non-array sequences
  - update QuadraturePoints.from_table() to store true order
  - simplify ProblemDefinition.__init__() - re-implements #212
  - clean up sfepy/fem/parseEq.py
  - clean up sfepy/fem/variables.py, add module docstring

- examples and tests:

  - new examples/linear_elasticity/elastic_contact_sphere.py + test
  - update docstrings of poisson.py, poisson_short_syntax.py examples

    - add mutual references, tutorial reference
    - remove comments in examples/diffusion/poisson.py as their text is in the
      tutorial

  - update docstring of poisson_periodic_boundary_condition.py example
  - update docstring of poisson_field_dependent_material.py example
  - new examples/diffusion/poisson_neumann.py + test - closes #195

- docs:

  - fix issue, closes #236: move downloads from http://code.google.com/p/sfepy/
    to sfepy.org/downloads
  - fix duplicated module index in html
  - remove link to GoogleWiki
  - add labels to tutorial

.. _2013.2-2013.3:

from 2013.2 to 2013.3
=====================

- new implementation of Mesh topology data structures in C:

  - merge branch 'cmesh'
  - rename mesh.pyx -> cmesh.pyx
  - new C mesh data structures:

    - new mesh.c, mesh.h
    - new MeshGeometry, MeshConnectivity, MeshTopology, Mesh, Region, Indices
    - mesh construction, printing
    - new MeshEntity, MeshEntityIterator + related functions
    - connectivity computation algorithms: mesh_transpose(), mesh_intersect(),
      mesh_build(), face and edge orientations

  - new CConnectivity, CMesh - cython wrappers of C mesh structures
  - new functions for getting entities incident to given entity:
    new Indices, me_get_incident(), me_get_incident2(), contains()
  - new functions for setting up local entities:

    - Mesh and CMesh: new .entities attribute, new LocalEntities struct

  - support DEBUG_MESH flag in sfepy/fem/extmods/
  - new debprintf(), use it in mesh connectivity setup functions
  - new tests/test_cmesh.py - test_cmesh_counts()
  - miscellaneous updates:

    - clean up common_python.c, common.h, fix size_t format specifiers
    - new mem_list_new(), mem_list_remove(), mem_check_ptr(), update
      mem_alloc_mem(), mem_free_mem()
    - new cmem_statistics() to access mem_statistics() from Python
    - new mem_get_cur_usage(), mem_get_max_usage(), get_cmem_usage() C functions
    - new get_cmem_usage()

- new implementation of regions based on CMesh, CMesh improvements
  (w.r.t. above):

  - merge branch 'regions'
  - Region (summary):

    - new attributes: kind, true_kind, entities, can, can_vertices, can_edges,
      can_faces, parent
    - removed attributes: true_cells, all_vertices
    - use CMesh
    - vertices, edges, faces, facets and cells are properties pointing to
      .entities attribute, automatically created when needed

  - refactor Region, part 1

    - use CMesh
    - new Region.setup_from_highest(), .setup_from_vertices()
    - vertices, edges, faces, facets and cells are properties pointing to
      .entities attribute, automatically created when needed

  - refactor Region, part 2 - update operators

    - remove old operators
    - new .eval_op_vertices(), .eval_op_edges(), .eval_op_faces(),
      .eval_op_facets(), .eval_op_cells()
    - simplify region_op()

  - refactor Region, part 3 - update special constructors

    - update .from_vertices(), update & rename .from_faces() -> .from_facets()

  - refactor Region, part 4 - group compatibility functions

    - new .igs property
    - update .update_shape(), .get_vertices_of_cells(), .get_vertices(),
      .get_edges(), .get_faces(), .get_surface_entities() -> .get_facets(),
      .get_cells()
    - remove .delete_zero_faces() body, .select_cells(),
      .select_cells_of_surface(), .has_cells_if_can()
    - update docstring

  - refactor Region, part 5 - all_vertices -> vertices

    - call region.update_shape() in Domain.create_region()
    - remove Region.select_cells_of_surface(), .setup_face_indices() calls
    - update FESurface.__init__()
    - update VolumeField._setup_geometry()
    - update Region.create_mapping(), .has_cells()

  - Region (other details):

    - new Region.__can, .__facet_kinds class attributes, add kind to
      .__init__() arguments
    - update Region.__init__() - new Region.set_kind()
    - replace Region.setup_face_indices() with new Region.get_facet_indices()
    - add group offset in Region.get_cells() by default
    - update Region.setup_mirror_region() - use parent region information,
      update term igs in Term.get_conn_info()
    - new Region.get_entities()
    - update Region: force creation of region entities to respect region kind

      - new Region._access()

    - warn about empty group indices in Region.igs()
    - update Region.iter_cells(), .get_charfun()
    - update Region.get_edge_graph()

  - C mesh:

    - new Mask struct
    - new mesh_select_complete() C function, new CMesh.get_complete()
    - new CMesh.get_surface_facets()
    - new mei_init_sub() C function - connectivity subset iterator
    - new CMesh.get_incident()- new mesh_count_incident(), mesh_get_incident()
      C functions
    - new CMesh.get_conn_as_graph()
    - new CMesh.get_from_cell_group(), .get_igs() adapter functions

      - for compatibility until new assembling is done
      - new CMesh.cell_groups attribute, update CMesh.from_mesh()

    - fix uninitialized memory error in mesh_build() - manifested on 32 bit
      platforms
    - update CMesh.get_incident() to optionally return offsets, update
      mesh_get_incident()
    - new mesh_get_local_ids() C function, CMesh.get_local_ids()
    - new CMesh.get_cell_coors(), mesh_get_cell_coors() C function
    - fix facet orientation computation in mesh_build()
    - fix array sizes in CMesh.setup_entities()
    - new CMesh.get_orientations()
    - improve debprintf() messages
    - new ind_print() C function
    - fix mei_init_conn(), conn_alloc() for no incident entities
    - check for zero length entities argument in CMesh - update
      .get_incident(), .get_local_ids(), .get_complete(), .get_igs()
    - fix facet orientation computation in mesh_build(), try 2: wrong last index
    - fix return type of CMesh.get_surface_facets()

  - Domain and Mesh:

    - rename sfepy/fem/parseReg.py -> sfepy/fem/parse_regions.py
    - update create_bnf() for new region selection operators
    - update Domain.create_region(), region_leaf(), region_op()

      - construct CMesh instance in Domain.__init__() - new Domain.cmesh
        attribute

    - remove old testing code in sfepy/fem/parse_regions.py
    - rename nodes -> vertices, elements -> cells in region selectors - update
      create_bnf(), region_leaf()
    - update transform_regions() - replace flags with kind, allow string only,
      update Domain.create_regions()
    - new Domain.get_cell_coors() - use in 'cells by function' region selector
    - remove Facets-related methods from Domain - remove .setup_facets(),
      .get_facets(), .surface_faces()
    - generalize mesh_get_cell_coors() -> mesh_get_centroids()

      - generalize CMesh.get_cell_coors() -> .get_centroids(),
        Domain.get_cell_coors() -> .get_centroids()

    - update transform_regions() for parent field, update
      Domain.create_region()
    - set region kind and definition in Domain.create_region()
    - call CMesh.setup_entities() in Domain.__init__()
    - remove Facets and most of sfepy/fem/facets.py contents, rename & update
      _build_orientation_map() -> build_orientation_map()
    - simplify FESurface.__init__() - use prepare_remap()
    - update FESurface.setup_mirror_connectivity() for CMesh, create .ori_map
      in .__init__()
    - new get_facet_dof_permutations(), permute_facets(),
      get_dof_orientation_maps() - updated and recombined code from removed
      Facets
    - remove for loop in FESurface.setup_mirror_connectivity()
    - update Domain.refine(), refine_2_3(), refine_2_4(), refine_3_4(),
      refine_3_8()
    - make FESurface.fis C-contiguous
    - update Mesh.from_region(), ._append_region_faces()

  - fields:

    - update H1NodalMixin._setup_facet_orientations()
    - update H1NodalMixin._setup_facet_dofs() - update ._setup_edge_dofs(),
      ._setup_face_dofs()
    - update Field._get_facet_dofs(), .get_dofs_in_region_group()
    - update H1HierarchicVolumeField._setup_facet_dofs() - update
      ._setup_edge_dofs(), ._setup_face_dofs()
    - remove unused facet_oris argument of H1NodalMixin._setup_facet_dofs()

  - terms:

    - adapt data types in ConnInfo.iter_igs(), Term.get_assembling_cells()
    - fix Term.get_physical_qps() for point integration
    - update NewTerm.integrate()

  - miscellaneous updates:

    - update define_box_regions()
    - update eval_in_els_and_qp()
    - update ConstantFunctionByRegion.get_constants()
    - new PhysicalQPs.__init__(), update get_physical_qps()
    - fix edit_filename() to use prefix
    - allow points and edges in MeditMeshIO, VTKMeshIO
    - fix Viewer.save_image() for file names with path
    - update smooth_mesh()

  - examples and tests:

    - update all examples and tests for new region syntax
    - update tests/test_poly_spaces.py
    - update tests/test_regions.py: add more selectors, number regions
    - update test_install.py
    - new examples/navier_stokes/stokes_slip_bc.py + test

  - scripts:

    - findSurf.py: update and clean up

      - new _get_facets(), get_surface_faces()
      - improve help message

    - script/gen_gallery.py:

      - fix figure paths in generate_images()
      - add custom view for stokes_slip_bc.py example

    - new script/plot_mesh.py - basic CMesh plotting

      - new sfepy/postprocess/plot_cmesh.py: new plot_wireframe()
      - new plot_entities(), label_global_entities(), label_local_entities()

  - configuration:

    - clean up sfepy/base/conf.py
    - fix key in transform_regions()

  - docs:

    - update primer and tutorial
    - update users guide

- new MultiProblem solver: allows "conjugate" solution of subproblems

  - merge branch 'subpb' - closes #226
  - new example that uses MultiProblem solver:

    - new examples/acoustics/vibro_acoustic3d.py,
      examples/acoustics/vibro_acoustic3d_mid.py
    - new mesh files meshes/2d/acoustic_wg_mid.vtk, meshes/3d/acoustic_wg.vtk

  - new tests/test_input_vibro_acoustic3d.py

- postprocessing and visualization:

  - update Viewer.build_mlab_pipeline() to display data in only_names order
  - new plot_points()

- homogenization:

  - update CoefNN, CoefN and CoefSym ("smart" set_variables_default); new
    CoefExprPar
  - fix CoefSym and CopyData coefficients

- terms:

  - fix SDLinearElasticTerm, new SDSufaceIntegrateTerm, remove
    SDSufaceNormalDotTerm
  - new dw_convect_v_grad_s (ConvectVGradSTerm, scalar gradient term with
    convective velocity)
  - update (adjoint) Navier-Stokes terms for 2D, update convect_build_vtg,
    convect_build_vtbg() for 2D
  - update LinearTractionTerm: material parameter is now optional (merge branch
    'ltr_opt_mat' - closes #228)
  - clean up: remove PermeabilityRTerm
  - support 'qp' evaluation mode in HyperElasticBase - update .integrate(),
    .get_fargs(), .get_eval_shape()

- miscellaneous updates:

  - update Mesh.localize() to remove vertices not used in any cell
  - update Probe.probe() to report non-finite probe points
  - fix ProblemDefinition.set_equations_instance() for keep_solvers=True
  - check topological dimension of cells in VolumeField._check_region()
  - update ensure_path() to raise error on failure
  - update HypermeshAsciiMeshIO: use component id as material id
  - fix tiled_mesh1d()
  - new dets_fast()
  - update Region.igs to return parent igs if parent exists

    - update get_dependency_graph() to include parent region in dependencies
    - raise proper error when dependency is missing

- examples and tests:

  - update poisson_functions.py example to use dw_convect_v_grad_s
  - new examples/navier_stokes/navier_stokes2d.py + test
  - update tests/test_poly_spaces.py for 2_3, 3_4 geometries, new 2_3_2z.mesh,
    3_4_2z.mesh
  - new tests/test_mesh_smoothing.py - mesh smoothing test
  - clean up and update tests/test_parsing.py for new regions

- scripts:

  - new script/gen_mesh_prev.py: generates mesh previews in a given directory
  - script/save_basis.py:

    - make sure that output directory exists
    - enable saving all geometry element permutations - new
      save_basis_on_mesh()

.. _2013.1-2013.2:

from 2013.1 to 2013.2
=====================

- automatic testing of term calls:

  - merge branch 'i154' - closes #154
  - define arg_shapes attributes of terms
  - use new get_arg_kinds() in Term.classify_args() - simplify code

    - new Term._check_variables()

  - update Term.classify_args() to support 'mode' attribute
  - update Term to support arg_shapes and geometries attributes
  - fix SurfaceLaplaceLayerTerm, SurfaceCoupleLayerTerm
  - fix dw_ul_volume() - fixes VolumeULTerm, BulkPressureULTerm
  - docs: update developer guide - describe new Term attributes
  - fix and update DeformationGradientTerm

    - rename dq_def_grad -> ev_def_grad
    - remove DeformationGradientTerm.__call__(), new .function(), .get_fargs(),
      .get_eval_shape()
    - update dq_def_grad()

  - fix dw_v_dot_grad_s_sw() (VectorDotGradScalarTerm)
  - new fmf_mulATC()
  - fix ScalarDotGradIScalarTerm - add vector mode
  - fix docstrings of VectorDotGradScalarTerm, TLMembraneTerm
  - fix SufaceNormalDotTerm - make material shape (D, 1) to be compatible with
    SDSufaceNormalDotTerm
  - fix Term.call_function() to clear C errors properly
  - fix SDSufaceNormalDotTerm
  - update Term.get() to accept integration argument
  - fix NSOFSurfMinDPressDiffTerm, DiffusionTerm, DiffusionCoupling,
    PiezoCouplingTerm
  - fix docstrings of SDLinearElasticTerm, StokesTerm, GradTerm, DivTerm
  - new tests/test_term_call_modes.py to automatically test (almost) all terms,
    skips terms with an empty arg_shapes attribute

- new elastic contact plane term:

  - new flag_points_in_polygon2d()
  - new ContactPlane - contact plane with polygonal boundary

    - new sfepy/mechanics/contact_planes.py
    - new plot_polygon(), plot_points()

  - new ContactPlaneTerm (dw_contact_plane)
  - add custom view for elastic_contact_planes.py example to
    script/gen_gallery.py

- terms:

  - rename DotTermsDotProductSurfaceTerm functions:

    - dw_surface_dot_vectornormscalar -> dw_surface_s_v_dot_n
    - dw_surface_dot_scalarnormvector -> dw_surface_v_dot_n_s

  - update DotProductSurfaceTerm: new 'scalar_norm_vector' mode

- postproc.py:

  - force re-read in VTKFileSource.set_filename()
  - add reload button to ViewerGUI, support watching of vtk files in Viewer

    - new ReloadSource
    - closes #217

  - fix default opacity in Viewer

- input-output:

  - update VTKMeshIO.write() to mark unfinished write by 'x' in 1. byte
  - new HDF5MeshIO.read_bounding_box()

- bases:

  - merge branch 'cbases'
  - move low-level function from bases.pyx to new lagrange.c

    - new sfepy/fem/extmods/lagrange.[ch]
    - get_barycentric_coors(), get_xi_simplex(), get_xi_tensor(),
      eval_lagrange_simplex(), eval_lagrange_tensor_product() ('_' prefix from
      name removed) translated from Cython to C
    - update bases.pyx, use fmf_pretend_nc(), import FMField directly

  - rename & update script/gen_lobatto_pyx.py -> script/gen_lobatto1d_c.py

    - generate lobatto1d.c and lobatto1h.c
    - new lobatto1d_template.c, lobatto1d_template.h

  - split lobatto_template.pyx into lobatto_bases.pyx and lobatto.c, low level
    functions in lobatto.c, as in lagrange.c

- miscellaneous updates:

  - clean up of many modules
  - new fmf_pretend_nc()
  - fix dependencies in sfepy/terms/extmods/setup.py - this fixes rebuilding
    terms.so even when an unrelated file was changed
  - fix Struct._str() for classes with overloaded .get() (e.g. Term)
  - fix Domain.save_regions_as_groups() for regions without true cells
  - fix iterative solvers for tolerances not given in conf - fix PyAMGSolver,
    PETScKrylovSolver, PETScParallelKrylovSolver
  - remove compatibility function sorted()
  - remove unused attribute of CharacteristicFunction
  - add f_tol tolerance option to ScipyBroyden nonlinear solver
  - add custom norms to RadialMesh

- tests and examples:

  - new examples/homogenization/perfusion_micro.py + test

    - homogenization of a double-porous medium
    - new 3D mesh: matrix with two disconnected channels

  - new examples/linear_elasticity/elastic_contact_planes.py + test -
    demonstrating use of contact plane term

- docs:

  - add development tab, new doc/development.rst
  - link wiki from development tab
  - add related projects section
  - new sections on term evaluation modes and postprocessing/probing

    - move description of term evaluation modes from developer to users guide
    - update users guide
    - closes #196

- gallery:

  - improve gen_gallery.py

    - add captions, contents and section titles
    - move _gallery_template contents to new doc/gallery_template.html

  - fix script/gen_gallery.py for new time stepping
  - fix output suffix for time-dependent problems in generate_images()
  - update docstring of gen_gallery.py - describe docs regeneration steps

.. _2012.4-2013.1:

from 2012.4 to 2013.1
=====================

- solvers:

  - move time stepping solvers to new sfepy/solvers/ts_solvers.py
  - redesign time stepping, unify use of stationary and evolutionary solvers:

    - new StationarySolver
    - update PDESolverApp.call()
    - update TimeSteppingSolver - change __init__(), __call__() arguments,
      remove .set_step_fun()
    - remove solve_stationary(), solve_evolutionary()
    - move prepare_matrix(), prepare_save_data(), make_implicit_step(),
      make_explicit_step() into sfepy/solvers/ts_solvers.py
    - new get_initial_state()
    - update SimpleTimeSteppingSolver.__call__() to implement fully the time
      stepping loop, new .solve_step()
    - update ExplicitTimeSteppingSolver, new .solve_step

  - new AdaptiveTimeSteppingSolver - implicit adaptive time stepping solver:

    - new get_min_dt(), adapt_time_step()
    - new VariableTimeStepper.from_conf(), .set_from_ts(),
      .set_n_digit_from_min_dt()
    - update VariableTimeStepper.set_step() to allow only step = 0
    - update VariableTimeStepper.__iter__() to include t1

- input-output:

  - fix ANSYSCDBMeshIO.read(), guess it with .inp suffix:

    - new .guess()
    - for mixed element meshes and coordinates without "solid" keyword

  - allow quadratic elements in ANSYSCDBMeshIO.read(), strip extra nodes:

    - allow more format variants
    - add remap argument to mesh_from_groups()
    - add .dat suffix

  - update ANSYSCDBMeshIO.read() to support nodes of boundary conditions
  - support node sets (nodes of boundary conditions) in HDF5MeshIO.read(),
    .write()

  - fix omitting nodes of boundary conditions in Mesh.write()
  - add nodal_bcs argument to Mesh.from_data()
  - update HDF5MeshIO.read() for old meshes without node_sets group

- mesh, domain, regions:

  - support mat_ids in merge_mesh()
  - new Mesh.__add__() to merge meshes
  - add verbose argument to gen_block_mesh(), gen_cylinder_mesh()
  - new gen_extended_block_mesh() mesh generator
  - fix smooth_mesh() - improve efficiency
  - support nodes and elements of set in create_bnf()
  - update region_leaf() and test for "nodes of set" selector:

    - prepare for "elements of set" selector
    - nodes of group support only int id

  - move bodies of functions for saving regions from ProblemDefinition to
    Domain:

    - new Domain.save_regions(), .save_regions_as_groups()
    - update ProblemDefinition.save_regions(), .save_regions_as_groups()

  - new Region.delete_zero_faces(), used in .complete_description()
  - use Region.get_cells() instead of direct access to cells

    - true cells are checked in Region.get_cells() if needed
    - true_cells attribute is properly initialized in
      Region.complete_description()

  - allow dot in set names in create_bnf()
  - allow dot in region names
  - fix Region.get_edge_graph() for degenerate edges

- fields, variables:

  - update setting of variables data (use Variable.set_data())
  - rename/update Variables.non_state_data_from_state() ->
    .set_data_from_state()
  - rename/update Variable.data_from_any() -> .set_data()
  - rename FieldVariable.data_from_qp() -> .set_data_from_qp()
  - new Variable.is_finite()

- new terms:

  - dw_tl_bulk_active (active bulk pressure term)

- miscellaneous updates:

  - make view button to print view and roll as arguments of postproc.py
  - set ts directly in ProblemDefinition.update_time_stepper()
  - add verbose argument to MyBar progress bar
  - allow to run python shell from debugger
  - more grammar elements available in parse_conf.py
  - radial mesh - fixes, integrals and derivatives
  - new get_face_areas() + test
  - new global option 'check_term_finiteness'
  - check finiteness in Term.evaluate()
  - fix compute_nodal_normals() for degenerate elements, check for zero normals
  - new configure_output()
  - remove many unused functions, code clean up

- scripts:

  - script/save_basis.py: plot nodes of selected dofs, new plot_nodes()
  - script/convert_mesh.py: new --center option

- examples and tests:

  - update linear_elastic_damping.py example to use adaptive time stepping
  - new tests/test_mesh_generators.py - test basic mesh generators

- docs:

  - add support section to main page

.. _2012.3-2012.4:

from 2012.3 to 2012.4
=====================

- initial implementation of elements with hierarchical basis:

  - merge hierarchic branch
  - fields and approximations:

    - split fields.py into fields_base.py, fields_nodal.py
    - new VolumeField, H1NodalMixin
    - new H1NodalVolumeField, H1DiscontinuousField, H1NodalSurfaceField
    - prepend underscore to non-API methods
    - new H1HierarchicVolumeField (fields_hierarchic.py)
    - move ._setup_esurface(), .create_mesh() into Field
    - new Field.from_args()
    - move setting DOFs from a function to a field method

      - new H1NodalMixin.set_dofs(), update EquationMap.map_equations()

    - add merge argument to Field.get_dofs_in_region_group()
    - new H1HierarchicVolumeField.set_dofs() - hack for constant values
    - define edge orientation (ap.ori) in H1HierarchicVolumeField
    - new code for face orientations

      - update H1HierarchicVolumeField._setup_facet_dofs()
      - new LobattoTensorProductPolySpace._get_counts(), ._get_face_axes_nodes()
      - update LobattoTensorProductPolySpace.__init__(), _define_nodes(),
        ._eval_base()

    - update Approximation for facet orientation, add ori attribute
    - new Approximation.clear_qp_base()
    - add iels argument to Approximation.get_base()
    - new Approximation.get_poly_space(), update .get_base()

  - bases and mappings:

    - update CVolumeMapping.describe() for basis per element
    - update VolumeMapping.get_mapping() for facet orientation argument
    - update PolySpace.eval_base() for facet orientation argument
    - support edge orientation in LobattoTensorProductPolySpace._eval_base()
    - add node_orders attribute to LobattoTensorProductPolySpace
    - add force_axis argument to PolySpace.eval_base()
    - use force_axis argument to fix Lagrange basis computations
    - fix gradient mode in eval_lobatto_tensor_product()
    - ensure c-contiguity in PolySpace.eval_base()
    - update PolySpace.any_from_args() for no simplex Lobatto class
    - raise error for surface mapping with hierarchical basis

  - facets and geometry elements:

    - new Facets.get_orientation()
    - fix axes orientation of 6th face of 3_8
    - new GeometryElement.get_grid()

      - new _get_grid_1_2(), _get_grid_2_3(), _get_grid_2_4(), _get_grid_3_4(),
        _get_grid_3_8()

  - linearization of higher order fields:

    - update linearizer.py for hierarchic basis

      - remove bf argument of DOF evaluation functions (use ps instead)
      - update get_eval_dofs() - add ps (poly. space) and ori (facet
        orientation) arguments
      - update create_output()
      - update Field.linearize(), eval_dofs() for updated linearizer.py

    - create linearized output for an expression

      - new create_expression_output()
      - new get_eval_expression(), eval_in_els_and_qp()

  - scripts:

    - update script/plot_condition_numbers.py
    - update script/save_basis.py for meshes (global basis)

      - new --lin-options, --dofs, --permutations, --plot-dofs options
      - support gradients in mesh mode
      - require output directory, print help by default

  - examples and tests:

    - new tests/test_normals.py - check surface normal orientations
    - new tests/test_poly_spaces.py - check hierarchic and nodal basis
      continuity
    - add two-element meshes
    - update examples/diffusion/sinbc.py to use hierarchical basis

  - miscellaneous updates:

    - new prepare_translate()
    - new GeometryElement.get_conn_permutations()
    - new make_h1_projection_data()
    - generalize sym2dim()
    - new sfepy/postprocess/plot_dofs.py

      - new _get_axes(), plot_mesh(), plot_global_dofs(), plot_local_dofs()

    - update filename attribute in UserMeshIO.read()

- unify C/Cython structures for volume and surface mappings:

  - merge surfvol branch
  - rename geometry.[ch] -> refmaps.[ch]
  - merge VolumeGeometry and SurfaceGeometry into Mapping
  - merge CVolumeMapping and CSurfaceMapping into CMapping
  - update sfepy/fem/extmods/setup.py
  - update map_print(), small fixes in map_describe(), map_integrate()
  - update terms for Mapping
  - update VolumeMapping, SurfaceMapping for CMapping
  - miscellaneous updates for CMapping
  - remove eval_real_extra(), dq_grad_extra()
  - remove dw_surface_dot_vector(), dw_surface_dot_scalar()

    - update dw_volume_dot_vector(), dw_volume_dot_scalar(),
      DotProductVolumeTerm
    - update docstring of DotProductSurfaceTerm

  - remove remaining references to CMapping.area
  - new CauchyStrainSTerm (ev_cauchy_strain_s)

- problem descriptions:

  - replace getattr() calls with ProblemConf.get_function()

    - fix for functions defined inside define() and not on the problem
      description module level

  - add allow_tuple, free_word arguments to create_bnf()
  - allow deferred setup of ProblemConf - add setup argument to .__init__()
  - update ProblemDefinition.setup_default_output() for missing options

- data probing:

  - new Probe.set_n_point(), update .__init__()
  - allow extrapolation in probes

    - new Probe.set_options() - set close limit, update .__init__(), .probe()

  - probe.py: new --close-limit option, update generate_probes(), clean up
  - update generate_probes() to support multiple figures per probe hook
  - add size_hint to Probe.set_options(), fix refinement for multiple calls

    - new Probe.reset_refinement() - called in .probe()

  - new write_results(), read_results(), update generate_probes(),
    postprocess()
  - move read_header(), get_data_name() to sfepy/fem/probes.py, update
  - new dict_from_options(), update ProblemConf.from_file_and_options()

    - make dict_from_string a standalone function

  - generalize override in ProblemConf.__init__()
  - update generate_probes() to support -1 as last step

- base:

  - new edit_tuple_strings(), add recur argument to edit_dict_strings()
  - update edit_dict_strings() for lists of replacements and tuple values
  - update assert_() for custom messages
  - add overwrite_by_none argument to update_dict_recursively()
  - remove sfepy/base/tasks.py
  - new edit_filename() utility function
  - new Struct._str(), .str_class(), update .__str__()
  - simplify dict-like methods of Struct

    - replace Struct.get_default_attr() with updated Struct.get()
    - rename .set_default_attr() -> set_default()

- logging:

  - split sfepy/base/log.py - new sfepy/base/log_plotter.py, clean up

    - rename ProcessPlotter -> LogPlotter

  - implement reading and plotting of Log text output

    - new read_log(), plot_log()

  - add log_name argument to get_logging_conf(), update docstring

- materials:

  - update Materials, Material for optional clearing of material data
  - improve updating of material parameters

    - update Materials.time_update() - force and clear arguments replaced by
      mode argument
    - update Material.__init__(), .update_special_data(),
      .update_special_constant_data(), .time_update()

- phononic materials:

  - new get_log_freqs(), update detect_band_gaps()
  - new get_ranges(), update cut_freq_range(), BandGaps.__call__()

    - opts.eig_range is a slice now

  - update transform_plot_data() to work with nan values

- linear combination boundary conditions:

  - rename _save_normals() -> _save_vectors(), update
  - new Region.get_edge_graph()
  - new get_edge_paths(), _get_edge_path()
  - new compute_nodal_edge_dirs()
  - new NormalDirectionOperator.get_vectors(), update .__init__()
  - new EdgeDirectionOperator, new 'edge_direction' LCBC

    - update LCBCOperators.add_from_bc()

  - fix and simplify LCBCOperator.treat_pbcs()

    - LCBCs combined with EPBCs in common nodes work correctly now,
      in particular for non-penetration conditions with general normals

- scripts:

  - remove obsolete scripts:

    - mesh_to_vtk.py, hfm3_mesh.py, convert.py, edit_neu.py, make_spkg.py,
      spymatrix.py

  - new script/sync_module_docs.py
  - script/save_basis.py: fix --plot-dofs for no --permutations, update help
  - update script/gen_gallery.py for custom visualizations
  - script/plot_condition_numbers.py: report number of quadrature points

- simple.py:

  - clean up, add docstring
  - new --save-ebc-nodes option, change meaning of --save-ebc

    - update solve_pde(), save_only(), PDESolverApp.call()
    - update generate_images(), test_hyperelastic_tlul.py, tests_basic.py

- postproc.py:

  - group options, clean up, use docstring for help
  - allow negative time steps in Viewer.call_mlab()

    - Python-like indexing from the last step
    - check that step is in the range available

  - update step help message

- solvers:

  - fix arguments in TimeStepper.from_conf()
  - fix ScipyFMinSolver.__call__() for old scipy (0.7.2)
  - add divergence tolerance option to PETSc solvers

- miscellaneous updates:

  - fix examples in script/config.py
  - new sfepy/postprocess/plot_facets.py

    - new plot_geometry(), plot_edges(), plot_faces(), draw_arrow()

  - new get_perpendiculars()
  - move .setup_coors() to Field, fix vertex indices
  - new Region.get_vertices_of_cells(), update .update_vertices()
  - fix VolumeField._setup_vertex_dofs() for superfluous vertices
  - fix VolumeField.average_qp_to_vertices() for superfluous vertices
  - new RadialVector, ExplicitRadialMesh, update RadialHyperbolicMesh,
    RadialMesh
  - new H1NodalSurfaceField.interp_v_vals_to_n_vals() - initial stub
  - fix ProblemDefinition.save_ebc() for EBCs needing ProblemDefinition
    instance
  - report (initial) data ranges in Viewer.build_mlab_pipeline()
  - allow context variables in ProblemDefinition.evaluate()

    - add strip_variables argument to .evaluate(), .create_evaluable()

  - rename State.get_scaled_norm() -> get_weighted_norm()

    - fix weighting by state parts, add weights, return_weights arguments

  - fix EquationMap._init_empty() to initialize n_epbc
  - clean up sfepy/fem/dof_info.py
  - fix gen_mesh_from_voxels()
  - remove unused function: reorder_dofs_on_mirror
  - fix integration over mirror surface
  - treat separately real and imaginary fill value parts in Field.extend_dofs()

    - new get_min_value() - fill value is 0 for vectors, minimum for scalars

  - fix dw_jump term, use the (optional) material parameter as a multiplier

- examples and tests:

  - update docstrings in examples/diffusion/poisson_functions.py
  - report number of failed tests in test_install.py
  - new example with field-dependent diffusion coefficient + test
    (examples/diffusion/poisson_field_dependent_material.py)
  - new diffusion example with periodic boundary conditions + test
    (examples/diffusion/poisson_periodic_boundary_condition.py)
  - new tests for field dependent and periodic diffusion examples
  - new acoustic 3d example (examples/acoustics/acoustics3d.py) + test
  - remove examples/diffusion/subdomains.py + test

- docs:

  - include scripts in developer guide
  - sync module index of developer guide with current sources

.. _2012.2-2012.3:

from 2012.2 to 2012.3
=====================

- terms:

  - fix BulkPressureULTerm and CompressibilityULTerm
  - move TLMembraneTerm.describe_membrane_geometry() into new
    describe_geometry()
  - fix TLMembraneTerm.eval_function() for updated transform_data()
  - new FibresActiveTLTerm.get_eval_shape()
  - update VolumeTLTerm.get_fargs() for 'el_avg' mode
  - fix term evaluation to clear properly the C error flag
  - new Term.call_get_fargs(), .call_function()
  - update Term.eval_real(), .eval_complex(), .evaluate()
  - new SurfaceNormalDotTerm (dw_surface_ndot)
  - new ScalarDotGradIScalarTerm (dw_s_dot_grad_i_s)
  - update DotProductSurfaceTerm: allow (vector * normal * scalar) product
  - new shape sensitivity terms: SDLinearElasticTerm (d_sd_lin_elastic),
    SDSufaceNormalDotTerm (d_sd_surface_ndot)
  - support 'el' mode in term evaluation:

    - update SurfaceFluxTerm for 'el' mode, update docstring

  - fix checking of term - integral compatibility

- merge pull request #186 from vlukes/meshtiler:

  - new tiled_mesh() function
  - move mesh generating functions to sfepy/mesh/mesh_generators.py
  - update genPerMesh.py
  - update scripts cylindergen.py and blockgen.py

- visualization:

  - postproc.py: allow finer control of opacity (closes Issue 194)
  - fix _get_scalars() for vectors

- base:

  - move try_imports() into sfepy.base.base, fix error printing, add docstring
  - update invert_dict() - add unique argument, add docstring
  - fix find_subclasses() for no name attribute case
  - add name_attr argument to find_subclasses(), load_classes()
  - update Output to accept a file descriptor instead of a file name
  - update Output for global options, new test_verbose_output()
  - support global options (closes Issue 190) - new sfepy/base/goptions.py:
    - new ValidatedDict class, goptions instance
    - default_goptions with 'verbose' option

- homogenization:

  - fix 'return_all' option in HomogenizationApp
  - fix CoefSymSym

- fem:

  - new get_ref_coors(), update Field.evaluate_at()
  - new make_l2_projection_data(), update make_l2_projection()
  - support 'space' and 'poly_space_base' items in transform_fields()
  - fix Region.delete_groups() - update also all_vertices etc.
  - remove unused qp_indx field in PhysicalQPs
  - fix stiffness_from_lame_mixed()
  - material parameters can be defined per region using region names:

    - update linear homogenization examples using new material definition
    - fix Material.__init__() for special and keyword values

  - new SplineBox geometry parametrization
  - update extend_cell_data() for surface regions
  - allow different base function values per element:

    - add base function attribute (bf) to VolumeGeometry, SurfaceGeometry
    - make bf of CVolumeMapping, CSurfaceMapping point to the above bf
    - add flag argument to {CVolumeMapping, CSurfaceMapping}.__cinit__()
    - update Approximation.describe_geometry()
    - new FMF_SetCellX1() macro
    - update terms for base functions per element

  - define integrals by order instead of quadrature name:

    - update transform_integrals(), Integrals.from_conf(), Integral.__init__()
    - update tests and examples for new integral definition

- input-output:

  - boundary condition id in Nastran format is used as the node group
  - update MeditMeshIO.read() to skip unsupported entities

- solvers:

  - update imports in PysparseEigenvalueSolver.__call__() for new Pysparse
  - update ScipyIterative solver to allow preconditioning
  - check solver consistency in ProblemDefinition.set_solvers_instances()
  - fix preconditioner arguments for qmr method in ScipyIterative
  - update ScipyIterative solver for iteration callback
  - new standard_call() decorator for linear solvers, checks also data shapes
  - update docstring of Newton (configuration options, small fixes)
  - update docstring of Solver (common solver configuration options)

- misc:

  - fix unique_rows() for non-contiguous arrays, add arguments of numpy.unique()

- examples and tests:

  - update ULF hyperelastic examples, change error tolerance, number of
    iterations
  - update hyperelastic tests: test |TL - UL| and |UL - UL_mixed|

- docs:

  - update for new integral definition
  - include special methods of nonlinear solvers
  - update User's Guide (Problem description file/Materials)
  - add notes on solving PDEs by FEM in the tutorial part of the documentation

.. _2012.1-2012.2:

from 2012.1 to 2012.2
=====================

- reimplement acoustic band gaps code using the homogenization engine:

  - merge band_gaps_he branch
  - rename eigen.py -> phonon.py
  - move AcousticBandGapsApp to new sfepy/homogenization/band_gaps_app.py
  - fix coding style in sfepy/homogenization/band_gaps_app.py
  - replace compute_density_volume_info() by DensityVolumeInfo
  - new MiniAppBase.process_options()
  - preserve order of requirements in HomogenizationEngine.call()
  - new SimpleEVP mini-application:

    - reimplement part of AcousticBandGapsApp.solve_eigen_problem() for
      'simple' problems

  - rename sfepy/homogenization/phono.py -> .../coefs_phononic.py
  - new Eigenmomenta mini-application:

    - update compute_eigenmomenta()
    - remove prepare_eigenmomenta()

  - new CoefDummy
  - update AcousticBandGapsApp.process_options(), .process_options_pv() -
    options moved to corrector/coefficient classes
  - update AcousticBandGapsApp.call() - compute ingredients for band gap
    detection via HomogenizationEngine
  - update frequency-dependent tensors:

    - update AcousticMassTensor, AcousticMassLiquidTensor, AppliedLoadTensor

  - update Eigenmomenta
  - update band gaps functions, new BandGaps mini-application:

    - update cut_freq_range(), split_chunks(), detect_band_gaps(),
      get_callback(), find_zero(), describe_gaps()
    - remove setup_band_gaps()

  - update try_set_defaults() for one-level recurrence
  - move plotting functions: coefs_phononic.py -> band_gaps_app.py
  - split SimpleEVP.__call__():

    - new SimpleEVP.prepare_matrices(), .post_process()

  - new SchurEVP mini-application
  - remove obsolete code in AcousticBandGapsApp:

    - remove make_save_hook()
    - remove AcousticBandGapsApp.fix_eig_range(), .solve_eigen_problem(),
      .eval_homogenized_coefs()

  - update BandGaps for custom detection function and rhs matrix
  - set incident wave dir. in AcousticBandGapsApp.call()
  - new dispersion mini-applications:

    - new ChristoffelAcousticTensor, PolarizationAngles, PhaseVelocity
      mini-applications
    - new compute_cat_sym_sym(), compute_cat_dim_sym(), compute_cat_dim_dim(),
    - remove compute_cat(), compute_polarization_angles()
    - new AcousticBandGapsApp.plot_dispersion()
    - remove AcousticBandGapsApp.compute_cat(), .compute_phase_velocity()
    - remove report_iw_cat()

  - new HomogenizationEngine option 'compute_only':

    - unify handling of dependencies (requirements and coefficients)
    - remove _get_parents()

  - fix de-duplication of names in HomogenizationEngine.call()
  - inherit AcousticBandGapsApp from HomogenizationApp (common options)
  - update Coefficients._save_dict():

    - support custom 'to_file_txt' callable attribute of coefficients
    - ignore unknown types

  - update saving of figures in AcousticBandGapsApp
  - allow saving band gaps logs:

    - new BandGaps.to_file_txt(), .save_log()
    - new application option 'log_save_name'

  - improve default plot resources
  - remove broken caching from AcousticBandGapsApp
  - update phononic examples for AcousticBandGapsApp:

    - new examples/phononic/band_gaps_conf.py
    - update examples/phononic/band_gaps.py,
      examples/phononic/band_gaps_rigid.py
    - remove coef_conf_elastic.py, gen_mesh.py, parametric.py, plot_gaps.py

  - update test_install.py

- homogenization:

  - new homogenization corrector function OnesDim - unit vector in the
    directions
  - proper coding style in coefs_base.py
  - update float formatting capabilities of Coefficients.to_file_latex()
  - allow passing custom application options to HomogenizationEngine
  - new get_volume_from_options(), update HomogenizationApp.call()
  - allow to add some additional information to coefficients files

- quadratures:

  - update table of 1_2 quadrature points

    - fix and update module docstring
    - update QuadraturePoints.__init__() - add symmetric argument
    - up to polynomial order 47

  - increase tolerance in test_quadratures()
  - create missing tensor product quadratures using 1_2 (line) quadratures:
    - order for tensor product geometries is given in terms of the 1D order
    - new QuadraturePoints.from_table(), _get_max_orders()
    - update QuadraturePoints.__init__()
    - update Integral.get_qp(), .integrate()
    - move & update Integral.get_actual_order() -> get_actual_order()

  - fix order in create_mass_matrix(), make_l2_projection()
  - test generated tensor product quadratures
  - cite PHAML in module docstring
  - update table of 2_3 quadrature points
  - fix polynomial order in get_poly()
  - report differences in test_weight_consistency(), test_quadratures()
  - create missing simplex (2_3, 3_4) quadratures

    - update QuadraturePoints.from_table()
    - new sfepy/fem/simplex_cubature.py

  - update test_quadratures()

- solvers:

  - rewrite stabilization material handling

    - create_stabil_mat() -> StabilizationFunction
    - stabilization material uses now a "regular" material function

  - update Oseen solver for new stabilization material handling
  - fix eig() to obey options from keyword arguments

    - change num -> n_eigs

  - fix PysparseEigenvalueSolver.__call__() for no r.h.s. matrix

    - set default n_eigs, clean up

  - new SimpleTimeSteppingSolver.process_conf()
  - new TimeSteppingSolver.set_step_fun()
  - new ExplicitTimeSteppingSolver
  - new MassOperator

- applications:

  - rename simple_app.py -> pde_solver_app.py
  - rename SimpleApp -> PDESolverApp
  - rename pde_solve() -> solve_pde(), add docstring

    - move it into sfepy/applications/pde_solver_app.py
    - remove sfepy/applications/top_level.py

  - clean up sfepy/applications/application.py
  - merge sfepy/solvers/generic.py with sfepy/applications/pde_solver_app.py

    - remove sfepy/solvers/generic.py
    - update solve_pde(), PDESolverApp.call()
    - rename solve_stationary_op() -> solve_stationary(),
      time_step_function() -> make_implicit_step(),
      solve_evolutionary_op() -> solve_evolutionary()

  - update PDESolverApp.call() for basic explicit time stepping

    - new make_explicit_step()

- input-output:

  - update Mesh._set_shape_info() - add element group dimensions
  - update MeshIO classes to read elements of lower dimension

    - update VTKMeshIO, ComsolMeshIO, AVSUCDMeshIO, HypermeshAsciiMeshIO,
      AbaqusMeshIO, NEUMeshIO, ANSYSCDBMeshIO
    - merge mesh_from_tetra_hexa(), mesh_from_tri_quad() -> mesh_from_groups()

  - change default value of omit_facets to False (Mesh.from_file() etc.)
  - update HDF5MeshIO.read_times() to return also time steps
  - update extract_times(), dump_to_vtk(), extractor.py for missing time steps

    - extraction now works with files where not all the time steps are saved

  - update VTKMeshIO.read() for pixels and voxels

- domain:

  - update Domain: store group dimension

    - update Domain.fix_element_orientation() to skip facet groups

  - update Facets for empty facet groups
  - update Region.get_n_cells() for facet groups

- scripts:

  - new script/plot_condition_numbers.py
  - update script/gen_lobatto_pyx.py to generate also derivatives
  - fix gen_lobatto() in script/gen_lobatto_pyx.py
  - new script/save_basis.py
  - clean up script/blockgen.py, new --2d option
  - clean up script/convert_mesh.py, new --refine option

- visualization:

  - use new _get_scalars() in domain specific plot functions
  - new plot_warp_scalar() domain specific plot function
  - linearizer:  update create_output() for custom evaluation functions

    - new get_eval_dofs(), get_eval_coors()
    - add coordinates argument to DOFs evaluation function in get_eval_dofs()
    - update create_output()

- new LobattoTensorProductPolySpace - initial implementation of hierarchic basis

- unify dot product and mass terms:

  - dw_mass, dw_mass_scalar, dw_surf_mass_scalar, dw_volume_wdot_scalar ->
    dw_volume_dot_vector, dw_surface_dot_vector, dw_volume_dot_scalar,
    dw_surface_dot_scalar
  - remove sfepy/terms/termsMass.py, sfepy/terms/extmods/termsMass.[ch]
  - remove MassVectorTerm, MassScalarTerm, MassScalarSurfaceTerm,
  - new sfepy/terms/terms_dot.py, sfepy/terms/extmods/terms_dot.[ch]
  - update examples

- terms:

  - update AdjDivGradTerm, AdjConvect1Term, AdjConvect2Term for new term
    evaluation
  - update objective function terms for new term evaluation

    - update SUPGCAdjStabilizationTerm (fix name from AdjSUPGCtabilizationTerm)
    - update SUPGPAdj1StabilizationTerm, SUPGPAdj2StabilizationTerm

  - update adjoint Navier-Stokes stabilization terms for new term evaluation

    - merge NSOFMinGrad1Term, NSOFMinGrad2Term -> NSOFMinGradTerm
    - update NSOFSurfMinDPressTerm, NSOFSurfMinDPressDiffTerm

  - update Navier-Stokes shape derivative terms for new term evaluation

    - update SDDivTerm, SDDivGradTerm, SDConvectTerm

  - update SDDotScalarTerm for new term evaluation (was TestPQTerm)

  - update stabilized Navier-Stokes shape derivative terms for new term
    evaluation

    - update SDGradDivStabilizationTerm, SDSUPGCStabilizationTerm,
      SDPSPGCStabilizationTerm, SDPSPGPStabilizationTerm
    - fix their docstrings

  - make material argument of DivGradTerm optional, implement evaluation mode
  - update PermeabilityRTerm for new term evaluation
  - update SDDotVolumeTerm to support also vectors (was SDDotScalarTerm)
  - new VectorDotGradScalarTerm (dw_v_dot_grad_s)
  - update DotProductVolumeTerm for matrix coefficient
  - add optional material argument to NonPenetrationTerm
  - fix term_mode in d_sd_... terms - set default value to 1

- examples:

  - update examples/navier_stokes/stabilized_navier_stokes.py
  - new examples/diffusion/time_poisson_explicit.py

- shaper.py:

  - fix --direct, --adjoint options
  - update solve_stokes(), solve_navier_stokes(), solve_generic_direct(),
    solve_direct(), solve_adjoint() for State
  - fix main()
  - update and clean up
  - update vec -> state in shape optimization functions and ShapeOptimFlowCase
  - update ShapeOptimFlowCase: handling of materials, term evaluation mode
  - update update_mesh(), solve_problem_for_design() for current code

- misc:

  - fix Equations.invalidate_term_caches() to invalidate all variables
  - further speed-up assemble.pyx by using pointer to iels, change int -> int32
  - major speed-up of some key functions

    - profiling using line_profiler -> use numpy.take() instead of fancy
      indexing
    - update Approximation.get_connectivity(), .describe_geometry()
      Field.setup_vertex_dofs(), .setup_coors(), Mapping.__init__(),
      create_adof_conn()

  - fix vg_getElementDiameters()
  - fix bf_actt_c1()
  - fix print_matrix_diff() for scipy with removed rowcol()
  - change default value of follow_epbc to False

    - update Equations.strip_state_vector(), State.get_reduced(),
      Variables.strip_state_vector(), FieldVariable.get_reduced()
  - new get_subdict()
  - new match_coors()
  - fix region comparison in Field.get_mapping()
  - fix package_check() for bogus versions
  - update Field.linearize()
  - new eval_lobatto_tensor_product()
  - make vertex_maps common, fix LagrangeTensorProductPolySpace._define_nodes()
  - update ProblemDefinition.get_integrals() - add names, kind arguments
  - remove unused MultiplierVariable, ConstantVariable

- docs:

  - add archlinux installation requirements
  - fix numpy docstring standard links
  - add coding style section to developer guide
  - refer to web pages using new doc/links.inc
  - add NTC logo to main page
  - update developer guide to reflect the above updates

.. _2011.4-2012.1:

from 2011.4 to 2012.1
=====================

- initial version of linearizer of higher order solutions:

  - merge linearizer branch
  - new Field.linearize()
  - clean up SimpleApp.process_options()
  - update output creation functions/methods
  - fix Field.get_true_order() for forced bubble DOFs
  - new Field.is_higher_order()
  - new application option 'linearization'
  - new FieldVariable.linearize()
  - silence Domain.setup_facets()
  - HDF5:

    - do not linearize when saving to 'h5' format
    - update HDF5MeshIO.write() to save field name in 'full' mode
    - update HDF5MeshIO.read_data()
    - allow 'strip' mode for 'h5' in ProblemDefinition.setup_output()
    - fix HDF5MeshIO.read_time_stepper() to close file on exception

  - move output creation to Field, new Field.create_output()

    - remove FieldVariable.extend_dofs(), .remove_extra_dofs(), .linearize()

  - update recover_bones(), recover_paraflow()
  - extractor.py: new --linearization option

    - support linearization when dumping to VTK
    - new create_problem(), parse_linearization()
    - catch ValueError in dump_to_vtk()

  - new Struct.update()
  - add linearization options to examples/navier_stokes/navier_stokes.py,
    examples/diffusion/sinbc.py
  - new tests/test_linearization.py
  - docs: update developer guide, new linearizer.rst

- solvers:

  - new PETScParallelKrylovSolver solver class, petsc_worker.py
  - add precond_side option to PETSc Krylov solvers
  - prevent unwanted importing of PETSc
  - allow changing properties of solvers via command line
  - fix PETSc Krylov solvers to obey nonzero initial guess
  - fix Newton, Oseen to pass current solution as initial guess to linear solver
  - add verbose option to Solver (base of all solvers)
  - use verbose option in Newton solver
  - new SchurGeneralized solver, update SchurComplement
  - new ScipyFMinSolver optimization solver
  - fix ScipyBroyden for non-array return value

- postprocessing and visualization:

  - allow reusing viewer scene in make_animation()
  - build mlab.pipeline after scene is activated

    - reason: some VTK objects or properties require a scene with a camera
      and interaction to be open to work properly

  - new plot_velocity() domain specific plot function
  - force source update in Viewer.render_scene()
  - always use Viewer source change hack
  - fix error reporting in Viewer.call_mlab()
  - catch only ImportError in mayavi imports
  - fix Viewer.build_mlab_pipeline() for no point scalars
  - postproc.py: new --fgcolor, --bgcolor options
  - fix default colors in Viewer.call_mlab()
  - fix mayavi imports in postprocess.utils

- homogenization:

  - new CopyData corrector
  - remove unused homog. coefficient functions
  - fix HomogenizationApp.call()

- input-output:

  - support tensors in vertices in VTKMeshIO.write()
  - fix supported_capabilities (nastran supports write())
  - update read_array() for unknown number of columns
  - fix NEUMeshIO.read() - problem with boundary conditions
  - new MeshIO.for_format(), update MeshIO.any_from_filename()
    - update script/convert_mesh.py
    - move output_writable_meshes() to sfepy/fem/meshio.py
  - blockgen.py, cylindergen.py: new --format option

- problem description:

  - allow extra arguments to define() in ProblemConf.from_file()
  - simple.py: new --define option

- schroedinger:

  - remove specialized DFT code (move to a separate project)
  - move code to new sfepy/physics/schroedinger_app.py
  - add result names arguments to SchroedingerApp.save_results()

- tests and examples:

  - update its2D_1.py example - use refine_mesh()
  - update Primer examples to use dw_point_load term
  - new examples/linear_elasticity/linear_viscoelastic.py + test
  - new examples/thermo_elasticity/thermo_elasticity.py + test
  - update examples/linear_elasticity/linear_viscoelastic.py

    - save in HDF5 format
    - new post_process() - compute strain, stresses
    - update linear_tension()
    - new main() - plot time histories

  - new examples/linear_elasticity/prestress_fibres.py + test

- material coefficients:

  - fix bulk_modulus_lame()
  - clean up sfepy/mechanics/matcoefs.py, add/update docstrings
  - improve testing of material parameter conversion functions
  - rename test_tensors() -> test_elastic_constants()
  - new test_conversion_functions(), test_stiffness_tensors()
  - use better names for material parameter conversion functions

    - youngpoisson_to_lame() -> lame_from_youngpoisson()
    - stiffness_tensor_*() -> stiffness_from_*()
    - bulk_modulus_lame() -> bulk_from_lame()
    - bulk_modulus_youngpoisson() -> bulk_from_youngpoisson()

- variable evaluation:

  - new find_ref_coors(), evaluate_in_rc() instead of evaluate_at()
  - update Field.evaluate_at() to use find_ref_coors(), evaluate_in_rc()
  - allow caching reference coordinates, add ret_ref_coors argument
  - change cache attribute ctree -> kdtree
  - print more detailed timings
  - update FieldVariable.evaluate_at()

- term evaluation:

  - fix Term.evaluate() to pass diff_var argument to .eval_*()
  - add step, time_derivative arguments to Term.get()
  - support preserving evaluate caches in State and Variable
  - fix ProblemDefinition.solve() to preserve evaluate caches
  - rewrite variable and evaluate cache history handling

    - Variable.data attribute is now a deque with length given by history
    - history is given by integer >= 0, not string
    - evaluate cache is split by time steps
      (evaluate_cache(mode) -> step_cache(step) -> cache(key))

  - update examples for new history syntax
  - remove sfepy/terms/cache.py, sfepy/terms/cachesHistory.py
  - remove code related to term caches
  - allow preserving evaluate caches in ProblemDefinition.evaluate()

- terms:

  - optional material parameter in LaplaceTerm
  - fix LaplaceTerm for test_msm_symbolic.py
  - fix arguments of dw_permeability_r()
  - new PointTermBase, update LinearPointSpringTerm
  - new concentrated point load term (dw_point_load)
  - fix acoustic terms
  - fix BulkPressureULTerm, DotProductVolumeTerm
  - new term assembling function mulATB_integrate()
  - update DiffusionCoupling term - simplify assembling
  - remove d_volume_dot term (subset of dw_volume_dot, same class name)
  - update DotProductVolumeTerm (dw_volume_dot)
  - fix TLMembraneTerm.describe_membrane_geometry()
  - fix parents of LinearElasticIsotropicTerm
  - update HyperElasticBase.get_family_data() for new evaluate cache handling
  - new THTerm, ETHTerm classes, support custom advance in Variable.advance()
  - update LinearElasticTHTerm, LinearElasticETHTerm for new term evaluation
  - fix and clean up compute_mean_decay()
  - update BiotTHTerm, BiotETHTerm for new term evaluation
  - update DotSProductVolumeOperatorWTHTerm, DotSProductVolumeOperatorWETHTerm
    for new term evaluation
  - update dw_volume_wdot_scalar()
  - update BiotStressTerm, BiotStressQTerm for new term evaluation
  - fix ConcentratedPointLoadTerm.check_shapes() for functions, update docstring
  - unify Cauchy strain, stress evaluation terms

    - dq_cauchy_strain, de_cauchy_strain -> ev_cauchy_strain
    - dq_cauchy_stress, de_cauchy_stress -> ev_cauchy_stress

  - unify Biot stress evaluation terms

    - de_biot_stress, dq_biot_stress -> ev_biot_stress

  - simplify Cauchy and Biot stress term evaluation in qp mode
  - unify diffusion velocity evaluation terms

    - de_diffusion_velocity, di_diffusion_integrate -> ev_diffusion_velocity

  - unify gradient and divergence evaluation terms

    - de_grad, dq_grad -> ev_grad
    - de_div, dq_div -> ev_div

  - rename di_surface_integrate -> ev_surface_integrate, update docstring
  - rename di_volume_integrate -> ev_volume_integrate, allow optional material
  - rename di_integrate_mat -> ev_integrate_mat, allow 'qp' mode
  - update CauchyStressTerm.function() for optional coefficient
  - new fading memory stress terms (ev_cauchy_stress_th, ev_cauchy_stress_eth)
  - simplify code of DiffusionVelocityTerm.function()
  - fix dw_volume_dot evaluation
  - update LinearPrestressTerm, LinearStrainFiberTerm for new term evaluation
  - remove obsolete sfepy/terms/terms_base.py, update docs

- logging:

  - update Log.plot_vlines() to add line also to text log file
  - move plotting-related imports in log.py into ProcessPlotter and Log
  - update Log.__init__(), .plot_vlines() to allow better plot reconstruction

- interactive:

  - move code from __init__.py to session.py in sfepy/interactive/
  - update isfepy for ipython 0.12 (adapted code from current isympy)
  - update isfepy docstring

- setup:
  - install scripts, examples and tests along sources, do not install docs
  - update setup.py

- misc:

  - new refine_mesh()
  - add conn attribute to GeometryElement
  - fix evaluate_at() for piecewise-constant approximations
  - speed-up assemble.pyx by using more pointer arithmetics
  - fix complex term evaluation
  - fix testing for string instances for Python 3
  - fix implicit function declaration warnings in generated terms.c
  - update test_install.py for updated its2D_3.py example
  - fix LagrangeNodes.append_tp_faces()
  - fix problem with connectivity in mirror surface
  - allow user defined base function in variable evaluation
  - add 2d von Mises stress calculation
  - update get_debug() for ipython 0.12
  - new get_mapping_data(), update get_jacobian(), new get_normals()
  - new mesh smoothing function smooth_mesh()
  - update Probe for ctree -> kdtree
  - update create_evaluable() to obey verbose argument in weak mode
  - add verbose argument to Variables.time_update(), Equations.time_update()
  - update ProblemDefinition.copy() to setup output of new problem
  - fix create_evaluable() to pass time stepper to terms
  - update Variable.init_data() to use .data_from_any()
  - remove unused evaluate()
  - cleanup in sfepy/geom, mesh and geometry tools moved to sfepy/mesh/
  - update term evaluation to use var_dict also with try_equations
  - update ProblemDefinition.create_evaluable()
  - fix section titles in docstrings
  - update formatting of term docstrings to remove sphinx warnings
  - many small docstring fixes to remove sphinx warnings

- docs:

  - update script/gen_gallery.py:

    - omit __init__.py files, unify link names
    - skip missing images

  - update many docstrings (especially all docstrings of terms to remove sphinx
    warnings)
  - update nodal stress description in Primer, mention refinement
  - update ubuntu installation section
  - update Primer for dw_point_load term
  - update for missing/removed/renamed modules
  - update users guide for new postproc.py options
  - cleanup whitespace, add table of contents to users guide
  - add basic info about solvers to users guide
  - fix and update Primer
  - update gitwash, include it to developer guide table of contents
  - update variables section in users guide
  - move release notes to doc/, add sphinx section labels
  - add news and archived news
  - update main html page, add links
  - add google analytics support
  - update release tasks
  - update developer guide for ev_* terms
  - remove termsLinElasticity_full.rst from developer guide
  - process_terms(): insert term call signature before arguments lists

.. _2011.3-2011.4:

from 2011.3 to 2011.4
=====================

- use cython instead of swig - merge cython branch:

  - wrap vxc() and vectorize it - new get_vxc():

    - update SchroedingerApp.iterate() to use get_vxc()
    - remove sfepy/physics/extmods/dft.i, update setup.py, __init__.py

  - new sfepy/fem/extmods/assemble.pyx:

    - new assemble_vector(), assemble_vector_complex()
    - new assemble_matrix(), assemble_matrix_complex()
    - use new assembling functions in Term.assemble_to()
    - remove old C assembling functions
    - update insert_sparse_to_csr() for assemble.pyx

  - remove unused functions in fem.c:

    - remove lagrange1*(), baseBiL*(), baseTriL*(),
      rezidual(), matrix(), inverse_element_mapping()

  - tests:

    - split tests/test_assembling.py

      - move test_eval_matrix(), test_vector_matrix(), test_surface_evaluate(),
        test_dq_de() into tests/test_term_consistency.py
      - move test_save_ebc() into new tests/test_ebcs.py
      - clean up and modernize

    - add new assembling tests to tests/test_assembling.py

      - new test_assemble_vector(), test_assemble_vector_complex(),
        test_assemble_matrix(), test_assemble_matrix_complex()

  - new sfepy/fem/extmods/bases.pyx:

    - new get_barycentric_coors(), _get_barycentric_coors()
    - new _eval_lagrange_simplex(), eval_lagrange_simplex()
    - new _eval_lagrange_tensor_product(), eval_lagrange_tensor_product()
    - new _get_xi_simplex(), _get_xi_tensor()
    - new evaluate_at()
    - update Field.evaluate_at() for new evaluate_at()

  - use new base evaluation functions:

    - update _eval_base() of LagrangeSimplexPolySpace,
      LagrangeSimplexBPolySpace, LagrangeTensorProductPolySpace
    - remove PolySpace.clear_c_errors()

  - new wrappers of orient_elements(), graph_components():

    - new sfepy/fem/extmods/mesh.pyx
    - add bases, mesh to __init__.py, remove meshutils
    - update Domain.fix_element_orientation()
    - update surface_components()
    - remove sfepy/fem/extmods/meshutils.i

  - new pyalloc(), pyfree() helper functions
  - move mesh_graph() and related functions to meshutils.c

    - remove raw_graph()
    - new create_mesh_graph()

  - define basic types in new types.pxd
  - new _fmfield.pyx, _fmfield.pxd - FMField-related wrappers:

    - new array2fmfield*()

  - remove fem module from setup.py
  - remove fem.c, fem.h, fem.i
  - remove geometry.i, update setup.py
  - new mappings.pyx:

    - new CVolumeMapping extension class
    - update setup.py for mappings.pyx, reuse _fmfield library
    - new mappings.pxd for mappings.pyx
    - new CSurfaceMapping extension class
    - use CVolumeMapping, CSurfaceMapping in VolumeMapping, SurfaceMapping

  - new Config.python_include()
  - use sfepy_common library for common C code
  - use sfepy_terms library for term C code
  - add ignore_errors argument to load_classes()
  - new array2pint2(), array2pint1() helper functions
  - simplify import of extmods.terms
  - new terms.pyx:

    - wrap dq_state_in_qp(), dq_grad(), dq_grad_extra(), dq_div_vector(),
      d_volume_surface(), di_surface_moment()
    - update eval_real(), eval_real_extra(), eval_complex()
    - wrap dq_finite_strain_tl(), dq_finite_strain_ul(),
      dq_tl_finite_strain_surface(), dq_tl_he_stress_bulk(),
      dq_ul_he_stress_bulk(), dq_tl_he_stress_neohook(),
      dq_ul_he_stress_neohook(), dq_tl_he_stress_mooney_rivlin(),
      dq_ul_he_stress_mooney_rivlin(), dq_tl_he_tan_mod_bulk(),
      dq_ul_he_tan_mod_bulk(), dq_tl_he_tan_mod_neohook(),
      dq_ul_he_tan_mod_neohook(), dq_tl_he_tan_mod_mooney_rivlin(),
      dq_ul_he_tan_mod_mooney_rivlin(), dw_he_rtm(), de_he_rtm(),
      dq_tl_stress_bulk_pressure(), dq_ul_stress_bulk_pressure(),
      dq_tl_tan_mod_bulk_pressure_u(), dq_ul_tan_mod_bulk_pressure_u(),
      dw_tl_volume(), dw_ul_volume(), dw_tl_diffusion(),
      dw_tl_surface_traction(), dq_def_grad(), he_residuum_from_mtx(),
      he_eval_from_mtx()
    - update hyperelastic term classes for terms.pyx:

      - update HyperElasticBase, HyperElastic{TL, UL}Base,
        SurfaceTractionTLTerm

    - wrap dw_volume_wdot_scalar(), dw_laplace(), d_laplace(),
      dw_diffusion(), d_diffusion(), dw_permeability_r(),
      dw_diffusion_coupling(), d_diffusion_coupling(),
      de_diffusion_velocity(), d_surface_flux()
    - wrap dw_lin_elastic_iso(), dw_lin_elastic(), d_lin_elastic(),
      dw_lin_prestress(), dw_lin_strain_fib(), de_cauchy_strain(),
      de_cauchy_stress(), dq_cauchy_strain(), dw_surface_ltr(),
      dw_volume_lvf(), dw_mass(), dw_mass_scalar(), d_mass_scalar(),
      dw_surf_mass_scalar()
    - wrap term_ns_asm_div_grad(), term_ns_asm_convect(),
      dw_lin_convect(), dw_div(), dw_grad(), dw_st_pspg_c(),
      dw_st_supg_p(), dw_st_supg_c(), dw_st_grad_div(), dw_biot_grad(),
      dw_biot_div(), d_biot_div(), dw_piezo_coupling(),
      d_piezo_coupling(), dw_electric_source()
    - wrap d_diffusion_sa(), dw_surf_laplace(), d_surf_laplace(),
      dw_surf_lcouple(), d_surf_lcouple()
    - add stubs for remaining functions
    - remove terms.i

  - new crcm.pyx:

    - rcm(), permute_in_place()
    - update setup.py and init files for crcm.pyx
    - remove rcm.i, array.i, common.i, fmfield.i

  - update .gitignore
  - update docs (swig -> cython)
  - fix script/gen_term_table.py for circular dependency
  - build:

    - update Clean.run() to remove cython-generated files
    - fix speed regression (numpy.distutils quirk)
    - make cython compulsory dependence in setup.py
    - set min. cython version to 0.14.1

- make proper class for physical quadrature points:

  - new PhysicalQPs, .get_merged_values(), .get_shape()
  - update get_physical_qps()
  - update Material for PhysicalQPs
  - update PhysicalQPs.get_shape() - make ig optional, check shape
    compatibility

    - also allow any raveled shape length >= 1

  - update get_physical_qps(), add n_total attribute to PhysicalQPs

- updated Lagrangian formulation:

  - fix evaluate.new_ulf_iteration()
  - update hyperelasticity (ULF) terms for new Term.evaluate()
  - remove cachesFiniteStrain.py:

    - functions moved to terms_hyperelastic_ul(tl).py

  - fix hyperelastic ULF terms: CompressibilityULTerm, VolumeULTerm
  - fix hyperelastic terms in ULF for mixed pressure-displacement formulation

- terms:

  - update DiffusionIntegrateTerm, renamed: d_diff... --> di_diff...
  - new arg_type: 'opt_material' - optional material

    - remove "_mass_scalar_w" terms, use "_mass_scalar" instead
    - remove "_w" terms in termsBasic.py and termNavierStokes.py

  - fix Term.classify_args()
  - fix SurfaceMomentTerm.get_eval_shape()
  - fix dw_volume_dot, SurfaceTerm, SumNodalValuesTerm, DiffusionSATerm
  - update DiffusionIntegrateTerm, DiffusionVelocityTerm,
    DiffusionRTerm, disable PermeabilityRTerm
  - update acoustic and diffusion terms

- problem description:

  - update ProblemConf.get_function() for passing functions directly
  - use ProblemConf.get_function() to get parametric hook
  - update ProblemConf.from_file() to accept define_args as tuple

- tensors:

  - update transform_data() for fourth order tensors
  - check numbers of points in transform_data()

    - do not meddle with data shape implicitly

  - fix transform_data() (second order case)
  - add mode argument to prepare_cylindrical_transform()
  - fix docstring of transform_data()

- polynomial spaces:

  - pass space, poly_space_base to Interpolant, SurfaceInterpolant:

    - update {SurfaceField, Field}.create_interpolant()

  - update PolySpace.eval_base() for array-like coors argument
  - new sfepy/fem/extmods/lobatto_template.pyx

    - new eval_lobatto()
    - update setup.py for lobatto.pyx
    - add generated sfepy/fem/extmods/lobatto.pyx to simplify building

- mesh refinement:

  - new Mesh.get_element_coors()
  - new refine_3_8()
  - update Domain.refine() for hexahedral meshes
  - new gen_misc_mesh(), gen_mesh_from_string()
  - new ProblemDefinition.refine_uniformly()
  - new refine_2_3(), refine_2_4()
  - update Domain.refine() for 2D meshes

- solvers:

  - simplify imports in ScipyIterative.__init__()
  - silence sparse efficiency warning
  - new linear solver - Schur complement

- input-output:

  - update MeditMeshIO.read() - add omit_facets argument
  - update Mesh.from_file(), MeshIO.read() - add omit_facets argument
  - remove *args from MeshIO.read(), .write()
  - new BDFMeshIO.write() function (Nastran data format)
  - fix VTKMeshIO.read()
  - clean up ioutils.py (imports, whitespace, coding style)
  - update supported_capabilities dict - prepare for boundary conditions
  - update Mesh._set_data() - add nodal_bcs argument
  - update NEUMeshIO.read() to read nodes of boundary conditions
  - update Mesh.from_file() - allow passing Mesh instance in place of file name

- fields, variables:

  - update Field.get_mapping() for initial time step
  - fix indx initialization in Variable.__init__()
  - update Field.setup_dof_conns() for string dc_type
  - report standalone vertices in evaluate_at()
  - check for standalone vertices in Field.evaluate_at()

- scripts:

  - turn script/gen_term_table.py into sphinx extension
  - new script/gen_gallery.py:

    - generate documentation and gallery of sfepy examples

  - new script/gen_lobatto_pyx.py
  - postproc.py:

    - update parse_view(), --view help message
    - new --no-offscreen option

  - add script/show_authors.py

- schroedinger:

  - more fixes for meshes with several element groups
  - new sfepy/physics/radial_mesh.py - RadialMesh, RadialHyperbolicMesh
  - fix for higher order approximations, clean up
  - update split of DFT and general options, update docstrings
  - add init_hook application option
  - add v_fun_name application option
  - update Potential for user arguments to evaluation function
  - new Potential.__len__()
  - rename --mesh option to -create-mesh
  - new --mesh, --mesh-dir options
  - use ensure_path(), allow only one of --mesh, --create-mesh

- homogenization:

  - update recover_micro_hook(), new parameter 'recovery_file_tag'
  - update homogenized coefficient classes for term_mode argument:

    - update MiniAppBase.__init__()
    - update {CoefSymSym, CoefFMSymSym, CoefDimSym, CoefNN, CoefN, CoefSym,
      CoefFMSym, CoefOne, CoefFMOne}.__call__()

  - new volume computation and selection
  - allow passing scalar volume to coefficient classes

    - new MiniAppBase._get_volume()

  - fix output key in CorrMiniApp.get_output()
  - update HomogenizationApp.process_options(), .call() - clean up
  - clean up Coefficients

- tests and examples:

  - new test_hyperelastic_tlul.py - compare TLF and ULF solutions
  - remove tests/test_input_hyperelastic.py

    - test duplicated by test_hyperelastic_tlul.py

  - add basic docstrings with equations to examples in gallery
  - fix linear_elastic_probes.py example
  - update test_install.py to test linear_elastic_probes.py example
  - add test for refine_3_8():

    - new refine(), test_refine_hexa()
    - rename test_refine() -> test_refine_tetra()

  - add test for fourth order tensors support in transform_data():

    - new test_transform_data4()

  - test mesh interpolation invariance (FieldVariable.evaluate_at()):

    - new gen_datas(), test_invariance(), test_invariance_qp()

  - add __init__.py files to allow package imports in examples/
  - fix importing utils.py in Navier-Stokes examples:

    - wrong file (from pytables) got imported in Windows

  - improve testing of uniform mesh refinement

- misc:

  - new get_jacobian() utility function
  - update set_mesh_coors() - add clear_all argument
  - fix argument types in gen_block_mesh(), gen_cylinder_mesh()
  - update FieldVariable.data_from_qp() for higher order approximations
  - fix test_install.py for windows
  - new ensure_path(), locate_files(), remove_files()
  - update Viewer.set_source_filename() - catch also AttributeError
  - fix syntax for python 2.5
  - remove FMField.stride (unused in cython wrappers)
  - update 'nodes of group' region selector for named groups
  - use NumPy C type definitions, fix pointer type and printf format warnings
  - remove unnecessary includes (fix _POSIX_C_SOURCE warnings)
  - remove unused vg_integrateChunk(), sg_integrateChunk()
  - remove caches attribute from Equations and Equation
  - improve getting default integrals in ProblemDefinition:

    - new ProblemDefinition.get_integrals()

  - fix common_python.c for 64bit systems
  - add reference element meshes
  - remove unused methods of GeometryElement:

    - remove .orient_edges(), .orient_faces3(), .orient_faces4()

  - remove sfepy/optimize/fluentutils.py
  - fix docstring of get_green_strain_sym3d() for LaTeX
  - fix centre argument in get_coors_in_tube()

- build:

  - update Clean.run() to clean also examples/, script/, tests/
  - include *.pxd files and lobatto_template.pyx in distribution
  - do not install version.h.in
  - add lobatto_template.pyx, version.h.in to MANIFEST.in

- documentation:

  - new examples.rst
  - update latex_preamble in doc/conf.py
  - add terms_constraints.rst to table of contents
  - new sfepy docs html layout
  - update region selection syntax in users guide
  - fix developer guide for renamed/deleted/new files
  - move installation docs into new doc/installation.rst, update
  - add remaining sfepy modules to developer guide
  - update script options in users guide
  - link examples, primer from index
  - add features to introduction
  - add applications to index
  - more on ebcs given by functions

.. _2011.2-2011.3:

from 2011.2 to 2011.3
=====================

- major update of terms - merge easier_terms branch:

  - aim: easier usage and definition of terms while retaining original C
    functions
  - approximations and fields:

    - new Approximation.get_connectivity()
    - update Approximation.describe_geometry():

      - new return_mapping argument
      - make volume mappings local to given region

    - include used polynomial space in geometry object:

      - update Approximation.describe_geometry()

    - manage mappings in Field:

      - new Field.clear_mappings(), .save_mappings(), .create_mapping(),
        .get_mapping()

  - equations and variables:

    - evaluate and cache quantities in FieldVariable:

      - new FieldVariable.evaluate(), .clear_evaluate_cache()

    - update Variable.__call__() history behaviour for step 0
    - update Variable.advance() to advance evaluate cache
    - copy equations variables in ProblemDefinition.create_evaluable()
    - invalidate evaluate cache in Variable.data_from_any()
    - factor out variable evaluation to new evaluate_variable.py:

      - update FieldVariable.evaluate()
      - new eval_real(), eval_complex()

    - update Equations.invalidate_term_caches() for evaluate caches:

      - new FieldVariable.invalidate_evaluate_cache()

  - move NewTerm.get_shape_kind() to get_shape_kind()
  - update Term:

    - new Term.get(), .get_mapping(), .get_data_shape()
    - refactor Term.evaluate():

      - new Term.check_shapes(), .eval_real(), Term.eval_complex(),
        split_complex_args()

    - new Term.get_assembling_cells()
    - remove Term.needs_local_chunk()
    - update Term.assemble_to()
    - fix Term.iter_groups() for no cells in a group

  - update many terms for new Term.evaluate():

    - dw_laplace, dw_diffusion, de_diffusion_velocity, dw_biot,
      dw_lin_elastic, de_cauchy_strain, de_cauchy_stress:

      - update both Python and C code, lots of simplifications

    - dw_mass_scalar, dw_mass_scalar_w
    - dw_lin_elastic_iso, dq_cauchy_strain, dq_cauchy_stress,
      dw_mass_vector, dw_piezo_coupling, dw_surface_ltr
    - remove dw_mass
    - dw_div_grad, dw_convect, dw_lin_convect, dw_stokes, dw_stokes_w
    - merge dq_lin_convect with dw_lin_convect
    - dq_grad, de_grad, dq_div, de_div (was d_div)
    - dw_point_lspring, dw_volume_lvf
    - di_volume_integrate (merged with de_average_variable,
      dq_state_in_volume_qp)
    - di_surface_integrate (was d_surface_integrate, merged with
      dq_state_in_surface_qp)
    - di_surface_integrate_w (was d_surface_integrate_w)
    - dw_volume_integrate, dw_volume_integrate_w,
      dw_surface_integrate, dw_surface_integrate_w, d_volume_dot,
      d_surface_dot, d_surface_dot_w, d_volume, d_surface, d_volume_surface,
      di_surface_moment, d_sum_vals
    - dw_volume_dot_w, allow different approximation of the arguments
    - di_integrate_mat (was de_volume_average_mat +
      di_volume_integrate_mat)
    - dw_jump, remove dw_jump()
    - dw_non_penetration term
    - dw_st_pspg_c, dw_st_supg_p, dw_st_supg_c, dw_st_grad_div
    - dw_surface_mass_scalar, dw_surface_mass_scalar_w, dw_bc_newton
    - remove dw_mass_scalar_fine_coarse
    - d_surface_flux:

      - rename d_hdpm_surfdvel -> d_surface_flux
      - remove termsHDPM.* files, move functions to termsLaplace.*

    - dw_tl_membrane
    - dw_electric_source
    - update basic hyperelastic TL terms for new Term.evaluate()

      - new HyperElasticBase.get_family_data()
      - new HyperElasticTLBase.integrate(), .function(), .compute_family_data(),
        .compute_stress(), .compute_tan_mod(), .get_fargs(), .get_eval_shape()
      - dw_tl_he_neohook, dw_tl_he_mooney_rivlin, dw_tl_bulk_penalty
      - dw_tl_fib_a, dw_tl_surface_traction

    - update TL perfusion terms for new Term.evaluate()

      - dw_tl_bulk_pressure, dw_tl_volume, dw_tl_diffusion terms

  - update examples for new Term.evaluate():

    - examples/biot/biot_npbc.py
    - examples/navier_stokes/navier_stokes.py
    - examples/linear_elasticity/material_nonlinearity.py
    - examples/biot/biot_npbc_lagrange.py
    - examples/homogenization/linear_elastic_mM.py

  - update tests:

    - update test_surface_evaluate()
    - update test_laplace_unit_*.py for d_surface_flux

  - update Region.select_cells_of_surface() for passing to C
  - averaging mode in vg_integrate(), sg_integrate():

    - use vg.integrate() in de_grad, de_div terms
    - remove de_integrate()

  - fix regions, groups for argument traces
  - remove obsolete term caches:

    - remove sfepy/terms/cachesBasic.py
    - remove FiniteStrainTLDataCache, FiniteStrainSurfaceTLDataCache
    - new eval_real_extra(), dq_grad_extra()

  - update sfepy/homogenization/recovery.py
  - update FESurface for mirror connectivities:

    - start FESurface.setup_mirror_connectivity()
    - update FESurface.get_connectivity() - add is_trace argument

  - update field and DOF connectivities for boundary traces:

    - is_trace is part of field DOF connectivity keys, active DOF
      connectivity keys and FieldVariable evaluate cache keys
    - update setup_dof_conns(), Field.setup_extra_data()
    - update Field.setup_surface_data(), .setup_dof_conns():

      - add is_trace argument

    - update Approximation.get_connectivity() - add is_trace argument
    - update Equations.get_graph_conns()
    - update FieldVariable.get_dof_conn(), .evaluate()
    - update Term.assemble_to()
    - update SurfaceField.setup_dof_conns() for is_trace argument

- docs:

  - add Primer tutorial
  - installation: new Python(x,y) instructions with umfpackpy
  - dev guide: describe directory structure
  - add local table of contents at several places
  - fix term signature table generation in process_terms()
  - fix typeset_term_table() for multi-equation definitions
  - dev guide: rewrite section on implementing new terms
  - new terms_new.rst
  - make tutorial images smaller
  - update and link release tasks to developer guide

- regions:

  - update define_box_regions() - add can_cells argument
  - add true_cells attribute to Region:

    - update Region.__init__(), .update_groups(), .set_faces(), .select_cells(),
      .select_cells_of_surface()
    - update Region.set_from_group(), .set_faces(), .complete_description()

  - update Region.get_n_cells() to return total count optionally
  - fix Region.set_faces()

- fields:

  - check region passed to field constructor:

    - new Field.check_region(), SurfaceField.check_region()

  - new SurfaceField.average_qp_to_vertices()
  - fix Approximation.describe_geometry() for surface fields
  - fix default fill value in Field.extend_dofs() for scalars
  - fix Field.get_dofs_in_region_group() for surface regions (no true cells)
  - remove unused Field.update_geometry()
  - replace FieldVariable.describe_geometry() by new .get_mapping()

- problem description:

  - update ProblemConf constructors - new override argument:

    - override given configuration items using override dict
    - new ProblemConf.dict_from_string()
    - new sfepy/base/parse_conf.py: create_bnf()

  - schroedinger.py: new --conf, --options options - allow override
  - simple.py: new --conf, --options options - allow override
  - new ProblemConf.get_function():

    - update assign_standard_hooks()
    - update schroedinger.py

  - remove obsolete 'fe' keyword

- input-output:

  - fix reading mat_id in VTKMeshIO, read node_groups (if it exists)
  - set TetgenMeshIO.getnodes(), .getele() verbose default to False
    so that it does not interfere with runTests.py
  - update UserMeshIO.read() for functions returning new mesh
  - fix blockgen.py, cylindergen.py scripts
  - update gen_block_mesh(), gen_cylinder_mesh() for non-array arguments
  - fix SimpleApp.setup_output_info() for MeshIO instances as filename_mesh

- implement geometrical surface groups in Domain:

  - new Domain.create_surface_group(), .clear_surface_groups()
  - update Approximation.describe_geometry() for surface groups
  - update Field.setup_extra_data()

- integrals:

  - update Integrals.get() to accept int as quadrature name
  - new Integral.get_key()
  - remove unused dim attribute from Integral

- terms:

  - fix Mooney-Rivlin stress term
  - remove VolumeDataCache, SurfaceDataCache:

    - simplify d_volume, d_surface terms

  - remove Term.describe_geometry(), use .get_mapping() instead

- solvers:

  - update processing of solver options
  - use new make_get_conf()
  - options in 'conf' can be overridden using kwargs
  - remove LinearSolver.set_tolerance()
  - remove SymeigEigenvalueSolver, clean and fix sfepy/solvers/eigen.py
  - remove symeig reference
  - merge common code in eigensolvers to a decorator
  - fix PETScKrylovSolver for nonlinear problems
  - fix output indentation in PysparseEigenvalueSolver.__call__()
  - always compute initial residual in time_step_function()

- schroedinger:

  - update schroedinger.py to follow style guide
  - schroedinger.py: new --save-restart, --load-restart options:

    - basic DFT restarting functionality

  - schroedinger.py: add save_dft_iterations application option
  - allow adding/subtracting zero to PotentialBase
  - make Potential iterable, numerical charge computation:

    - new Potential.__iter__(), .get_distance(), .get_charge()

  - update getvxc() to compute also energy
  - fix schroedinger.py for meshes with several element groups
  - schroedinger.py: move conf code to new ProblemConf.from_file_and_options()

- scripts:

  - extractor.py: improve help message
  - remove sfepy_gui.py
  - update test_install.py to report success/failure in log of times
  - update test_install.py to test --config option

- base:

  - new print_array_info(), update_dict_recursively(), edit_dict_strings()

- misc:

  - new SimpleApp.save_dict(), .load_dict()
  - make PolySpace.node_coors C-contiguous
  - fix ElasticConstants for SymPy 0.7.0
  - remove unused Domain.get_orientation()
  - update evaluate_at() to be more robust for tensor product elements
  - update EquationMap.map_equations() for array values
  - fix mayavi imports for version 4.0.0
  - fix output shape in dot_sequences()
  - fix State.get_scaled_norm() for zero scaling norm
  - fix set_mesh_coors(), ProblemDefinition.set_mesh_coors()
  - remove geometries and geometries0 attributes

- implement updated Lagrangian formulation for finite strain elasticity

  - example
  - setup reference state for nonlinear homogenization
  - fix and add terms, support for nonlinear homogenization
  - updating reference geometry
  - new "ulf" options, setup iter_hook
  - setup actual and initial coordinates
  - fix computation of the deformation gradient
  - new example hyperelastic_ul_up.py - ULF, displacement-pressure formulation
  - new dw_ul_volume and dw_ul_compressible terms

- new terms:

  - dw_div_w (weighted divergence term of a test function)

- examples:

  - new examples/standalone/thermal_electric/thermal_electric.py:

    - update test_install.py to test it

  - new examples/linear_elasticity/its2D_1.py
  - new examples/linear_elasticity/its2D_2.py
  - new examples/linear_elasticity/its2D_3.py

    - update test_install.py to test it

  - new examples/linear_elasticity/its2D_4.py

- setup:

  - fix setup.py files to correctly define DEBUG_FMF flag
  - allow calling "python setup.py" via Makefile
  - update site_cfg_template.py
  - fix matplotlib version for ubuntu 10.04
  - update package_check() for alternative names
  - update mayavi version check

.. _2011.1-2011.2:

from 2011.1 to 2011.2
=====================

- experimental implementation of terms - merge new_terms branch:

  - aim: easier usage and definition of terms
  - new NewTerm class, to coexist with Term for some time
  - update FieldVariable:

    - new .set_current_group(), .clear_current_group()
    - new .val(), .val_qp(), .grad(), .grad_qp(), .iter_dofs()
    - new .assign_geometries(), .clear_bases(), .setup_bases()
    - new .get_element_zeros(), .get_component_indices()
    - add spatial dimension attribute
    - cache base functions and geometries

  - new get_range_indices()
  - new examples:

    - examples/miscellaneous/compare_scalar_terms.py
    - examples/miscellaneous/compare_vector_terms.py

  - new new terms:

    - dw_new_diffusion
    - dw_new_mass_scalar
    - dw_new_mass
    - dw_new_lin_elastic

- implement basic membrane elements - merge tl_membrane_term branch:

  - new dw_tl_membrane (Mooney-Rivlin membrane term)

- update build system to use exclusively setup.py:

  - update setup.py to check dependencies, prepare for using Cython
  - update setup.py clean to remove files generated by in-place build
  - new setup.py options: htmldocs, pdfdocs, doxygendocs
  - remove all Makefiles
  - new MANIFEST.in for excluding paths/files from distribution
  - site configuration:

    - update options in site_cfg_template.py, update Config
    - remove archlib, numpy_include
    - rename opt_flags -> compile_flags

  - set options for building extension modules via Config
  - update and clean up setup.py files

- docs:

  - update for no Makefiles, add sympy as dependency
  - add sphinx build file for windows
  - fix doc/doxygen.config to exclude generated files
  - update tutorial:

    - add long, short syntax sections, add TOC
    - add basic notions, sneak peek sections

- boundary condition:

  - allow switching boundary conditions on/off depending on time:

    - update boundary conditions to have times attribute
    - update Conditions, EssentialBC, PeriodicBC, LinearCombinationBC
    - use new is_active_bc(), update matrix graph as needed

  - pass problem definition to user EBC functions:

    - update tests and examples

- postprocessing and visualization:

  - fix Viewer and add_subdomains_surface() for animations:

    - new FileSource.setup_mat_id(), set mat_id in .create_source()

  - update ViewerGUI - 'make snapshots' button
  - fix animation view setting in Viewer, ViewerGUI:

    - use all view components including distance and focal point

  - postproc.py:

    - new --opacity option
    - new --domain-specific option

  - new parse_domain_specific()
  - new sfepy/postprocess/domain_specific.py
  - new DomainSpecificPlot class, plot_displacements()

- linalg:

  - prevent nans in cg_eigs()
  - fix cg_eigs() for array input
  - new normalize_vectors()
  - update dot_sequences() to use numpy.core.umath_tests.matrix_multiply():

    - more than ten times faster, if available!

- input-output:

  - extractor.py: new --times option
  - support for variable time steps:

    - new extract_times()
    - new MeshIO.read_times(), HDF5MeshIO.read_times()
    - update VTKMeshIO, HDF5MeshIO to save true time step data

- solvers:

  - update Newton:

    - new 'give_up_warp' configuration option
    - return number of iterations in status
    - support iter_hook
    - improve bad linear solver convergence report
    - allow linear solver precision relative to residual norm:

      - new 'lin_precision' solver option

  - new VariableTimeStepper - time stepper with a variable time step
  - update ProblemDefinition for nonlinear solver iter_hook
  - improve convergence reporting of ScipyIterative, PETScKrylovSolver
  - add tolerance arguments to linear solvers

- homogenization:

  - new CorrSetBCS corrector - "zero" state with applied boundary conditions

- misc:

  - new invert_remap()
  - new Region.from_faces(), .set_faces()
  - new guess_time_units()
  - add test_install.py

- new terms:

  - dw_lin_strain_fib (linear prestrain term - defined by direct. fibers)

- removed terms:

  - remove d_surf_diffusion_integrate term (same as d_hdpm_surfdvel)

- examples:

  - new examples/diffusion/poisson_short_syntax.py
  - rearrange examples/diffusion/poisson.py

- many bug fixes =:)

.. _2010.4-2011.1:

from 2010.4 to 2011.1
=====================

- implement discontinuous fields - merge 'discontinuous' branch

  - use mesh connectivity to construct reference maps

    - independently from approximations of variables (fields)
    - always P1 or Q1, based on element geometry
    - update Approximation.describe_geometry() and related functions

  - new DiscontinuousField, DiscontinuousApproximation classes
  - use DiscontinuousField for P0, Q0 approximations
  - new eval_nodal_coors()
  - update Approximation.eval_extra_coor()
  - new Field.average_qp_to_vertices(), .interp_to_qp()

- update surface fields:

  - new SurfaceApproximation class, used by SurfaceField
  - new SurfaceInterpolant class, used by SurfaceField

- fields:

  - new Field.get_true_order(), .get_vertices()
  - new Field.evaluate_at()

    - allows different behaviour for Field subclasses

  - new Field.get_output_approx_order() - correct output order
  - remove Approximations class, move its functionality to Field

    - simplification to remove a layer of code that was not needed
    - Field.aps is an ordinary dict
    - methods using is_surface split between Field and SurfaceField

- state, variables, degrees of freedom:

  - new DofInfo.get_n_dof_total()
  - new Linear Combination BC operator: IntegralMeanValueOperator
  - new EquationMap.get_operator()
  - new State.from_variables(), .set_parts()
  - add force argument to State.set_full()
  - new Variables.check_vector_size(), use it to check DOF vectors
  - fix Variables.state_to_output() for saving by parts
  - fix Variable.advance() - prevent history modification
  - new FieldVariable.apply_ebc(), .apply_ic()

    - update Variables.apply_ebc(), .apply_ic()

  - new FieldVariable.get_full()

    - update Variables.make_full_vec(), remove var_name argument

  - new FieldVariable.get_reduced()

    - update Variables.strip_state_vector(), fix for non-first variables

  - remove Variable.get_indx(), update State accordingly

    - the indx attribute of a variable is local only, it does not index the
      state vector - dangerous to expose it

- materials: rewrite updating of material parameters

  - allow material nonlinearity, i.e. parameters dependent on state
  - ProblemDefinition instance needs to be passed into
    Materials.time_update() and related functions
  - material user function syntax changed

    - takes ts, coors, mode, equations, term, problem and group_indx
    - the coors argument are the QP coordinates for all element groups
      merged

- equations and evaluation:

  - split Equations.time_update() - new Equations.time_update_materials()
  - fix term evaluation with complex variables
  - update Equations.eval_tangent_matrices() - names argument
  - fix Equations.eval_tangent_matrices() for multi-variate terms

    - clear the matrix there, not in BasicEvaluator.eval_tangent_matrix()

  - update Equations.eval_residuals() - by_blocks, names arguments
  - new Equations.print_terms()
  - add regions argument to create_evaluable()

- terms:

  - new register_term()
  - ensure that each Term instance has caches attribute
  - ensure that all terms in equations share the same DataCaches instance

    - new Equations.setup_caches()

  - update Term.get_arg_name() - add docstring, join argument
  - fix Term.assign_caches(), .get_cache() for material arguments
  - update cachesBasic for complex values

- mesh, domain, regions:

  - generate meshes using: 2D - triangle, 3D - tetgen
  - speed-up mesh reading by using numpy.fromfile()

    - update read_array()

  - update VTKMeshIO, ComsolMeshIO, VTKMeshIO, MeditMeshIO
  - update skip_read_line(), read_token()
  - coordinate transformation matrix can be defined in options
  - allows translation in coordinate transformation
  - new Domain.refine() - uniform tetrahedral mesh refinement

    - new sfepy/fem/refine.py: refine_3_4()

  - new Region.from_vertices()
  - new Region.update_shape(), update Region.complete_description()
  - new Facets.get_complete_facets()

- problem definition:

  - new ProblemDefinition.setup_hooks()
  - fix ProblemDefinition.solve() for LCBCs
  - new ProblemDefinition.set_output_dir(), .set_solvers_instances()
  - update ProblemDefinition.init_solvers() to report presolve time
  - new ProblemDefinition.set_equations_instance()

- solvers:

  - fix time_step_function() to allow presolve if possible
  - fix Newton, SemismoothNewton exception handling (define ok)
  - update/fix TimeStepper construction arguments
  - new pre_process_hook - called in solve_direct
  - fix TimeStepper, get_print_info() for n_step set to 1
  - fix time stepping solver for quasistatic linear problems

    - new prepare_matrix()

- sfepy.linalg:

  - new get_coors_in_ball()
  - new assemble1d(), unique_rows(), infinity_norm()
  - new insert_strided_axis()
  - update cg_eigs(), sym_tri_eigen(), allow eigenvalue selection
  - fix and update dot_sequences() - general mode argument

- large deformation:

  - fix HyperElasticBase for several element groups
  - fix HyperElasticBase for matrix-only assembling
  - compute stress in matrix mode if no previous residual mode call
  - fix BulkPressureTLTerm (dw_tl_bulk_pressure) for several element groups
  - update VolumeTLTerm to work in initial step
  - fix error handling in hyperelastic term caches
  - new sfepy/mechanics/membranes.py: functions for membranes

    - new describe_deformation(), get_tangent_stress_matrix(),
      create_transformation_matrix(), create_mapping(), get_invariants()

- schroedinger: update and clean up

  - update for recent changes (on the fly geometries, Materials under
    Equations, for material nonlinearity, ...)
  - add iter_hook_final application option
  - fix getting parametric hook for inputs with define()

- homogenization:

  - update plotPerfusionCoefs.py
  - new CorrEqPar - parametrized equation via 'eq_pars'
  - update TCorrectorsViaPressureEVP for Variable.get_full()
  - allow complex coefficients
  - fix CorrectorsPermeability
  - update PressureEigenvalueProblem
  - remove obsolete fading memory coefficients

    - remove ViscousFMCoef, BiotFMCoef, BiotFM2Coef, FMRBiotModulus

  - update fading memory coefficients CoefFMSymSym, CoefFMSym, CoefFMOne
  - fading corrector file names obtained by set_variables()
  - update time dependent pressure eigenvalue problem based correctors
  - update TCorrectorsViaPressureEVP
  - update TCorrectorsRSViaPressureEVP, TCorrectorsPressureViaPressureEVP
  - raise original exception in MiniAppBase.init_solvers()
  - update recover_bones()
  - fix variable names in CorrMiniApp.get_output() for no components

- genPerMesh.py:

  - new --repeat option
  - split and move functionality into sfepy.fem
  - move fix_double_nodes(), get_min_edge_size(), get_min_vertex_distance(),
    get_min_vertex_distance_naive() into sfepy/fem/mesh.py
  - new compose_periodic_mesh()
  - remove broken --test option

- new terms:

  - d_sum_vals (sum nodal values, for use in postprocessing)
  - d_diffusion_integrate (diffusion integral term)
  - d_surf_diffusion_integrate (diffusion surface integral term)
  - dw_diffusion_coupling (diffusion copupling term)
  - new d_div term (evaluate divergence term)
  - d_surface_integrate_w (integrate a variable over a surface)
  - d_surface_dot_w (surface dot product for both scalar and vector fields)

- clean up in acoustic terms

  - some of them replaced by more general diffusion terms

- simplify acoustic/diffusion sensitivity terms

  - d_llaplace_p_sa1, d_llaplace_p_sa2,  d_llaplace_t_sa2 -> d_diffusion_sa
  - dw_surface_llaplace -> dw_surface_laplace

- remove obsolete code, clean up:

  - BasicEvaluator.strip_state_vector(), LCBCEvaluator.strip_state_vector()
  - remove obsolete function and code (_fix_scalar_dc())
  - remove Field.get_extra_nodes_as_simplices(), .write_mesh()
  - simple.py: remove --save-region-field-meshes option
  - remove code depending on removed Field.get_extra_nodes_as_simplices()

    - Mesh.from_region_and_field()
    - ProblemDefinition.save_region_field_meshes()

  - remove Field.interp_c_vals_to_n_vals()
  - remove parameter 'shape' from term di_volume_integrate_mat
  - remove read_tuple()

- docs:

  - new projections.rst, fields.rst, variables.rst

- misc:

  - remove star imports
  - fix Output.__init__() arguments
  - new Container.extend()
  - allow construction of OneTypeList from sequence
  - fix acoustic band gaps code
  - new Interpolant.get_geom_poly_space()
  - new make_l2_projection() for scalar field variables
  - add tetrahedron quadratures of order 4 and 6
  - update get_physical_qps() - use slices for efficiency
  - update Viewer - plot scalar cell data as cell data if possible
  - isfepy: update startup message

- tests and examples:

  - tests/test_projections.py: new test_projection_tri_quad()
  - new examples/linear_elasticity/material_nonlinearity.py + test
  - fix, update examples/diffusion/poisson_parametric_study.py
  - update tests/test_tensors.py

    - new test_transform_data(), test_stress_transform()

  - new tests/test_linalg.py - test dot_sequences(), insert_strided_axis()
  - update tests/test_linalg.py - new test_unique_rows(), test_assemble1d()
  - new tests/test_domain.py - very basic tests of facets and refinement

- many bug fixes

.. _2010.3-2010.4:

from 2010.3 to 2010.4
=====================

- base:

  - better printing formatting for basic data types

- docs:

  - use viewcode Sphinx extension
  - add gitwash tutorial (adapted from Numpy)

- sfepy.linalg:

  - new insert_sparse_to_csr() - insert a sparse matrix into a CSR matrix
  - new compose_sparse()
  - improve structuring:

    - move some functions from sfepy.linalg.utils to sfepy.linalg.geometry
    - remove unneeded functions

  - simplify importing:

    - import all its contents into sfepy.linalg namespace

  - new sfepy/linalg/eigen.py - eigenvalues utility functions
  - new sfepy/linalg/geometry.py - barycentic coordinates and simplex utilities

- conditions:

  - make active LCBC-constrained DOF information always defined
  - update make_global_lcbc_operator() to preserve matrix blocks

    - also create and return global active LCBC-constrained DOF information

  - new NormalDirectionOperator class

- solvers:

  - solvers: provide default name and kind for any type of conf
  - allow (re)setting data of an existing TimeStepper instance
  - use a single time stepper instance in ProblemDefinition

    - pass the instance to .set_equations() as user data to satisfy
      time-dependent term arguments

  - update Newton, SemismoothNewton - raise original residual/matrix exceptions
  - update SemismoothNewton to use compose_sparse()

    - the Jacobian needs no longer to have the non-smooth part preallocated

- refactoring of geometries (reference mappings) - merge 'geo' branch

  - create geometries as needed on the fly, similarly to term caches
  - equations only assign container for geometries to terms
  - geometries no longer stored in Approximations instances

    - greatly simplify Approximations.describe_geometry()

- new sfepy/fem/mappings.py:

  - handle reference element mappings by new Mapping, VolumeMapping,
    SurfaceMapping classes

- update Equations to create, hold and update Materials:

  - only materials actually present in equations are updated during
    ProblemDefinition.time_update() call now
  - update materials in ProblemDefinition to be created on demand
  - similar to creating variables

- DOF vector synchronization with variables - merge 'state' branch

  - new sfepy/fem/state.py
  - new State class for handling state Variables

- Domain and Mesh:

  - new Facets class for handling edges and faces
  - remove C code superseded by Facets
  - remove unused code superseded by scipy.spatial
  - new Mesh.explode_groups()

- update Field:

  - simplify Field, Approximations - assume single base and region
  - new SurfaceField - subclass of Field

    - enrich the field region syntax - allow (region, 'surface') tuple
    - add is_surface attribute to Approximations, Approximation
    - update Mesh.from_region() for surface field regions

      - useful for saving SurfaceField variables with file_per_var option

  - simplify setting Field approximation order and Interpolant construction
  - move code for getting DOFs in a region to Field
  - move DOF manipulation functions to Field

- update Equations:

  - allow passing additional connectivities to Equations.create_matrix_graph()
  - allow passing single Term to Equation.__init__()
  - update Equations.eval_tangent_matrices() - block assembling mode

- update Variables:

  - set _variables attribute in Variables.__setitem__()
    so that any Variable has it once it is added to Variables
  - new MultiplierVariable - subclass of FieldVariable

- update Terms:

  - allow different Term integration per call mode
  - simplify setting of term geometry and connectivity types:

    - new Term.integration attribute
    - new Term.setup_integration() to determine geometry and connectivity
      types according to the integration attribute
    - remove Volume, Surface, Edge, Point, SurfaceExtra constants
    - geometry types lower-cased

- expression evaluation:

  - pass integral instances instead of integral names where applicable

    - pass Term instance to .init_data() of DataCache subclasses
    - update all affected terms and term caches

  - enable calling user functions on tangent matrix in evaluators
  - check argument names consistency in ProblemDefinition.create_evaluable()

- implement higher order elements - merge 'ori' branch:

  - new NodeDescription class
  - prepare all possible facet DOF permutations
  - update Facets to store raw orientation
  - reimplement Approximations.setup_global_base():

    - for any polynomial degrees (including extra face DOFs), no C

  - update computation of extra node coordinates, no C
  - remove obsolete/unused methods and C functions
  - prepare remap vectors and DOF indices for all DOF kinds

- new sfepy/fem/projections.py

  - start projections between FE spaces

- homogenization:

  - remove unused (obsolete) correctors and coefficients
  - remove 'auxiliary' coefficients

- new sfepy/mechanics/friction.py, DualMesh class
- problem description file:

  - allow optional arguments to define()
  - update field keywords to match Field constructor arguments
  - new ANSYS CDB file reader

- output:

  - new FieldVariable.create_output(), simplify Variables.state_to_output()
  - update Variables.state_to_output() - allow skipping variables

- new terms:

  - dw_non_penetration (non-penetration condition term)
  - dw_surface_lcouple (acoustic term - derivatives in surface directions)
  - dw_surface_llaplace (acoustic term - derivatives in surface directions)
  - dq_div (new divergence in QP term)

- scripts:

  - new friction_slip.py (work in progress)
  - compare_elastic_materials.py: new --no-plot option
  - postproc.py:

    - new --subdomains option
    - update Viewer - new vector cut plane plotting mode

- tests and examples:

  - new examples/biot/biot_npbc_lagrange.py + test

    - uses dw_non_penetration term

  - update tests/test_volume.py to report volumes
  - update examples/navier_stokes/navier_stokes.py

    - check divergence-free solution

  - new tests/test_sparse.py - test compose_sparse()
  - new 'linear_elastic_up.py' example + test

    - linear elasticity, mixed formulation

  - new test_eval_matrix()
  - tests/test_meshio.py: new _compare_meshes(), test_write_read_meshes()
  - new tests/test_projections.py

- many bug fixes

.. _2010.2-2010.3:

from 2010.2 to 2010.3
=====================

- refactor for interactive use, making things simpler:

  - redesign term evaluation: non-assembling modes, hierarchy of calls:

    - hierarchy: ProblemDefinition.evaluate() - evaluate() -
      Equation.evaluate() - Term.evaluate()
    - each level can be used by itself
    - 'eval', 'el_avg', 'qp' and 'weak' modes
    - split call_mode into (call_)mode and new term_mode
    - split evaluate() into create_evaluable() and eval_equations()

  - new Domain methods to access underlying mesh
  - refactor Field, remove Fields:

    - update Field construction (remove bases)
    - move DOF connectivity setup to fields

  - refactor construction of Variables

    - move field-specific methods into FieldVariable

  - refactor Materials, Material:

    - remove regions from Material definition:

      - a Material instance is now really just a collection of values
      - region is given by a term using the particular Material

    - split material update code into several functions
    - allow mixing special constants with parameters given by user-defined
      function by passing mode='special_constant' to the function

  - refactor construction of Equations, Equation:

    - Equation.__init__() accepts Terms instance directly
    - make parse_definition() a regular function
    - update Equations to create and hold Variables
    - variables collected from individual terms
    - Equations now hold geometries instead of ProblemDefinition
    - remove term prefixes (namespaces) from description of equations
    - move setup of equations from ProblemDefinition to Equations
    - move mirror region handling to Region
    - move creation of ConnInfo into Term
    - move assembling to Equations (and Term)

  - refactor Terms, Term:

    - allow override of term arguments in Term.get_args()
    - new Term.new() factory constructor
    - simplified equation parser (full argument parsing now in
      create_arg_parser())
    - support basic arithmetics
    - set term integral at time of term construction

  - new basic boundary condition classes:
    BoundaryConditions, BoundaryCondition, EssentialBC, PeriodicBC,
    LinearCombinationBC

    - allow Function instances in conditions

  - refactor linear combination BC

    - new LCBCOperator, RigidOperator, NoPenetrationOperator,
      LCBCOperators, make_global_lcbc_operator()

  - refactor DofInfo into proper class (and module)
  - refactor equation mapping into EquationMap class
  - implement simplified integral specification in equations

    - the integral can either be a string representation of a non-negative
      integer (the integral order) or 'a' (automatic order) or a string
      beginning with 'i' (existing custom integral name)
    - integrals are created on demand

  - ConnInfo now stores directly variables instead of their names
  - update ProblemDefinition for interactive use:

    - evaluators do not hold the tangent matrix
    - split and update ProblemDefinition.time_update()

  - remove unnecessary arguments from evaluators and generic solvers
  - remove historical cruft, obsolete code
  - update all examples
  - update applications for new term evaluation:

    - schroedinger.py
    - shaper.py

- simplify interactive construction of solvers:

  - when a non-abstract class is used, name and kind are inferred
    automatically

- improve tests of examples:

  - update TestInput to call hook functions and to use solve_direct()
  - simplify TestInputEvolutionary
  - check nonlinear solver stopping conditions also for evolutionary
    problems

- homogenization:

  - new CoefSum and CorrSum
  - new CoefEval - evaluate expression (e.g. 'c.A/2 + c.B*c.C')
  - update for new evaluation code
  - simplify saving/dumping of correctors by new CorrSolution class
  - correctors stored by variables, not as the whole state vector
  - user should provide set_variables() functions for all required
    correctors/coefficients
  - pass only the direct dependencies to coefficient and corrector mini_apps

- mesh readers:

  - add support for 2d abaqus quad and tri elements
  - add full read and write support for comsol mesh format for sfepy
    supported types

- examples:

  - update examples/quantum:

    - unify 2D and 3D versions
    - remove broken DFT examples

  - new example + test (linear_elastic_tractions.py):

    - employs simplified integral definition

  - new examples/standalone/interactive/linear_elasticity.py

- tests:

  - new tests/test_high_level.py

- documentation:

  - improve docstrings:

    - add argument description for all terms
    - prepend term call signature(s) into term docstrings

  - new tutorial "Interactive Example: Linear Elasticity"

- many bug fixes
- base:

  - update Container class to be more dict-like

- new AUTHORS file

.. _2010.1-2010.2:

from 2010.1 to 2010.2
=====================

- new mesh readers:

  - MED (Salome, PythonOCC) format
  - Gambit NEU mesh format
  - UserMeshIO class:

    - creating, writing meshes by user-supplied functions

- mechanics:

  - ElasticConstants class - conversion formulas for elastic constants
  - StressTransform class to convert various stress tensors
  - basic tensor transformations

- updated documentation:

  - new sections in developer guide
  - updated tutorial
  - many new docstrings

- solvers:

  - semi-smooth Newton method
  - allow registering custom solvers

- examples:

  - usage of functions to define various parameter
  - usage of probes

- scripts:

  - simple.py: new --log, --quiet options
  - postproc.py: new --wireframe, --group-names options
  - extractor.py: new --same-dir, --to, --step options
  - split homogen.py:

    - HomogenizationApp moved to sfepy/homogenization/homogen_app.py

- new tests:

  - test region construction
  - test quadratures using symbolic integration
  - test semi-smooth Newton solver

- miscellaneous updates:

  - automatic order of variables
  - refactor integrals and quadratures
  - improve printing of Struct instances
  - IPython-enabled debug()
  - fixed probes in 2D
  - split Material.time_update() to allow easier setting of data
  - region selection of several nodes or elements by their ids
  - update dump_to_vtk() for stationary results (no time stepper)
  - update import_file(), load_classes() for package namespaces
  - update all terms for the new Term constructor
  - refactor dof connectivity setup, info, active dof info
  - refactor term argument checking
  - update equations and terms construction
  - update HomogenizationEngine to allow inter-coefficient dependencies
  - update term and cache table generation
  - run tests in alphabetical order
  - fix paths to meshes and other data in system-wide installation
  - new get_lattice_volume()
  - many small bug fixes

- new terms:

  - dw_stokes_w (Stokes term weighted by scalar function)
  - dq_biot_stress (Biot stress term in QP)
  - dq_cauchy_strain (Cauchy strain term in QP)
  - dq_cauchy_stress (Cauchy stress term in QP)
  - dq_def_grad (deformation gradient term)
  - dw_lin_prestress (linear prestress term)
  - dw_surface_mass_scalar_w (weighted surface scalar mass term)
  - de_biot_stress (averaged Biot stress term)
  - di_surface_moment (surface moment term)

.. _2009.4-2010.1:

from 2009.4 to 2010.1
=====================

- new sphinx-based documentation
- major branches merged:

  - 'interp' branch: interpolation between different meshes
  - 'shape' branch: shape optimization in optimal flow problems

- fast evaluation (in C) of Lagrange base functions:

  - new sfepy/fem/poly_spaces.py, tests/test_fem.py

- new GeometryElement class:

  - tensor product geometry now in [0, 1] instead of [-1, 1]
  - remove sfepy/eldesc/*

- clean-up of examples and meshes
- examples:

  - perfusion in the total Lagrangian (TL) formulation
  - active fibres in the TL formulation

- homogenization:

  - new examples:

    - linear elasticity, micro-macro coupling + test, micro-recovery

  - updated homogenization engine:

    - support for coefficients summing
    - improved saving of correctors

  - new acoustic and perfusion homogenized coefficients

- data probing:

  - automatic refinement of probe points
  - speed-up:

    - point caching, use cKDTree for speed
    - generate_probes() can reuse problem, probes, etc.

  - new PointsProbe data probe
  - update generate_probes() for multiple probe hooks

- postprocessing and visualization:

  - VTK source construction for any format supported by MeshIO classes
  - HDF5FileSource -> GenericFileSource
  - new GenericSequenceFileSource

- graphical logging:

  - support logging to a text file, vertical line plot (see live_plot.py)
  - update Log and ProcessPlotter for several Log instances
  - Log class: wait until figure save is acknowledged
  - convergence log support in Newton and Oseen solvers

- schroedinger: components of V evaluated point-wise in QPs
- miscellaneous updates:

  - new --save-regions-as-groups option in simple.py
  - move and update functions from extractor.py into time_history.py
  - Oseen solver: leave setup of stabilization parameters to user
  - allow also 'dq', 'de' call modes in InstantaneousBase._call()
  - split termsHyperElasticity.py to base, TL and UL parts
  - utilities for work with units of physical quantities:

    - new sfepy/mechanics/units.py

  - functions to compute tensor-related quantities usual in continuum mechanics:

    - new sfepy/mechanics/tensors.py

  - many bug fixes

- new terms:

  - d_surface (surface of a subdomain)
  - dw_volume_integrate_variable (volume integration a variable coefficient)
  - dw_diffusion_r (diffusion-like term)
  - TL formulation terms:

    - dw_tl_fib_a (hyperelastic active fibres)
    - dw_tl_bulk_pressure (hyperelastic bulk pressure)
    - dw_tl_volume (volume)
    - dw_tl_diffusion (diffusion with deformation-dependent permeability)
    - dw_tl_surface_traction (surface traction)

  - acoustic terms:

    - dw_acoustic (acoustic term)
    - d_acoustic_surface (acoustic surface term (in-plane directions))
    - d_acoustic_alpha (evaluation of acoustic term (in-plane directions))
    - dw_acoustic_integrate (integration of acoustic term (in-plane directions))
    - terms for sensitivity analysis:

      - d_sa_acoustic_alpha, d_sa_acoustic_alpha2, d_sa_acoustic_z,
        d_sa_acoustic_z2

.. _2009.3-2009.4:

from 2009.3 to 2009.4
=====================

- major branches merged:

  - 'ulf' branch: updated Lagrangian (UL) formulation
  - 'functions' branch:

    - unified passing extra arguments to boundary condition, material, and region
      functions
    - physical quadrature point generation
    - unified/improved handling of material parameters in terms:

      - all material parameters defined in physical quadrature points
      - all terms updated, some terms were coalesced into one

  - 'porous' branch: homogenized porous media

- input file keywords:

  - new 'functions' keyword

- simplifications & unifications:

  - results of all time steps of an evolutionary problem can be saved to a
    single HDF5 file
  - enable passing variables data to ProblemDefinition.solve()
  - runTests.py: allow multiple test files as command line arguments
  - Viewer.call_mlab() split and refactored
  - short syntax for periodic boundary conditions
  - simplified input file syntax of materials

- postprocessing and visualization:

  - using FileSource class abstracts the particular format for storing results:

    - VTK, HDF5 supported now

  - support for file sequences (evolutionary simulations)

    - time step selection for HDF5 (single) and VTK (sequence) files

  - animations (using ffmpeg)
  - minimalistic ViewerGUI
  - show scalar bars
  - various vector plotting modes
  - watch results file (HDF5) and add time steps as they are saved
  - listing data ranges works offscreen, summary for file sequence
  - sfepy_gui.py:  Mayavi2-based GUI to launch simulations

- changes aimed at interactive work:

  - Domain, Region, Field creation refactoring

- data probing - postprocessing mode:

  - read a previously probed data from the probe text file, re-plot them, and
    integrate them along the probe

- graphical logging:

  - dynamic adding of data groups (new axes) to Log and ProcessPlotter

- many bug fixes, namely:

  - fix import_file() for multiple imports
  - fix saving results with piece-wise constant (Q0) approximation

- miscellaneous updates:

  - quasistatic time stepping
  - new zero-order elements: 3_4_P0, 3_8_Q0
  - more elastic tensor construction functions:

    - elastic tensor from Young's modulus and Poisson's ratio
    - elastic tensors for use in mixed formulation

  - setting of parameter variables by a user-defined function
  - gen_block_mesh() can generate also 2D meshes
  - reversed Cuthill-McKee permutation algorithm, graph in-place permutation

- new terms:

  - dw_volume_wdot_scalar_eth (exponential decay dot product convolution term)
  - dw_biot_eth (exponential decay Biot convolution term)
  - dw_lin_elastic_eth (exponential decay elastic convolution term)
  - updated Lagrangian (UL) formulation terms:

    - dw_ul_bulk_penalty, dw_ul_he_neohook, dw_ul_he_mooney_rivlin

.. _2009.2-2009.3:

from 2009.2 to 2009.3
=====================

- basic support for Windows installation via numpy distutils (finally!):

  - installation using standard "python setup.py install"...

- postproc.py:

  - quite usable now for fast first glance at the results
  - plots point, cell data of all kinds (scalar, vector, tensor)
  - draw iso-surface in 3D mode
  - fixed filename in Viewer for femhub notebook
  - new options: --scalar-mode, --list-names, --only-names,
    --rel-text-width, --no-show, --roll, --view, --all, --layout, -o

- cylindergen.py:

  - cylindrical mesh generator

- probe.py:

  - can probe selected quantities only
  - new options: --only-names, --auto-dir, --same-dir

- isfepy:

  - new options: --no-wx, --no-viewer

- phono: basic support for liquid inclusions

  - support for inner band gaps detection (brute force) and plotting

- homogenization: added new-style piezo-elastic corrector and coefficient classes
- schroedinger: fixed charge density computation
- solvers:

  - added SciPy direct sparse solvers (ls.scipy_direct) - unified
    umfpack, superlu

- new terms:

  - de_grad (element average of gradient)
  - d_volume_surface (compute volume using surface integral)
  - dw_bc_newton (Newton boundary condition)
  - dq_state_in_volume_qp, dq_state_in_surface_qp (interpolating state
    into quadrature points)
  - dw_surface_integrate_variable (weak surface term with variable coefficient)

.. _2009.1-2009.2:

from 2009.1 to 2009.2:
======================

- scripts:

  - added probe.py - a script to probe and plot results saved in result files

    - data probes along geometrical objects (e.g. lines, rays) intersecting the
      mesh

  - added postproc.py - a script to visualize results saved in result files

    - added Viewer class - 3D plots using mayavi2
    - rudimentary automatic mode only

  - added isfepy (interactive sfepy) IPython shell

    - uses new pde_solve(), pre-imports mayavi2 based Viewer

- short input syntax for LCBC conditions, fields, integrals, materials and
  solvers
- automatic html documentation generation via doxygen
- new mesh readers:

  - Nastran (.bdf) format
  - Abaqus ascii (.inp)

- new example problems:

  - subdomains.py + test - test dw_jump interface term
  - stabilized_navier_stokes.py input + test - test Oseen solver
  - acoustics.py + test - compute complex acoustic pressure

- solvers:

  - changed API of nonlinear solvers so that scipy solvers can be used
  - added Broyden and Anderson nonlinear solvers (SciPy implementation)
  - updated Oseen solver

- major rewrite of handling of dof connectivities, matrix graph and term
  geometries:

  - lots of dof connectivity related code was simplified/removed
  - extra connectivity data (surface, point) computed on demand in
    Variables.setup_dof_conns()
  - support for terms with traces of variables on interface regions
  - surface data computation for terms of volume dof_conn and Surface geometry

- extended syntax of equations to allow boundary traces of variables:

  - to use when a field value at an interface boundary is needed from the
    neighbouring subdomain side and the field is not defined there
  - example: dw_jump.isurf.Gamma12_1( jump1.val, q1, p1, tr(p2) )

- refactored linear combination boundary conditions (LCBC) code:

  - fixed rigid LCBC for multi-field problems
  - added no penetration LCBC
  - major speed-up (several orders) of LCBC operator construction

    - assembled via the coo_matrix instead of the lil_matrix

  - fixed check_tangent_matrix() for LCBC

- applications:

  - homogenization:

    - prefactorize the matrix for linear corrector problems - major speed-up

  - phononic materials:

    - plot also middle eigenvalues in 3D, fixed plot labels, polarization angles
    - caching of eigenvalue problem solution and Christoffel acoustic tensor

  - schroedinger.py:

    - choose and call DFT solver via solver interface

- general:

  - fixed boundary quadrature points for multi-field problems
  - fixed complex assembling
  - fixed live plotting (ProcessPlotter) for multi-core machines
  - improved Output class, simplified its usage
  - Struct.__str__() prints in alphabetical order
  - unified version information by introducing sfepy.__version__
  - polished MeshIO class
  - implemented region selection by node groups
  - refactored Mesh nodes, lots of simplifications
  - many small fixes and updates

- new terms:

  - dw_jump (scalar interface jump term)
  - dw_surface_mass_scalar (scalar mass on a surface boundary)

.. _2008.4-2009.1:

from 2008.4 to 2009.1:
======================

- new solvers:

  - simple backtracking steepest descent optimization solver
  - PETSc Krylov solvers via petsc4py, sequential mode
  - LOBPCG eigenvalue solver (SciPy implementation)

- new mesh readers:

  - mesh3d (hermes3d)
  - AVS UCD ascii mesh
  - Hypermesh ascii mesh

- homogenization:

  - MiniAppBase base class for "mini-applications": micro-problem correctors,
    homogenized coefficients
  - unified approach to resolve data dependencies: HomogenizationEngine class

- applications:

  - phononic materials:

    - dispersion analysis, phase velocity computation for phononic materials
    - homogenized coefficients computed via the HomogenizationEngine
    - caching of coefficients to speed up parametric runs

  - schroedinger.py:

    - all functionality moved into SchroedingerApp class
    - inherits from SimpleApp -> can be parametrized
    - fixed DFT iterations, iteration plot saving
    - basic smearing around Fermi limit

- scripts:

  - convert_mesh.py:

    - --scale option, support different scaling for each axis

- general:

  - terms, caches now imported dynamically by load_classes()

    - to add a new term/cache module just put it into sfepy/terms/

  - better setup of Application options

    - automatic option update in parametric studies

  - default configuration options for the linear, nonlinear and eigen- solvers
  - various 64bit fixes
  - allow empty output prefix, combined output to file and terminal

- new terms:

  - dw_electric_source (electric source term)

.. _00.50.00-2008.4:

from 00.50.00 to 2008.4:
========================

- framework for running parametric studies
- allow time derivatives of variables as term arguments

  - example (transient diffusion):
    """dw_mass_scalar.i1.Omega( s, dT/dt )
    + dw_laplace.i1.Omega( coef.val, s, T ) = 0"""

- initial conditions via ics, ic_<number> keywords
- enhanced acoustic band gaps code

  - dispersion analysis (polarization angle calculation)
  - applied load tensor computation
  - phase velocity computation for periodic perforated media with empty holes

- term base classes

  - actual term code reduced significantly
  - adding new terms is even easier now

- type of term arguments determined fully at run-time

  - many terms were unified
  - the same term can be used both for the finite element assembling and the
    evaluation of the weak form for known fields (dw_ = d_)

- live plotting using multiprocessing module

  - assumes GTKAgg matplotlib backend
  - support for setting x axis values and labels and y labels
  - figure saving

- printing messages: Output class
- homogenized coefficients classes prototypes
- improved schroedinger.py

  - plotting DFT iterations

- created sfepy/mechanics

  - conversions of elastic constants and transformations to plane

- float format used for saving results can be set by the 'float_format' option
- new terms:

  - dw_piezo_coupling (piezo-electric coupling term)
  - dw_biot (Biot coupling term, former dw_biot_div, dw_biot_grad, ...)
  - dw_stokes (Stokes coupling term, former dw_div, dw_grad, ...)
  - dw_lin_elastic_th (linear elasticity fading memory, former dw_lin_viscous_th)
  - dw_biot_th (Biot fading memory terms unified)

.. _00.46.02-00.50.00:

from 00.46.02 to 00.50.00:
==========================

- finite strain elasticity: neo-Hookean, Mooney-Rivlin materials

  - total Lagrangian (TL) formulation
  - geometric data via finite_strain_tl DataCache

- solving problems in complex numbers
- generalized equations to allow linear combination of terms

  - example: """2 * aterm.i1.Omega( v, u ) = - 3.0 * bterm.i1.Omega2( v, u )"""

- run-time type of state term arguments

  - removed all *_r terms, now useless

- 'elements by function( domain )' region selector
- refactoring to follow Python coding style guidelines
- history support in variables
- MeshIO.read_dimension() to quickly get dimension in an input file
- improved site_cfg_template.py
- improved schroedinger.py
- new terms:

  - de_average_variable (average a variable in elements)
  - dw_surface_integrate (integrate over surface operator)
  - dw_tl_bulk_penalty (bulk penalty in TL formulation)
  - dw_tl_he_neohook (neo-Hooekan term in TL formulation)
  - dw_tl_he_mooney_rivlin (Mooney-Rivlin term in TL formulation)

.. _00.41.03-00.46.02:

from 00.41.03 to 00.46.02:
==========================

- alternative short syntax for specifying essential boundary conditions,
  variables and  regions
- saving results per variable (useful when variables defined in different
  subdomains)
- manufactured solutions tests:

  - SymPy support

- new eigenvalue solvers:

  - removed symeig dependence

- linear solvers based on PyAMG
- simple block mesh generator
- unified HDF5 mesh/solution reading/writing
- site configuration now via script/config.py + site_cfg.py
- example: computing homogenized elastic coefficients
- new terms and lots of reorganization:

  - Biot terms
  - some fading memory terms

.. _00.35.01-00.41.03:

from 00.35.01 to 00.41.03:
==========================

- works on 64 bits
- support for various mesh formats:

  - medit: .mesh
  - text VTK: .vtk
  - tetgen: .node + .ele
  - comsol: .txt

- Schroedinger equation solver

  - run via 'schroedinger.py'

- input files:

  - new syntax for variables and boundary conditions
  - improved handling of degrees of freedom

- more descriptive (and less) test and simulation messages
- new handling of approximations (-> lots of thing simplified)
- material parameters can be defined in mesh vertices
- simple.py: allow user-specified postProcessHook function
- documentation generation via prettydoc
- new solvers:

  - generic time-dependent problem solver
  - pysparse, symeig, scipy-based eigenproblem solvers
  - scipy-based iterative solvers

- new terms:

  - dw_volume_integrate (volume integral operator)
  - dw_mass_scalar_r (rhs for time-dependent Poisson problem)
  - di_volume_integrate_mat (integrate material parameters)
  - dw_volume_wdot and related terms (weighted dot product)
  - dw_mass_scalar_variable (scalar mass term with variable coefficients)
  - dw_lin_elastic and related terms (anisotropic linear elasticity)
  - dw_lin_viscous (linear viscosity)
  - de_cauchy_stress (element-averaged Cauchy stress)

from 00.31.06 to 00.35.01:
==========================

- per term integration, major rewrite of sfe.fem and related:

  - term.integral.domain( arguments ) syntax
  - 'integral_*' keyword for input files
  - each term can use its own quadrature points
  - 'field' keyword syntax changed

- minor:

  - genDocs.py: PDF term documentation generator

from 00.26.01 to 00.31.06:
==========================

- acoustic band gaps determination:

  - zones of frequencies where elastic waves do not propagate
  - based on homogenization of material made of small inclusions periodically
    embedded in an elastic matrix
  - run via 'eigen.py'

- general linear combination boundary conditions - 'lcbc' keyword:

  - rigid body motion constraint imposed on regions

- new Solver classes, solver reorganization:

  - all solvers now configured in a uniform way ('solver_[0-9]+' keywords)
  - specify solvers using 'options' keyword

- new terms:

  - dot product in a volume or on a surface region
  - vector field mass term
  - scalar field "mass", fine-coarse scalar "mass" terms:

    - used for coarse mesh -> fine mesh interpolation of scalar fields

- minor:

  - added updated findSurf.py into distribution - extract surface from a mesh
  - script/kill_*
  - script/writeMesh2D.m
  - script/writeSparseMatrixHDF5.m

from 00.22.02 to 00.26.01:
==========================

- testing framework (in the spirit of unit tests):

  - particularly tests that standard input files work
  - runTests.py: output filtering

- linear spring term (kind of a relaxed Dirichlet BC on node displacements)
- volume term
- Laplace term in 2D
- chained periodic boundary conditions resolving
- new options for simple.py: --save-field-meshes, --solve-not
- periodic mesh merger (genPerMesh.py)
- minor:

  - improved region saving
  - growing term data cache
  - sparse matrix saving into HDF5
  - point dof connectivity and geometry
  - region handling improvements (canCells flag)
  - nonlinear solver status reporting
  - distribution: test and example meshes included in the release tarball
