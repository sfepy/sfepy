# created: 20.07.2007 (-1)

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
