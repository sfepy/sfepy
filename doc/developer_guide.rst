Developer Guide
===============

This section purports to document the *SfePy* internals. It is mainly useful for those who wish to develop *SfePy* and understand the inner workings of the code.

Module Index
------------

sfepy.applications package
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/applications/module_sfepy_applications_application
   src/sfepy/applications/module_sfepy_applications_simple_app
   src/sfepy/applications/module_sfepy_applications_top_level

sfepy.base package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/base/module_sfepy_base_base
   src/sfepy/base/module_sfepy_base_conf
   src/sfepy/base/module_sfepy_base_getch
   src/sfepy/base/module_sfepy_base_ioutils
   src/sfepy/base/module_sfepy_base_la
   src/sfepy/base/module_sfepy_base_log
   src/sfepy/base/module_sfepy_base_plotutils
   src/sfepy/base/module_sfepy_base_progressbar
   src/sfepy/base/module_sfepy_base_reader
   src/sfepy/base/module_sfepy_base_setup
   src/sfepy/base/module_sfepy_base_tasks
   src/sfepy/base/module_sfepy_base_testing

sfepy.eldesc package
^^^^^^^^^^^^^^^^^^^^

This package contains the element description files.

sfepy.fem package
^^^^^^^^^^^^^^^^^

WARNING: The code in the fem package is undergoing rapid change. It is best to refer directly to the code base until the code stabilizes.

sfepy.geom package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/geom/module_sfepy_geom_femlab
   src/sfepy/geom/module_sfepy_geom_geometry
   src/sfepy/geom/module_sfepy_geom_gmsh
   src/sfepy/geom/module_sfepy_geom_meshutils
   src/sfepy/geom/module_sfepy_geom_tetgen

sfepy.homogenization package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/homogenization/module_sfepy_homogenization_coefficients
   src/sfepy/homogenization/module_sfepy_homogenization_coefs_base
   src/sfepy/homogenization/module_sfepy_homogenization_coefs_elastic
   src/sfepy/homogenization/module_sfepy_homogenization_coefs_piezo
   src/sfepy/homogenization/module_sfepy_homogenization_engine
   src/sfepy/homogenization/module_sfepy_homogenization_pfdpm
   src/sfepy/homogenization/module_sfepy_homogenization_phono
   src/sfepy/homogenization/module_sfepy_homogenization_prolong
   src/sfepy/homogenization/module_sfepy_homogenization_recovery
   src/sfepy/homogenization/module_sfepy_homogenization_utils

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mechanics/module_sfepy_mechanics_matcoefs

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/physics/module_sfepy_physics_extmods_dft

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/postprocess/module_sfepy_postprocess_dataset_manager
   src/sfepy/postprocess/module_sfepy_postprocess_sources
   src/sfepy/postprocess/module_sfepy_postprocess_time_history
   src/sfepy/postprocess/module_sfepy_postprocess_utils
   src/sfepy/postprocess/module_sfepy_postprocess_viewer

sfepy.solvers package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/solvers/module_sfepy_solvers_eigen
   src/sfepy/solvers/module_sfepy_solvers_generic
   src/sfepy/solvers/module_sfepy_solvers_ls
   src/sfepy/solvers/module_sfepy_solvers_nls
   src/sfepy/solvers/module_sfepy_solvers_optimize
   src/sfepy/solvers/module_sfepy_solvers_oseen
   src/sfepy/solvers/module_sfepy_solvers_setup
   src/sfepy/solvers/module_sfepy_solvers_solvers
   src/sfepy/solvers/module_sfepy_solvers_ts

sfepy.terms package
^^^^^^^^^^^^^^^^^^^

The terms documentation is incomplete here as yet. Please refer to the terms
documentation in sfepy_manual.pdf until it is merged with Sphinx documentation.

.. toctree::
   :maxdepth: 2

   src/sfepy/terms/module_sfepy_terms_termsAdjointNavierStokes
   src/sfepy/terms/module_sfepy_terms_termsBasic
   src/sfepy/terms/module_sfepy_terms_termsBiot
   src/sfepy/terms/module_sfepy_terms_termsElectric
   src/sfepy/terms/module_sfepy_terms_termsLaplace
   src/sfepy/terms/module_sfepy_terms_termsLinElasticity
   src/sfepy/terms/module_sfepy_terms_termsMass
   src/sfepy/terms/module_sfepy_terms_termsNavierStokes
   src/sfepy/terms/module_sfepy_terms_termsPiezo
   src/sfepy/terms/module_sfepy_terms_termsPoint
   src/sfepy/terms/module_sfepy_terms_termsSurface
   src/sfepy/terms/module_sfepy_terms_termsVolume
   src/sfepy/terms/module_sfepy_terms_terms_fibres
   src/sfepy/terms/module_sfepy_terms_terms_hyperelastic_base
   src/sfepy/terms/module_sfepy_terms_terms_hyperelastic_tl
   src/sfepy/terms/module_sfepy_terms_terms_hyperelastic_ul
