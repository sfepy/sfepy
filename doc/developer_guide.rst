Developer Guide
===============

This section purports to document the SfePy internals. It is mainly useful for those who which to develop SfePy and understand the inner workings of the code.

Module Index
------------

sfepy.applications package
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   module_sfepy_applications_application
   module_sfepy_applications_simple_app
   module_sfepy_applications_top_level

sfepy.base package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_base_base
   module_sfepy_base_conf
   module_sfepy_base_getch
   module_sfepy_base_ioutils
   module_sfepy_base_la
   module_sfepy_base_log
   module_sfepy_base_plotutils
   module_sfepy_base_progressbar
   module_sfepy_base_reader
   module_sfepy_base_setup
   module_sfepy_base_tasks
   module_sfepy_base_testing

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

   module_sfepy_geom_femlab
   module_sfepy_geom_geometry
   module_sfepy_geom_gmsh
   module_sfepy_geom_meshutils
   module_sfepy_geom_tetgen

sfepy.homogenization package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_homogenization_coefficients
   module_sfepy_homogenization_coefs_base
   module_sfepy_homogenization_coefs_elastic
   module_sfepy_homogenization_coefs_piezo
   module_sfepy_homogenization_engine
   module_sfepy_homogenization_pfdpm
   module_sfepy_homogenization_phono
   module_sfepy_homogenization_prolong
   module_sfepy_homogenization_recovery
   module_sfepy_homogenization_utils

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_mechanics_matcoefs

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_physics_extmods_dft

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_postprocess_dataset_manager
   module_sfepy_postprocess_sources
   module_sfepy_postprocess_utils
   module_sfepy_postprocess_viewer

sfepy.solvers package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   module_sfepy_solvers_eigen
   module_sfepy_solvers_generic
   module_sfepy_solvers_ls
   module_sfepy_solvers_nls
   module_sfepy_solvers_optimize
   module_sfepy_solvers_oseen
   module_sfepy_solvers_setup
   module_sfepy_solvers_solvers
   module_sfepy_solvers_ts

sfepy.terms package
^^^^^^^^^^^^^^^^^^^

The terms documentation is incomplete here as yet. Please refer to the terms documentation until it is merged with Sphinx documentation.

