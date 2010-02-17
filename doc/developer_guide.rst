Developer Guide
===============

This section purports to document the *SfePy* internals. It is mainly useful for those who wish to develop *SfePy* and understand the inner workings of the code.

Module Index
------------

sfepy.applications package
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 2

   src/sfepy/applications/application
   src/sfepy/applications/simple_app
   src/sfepy/applications/top_level

sfepy.base package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/base/base
   src/sfepy/base/conf
   src/sfepy/base/getch
   src/sfepy/base/ioutils
   src/sfepy/base/la
   src/sfepy/base/log
   src/sfepy/base/plotutils
   src/sfepy/base/progressbar
   src/sfepy/base/reader
   src/sfepy/base/tasks
   src/sfepy/base/testing

sfepy.fem package
^^^^^^^^^^^^^^^^^

WARNING: The code in the fem package is undergoing rapid change. It is best to refer directly to the code base until the code stabilizes.

.. toctree::
   :maxdepth: 2

   src/sfepy/fem/geometry_element
   src/sfepy/fem/poly_spaces
   src/sfepy/fem/probes

sfepy.geom package
^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/geom/femlab
   src/sfepy/geom/geometry
   src/sfepy/geom/gmsh
   src/sfepy/geom/meshutils
   src/sfepy/geom/tetgen

sfepy.homogenization package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/homogenization/coefficients
   src/sfepy/homogenization/coefs_base
   src/sfepy/homogenization/coefs_elastic
   src/sfepy/homogenization/coefs_piezo
   src/sfepy/homogenization/convolutions
   src/sfepy/homogenization/engine
   src/sfepy/homogenization/pfdpm
   src/sfepy/homogenization/phono
   src/sfepy/homogenization/prolong
   src/sfepy/homogenization/recovery
   src/sfepy/homogenization/utils

sfepy.interactive package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/interactive

sfepy.mechanics package
^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/mechanics/matcoefs
   src/sfepy/mechanics/tensors
   src/sfepy/mechanics/units

sfepy.physics package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/physics/dft

sfepy.postprocess package
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/postprocess/dataset_manager
   src/sfepy/postprocess/sources
   src/sfepy/postprocess/time_history
   src/sfepy/postprocess/utils
   src/sfepy/postprocess/viewer

sfepy.solvers package
^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   src/sfepy/solvers/eigen
   src/sfepy/solvers/generic
   src/sfepy/solvers/ls
   src/sfepy/solvers/nls
   src/sfepy/solvers/optimize
   src/sfepy/solvers/oseen
   src/sfepy/solvers/solvers
   src/sfepy/solvers/ts

sfepy.terms package
^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   term_table
   src/sfepy/terms/termsAcoustic
   src/sfepy/terms/termsAdjointNavierStokes
   src/sfepy/terms/termsBasic
   src/sfepy/terms/termsBiot
   src/sfepy/terms/termsElectric
   src/sfepy/terms/termsLaplace
   src/sfepy/terms/termsLinElasticity
   src/sfepy/terms/termsMass
   src/sfepy/terms/termsNavierStokes
   src/sfepy/terms/termsPiezo
   src/sfepy/terms/termsPoint
   src/sfepy/terms/termsSurface
   src/sfepy/terms/termsVolume
   src/sfepy/terms/terms_fibres
   src/sfepy/terms/terms_hyperelastic_base
   src/sfepy/terms/terms_hyperelastic_tl
   src/sfepy/terms/terms_hyperelastic_ul

sfepy.terms package - full inheritance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The terms documentation is incomplete here as yet. Please refer to the terms
documentation in sfepy_manual.pdf until it is merged with Sphinx documentation.

This section is the same as the previous, but shows the full inheritance for the
terms classes.

.. toctree::
   :maxdepth: 2
   
   src/sfepy/terms/termsLinElasticity_full
