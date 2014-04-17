.. SfePy documentation master file, created by
   sphinx-quickstart on Wed Oct 14 00:02:22 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. only:: html

   SfePy: Simple Finite Elements in Python
   =======================================

   SfePy is a software for solving systems of coupled partial differential
   equations (PDEs) by the finite element method in 2D and 3D. It can be viewed
   both as black-box PDE solver, and as a Python package which can be used for
   building custom applications. The word "simple" means that complex FEM
   problems can be coded very easily and rapidly.

   SfePy can use many terms to build systems of partial differential
   equations (PDEs) to be solved - follow this quick link to the table of
   all available terms to see the full list: :ref:`term_overview`. SfePy
   comes with a number of examples that can get you started, check
   :ref:`examples`. Some more advanced features are discussed in :doc:`primer`.

   **License:** :doc:`BSD <license>`

   Links
   -----

      .. list-table::

         * - Development version documentation
           - http://sfepy.org/doc-devel
         * - Latest release version documentation
           - http://sfepy.org/doc
         * - Automatically generated gallery
           - http://sfepy.org/gallery/gallery
         * - Mailing list (both user and developer)
           - http://groups.google.com/group/sfepy-devel
         * - Github organization
           - http://github.com/sfepy
         * - Source code (main git repository)
           - https://github.com/sfepy/sfepy
         * - Bug/issue tracking
           - https://github.com/sfepy/sfepy/issues

   Applications
   ------------

   Here we list some of the applications SfePy is developed for.

   - homogenization of porous media (parallel flows in a deformable porous
     medium)
   - acoustic band gaps (homogenization of a strongly heterogenous
     elastic structure: phononic materials)
   - acoustic waves in thin perforated layers
   - shape optimization in incompressible flow problems
   - finite element formulation of Schroedinger equation

   Featured Applications
   ---------------------

   - Fish heart model: http://sfepy.org/fish_heart
   - Phononic materials: http://sfepy.org/phononic

   Related Projects
   ----------------

   - Semi-automatic generation of finite element meshes from CT/MR scans stored
     in the DICOM file format: http://sfepy.org/dicom2fem

   Support
   -------

   Work on SfePy is partially supported by the following projects:

   - project GAP108/11/0853 (Nanostructures with transition metals: Towards
     ab-initio material design) of Czech Science Foundation, since 2011;
   - project NT13326 of Ministry of Health of the Czech Republic, since 2012;
   - project GAP101/12/2315 (Modelling of acoustic wave propagation in strongly
     heterogeneous media; multi-scale numerical and analytical approaches) of
     Czech Science Foundation, since 2012.

   In past, work on SfePy was partially supported by the following projects:

   - project MSM4977751303 (Failure prediction of heterogeneous materials,
     components of mechanical and biomechanical systems) of Ministry of
     Education, Youth and Sports of the Czech Republic, in 2005-2011;
   - project GA101/07/1471 (Finite element modelling of linear, non-linear and
     multiscale effects in wave propagation in solids and heterogeneous media)
     of Czech Science Foundation, in 2007-2011;
   - project GA106/09/0740 (Microstructure oriented hierarchical modeling of
     brain perfusion for CT based cerebral blood flow evaluation) of Czech
     Science Foundation, in 2009-2012.

.. _documentation:

Documentation
=============

.. toctree::
   :maxdepth: 2

   introduction
   installation
   tutorial
   primer
   users_guide
   examples
   developer_guide
   notes

PDF version of the documentation: :download:`sfepy_manual.pdf`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

