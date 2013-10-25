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

   **License:** `BSD <http://www.opensource.org/licenses/bsd-license.php>`_

   Links
   -----

      .. list-table::

         * - Development version documentation
           - http://sfepy.org/doc-devel
         * - Latest release version documentation
           - http://sfepy.org/doc
         * - Automatically generated gallery
           - http://sfepy.org/gallery/gallery
         * - Wiki pages, etc.
           - http://code.google.com/p/sfepy
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

   - project GA108/11/0853 (Nanostructures with transition metals: Towards
     ab-initio material design) of Czech Science Foundation,
   - project GA106/09/0740 (Microstructure oriented hierarchical modeling of
     brain perfusion for CT based cerebral blood flow evaluation) of Czech
     Science Foundation,
   - project NT13326 of Ministry of Health of the Czech Republic.

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

