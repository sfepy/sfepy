.. include:: links.inc

.. SfePy documentation master file, created by
   sphinx-quickstart on Wed Oct 14 00:02:22 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. only:: html

   SfePy: Simple Finite Elements in Python
   =======================================

   SfePy is a software for solving systems of coupled partial differential
   equations (PDEs) by the finite element method in 1D, 2D and 3D. It can be
   viewed both as black-box PDE solver, and as a Python package which can be
   used for building custom applications. The word "simple" means that complex
   FEM problems can be coded very easily and rapidly.

   SfePy can use many terms to build the PDEs to be solved, see
   :ref:`term_overview`. SfePy comes also with a number of examples that can
   get you started, check :ref:`sfepy-gallery-examples-index` and
   :doc:`tutorial`. Some more advanced features are discussed in :doc:`primer`.

   SfePy can be used in parallel (work in progress), see
   :ref:`solving_problems_in_parallel`.
   There is also a preliminary support for the isogeometric analysis,
   outlined in :ref:`isogeometric_analysis`.

   The easiest way to install SfePy is to use `Anaconda`_, see
   :ref:`multi_platform_distributions_notes` in
   :ref:`introduction_installation`.

   **License:** :doc:`BSD <license>`

   Links
   -----

      .. list-table::

         * - Development version documentation
           - https://sfepy.org/doc-devel
         * - Latest release version documentation
           - https://sfepy.org/doc
         * - Automatically generated gallery
           - https://sfepy.org/gallery
         * - Mailing list (both user and developer)
           - `sfepy(at)python.org
             <https://mail.python.org/mm3/mailman3/lists/sfepy.python.org>`_
         * - Github organization
           - https://github.com/sfepy
         * - Source code (main git repository)
           - https://github.com/sfepy/sfepy
         * - Bug/issue tracking
           - https://github.com/sfepy/sfepy/issues

   Applications
   ------------

   Here we list some of the applications SfePy is/was developed for.

   - homogenization of porous media - parallel flows in a deformable porous
     medium
   - acoustic band gaps, homogenization of a strongly heterogenous
     elastic structure: phononic materials
   - acoustic waves in thin perforated layers
   - finite element formulation of Schroedinger equation
   - flow of a two-phase non-Newtonian fluid medium in a general domain - oil
     expression in screw presses/extruders

   Related Projects
   ----------------

   - Semi-automatic generation of finite element meshes from CT/MR scans stored
     in the DICOM file format: http://sfepy.org/dicom2fem

   Citing
   ------

   If you would like to cite the SfePy package in a paper or presentation, the
   following can be used:

   - General use:

     - Plain text:

       R. Cimrman. SfePy - write your own FE application. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 6th European Con- ference
       on Python in Science (EuroSciPy 2013), pages 65–70, 2014.
       http://arxiv.org/abs/1404.6391.

     - BibTeX::

         @InProceedings{cimrman14:_sfepy_write_your_own_fe_applic,
           author =       {Robert Cimrman},
           title =        {{SfePy} - Write Your Own {FE} Application},
           booktitle =    {Proceedings of the 6th European Conference on
                           Python in Science (EuroSciPy 2013)},
           pages =        {65--70},
           year =         2014,
           editor =       {Pierre de Buyl and Nelle Varoquaux},
           note =         {http://arxiv.org/abs/1404.6391},
         }

   - IGA-specific use:

     - Plain text:

       R. Cimrman. Enhancing SfePy with isogeometric analysis. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 7th European Conference on
       Python in Science (EuroSciPy 2014), pages 65–72, 2014.
       http://arxiv.org/abs/1412.6407.

     - BibTeX::

         @InProceedings{cimrman14:_enhan_sfepy_isogeom_analy,
           author =       {Robert Cimrman},
           title =        {Enhancing {SfePy} with Isogeometric Analysis},
           booktitle =    {Proceedings of the 7th European Conference on
                           Python in Science (EuroSciPy 2014)},
           pages =        {65--72},
           year =         2014,
           editor =       {Pierre de Buyl and Nelle Varoquaux},
           note =         {http://arxiv.org/abs/1412.6407},
         }

   Support
   -------

   Work on SfePy is partially supported by the following ongoing projects:

   - project GA17-12925S (Strength of materials and mechanical components based
     on iron: Multi-scale approach) of Czech Science Foundation, since 2017;
   - project GA17-01618S (Fluid Acoustics in Periodic Micro-Architectures) of
     Czech Science Foundation, since 2017;
   - project GA16-03823S (Homogenization and multi-scale computational modelling
     of flow and nonlinear interactions in porous smart structures) of Czech
     Science Foundation, since 2016;

   In past, work on SfePy was partially supported by the following projects:

   - project TA04010992 (Research and development of technologies for obtaining
     vegetable oils and cakes with an emphasis on quality of cakes as feed and
     the appropriateness of the used construction materials) of The
     Technology Agency of the Czech Republic, in 2014-2017.
   - project GAP101/12/2315 (Modelling of acoustic wave propagation in strongly
     heterogeneous media; multi-scale numerical and analytical approaches) of
     Czech Science Foundation, 2012-2016.
   - project NT13326 of Ministry of Health of the Czech Republic, in 2012-2015;
   - project GAP108/11/0853 (Nanostructures with transition metals: Towards
     ab-initio material design) of Czech Science Foundation, in 2011-2015;
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
   users_guide
   examples
   theory
   developer_guide

PDF version of the documentation: :download:`sfepy_manual.pdf`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
