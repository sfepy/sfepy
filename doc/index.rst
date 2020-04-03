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
   get you started, check `Gallery`_, :doc:`examples` and
   :doc:`tutorial`. Some more advanced features are discussed in :doc:`primer`.

   SfePy can be used in parallel (work in progress), see
   :ref:`solving_problems_in_parallel`.
   There is also a preliminary support for the isogeometric analysis,
   outlined in :ref:`isogeometric_analysis`.

   The easiest way to install SfePy is to use `Anaconda`_, see
   :ref:`multi_platform_distributions_notes` in
   :ref:`introduction_installation`.

   **License:** :doc:`BSD <license>`

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
   - PUCGen (Periodic Unit Cell Generator): https://github.com/sfepy/pucgen

   Citing
   ------

   If you would like to cite the SfePy package in a paper or presentation, the
   following reference (`full-text link <https://rdcu.be/bCn16>`_) can be used:

   - Plain text:

     Cimrman, R., Lukeš, V., Rohan, E., 2019. Multiscale finite element
     calculations in Python using SfePy. Adv Comput Math.
     https://doi.org/10.1007/s10444-019-09666-0

   - BibTeX::

       @article{Cimrman_Lukes_Rohan_2019,
         title =        {Multiscale finite element calculations in Python using SfePy},
         ISSN =         {1572-9044},
         url =          {https://doi.org/10.1007/s10444-019-09666-0},
         DOI =          {10.1007/s10444-019-09666-0},
         journal =      {Advances in Computational Mathematics},
         author =       {Cimrman, Robert and Lukeš, Vladimír and Rohan, Eduard},
         year =         2019,
       }

   - Other references:

     - R. Cimrman. SfePy - write your own FE application. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 6th European Con- ference
       on Python in Science (EuroSciPy 2013), pages 65–70, 2014.
       http://arxiv.org/abs/1404.6391.

     - R. Cimrman. Enhancing SfePy with isogeometric analysis. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 7th European Conference on
       Python in Science (EuroSciPy 2014), pages 65–72, 2014.
       http://arxiv.org/abs/1412.6407.

   Support
   -------

   Work on SfePy is partially supported by the following ongoing projects:

   - project GA19-04956S (Dynamic and nonlinear behaviour of smart structures;
     modelling and optimization) of the Czech Science Foundation, since 2019;
   - the European Regional Development Fund-Project "Application of Modern
     Technologies in Medicine and Industry"
     (No. CZ.02.1.01/0.0/0.0/17_048/0007280) of the Czech Ministry of
     Education, Youth and Sports, since 2018.

   See also :doc:`archived_support`.

   Manual in PDF
   --------------

   PDF version of the documentation: :download:`sfepy_manual.pdf`

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :hidden:
   :maxdepth: 3

   About <self>
   documentation
   development
   examples/gallery
