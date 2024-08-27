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
   get you started, check :doc:`examples/gallery`, :doc:`examples` and
   :doc:`tutorial`. Some more advanced features are discussed in :doc:`primer`.
   Examples of scientific results obtained with the help of SfePy are listed in
   :ref:`example_applications`. See also :doc:`faq`.

   SfePy can be used in parallel (work in progress), see
   :ref:`solving_problems_in_parallel`.
   There is also a preliminary support for the isogeometric analysis,
   outlined in :ref:`isogeometric_analysis`.

   The easiest way to install SfePy is to use `pip`_ or `conda`_, see
   :ref:`introduction_installation`.

   **License:** :doc:`BSD <license>`

   Related Projects
   ----------------

   - Semi-automatic generation of finite element meshes from CT/MR scans stored
     in the DICOM file format: http://sfepy.org/dicom2fem
   - PUCGen (Periodic Unit Cell Generator): https://github.com/sfepy/pucgen
   - Utilities to run parametric studies in parallel, and to scoop the output
     files produced by the studies into a dataframe:
     https://github.com/rc/soops

   Citing
   ------

   If you would like to cite the SfePy package in a paper or presentation, the
   following references can be used:

   - General article:

     - Plain text:

       Cimrman, R., Lukeš, V., Rohan, E., 2019. Multiscale finite element
       calculations in Python using SfePy. Advances in Computational
       Mathematics 45, 1897-1921. https://doi.org/10.1007/s10444-019-09666-0

       (preprint: https://arxiv.org/abs/1810.00674)

     - BibTeX::

         @article{Cimrman_Lukes_Rohan_2019,
           title =        {Multiscale finite element calculations in Python using SfePy},
           author =       {Cimrman, Robert and Lukeš, Vladimír and Rohan, Eduard},
           issn =         {1572-9044},
           doi =          {10.1007/s10444-019-09666-0},
           journal =      {Advances in Computational Mathematics},
           year =         2019,
         }

   - Performance related data for version 2021.1 and a description of the
     multi-linear terms implementation are given in:

     - Plain text:

       Cimrman, R., 2021. Fast evaluation of finite element weak forms using
       python tensor contraction packages. Advances in Engineering Software
       159, 103033. https://doi.org/10.1016/j.advengsoft.2021.103033

       (preprint: https://arxiv.org/abs/2107.04121)

     - BibTeX::

         @article{Cimrman_2021,
           title =        {Fast Evaluation of Finite Element Weak Forms Using Python Tensor Contraction Packages},
           author =       {Cimrman, Robert},
           issn =         {0965-9978},
           doi =          {10.1016/j.advengsoft.2021.103033},
           journal =      {Advances in Engineering Software},
           volume =       159,
           year =         2021,
         }

   - :ref:`example_applications` have links to related scientific articles.

   - Other references:

     - R. Cimrman. SfePy - write your own FE application. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 6th European Conference
       on Python in Science (EuroSciPy 2013), pages 65–70, 2014.
       http://arxiv.org/abs/1404.6391.

     - R. Cimrman. Enhancing SfePy with isogeometric analysis. In P. de Buyl
       and N. Varoquaux, editors, Proceedings of the 7th European Conference on
       Python in Science (EuroSciPy 2014), pages 65–72, 2014.
       http://arxiv.org/abs/1412.6407.

   Support
   -------

   Work on SfePy is partially supported by the following ongoing projects:

   - project `GA24-12291S`_ (Multiscale modelling of acoustics-driven fluid
     suspensions flows in adaptive porous structures) of the Czech Science
     Foundation, since 2024
   - project `GA23-06220S`_ (Flexoelectric periodic structures for fluid
     transport and energy harvesting) of the Czech Science Foundation, since
     2023
   - project `GF22-00863K`_ (Controllable metamaterials and smart structures:
     Nonlinear problems, modelling and experiments) of the Czech Science
     Foundation (LEAD Agency), since 2022

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
