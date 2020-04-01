.. include:: links.inc

Introduction
============

.. only:: html

   .. contents:: Table of Contents
      :local:
      :backlinks: top

*SfePy* (http://sfepy.org) is a software for solving systems of coupled partial
differential equations (PDEs) by the finite element method in 1D, 2D and 3D. It
can be viewed both as black-box PDE solver, and as a Python package which can
be used for building custom applications. The word "simple" means that complex
FEM problems can be coded very easily and rapidly.

There is also a preliminary support for the isogeometric analysis, outlined in
:ref:`isogeometric_analysis`.

The code is written almost entirely in `Python`_, with exception of the most
time demanding routines - those are written in C and wrapped by `Cython`_ or
written directly in Cython.

*SfePy* is a free software released under the `New BSD License`_. It relies on
`NumPy`_ and `SciPy`_ (an excellent collection of tools for scientific
computations in Python). It is a multi-platform software that should work on
Linux, Mac OS X and Windows.

*SfePy* was originally developed as a flexible framework to quickly implement
and test the mathematical models developed during our various research
projects. It has evolved, however, to a rather full-featured (yet small) finite
element code. Many terms have been implemented that can be used to build the
PDEs, see :ref:`term_overview`. SfePy comes also with a number of examples that
can get you started, check :ref:`sfepy-gallery-examples-index` and
:doc:`tutorial`. Some more advanced features are discussed in :doc:`primer`.
