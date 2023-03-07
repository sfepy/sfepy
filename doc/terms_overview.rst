.. _term_overview:

Term Overview
=============

Term Syntax
-----------

In general, the syntax of a term call is:

.. centered::
   ``<term name>.<i>.<r>( <arg1>, <arg2>, ... )``,

where ``<i>`` denotes an integral name (i.e. a name of numerical quadrature
to use) and ``<r>`` marks a region (domain of the integral).

The following notation is used:

.. list-table:: Notation.
   :widths: 20 80
   :header-rows: 1

   * - symbol
     - meaning
   * - :math:`\Omega`
     - cell (volume) (sub)domain
   * - :math:`\Gamma`
     - facet (surface) (sub)domain
   * - :math:`\cal{D}`
     - cell or facet (sub)domain
   * - :math:`d`
     - dimension of space
   * - :math:`t`
     - time
   * - :math:`y`
     - any function
   * - :math:`\ul{y}`
     - any vector function
   * - :math:`\ul{n}`
     - unit outward normal
   * - :math:`q`
     - scalar test or parameter function
   * - :math:`p`
     - scalar unknown or parameter function
   * - :math:`\ul{v}`
     - vector test or parameter function
   * - :math:`\ul{w}`, :math:`\ul{u}`
     - vector unknown or parameter function
   * - :math:`\ull{e}(\ul{u})`
     - Cauchy strain tensor (:math:`\frac{1}{2}((\nabla u) + (\nabla u)^T)`)
   * - :math:`\ull{F}`
     - deformation gradient :math:`F_{ij} = \pdiff{x_i}{X_j}`
   * - :math:`J`
     - :math:`\det(F)`
   * - :math:`\ull{C}`
     -  right Cauchy-Green deformation tensor :math:`C = F^T F`
   * - :math:`\ull{E}(\ul{u})`
     - Green strain tensor :math:`E_{ij} = \frac{1}{2}(\pdiff{u_i}{X_j} +
       \pdiff{u_j}{X_i} + \pdiff{u_m}{X_i}\pdiff{u_m}{X_j})`
   * - :math:`\ull{S}`
     -  second Piola-Kirchhoff stress tensor
   * - :math:`\ul{f}`
     - vector volume forces
   * - :math:`f`
     - scalar volume force (source)
   * - :math:`\rho`
     - density
   * - :math:`\nu`
     - kinematic viscosity
   * - :math:`c`, :math:`\ul{c}`, :math:`\ull{c}`
     - any constant
   * - :math:`\delta_{ij}, \ull{I}`
     - Kronecker delta, identity matrix
   * - :math:`\tr{\ull{\bullet}}`
     - trace of a second order tensor (:math:`\sum_{i=1}^d \bullet_{ii}`)
   * - :math:`\dev{\ull{\bullet}}`
     - deviator of a second order tensor
       (:math:`\ull{\bullet} - \frac{1}{d}\tr{\ull{\bullet}}`)
   * - :math:`T_K \in \Tcal_h`
     - :math:`K`-th element of triangulation (= mesh) :math:`\Tcal_h` of
       domain :math:`\Omega`
   * - :math:`K \from \Ical_h`
     - :math:`K` is assigned values from :math:`\{0, 1, \dots, N_h-1\}
       \equiv \Ical_h` in ascending order

The suffix ":math:`_0`" denotes a quantity related to a previous time step.

Term names are (usually) prefixed according to the following conventions:

.. list-table:: Term name prefixes.
   :widths: 5 20 25 50
   :header-rows: 1

   * - prefix
     - meaning
     - evaluation modes
     - meaning
   * - dw
     - discrete weak
     - `'weak'`
     - terms having a virtual (test) argument and zero or more unknown
       arguments, used for FE assembling
   * - ev
     - evaluate
     - `'eval'`, `'el_eval'`, `'el_avg'`, `'qp'`
     - terms having all arguments known, modes `'el_avg'`, `'qp'` are not
       supported by all `ev_` terms
   * - de
     - discrete einsum
     - any (work in progress)
     - multi-linear terms defined using an enriched einsum notation

Evaluation modes 'eval', 'el_avg' and 'qp' are defined as follows:

.. list-table:: Evaluation modes.
   :widths: 20 80
   :header-rows: 1

   * - mode
     - definition
   * - 'eval'
     - :math:`\int_{\cal{D}} (\cdot)`
   * - 'el_avg'
     - vector for :math:`K \from \Ical_h: \int_{T_K} (\cdot) / \int_{T_K} 1`
   * - 'qp'
     - :math:`(\cdot)|_{qp}`


.. _term_table:

Term Table
----------

Below we list all the terms available in automatically generated tables. The
first column lists the name, the second column the argument lists and the third
column the mathematical definition of each term. The terms are devided into the
following tables:

* `Table of basic terms`_

* `Table of large deformation terms`_ (total/updated Lagrangian formulation)

* `Table of sensitivity terms`_

* `Table of special terms`_

* `Table of multi-linear terms`_

The notation ``<virtual>`` corresponds to a test function,
``<state>`` to a unknown function and ``<parameter>`` to a known function. By
``<material>`` we denote material (constitutive) parameters, or, in general, any
given function of space and time that parameterizes a term, for example
a given traction force vector.

.. include:: term_table.rst
