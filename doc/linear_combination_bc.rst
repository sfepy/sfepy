Linear Combination Boundary Conditions
======================================

By linear combination boundary conditions (LCBCs) we mean conditions of the
following type:

.. math::
    :label: eq_lcbc

    \sum_{i = 1}^{L} a_i u_i(\ul{x}) = 0 \;, \quad \forall \ul{x} \in
    \omega \;,

where :math:`a_i` are given coefficients, :math:`u_i(\ul{x})` are some
components of unknown fields evaluated point-wise in points
:math:`\ul{x} \in \omega`, and :math:`\omega` is a subset of the entire
domain :math:`\Omega` (e.g. a part of its boundary). Note that the
coefficients :math:`a_i` can also depend on :math:`\ul{x}`.

A typical example is the *no penetration* condition
:math:`\ul{u}(\ul{x}) \cdot \ul{n}(\ul{x}) = 0`, where :math:`\ul{u}` are the
velocity (or displacement) components, and :math:`\ul{n}` is the unit
normal outward to the domain boundary.

Enforcing LCBCs
---------------

There are several methods to enforce the conditions:

* penalty method
* substitution method

We use the substitution method, e.i. we choose :math:`j` such that
:math:`a_j \neq 0` and substitute

.. math::
    :label: eq_lcbc_subs0

    u_j(\ul{x}) = - \frac{1}{a_j} \sum_{i = 1, i \neq j}^{L} a_i u_i(\ul{x}) \;,
    \quad \forall \ul{x} \in \omega \;,

into the equations. This is done, however, after the discretization by
the finite element method, as explained below.

Let us denote :math:`c_i = \frac{a_i}{a_j}` (:math:`j` is fixed). Then

.. math::
    :label: eq_lcbc_subs

    u_j(\ul{x}) = - \sum_{i = 1, i \neq j}^{L} c_i u_i(\ul{x}) \;,
    \quad \forall \ul{x} \in \omega \;.


Weak Formulation
----------------

We multiply :eq:`eq_lcbc_subs` by a test function :math:`v_j` and
integrate the equation over :math:`\omega` to obtain

.. math::
    :label: eq_lcbc_weak

    \int_{\omega} v_j u_j(\ul{x})
    = - \int_{\omega} v_j \sum_{i = 1, i \neq j}^{L} c_i u_i(\ul{x}) \;,
    \quad \forall v_j \in H(\omega) \;,

where :math:`H(\omega)` is some suitable function space (e.g. the same space
which :math:`u_j` belongs to).

Finite Element Approximation
----------------------------

On a finite element :math:`T_K` (facet or cell) we have :math:`u_i(\ul{x}) =
\sum_{k=1}^{N} u_i^k \phi^k (\ul{x})`, where :math:`\phi^k` are the
local (element) base functions. Using the more compact matrix notation
:math:`\ub_i = [u_i^1, \dots, u_i^N]`, :math:`\vphib = [\vphib^1, \dots,
\vphib^N]^T` we have :math:`u_i(\ul{x}) = \vphib(\ul{x}) \ub_i` and similarly
:math:`v_i(\ul{x}) = \vphib(\ul{x}) \vb_i`.

The relation :eq:`eq_lcbc_subs`, restricted to :math:`T_K`, can be
then written (we omit the :math:`\ul{x}` arguments) as

.. math::
    :label: eq_lcbc_weak_fe

    \int_{T_K} \vb_j^T \vphib^T \vphib \ub_j
    = - \int_{T_K} \vb_j^T \vphib^T\sum_{i = 1, i \neq j}^{L} c_i
    \vphib\ub_i \;, \quad \forall \vb_j

As :eq:`eq_lcbc_weak_fe` holds for any :math:`\vb_j`, we have a linear
system to solve. After denoting the "mass" matrices :math:`\Mb =
\int_{T_K} \vphib^T \vphib`, :math:`\Mb_i = \int_{T_K} c_i \vphib^T
\vphib` the linear system is

.. math::
    :label: eq_lcbc_weak_fe_m

    \Mb \ub_j = - \sum_{i = 1, i \neq j}^{L} \Mb_i \ub_i \;.

Then the individual coefficients :math:`\ub_j` can be expressed as

.. math::
    :label: eq_lcbc_weak_fe_m_s

    \ub_j = - \Mb^{-1} \sum_{i = 1, i \neq j}^{L} \Mb_i \ub_i \;.

Implementation
--------------

Above is the general treatment. The code uses its somewhat simplified
version described here. If the coefficients :math:`c_i` are constant in
the element :math:`T_K`, i.e. :math:`c_i(\ul{x}) = \bar c_i` for
:math:`x \in T_K`, we can readily see that :math:`\Mb_i = \bar c_i
\Mb`. The relation :eq:`eq_lcbc_weak_fe_m_s` then reduces to

.. math::
    :label: eq_lcbc_weak_fe_impl

    \ub_j = - \Mb^{-1} \sum_{i = 1, i \neq j}^{L} \bar c_i \Mb \ub_i
    = \sum_{i = 1, i \neq j}^{L} \bar c_i \ub_i \;,

hence we can work with the individual components of the coefficient
vectors (= degrees of freedom) only, as the above relation means, that
:math:`u_j^k = \bar c_i u_i^k` for :math:`k = 1, \dots, N`.
