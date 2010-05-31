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
    :label: eq_lcbc_subs

    u_j(\ul{x}) = - \frac{1}{a_j} \sum_{i = 1, i \neq j}^{L} a_i u_i(\ul{x}) \;,
    \quad \forall \ul{x} \in \omega \;,

into the equations. This is done, however, after the discretization by
the finite element method.

Finite Element Approximation
----------------------------

**NOTE** Below is the correct treatment. The code uses its somewhat
simplified version, as it employs the fact that we use only Lagrange
base functions.

On a finite element :math:`T_K` (surface or volume) we have :math:`u_i(\ul{x}) =
\sum_{k=1}^{N} u_i^k \phi^k (\ul{x})`, where :math:`\phi^k` are the
local (element) base functions. The relation :eq:`eq_lcbc_subs` can be
then written as

.. math::
    :label: eq_lcbc_subs_fe

    \sum_{k=1}^{N} u_j^k \phi^k(\ul{x}) = - \frac{1}{a_j} \sum_{i=1,
    i \neq j}^{L} a_i \sum_{k=1}^{N} u_i^k \phi^k(\ul{x}) \;,
    \quad \forall \ul{x} \in T_K \;.

This relation cannot be used directly - we need to express the
individual coefficients $u_j^k$. As we have :math:`N` coefficients in
the element $T_K$, we evaluate :eq:`eq_lcbc_subs_fe` in :math:`N` points
:math:`\ul{x}^k \in T_K` and obtain the following linear system to
solve:


