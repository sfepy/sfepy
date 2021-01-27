Implementation of Essential Boundary Conditions
===============================================

The essential boundary conditions can be applied in several ways. Here we
describe the implementation used in SfePy.

Motivation
----------

Let us solve a linear system :math:`A x = b` with :math:`n \times n` matrix
:math:`A` with :math:`n_f` values in the :math:`x` vector known. The known
values can be for example EBC values on a boundary, if :math:`A` comes from a
PDE discretization. If we put the known fixed values into a vector :math:`x_f`,
that has the same size as :math:`x`, and has zeros in positions that are not
fixed, we can easily construct a :math:`n \times n_r` matrix :math:`T` that
maps the reduced vector :math:`x_r` of size :math:`n_r = n - n_f`, where the
fixed values are removed, to the full vector :math:`x`:

.. math::

   x = T x_r + x_f \;.

With that the reduced linear system with a :math:`n_r \times n_r` can be
formed:

.. math::

   T^T A T x_r = T^T (b - A x_f)

that can be solved by a linear solver. We can see, that the (non-zero) known
values are now on the right-hand side of the linear system. When the known
values are all zero, we have simply

.. math::

   T^T A T x_r = T^T b \;,

which is convenient, as it allows simply throwing away the A and b entries
corresponding to the known values already during the finite element assembling.

Implementation
--------------

All PDEs in SfePy are solved in a uniform way as a system of non-linear
equations

.. math::

   f(u) = 0 \;,

where :math:`f` is the nonlinear function and :math:`u` the vector of unknown
DOFs. This system is solved iteratively by the Newton method

.. math::

   u^{new} = u^{old} - (\tdiff{f}{u^{old}})^{-1} f(u^{old})

until a convergence criterion is met. Each iteration involves solution of the
system of linear equations

.. math::

   K \Delta u = r \;,

where the tangent matrix :math:`K` and the residual :math:`r` are

.. math::

   K \equiv \tdiff{f}{u^{old}} \;,

   r \equiv f(u^{old}) \;.

Then

.. math::

   u^{new} = u^{old} - \Delta u \;.

If the initial (old) vector :math:`u^{old}` contains the values of EBCs at
correct positions, the increment :math:`\Delta u` is zero at those
positions. This allows us to assemble directly the reduced matrix :math:`T^T K
T`, the right-hand side :math:`T^T r`, and ignore the values of EBCs during
assembling. The EBCs are satisfied automatically by applying them to the
initial guess :math:`u^{0}`, that is given to the Newton solver.

Linear Problems
^^^^^^^^^^^^^^^

For linear problems we have

.. math::

   f(u) \equiv A u - b = 0 \;,

   \tdiff{f}{u} = A \;,

and so the Newton method converges in a single iteration:

.. math::

   u^{new} = u^{old} - A^{-1} (A u^{old} - b) = A^{-1} b \;.

Evaluation of Residual and Tangent Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The evaluation of the residual :math:`f` as well as the tangent matrix
:math:`K` within the Newton solver proceeds in the following steps:

- The EBCs are applied to the full DOF vector :math:`u`.
- The reduced vector :math:`u_r` is passed to the Newton solver.
- Newton iteration loop:

  - Evaluation of :math:`f_r` or :math:`K_r`:

    #. :math:`u` is reconstructed from :math:`u_r`;
    #. local element contributions are evaluated using :math:`u`;
    #. local element contributions are assembled into :math:`f_r` or
       :math:`K_r` - values corresponding to fixed DOF positions are thrown
       away.

  - The reduced system :math:`K_r \Delta u_r = r_r` is solved.
  - Solution is updated: :math:`u_r \leftarrow u_r - \Delta u_r`.
  - The loop is terminated if a stopping condition is satisfied, the solver
    returns the final :math:`u_r`.

- The final :math:`u` is reconstructed from :math:`u_r`.
