.. include:: links.inc

.. _sec-solving-pdes-fem:

Notes on solving PDEs by the Finite Element Method
--------------------------------------------------

The Finite Element Method (FEM) is the numerical method for solving Partial
Differential Equations (PDEs). FEM was developed in the middle of XX. century
and now it is widely used in different areas of science and engineering,
including mechanical and structural design, biomedicine, electrical and power
design, fluid dynamics and other. FEM is based on a very elegant mathematical
theory of weak solution of PDEs. In this section we will briefly discuss basic
ideas underlying FEM.

.. _poisson-weak-form-tutorial:

Strong form of Poisson's equation and its integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us start our discussion about FEM with the strong form of Poisson's
equation

.. math::
   :label: eq_spoisson

   \Delta T = f(x), \quad x \in \Omega,

.. math::
   :label: eq_spoisson_dbc

   T = u(x), \quad x \in \Gamma_D,

.. math::
   :label: eq_spoisson_nbc

   \nabla T \cdot \mathbf{n} = g(x), \quad x \in \Gamma_N,

where :math:`\Omega \subset \mathbb{R}^n` is the solution domain with the
boundary :math:`\partial \Omega`, :math:`\Gamma_D` is the part of the boundary
where Dirichlet boundary conditions are given, :math:`\Gamma_N` is the part of
the boundary where Neumann boundary conditions are given, :math:`T(x)` is the
unknown function to be found, :math:`f(x), u(x), g(x)` are known functions.

FEM is based on a weak formulation. The weak form of the
equation :eq:`eq_spoisson` is

.. math::
    \int\limits_{\Omega} (\Delta T - f) \cdot s \, \mathrm{d}\Omega = 0,

where :math:`s` is a **test** function. Integrating this equation by parts

.. math::
  0 = \int\limits_{\Omega} (\Delta T - f) \cdot s \, \mathrm{d}\Omega
    = \int\limits_{\Omega} \nabla \cdot (\nabla T) \cdot s \, \mathrm{d}\Omega
      - \int_{\Omega} f \cdot s \, \mathrm{d}\Omega =

.. math::
    = - \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      + \int\limits_{\Omega} \nabla \cdot (\nabla T \cdot s) \, \mathrm{d}\Omega
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega

and applying Gauss theorem we obtain:

.. math::
  0 = - \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      + \int\limits_{\Gamma_D \cup \Gamma_N} \!\!\!\!
            s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega

or

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
  =
    \int\limits_{\Gamma_D \cup \Gamma_N} \!\!\!\!
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
  - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega.

The surface integral term can be split into two integrals, one over the
Dirichlet part of the surface and second over the Neumann part

.. math::
    :label: eq_wpoisson_full

    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_D}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    + \int\limits_{\Gamma_N}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega.

The equation :eq:`eq_wpoisson_full` is the initial weak form of the
Poisson's problem :eq:`eq_spoisson`--:eq:`eq_spoisson_nbc`. But we can not work
with it without applying the boundary conditions. So it is time to talk about
the boundary conditions.

Dirichlet Boundary Conditions
"""""""""""""""""""""""""""""

On the Dirichlet part of the surface we have two restrictions. One is the Dirichlet
boundary conditions :math:`T(x) = u(x)` as they are, and the second is the
integral term over :math:`\Gamma_D` in equation :eq:`eq_wpoisson_full`. To be
consistent we have to use only the Dirichlet conditions and avoid the integral
term. To implement this we can take the function :math:`T \in V(\Omega)` and
the test function :math:`s \in V_0(\Omega)`, where

.. math::
    V(\Omega) = \{v(x) \in H^1(\Omega)\},

.. math::
  V_0(\Omega) = \{v(x) \in H^1(\Omega); v(x) = 0, x \in \Gamma_D\}.

In other words the unknown function :math:`T` must be continuous together with
its gradient in the domain. In contrast the test function :math:`s` must be
also continuous together with its gradient in the domain but it should be zero
on the surface :math:`\Gamma_D`.

With this requirement the integral term over Dirichlet part of the surface
is vanishing and the weak form of the Poisson equation for
:math:`T \in V(\Omega)` and :math:`s \in V_0(\Omega)` becomes

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_N}
        s \cdot (\nabla T \cdot \mathbf{n}) \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega,

    T(x) = u(x), \quad x \in \Gamma_D.

That is why Dirichlet conditions in FEM terminology are called
**Essential Boundary Conditions**. These conditions are not a part of the weak
form and they are used as they are.

Neumann Boundary Conditions
"""""""""""""""""""""""""""

The Neumann boundary conditions correspond to the known flux
:math:`g(x) = \nabla T \cdot \mathbf{n}`. The integral term over the Neumann
surface in the equation :eq:`eq_wpoisson_full` contains exactly the same flux.
So we can use the known function :math:`g(x)` in the integral term:

.. math::
    \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
    =
    \int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma
    - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega,

where test function :math:`s` also belongs to the space :math:`V_0`.

That is why Neumann conditions in FEM terminology are called
**Natural Boundary Conditions**. These conditions are a part of weak form terms.

The weak form of the Poisson's equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can write the resulting weak form for the Poisson's problem
:eq:`eq_spoisson`--:eq:`eq_spoisson_nbc`. For any test function
:math:`s \in V_0(\Omega)` find :math:`T \in V(\Omega)`
such that

.. math::
    :label: eq_wpoisson_final

    \boxed {
    \begin{split}
      \int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega
      & =
      \int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma
      - \int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega, \quad \mbox{and}\\
      T(x) & = u(x), \quad x \in \Gamma_D.
    \end{split}
    }

Discussion of discretization and meshing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is planned to have an example of the discretization based on the Poisson's
equation weak form :eq:`eq_wpoisson_final`. For now, please refer to the
wikipedia page `Finite Element Method`_ for a basic description of the
disretization and meshing.


Numerical solution of the problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To solve numerically given problem based on the weak form
:eq:`eq_wpoisson_final` we have to go through 5 steps:

#. Define geometry of the domain :math:`\Omega` and
   surfaces :math:`\Gamma_D` and :math:`\Gamma_N`.
#. Define the known functions :math:`f`, :math:`u` and :math:`g`.
#. Define the unknown function :math:`T` and the test functions :math:`s`.
#. Define essential boundary conditions (Dirichlet conditions)
   :math:`T(x) = u(x), x \in \Gamma_D.`
#. Define equation and natural boundary conditions (Neumann conditions)
   as the set of all integral terms
   :math:`\int\limits_{\Omega} \nabla T \cdot \nabla s \, \mathrm{d}\Omega`,
   :math:`\int\limits_{\Gamma_N} g \cdot s \, \mathrm{d}\Gamma`,
   :math:`\int\limits_{\Omega} f \cdot s \, \mathrm{d}\Omega`.
