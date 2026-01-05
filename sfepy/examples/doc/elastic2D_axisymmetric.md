# Axisymmetric Linear Elasticity (2D Formulation)

This document presents a two-dimensional axisymmetric formulation of
the three-dimensional linear elasticity equations.
The formulation is derived under standard axisymmetry assumptions
and is equivalent to the original 3D problem.

---

## 1. Three-Dimensional Linear Elasticity

Let $ \Omega_{3D} \subset \mathbb{R}^3 $ denote an elastic body.
The static equilibrium equations read

$$
\nabla \cdot \boldsymbol{\sigma} + \boldsymbol{b} = 0
\quad \text{in } \Omega_{3D}.
$$

In the absence of body forces ($ \boldsymbol{b} = 0 $),

$$
\nabla \cdot \boldsymbol{\sigma} = 0.
$$

The material is assumed to be linear, isotropic, and homogeneous, with
the constitutive relation

$$
\boldsymbol{\sigma}
= 2\mu\, \boldsymbol{\varepsilon}(\boldsymbol{u})
+ \lambda \, \mathrm{tr}(\boldsymbol{\varepsilon}(\boldsymbol{u})) \mathbf{I},
$$

where the infinitesimal strain tensor is

$$
\boldsymbol{\varepsilon}(\boldsymbol{u})
= \tfrac{1}{2} \left( \nabla \boldsymbol{u}
+ (\nabla \boldsymbol{u})^{T} \right).
$$

---

## 2. Axisymmetric Assumptions

Cylindrical coordinates $ (r, \theta, z) $ are introduced.
The following assumptions are made:

- all fields are independent of $ \theta $,
- the circumferential displacement vanishes: $ u_\theta = 0 $.

The displacement field therefore reduces to

$$
\boldsymbol{u}(r,z) = (u_r(r,z),\, 0,\, u_z(r,z)).
$$

The computational domain is the meridian cross-section

$$
\Omega \subset \mathbb{R}^2_{(r,z)}.
$$

---

## 3. Axisymmetric Strong Form

Under the above assumptions, the equilibrium equations become

**Radial equilibrium**
$$
\partial_r \sigma_{rr}
+ \partial_z \sigma_{rz}
+ \frac{\sigma_{rr} - \sigma_{\theta\theta}}{r}
= 0
\quad \text{in } \Omega,
$$

**Axial equilibrium**
$$
\partial_r \sigma_{rz}
+ \partial_z \sigma_{zz}
+ \frac{\sigma_{rz}}{r}
= 0
\quad \text{in } \Omega.
$$

---

## 4. Axisymmetric Weak Formulation

Let $ v_r $ and $ v_z $ denote admissible test functions.
Assuming zero traction boundary conditions,
integration over the circumferential direction yields the volume element

$$
d\Omega_{3D} = 2\pi r \, d\Omega.
$$

### 4.1 Radial equilibrium (weak form)

$$
\int_\Omega
\left(
\sigma_{rr} \, \partial_r v_r
+ \sigma_{rz} \, \partial_z v_r
\right) (2\pi r)\, d\Omega
+
\int_\Omega
\sigma_{\theta\theta} \, v_r \, (2\pi)\, d\Omega
= 0.
$$

### 4.2 Axial equilibrium (weak form)

$$
\int_\Omega
\left(
\sigma_{rz} \, \partial_r v_z
+ \sigma_{zz} \, \partial_z v_z
\right) (2\pi r)\, d\Omega
= 0.
$$

---

## 5. Stress–Strain Relations in Axisymmetry

The non-zero stress components are given by

$$
\begin{aligned}
\sigma_{rr} &= 2\mu \partial_r u_r
+ \lambda \left(\partial_r u_r + \frac{u_r}{r} + \partial_z u_z \right), \\
\sigma_{\theta\theta} &= 2\mu \frac{u_r}{r}
+ \lambda \left(\partial_r u_r + \frac{u_r}{r} + \partial_z u_z \right), \\
\sigma_{zz} &= 2\mu \partial_z u_z
+ \lambda \left(\partial_r u_r + \frac{u_r}{r} + \partial_z u_z \right), \\
\sigma_{rz} &= \mu \left(\partial_z u_r + \partial_r u_z \right).
\end{aligned}
$$

---

## 6. Compact Bilinear Form

The weak formulation can be written in compact matrix form as follows.

### Radial test function

$$
\begin{aligned}
\int_\Omega (\nabla v_r)^T K_{rr}(r) (\nabla u_r)\, d\Omega
+ \int_\Omega (\nabla v_r)^T K_{rz}(r) (\nabla u_z)\, d\Omega \\
- \lambda \int_\Omega u_z \, \partial_z v_r \, d\Omega
+ \int_\Omega \frac{\lambda + 2\mu}{r} u_r v_r \, d\Omega
= 0.
\end{aligned}
$$

### Axial test function

$$
\begin{aligned}
\int_\Omega (\nabla v_z)^T K_{zr}(r) (\nabla u_r)\, d\Omega
+ \int_\Omega (\nabla v_z)^T K_{zz}(r) (\nabla u_z)\, d\Omega \\
+ \lambda \int_\Omega u_r \, \partial_z v_z \, d\Omega
= 0.
\end{aligned}
$$

### Axisymmetric material matrices

$$
K_{rr}(r) = r
\begin{pmatrix}
\lambda + 2\mu & 0 \\
0 & \mu
\end{pmatrix},
\qquad
K_{rz}(r) = r
\begin{pmatrix}
0 & \lambda \\
\mu & 0
\end{pmatrix},
$$

$$
K_{zr}(r) = r
\begin{pmatrix}
0 & \mu \\
\lambda & 0
\end{pmatrix},
\qquad
K_{zz}(r) = r
\begin{pmatrix}
\mu & 0 \\
0 & \lambda + 2\mu
\end{pmatrix}.
$$

---

## 7. Remark

This formulation is mathematically equivalent to the original
three-dimensional linear elasticity problem under the stated
axisymmetry assumptions.
It is directly suitable for two-dimensional finite element
discretizations using axisymmetric weighting.
