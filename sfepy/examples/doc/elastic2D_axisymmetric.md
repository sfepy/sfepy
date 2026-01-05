# Axisymmetric Linear Elasticity (2D Formulation)

This document presents a two-dimensional axisymmetric formulation of
the three-dimensional linear elasticity equations.
The formulation is derived under standard axisymmetry assumptions
and is equivalent to the original three-dimensional problem.

---

## 1. Three-Dimensional Linear Elasticity

Let $\Omega_{3D}\subset\mathbb{R}^3$ denote an elastic body.
The static equilibrium equations read

$$
\nabla\cdot\boldsymbol{\sigma}+\boldsymbol{b}=0
\quad\text{in }\Omega_{3D}.
$$

In the absence of body forces ($\boldsymbol{b}=0$),

$$
\nabla\cdot\boldsymbol{\sigma}=0.
$$

The material is assumed to be linear, isotropic, and homogeneous, with
the constitutive relation

$$
\boldsymbol{\sigma}
=2\mu\,\boldsymbol{\varepsilon}(\boldsymbol{u})
+\lambda\,\mathrm{tr}\!\left(\boldsymbol{\varepsilon}(\boldsymbol{u})\right)\mathbf{I},
$$

where the infinitesimal strain tensor is defined as

$$
\boldsymbol{\varepsilon}(\boldsymbol{u})
=\frac{1}{2}\left(
\nabla\boldsymbol{u}
+(\nabla\boldsymbol{u})^{T}
\right).
$$

---

## 2. Axisymmetric Assumptions

Cylindrical coordinates $(r,\theta,z)$ are introduced.
The following assumptions are made:

- all fields are independent of $\theta$,
- the circumferential displacement vanishes: $u_\theta=0$.

The displacement field therefore reduces to

$$
\boldsymbol{u}(r,z)=(u_r(r,z),\,0,\,u_z(r,z)).
$$

The computational domain is the meridian cross-section

$$
\Omega\subset\mathbb{R}^2_{(r,z)}.
$$

---

## 3. Axisymmetric Strong Form

**Radial equilibrium**

$$
\partial_r\sigma_{rr}
+\partial_z\sigma_{rz}
+\frac{\sigma_{rr}-\sigma_{\theta\theta}}{r}
=0
\quad\text{in }\Omega.
$$

**Axial equilibrium**

$$
\partial_r\sigma_{rz}
+\partial_z\sigma_{zz}
+\frac{\sigma_{rz}}{r}
=0
\quad\text{in }\Omega.
$$

---

## 4. Axisymmetric Weak Formulation

Let $v_r$ and $v_z$ denote admissible test functions.
Assuming zero traction boundary conditions,
integration over the circumferential direction yields

$$
\mathrm{d}\Omega_{3D}=2\pi r\,\mathrm{d}\Omega.
$$

### Radial equilibrium (weak form)

$$
\int_\Omega
\left(
\sigma_{rr}\,\partial_r v_r
+\sigma_{rz}\,\partial_z v_r
\right)(2\pi r)\,\mathrm{d}\Omega
+\int_\Omega
\sigma_{\theta\theta}\,v_r\,(2\pi)\,\mathrm{d}\Omega
=0.
$$

### Axial equilibrium (weak form)

$$
\int_\Omega
\left(
\sigma_{rz}\,\partial_r v_z
+\sigma_{zz}\,\partial_z v_z
\right)(2\pi r)\,\mathrm{d}\Omega
=0.
$$

---

## 5. Stress–Strain Relations in Axisymmetry

$$
\begin{aligned}
\sigma_{rr}
&=2\mu\,\partial_r u_r
+\lambda\left(
\partial_r u_r+\frac{u_r}{r}+\partial_z u_z
\right),\\
\sigma_{\theta\theta}
&=2\mu\,\frac{u_r}{r}
+\lambda\left(
\partial_r u_r+\frac{u_r}{r}+\partial_z u_z
\right),\\
\sigma_{zz}
&=2\mu\,\partial_z u_z
+\lambda\left(
\partial_r u_r+\frac{u_r}{r}+\partial_z u_z
\right),\\
\sigma_{rz}
&=\mu\left(
\partial_z u_r+\partial_r u_z
\right).
\end{aligned}
$$

---

## 6. Compact Bilinear Form

### Radial test function

$$
\begin{aligned}
\int_\Omega (\nabla v_r)^{T}K_{rr}(r)(\nabla u_r)\,\mathrm{d}\Omega
&+\int_\Omega (\nabla v_r)^{T}K_{rz}(r)(\nabla u_z)\,\mathrm{d}\Omega\\
&-\lambda\int_\Omega u_z\,\partial_z v_r\,\mathrm{d}\Omega
+\int_\Omega \frac{\lambda+2\mu}{r}u_r v_r\,\mathrm{d}\Omega
=0.
\end{aligned}
$$

### Axial test function

$$
\begin{aligned}
\int_\Omega (\nabla v_z)^{T}K_{zr}(r)(\nabla u_r)\,\mathrm{d}\Omega
&+\int_\Omega (\nabla v_z)^{T}K_{zz}(r)(\nabla u_z)\,\mathrm{d}\Omega\\
&+\lambda\int_\Omega u_r\,\partial_z v_z\,\mathrm{d}\Omega
=0.
\end{aligned}
$$

---

## 7. Remark

This formulation is mathematically equivalent to the original
three-dimensional linear elasticity problem under the stated
axisymmetry assumptions.
It is directly suitable for two-dimensional finite element
discretizations using axisymmetric weighting.
