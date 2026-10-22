<h1 align="center">Full formulation</h1>

<p align="center"><b>English</b> · <a href="FORMULACAO.pt-BR.md">Português</a> · <a href="../README.md">back to the README</a></p>

Formulation of Ferreira (2024, *Physica B* 695, 416494; MRS Meeting 2026) of the thermal-field tensor in nucleation,
with the Gurtin–Murdoch surface extended to elasto-plasticity in this repository, and its propagation to the secondary
dendrite arm spacing (SDAS) through the modified Rappaz–Boettinger (MRB) model. Sections 1–4 follow the form given by
the author; sections 5–8 transcribe what the scripts in `model/` compute.

## Contents

1. [Thermal field from the first law](#1-thermal-field-from-the-first-law)
2. [Thermal-field tensor](#2-thermal-field-tensor)
3. [Nucleation: free energy, critical radii and tensors](#3-nucleation-free-energy-critical-radii-and-tensors)
4. [Nucleation rate](#4-nucleation-rate)
5. [Surface energy: Gurtin–Murdoch](#5-surface-energy-gurtinmurdoch)
6. [Elasto-plastic surface](#6-elasto-plastic-surface)
7. [Numerical solution against the thermal gradient](#7-numerical-solution-against-the-thermal-gradient)
8. [Secondary dendrite arm spacing (MRB)](#8-secondary-dendrite-arm-spacing-mrb)
9. [Symbols](#9-symbols)

**Notation.** $`\mathbf{T}`$ (bold) is the thermal field, in K, and $`\nabla\mathbf{T}`$ its gradient; $`T`$ is the temperature
and $`\Delta T`$ the undercooling. $`\boldsymbol{\Gamma}`$ is the thermal-field tensor (generalized Gibbs–Thomson coefficient)
and $`\mathbf{A}`$ the area of the new surface, normal to the field $`\mathbf{T}`$.

---

## 1. Thermal field from the first law

By the first law the energy $`E`$ is a point function, while $`Q`$ and $`W`$ depend on the path:

```math
\delta Q - \delta W = dE
```

Taking $`dE \equiv \partial E_{,i}`$ with respect to the primitive variable $`T`$ defines two thermal-field gradients:

```math
\nabla\mathbf{T}_{\mathrm{pi}} = \nabla\left[\frac{\partial T}{\delta Q}\,(Q-W)\right] = \nabla\left[\frac{\partial T}{\delta Q}\,E\right]
\qquad\text{(path independent)}
```

```math
\nabla\mathbf{T}_{\mathrm{pd}} = \nabla\left[\frac{\partial T}{\delta Q}\,Q\right]
= \nabla\mathbf{T}_{\mathrm{pi}} + \nabla\left[\frac{\partial T}{\delta Q}\,W\right]
\qquad\text{(path dependent)}
```

Since $`E`$ and $`W`$ are independent, the cross terms are neglected:

```math
\nabla\mathbf{T} = \nabla\left[\frac{\partial T}{\delta Q}\,(E+W)\right]
\cong \nabla\left[\frac{\partial T}{\partial E}\,E\right] + \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

With $`E = U + K + E_P + \Sigma_{other}`$ and the internal energy split into its forms,

```math
\nabla\mathbf{T} \cong
\nabla\left[\frac{\partial T}{\partial U}\left(U_{transl} + U_{vib} + U_{rot} + U_{elect} + U_{mag} + U_{other}\right)\right]
+ \nabla\left[\frac{\partial T}{\partial K}\,K\right]
+ \nabla\left[\frac{\partial T}{\partial E_P}\,E_P\right]
+ \nabla\left[\frac{\partial T}{\partial \Sigma_{other}}\,\Sigma_{other}\right]
+ \nabla\left[\frac{\partial T}{\delta W}\,W\right]
```

Only the gradient of the thermal field displaces the solid–liquid equilibrium. In a transient solidification
experiment, $`\nabla\mathbf{T}`$ at the liquidus isotherm is $`G_L = \dot T / V_L = T_R / V_L`$.

## 2. Thermal-field tensor

Every newly created or deformed surface, of area $`\mathbf{A}`$, is normal to the thermal field:

```math
\boldsymbol{\Gamma} = \mathbf{A}\cdot\nabla\mathbf{T} \cong \mathbf{A}\cdot\nabla\left[\frac{\partial T}{\partial E}\,E\right]
```

## 3. Nucleation: free energy, critical radii and tensors

Total free energy of the nucleus (cap of angle $`\theta`$), with $`\Delta G_V = \Delta S_V\,\Delta T`$ and
$`\Delta G_S = \gamma_{SL}`$, where $`\Delta S_V`$, $`\Delta T`$, $`\gamma_{SL}`$ and $`\theta`$ depend on the radius:

```math
\Delta G = \left(\frac{1}{3}\pi r^3\,\Delta S_V\,\Delta T + \pi r^2\,\gamma_{SL}\right) f(\theta),
\qquad f(\theta) = 2 - 3\cos\theta + \cos^3\theta
```

The maximum with respect to $`r`$, $`\partial\Delta G/\partial r = 0`$, gives the quadratic equation

```math
\left[\left(\frac{\partial \Delta S_V}{\partial r}\Delta T + \Delta S_V\frac{\partial \Delta T}{\partial r}\right) f(\theta)
+ \Delta S_V\,\Delta T\,\frac{\partial f(\theta)}{\partial r}\right] r_C^2
+ 3\left[\left(\Delta S_V\,\Delta T + \frac{\partial \gamma_{SL}}{\partial r}\right) f(\theta)
+ \gamma_{SL}\frac{\partial f(\theta)}{\partial r}\right] r_C + 6\gamma_{SL}\, f(\theta) = 0
```

**Second-order critical radius**, in terms of the bulk, surface and configurational free energies
($`\Delta G_C = \gamma_{SL}\,(\partial f/\partial r)/f`$):

```math
r_C^{2nd} = -\frac{3}{2}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{2}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta G_V}{\partial r}}
```

**First-order critical radius** and classical nucleation theory (CNT), recovered when
$`\partial\Delta G_S/\partial r = \Delta G_C = 0`$:

```math
r_C^{1st} = -\frac{2\,\Delta G_S}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C},
\qquad r_C^{CNT} = -\frac{2\,\Delta G_S}{\Delta G_V}
```

**Second-order thermal-field tensor:**

```math
\boldsymbol{\Gamma}^{2nd} = -\frac{3}{4}\,
\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial \gamma_{SL}}{\partial r}\right) f(\theta) + \gamma_{SL}\dfrac{\partial f(\theta)}{\partial r}}
{\dfrac{1}{\Delta T}\dfrac{\partial\left(\Delta S_V\,\Delta T\, f(\theta)\right)}{\partial r}}
= -\frac{3}{4}\,\frac{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}{\dfrac{\partial \Delta S_V}{\partial r}}
= \mathbf{A}\cdot\nabla\mathbf{T}
```

**First-order thermal-field tensor**, with entropies $`S \equiv G/\Delta T`$:

```math
\boldsymbol{\Gamma}^{1st} = -\frac{\Delta G_S\,\Delta T}{\Delta G_V + \dfrac{\partial \Delta G_S}{\partial r} + \Delta G_C}
= -\frac{\Delta G_S}{\Delta S_V + \dfrac{\partial \Delta S_S}{\partial r} + \Delta S_C}
```

Both have the Gibbs–Thomson form:

```math
\boldsymbol{\Gamma}^{1st} = \tfrac{1}{2}\,\Delta T\, r_C^{1st}, \qquad \boldsymbol{\Gamma}^{2nd} = \tfrac{1}{2}\,\Delta T\, r_C^{2nd}
```

$`\Gamma^{1st}`$ contains the surface entropy $`\Delta S_S = (1/\Delta T)\,\partial\gamma_{SL}/\partial r`$, i.e. the derivative of
$`\gamma_{SL}`$; $`\Gamma^{2nd}`$ depends on $`\gamma_{SL}`$.

## 4. Nucleation rate

```math
I = \frac{D\,A}{\lambda^4}\,\frac{N}{V}\,\exp\left(\frac{\Delta G_C}{\Delta G_{C,Eq}}\right)
```

with $`D`$ the self-diffusion coefficient, $`N/V`$ the density of modes of the new phase, $`\lambda`$ the jump distance to
the new surface and $`A`$ the area of the new phase.

## 5. Surface energy: Gurtin–Murdoch

On the nucleus surface $`\mathcal{S}`$, with surface gradient $`\nabla_S`$ and tangent projector
$`\mathbf{P} = \mathbf{I} - \mathbf{n}\otimes\mathbf{n}`$ (Gurtin and Murdoch, 1978; Ferreira, 2022):

```math
\nabla_S \cdot \boldsymbol{\Sigma} = \boldsymbol{\Theta}\,\mathbf{n} + \rho_0\,\ddot{\mathbf{u}}
```

```math
\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E} + (\lambda_0 - \sigma)\,(\mathrm{tr}\,\mathbf{E})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u}
\qquad\text{on } \mathcal{S},
\qquad \mathbf{E} = \tfrac{1}{2}\left(D\mathbf{u} + D\mathbf{u}^{T}\right)
```

with $`\sigma`$ the residual surface tension, $`\rho_0`$ the surface mass density and $`\lambda_0`$, $`\mu_0`$ the surface
Lamé moduli. Surface stress and surface energy are related by the Shuttleworth equation,

```math
\Sigma_{\alpha\beta} = \gamma\,\delta_{\alpha\beta} + \frac{\partial \gamma}{\partial \varepsilon_{\alpha\beta}},
```

and the surface free energy in metric form (Ricci scalar $`R = 2/r_0^2`$ on the 2-sphere) is

```math
\Delta G_S = \gamma_0\,\Phi\!\left(r_0^2 R\right) - \frac{1}{2A}\int_{\mathcal{S}} \Sigma_{ij}\,g^{ij}\,dA .
```

**Spherical solution used in the scripts.** With the bulk Lamé constants $`\lambda`$, $`\mu`$:

```math
y(r) = (3\lambda + 2\mu)\,r + 2(\lambda_0 + \mu_0), \qquad y_0 = y(r_{eq})
```

```math
s(r) = -2\,\sigma_0 \ln\frac{y_0}{y}, \qquad \tau(r) = -\sigma_0 + s(r), \qquad
\gamma_{SL}(r) = \gamma_0\left(\frac{r_{eq}}{r}\right)^{2} - s(r)
```

```math
\frac{\partial \gamma_{SL}}{\partial r} = -\frac{2\gamma_0\,r_{eq}^2}{r^3} - \frac{2\sigma_0}{r\,(1+\alpha)},
\qquad \alpha = \frac{2\left(\sigma_0 + 2(\lambda_0+\mu_0)\right)}{r\,(3\lambda + 2\mu)}
```

$`\sigma_0 < 0`$ is a compressive force; in the scripts `sigma_0` holds its magnitude, so that $`\tau(r_{eq}) = -\sigma_0`$.
$`\gamma_0`$, $`\sigma_0`$ and $`r_{eq}`$ are equilibrium values (no deformation of the created area). With
$`E_A = 2\ln(y/y_0)`$, the areal strain of the surface, the elastic solution is $`s = C_A E_A`$ with $`C_A = \sigma_0`$.

## 6. Elasto-plastic surface

**Kinematics** (tangent plane) and the Gurtin–Murdoch law with the elastic part:

```math
\mathbf{E} = \mathbf{E}^{e} + \mathbf{E}^{p}, \qquad \mathbf{E}^{p} = \mathbf{P}\,\mathbf{E}^{p}\,\mathbf{P},
\qquad
\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E}^{e} + (\lambda_0 - \sigma)\,(\mathrm{tr}\,\mathbf{E}^{e})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u}
```

**Yield on the excess stress** $`\tilde{\boldsymbol{\Sigma}} = \boldsymbol{\Sigma} - \sigma\,\mathbf{P}`$ (the residual stress alone
never yields):

```math
f_A = \left|\mathrm{tr}\,\tilde{\boldsymbol{\Sigma}}\right| - 2\,\sigma_y^{s}(\alpha_A) \le 0, \qquad \sigma_y^{s} = \sigma_{y0}^{s} + H_A\,\alpha_A
\qquad\text{(areal: coherent → semicoherent interface)}
```

```math
f_D = \left\lVert \mathrm{dev}_S\, \tilde{\boldsymbol{\Sigma}} \right\rVert - \sigma_{yD}(\alpha_D) \le 0, \qquad \sigma_{yD} = \sigma_{yD0} + H_D\,\alpha_D
\qquad\text{(J}_2\text{ in the tangent plane)}
```

**Associative flow rule and Kuhn–Tucker conditions:**

```math
\dot{\mathbf{E}}^{p} = \dot{\lambda}_A\,\mathrm{sign}\!\left(\mathrm{tr}\,\tilde{\boldsymbol{\Sigma}}\right)\tfrac{1}{2}\mathbf{P}
+ \dot{\lambda}_D\,\frac{\mathrm{dev}_S\,\tilde{\boldsymbol{\Sigma}}}{\lVert \mathrm{dev}_S\,\tilde{\boldsymbol{\Sigma}} \rVert},
\qquad \dot{\alpha}_k = \dot{\lambda}_k, \qquad \dot{\lambda}_k \ge 0,\quad f_k \le 0,\quad \dot{\lambda}_k f_k = 0
```

On a sphere $`\mathbf{E} = \varepsilon\,\mathbf{P}`$ is equibiaxial, $`\mathrm{dev}_S\,\mathbf{E} = \mathbf{0}`$, and only the areal mechanism acts.

**Return mapping** (step $`n \rightarrow n+1`$, $`r`$ decreasing from $`r_{eq}`$):

```math
s^{tr} = C_A\left(E_A - E_{A,n}^{p}\right), \qquad f^{tr} = |s^{tr}| - \left(\sigma_{y0}^{s} + H_A\,\alpha_{A,n}\right)
```

```math
f^{tr} > 0:\quad \Delta\lambda = \frac{f^{tr}}{C_A + H_A},\quad s = s^{tr} - C_A\,\Delta\lambda\,\mathrm{sign}(s^{tr}),\quad
E_{A,n+1}^{p} = E_{A,n}^{p} + \Delta\lambda\,\mathrm{sign}(s^{tr}),\quad \alpha_{A,n+1} = \alpha_{A,n} + \Delta\lambda
```

with consistent tangent $`\partial s/\partial E_A = C_A H_A/(C_A + H_A)`$ and the closed form for monotonic loading

```math
s = \mathrm{sign}(E_A)\left[\sigma_{y0}^{s} + \frac{H_A\,C_A}{C_A + H_A}\left(|E_A| - \frac{\sigma_{y0}^{s}}{C_A}\right)\right],
\qquad |E_A| > \sigma_{y0}^{s}/C_A .
```

**Plastic area and the equilibrium reference.** The plastic area moves the reference at which $`\gamma_0`$, $`\sigma_0`$
and $`r_{eq}`$ were defined:

```math
\gamma_{SL} = \gamma_0\left(\frac{r_{eq,p}}{r}\right)^{2} - s, \qquad r_{eq,p} = r_{eq}\,\exp\!\left(E_A^{p}/2\right),
\qquad
\frac{\partial \gamma_{SL}}{\partial r} = -\frac{2\gamma_0\,r_{eq,p}^2}{r^3} + \frac{\gamma_0\,r_{eq,p}^2}{r^2}\frac{\partial E_A^{p}}{\partial r} - \frac{\partial s}{\partial E_A}\frac{\partial E_A}{\partial r}
```

With perfect plasticity ($`H_A = 0`$), after yield:

```math
\gamma_{SL} = \gamma_0\,e^{\sigma_y^{s}/C_A} + \sigma_y^{s}, \qquad \Gamma^{2nd} = \frac{\gamma_{SL}}{|\Delta S_V|}
\quad\text{(independent of the radius and of the thermal gradient)}
```

**Energy of the plastic surface.** Free energy $`\psi = \psi^{e}(\mathbf{E}^{e}) + \psi^{p}(\alpha_A)`$, with
$`\psi^{p} = \tfrac{1}{2}H_A\alpha_A^2`$ and hardening force $`K = \partial\psi^{p}/\partial\alpha_A = H_A\alpha_A`$. By the
Clausius–Duhem inequality, $`D = s\,\dot{E}_A^{p} - K\dot{\alpha}_A = \sigma_{y0}^{s}\,\dot{\alpha}_A \ge 0`$, and

```math
w_p = \int \left(\sigma_{y0}^{s} + H_A\,\alpha_A\right) d\alpha_A = \psi^{p} + q_{diss}, \qquad q_{diss} = \sigma_{y0}^{s}\,\alpha_A .
```

| Surface state (`GM_DISSIPATION_MODE`) | $`\gamma_{SL}`$ |
|---|---|
| loaded, elastic part (`stored`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s`$ |
| loaded + stored energy (`thermo`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s + \psi^{p}`$ |
| loaded, total (`barrier`) | $`\gamma_0\,(r_{eq,p}/r)^2 - s + w_p`$ |
| after the elastic restoration (`restored`) | $`\gamma_0`$ |
| after the restoration + plastic work (`restored_barrier`) | $`\gamma_0 + w_p`$ |

## 7. Numerical solution against the thermal gradient

The script steps the radius from equilibrium, $`r_i = r_{eq} - i\,\Delta r`$ with $`\Delta r = 0.005\,r_{eq}`$, and at each
step computes $`\gamma_{SL}`$ and $`\partial\gamma_{SL}/\partial r`$ (sections 5–6) and then:

**Undercooling.** The root the script obtains from `f_DT` coincides with the capillary relation

```math
\Delta T = \frac{2\,\gamma_{SL}}{r\,|\Delta S_V|}, \qquad \Delta S_V = -\frac{\Delta H\,\rho_s}{T_L}
```

**Thermal gradient** (homogeneous nucleus, $`\theta = \pi`$, $`f = 4`$):

```math
\nabla T = \frac{\Delta T}{8\pi r}
```

**Second-order angle** $`\theta_2`$: root of

```math
-\frac{3}{2}\,\frac{\left(\Delta S_V\,\Delta T + \dfrac{\partial\gamma_{SL}}{\partial r}\right)\dfrac{f(\theta)}{4} + \gamma_{SL}\,f'_r(\theta)}
{\dfrac{\partial(\Delta S_V\Delta T)}{\partial r}\dfrac{f(\theta)}{4} + \Delta S_V\,\Delta T\,f'_r(\theta)} - r_{het}(r,\theta) = 0,
\qquad
\frac{\partial(\Delta S_V\Delta T)}{\partial r} = -\frac{3}{2r}\left(\Delta S_V\,\Delta T + \frac{\partial\gamma_{SL}}{\partial r}\right)
```

with $`f'_r(\theta) = -3f(\theta) - (1-\cos\theta)(2-\cos\theta-\cos^2\theta)`$ (function `dfthetadr` of the script) and the cap radius

```math
r_{het}(r,\theta) = r\,\exp\!\left[\int_{\pi}^{\theta}\frac{\sin\theta'\,(1+\cos\theta')}{2-\cos\theta'-\cos^2\theta'}\,d\theta'\right]
```

**Heterogeneous second-order tensor** (column `GT_het_2nd`):

```math
\nabla T_{het} = \frac{\Delta T}{4\pi\,r_{het}\,(1-\cos\theta_2)}, \qquad
\Gamma^{2nd}_{het} = 2\pi r^2\,(1-\cos\theta_2)\,\nabla T_{het} = \frac{r^2}{r_{het}}\,\frac{\Delta T}{2}
```

Since $`\theta_2 \approx \pi`$, $`r_{het} \approx r`$ and $`\Gamma^{2nd}_{het} \approx r\,\Delta T/2 = \gamma_{SL}/|\Delta S_V| = 4\pi r^2\,\nabla T`$.

**State at an experimental gradient.** The quantities are interpolated linearly on the curve $`\nabla T(r)`$ at the
measured gradient (mean or local).

## 8. Secondary dendrite arm spacing (MRB)

For $`n`$ solutes (Al-3Cu-5Nb-0.1Fe: Cu, Nb, Fe), the modified Rappaz–Boettinger model is solved by fixed-point
iteration on SDAS:

```math
\mathrm{SDAS} = 5.5\,\left(M\,t_{SL}\right)^{1/3}, \qquad
M = \bar{\beta}\left(-\frac{\Gamma}{p_1}\right)\ln\frac{s_2}{t_3}
```

```math
p_1 = \sum_i \frac{m_{L,i}(1-k_{ef,i})(C_{f,i}-C_{0,i})}{D_{L,i}}, \qquad
s_2 = \sum_i \frac{m_{L,i}(1-k_{ef,i})\,C_{f,i}}{D_{L,i}}, \qquad
t_3 = \sum_i \frac{m_{L,i}(1-k_{ef,i})\,C_{0,i}}{D_{L,i}}
```

**Back-diffusion in the solid:**

```math
Fo_i = \frac{D_{S,i}\,t_{SL}}{(\mathrm{SDAS}/2)^2}, \qquad
\gamma_i = \frac{4\,Fo_i\,k_{0,i}}{4\,Fo_i\,k_{0,i} + 1}, \qquad
\beta_i = \frac{2\,\gamma_i\,Fo_i}{2\,Fo_i + \gamma_i}, \qquad
\bar{\beta} = \sum_i w_i\,\beta_i, \quad w_i = \frac{C_{0,i}}{\sum_j C_{0,j}}
```

(in the script $`\beta_3`$ uses the Fourier number of solute 2: $`\beta_3 = 2\gamma_3 Fo_2/(2Fo_3 + \gamma_3)`$).

**Effective partition** with the mean diffusion length $`\bar\delta_i`$ (mean of $`D_{L,i}/V`$ over the velocity range):

```math
k_{ef,i} = \frac{k_{0,i}}{k_{0,i} + (1-k_{0,i})\exp\!\left(-V\,\bar\delta_i / D_{L,i}\right)}, \qquad
D_{L,i} = D_{0L,i}\,e^{-Q_{L,i}/RT_L}, \quad D_{S,i} = D_{0S,i}\,e^{-Q_{S,i}/RT_S}
```

The original RB model is recovered with $`\bar\beta = 1`$ and $`k_{ef} = k_0`$.

**Casting laws** (Al-3Cu-5Nb-0.1Fe, horizontal; P in mm):

```math
V_L = 3.4\,P^{-0.56}\ \mathrm{mm\,s^{-1}}, \qquad t_{SL} = 4.93116\,P^{0.87606}\ \mathrm{s}, \qquad
G_L(P) = 10405.54 - 200.22\,P + 1.04692\,P^2\ \mathrm{K\,m^{-1}}
```

**SDAS explicit in the thermal gradient.** Inverting $`G_L(P)`$ (monotonic for $`P < 95.6`$ mm):

```math
P(\nabla T_L) = \frac{200.22 - \sqrt{40088.05 - 4.18768\,(10405.54 - \nabla T_L)}}{2.09384}, \qquad
t_{SL}(\nabla T_L) = 4.93116\,P(\nabla T_L)^{0.87606}
```

```math
\mathrm{SDAS}(\nabla T_L) = 5.5\left[M\,4.93116\,P(\nabla T_L)^{0.87606}\right]^{1/3}
```

and, with $`t_{SL} = \Delta T_{SL}/(\nabla T_L\,V_L)`$, the form with exponent $`-1/3`$:

```math
\mathrm{SDAS} = 5.5\left[\frac{M\,\Delta T_{SL}}{\nabla T_L\,V_L}\right]^{1/3} = K\left(\nabla T_L\,V_L\right)^{-1/3}
```

**Local tensor.** With $`\Gamma`$ computed at the gradient measured at each position, $`\Gamma(P) = \Gamma^{2nd}(\nabla T_L(P))`$
enters $`M`$ point by point. On the elastic surface $`\Gamma^{2nd} \approx 4\pi r^2\,\nabla T`$ is proportional to the
gradient; on the elasto-plastic surface with the updated reference and $`H_A = 0`$, $`\Gamma^{2nd}`$ is constant after yield.

## 9. Symbols

| Symbol | Meaning | Unit |
|---|---|---|
| $`\mathbf{T}`$, $`\nabla\mathbf{T}`$, $`G_L`$ | thermal field, its gradient, gradient at the liquidus isotherm | K, K m⁻¹ |
| $`\Delta T`$ | undercooling | K |
| $`\Delta S_V`$ | volumetric entropy of fusion | J m⁻³ K⁻¹ |
| $`\gamma_{SL}`$, $`\gamma_0`$ | solid–liquid surface energy and its equilibrium value | J m⁻² |
| $`\sigma_0`$, $`\tau`$, $`s`$ | residual surface tension, surface tension, excess surface stress | N m⁻¹ |
| $`\lambda_0`$, $`\mu_0`$; $`\lambda`$, $`\mu`$ | Lamé moduli of the surface; of the bulk | N m⁻¹; Pa |
| $`r`$, $`r_C`$, $`r_{eq}`$, $`r_{het}`$ | radius, critical radius, equilibrium radius, cap radius | m |
| $`\theta`$, $`f(\theta)`$ | contact angle and geometric factor | rad, – |
| $`\Gamma^{1st}`$, $`\Gamma^{2nd}`$ | first- and second-order thermal-field tensors | m K |
| $`E_A`$, $`E_A^p`$, $`\alpha_A`$ | total and plastic areal strain, accumulated plastic strain | – |
| $`\sigma_{y0}^s`$, $`H_A`$, $`C_A`$ | areal yield stress, hardening, areal modulus | N m⁻¹ |
| $`w_p`$, $`\psi^p`$, $`q_{diss}`$ | plastic work, stored energy, dissipated heat | J m⁻² |
| $`t_{SL}`$, $`V_L`$, $`T_R`$ | local solidification time, tip velocity, cooling rate | s, m s⁻¹, K s⁻¹ |
| $`m_{L,i}`$, $`k_{0,i}`$, $`k_{ef,i}`$ | liquidus slope, equilibrium and effective partition coefficients | K wt.%⁻¹, – |
| $`C_{0,i}`$, $`C_{f,i}`$ | initial and final liquid concentration | wt.% |
| $`D_{L,i}`$, $`D_{S,i}`$ | diffusivity in the liquid and in the solid | m² s⁻¹ |
| $`Fo_i`$, $`\beta_i`$, $`\bar\beta`$ | Fourier number, back-diffusion parameter and its weighted mean | – |
