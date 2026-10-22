# gm-elastoplastic-nucleation

Elasto-plastic extension of the Gurtin–Murdoch (GM) surface stress in the thermal-field-tensor
nucleation model of Ferreira (Physica B 695, 416494, 2024), applied to the FCC_A1 phase of the
Al-3Cu-5Nb-0.1Fe alloy and propagated to the secondary dendrite arm spacing (SDAS) through the
modified Rappaz–Boettinger (MRB) model. A small CPU molecular-dynamics (MD) study of the
solid–liquid interface stress of Al is included as a check of the surface parameters.

Author: I.L. Ferreira, Faculty of Mechanical Engineering, Federal University of Pará (UFPA), Brazil.

## Contents

| Folder | What it holds |
|---|---|
| `model/` | Nucleation model with the elasto-plastic GM surface (`..._ElastoPlastic_GM.py`) and the two SDAS scripts (`..._GT_point_to_point.py`, `..._SDAS_vs_GL.py`); helpers `derivative_shim.py`, `midpoint.py` |
| `md_cpu/` | LAMMPS inputs and analysis for the planar interface stress and the Laplace pressure of a critical seed (`E5_results.md`) |
| `results/` | CSV tables and figures produced by the scripts and the MD runs |

## Theory

### 1. Gurtin–Murdoch surface (Gurtin and Murdoch, 1978; Ferreira, 2022, 2024)

On the nucleus surface $\mathcal{S}$, with surface gradient $\nabla_S$ and tangent projector
$\mathbf{P} = \mathbf{I} - \mathbf{n}\otimes\mathbf{n}$:

$$\nabla_S \cdot \boldsymbol{\Sigma} = \boldsymbol{\Theta}\,\mathbf{n} + \rho_0\,\ddot{\mathbf{u}}$$

$$\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E} + (\lambda_0 - \sigma)\,(\operatorname{tr}\mathbf{E})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u} \qquad \text{on } \mathcal{S}$$

$$\mathbf{E} = \tfrac{1}{2}\left(D\mathbf{u} + D\mathbf{u}^{T}\right)$$

where $\boldsymbol{\Theta}$, $\mathbf{E}$ and $\mathbf{u}$ are the stress, strain and displacement fields,
$\sigma$ and $\rho_0$ the residual surface tension and surface mass density, and $\lambda_0$, $\mu_0$ the
surface Lamé moduli. Surface stress and surface energy are related by the Shuttleworth equation,

$$\Sigma_{\alpha\beta} = \gamma\,\delta_{\alpha\beta} + \frac{\partial \gamma}{\partial \varepsilon_{\alpha\beta}},$$

and the surface free energy of the nucleus is

$$\Delta G_S = \gamma_0\,\Phi\!\left(r_0^2 R\right) - \frac{1}{2A}\int_{\mathcal{S}} \Sigma_{ij}\,g^{ij}\,dA ,$$

with $R$ the Ricci scalar ($R = 2/r_0^2$ for a 2-sphere), $g^{ij}$ the metric of the deformed surface and
$r_0$ the radius at the equilibrium (undeformed) state.

**Spherical solution used in the scripts.** With bulk Lamé constants $\lambda$, $\mu$,

$$y(r) = (3\lambda + 2\mu)\,r + 2(\lambda_0 + \mu_0), \qquad y_0 = y(r_{eq}),$$

the excess surface stress, surface tension and surface energy are

$$s(r) = -2\,\sigma_0 \ln\frac{y_0}{y}, \qquad \tau(r) = -\sigma_0 + s(r), \qquad \gamma_{SL}(r) = \gamma_0\left(\frac{r_{eq}}{r}\right)^{2} - s(r).$$

$\sigma_0 < 0$ is a compressive force; in the scripts `sigma_0` holds its magnitude, so that
$\tau(r_{eq}) = -\sigma_0$. $\gamma_0$, $\sigma_0$ and $r_{eq}$ are equilibrium values: no deformation of the
created area. Writing $E_A = 2\ln(y/y_0)$ for the areal strain of the surface (partitioned with the bulk),
the elastic solution is $s = C_A E_A$ with areal modulus $C_A = \sigma_0$.

### 2. Elasto-plastic proposal

**Kinematics.** Additive split of the surface strain in the tangent plane,

$$\mathbf{E} = \mathbf{E}^{e} + \mathbf{E}^{p}, \qquad \mathbf{E}^{p} = \mathbf{P}\,\mathbf{E}^{p}\,\mathbf{P},$$

and the Gurtin–Murdoch law is evaluated with the elastic part:

$$\boldsymbol{\Sigma} = \sigma\,\mathbf{I} + 2\,(\mu_0 - \sigma)\,\mathbf{E}^{e} + (\lambda_0 - \sigma)\,(\operatorname{tr}\mathbf{E}^{e})\,\mathbf{I} + \sigma\,\nabla_S\mathbf{u}.$$

**Yield on the excess stress.** With $\tilde{\boldsymbol{\Sigma}} = \boldsymbol{\Sigma} - \sigma\,\mathbf{P}$, so that the
residual tension alone never yields, two decoupled mechanisms:

$$f_A = \left|\operatorname{tr}\tilde{\boldsymbol{\Sigma}}\right| - 2\,\sigma_y^{s}(\alpha_A) \le 0, \qquad \sigma_y^{s} = \sigma_{y0}^{s} + H_A\,\alpha_A \quad \text{(areal: coherent} \rightarrow \text{semicoherent interface)}$$

$$f_D = \left\lVert \operatorname{dev}_S \tilde{\boldsymbol{\Sigma}} \right\rVert - \sigma_{yD}(\alpha_D) \le 0, \qquad \sigma_{yD} = \sigma_{yD0} + H_D\,\alpha_D \quad \text{(J}_2\text{ in the tangent plane)}$$

**Flow rule and hardening (associative), with Kuhn–Tucker conditions:**

$$\dot{\mathbf{E}}^{p} = \dot{\lambda}_A\,\operatorname{sign}\!\left(\operatorname{tr}\tilde{\boldsymbol{\Sigma}}\right)\tfrac{1}{2}\mathbf{P} + \dot{\lambda}_D\,\frac{\operatorname{dev}_S\tilde{\boldsymbol{\Sigma}}}{\lVert \operatorname{dev}_S\tilde{\boldsymbol{\Sigma}} \rVert}, \qquad \dot{\alpha}_k = \dot{\lambda}_k, \qquad \dot{\lambda}_k \ge 0,\; f_k \le 0,\; \dot{\lambda}_k f_k = 0 .$$

For a spherical nucleus the surface strain is equibiaxial, $\mathbf{E} = \varepsilon\,\mathbf{P}$, so
$\operatorname{dev}_S\mathbf{E} = \mathbf{0}$ and the J$_2$ mechanism never activates (checked at run time);
only the areal mechanism changes the result.

**Return mapping (sphere, step $n \rightarrow n+1$, $r$ decreasing from $r_{eq}$):**

$$s^{tr} = C_A\left(E_A - E_{A,n}^{p}\right), \qquad f^{tr} = |s^{tr}| - \left(\sigma_{y0}^{s} + H_A\,\alpha_{A,n}\right)$$

$$\text{if } f^{tr} > 0:\quad \Delta\lambda = \frac{f^{tr}}{C_A + H_A},\quad s = s^{tr} - C_A\,\Delta\lambda\,\operatorname{sign}(s^{tr}),\quad E_A^{p} \mathrel{+}= \Delta\lambda\,\operatorname{sign}(s^{tr}),\quad \alpha_A \mathrel{+}= \Delta\lambda$$

with consistent tangent $\partial s/\partial E_A = C_A H_A/(C_A + H_A)$ and the closed form, used as a check,

$$s = \operatorname{sign}(E_A)\left[\sigma_{y0}^{s} + \frac{H_A\,C_A}{C_A + H_A}\left(|E_A| - \frac{\sigma_{y0}^{s}}{C_A}\right)\right] \quad \text{for } |E_A| > \sigma_{y0}^{s}/C_A .$$

**Plastic area and the equilibrium reference.** The plastic area does not compress the surface energy; it
moves the reference at which $\gamma_0$, $\sigma_0$ and $r_{eq}$ were defined:

$$\gamma_{SL} = \gamma_0\left(\frac{r_{eq,p}}{r}\right)^{2} - s, \qquad r_{eq,p} = r_{eq}\,\exp\!\left(E_A^{p}/2\right)$$

(option `GM_PLASTIC_REFERENCE=1`; with $E_A^p = 0$ the elastic solution is recovered). For perfect
plasticity, $H_A = 0$, the surface energy after yield is independent of the radius and of the thermal
gradient:

$$\gamma_{SL} = \gamma_0\,e^{\sigma_y^{s}/C_A} + \sigma_y^{s}, \qquad \Gamma^{2nd} = \frac{\gamma_{SL}}{|\Delta S_V|} .$$

**Dissipated work.** $w_p = \int \sigma_y^{s}\,d\alpha_A$ [J m⁻²] is either kept out of the barrier and
reported as $\Delta G_{diss} = 4\pi r^2 w_p$ (`stored`, default) or added to $\gamma_{SL}$ (`barrier`, upper bound).
The matrix is a liquid (no shear), so plasticity is restricted to the surface.

## Main results (Al-3Cu-5Nb-0.1Fe, mean gradient 3294 K/m)

The continuum formulation reproduces the experimental SDAS with both orders of the thermal-field
tensor: MRB with Gamma 1st order (3.887e-7 m K) and 2nd order (1.0685e-6 m K) fall within the
experimental scatter (Ferreira, JMEP manuscript, Table 6). The elasto-plastic surface of this
repository is an extension of that formulation, not a correction of it.

* Elastic branch: `GM_SURFACE_MODEL=elastic` reproduces the original solution exactly.
* Elasto-plastic surface with the plastic area updating the reference (`sigma_y` = 0.0695 N/m, `H` = 0):
  Gamma 2nd order = 3.887e-7 m K, the same at every experimental gradient (1.1–10.8 × 10³ K/m), and
  MRB SDAS = 7.67 t_SL^0.377 µm. `sigma_y` was obtained by matching the SDAS data; it is not an
  independent measurement.
* MD values of the GM surface parameters (Borovikov EAM, pure Al): `sigma_0` = −0.030 ± 0.007 N/m
  (24 000 atoms, 1 ns) and, in this repository, −0.040 ± 0.004 N/m (8192 atoms, 100 ps). Inputs taken
  from MD are an open issue on the MD side (pure Al, no solute, nanometre nuclei), not of the
  continuum theory.
* SDAS written explicitly in the thermal gradient: `SDAS = K (G_L V_L)^(-1/3)`; the free exponent
  against the cooling rate is −0.280 ± 0.059 (95 %), which contains −1/3.

## Running

Python 3.11 with NumPy, SciPy and Matplotlib. From `model/`:

```bash
# elastic GM (reproduces the original script exactly) and elasto-plastic GM
GM_SURFACE_MODEL=elastic python Nucleation_model_..._ElastoPlastic_GM.py
python Nucleation_model_..._ElastoPlastic_GM.py
# elasto-plastic with the plastic area updating the reference
GM_PLASTIC_REFERENCE=1 GM_SIGMA_YA0=0.0695 GM_H_A=0 python Nucleation_model_..._ElastoPlastic_GM.py
# SDAS (reads the GT_point_to_point_*.csv files written above)
python scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026_GT_point_to_point.py
python scriptsdasmodelAl3Cu5Nb01Fe_zoom_horizontal_B_MRS_Meeting_2026_SDAS_vs_GL.py
```

MD (LAMMPS 10Sep2025; on Windows the OpenMP build: `lmp -sf omp -pk omp 8 -in ...`). Download
`Ni-Al-Nb.eam.fs` (Borovikov et al., Int. J. Plast. 178, 104004, 2024) from the NIST Interatomic
Potentials Repository into `potentials/`. `in.planar_stress` starts from a solid–liquid coexistence slab
(`-var coexist <file>`).

## Data sources

Thermal gradients (`G_L = T_R / V_L`) and SDAS: horizontal transient solidification of
Al-3Cu-5Nb-0.1Fe, Mendes et al. (2023). Alloy, diffusion and partition data as in the original
scripts of the author.

## Status

Research code under development. The yield stress of the interface still needs an independent value
(the series of seeds of 1.5–3.0 nm in `md_cpu/` would give the radius dependence of the interface stress).
