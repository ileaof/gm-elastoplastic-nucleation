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

## Elasto-plastic surface model

Small-strain additive split in the tangent plane, `eps_s = eps_s^e + eps_s^p`, with yield acting on
the excess stress `s~ = sigma_s - tau0 P`, so that the residual stress alone never yields:

* areal (dilatational) mechanism, `|tr s~| <= 2 sigma_y^s(alpha)`: loss of coherency of the interface;
* J2 on the tangent-plane deviator: inactive for a spherical nucleus (equibiaxial strain), checked at run time.

`gamma_0`, `sigma_0` and `r_eq` are equilibrium values (undeformed created area). With the option
`GM_PLASTIC_REFERENCE=1` the plastic area moves the reference, `gamma = gamma_0 (r_eq,p / r)^2 - s` with
`r_eq,p = r_eq exp(E_A^p / 2)`. For perfect plasticity (`H = 0`) the surface energy after yield is
`gamma = gamma_0 exp(sigma_y / C_A) + sigma_y`, independent of the radius and of the thermal gradient.

Sign convention: `sigma_0 < 0` (compressive force). In the scripts `sigma_0` is its magnitude, since
`tau = -sigma_0 - 2 sigma_0 ln(y0/y)` gives `tau(r_eq) = -sigma_0`.

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
