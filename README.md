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

| Surface | Gamma_het 2nd order [m K] | SDAS (MRB, prefactor 5.5) | MAPE vs experiment |
|---|---|---|---|
| Elastic GM (original, `sigma_0` = 1.3761 N/m) | 1.0685e-6 | 9.38 t_SL^0.379 µm | 23.8 % |
| Elasto-plastic, fixed reference (`sigma_y` = 0.5 N/m) | 1.0146e-6 | 9.28 t_SL^0.379 µm | 22.5 % |
| Elasto-plastic, reference updated (`sigma_y` = 0.0695 N/m, `H` = 0) | 3.887e-7 | 7.67 t_SL^0.377 µm | 4.6 % |

* With the updated reference the 2nd-order tensor is the same at every experimental gradient
  (1.1–10.8 × 10³ K/m), so point-to-point and mean-gradient calculations coincide.
* `sigma_y` = 0.0695 N/m was obtained by matching the SDAS data; it is not an independent measurement.
* MD values of the GM parameters (Borovikov EAM, pure Al) give `sigma_0` = −0.030 ± 0.007 N/m
  (24 000 atoms, 1 ns) and, in this repository, −0.040 ± 0.004 N/m (8192 atoms, 100 ps). With the MD
  `sigma_0` the elastic GM surface already gives a 2nd-order tensor of 2.97e-7 m K, consistent with the
  data: the excess surface stress of the nucleus is of order 0.03–0.08 N/m, not ~1 N/m.
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
