# Continuity dossier — elasto-plastic Gurtin–Murdoch surface

State on 2026-10-23. For an agent who takes over the work without the conversation history.
Portuguese version: `HANDOFF.md`.

## 1. Who and what

* **Author:** Prof. Ivaldo Leão Ferreira (UFPA, ileao@ufpa.br). Reply in **Portuguese**; papers in **English**.
* **The author's theory:** the Gibbs–Thomson coefficient is a thermal-field tensor, Γ = A·∇T. Maximizing ΔG gives a
  quadratic in r. The second-order solution is the complete one ("precise"); the first-order one is an approximation.
  The surface energy comes from the Gurtin–Murdoch (GM) surface, with σ₀ < 0 (compressive; the scripts store the
  magnitude |σ₀|).
* **Goal of this work:**
  - generalize the GM surface to elasto-plasticity, within non-equilibrium thermodynamics;
  - apply it to all nucleation models and to the SDAS (MRB) models;
  - write the paper.

## 2. Constraints and preferences (mandatory)

1. **Do not change** anything in `Sim Luane Final/`; work only in `modelos_elasto-plastos/`.
2. **Do not push** to the GPU repository `nucleation-molecular-dynamics`. The public repository
   `github.com/ileaof/gm-elastoplastic-nucleation` (local: `C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation`,
   commit `ddb3d00`, 2026-10-24) **was synced** at the author's request (it will be cited in the paper): published
   `nucleacao/` (the `stored_total` copies, `gm_elastoplastic.py`, figures, tables), `sdas/`, `ferramentas/` and
   `paper/` (including `nucleation_elastoplastic_AlSiMgFe.py` and its outputs), with both READMEs updated to
   `stored_total`. Deliberately not in the repo: the `.docx` manuscript, `HANDOFF.md` (internal), `Sim Luane Final`
   and the originals that exist only in the author's workspace.
3. **Physical values:** ask before assuming. Describe facts and ask; do not present the author's code as an "error".
4. **SDAS:**
   - the author says no label is swapped (1st order approximate, 2nd order precise);
   - the most precise curve is the one with the largest Γ (2nd order);
   - keep the original curves and legends.
   The only change in the copies is the legend *handles* (see §5).
5. **Nucleation rate:** keep exp(+ΔG*/ΔG*_ref) as in the original.
6. **Physics adopted by the author:**
   - the plastic deformation of the surface is energy stored for phase nucleation;
   - the elasto-plastic Γ must be very close to the elastic one and grow in the same way;
   - the plastically deformed Γ is larger than the equilibrium one, which allows the deformation-induced nucleation
     of another phase;
   - only the elastic part is restored;
   - the phase change is governed by the rate of density change;
   - ΔT_SL must be the mean integral (as the "Mean Integral" values of the `Thermal_data_*` scripts);
   - solidification experiments have high gradients: figures must cover the whole range of ∇T (up to ≥ 10⁵ K/m).
7. **Environment:**
   - Windows, Python 3.11, LAMMPS without MPI (use OpenMP);
   - **never** `taskkill //IM python.exe` (it kills everything);
   - in bash heredocs, the `\n` of `\nabla` becomes a line break: edit LaTeX with the Edit tool;
   - Word files open in Word are locked: save under another name.

## 3. Current model (decisions in order)

Core file: `nucleacao/gm_elastoplastic.py`, class `GMElastoPlasticSurface`, `step(r)` → (τ, s, γ, ∂γ/∂r).

* **Areal yield:**
  - criterion |s| ≤ σ_y + Hα, with s = |σ₀|E_A^e;
  - return mapping and consistent tangent C·H/(C + H);
  - H = 0 (no hardening at the liquidus).
* **Adopted mode `stored_total`** (author's choice, option "c"):
  - formula: γ = γ₀(r_eq/r)² + |σ₀|(|E^e| + |E^p|) − q_diss = γ_el − σ_y·α;
  - fixed reference (r_eq does not change);
  - effect: Γ is 0.5–4.2 % below the elastic one.
* **Rejected modes** (kept for comparison):
  - `complement` (updated reference, Γ down to −29 %);
  - `barrier`, `stored` (the rate decreased monotonically — "physically strange" for the author);
  - `restored*`.
* **Elastic restoration** (`ferramentas/ajustar_restauracao.py`):
  - ε̇_v = β_s·Ṫ/ΔT_SL;
  - f = ε̇_vτ_r/(1 + ε̇_vτ_r);
  - Γ = Γ_EP − f·|s|/|ΔS_V|, with |s| ≤ σ_y.
  τ_r was not defined; f = 1 (bound) is used. The old version (restoring everything down to γ₀, τ_r = 7.84e3 s fitted
  to the SDAS, error 3.6 %) was **rejected** by the author because it brought Γ to equilibrium.
* **σ_y** (`ferramentas/sigma_y_ligas.txt`): 0.0386 / 0.1572 / 0.0423 N/m. They come from an old fit in `barrier`
  mode and are not a fit to the SDAS in the current mode.
* **Sweep:** `GM_DT_MAX` defaults to 30 K, so every copy covers the whole range up to ~10⁶ K/m.

## 4. Current numbers (`stored_total` mode)

| Alloy | yield | Γ²ⁿᵈ el./EP (mean ∇T) | SDAS 2nd el. / EP | γ at 10⁵ K/m el./EP | Γ²ⁿᵈ at 10⁵ K/m el./EP |
|---|---:|---|---|---|---|
| Al-3Cu-5Nb-0.1Fe (3294 K/m) | 731 K/m | 1.0685e-6 / 1.0627e-6 | 9.38 / 9.37·t^0.379 | 5.58 / 5.53 J/m² | 5.693e-6 / 5.644e-6 |
| Al-0.8Si-0.6Mg-0.2Fe (7902) | 925 K/m | 1.5745e-6 / 1.5083e-6 | 6.71 / 6.66·t^0.387 | 4.45 / 4.24 | 4.805e-6 / 4.581e-6 |
| Al-3Cu-0.5Mg-0.2Fe (3708) | 562 K/m | 1.0973e-6 / 1.0878e-6 | — / 6.93·t^0.384 | 4.35 / 4.29 | 5.210e-6 / 5.142e-6 |

**Point the author may come back to:** he stated that "with the elastic one 10⁵ K/m is not possible". The facts are:
* in the copies (closed-form ΔT), elastic and EP both reach ~1.1×10⁶ K/m;
* what does not reach it is the original script, because `brentq(f_DT, 0.0001, 1)` limits ΔT ≤ 1 K and the number of
  points is fixed.
This was told to him and written as a fact in the paper; he has not replied yet.

## 5. What was done in this session (files)

* **Nucleation copies:**
  - `nucleacao/*_EP.py` are generated by `ferramentas/aplicar_elastoplastico.py`, with `ENERGY_MODE = 'stored_total'`;
  - `nucleacao/figure_layout.py` moves texts and legends inside the axes in all figures of the scripts, and never
    blocks the display (a figure that cannot be adjusted is left as it is, with a warning).
* **Figures generated by every nucleation copy** (the original code of the rate against ΔT is kept):
  - rate against ΔT (figure 13): the elastic rate of the same script in blue, computed by
    `figure_layout.elastic_rate()`, which runs the script in elastic mode in a separate process with `GM_DUMP_RATE`
    (the script saves the arrays and stops before the figure), plus the elastic range and the yield;
  - new figures against ∇T over the whole range (log scale, up to ~10⁶ K/m): nucleation rate (block
    `RATE_VS_GRADIENT`), second-order Γ, surface energy γ and surface tension τ (`figure_layout.vs_gradient_figure`),
    each with the elastic curve, the yield (shaded elastic range), the mean experimental ∇T and the values at 10⁵ K/m;
  - only the alloy name of the copy stays in the rate figure (`ALLOY_TEXT`). The Al-Si-Mg-Fe original also contains
    "Al-6wt%Cu-4wt%Si-0.11wt%Fe", which was the only one visible in the original.
* **SDAS copies:**
  - `sdas/*_EP.py` are generated by `ferramentas/aplicar_sdas.py`;
  - `legend_handles()` ties each legend text to its own line (in the original, the 2nd order is drawn before the 1st
    and `plt.legend` with texts only swapped the line styles);
  - the grey curve (Al-Cu-Nb-Fe) is the constant Γ at the mean gradient minus the restored elastic part (f = 1):
    9.30·t^0.379. The version with the local point-to-point Γ was judged wrong by the author.
* **Composed surface tension** (author, 2026-10-23): τ = (σ₀ + s^e) + |σ₀|E_A^p + σ_y = τ_el + σ_y (H = 0), in
  `gm_elastoplastic.py` (`hist`, `tau_total`); figure `tensao_superficial_alto_gradiente.py`.
* **Comparison figures in `nucleacao/figuras/`:** `comparar_curvas_EP.py`, `comparar_taxa_nucleacao.py`,
  `energia_superficie_alto_gradiente.py` and `gibbs_thomson_alto_gradiente.py`, all with the yield marked.
* **Fix:** `rodar_modelos_EP.py` now uses the newest table of each folder (it used to read CSVs of old modes).
* **Open report:** the author said the nucleation copies "no longer generate figures". The test here (qtagg) created
  all figures without error. It is still unknown how he runs them (IDE or backend) and which message appears.
* **Paper:**
  - file: `C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v6_final.docx` (final version for publication, identical copy of v6: the four formulations of the presentation + spherical solutions of ΔG_S and its potential; v4 = general notation of slides 12–16 in elasto-plastic form, v5 = + derivation from the metric with area variation);
  - 19 pages, 28 native OMML equations, 6 figures, 7 tables, in English; includes the composed surface tension
    (Eq. 22) and Section 5.4 with Γ, γ and τ of the three alloys over the whole ∇T range (up to 10⁶ K/m) and Table 6 at 10⁵ K/m;
  - v1 to v4 are **outdated**;
  - sources in `paper/`: `build_paper.py` (python-docx + latex2mathml + Office MML2OMML.XSL) and `make_figures.py`;
  - verification render: Word COM → PDF → PyMuPDF;
  - **ScholarOne rejected v6_final** ("Caught XML parse error / Validation error"): the plain text carried
    105 non-breaking hyphens (U+2011, in the alloy designations) and ~200 Unicode superscripts (U+207x:
    the "10^-6 K/m"). `sanitize_text()` in `build_paper.py` (commit `6339354` in the repo) converts U+2011
    to "-" and Unicode superscripts to real superscript runs (ASCII digits); a verification rebuild came
    out with 0 offending characters and the 303 equations untouched. For resubmission:
    `Manuscript_..._v6_final_clean.docx` (sanitized v6_final, same content). Figures 26 and 27 of the
    tweets come from `paper/nucleation_elastoplastic_AlSiMgFe.py` (window 26 = 3x3 curves panel, 27 = 2x3
    energies/entropies panel).
* **README.md / README.en.md** updated for `stored_total`.

## 6. Pending items and suggestions

0. **Author decision (2026-10-23):** publish with γ = γ₀(r_eq/r)² − s (script form) — the form of all results of paper
   v5. The spherical metric form (Eq. 89 of the Materials paper, with the area integral) and the reading of Σ in the
   integral from A₀ to A₁ are left for later (tests in `ferramentas/formas_energia_gt_sdas.py` and
   `sdas_metrica_area_EP.py`). Also open: which Γ the author considers 2nd order in the SDAS script (the comments he
   pasted call `z3n_1st_Eq`, with Γ = 3.887e-7, 2nd order; the script and the nucleation model give 1.0685e-6 as 2nd).

1. **SDAS:** the 2nd-order EP curve (9.37) and the restored one (9.30) coincide with the elastic 2nd order (9.38),
   which lies above the points. The 1st order (7.67) "fits", but the author considers it only an approximation. He
   said he will look at the SDAS later. Do not fit τ_r or invent laws without asking.
2. **Paper v2:**
   - review with the author;
   - the Mendes et al. (2023) reference is incomplete (placeholder);
   - the SDAS is not in the paper;
   - figures: Fig. 1 = `curvas_EP_Al3Cu5Nb01Fe.png`; Figs. 3–5 = `gibbs_thomson_1e5.png`, `energia_superficie_1e5.png`,
     `tensao_superficial_1e5.png` (from `nucleacao/figuras/`).
3. **GitHub repository:** update with the `stored_total` mode and the paper only if the author asks.
4. **Wrong text in the original:** figure 16 of Al-Cu-Nb-Fe shows the name of another alloy (from the original). The
   author was asked and did not reply.
5. **H₂O:** provisional σ_y (0.0695 N/m).

## 7. How to reproduce everything

See "How to run" in `README.en.md`. Order:
1. `aplicar_elastoplastico.py`;
2. `rodar_modelos_EP.py`;
3. `ajustar_restauracao.py`;
4. `aplicar_sdas.py`;
5. `rodar_sdas_EP.py`;
6. figure scripts.

The paper: `python -I paper/build_paper.py`, after `paper/make_figures.py` and after copying
`curvas_EP_Al3Cu5Nb01Fe.png` to `paper/Fig1_curves.png`. Output paths are at the top of the scripts.

## 8. E6 (MD recrystallization) — what to run on the GPU machine and what to send back

Context: the professor wants to prove that the stored plastic deformation (`stored_total`, sigma_y*alpha)
reappears as a potential in recrystallization. This laptop (i7-1360P, no NVIDIA GPU) already completed the
pilots and the small series (2 048 atoms: 500-700 K Jelinek alloy, 500-900 K Mendelev control) — everything
is in the repo, commit `e98d66b`, folder `recristalizacao/` (PLANO.md has the numbers and the night log).
The remaining runs are heavy and belong on the **NVIDIA machine**. Full procedure: `recristalizacao/GPU.md`
in the repo.

**On the GPU machine** (LAMMPS; MEAM runs on MPI+OpenMP, EAM takes `-sf gpu`):
1. `git clone https://github.com/ileaof/gm-elastoplastic-nucleation.git`
   `cd gm-elastoplastic-nucleation/recristalizacao`
2. Set `LMP` / `NP` (MPI ranks) / `NT` (threads) at the top of `rodar_gpu.sh`.
3. `bash rodar_gpu.sh` — G0 validation (~1 min); G1 big cell 16 384 atoms (~1-3 h);
   G2 series 500-900 K on the big cell (~10-16 h); G3 fine frontier 750/775/800 K x 3 thermal seeds
   (~4-8 h); G4 long anneal 1 ns at 750 K (~10-20 h). Per section: `SECTIONS="G1 G3" bash rodar_gpu.sh`.
4. `bash juntar_resultados.sh` -> creates `results_e6.tar.gz` (logs, stored/perfil .dat, figures; the
   dumps stay there — they are big).

**What to send back to THIS machine**
(`C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\modelos_elasto-plastos\recristalizacao\`):
* **Required:** the `results_e6.tar.gz` (OneDrive/transfer/USB), unzipped into the folder above. If the GPU
  machine has the GitHub credentials, `git push` (the `results_e6/` folder) and pull here also works.
* **Optional (for OVITO):** the `recozido_*.dump` of interest and the large `deformed_jelinek.dump`.
* **Not needed:** the `data_*.data` (G1 recreates its own) nor the large raw logs.

After it comes back: run `python analyze_recristalizacao.py` in the folder and report to the professor
(E_st(eps), delta, the statistical recrystallization threshold, comparison with sigma_y*alpha/delta).
Questions for the professor on return: (a) new design for 800-900 K (finer dt, barostat, bigger cell);
(b) thermal term (Gamma) in the anneal; (c) delta measured or imposed; (d) the "recrystallized" criterion.

Agent's persistent memory:
`C:\Users\ileao\.claude\projects\c--MacbookPro-ufpa-PPGEM-SDAS-Al-Si-Mg-Fe\memory\` (index `MEMORY.md`).
