# E5 (CPU) — tensão de interface: placa plana e piloto de núcleo esférico (Laplace)

Máquina: i7-1360P, LAMMPS 10Sep2025 do Windows com OpenMP (`lmp -sf omp -pk omp 8`; o binário do Windows não
usa MPI). Desempenho: ~1,0 h/ns para 8192 átomos (Borovikov EAM). Data: 2026-10-20/22.

Objetivo: obter o excesso de tensão de Gurtin–Murdoch do núcleo, s = f(r) − τ0, para testar o modelo
elasto-plástico de superfície (escoamento de área σ_y^s) do script `..._ElastoPlastic_GM.py`.

## Placa plana (`in.planar_stress`, `analyze_planar.py`)
Coexistência da E2 (`../e2_melting/coexist_npt_borovikov_930.data`), NPH z, 100 ps, perfil σ_αα(z) por átomo,
excesso sobre as fases bulk com divisória equimolar.

**f_plana = τ0 = −0,040 ± 0,004 N/m** (blocos de 20 ps; sistemático maior). Concorda com a E5 da GPU
(`Nucleation_MD_GPU/md/e5_surface_stress/`: σ0 = −0,030 ± 0,007 N/m, 24 000 át., 1 ns).

## Piloto de núcleo (`in.seed_laplace`, `analyze_laplace.py`)
Protocolo do JMEP Sec. 3.4, r0 = 2,0 nm, T* = 809,3 K (Tabela 4), caixa 24³ (55 296 át.); f = ΔP·r/2.

| | valor |
|---|---|
| a(T*) | 4,108 Å |
| ΔP (janela 42–68 ps) | 349 ± 102 bar |
| f (r = r0) | +0,035 ± 0,010 N/m (≈ 0,039–0,041 com o raio pela contagem de sólidos) |
| γ = Γ·\|ΔS_V\| (Tabela 4) | 0,108 J/m² |
| s = f − τ0 | +0,075 ± 0,011 N/m |

Ressalvas: o núcleo cresceu na produção (N_sólido 2538 → 3379; caixa menor que a do paper) e há um só seed.

## Decisões do autor
- **σ0 < 0 (força compressiva).** Nos scripts `sigma_0` é o módulo (τ(r_eq) = −sigma_0). σ0 do MD entra como +0,030.
- Série de seeds (1,5–3,0 nm, caixa 32³, ~6–7 h CPU): não lançada; só é necessária para testar a tendência τ(r).
