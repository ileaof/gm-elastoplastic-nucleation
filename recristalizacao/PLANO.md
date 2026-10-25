# E6 — Recristalização por MD com a deformação plástica armazenada como potencial

Proposta de 2026-10-24/25 (autor: Prof. Ivaldo Leão Ferreira). O material tem deformação plástica
armazenada (modo `stored_total`, energia σ_y·α por área na interface); elevando a temperatura até ativar
a difusão atômica, a recristalização deve iniciar, liberando essa energia como força motriz. O MD deve
provar a existência, medir a magnitude e comparar com o potencial previsto pelas mesmas equações.

## 1. A física a provar (a mesma estrutura do paper)

No modo `stored_total` a interface carrega energia plástica por área, ω_st = σ_y·α (H = 0; o sweep dá
γ_el − γ_EP por ∇T). A recristalização é a nucleação do grão isento de deformação:

	ΔG_rec(r) = 4π r² γ_gb f(θ)/4 − (4/3)π r³ · (ω_st/δ)          →  r* ≈ 2 γ_gb / (ω_st/δ)

a mesma quadrática do paper, com o papel do super-resfriamento (ΔS_V·ΔT) trocado pela energia
armazenada por volume, ω_st/δ — δ = espessura da zona deformada (a converter J m⁻² → J m⁻³).
Fecha coerente com as regras do professor: a parte elástica é restaurável por si (`restoration`,
f = 1); a permanente só sai consumida pela recristalização (α → 0 devolve a interface ao equilíbrio).

## 2. Números contínuos previstos (do sweep, fonte: HANDOFF.md §4)

| Liga | ω_st = γ_el − γ_EP em 10⁵ K/m [J m⁻²] | ω_st/δ com δ = 5 nm [J m⁻³] | escoamento ∇T [K m⁻¹] |
|---|---:|---:|---:|
| Al-3Cu-5Nb-0,1Fe | 0,049 | ~1×10⁷ | 731 |
| Al-0,8Si-0,6Mg-0,2Fe | 0,208 | ~4×10⁷ | 925 |
| Al-3Cu-0,5Mg-0,2Fe | 0,060 | ~1×10⁷ | 562 |

Faixa de referência de metais deformados: 10⁵–10⁷ J m⁻³. O potencial previsto (com δ de poucos nm)
está no topo dela — mensurável em MD.

## 3. Capacidades desta máquina (verificado 2026-10-24)

* **CPU:** Intel i7-1360P, 12 núcleos / 16 threads (4P+8E), 2,2 GHz base.
* **RAM:** 32 GB (3,8 GB livres no momento; MD de ≤10⁵ átomos usa ≤1 GB).
* **GPU:** Intel Iris Xe (2 GB, integrada) — **não há NVIDIA**; rota LAMMPS **OpenMP** (o estudo GPU
  continua na outra máquina, `nucleation-molecular-dynamics` — sem push, regra do dossier).
* **Disco:** 157 GB livres em C.
* **LAMMPS:** 10Sep2025 (`C:\Users\ileao\AppData\Local\LAMMPS 64-bit 10Sep2025\bin\lmp`), pacote OPENMP
  presente (GPU/KOKKOS sem uso aqui). Convenções igual ao `e2/e5` do `Nucleation_MD` (units metal,
  dt = 0,002 ps, pot_*.lmp, comentários no topo).
* **Veredicto:** suficiente para prova de princípio — sistemas de até ~50×10³ átomos; passo 1 (medir)
  em minutos; recozimentos de 0,1–1 ns em horas por temperatura. Série pesada fica na outra máquina.

**Cronometragem medida aqui (bench real, 1000 passos NVT, 2 048 átomos; logs `bench_*.log`):**

| Configuração | Tempo | Custo |
|---|---:|---|
| Jelinek MEAM liga, 8 threads | 51,2 s | 25 μs/átomo/passo |
| Jelinek MEAM liga, 16 threads | 45,3 s | 22 μs/átomo/passo (só 12 % melhor: satura) |
| Mendelev EAM Al puro, 8 threads | 4,4 s | 2,1 μs/átomo/passo (12× mais rápido) |

**Tempos projetados das corridas do plano (Jelinek, 16 threads):**

| Corrida | Passos | Tempo |
|---|---:|---|
| E6.1 piloto (2 048 átomos, 85 000 passos) | 85 k | **~1 h** |
| E6.1 com Mendelev (idem) | 85 k | **~6 min** |
| E6.2 recozimento, 1 temperatura (60 000 passos, 2 048 átomos) | 60 k | **~45 min** |
| E6.2 série das 6 temperaturas (500–900 K) | 360 k | **~4,5 h** (uma noite) |
| E6.1 escala média (nx=10, 4 000 átomos) | 85 k | ~2 h |
| E6.1 grande (nx=16, 16 384 átomos) | 85 k | ~8 h → série larga ir para a máquina GPU |

Recomendação: `-pk omp 8` (16 threads quase não ganha); usar o Mendelev para validar a pipeline e o
Jelinek para os números da liga.

## 4. Potenciais (já no repositório MD, copiados para cá)

| Potencial | Arquivo | Papel |
|---|---|---|
| **Jelinek 2012, MEAM Al-Si-Mg-Cu-Fe** | `AlSiMgCuFe_Jelinek2012/` | **proposta principal**: simula a liga real (Al-0,8Si-0,6Mg-0,2Fe, sem o Cu). T_m já cravada no e2: **936 ± 3 K** (DTA 936,58) |
| Mendelev 2008, EAM Al puro | `Al_Mendelev2008/Al_mm.eam.fs` | alternativa de controle (mecanismo; T_m 932,3 K) |

Composição da liga proposta (at% de 99,4% Al da liga): Al 98,46 / Si 0,77 / Mg 0,67 / Fe 0,096 %.

## 5. As três corridas propostas

1. **E6.1 deformar e medir** (`in.deformar.lmp`): FCC com solutos, equilíbrio a 300 K, cisalhamento
   simples em degraus (ε = 0,0–0,30), relaxação a cada degrau; grava `stored.dat` com
   E(ε) por átomo e escreve `data_<ε>.data` para o recozimento. Minutos.
2. **E6.2 recozimento** (`in.recozer.lmp`): degraus de T = 500–900 K (T_m Jelinek 936 K), NPT;
   termo (pe, temp, press) + CNA (r_c = 0,854·a₀) e dumps para o OVITO. O início da recristalização
   = queda de PE acima do relaxamento (recuperação) acompanhada de FCC novo.
3. **E6.3 comparação** (python): T_onset vs ω_st^MD; E_st(ε) vs σ_y·α (o potencial previsto);
   δ medida pela espessura da zona rica em energia no perfil do degrau.

## 6. Decisões do professor antes de rodar

1. Jelinek (Al-Si-Mg-Fe, proposto) ou Al puro de Mendelev como primeiro alvo?
2. Recozimento só com o potencial plástico, ou também com termo térmico (Γ)?
3. δ: medida no MD dos perfis (proposto) ou valor imposto?
4. Critério de "recristalizado": queda de PE + fração FCC/CNA (proposto) ou orientação nova (PTM/OVITO)?
5. Escala da primeira série: piloto 2 048–16 384 átomos (esta máquina) ou já a série grande?

## 7. Como rodar (depois do ok)

	cd "C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\modelos_elasto-plastos\recristalizacao"
	lmp -sf omp -pk omp 16 -in in.deformar.lmp -var pot jelinek -var nx 8 -var eps_list 0.05 0.10 0.14 0.20 0.30
	lmp -sf omp -pk omp 16 -in in.recozer.lmp   -var pot jelinek -var data data_0.14.data -var T 800
## 8. Registro da noite (25/10/2026)

| Hora | Evento |
|---|---|
| 23:15-23:20 | 4 erros de sintaxe/protocolo pegos pela validação (ave/chunk, triclinic, unidades reduzidas, erate) |
| 23:36 | Piloto Mendelev ✓: pico de E_st exatamente em ε=0,14; ~75% liberados após o escoamento |
| 00:26 | Piloto Jelinek ✓ (mesma assinatura); evento anomal no relax de ε=0,30 (T→0) — re-corrida pendente |
| 01:51 | T700: lost atoms (2 ejetados) → protocolo endurecido (minimize + neigh_modify every 1) |
| 03:02 | T500/T600/T700 ✓ com protocolo endurecido (uniformizada a série) |
| 03:44 | T800: ejeção de novo — conclusão: a fronteira 700-800 K É o início violento |
| 03:50 | Lançada a noite 3: 800/850/900 K com subira blindada (nve/limit + berendsen) + série Mendelev |
| 03:58 | T800 ejetou mesmo blindado (oscilou 40-1380 K na subida): regime 800-900 K e intrinsecamente catárstico nesta célula — redesenho (dt, barostato, célula maior) fica como decisão do professor; série Jelinek válida: 500-700 K |
| 04:03 | Lançada noite 4: só a série Mendelev (500-900 K) para completar a matriz |
