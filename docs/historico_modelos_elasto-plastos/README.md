# Modelos elasto-plásticos

Cópias dos modelos de nucleação (superfície de Gurtin–Murdoch elástica) com a superfície de
Gurtin–Murdoch **elasto-plástica**, e cópias dos modelos de SDAS (MRB) com o tensor de Gibbs–Thomson
elástico de 1ª e 2ª ordem (como no original) e o elasto-plástico de 2ª ordem.
**Os arquivos originais (`Sim Luane Final/`) não foram alterados**: as ferramentas só os leem.

Formulação completa: paper `paper/` (Word em `OneDrive/Documentos/Paper Nucleation/
Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v2.docx`). Para retomar o trabalho, ver `HANDOFF.md`.
Versão em inglês: `README.en.md` (dossier: `HANDOFF.en.md`).

## Organização

| Pasta | Conteúdo |
|---|---|
| `nucleacao/` | 11 modelos de nucleação `*_EP.py`, o módulo comum `gm_elastoplastic.py`, `figure_layout.py` (textos e legendas dentro dos eixos) e os auxiliares (`derivative_shim.py`, `midpoint.py`, `mpl_textguard.py`) |
| `nucleacao/figuras/` | comparações elástico x elasto-plástico contra ∇T (`curvas_EP_*`, `taxa_nucleacao_*`, `energia_superficie_1e5.png`) |
| `nucleacao/resultados/` | uma pasta por modelo, com as tabelas das rodadas (a mais recente vale) |
| `nucleacao/nao_convertidos/` | `Nucleation_model_..._proposal_luana_GL_data_indev.py`, copiado sem conversão (ver abaixo) |
| `sdas/` | 5 modelos de SDAS `*_EP.py` e `figuras/` |
| `ferramentas/` | scripts que fizeram as cópias, as rodadas, as figuras e as verificações |
| `paper/` | fontes do paper (`build_paper.py`, `make_figures.py`, figuras) |
| `estudo_AlCuNbFe/` | estudo inicial em Al-3Cu-5Nb-0,1Fe (modos antigos, Γ ponto a ponto, SDAS x ∇T); histórico |

## Superfície elasto-plástica (`nucleacao/gm_elastoplastic.py`)

O módulo faz o papel de `sigma_func`, `gamma_func_r` e `dgammadr_func_r` no laço principal: escoamento de
área sobre a tensão em excesso s = |σ₀|·E_A^e, |s| ≤ σ_y + Hα (return mapping, r diminuindo a partir de
r_eq). O escoamento marca **onde termina o regime elástico e começa o plástico** (560–925 K/m nas três ligas).

**Energia de superfície adotada (modo `stored_total`, padrão das cópias, decisão do autor em 2026-10-22):**
a deformação permanente complementa a elástica, a sua energia fica armazenada na interface e só o calor
dissipado pelo fluxo plástico sai; a referência de equilíbrio não muda (r_eq fixo):

    γ_SL = γ₀·(r_eq/r)² + |σ₀|·(|E_A^e| + |E_A^p|) − q_diss = γ_SL,elástico − σ_y·α        (H = 0)

Por isso o Γ elasto-plástico acompanha o elástico, cresce do mesmo jeito com ∇T e fica sempre acima do
equilíbrio γ₀/|ΔS_V| (energia armazenada disponível para nucleação induzida por deformação).
Derivada: ∂γ/∂r = −2γ₀r_eq²/r³ − (|σ₀| − σ_y)·∂E_A/∂r depois do escoamento (H = 0).

Outros modos (comparação): `complement` (referência atualizada r_eq,p = r_eq·exp(E_A^p/2); Γ até 29 % abaixo
em 8×10⁴ K/m), `barrier` (+w_p), `stored` (só o elástico), `thermo`, `restored`, `restored_barrier`.

Configuração (variáveis de ambiente): `GM_SURFACE_MODEL` (`elastic` | `elastoplastic`), `GM_SIGMA_YA0`,
`GM_H_A` (padrão 0), `GM_DISSIPATION_MODE` (padrão `stored_total`), `GM_PLASTIC_REFERENCE` (padrão 0 em
`stored_total`, 1 nos outros), `GM_DT_MAX` (padrão 30 K: a varredura cobre toda a faixa, até ~10⁶ K/m).

Outras mudanças nas cópias:

* ΔT em forma fechada, ΔT = 2γ/(r|ΔS_V|) (a raiz de `f_DT` dos originais); ∇T = ΔT/(8πr). O
  `brentq(f_DT, 0.0001, 1)` dos originais só encontra ΔT ≤ 1 K: é por isso que o script original não chega
  a 10⁵ K/m (ΔT ≈ 5 K). Nas cópias, elástico e elasto-plástico vão até ~1,1×10⁶ K/m (fim: raiz do ângulo de
  2ª ordem);
* varredura até ΔT = `GM_DT_MAX`, vetores com 400 posições; θ de 1ª ordem NaN quando não tem raiz;
* índice `j` do gradiente experimental achado na varredura; CSVs com o modelo de superfície no nome;
* `*_dev_gam`: forma métrica do original, 3γ₀/(1 + 2r²/r_eq²);
* figura da taxa contra ΔT (mantida) com eixos automáticos, a curva elástica do mesmo script e o escoamento;
  todas as figuras: textos e legendas que saíam dos eixos vão para dentro (`figure_layout.py`). A expressão da
  taxa, exp(+ΔG*/ΔG*_ref), não mudou;
* **figuras novas contra ∇T em toda a faixa (log, até ~10⁶ K/m)**, geradas pelo próprio script: taxa de
  nucleação, Γ de 2ª ordem, energia de superfície γ e tensão superficial τ, cada uma com a curva elástica do
  mesmo script (`figure_layout.elastic_rate`, processo separado), o escoamento, o ∇T experimental médio e os
  valores em 10⁵ K/m (`figure_layout.vs_gradient_figure`).

### Resultados (modo `stored_total`)

| Liga | escoamento ∇T | Γ²ⁿᵈ elástico / EP no ∇T médio | Γ²ⁿᵈ elástico / EP em 10⁵ K/m | γ elástico / EP em 10⁵ K/m |
|---|---:|---|---|---|
| Al-3Cu-5Nb-0,1Fe (3294 K/m) | 731 K/m | 1,0685e-6 / 1,0627e-6 (−0,5 %) | 5,693e-6 / 5,644e-6 | 5,58 / 5,53 J/m² |
| Al-0,8Si-0,6Mg-0,2Fe (7902 K/m) | 925 K/m | 1,5745e-6 / 1,5083e-6 (−4,2 %) | 4,805e-6 / 4,581e-6 | 4,45 / 4,24 J/m² |
| Al-3Cu-0,5Mg-0,2Fe (3708 K/m) | 562 K/m | 1,0973e-6 / 1,0878e-6 (−0,9 %) | 5,210e-6 / 5,142e-6 | 4,35 / 4,29 J/m² |

Pico da taxa de nucleação (elástico → EP): 6,62e49 → 6,37e49 (Al-Cu-Nb-Fe), 1,01e50 → 8,19e49
(Al-Si-Mg-Fe), 1,43e50 → 1,33e50 (Al-Cu-Mg-Fe). Figuras: `comparar_curvas_EP.py`,
`comparar_taxa_nucleacao.py`, `energia_superficie_alto_gradiente.py` (todas com o escoamento marcado:
faixa sombreada = regime elástico, linha traço-ponto = escoamento). Com H = 0, depois do escoamento
τ = σ₀ − σ_y fica constante.

Γ²ⁿᵈ de todas as cópias: `ferramentas/gamma_2nd_EP_ligas.txt` (`rodar_modelos_EP.py`; dev_gam 1,0611e-6;
H₂O 9,8668e-7).

### Verificação

* `GM_SURFACE_MODEL=elastic`: cada cópia reproduz o original **exatamente** (`verificar_copias_elasticas.py`
  → `verificacao_elastica.txt`);
* derivada analítica de `stored_total` x diferenças finitas: 3×10⁻⁶; sem q_diss, igual ao elástico a 2,5×10⁻¹⁶.

## σ_y por liga

`ferramentas/sigma_y_ligas.txt`: 0,0386 (Al-Cu-Nb-Fe), 0,1572 (Al-Si-Mg-Fe), 0,0423 (Al-Cu-Mg-Fe),
0,0695 N/m (H₂O, provisório). Vêm de um ajuste antigo (modo `barrier`, `ajustar_sigma_y.py`); no modo
`stored_total` não são ajuste ao SDAS. Como q_diss = σ_y·α some para σ_y → 0 e σ_y → ∞, o desvio do Γ
elástico é máximo em σ_y intermediário.

## Cópias de SDAS (`ferramentas/aplicar_sdas.py`)

Cada cópia mantém as curvas e os rótulos do original e acrescenta o MRB com o Γ_het de 2ª ordem
elasto-plástico (pontilhada grossa). A legenda recebe as linhas explicitamente (handles), para que cada
texto fique ao lado da linha da sua curva: tracejada = 1ª ordem, traço-ponto = 2ª ordem (no original o
plt.legend só com textos associava as linhas pela ordem de desenho, 2ª ordem desenhada antes da 1ª).
Figuras em `sdas/figuras/` (`rodar_sdas_EP.py`).

| Cópia de SDAS | Γ²ⁿᵈ EP [m K] | SDAS (MRB) EP | SDAS (MRB) 2ª ordem elástico |
|---|---:|---|---|
| `Al3Cu5Nb01Fe_..._B_MRS_Meeting_2026_EP` | 1,0627e-6 | 9,37·t_SL^0,379 µm | 9,38·t_SL^0,379 µm |
| `Al08Si06Mg02Fe_..._new_script_EP` e `..._Ferreira_2024_EP` | 1,5083e-6 | 6,66·t_SL^0,387 µm | 6,71·t_SL^0,387 µm |
| `Al3Cu05Mg02Fe_..._Ferreira_2025_EP` | 1,0878e-6 | 6,93·t_SL^0,384 µm | — |
| `Al3Cu05Mg02Fe_..._Ferreira_2025_dev_gam_EP` | 1,0611e-6 | 6,90·t_SL^0,384 µm | — |

O SDAS elasto-plástico de 2ª ordem coincide com o elástico de 2ª ordem.

### Restauração elástica pela taxa de variação de densidade (Al-3Cu-5Nb-0,1Fe)

Só a parte elástica (recuperável) é restaurada; a energia da deformação permanente fica
(`ferramentas/ajustar_restauracao.py`):

    ε̇_v = β_s·Ṫ/ΔT_SL,   β_s = (ρ_S − ρ_L)/ρ_S = 0,0767,   Ṫ = ∇T_L(P)·V_L(P)
    ΔT_SL = 1/(90 − 5)·∫₅⁹⁰ Ṫ(P)·t_SL(P) dP = 185,3 K   (médio integral na faixa medida)
    f = ε̇_v·τ_r/(1 + ε̇_v·τ_r),   Γ(P) = Γ_EP(∇T_L) − f·|s|/|ΔS_V|,   |s| ≤ σ_y

Sem τ_r definido usa-se f = 1 (restauração elástica máxima; `RESTORATION_TAU` dá τ_r em s). A redução máxima é
σ_y/|ΔS_V| = 3,94e-8 m K (2–10 % de Γ); Γ fica 1,2–6,3 vezes o equilíbrio ao longo do lingote. A curva cinza
da Fig. 1 da cópia `Al3Cu5Nb01Fe_..._EP` usa, como as outras curvas MRB, o Γ constante no gradiente médio menos
a parte elástica restaurada (f = 1): 9,30·t_SL^0,379 µm. (Uma versão com Γ local ponto a ponto ficava acima
dos pontos nos gradientes altos e foi considerada errada pelo autor.)

## Arquivo não convertido

`Nucleation_model_FCC_A1_new_paper_proposal_luana_GL_data_indev.py` tem outra estrutura: calcula r_hom
resolvendo `func_r_new` para um gradiente dado (com a superfície elástica dentro dessa função) e para em
`quit()` antes do laço. Aplicar a superfície elasto-plástica exigiria reescrever `func_r_new`; a cópia
está em `nucleacao/nao_convertidos/` sem alteração.

## Observação sobre uma tabela original

No `..._Sim_Luane.py`, o cabeçalho de `tabela_paper.csv` tem 10 rótulos e as linhas têm 11 valores (não há
o rótulo `GT_het`): a coluna rotulada `GT_het_2nd` contém `GT_het_v` e o Γ de 2ª ordem é a última coluna.
A extração em `rodar_modelos_EP.py` leva isso em conta; o script original não foi alterado.

## Como rodar

```bash
cd ferramentas
python aplicar_elastoplastico.py        # recria as cópias de nucleação (ENERGY_MODE, lê sigma_y_ligas.txt)
python verificar_copias_elasticas.py <pasta temporária>   # modo elástico x original
python rodar_modelos_EP.py              # rodadas elasto-plásticas e Γ de 2ª ordem (--extract só extrai)
python ajustar_restauracao.py           # restauração elástica (Al-Cu-Nb-Fe), antes de aplicar_sdas.py
python aplicar_sdas.py                  # recria as cópias de SDAS
python rodar_sdas_EP.py                 # figuras de SDAS
python comparar_taxa_nucleacao.py       # taxa de nucleação: elástica x elasto-plástica
python comparar_curvas_EP.py            # demais curvas contra ∇T
python energia_superficie_alto_gradiente.py   # γ até 10⁶ K/m, valores em 10⁵ K/m
python gibbs_thomson_alto_gradiente.py        # Γ²ⁿᵈ até 10⁶ K/m, valores em 10⁵ K/m
```

Um modelo isolado, por exemplo com outro σ_y ou no modo elástico:

```bash
cd nucleacao
GM_SIGMA_YA0=0.05 python Nucleation_model_..._EP.py
GM_SURFACE_MODEL=elastic python Nucleation_model_..._EP.py
```
