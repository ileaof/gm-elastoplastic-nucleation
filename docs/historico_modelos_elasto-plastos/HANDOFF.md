# Dossier de continuidade — superfície de Gurtin–Murdoch elasto-plástica

Estado em 2026-10-22. Para um agente que assuma o trabalho sem o histórico da conversa.
Versão em inglês: `HANDOFF.en.md`.

## 1. Quem e o quê

* **Autor:** Prof. Ivaldo Leão Ferreira (UFPA, ileao@ufpa.br). Responder em **português**; papers em **inglês**.
* **Teoria do autor:** o coeficiente de Gibbs–Thomson é um tensor de campo térmico, Γ = A·∇T. A maximização de ΔG
  dá uma quadrática em r. A solução de 2ª ordem é a completa ("precisa"); a de 1ª ordem é uma aproximação.
  A energia de superfície vem da superfície de Gurtin–Murdoch (GM), com σ₀ < 0 (compressiva; os scripts guardam o
  módulo |σ₀|).
* **Objetivo deste trabalho:**
  - generalizar a superfície GM para elasto-plasticidade, em termodinâmica fora do equilíbrio;
  - aplicar a todos os modelos de nucleação e aos de SDAS (MRB);
  - escrever o paper.

## 2. Restrições e preferências (obrigatórias)

1. **Não alterar** nada em `Sim Luane Final/`; trabalhar só em `modelos_elasto-plastos/`.
2. **Não fazer push** no repositório GPU `nucleation-molecular-dynamics`. O repositório público
   `github.com/ileaof/gm-elastoplastic-nucleation` (local: `C:\Users\ileao\OneDrive\Documentos\gm-elastoplastic-nucleation`,
   commit `ddb3d00`, 24/10/2026) **foi sincronizado** a pedido do autor (será citado no paper): publicados
   `nucleacao/` (cópias `stored_total`, `gm_elastoplastic.py`, figuras, tabelas), `sdas/`, `ferramentas/` e
   `paper/` (inclui `nucleation_elastoplastic_AlSiMgFe.py` e os outputs), com os dois READMEs atualizados para
   `stored_total`. Fora do repo de propósito: o manuscrito `.docx`, o `HANDOFF.md` (interno), `Sim Luane Final`
   e os originais que só existem no workspace do autor.
3. **Valores físicos:** perguntar antes de assumir. Descrever fatos e perguntar; não apresentar o código do autor como
   "erro".
4. **SDAS:**
   - o autor diz que nenhum rótulo está trocado (1ª ordem aproximada, 2ª precisa);
   - a curva mais precisa é a de maior Γ (2ª ordem);
   - manter curvas e legendas originais.
   A correção feita nas cópias é só de *handles* (ver §5).
5. **Taxa de nucleação:** manter exp(+ΔG*/ΔG*_ref) como no original.
6. **Física adotada pelo autor:**
   - a deformação plástica da superfície é energia armazenada para a nucleação de fases;
   - o Γ elasto-plástico deve ser muito parecido com o elástico e crescer do mesmo jeito;
   - o Γ deformado plasticamente é maior que o de equilíbrio, o que permite a nucleação de outra fase induzida por
     deformação;
   - só a parte elástica é restaurada;
   - a mudança de fase é governada pela taxa de variação de densidade;
   - ΔT_SL deve ser o médio integral (como os "Mean Integral" dos scripts `Thermal_data_*`).
7. **Ambiente:**
   - Windows, Python 3.11, LAMMPS sem MPI (usar OpenMP);
   - **nunca** `taskkill //IM python.exe` (mata tudo);
   - em heredocs bash, `\n` em `\nabla` vira quebra de linha: editar LaTeX com a ferramenta Edit;
   - arquivos Word abertos no Word ficam travados: salvar com outro nome.

## 3. Modelo atual (decisões em ordem)

Arquivo central: `nucleacao/gm_elastoplastic.py`, classe `GMElastoPlasticSurface`, `step(r)` → (τ, s, γ, ∂γ/∂r).

* **Escoamento de área:**
  - critério |s| ≤ σ_y + Hα, com s = |σ₀|E_A^e;
  - return mapping e tangente consistente C·H/(C + H);
  - H = 0 (sem encruamento na liquidus).
* **Modo adotado `stored_total`** (escolha do autor, opção "c"):
  - fórmula: γ = γ₀(r_eq/r)² + |σ₀|(|E^e| + |E^p|) − q_diss = γ_el − σ_y·α;
  - referência fixa (r_eq não muda);
  - efeito: Γ fica 0,5–4,2 % abaixo do elástico.
* **Modos rejeitados** (mantidos para comparação):
  - `complement` (referência atualizada, Γ até −29 %);
  - `barrier`, `stored` (taxa caía monotonamente — "fisicamente estranho" para o autor);
  - `restored*`.
* **Restauração elástica** (`ferramentas/ajustar_restauracao.py`):
  - ε̇_v = β_s·Ṫ/ΔT_SL;
  - f = ε̇_vτ_r/(1 + ε̇_vτ_r);
  - Γ = Γ_EP − f·|s|/|ΔS_V|, com |s| ≤ σ_y.
  τ_r não foi definido; usa-se f = 1 (limite). A versão antiga (restaurar tudo até γ₀, τ_r = 7,84e3 s ajustado ao
  SDAS, erro 3,6 %) foi **rejeitada** pelo autor porque levava Γ ao equilíbrio.
* **σ_y** (`ferramentas/sigma_y_ligas.txt`): 0,0386 / 0,1572 / 0,0423 N/m. Vêm de um ajuste antigo no modo `barrier`
  e não são ajuste ao SDAS no modo atual.

## 4. Números atuais (modo `stored_total`)

| Liga | escoamento | Γ²ⁿᵈ el./EP (∇T médio) | SDAS 2ª el. / EP | γ em 10⁵ K/m el./EP |
|---|---:|---|---|---|
| Al-3Cu-5Nb-0,1Fe (3294 K/m) | 731 K/m | 1,0685e-6 / 1,0627e-6 | 9,38 / 9,37·t^0,379 | 5,58 / 5,53 J/m² |
| Al-0,8Si-0,6Mg-0,2Fe (7902) | 925 K/m | 1,5745e-6 / 1,5083e-6 | 6,71 / 6,66·t^0,387 | 4,45 / 4,24 |
| Al-3Cu-0,5Mg-0,2Fe (3708) | 562 K/m | 1,0973e-6 / 1,0878e-6 | — / 6,93·t^0,384 | 4,35 / 4,29 |

**Ponto que o autor pode retomar:** ele afirmou que "com o elástico não dá" 10⁵ K/m. Os fatos são:
* nas cópias (ΔT em forma fechada), o elástico e o EP vão até ~1,1×10⁶ K/m;
* quem não chega é o script original, porque `brentq(f_DT, 0.0001, 1)` limita ΔT ≤ 1 K e o número de pontos é fixo.
Isso foi dito a ele e escrito como fato no paper; ele ainda não respondeu.

## 5. O que foi feito nesta sessão (arquivos)

* **Cópias de nucleação:**
  - `nucleacao/*_EP.py` são geradas por `ferramentas/aplicar_elastoplastico.py`, com `ENERGY_MODE = 'stored_total'`;
  - `nucleacao/figure_layout.py` põe textos e legendas dentro dos eixos em todas as figuras dos scripts;
  - o gráfico 13 (taxa) tem os textos no canto superior direito.
* **Cópias de SDAS:**
  - `sdas/*_EP.py` são geradas por `ferramentas/aplicar_sdas.py`;
  - `legend_handles()` liga cada texto da legenda à sua linha (no original, a 2ª ordem é desenhada antes da 1ª e o
    `plt.legend` só com textos trocava os traços);
  - a curva cinza (Al-Cu-Nb-Fe) é o Γ constante no gradiente médio menos a parte elástica restaurada (f = 1):
    9,30·t^0,379. A versão com Γ local ponto a ponto foi considerada errada pelo autor.
* **Tensão superficial composta** (autor, 2026-10-23): τ = (σ₀ + s^e) + |σ₀|E_A^p + σ_y = τ_el + σ_y (H = 0), em
  `gm_elastoplastic.py` (`hist`, `tau_total`); figura `tensao_superficial_alto_gradiente.py`.
* `nucleacao/figure_layout.py` não impede a exibição: se falhar numa figura, avisa e a deixa como está. Depois de o
  autor relatar que as cópias de nucleação "não geram mais figuras", o teste aqui com qtagg mostrou 26 figuras criadas
  e entregues ao `plt.show()`, sem erro. Falta saber como ele roda (IDE ou backend) e qual mensagem aparece.
* **Figuras contra ∇T,** com o escoamento marcado (faixa sombreada = elástico, traço-ponto = escoamento):
  `comparar_curvas_EP.py`, `comparar_taxa_nucleacao.py`, `energia_superficie_alto_gradiente.py` e
  `gibbs_thomson_alto_gradiente.py` (Γ²ⁿᵈ até 10⁶ K/m), todas gerando em
  `nucleacao/figuras/`.
* **Correção:** `rodar_modelos_EP.py` agora usa a tabela mais recente de cada pasta. Antes lia CSVs de modos antigos.
* **Figuras da taxa nas 11 cópias de nucleação** (em `aplicar_elastoplastico.py`; o código das figuras originais da
  taxa contra ΔT foi mantido):
  - na figura contra ΔT (gráfico 13), a taxa elástica do mesmo script aparece junto, em azul. Ela é calculada por
    `figure_layout.elastic_rate()`, que roda o script em modo elástico num processo separado com `GM_DUMP_RATE` (o
    script grava a taxa e para antes da figura). Também aparece a faixa elástica com o escoamento;
  - foi criada uma figura nova da taxa contra ∇T (`figG`, escala log), com o elástico, o escoamento e o ∇T
    experimental médio. Ela fica no bloco `RATE_VS_GRADIENT`, logo antes do `plt.show()`;
  - figuras novas contra ∇T em toda a faixa: Γ²ⁿᵈ, γ e τ (`figure_layout.vs_gradient_figure`), com o elástico, o
    escoamento, o ∇T médio e os valores em 10⁵ K/m. `GM_DT_MAX` passou a ter padrão 30 K (varredura até ~10⁶ K/m);
  - na figura da taxa fica só o nome da liga da cópia (`ALLOY_TEXT`). O original da Al-Si-Mg-Fe tem também
    "Al-6wt%Cu-4wt%Si-0.11wt%Fe", que era a única visível no original.
* **Paper:**
  - arquivo: `C:\Users\ileao\OneDrive\Documentos\Paper Nucleation\Manuscript_ElastoPlastic_Surface_NonEquilibrium_Nucleation_v6_final.docx` (versão final para publicação, cópia idêntica da v6: quatro formulações da apresentação + soluções esféricas de ΔG_S e do seu potencial; v4 = notação geral dos slides 12–16 em forma elasto-plástica, v5 = + dedução pela métrica com variação da área);
  - 19 páginas, 28 equações OMML nativas, 6 figuras, 7 tabelas, em inglês; inclui a tensão superficial composta
    (Eq. 22) e a seção 5.4 com Γ, γ e τ das três ligas em toda a faixa de ∇T (até 10⁶ K/m) e a Tabela 6 em 10⁵ K/m;
  - v1 a v4 estão **desatualizadas**;
  - fontes em `paper/`: `build_paper.py` (python-docx + latex2mathml + MML2OMML.XSL do Office) e `make_figures.py`;
  - render de verificação: Word COM → PDF → PyMuPDF;
  - **ScholarOne rejeitou o v6_final** ("Caught XML parse error / Validation error"): o texto puro tinha
    105 hifens não-quebráveis (U+2011, nas designações `Al‑3Cu‑…`) e ~200 sobrescritos Unicode (U+207x:
    os "10⁻⁶ K/m"). O `sanitize_text()` no `build_paper.py` (commit `6339354` no repo) converte U+2011 em "-"
    e sobrescritos Unicode em runs de sobrescrito real (com dígitos ASCII); rebuild de verificação saiu
    com 0 caracteres ofensivos e as 303 equações intactas. Para reenvio:
    `Manuscript_..._v6_final_clean.docx` (v6_final saneado, mesmo conteúdo). As figuras 26 e 27 dos tweets
    vêm de `paper/nucleation_elastoplastic_AlSiMgFe.py` (janelas 26 = painel 3×3 de curvas, 27 = painel 2×3
    de energias e entropias).
* **README.md** atualizado para `stored_total`.

## 6. Pendências e sugestões

0. **Decisão do autor (2026-10-23):** publicar com γ = γ₀(r_eq/r)² − s (forma dos scripts) — é a forma de todos os
   resultados do paper v5. A forma esférica da métrica (Eq. 89 do paper da Materials, com integral na área) e a leitura
   de Σ na integral de A₀ a A₁ ficam para depois (testes em `ferramentas/formas_energia_gt_sdas.py` e
   `sdas_metrica_area_EP.py`). Também em aberto: qual Γ o autor considera de 2ª ordem no script de SDAS (os comentários
   que ele colou chamam `z3n_1st_Eq`, com Γ = 3,887e-7, de 2ª ordem; o script e a nucleação dão 1,0685e-6 como 2ª).

1. **SDAS:** a curva EP de 2ª ordem (9,37) e a restaurada (9,30) coincidem com a 2ª ordem elástica (9,38), que fica
   acima dos pontos. A 1ª ordem (7,67) "bate", mas o autor a considera só uma aproximação. Ele disse que vê o SDAS
   depois. Não ajustar τ_r nem inventar leis sem perguntar.
2. **Paper v2:**
   - revisar com o autor;
   - a referência Mendes et al. (2023) está incompleta (placeholder);
   - o SDAS não está no paper;
   - figuras: Fig. 1 = `curvas_EP_Al3Cu5Nb01Fe.png`; Figs. 3–5 = `gibbs_thomson_1e5.png`, `energia_superficie_1e5.png`,
     `tensao_superficial_1e5.png` (de `nucleacao/figuras/`).
3. **Repositório GitHub:** atualizar com o modo `stored_total` e o paper só se o autor pedir.
4. **Texto errado no original:** o gráfico 16 da Al-Cu-Nb-Fe tem o nome de outra liga no texto (vem do original). O
   autor foi perguntado e não respondeu.
5. **H₂O:** σ_y provisório (0,0695 N/m).

## 7. Como reproduzir tudo

Ver "Como rodar" no `README.md`. Ordem:
1. `aplicar_elastoplastico.py`;
2. `rodar_modelos_EP.py`;
3. `ajustar_restauracao.py`;
4. `aplicar_sdas.py`;
5. `rodar_sdas_EP.py`;
6. scripts de figuras.

O paper: `python -I paper/build_paper.py`, depois de `paper/make_figures.py` e de copiar `curvas_EP_Al3Cu5Nb01Fe.png`
para `paper/Fig1_curves.png`. Os caminhos de saída estão no topo dos scripts.

## 8. E6 (recristalização por MD) — o que rodar na máquina GPU e o que trazer de volta

Contexto: o professor quer provar que a deformação plástica armazenada (`stored_total`, σ_y·α) reaparece como
potencial na recristalização. Este laptop (i7-1360P, sem NVIDIA) já completou os pilotos e a série pequena
(2 048 átomos, 500–700 K da liga; 500–900 K do controle) — tudo no repo, commit `e98d66b`, pasta
`recristalizacao/` (o `PLANO.md` tem os números e o registro da noite). O que falta é pesado e roda na
**máquina com NVIDIA**. Procedimento detalhado: `recristalizacao/GPU.md` no repo.

**Na máquina GPU** (LAMMPS; MEAM roda em MPI+OpenMP, EAM aceita `-sf gpu`):
1. `git clone https://github.com/ileaof/gm-elastoplastic-nucleation.git`
   `cd gm-elastoplastic-nucleation/recristalizacao`
2. Ajustar `LMP` / `NP` (ranks MPI) / `NT` (threads) no topo do `rodar_gpu.sh`.
3. `bash rodar_gpu.sh` — G0 validação (~1 min); G1 célula grande 16 384 átomos (~1–3 h);
   G2 série 500–900 K na célula grande (~10–16 h); G3 fronteira fina 750/775/800 K × 3 sementes
   (~4–8 h); G4 recozimento longo 1 ns em 750 K (~10–20 h). Por seção:
   `SECTIONS="G1 G3" bash rodar_gpu.sh`.
4. `bash juntar_resultados.sh` → cria `results_e6.tar.gz` (logs, `stored/perfil .dat`, figuras; os
   dumps ficam lá — grandes).

**O que mandar de volta para ESTA máquina** (`C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\modelos_elasto-plastos\recristalizacao\`):
* **Obrigatório:** o `results_e6.tar.gz` (Any OneDrive/transferência/USB), descompactar na pasta acima. Se a máquina
  GPU tiver as credenciais do GitHub, pode ser por `git push` (pasta `results_e6/`) e pull daqui.
* **Opcional (para o OVITO):** os `recozido_*.dump` de interesse e o `deformed_jelinek.dump` grande.
* **Não precisa mandar:** os `data_*.data` (o G1 recria os deles), nem os logs brutos grandes.

Depois de voltar: rodar `python analyze_recristalizacao.py` na pasta e reportar ao professor
(E_st(ε), δ, limiar de recristalização limiar/estatístico, comparação com σ_y·α/δ).
Perguntas ao professor na volta: (a) novo desenho para 800–900 K (dt menor, barostato, célula);
(b) termo térmico (Γ) no recozimento; (c) δ medida ou imposta; (d) critério de "recristalizado".

## 9. Estado do paper e do repositório (28/10/2026) — LEIA ANTES DE TOCAR NO MANUSCRITO

**O manuscrito é o `Manuscript_..._v7_total.docx`** (OneDrive/Paper Nucleation; a linhagem que o
professor fixou: "use a versão 7 para o build total, pois a v7 é a única que fala sobre recristalização").
Contém 12 figuras, 11 tabelas, 45 equações, 132 runs saneados, e as seções 5.7–5.10 com a campanha MD:
5.7/5.8 (E6: pilotos, série pequena, descarregamento E6.4, fronteira fina — Figs 7–10, Tabelas 8–9),
5.9 (E6.5: célula de 62 500, as séries 800/850/900 NVT — Fig 11, Tabela 10), 5.10 (E6.6: descarregamentos
γ=0,10/0,30 + anneals NPT/NVT, platôs — Fig 12, Tabela 11).

**Regras de continuidade (não negociadas):**

1. **Os builds nunca sobrescrevem a v7_final** (25/10, o marco) — o build atual é o `v7_total.docx`.
   A cópia de rastro do passe-final: `..._v7_total_2710a.docx`.
2. **A `..._v8_final.docx` (26/10) é sobressalente** (o meu build intermediário antes da ordem da
   linhagem v7); a v7_total contém tudo dela + E6.6 + o passe final. Apagar só com ok do professor.
3. **Edições diretas no Word são PERDIDAS no próximo build** (o conteúdo vem do código): o professor
   abriu/re-salvou o v7_total em 27/10 21:44 SEM mudança de conteúdo (verificado por diff); se ele
   editar de verdade, PORTAR a mudança para o `build_paper.py` antes de qualquer rebuild.
4. **Física adotada e publicada**: o γ=0,10 na célula grande é elástica (unload ε=0,00, resto ~0) —
   o excedente da Tabela 10 é energia do estado carregado; o plástico entra com γ=0,30 (desliza,
   descarrega em ε=0,252, resto 0,0119 eV/átomo = 1,17 J/m²); a liberação é gradual (47–67 % em
   100–300 ps) sobre platô de +3–6 meV/átomo (3–6×10⁷ J/m³ ≈ a ordem de σ_y^s·α/δ com δ = 0,19 nm);
   a catástrofe das células pequenas não existe na grande; TODA fração de defeito em dump quente =
   excesso sobre a régua do mesmo barostato (NVT 1,1→28,6 %: 500–900 K; NPT 42,6/52,7/62,2 %).
5. As **siglas definidas na primeira citação** (FCC/MEAM/EAM/LAMMPS na 5.7; NVT/NPT/CNA na 5.9; ΔPE e
   ΔE_st na legenda da Tabela 10; os códigos E6.5/E6.6 ancorados: "E6.5 in the records of the
   repository cited here"). Manter assim em qualquer texto novo.
6. O evento das tabelas do paper = o passo de INÍCIO do par do maior salto, com o primeiro diff do
   thermo saltado (aquecimento); o recuo longo = "1-2 ps" (a linha antiga "1 ps" e o artefato "0 ps"
   do 775/s11 foram corrigidos em 7184011).
7. **Pushes**: `f3e5ce0` (as simulações consolidadas, 109 arquivos), `7184011` (o passe final:
   Fig 10 em inglês a partir do parse da própria tabela + abstract/conclusões), `36e3e8e` (as siglas).
   O abstract passou de 1461 para 1754 caracteres COM o ok do passe — o professor pode vetar/reredigir.
8. Fora do repo por desenho: o `.docx`, este dossiê, o `E65_registro.md` (interno; espelhado na caixa
   como REGISTRO_E65.md), a galeria `galeria_recristalizacao.html` (local). Os dumps grandes (53 MB × N)
   ficam na caixa (~/e6_rec/) para o OVITO.
9. **As duas casas do pipeline do paper (amb nesta máquina):**
   - ORIGINAL: `modelos_elasto-plastos\paper\` (o workspace; as Figs 1–6 de 23/10; por lá saíram a
     v6-final e a v7-final de 25/10, com o `make_e6_figures.py` = o gerador das Figs 7–10);
     **sincronizada em 28/10** com o `build_paper.py` do repo (OUT = v7_total nas duas).
   - ATUAL: `OneDrive\gm-elastoplastic-nucleation\paper\` (o repo; por lá correram os builds do
     E6.5/E6.6: v8 (residual), v7_total e o passe final). O `.docx` de SAÍDA sempre foi
     `OneDrive\Documentos\Paper Nucleation\` (o OUT de todas as versões).
   - **Buraco de reprodução**: `make_e6_figures.py` (o gerador das Figs 7–10) só existe no workspace —
     NÃO foi publicado no repo (as PNG sim). Incluir no próximo push com ok do professor.
10. **A NOTE1 ficou fora do ar desde 28/10 (desligada pelo professor, com tudo entregue)**: ssh
    `labmo@192.168.15.96 -p 2222` só volta quando ela for ligada. Nada depende dela.

Memória persistente do agente:
`C:\Users\ileao\.claude\projects\c--MacbookPro-ufpa-PPGEM-SDAS-Al-Si-Mg-Fe\memory\` (índice `MEMORY.md`).
