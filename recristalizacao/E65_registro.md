# Registro E6.5 — análise do recozimento a 900 K na caixa GPU (26/10/2026)

Estado: E6 fechado (25/10, repo commits até `d109af0`). Na noite de 25/10 a caixa (NOTE1) rodou
um experimento NOVO, desenhado por ela sob o professor: **E6.5** — controle Mendelev (EAM Al puro),
célula grande `nx=25` = **62 500 átomos**, FCC criado + equilíbrio 10 ps a 300 K + cisalhamento em
linha (erigido dentro do script, `erate 0.010`, 10 ps) + rampa blindada 300→900 K (5 000 passos) +
recoximento **NVT 900 K × 1 ns**, dumps a cada 25 000 passos (50 ps). Terminou 21:29 de 25/10
("E6.5 completa", 33:55 de relógio), sem átomos perdidos. A caixa estava ociosa esta manhã
(nenhum lmp rodando, carga 0,21). O dump não tinha sido analisado — **esta é a análise que faltava**.

## 1. Fatos da análise (dump `reciz_mendelev_T900.dump`: 20 frames, 20→970 ps)

1. **O rótulo "5" do `cna/atom` a 900 K é artefato térmico, não estrutura.** Prova por controle
   (`in.ctrl_fcc.lmp`, mesmo setup, FCC puro SEM cisalhamento, 30 000 passos, 1:54):
   - 300 K: 99,99% FCC (62 494 de 62 500) — o CNA está limpo;
   - 900 K: **71,5% FCC / 28,5% "5"**, estável (+40 ps: 71,3/28,7).
   O rótulo "5" num FCC puro a 0,97·T_m vem da vibração térmica com `rc = 0,854·a₀`.
   A 300 K o mesmíssimo compute não marca nada ⇒ qualquer fração de defeito a 900 K
   deve ser lida como **excesso sobre o baseline de 28,6%**.
2. **Excesso real do cisalhamento: ~21,8% dos sítios** (bruto 50,4% − baseline 28,6%),
   estável nos 970 ps (mistura de rótulos nos vizinhos: 65/35 — leve correlação, sem manchas grandes).
3. **PE (média/frame, do CSV):**
   - equilíbrio 300 K: −3,4107 → −3,3733 eV/átomo (FCC puro a 300 K no controle: −3,3733 ✓);
   - fim do cisalhamento (γ = 0,10, 300 K): −3,3516 ⇒ **+0,0217 eV/átomo armazenado na deformação**;
   - hold 900 K: −3,2598 constante de 20 a 970 ps; puro a 900 K: −3,2867 (controle);
     ⇒ **+0,0267 eV/átomo remanescente ≈ +2,6×10⁸ J/m³** (célula 1 034 131 Å³).
     (Faixa de referência do PLANO.md §2 para metais deformados: 10⁵–10⁷ J/m³ — o remanescente do MD
     fica ~26× acima do topo; ressalvas: cisalhamento a 10⁹ s⁻¹, **sem descarregamento elástico**, NVT.)

## 2. Fatos da história termodinâmica (thermo do `rec_gpu.log`)

| passo | fase | T (K) | press (kbar) | PE/átomo (eV) | Pxy (kbar) |
|---:|---|---:|---:|---:|---:|
| 0 | FCC criado | 300 | 2,6 | −3,4107 | ~0 |
| 5000 | fim do equilíbrio | 297 | 17,5 | −3,3733 | ~0 |
| 10000 | fim do cisalhamento | 301 | 25,8 | −3,3516 | **−42,4** |
| 15000 | fim da rampa | 893 | 53,1 | −3,2610 | **−33,0** |
| 515000 | fim do hold | 895 | 53,4 | −3,2598 | **−33,0** |

* No cisalhamento, Pxy sobe com saturação (−8 → −42 kbar em 10 ps): curva de escoamento a 300 K
  (taxa 10⁹ s⁻¹); PE sobe +0,022 eV/átomo.
* Na rampa, Pxy relaxa de −42,4 para −33,0 kbar (−22%) — recuperação durante a subida.
* No hold, **nada evolui em 970 ps**: PE constante, Pxy constante, frações CNA constantes.
  O material cisalhado resiste (Pxy −33 kbar não desliza); **sem recuperação, sem
  recristinação, sem evento** — estado meta-estável.
* Confinamento NVT: pressão 47–53 kbar (célula rígida no volume de a₀ = 4,045 construído a 300 K,
  sem a expansão térmica de ~1,4% linear até 900 K).

## 3. Discrepância no script (fato, não decisão)

O comentário do `in.rec_gpu.lmp` diz "cisalhamento a 0.20 (10 ps)", mas `erate 0.010 × 10 ps = 0,10`
e o tilt do box confirma: 10,1125/101,125 = **γ = 0,10**. Se o professor queria 0,20, a corrida
ficou com metade da deformação.

## 4. O que o E6.5 responde e o que fica aberto

* **Responde (a) do HANDOFF_E6_GPU §4**: o desenho novo (célula 62 500 + rampa blindada) **não ejeta**
  a 900 K — a catástrofe 800–900 K da célula pequena (2 048) e da 16 384 (2 ps) foi desenhada embora.
* **Fica aberto**: com γ = 0,10, sem descarregamento e NVT, não há sinal de recristalização em 1 ns
  a 900 K — a energia armazenada (+2,6×10⁸ J/m³) permanece retida; o hold não a libera.

## 5. Arquivos

* **Caixa** (`labmo@NOTE1:~/e6_rec/`, ssh -p 2222): `in.rec_gpu.lmp`, `rec_gpu.log`,
  `reciz_mendelev_T900.dump` (53 MB — para OVITO, fica lá), `in.ctrl_fcc.lmp`, `ctrl_fcc.log`,
  `ctrl_fcc_T{300,900,900_late}.dump`, `analyze_e65.py`, `check_e65_geom.py`, `e65_analise/`
  (`e65_resumo.csv`, `fig_E65_900K.png`), `REGISTRO_E65.md` (este registro).
* **Laptop** (`modelos_elasto-plastos/recristalizacao/`): idem os pequenos + `e65_analise/` +
  `fig_E65_resumo.png` (3 painéis: PE com baseline, fração "5" com excesso, Pxy residual) +
  `fig_e65_final.py`, `E65_registro.md`.
* O clone da caixa (`~/e6_rec/repo`) está em `9944458` (4 commits atrás do laptop: `a4a5332…d109af0`,
  que trazem a Tabela 9/Fig. 10). Nada do E6.5 foi commitado — **só com pedido do autor**.

## 8. Cadeia E6.6 lançada (26/10 18:56, autorização do professor: "faça, mas deixe em background...")

Ao autorizar a sequência 5→2→3, o professor mandou deixar o cálculo rodando na caixa mesmo com o
ssh desconectado. Lançada com `nohup ... < /dev/null` (desanexada de fato; testado ontem). A cadeia
(`ser_e66.sh`, 21 corridas, todas `-sf gpu`, ordem):

1. **Controles NVT frios** 500/600/700 K (~2 min cada) — fecha a régua do baseline na faixa toda;
2. **Controles NPT** 800/850/900 K (~42 min cada, 1 ns + trajetória) — a régua nova para as corridas
   a zero de pressão (as NVT têm a régua do controle NVT; comparar estados de volume livre contra
   a régua NVT confundiria com o efeito de volume);
3. **Item 2 (discriminante)**: γ = 0,10 *descarregado* (reverso em degraus até pxy ≥ 0, protocolo
   E6.4, o cruzamento impresso como `DESCA_CROSS` no log) + **recrozimento NPT P=0 de 1 ns**,
   800/850/900 K × sementes 101/202 (6 corridas);
4. **Item 2b**: γ = 0,30 descarregado + NPT, 1 semente (o estado que decaiu nas Figs 8/10 na
   célula grande de 16 384 — agora na de 62 500, para testar se a catástrofe sobrevive ao redesenho);
5. **Item 3**: γ = 0,20 e 0,30 *carregados* + NVT 1 ns (a variação em γ do desenho E6.5), 6 corridas.

Nova entrada de LAMMPS: `in.descar_npt.lmp` (skeleton do E6.5 + loop de descarregamento com
`variable pnow equal pxy` + rampa/hold NPT), `in.ctrl_npt.lmp` (régua NPT), `in.rec_largo.lmp`
(E6.5 com `gam` variável). Tempo previsto pela taxa medida: **~13–14 h** (termina na manhã de 27/10).
Entrega prevista: análise (excesso sobre baseline **por régua certa** — NVT para as largas, NPT para
as descarregadas), as 21 corridas em registro, e a extensão da seção 5.9 do paper só depois de os
números existirem e o professor aprovar.

## 9. Resultados da cadeia E6.6 (27/10, 12:08; análise por analyze_e66.py + fig_e66_serie.py)

21/21 corridas ok, **zero átomos perdidos em qualquer um** (nem o γ=0,30 descarregado + NPT). Resultados:

1. **O cisalhamento a γ=0,10 na célula de 62 500 é essencialmente elástico** (taxa 10¹⁰ s⁻¹):
   o descarregamento reverso (protocolo E6.4) volta ao cruzamento pxy ≥ 0 em **ε_cruz = 0,000**
   (as 2 sementes × 3 T idênticas), com resto plástico **0,02–0,03 J m⁻² = ~zero** (nível do ruído de
   equilíbrio; = o E6.4 do estado elástico de pico, −0,007 J/m²).
   ⇒ **Correção de leitura: o "+0,026 eV/átomo" da série E6.5 (Tabela 10) é a energia do estado
   CARREGADO (tilt fixado pela geometria, σxy −33 a −35 kbar), não armazenado plástico.**
2. **γ=0,30 descarregado**: cruzamento em **ε = 0,252** (= a Tabela 8 do E6.4: 0,250 ✓) com resto
   plástico **+0,0119 eV/átomo = 1,17 J m⁻²** (o mesmo estado na célula pequena, protocolo em degraus:
   0,184 J m⁻² — o cisalhamento contínuo único armazena ~2× por átomo do que os degraus; 0,0059 eV/at).
3. **A catástrofe NÃO volta na célula grande**: o estado 0,30-descarregado sob NPT (P=0) em 800/850/900 K
   não perde átomo, não tem evento (ΔPE por passo máximo = ruído 0,002–0,003 eV/at) e **libera
   gradualmente ~0,0045–0,0056 eV/átomo em 100–300 ps** — recuperação por difusão, não o colapso de 2 ps
   do Jelinek 16 384 nem o de 55–60 ps do Mendelev 2 048: **a catástrofe era fenômeno da célula pequena
   (2 048) e da fronteira violenta 750–800 K nela**.
4. **Platô comum**: os três estados com armazenagem plástica real (0,20 carregado, 0,30 carregado,
   0,30 descarregado) liberam **47–67 % do excesso inicial** sobre a régua e estabilizam em
   **+0,003 a +0,006 eV/átomo acima da régua**, que persiste até os 1 ns — o resto plástico que não
   sai em 1 ns; nenhuma recristalização completa (a fração de excesso CNA sobre a régua correta cai
   de +9/+11/+6 % para +4 a +6,5 % e fica).
5. **Régua completa do artefato CNA (controles puros, 62 500 átomos):**
   NVT: 1,1 % (500 K) / 3,7 / 9,0 / 17,7 / 22,8 / 28,6 % (900 K); **NPT (P=0): 42,6 / 52,7 / 62,2 %**
   (800/850/900 — a vibração a volume livre desclassifica MAIS; press ≈ 0 verificado). Toda fração de
   defeito em dump quente agora tem régua para os dois ensembles.
6. As sementes 101/202 do descarregado-γ0,10 replicam a régua NPT (PE 0,0009 eV/at — replicação ✓).

Figura: `e66_analise/fig_E66_serie.png` (2×3: PE com as duas réguas; excesso CNA sobre a régua certa).
Saídas de dados: `e66_analise/{reguas.csv, descarrega.csv, largo.csv, frames_*.csv}` (box + laptop).

## 11. Build total na linhagem v7 (27/10, pedido do professor: "Use a versão 7 para o build total")

O professor definiu a linhagem: o build TOTAL (E6 + E6.5 + E6.6) sai como **v7** — e como builds nunca
sobrescrevem, saiu como **`Manuscript_..._v7_total.docx`** (4,59 MB; o `v7_final` de 25/10 ficou
intacto como registro). Conteúdo do passe:

1. **Emenda da 5.9**: a frase do "stored deformation survives" agora diz que o estado CARREGADO
   γ=0,10 permanece e que a sobra é a energia da geometria carregada (a prova vem da 5.10);
2. **Subseção 5.10 nova** (3 parágrafos): (i) o descarregamento na célula de 62 500 — γ=0,10 volta ao
   ε=0,00 elástico (resto 0,02–0,03 J/m², ruído), γ=0,30 desliza (pico 7,3 GPa perto de ε≈0,18, cai para
   0,7–1,2 GPa) e cruza em **ε=0,252** (= 0,250 da Tabela 8 ✓) com resto **0,0119 eV/átomo (1,17 J/m²)**,
   ~2× o por-átomo da pequena (0,0059); (ii) as anneais contra a régua do PRÓPRIO barostato (NPT contra
   régua NPT, NVT contra NVT — o teste de volume livre), zero átomos perdidos nas 21, liberação gradual
   de 47–67 % em 100–300 ps, platô comum de +3 a +6 meV/átomo que os 1 ns não consomem, a catástrofe
   da Fig 10 ausente da célula grande, as sementes 101/202 replicando a régua NPT (±1,7 meV/átomo);
   (iii) o fecho da escala: o platô = 3–6×10⁷ J/m³ ≈ a ordem do potencial σ_y^s·α/δ com o δ medido
   (3,2×10⁷ com os parâmetros do controle) — **a leitura desse acordo fica para o professor** (o σ_y
   é do fit antigo; δ = 0,19 nm);
3. **Fig 12 nova** (`paper/Fig12_E66_unload_release.png`): (a) as curvas de ida/volta σxy(ε)
   [elástica re-fazendo o caminho; 0,30 com glide e cruzamento marcado] + (b–d) os platôs do excesso
   sobre a régua em 800/850/900 K com a % liberada em caixa;
4. **Tabela 11 nova** (15 linhas, computada do build): as 12 corridas de estado + os 3 restos do
   descarregado — cruzamentos, restos, excessos no início/1 ns, % liberada, σxy médio.

Verificações: 12 imagens inline (Figs 1–12 ✓), 45 equações, 131 runs saneados, ordem
5.7→5.8→5.9→5.10→Conclusões correta; a fidelidade v8↔v7_total mostra SÓ a frase emendada + o bloco
5.10. Dados no repo para o build: `recristalizacao/e66_data/` (18 logs) + `e65_serie26/` (réguas).

**Ficam com o professor**: (a) eixos da Fig 10 em português (não mexi); (b) a discussão do platô contra
σ_y^s (o texto diz que a leitura é dele, pela calibração de σ_y); (c) o abstract/Conclusões seguem sem
os resultados novos; (d) commit e push feitos no pedido dele ("simulações consolidadas"): **f3e5ce0**
(109 arquivos, d109af0→f3e5ce0, main), working tree limpa, e a caixa puxou o mesmo hash — o repo
público está na citação do paper atualizado.

## 13. Passe-final do build (27/10, autorizado: "então deixe tudo pronto pode fazer esse build")

Os três itens abertos entraram e o `v7_total.docx` foi reconstruído no lugar (o estado anterior
estacionado como `..._v7_total_2710a.docx`; o push do passe = **7184011**):

1. **Fig 10 em inglês** (`make_fig10_en.py`) e derivada do MESMO parse que constrói a Tabela 9
   (fonte única). Com a correção, o salto de primeira linha do thermo deixa de ser lido como evento:
   o 775/s11 volta ao evento de 19 ps (o +88 eV no passo 0 era artefato) e o recuo longo fecha em
   **"1-2 ps"**, no lugar do "1" que discordava do texto ("2 ps") — a inconsistência que eu havia
   sinalizado. Tabela 9 agora: 750: 41/64/17; 775: 19/limpo/66; 800: limpo/limpo/56; longo 1-2.
2. **Abstract** (1754 caracteres, de 1461): uma frase do MD (liberação gradual 47-67 % em 100-300 ps,
   platô persistente 3-6×10⁷ J m⁻³ = a ordem de σ_y^s·α/δ com δ medido). As frases são minhas, com o
   ok do passe — o professor pode vetar/reredigir.
3. **Conclusões**: bula nova da separação de protocolos na célula de 62 500 (γ=0,10 elástico → ε=0,00;
   γ=0,30 desliza e descarrega em ε=0,252, resto 0,0119 eV/átomo = 1,17 J/m²; liberação 47-67 %, platô
   +3-6 meV/átomo que 1 ns não consome; a catástrofe ausente da célula grande). **5.10** fecha com a
   concordância explícita (fator ~2, a mesma ordem da §5.7).

Verificado no .docx final: 12 figuras (Figs1-12), 45 equações, 132 runs saneados, abstract com
'plateau' ✓, a bula da célula grande ✓, a ordem das figuras/tabelas correta.

## 14. Siglas definidas (27/10, "faltaram definir as siglas como NPT e outras") + push 36e3e8e

* §5.7 (texto pre-existente): FCC = face-centered cubic; MEAM = modified embedded-atom method (Jelinek);
  EAM = embedded-atom method (Mendelev); LAMMPS = Large-scale Atomic/Molecular Massively Parallel
  Simulator (o "code" soletrado na primeira citação; a sigla nua saiu do texto).
* §5.9: NVT soletrado ("the number of atoms, volume and temperature held (NVT)"); NPT soletrado na
  primeira citação ("NPT, the ensemble in which the number of atoms, the pressure and the temperature
  are held, here at zero pressure"); CNA linkado ("the common-neighbor analysis (CNA)"); os códigos
  E6.5/E6.6 ancorados ("E6.5 in the records of the repository cited here" — o paper cita o repo e os
  códigos apontam para os registros dele).
* Tabela 10: a legenda define ΔPE (o excedente de energia potencial sobre o controle da própria T) e
  ΔE_st (o mesmo como densidade volumétrica). Fig 12 regenerada com "reference" no lugar de "ruler"
  (painéis e legenda; o registro de linguagem informal saiu do texto e da figura).
* Verificação no docx construído: zero 'ruler' em resto; (FCC)/(MEAM)/(EAM)/(CNA)/(NVT)/NPT/LAMMPS
  todos presentes; 12 figuras; 132 runs saneados. push: 7184011→36e3e8e.

## 15. Registro no HANDOFF.md e README.md (28/10) + push a310587

Pedido: "registre isso no Handoff.md e README.md". Feito:

* **HANDOFF.md** (interno, `modelos_elasto-plastos/`): nova §9 "Estado do paper e do repositório"
  com o manuscrito corrente (v7_total: 12 figuras/11 tabelas/45 equações, seções 5.7–5.10), as 8
  regras de continuidade (builds nunca sobrescrevem a v7_final; a v8 = sobressalente com apagamento
  pendente; **edição direta no Word se perde no build — portar para o build_paper.py**; a física
  adotada publicada; as siglas na primeira citação; a convenção do evento; o histórico de pushes;
  o que fica fora do repo).
* **README.md + README.pt-BR.md** (o par público): a linha nova do `recristalizacao/` na tabela de
  layout (as séries E6/E6.5/E6.6, os logs crus das Tabelas 8/9/11, as réguas, os scripts); a linha
  do `paper/` ampliada (5.7–5.10, Figs 7–12, Tabelas 8–11, a saída v7_total); o bloco de "Como
  rodar" das cadeias (ser_e65/ser_e66 + analyze_e66 na caixa) e a convenção dos códigos de série
  (E5, E6.1–E6.6 = os registros, ancorados ao repo); o Estado com o resultado MD (réguas,
  descarregamentos, platô 3–6×10⁷ J/m³).
* push: 36e3e8e → **a310587** (2 arquivos, 42 linhas). A caixa segue sincronizada pelo pull.

## 16. A origem do paper (28/10, pergunta: "originalmente, em que diretório foi gerado?")

Duas casas do pipeline do paper, as duas NESTA máquina:
* **ORIGINAL: `modelos_elasto-plastos\paper\`** — as Figs 1–6 com data de 23/10 12:06; por lá saíram
  a v6-final e a **v7-final (25/10, 20:48)**; o `make_e6_figures.py` (25/10 20:13) = o gerador das
  Figs 7–10 **e não está no repositório** (as PNG sim — buraco de reprodução; incluir no próximo
  push com ok do professor);
* **ATUAL: `gm-elastoplastic-nucleation\paper\` (OneDrive, o repo)** — os builds do E6.5/E6.6 (a v8
  residual, a v7_total e o passe final) correram por aqui, com o OUT = v7_total.
* O **.docx de saída SEMPRE** = `OneDrive\Documentos\Paper Nucleation\` (o OUT nunca mudou; OneDrive
  sincroniza a saída).
* **Ação**: o `build_paper.py` do repo foi copiado para a casa original (28/10) — OUT = v7_total nas
  duas; nenhum agente monta por um código velho.
* A **NOTE1 está desligada** desde hoje (o professor a desligou com tudo entregue): o mirroring do
  registro fica pendente para quando voltar; nada depende dela.

## 12. Commit e push (27/10, a pedido: "é hora de um commit and push")

**f3e5ce0** — "E6.5/E6.6: the 62 500-atom series, the hot-CNA baselines, the unloadings and the release
plateaus" — 109 arquivos, 51 986 linhas: os inputs (in.rec_gpu, in.ctrl_fcc/npt, in.descar_npt,
in.rec_largo, ser_e65/66.sh), as réguas e dados (e65_serie26/, e66_analise/ com régua/descarrega/largo/
frames + fig_E66_serie.png, e66_data/ com os 18 logs, baseline_puro.dat), os crus das Tabelas 8/9
(descarga/deform logs e gpu2/ — o gap de reproducibilidade do build fechado), os scripts de análise
(analyze_e65/e66, check_e65_geom, escreve_baseline, make_fig_E65/E66, fig_e65_serie) e o paper/
(build_paper.py emendado + Fig 11/12). PLANO.md ganhou §9 (o registro público E6.5/E6.6) e GPU.md
ganhou o capítulo das cadeias novas. FORA por desenho: o .docx, o E65_registro.md (interno), a galeria,
os dumps grandes (ficam na caixa para OVITO). A caixa puxou (repo em f3e5ce0).

## 10. Alimentação do paper (26/10, a pedido do professor)

* **Criação da Fig 11 do paper** (`paper/Fig11_E65_series.png`, repo): a série E6.5 em 3 painéis
  (a: PE com a linha do controle puro por T; b: rótulo "5" com baseline e o excesso sombreado;
  c: Pxy residual), rótulos em inglês, matplotlib padrão — `make_fig_E65.py` (laptop, recristalizacao/).
* **Subseção 5.9 nova no `build_paper.py`** (depois da 5.8, antes das Conclusões), com 4 parágrafos:
  o desenho (62 500, cisalhamento em linha a 10¹⁰ s⁻¹, rampa cega, NVT — sem átomo perdido), a
  descoberta do **baseline do CNA quente** (17,7/22,8/28,6% no FCC puro a 800/850/900 K; 99,99% a 300 K;
  frações de defeito em dump quente só com excesso sobre baseline; energias não afetadas), os resultados
  (ΔPE +0,025–0,027 eV/átomo = 2,4–2,6×10⁸ J/m³; σxy −33/−35 kbar residual; nada de liberação em 1 ns) e
  a leitura factual do contraste com as Figs 8/10 (estado descarregado ε=0,30 × carregado γ=0,10; NPT × NVT).
* **Tabela 10** computada do build, os dados pequenos copiados para o repo em
  `recristalizacao/e65_serie26/` (`e65_resumo_T{800,850,900}.csv`, `rec_gpu_T{800,850,900}.log`,
  `baseline_puro.dat` de `escreve_baseline.py`).
* **Saída**: `Manuscript_..._v8_final.docx` (novo arquivo, 4,05 MB; o v7_final ficou intocado).
  Verificação de fidelidade (`fidelidade_v8.py`): v7 = 225 parágrafos, v8 = 234; os 9 novos são
  TODOS do bloco 5.9 — nenhum texto do v7 foi removido ou alterado; 11 figuras inline;
  Tabela 10 = os números da análise.
* **Aviso de processo descoberto ao construir**: o build lê `descarga_*.log`, `deform_*.log` e
  `gpu2/recozer_jelinek_*.log` de dentro do REPO e descarta em silêncio a Tabela 8, a 9 e a Fig 10 se
  não estiverem (eram gitignored e a pasta limpa). Restaurei as cópias no repo (laptop → repo) antes
  do rebuild. Em builds futuros conferir no output do `verifica_v8.py` que Tables 8/9/Fig 10 entraram.
* **Bandeiras para o professor**: (i) Tabela 9/linha do recuo longo: início do evento = 1 ps na linha,
  o texto do paper diz "2 ps" — inconsistência pré-existente do v7 (não mexi: é texto dele); (ii) eixos
  da Fig 10 em português num paper em inglês; (iii) abstract e Conclusões NÃO mudaram (o E6.5 fica
  fora deles até ele decidir).

## 6. Perguntas ao professor (decisões, não suposições)

1. Cisalhamento: manter γ = 0,10 (o que rodou) ou refazer com 0,20/0,30 (o comentário sugeriu 0,20)?
2. Descarregar antes do recozimento (análogo do E6.4: reverter até Pxy = 0), para casar com o quadro
   `stored_total` (estado descarregado, só plástico remanescente)? Hoje o estado fica carregado (−33 kbar).
3. NVT (47–53 kbar, sem expansão térmica) ou NPT a P = 0?
4. Critério de "recristalizado": excesso de CNA com baseline (o que fiz) + PE, ou PTM/OVITO nos dumps?
5. Série: 800/850 K no desenho E6.5? A versão liga (Jelinek MEAM) na célula de 62 500 é impraticável
   nesta caixa (~1,6 s/passo ⇒ ~9 dias por 1 ns; o Mendelev acelera na GPU, o MEAM não).

## 7. Série executada a pedido do professor (26/10, somente GPU)

O professor autorizou: "utilize somente a placa gráfica do outro micro é mais eficiente. Prossiga!"
Rodou em cadeia (`ser_e65.sh`, tudo com `-sf gpu`, ordem: controles depois corridas):

| corrida | início → fim | relógio |
|---|---|---:|
| controle FCC puro 850 K | 09:36:23 → 09:38:33 | **2:10** |
| controle FCC puro 800 K | 09:38:33 → 09:40:39 | **2:06** |
| E6.5 850 K (1 ns) | 09:40:39 → 10:20:19 | **39:40** |
| E6.5 800 K (1 ns) | 10:20:19 → 10:58:45 | **38:26** |

(Série completa: 09:36 → 10:58, 1 h 22 min, ~40 min por 1 ns — 17% acima dos 33:55 medidos no 900 K.
Nenhum átomos perdido; nenhuma interação.)

### Resultados da série (excesso de CNA = bruto − baseline do puro na MESMA temperatura)

| T | f"5" bruto | baseline puro | **excesso** | PE cisalhado | PE puro | **ΔPE** | densidade |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 800 K | 37,6% | 17,7% | **20,0%** | −3,2765 | −3,3018 | **+0,0253 eV/átomo** | 2,45×10⁸ J/m³ |
| 850 K | 44,2% | 22,8% | **21,5%** | −3,2683 | −3,2940 | **+0,0257 eV/átomo** | 2,49×10⁸ J/m³ |
| 900 K | 50,4% | 28,6% | **21,9%** | −3,2600 | −3,2864 | **+0,0265 eV/átomo** | 2,56×10⁸ J/m³ |

* Pxy residual no hold: −34,9 / −34,0 / −33,1 kbar (800/850/900) — ~constante, não desliza.
* **O armazenado remanescente é praticamente independente da temperatura** no desenho E6.5
  (ΔPE ≈ +0,026 eV/átomo nas três T) e NÃO se libera em 1 ns nas três — estado meta-estável em toda a faixa.
* O excesso de CNA ~21% é estável e o bruto sobe com T só porque o artefato térmico sobe (17,7→28,6%).
* Contraste com os dados antigos (célula 16 384, pré-deformada ε=0,30, DESCARREGADA, NPT): decaimento
  violento em ~2 ps a 800 K. O desenho novo (γ=0,10, carregado, NVT) não decai. A diferença candidata
  a separar os regimes: magnitude do armazenado (0,30 descarregado ⇒ muito mais energia) e/ou o
  barostato. Ficam para a decisão do professor: (b) descarregar até Pxy=0? (c) NVT ou NPT?
* Figuras: `e65_analise/fig_E65_serie.png` (3 painéis das três T), `e65_analise_T{800,850}/fig_E65_T{800,850}.png`,
  CSVs ao lado. Controles: `ctrl_fcc_T{800,850}_{300,hold0,late}.dump/.log` (laptop e caixa).