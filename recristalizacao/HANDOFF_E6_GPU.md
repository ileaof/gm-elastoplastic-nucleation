# E6 na máquina GPU — procedimento de continuidade (para o operador/agente que assume lá)

Objetivo (prof. Ivaldo Leão, 24-25/10/2026): provar por MD que a deformação plástica armazenada da
superfície elasto-plástica (`stored_total`, energia σ_y·α) reaparece como potencial na recristalização:
aquecendo o sólido deformado até ativar a difusão, a recristalização inicia e libera essa energia.
O laptop principal (sem NVIDIA) já fez os pilotos e a série pequena; o que falta é pesado e roda aqui na
máquina GPU. Contexto completo: `PLANO.md` (a física, a cronometragem, o registro da noite de 25/10) e
`GPU.md` (tabela das seções e tempos).

## 1. O que rodar (na máquina GPU)

    git clone https://github.com/ileaof/gm-elastoplastic-nucleation.git
    cd gm-elastoplastic-nucleation/recristalizacao
    bash rodar_gpu.sh                   # a sequência inteira (G0..G4), uma vez

Configurações no topo do `rodar_gpu.sh` (ou por variável de ambiente):
`LMP` (caminho do lmp), `NP` (ranks MPI), `NT` (threads OpenMP), `SEEDS` (sementes térmicas do G3),
`SECTIONS` (subconjunto, ex.: `SECTIONS="G0 G1 G3"`).

| Seção | Corrida | Tempo grosso |
|---|---|---:|
| G0 | validação (500 passos, MEAM) | ~1 min |
| G1 | E6.1 grande: 16 384 átomos da liga deformados (Jelinek) | ~1–3 h |
| G2 | recozimento da célula grande: 500–900 K | ~10–16 h |
| G3 | **fronteira fina**: 750/775/800 K × 3 sementes | ~4–8 h |
| G4 | **recoximento longo**: 1 ns em 750 K | ~10–20 h |

Detalhes de física: MEAM (Jelinek/liga real) não acelera em GPU — MPI+OpenMP; Mendelev (EAM Al puro,
controle) aceita `-sf gpu`. Os eventos violentos ocorrem nos primeiros 10 ps — é o G3 que localiza o
limiar com estatística. Se alguma seção der `Lost atoms`/ERROR, o `rodar_gpu.sh` aborta e a seção pode
ser refeita depois de revisada a configuração (veja o registro da noite em `PLANO.md` §8: em 2 048
átomos a fronteira 700–800 K ejeta átomos mesmo com subida blindada — no G3, célula maior, o regime
pode ser estável; não edite a física das entradas sem o professor).

## 2. O que mandar de volta para o laptop principal

Na máquina GPU, no fim:

    bash juntar_resultados.sh      # cria results_e6.tar.gz com os arquivos pequenos

Depois, uma das rotas:

* **Arquivo (sempre funciona):** enviar o `results_e6.tar.gz` ao laptop (OneDrive, transferência ou USB)
  e descompactar dentro de
  `C:\MacbookPro\ufpa\PPGEM\SDAS Al_Si_Mg_Fe\modelos_elasto-plastos\recristalizacao\`
  (mantendo a pasta `results_e6/` ou extraindo direto — o agente local organiza).
* **Git (se a máquina tiver credenciais de push):**
  `git add . && git commit -m "E6 GPU results" && git push` — o agente do laptop faz pull e organiza.
* **Opcional (só se for renderizar no OVITO):** os `recozido_*.dump` de interesse e o
  `deformed_jelinek.dump` — são grandes; a rota de arquivo cobre.

## 3. O que NÃO precisa mandar

* `data_*.data` (o G1 recria as bases dele localmente);
* os logs brutos completos do LAMMPS dentro do tar (o `juntar_resultados.sh` já filtra o essencial).

## 4. No laptop, depois dos resultados chegarem

Rodar `python analyze_recristalizacao.py` na pasta e reportar ao professor: E_st(ε) na célula grande, a
espessura δ medida, o limiar de recristalização (estatístico, das 3 sementes do G3) e a comparação com o
potencial previsto σ_y·α/δ. Perguntas pendentes do professor para a discussão: (a) novo desenho para
800–900 K (dt menor, barostato, célula); (b) termo térmico (Γ) no recozimento; (c) δ medida ou imposta;
(d) critério final de "recristalizado" (CNA + queda de PE: o adotado nos pilotos).