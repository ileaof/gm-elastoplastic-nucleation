# E6 na máquina GPU — procedimento completo

## Baixar

    git clone https://github.com/ileaof/gm-elastoplastic-nucleation.git
    cd gm-elastoplastic-nucleation/recristalizacao

## Rodar (só uma vez, sequência inteira: ~1-2 dias)

    bash rodar_gpu.sh

Por seção (editável no topo do `rodar_gpu.sh`):

    LMP=/caminho/do/lmp NP=16 NT=8 SECTIONS="G0 G1 G3" bash rodar_gpu.sh

| Seção | O que faz | Tempo grosso (16 ranks) |
|---|---|---:|
| G0 | valida o ambiente (500 passos, MEAM) | ~1 min |
| G1 | E6.1 **grande**: 16 384 átomos deformados (Jelinek, liga) | ~1-3 h |
| G2 | recozimento da célula grande: 500-900 K (6 patamares) | ~10-16 h |
| G3 | **fronteira fina**: 750 / 775 / 800 K × 3 sementes térmicas | ~4-8 h |
| G4 | **recoximento longo**: 1 ns em 750 K | ~10-20 h |

Notas: MEAM (Jelinek) não tem versão GPU/omp-acelerada no LAMMPS — roda em MPI+OpenMP
(`NP`/`NT`); Mendelev (EAM) acelera em GPU (`-sf gpu`) e é o controle. Os eventos violentos
acontecem nos primeiros 10 ps — o G3 é que localiza o limiar com estatística.

## Trazer de volta ao laptop

    bash juntar_resultados.sh      # cria results_e6.tar.gz com logs/termos/perfis/figuras

Depois, qualquer um dos dois caminhos:

1. **Git** (se a máquina tiver as credenciais): `git add results_e6 && git commit -m "E6 GPU results" && git push` — no laptop: `git pull`.
2. **Arquivo**: copie o `results_e6.tar.gz` (OneDrive/pendrive/e-mail) e descompacte em
   `modelos_elasto-plastos/recristalizacao/` — os dumps ficarão na GPU; se o OVITO for renderí-los,
   copie também os `recozido_*.dump` de interesse.

De volta ao laptop, o agente roda `analyze_recristalizacao.py` sobre os resultados completos.