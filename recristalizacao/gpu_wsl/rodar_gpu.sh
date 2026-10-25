#!/bin/bash
# E6-GPU — as corridas pesadas da recristalização, para a máquina com NVIDIA:
#   G1 E6.1 grande (16 384 átomos, Jelinek MEAM) -> produz a base deformada data_jelinek_0.14.data
#   G2 série de recozimento na célula grande (500-900 K)
#   G3 fronteira fina (750/775/800 K) x 3 sementes térmicas
#   G4 recozimento longo: 1 ns em 750 K
# MEAM (Jelinek) nao acelera em GPU: roda em MPI/OpenMP. EAM (Mendelev) acelera: -sf gpu.
# Configurar no topo (ou por variável de ambiente). Baixar com:
#   git clone https://github.com/ileaof/gm-elastoplastic-nucleation.git
#   cd gm-elastoplastic-nucleation/recristalizacao
#   bash rodar_gpu.sh            (tudo)         ou:  SECTIONS="G1 G3" bash rodar_gpu.sh
set -u
cd "$(dirname "$0")"
LMP="${LMP:-lmp}"
NP="${NP:-16}"                    # ranks MPI para o MEAM (ajuste aos núcleos da máquina)
NT="${NT:-8}"                     # threads OpenMP por rank
SEEDS="${SEEDS:-11 22 33}"
SECTIONS="${SECTIONS:-G0 G1 G2 G3 G4}"
MPIRUN="$(command -v mpirun || true)"
[ -n "$MPIRUN" ] && MPIRUN="mpirun -np $NP --oversubscribe -bind-to none"

run() { # run <log> <in> <var...>
  local log="$1" in="$2"; shift 2
  echo "GPU STAGE inicio $(date +%H:%M) $log"
  set -- -in "$in" -log "$log" "$@"
  if [ -n "$MPIRUN" ]; then
    $MPIRUN "$LMP" -sf omp -pk omp ${NT} "$@" >> gpu_console.log 2>&1
  else
    "$LMP" -sf omp -pk omp ${NT} "$@" >> gpu_console.log 2>&1
  fi
  local rc=$?
  if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
    echo "GPU STAGE FALHOU $(date +%H:%M) rc=$rc $log"
    exit 1
  fi
  echo "GPU STAGE ok $(date +%H:%M) $log"
}

run_gpu() { # run_gpu <log> <in> <var...>  (Mendelev, EAM acelerado por GPU)
  local log="$1" in="$2"; shift 2
  echo "GPU STAGE inicio $(date +%H:%M) $log"
  "$LMP" -sf gpu -pk gpu 1 -in "$in" -log "$log" "$@" >> gpu_console.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
    echo "GPU STAGE FALHOU $(date +%H:%M) rc=$rc $log"
    exit 1
  fi
  echo "GPU STAGE ok $(date +%H:%M) $log"
}

case " $SECTIONS " in *" G0 "*) # G0 — checagem do ambiente: 1 000 passos MEAM
  run validacao_gpu.log in.bench.lmp -var pot jelinek -var alloy 1 -var nx 4 -var nsteps 500;; esac

case " $SECTIONS " in *" G1 "*) # G1 — célula grande deformada (Jelinek, ~1-3 h)
  run deform_g.log in.deformar.lmp -var pot jelinek -var alloy 1 -var nx 16;; esac

case " $SECTIONS " in *" G2 "*) # G2 — recozimento na célula grande (~1,5-3 h por patamar)
  for T in 500 600 700 800 850 900; do
    run recozer_jelinek_T$T.log in.recozer.lmp -var pot jelinek -var alloy 1 \
        -var datafile data_jelinek_0.14.data -var T $T -var dumpfreq 5000
  done;; esac

case " $SECTIONS " in *" G3 "*) # G3 — fronteira fina x sementes térmicas (~45-90 min cada)
  for T in 750 775 800; do
    for S in $SEEDS; do
      run recozer_jelinek_T${T}_s$S.log in.recozer.lmp -var pot jelinek -var alloy 1 \
          -var datafile data_jelinek_0.14.data -var T $T -var seed $S \
          -var tag _s$S -var dumpfreq 5000
    done
  done;; esac

case " $SECTIONS " in *" G4 "*) # G4 — recozimento longo: 1 ns na fronteira (~10-20 h)
  run recozer_jelinek_longo_T750.log in.recozer.lmp -var pot jelinek -var alloy 1 \
      -var datafile data_jelinek_0.14.data -var T 750 \
      -var nrampe 5000 -var nprod 500000 -var tag _longo -var dumpfreq 10000;; esac

echo "GPU fim $(date +%H:%M): rode bash juntar_resultados.sh para levar os resultados" | tee -a gpu_console.log