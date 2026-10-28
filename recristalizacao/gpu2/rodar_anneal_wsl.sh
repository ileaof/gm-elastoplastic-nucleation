#!/bin/bash
# E6.2c — recozimentos blindados nesta máquina (RTX 3050 + LAMMPS WSL compilado).
# Jelinek (MEAM, CPU/omp) e Mendelev (EAM, GPU).
# Uso (dentro do WSL, no directorio recristalizacao):  bash rodar_anneal_wsl.sh [[T ...]]
set -u
export PATH="/home/labmo/lammps-gpu/bin:/usr/local/cuda-13.2/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/cuda-13.2/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=8

T_LIST=${*:-700 800}
for T in $T_LIST; do
  for POT in jelinek mendelev; do
    SF="-sf omp"; PK="-pk omp 8"
    if [ "$POT" = mendelev ]; then SF="-sf gpu"; PK="-pk gpu 1"; fi
    echo "--- E6.2c: $POT T=$T ---"
    lmp $SF $PK -log "recozer_${POT}_T$T.log" -in in.recozer_cap.lmp \
        -var pot $POT -var datafile "data_${POT}_0.14.data" -var T $T
    rc=$?
    if [ $rc -ne 0 ] || grep -q "ERROR:" "recozer_${POT}_T$T.log" 2>/dev/null; then
      echo "STAGE-FALHOU rc=$rc $POT T=$T"
    else
      echo "STAGE-ok $POT T=$T"
    fi
  done
done
echo "recozimentos concluidos"