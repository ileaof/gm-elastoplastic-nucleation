#!/bin/bash
# E6.4 nesta máquina (WSL + LAMMPS compilado): descarregamento elástico nos estados saudáveis.
# Jelinek (MEAM, CPU/omp); Mendelev (EAM) na GPU. Igual ao rodar_e64.sh remoto, caminhos daqui.
set -u
cd "$(dirname "$0")" || exit 1
export PATH="/home/labmo/lammps-gpu/bin:/usr/local/cuda-13.2/bin:$PATH"
export LD_LIBRARY_PATH="/usr/local/cuda-13.2/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=8
LMP=lmp

run() {
  local log="$1"; shift
  echo "E64 STAGE inicio $(date +%H:%M) $log"
  "$LMP" -log "$log" "$@" > /dev/null 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
    echo "E64 STAGE FALHOU rc=$rc $log"
  else
    echo "E64 STAGE ok $log"
  fi
}

run descarga_jelinek_0.20.log   -sf omp -pk omp 8  -in in.descarregar.lmp -var pot jelinek   -var datafile data_jelinek_0.20.data
run descarga_mendelev_0.14.log  -sf gpu -pk gpu 1  -in in.descarregar.lmp -var pot mendelev  -var datafile data_mendelev_0.14.data
run descarga_mendelev_0.20.log  -sf gpu -pk gpu 1  -in in.descarregar.lmp -var pot mendelev  -var datafile data_mendelev_0.20.data
echo "E64 completa"