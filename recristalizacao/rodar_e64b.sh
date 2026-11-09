#!/bin/bash
# E6.4b: descarregamento dos estados POS-escoamento (eps=0.30 limpos, da caixa GPU)
set -u
cd "/c/MacbookPro/ufpa/PPGEM/SDAS Al_Si_Mg_Fe/modelos_elasto-plastos/recristalizacao" || exit 1
LMP="/c/Users/ileao/AppData/Local/LAMMPS 64-bit 10Sep2025/bin/lmp"
run() {
  local log="$1"; shift
  echo "STAGE inicio $(date +%H:%M) $log" | tee -a noite_console.log
  "$LMP" -sf omp -pk omp 8 -log "$log" "$@" >> noite_console.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
    echo "STAGE FALHOU $(date +%H:%M) rc=$rc $log" | tee -a noite_console.log
    exit 1
  fi
  echo "STAGE ok $(date +%H:%M) $log" | tee -a noite_console.log
}
run descarga_jelinek_0.30.log -in in.descarregar.lmp -var pot jelinek -var alloy 1 -var datafile data_jelinek_0.30.data
run descarga_mendelev_0.30.log -in in.descarregar.lmp -var pot mendelev -var alloy 0 -var datafile data_mendelev_0.30.data
echo "STAGE fim $(date +%H:%M) E6.4b completa" | tee -a noite_console.log
