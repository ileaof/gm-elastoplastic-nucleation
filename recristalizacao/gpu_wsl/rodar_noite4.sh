#!/bin/bash
# Noite do E6, parte 4 (25/10/2026, 04:03): série Mendelev completa (500-900 K, protocolo endurecido).
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

for T in 500 600 700 800 850 900; do
  run recozer_mendelev_T$T.log -in in.recozer.lmp -var pot mendelev -var alloy 0 \
      -var datafile data_mendelev_0.14.data -var T $T
done
echo "STAGE fim $(date +%H:%M) noite 4 completa" | tee -a noite_console.log