#!/bin/bash
# Mendelev EAM na RTX 3050 de verdade (-sf gpu): a série de recozimentos 500-900 K
# e o descarregamento do estado 0.30. Rápido no cartão (minutos).
set -u
cd ~/e6_rec
LMP="$HOME/lammps-gpu/bin/lmp"
P=progresso.txt
export OMP_NUM_THREADS=8
cp -n /mnt/c/Users/labmo_agmg0g1/OneDrive/Documentos/gm-elastoplastic-nucleation/recristalizacao/data_mendelev_0.30.data . || true

for T in 500 600 700 800 850 900; do
  log="gpu_mendelev_T${T}.log"
  echo "STAGE inicio $(date +%H:%M) mendelev_gpu T=$T" >> "$P"
  "$LMP" -sf gpu -pk gpu 1 -in in.recozer.lmp -var pot mendelev -var alloy 0 \
         -var datafile data_mendelev_0.14.data -var T $T -var seed $((T + 7)) -var dumpfreq 5000 \
         -log "$log" >> cadeia_console.log 2>&1
  rc=$?
  if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
    echo "STAGE FALHOU rc=$rc mendelev_gpu T=$T" >> "$P"
  else
    echo "STAGE ok $(date +%H:%M) mendelev_gpu T=$T" >> "$P"
  fi
done

log="gpu_descarga_mendelev_0.30.log"
echo "STAGE inicio $(date +%H:%M) descarga mend 0.30 gpu" >> "$P"
"$LMP" -sf gpu -pk gpu 1 -in in.descarregar.lmp -var pot mendelev -var alloy 0 \
       -var datafile data_mendelev_0.30.data -log "$log" >> cadeia_console.log 2>&1
rc=$?
if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
  echo "STAGE FALHOU rc=$rc descarga_gpu" >> "$P"
else
  echo "STAGE ok $(date +%H:%M) descarga_gpu" >> "$P"
fi
echo "GPU-MENDELEV fim $(date)" >> "$P"