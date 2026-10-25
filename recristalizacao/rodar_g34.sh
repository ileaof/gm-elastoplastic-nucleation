#!/bin/bash
# G3 + G4 refazidos (com o tag corrigido) enquanto o nx=16 segue
set -u
cd ~/e6_rec
LMP=$HOME/lammps-gpu/bin/lmp
export OMP_NUM_THREADS=8
P=progresso.txt
for S in 11 22 33; do
  for T in 750 775 800; do
    log="recozer_jelinek_T${T}_s${S}.log"
    echo "STAGE inicio $(date +%H:%M) anneal T=$T s=$S" >> $P
    $LMP -sf omp -pk omp 8 -log $log -in in.recozer.lmp -in in.recozer.lmp -var pot jelinek -var alloy 1          -var datafile data_jelinek_0.14.data -var T $T          -var seed $((S * 1000 + T)) -var tag _s${S} -var dumpfreq 5000 >> cadeia_console.log 2>&1
    rc=$?
    if [ $rc -ne 0 ] || grep -q 'ERROR:' $log 2>/dev/null; then
      echo "STAGE FALHOU rc=$rc anneal T=$T s=$S" >> $P
    else
      echo "STAGE ok $(date +%H:%M) anneal T=$T s=$S" >> $P
    fi
  done
done
log=longo_T750_5.log
echo "STAGE inicio $(date +%H:%M) longo 750" >> $P
$LMP -sf omp -pk omp 8 -log $log -in in.recozer.lmp -var pot jelinek -var alloy 1      -var datafile data_jelinek_0.14.data -var T 750 -var seed 777      -var nrampe 5000 -var nprod 500000 -var tag _longo -var dumpfreq 10000 >> cadeia_console.log 2>&1
rc=$?
if [ $rc -ne 0 ] || grep -q 'ERROR:' $log 2>/dev/null; then echo "STAGE FALHOU rc=$rc longo" >> $P; else echo "STAGE ok $(date +%H:%M) longo" >> $P; fi
echo "CADEIA2 fim $(date)" >> $P
