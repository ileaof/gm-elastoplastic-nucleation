#!/bin/bash
# Cadeia pesada do E6 na caixa GPU (NOTE1, RTX 3050): G3 -> G4 -> G1.
# Lançamento:  nohup bash rodar_cadeia_gpu.sh > cadeia_nohup.log 2>&1 &
# Progresso:   progresso.txt (linhas STAGE); console: cadeia_console.log
set -u
cd "$(dirname "$0")"
LMP="$HOME/lammps-gpu/bin/lmp"
export OMP_NUM_THREADS=8
P=progresso.txt

echo "CADEIA inicio $(date)" >> "$P"

# base do estado deformado (do clone compartilhado, se faltar)
cp -n /mnt/c/Users/labmo_agmg0g1/OneDrive/Documentos/gm-elastoplastic-nucleation/recristalizacao/data_jelinek_0.14.data . || true
[ -f data_jelinek_0.14.data ] || { echo "CADEIA FALHOU: base data_jelinek_0.14.data ausente $(date)" >> "$P"; exit 1; }

# G3 — fronteira fina: 750/775/800 K × 3 sementes térmicas (9 recozimentos, dumps espaçados)
for S in 11 22 33; do
  for T in 750 775 800; do
    log="recozer_jelinek_T${T}_s${S}.log"
    echo "STAGE inicio $(date +%H:%M) anneal T=$T s=$S" >> "$P"
    "$LMP" -sf omp -pk omp 8 -log "$log" -in in.recozer.lmp -var pot jelinek -var alloy 1 \
           -var datafile data_jelinek_0.14.data -var T $T \
           -var seed $((S * 1000 + T)) -var tag _s${S} -var dumpfreq 5000 >> cadeia_console.log 2>&1
    rc=$?
    if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
      echo "STAGE FALHOU rc=$rc anneal T=$T s=$S" >> "$P"
    else
      echo "STAGE ok $(date +%H:%M) anneal T=$T s=$S" >> "$P"
    fi
  done
done

# G4 — recozimento longo: 1 ns em 750 K (a 5ª semente térmica)
log="recozer_jelinek_longo_T750.log"
echo "STAGE inicio $(date +%H:%M) longo 750 K" >> "$P"
"$LMP" -sf omp -pk omp 8 -log "$log" -in in.recozer.lmp -var pot jelinek -var alloy 1 \
       -var datafile data_jelinek_0.14.data -var T 750 -var seed 777 \
       -var nrampe 5000 -var nprod 500000 -var tag _longo -var dumpfreq 10000 >> cadeia_console.log 2>&1
rc=$?
if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
  echo "STAGE FALHOU rc=$rc longo" >> "$P"
else
  echo "STAGE ok $(date +%H:%M) longo" >> "$P"
fi

# G1 — célula grande: 16 384 átomos deformados (base para a série da célula grande, rodar depois)
log="deform_g16.log"
echo "STAGE inicio $(date +%H:%M) deformar nx=16" >> "$P"
"$LMP" -sf omp -pk omp 8 -log "$log" -in in.deformar.lmp -var pot jelinek -var alloy 1 -var nx 16 >> cadeia_console.log 2>&1
rc=$?
if [ $rc -ne 0 ] || grep -q "ERROR:" "$log" 2>/dev/null; then
  echo "STAGE FALHOU rc=$rc nx16" >> "$P"
else
  echo "STAGE ok $(date +%H:%M) nx16" >> "$P"
fi

echo "CADEIA fim $(date)" >> "$P"