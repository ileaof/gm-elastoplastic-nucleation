#!/bin/bash
cd ~/e6_rec
LMP=$HOME/lammps-gpu/bin/lmp
echo "SERIE E6.5 inicio $(date) - somente GPU (-sf gpu)"
for T in 850 800; do
  echo "ctrl T=$T inicio $(date)"
  $LMP -sf gpu -var T2 $T -in in.ctrl_fcc.lmp -log ctrl_fcc_T$T.log > nohup_ctrl_T$T.log 2>&1
  if [ $? -eq 0 ]; then touch ctrl_$T.ok; echo "ctrl T=$T OK $(date)"; else echo "ctrl T=$T FALHOU $(date)"; fi
done
for T in 850 800; do
  echo "rec T=$T inicio $(date)"
  $LMP -sf gpu -var T2 $T -in in.rec_gpu.lmp -log rec_gpu_T$T.log > nohup_rec_T$T.log 2>&1
  if [ $? -eq 0 ]; then touch rec_$T.ok; echo "rec T=$T OK $(date)"; else echo "rec T=$T FALHOU $(date)"; fi
done
echo "SERIE_E65_COMPLETA $(date)"
