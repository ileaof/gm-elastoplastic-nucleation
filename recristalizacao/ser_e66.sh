#!/bin/bash
# E6.6 — cadeia autorizada pelo professor (26/10, "faça, deixe em background"):
# item 5 (controles NVT frios 500/600/700 K) -> controles NPT 800/850/900 K (a régua)
# -> item 2 (γ=0.10 descarregado + NPT, 2 sementes) -> item 2b (γ=0.30 descarregado + NPT)
# -> item 3 (γ=0.20 e 0.30 carregados, NVT). Tudo somente na placa (-sf gpu).
# Rodada com nohup: sobrevive a desconexões do ssh.
cd ~/e6_rec
LMP=$HOME/lammps-gpu/bin/lmp
echo "E6.6 inicio $(date) - somente GPU"

stage () {   # stage <tag> <args...>
  local tag="$1"; shift
  echo "STAGE inicio $tag $(date)"
  "$@" > nohup_$tag.log 2>&1
  if [ $? -eq 0 ]; then touch ok_$tag.flag; echo "STAGE ok $tag $(date)"
  else echo "STAGE FALHOU $tag $(date)"; fi
}

# ---- item 5: controles NVT frios (bem rápidos: ~2 min; usa os dumps de 1 ns? não: o ctrl_fcc roda 30k passos)
for T in 500 600 700; do
  stage ctrl_T$T $LMP -sf gpu -var T2 $T -in in.ctrl_fcc.lmp -log ctrl_fcc_T$T.log
done

# ---- régua NPT dos controles puros (1 ns + dumps)
for T in 800 850 900; do
  stage ctrlnpt_T$T $LMP -sf gpu -var T2 $T -in in.ctrl_npt.lmp -log ctrl_npt_T$T.log
done

# ---- item 2: γ = 0.10 descarregado + NPT, 2 sementes (a medição do resto plástico + o teste P=0)
for T in 800 850 900; do for S in 101 202; do
  stage descarnpt_T${T}_s$S $LMP -sf gpu -var T2 $T -var seed $S -var gam 0.10 \
        -in in.descar_npt.lmp -log desca_npt_T${T}_s$S.log
done; done

# ---- item 2b: γ = 0.30 descarregado + NPT (o estado que decaiu nas Figs 8/10, na célula grande)
for T in 800 850 900; do
  stage desca30_T$T $LMP -sf gpu -var T2 $T -var seed 101 -var gam 0.30 \
        -in in.descar_npt.lmp -log desca_npt30_T$T.log
done

# ---- item 3: γ = 0.20 e 0.30 carregados + NVT (a variação em γ do desenho E6.5)
for T in 800 850 900; do
  stage largo20_T$T $LMP -sf gpu -var T2 $T -var gam 0.20 -in in.rec_largo.lmp -log rec_largo20_T$T.log
done
for T in 800 850 900; do
  stage largo30_T$T $LMP -sf gpu -var T2 $T -var gam 0.30 -in in.rec_largo.lmp -log rec_largo30_T$T.log
done

echo "E6.6_COMPLETA $(date)"
df -h ~ | tail -1