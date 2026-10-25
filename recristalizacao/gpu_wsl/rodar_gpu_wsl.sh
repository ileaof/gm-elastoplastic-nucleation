#!/bin/bash

# ==============================================================================
# Script para rodar Recristalização (E6) via GPU NVIDIA
# Baseado no GPU.md e PLANO.md
# ==============================================================================

# --- Configurações de Hardware ---
# LMP: executável com pacote GPU+MEAM+OPENMP (compilado pelo compilar_lammps_gpu_wsl.sh)
LMP="$HOME/lammps-gpu/bin/lmp"
# NP: MPI (o build aqui é serial); NT: threads OpenMP (12 cores do WSL)
NP=1
NT=8

echo "Iniciando simulações E6 via GPU..."
echo "Configuração: LMP=$LMP, NP=$NP, NT=$NT"

# Jelinek (MEAM): sem versão GPU no LAMMPS — roda no CPU via meam/omp (GPU.md);
# Mendelev (EAM): acelera de verdade na RTX 3050 com -sf gpu.
# -v é usado em vez de -var (forma atual do LAMMPS); -var ainda funciona e é mantido.

# 1. DEFORMAÇÃO E MEDIÇÃO (E6.1)
# Parametros: pot (jelinek/mendelev), nx (escala), eps_list (deformações)
echo "--- Executando E6.1: Deformação e Medição ---"
$LMP -sf omp -pk omp $NT -log deform_jelinek.log -in in.deformar.lmp -var pot jelinek -var nx 8 -var eps_list 0.05 0.10 0.14 0.20 0.30 || exit 1

echo "--- Executando E6.1 (controle Mendelev, GPU) ---"
$LMP -sf gpu -pk gpu 1 -log deform_mendelev.log -in in.deformar.lmp -var pot mendelev -var nx 8 -var eps_list 0.05 0.10 0.14 0.20 0.30 || exit 1

# 2. RECOZIMENTO (E6.2)
# Parametros: data (arquivo de entrada), T (temperatura)
echo "--- Executando E6.2: Recozimento ---"
# Exemplo para epsilon 0.14 e T = 800K
$LMP -sf omp -pk omp $NT -log recozer_jelinek_T800.log -in in.recozer.lmp -var pot jelinek -var datafile data_jelinek_0.14.data -var T 800 || exit 1

echo "--- Executando E6.2 (controle Mendelev, GPU) ---"
$LMP -sf gpu -pk gpu 1 -log recozer_mendelev_T800.log -in in.recozer.lmp -var pot mendelev -var datafile data_mendelev_0.14.data -var T 800 || exit 1

echo "Simulações concluídas com sucesso na GPU."
echo "Agora execute: bash juntar_resultados.sh"
