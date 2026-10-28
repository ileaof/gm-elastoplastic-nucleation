#!/bin/bash

# ==============================================================================
# Script para rodar Recristalização (E6) via CPU (OpenMP)
# Baseado no PLANO.md e GPU.md
# ==============================================================================

# --- Configurações de Hardware ---
# NP: Número de processos MPI (manter 1 para rodar puramente com OpenMP no laptop)
NP=1
# NT: Número de threads OpenMP (ajuste conforme seu processador, ex: 8 ou 12)
NT=16
# Comando do LAMMPS (se 'lmp' não estiver no PATH, coloque o caminho completo aqui)
LMP=lmp

echo "Iniciando simulações E6 via CPU..."
echo "Configuração: NP=$NP, NT=$NT"

# 1. DEFORMAÇÃO E MEDIÇÃO (E6.1)
# Parâmetros:
#   pot: jelinek (liga) ou mendelev (Al puro)
#   nx: escala da célula (8 = piloto, 16 = grande)
#   eps_list: lista de deformações para medir a energia armazenada
echo "--- Executando E6.1: Deformação e Medição ---"
$LMP -sf omp -pk omp $NT -in in.deformar.lmp -var pot jelinek -var nx 8 -var eps_list 0.05 0.10 0.14 0.20 0.30

# 2. RECOZIMENTO (E6.2)
# Parâmetros:
#   data: o arquivo de dados gerado no passo anterior (ex: data_0.14.data para epsilon=0.14)
#   T: Temperatura de recozimento (ex: 800 K)
echo "--- Executando E6.2: Recozimento ---"
# Exemplo para a deformação de 0.14 e temperatura de 800K
$LMP -sf omp -pk omp $NT -in in.recozer.lmp -var pot jelinek -var data data_0.14.data -var T 800

echo "Simulações concluídas. Agora você pode rodar o 'bash juntar_resultados.sh' para coletar os dados."
