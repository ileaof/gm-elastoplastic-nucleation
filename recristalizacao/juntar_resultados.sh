#!/bin/bash
# E6 — coleta dos resultados na máquina GPU para trazer de volta ao laptop:
#   bash juntar_resultados.sh
# Cria results_e6.tar.gz com os arquivos pequenos (logs, termos, perfis, figuras) — os
# dumps ficam no lugar (grandes); se o OVITO for usar, copiem os *.dump à parte.
set -u
cd "$(dirname "$0")"
mkdir -p results_e6
cp -v stored_*.dat perfil_Es_*.dat recozer_*.log deform_*.log bench_*.log \
      validacao_*.log fig_E6_*.png noite_console.log gpu_console.log 2>/dev/null | tail -1
echo "--- dumps (grandes, ficam fora do tar):"; ls -lh *.dump *.data 2>/dev/null | awk '{print $5, $9}'
tar -czf results_e6.tar.gz results_e6
echo "resultado: $(ls -lh results_e6.tar.gz | awk '{print $5, $9}')"
echo "de volta ao laptop: descompacte em modelos_elasto-plastos/recristalizacao/ e avise o agente"