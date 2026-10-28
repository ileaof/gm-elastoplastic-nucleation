#!/bin/bash
# Compila o LAMMPS (stable_30Sep2026) dentro do WSL Ubuntu 22.04 com os pacotes
# que o E6 precisa e que faltam no lmp do apt: MEAM (Jelinek), OPENMP (meam/omp) e
# GPU (CUDA — a RTX 3050 acelera o Mendelev/EAM; o MEAM roda em OpenMP).
# Sem sudo: tudo em $HOME (fonte em ~/src, build em ~/build, binário em ~/lammps-gpu).
# Uso:  bash compilar_lammps_gpu_wsl.sh   (dentro do WSL)
set -eu

TAG=stable_30Sep2026
export CUDA_HOME=/usr/local/cuda-13.2
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

SRC=$HOME/src/lammps-$TAG
BUILD=$HOME/build/lammps-$TAG

mkdir -p ~/src ~/build
cd ~/src
if [ ! -d "$SRC" ]; then
  echo "== baixando $TAG ..."
  curl -sL "https://github.com/lammps/lammps/archive/refs/tags/$TAG.tar.gz" -o lammps-$TAG.tar.gz
  tar xzf lammps-$TAG.tar.gz
fi

mkdir -p "$BUILD"
cd "$BUILD"
cmake "$SRC/cmake" \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=no \
  -D PKG_MEAM=yes \
  -D PKG_MANYBODY=yes \
  -D PKG_OPENMP=yes \
  -D PKG_GPU=yes \
  -D GPU_API=cuda \
  -D GPU_ARCH=86 \
  -D CMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" \
  -D CMAKE_INSTALL_PREFIX="$HOME/lammps-gpu"
make -j6
make install

echo "== pronto! =="
"$HOME/lammps-gpu/bin/lmp" -h 2>/dev/null | grep -A6 'Installed packages'