#!/bin/bash
# Compiles LAMMPS stable_30Sep2026 (the version of every E6 log) inside the Linux VM.
# Packages: MEAM (Jelinek 2012), MANYBODY (Mendelev 2008 eam/fs), OPENMP.
# No GPU package: Apple Silicon has no CUDA; the former "-sf gpu" runs use "-sf omp -pk omp N".
set -eu
TAG=stable_30Sep2026
SRC=$HOME/src/lammps-$TAG
BUILD=$HOME/build/lammps-$TAG
PREFIX=$HOME/lammps

sudo apt-get update -y
sudo apt-get install -y build-essential cmake git curl python3 python3-venv

mkdir -p ~/src ~/build
if [ ! -d "$SRC" ]; then
  curl -sL "https://github.com/lammps/lammps/archive/refs/tags/$TAG.tar.gz" -o ~/src/lammps-$TAG.tar.gz
  tar xzf ~/src/lammps-$TAG.tar.gz -C ~/src
fi
mkdir -p "$BUILD" && cd "$BUILD"
cmake "$SRC/cmake" \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=no \
  -D PKG_MEAM=yes \
  -D PKG_MANYBODY=yes \
  -D PKG_OPENMP=yes \
  -D CMAKE_INSTALL_PREFIX="$PREFIX"
make -j"$(nproc)"
make install
grep -q 'lammps/bin' ~/.bashrc || echo 'export PATH="$HOME/lammps/bin:$PATH"' >> ~/.bashrc
"$PREFIX/bin/lmp" -h | sed -n '1,3p;/Installed packages/,/^$/p'
echo "== LAMMPS pronto em $PREFIX/bin/lmp =="
