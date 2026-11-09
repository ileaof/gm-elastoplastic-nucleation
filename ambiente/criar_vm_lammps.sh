#!/bin/bash
# Creates the Multipass VM "lammps" (Ubuntu 22.04), mounts this repository inside it
# and compiles LAMMPS there with ambiente/compilar_lammps_vm.sh.
#
# Run it ON THE MAC only when no heavy computation is running: the build uses every vCPU of the VM.
# Usage:  bash ambiente/criar_vm_lammps.sh  [CPUS]  [MEMORIA]
set -eu
CPUS=${1:-12}
MEM=${2:-16G}
REPO="$(cd "$(dirname "$0")/.." && pwd)"

multipass launch 22.04 --name lammps --cpus "$CPUS" --memory "$MEM" --disk 20G
multipass mount "$REPO" lammps:/home/ubuntu/gm-elastoplastic-nucleation
multipass exec lammps -- bash /home/ubuntu/gm-elastoplastic-nucleation/ambiente/compilar_lammps_vm.sh
echo "Entrar na VM:  multipass shell lammps   (repo em ~/gm-elastoplastic-nucleation)"
