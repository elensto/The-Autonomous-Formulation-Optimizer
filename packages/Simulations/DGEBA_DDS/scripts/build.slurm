#!/usr/bin/env bash
#SBATCH --job-name=dds-build
#SBATCH --partition=grafite
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=reports/build-%j.out
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
REFERENCE=$HOME/formulationlab-upet-20261020
# Reuse already downloaded LAMMPS sources and compiler, with an isolated build.
source "$REFERENCE/activate.sh"
cd "$ROOT"
cmake -S "$REFERENCE/lammps-metatomic/cmake" -B build-lammps -G Ninja \
 -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
 -DCMAKE_INSTALL_PREFIX="$ROOT/lammps-dds" \
 -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DLAMMPS_INSTALL_RPATH=ON \
 -DBUILD_MPI=OFF -DBUILD_OMP=ON -DPKG_OPENMP=ON \
 -DPKG_MOLECULE=ON -DPKG_KSPACE=ON -DPKG_REACTION=ON -DFFT=KISS
cmake --build build-lammps --parallel "$SLURM_CPUS_PER_TASK"
cmake --install build-lammps
"$ROOT/lammps-dds/bin/lmp" -help > reports/lammps-help.txt
grep -q 'bond/react' reports/lammps-help.txt
touch reports/BUILD_COMPLETED
