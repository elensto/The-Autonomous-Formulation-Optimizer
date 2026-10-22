#!/usr/bin/env bash
#SBATCH --job-name=dds-pilot
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
date -u
hostname
"$ROOT/lammps-dds/bin/lmp" -help > lammps-help.txt
grep -q 'bond/react' lammps-help.txt
cp "$ROOT"/inputs/* .
"$ROOT/lammps-dds/bin/lmp" -sf omp -pk omp "$OMP_NUM_THREADS" -in in.pilot -log log.lammps
python3 "$ROOT/scripts/analyze_pilot.py" .
date -u
