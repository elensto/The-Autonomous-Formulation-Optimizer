#!/usr/bin/env bash
#SBATCH --job-name=dds-validate
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768M
#SBATCH --time=00:10:00
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
LMP="$ROOT/lammps-dds/bin/lmp"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
cp "$ROOT"/inputs/* .
python3 "$ROOT/scripts/export_topology.py" "$ROOT/runs/pilot-20261021T164259/pilot.data" data.start
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.cure_campaign \
 -var chunk_steps 100 -var blocks 1 -var cycles 2 \
 -var ramp_steps 100 -var anneal_steps 100 -var postanneal_steps 100 -log log.cure.lammps
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.equilibrate -var equil_steps 100 -log log.equilibrium.lammps
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.tg_screen \
 -var heat_steps 100 -var hold_steps 100 -var plateau_equil_steps 100 -var plateau_prod_steps 100 -log log.tg.lammps
echo PIPELINE_SYNTAX_AND_EXECUTION_VALIDATED_NOT_PHYSICAL_PRODUCTION
