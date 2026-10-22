#!/usr/bin/env bash
#SBATCH --job-name=dds-cure90
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768M
#SBATCH --time=24:00:00
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
PILOT="$ROOT/runs/pilot-20261021T164259"
LMP="$ROOT/lammps-dds/bin/lmp"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
date -u
hostname
cp "$ROOT"/inputs/* .
python3 "$ROOT/scripts/export_topology.py" "$PILOT/pilot.data" data.start
# A bounded smoke tests both reactive cycles and the intervening anneal before
# the long run. It uses a separate directory and never changes production input.
mkdir validation
cp data.start *Top.txt map*.txt in.cure_campaign validation/
(
 cd validation
 "$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.cure_campaign \
   -var chunk_steps 100 -var blocks 1 -var cycles 2 \
   -var ramp_steps 100 -var anneal_steps 100 -var postanneal_steps 100 \
   -log log.lammps
)
echo CAMPAIGN_VALIDATION_PASSED
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.cure_campaign -log log.lammps
python3 "$ROOT/scripts/analyze_campaign.py" .
date -u
