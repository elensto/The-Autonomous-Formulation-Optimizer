#!/usr/bin/env bash
#SBATCH --job-name=dds-properties
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768M
#SBATCH --time=12:00:00
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
CURE_DIR="$1"
LMP="$ROOT/lammps-dds/bin/lmp"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
date -u
hostname
if ! python3 -c 'import json,sys;sys.exit(0 if json.load(open(sys.argv[1]))["structure_gate_passed"] else 1)' "$CURE_DIR/campaign_summary.json"; then
 echo "PROPERTIES_SKIPPED: cure target/connectivity gate not passed" | tee properties_status.txt
 exit 0
fi
cp "$ROOT/inputs/in.equilibrate" "$ROOT/inputs/in.tg_screen" .
cp "$CURE_DIR/final.cure.data" .
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.equilibrate -log log.equilibrium.lammps
if ! python3 "$ROOT/scripts/analyze_properties.py" equilibrium .; then
 echo "TG_SKIPPED: NPT stationarity gate not passed" | tee properties_status.txt
 exit 0
fi
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.tg_screen -log log.tg.lammps
python3 "$ROOT/scripts/analyze_properties.py" tg .
echo "PRELIMINARY_PROPERTIES_PIPELINE_COMPLETED" | tee properties_status.txt
date -u
