#!/usr/bin/env bash
#SBATCH --job-name=dds-formulations
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768M
#SBATCH --time=1-12:00:00
set -euo pipefail
BATCH="$1"
test -f "$BATCH/VALIDATION_PASSED"
CASE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BATCH/case_ids.txt")
test -n "$CASE"
bash "$BATCH/common/parallel_runner.sh" "$BATCH" "$CASE" production
