#!/usr/bin/env bash
#SBATCH --job-name=dds-batch-check
#SBATCH --partition=grafeno
#SBATCH --account=users
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=768M
#SBATCH --time=00:20:00
set -euo pipefail
BATCH="$1"
# Exercise deficient, near-stoichiometric and excess-hardener compositions.
for CASE in D100_H040 D100_H049 D087_H050; do
 bash "$BATCH/common/parallel_runner.sh" "$BATCH" "$CASE" smoke
done
touch "$BATCH/VALIDATION_PASSED"
echo BATCH_EXECUTION_VALIDATED_NOT_PHYSICAL_PRODUCTION
