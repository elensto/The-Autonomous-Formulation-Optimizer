#!/usr/bin/env bash
# All inputs, data, scripts and results live in this isolated batch snapshot.
set -euo pipefail
BATCH="$1"
CASE="$2"
MODE="${3:-production}"
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
LMP="$ROOT/lammps-dds/bin/lmp"
if [[ "$MODE" == smoke ]]; then WORK="$BATCH/validation/$CASE"; else WORK="$BATCH/runs/$CASE"; fi
mkdir -p "$WORK"
cd "$WORK"
if [[ -f STARTED ]]; then echo "Refusing to overwrite an existing case: $WORK"; exit 2; fi
touch STARTED
cp "$BATCH/common/"* .
cp "$BATCH/cases/$CASE/data.initial" "$BATCH/cases/$CASE/case.json" .
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_PROC_BIND=spread OMP_PLACES=cores
read -r EPOXIES TARGET SEED < <(python3 -c 'import json;d=json.load(open("case.json"));print(d["epoxy_groups"],d["target_epoxy_conversion_percent"],d["seed"])')
date -u
hostname
echo PREPARATION | tee phase.txt
PREP=();CURE=();EQ=();TG=()
if [[ "$MODE" == smoke ]]; then
 PREP=(-var heat_steps 100 -var equil_steps 100 -var cure_steps 100)
 CURE=(-var chunk_steps 100 -var blocks 1 -var cycles 2 -var ramp_steps 100 -var anneal_steps 100 -var postanneal_steps 100)
 EQ=(-var equil_steps 100)
 TG=(-var heat_steps 100 -var hold_steps 100 -var plateau_equil_steps 100 -var plateau_prod_steps 100)
fi
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.prepare -var epoxies "$EPOXIES" -var seed "$SEED" "${PREP[@]}" -log log.prepare.lammps
python3 export_topology.py pilot.data data.start
echo CURING | tee phase.txt
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.cure_campaign -var epoxies "$EPOXIES" -var target "$TARGET" "${CURE[@]}" -log log.cure.lammps
python3 parallel_analysis.py .
if [[ "$MODE" != smoke ]] && ! python3 -c 'import json,sys;sys.exit(0 if json.load(open("campaign_summary.json"))["structure_gate_passed"] else 1)'; then
 echo CURE_GATE_NOT_PASSED | tee phase.txt
 exit 0
fi
echo EQUILIBRATING | tee phase.txt
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.equilibrate "${EQ[@]}" -log log.equilibrium.lammps
if [[ "$MODE" != smoke ]] && ! python3 analyze_properties.py equilibrium .; then
 echo EQUILIBRIUM_GATE_NOT_PASSED | tee phase.txt
 exit 0
fi
echo TG_SCAN | tee phase.txt
"$LMP" -nonbuf -sf omp -pk omp "$OMP_NUM_THREADS" -in in.tg_screen "${TG[@]}" -log log.tg.lammps
if [[ "$MODE" == smoke ]]; then
 echo SMOKE_ONLY_NOT_PHYSICAL_RESULT | tee phase.txt
else
 python3 analyze_properties.py tg .
 echo PRELIMINARY_PROPERTIES_COMPLETED | tee phase.txt
fi
date -u
