#!/usr/bin/env bash
set -euo pipefail
ROOT=$HOME/formulationlab-dds-livraghi-pilot-20261021
export PATH="$ROOT/lammps-dds/bin:$PATH"
date -u
hostname
git -C $HOME/formulationlab-upet-20261020/lammps-metatomic rev-parse HEAD || true
command -v lmp
lmp -help > reports/lammps-help.txt
grep -E 'REACTION|MOLECULE|KSPACE|bond/react|lj/cut/coul/long' reports/lammps-help.txt
grep -q 'bond/react' reports/lammps-help.txt
python3 --version
sha256sum inputs/data.initial inputs/*Top.txt inputs/map*.txt
