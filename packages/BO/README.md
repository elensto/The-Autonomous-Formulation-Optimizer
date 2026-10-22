# BO DGEBA / 4,4'-DDS

## Scope and truthfulness

This package adapts the GP target-window approach in Nannan's supplied **DGEBA/PGE/EDA** notebook to the team's **DGEBA/4,4'-DDS** pathway. It preserves Nannan's original notebook/script unchanged in `reference/`, imports her Matérn-5/2 `GPSurrogate` and probability-of-target-interval method, and replaces the old candidate chemistry. **No newly calculated DDS Tg values were provided. No trained GP model based on real DDS Tg results is claimed.** The only real simulation evidence in the source report is an ~8% conversion pilot and a reported 90%-target curing campaign, with no approved Tg observation at the time of the update.

## Two optimizer modes

**Mode A: `conversion`** selects *real, accessible* LAMMPS cured-network snapshots of a **fixed DGEBA:DDS composition**, varying achieved cure conversion. It only chooses a structure if its manifest row states `CURED_NETWORK_VALIDATED` and `physcheck_decision=ACCEPT`, and an actual file exists. The inherited GP requires ≥4 **distinct accepted conversion values**. Otherwise choose maximally diverse eligible cure states, without fake predictions. No available validated DDS cured structures were supplied, so initial mode A returns `NO_ELIGIBLE_CANDIDATES`. The current 8%-cure diagnostic is not an eligible cured-glass structure.

**Mode B: `composition`** searches **integer DGEBA/DDS molecule counts** plus requested cure conversion. It proposes **network-generation and curing jobs**, never pretends a new network is already available. Default DGEBA=100 with DDS count 40,45,50,55,60 (ratio 0.8,0.9,1.0,1.1,1.2) and target cure 60%,75%,85%,90%, with impossible 90% epoxy conversion at ratio 0.8 excluded. The candidate grid is a plan, not simulation results. Some off-stoichiometric compositions may be difficult to cure/gel and require PhysCheck. The default 150–220 °C target window is a research target, **not an achieved Tg**. No new structures or calculations are produced by BO itself.

Stoichiometry is the **amine reactive H / initial epoxy equivalent** ratio (not a DGEBA mole fraction):

\[r=\frac{4 N_{DDS}}{2 N_{DGEBA}}=\frac{2N_{DDS}}{N_{DGEBA}}.\]

For 100 DGEBA + 50 DDS: `r = 1.0`, initial 200 epoxides and 200 amine reactive hydrogens. For ratio 0.8, up to 80% **of initial epoxy sites** could react under ideal solely amine/epoxy stoichiometry (a theoretical ceiling, not a prediction of achieved conversion). Cure conversion must be measured from reaction outputs. In composition mode the optimizer should receive realized conversion, not simply the requested target, as the input for any accepted Tg label.

## Run it

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v

# Real chemistry: formulation proposals, no Tg labels yet -> diversity/exploration
python src/dds_optimizer.py --mode composition --generate-pool \
  --observations data/md_tg_observations.csv --output outputs/initial_composition.csv

# Real chemistry: cure-conversion selection requires actual validated cured files
python src/dds_optimizer.py --mode conversion \
  --candidates data/dds_network_manifest.csv \
  --observations data/md_tg_observations.csv --output outputs/initial_conversion.csv

# Synthetic software fixture only! This is NOT a real MD study.
python src/dds_optimizer.py --mode conversion \
  --candidates examples/DEMO_SYNTHETIC_ONLY/DEMO_candidates.csv \
  --observations examples/DEMO_SYNTHETIC_ONLY/DEMO_observations.csv \
  --allow-synthetic-demo --output outputs/DEMO_not_MD_proposals.csv
```

For a notebook experience run `jupyter notebook DDS_dual_mode_BO.ipynb`, or open the notebook in VS Code and run its cells. It calls the tested source module; this avoids divergences between a copied 500-line notebook and production code.

## Data contracts (CSV)

`data/dds_network_manifest.csv` / `data/dds_composition_candidate_pool.csv`:
`candidate_id,structure_id,n_DGEBA,n_DDS,ratio,conversion_pct,structure_status,file_path,physcheck_decision,job_status` plus optional fields. For conversion mode, `structure_status=CURED_NETWORK_VALIDATED` and `physcheck_decision=ACCEPT` are required, along with an existing file path. A *cured* artifact only becomes ready once the curing agent, LAMMPS, topology audit, charges, achieved conversion and equilibration have passed review.

`data/md_tg_observations.csv`:
`observation_id,candidate_id,n_DGEBA,n_DDS,ratio,conversion_pct,tg_C,tg_se_C,physcheck_decision,bo_trainable,source_type,experiment_status`.

Only `source_type=real_classical_MD`, `physcheck_decision=ACCEPT`, `bo_trainable=true`, and `experiment_status=completed` can train the real model. In practice the team should also record network replica, file SHA-256, force-field version, temperature scan, cooling rate, Tg fit quality, and provenance. Those are external PhysCheck obligations at present, not comprehensively enforced by this BO code.

The *same candidate* cannot be re-proposed when completed, pending or running. In `conversion` mode, the default policy also excludes already measured conversion levels across replicas, prioritizing distinct experimental design locations. If replicates are desired, implement a deliberate replicate policy.

## How the BO chooses

- **0–3 distinct approved points:** deterministic farthest-first initial exploration in normalized feature space, no numerical Tg predictions.
- **≥4 distinct approved points:** reuse Nannan's Matérn-5/2 Gaussian Process with uncertainty-aware training and the probability of satisfying the target Tg interval. Batch diversity is implemented with minimum normalized feature separation. With very sparse training data, GP predictions are exploratory and should be reviewed by a BO expert.
- **No eligible real conversion structures:** stop with `NO_ELIGIBLE_CANDIDATES` rather than proposing nonexistent LAMMPS files.

## Next integration tasks

1. **Curing teammate:** after Carbono job 175498, export one or more actual cured network files, exact verified conversion, charge conservation and connectivity. Verify the source 100:50 stoichiometry and distinguish requested from achieved conversion. New DDS counts need adapted packing and reaction topologies; do not change a count in a metadata file and call that a cured network.
2. **Experiment Agent:** consume proposed rows. If `proposed_job_type=existing_network_MD`, route to approved LAMMPS Tg job. If `generate_and_cure_network`, build/equilibrate network first, then evaluate Tg. Store job IDs and prevent duplicate in-flight submissions.
3. **PhysCheck:** confirm valid AMBER-style Livraghi chemistry/charges (not Fan CHARMM), equilibrium, Tg fit quality and comparable cooling rates. Only then mark training observations approved.
4. **BO reviewer:** validate small-data GP priors, interval acquisition performance, and treatment of independent replica uncertainties. Consider a simpler 1D monotonic/uncertainty-aware approach for cure conversion when only a handful of approved points exist.
5. **Integrator:** connect observation files to real experiment IDs, add persistent job reservations, status transitions and result storage. These BO scripts are *not* the complete autonomous agent pipeline.

## Known limitations

The composition grid changes stoichiometric ratio, but it does **not** vary resin/hardener chemical identities or introduce a diluent. Varying DDS amount changes both crosslink limits and network chemistry; not all requested cure levels may be reachable. Composition and conversion should remain distinct predictors. The GP does not correct different cooling-rate conditions automatically. The example dataset is emphatically **synthetic**, with placeholder `.data` files that are not LAMMPS structures. Never submit those examples to MD. The original Nannan notebook is archived in `reference/` unchanged.
