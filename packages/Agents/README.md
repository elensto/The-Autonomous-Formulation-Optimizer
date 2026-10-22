# FormulationLab — Dual-path Experiment + PhysCheck Agents

**2-day hackathon starter package, 21 October 2026.** Python standard library only. It is a real, executable *control-plane prototype*, not a ready-to-run cured-network molecular dynamics simulation or a network curing algorithm.

## What it supports

| Route | BO proposal | Experiment Agent | PhysCheck | Status |
|---|---|---|---|---|
| A: DGEBA/Jeffamine D230 | `mode=cure_conversion`, exact conversion or `structure_id` | Search verified Fan network, prepare screened LAMMPS cooling job | Checks source topology, metadata, force-field approval, relaxed restart and protocol | **Works once exact vetted archive files/restarts are imported** |
| B1: DGEBA/DDS cure-conversion | `mode=cure_conversion`, exact conversion or `structure_id` | Select DDS cured networks if they exist | Same checks but DDS-specific force-field ID; uncured Packmol inputs rejected | **No validated cured DDS networks currently supplied** |
| B2: DGEBA/DDS formulation | `mode=formulation`, `formulation_id`, `components`, conversion/structure | Search for a matching, already generated cured network | Confirm resin/hardener identity, validated cure topology/provenance and FF | **Interface ready, needs external curing backend** |

**No automatic chemistry conversion:** Even for DDS, BO proposals never create a new network without a validated curing pipeline. The agent responds with `MODIFY: NO_CURED_NETWORK_FOR_FORMULATION` and recommends external curing/structure import. It does not pretend a Packmol mixture is a cured network.

## Scientific status of supplied data

`data/network_manifest.csv` reconstructs the previous manifest schema from conversation records: **27 D230 archival placeholders, five DDS uncured structure placeholders, and one *unlinked* teammate benchmark evidence record** (which is *not* an additional unique cured network). Names, conversions, replicas and hashes for the 27 original archive files were not available in this session and are intentionally blank. **None of these records can start MD yet.**

The teammate reported classical MD smoke testing on a 7,452-atom approximately 80%-cured DGEBA/D230 network, **not** a completed/equilibrated Tg scan. The DDS paper provides starting configurations and scripts, not five ready-cured Tg networks. The five records in `data/literature_Tg_NOT_BO.csv` include Orselly (2022) Figure 8; the literature records must never be fed to the real BO as new Tg observations.

The original uncured composition library is retained as `data/original_uncured_components_REFERENCE_ONLY.json`, not used as cured Tg network data.

## Repository organization

```
FormulationLab_DualPath_Agents/
  agents/
    core.py                 Unified IDs, CSV readers, structure lookup, PhysCheck rules
    experiment_agent.py     Proposal -> verified structure -> job.json + guarded in.tg_scan.lmp
    physcheck_agent.py      Standalone pre/post decisions
  config/
    protocol_D230.json      D230 thermal screening protocol (not yet scientist-reviewed)
    protocol_DDS.json       DDS thermal screening protocol (not yet scientist-reviewed)
    lammps_restart_header_*.REVIEW_REQUIRED.lmp  Source FF template placeholders
  data/
    network_manifest.csv    Status-aware manifest, not invented files/conversions
    literature_Tg_NOT_BO.csv
    approved_Tg_observations_EMPTY.csv
    simulation_feasibility_teammate.md
    original_uncured_components_REFERENCE_ONLY.json
  examples/
    bo_proposal_D230.json
    bo_proposal_DDS.json
    bo_proposal_DDS_composition.json
    md_result_PENDING.json
  raw/FAN_D230/, raw/LIVRAGHI_DDS/   **empty source archive locations**
  tests/test_agents.py
  outputs/                   Generated decisions and future jobs
```

## Try it (no extra packages)

From repository root:

```bash
python -m unittest discover -s tests -v
python agents/experiment_agent.py \
  --proposal examples/bo_proposal_D230.json \
  --manifest data/network_manifest.csv \
  --protocol config/protocol_D230.json \
  --out outputs/d230_decision.json
python agents/experiment_agent.py \
  --proposal examples/bo_proposal_DDS_composition.json \
  --manifest data/network_manifest.csv \
  --protocol config/protocol_DDS.json \
  --out outputs/dds_decision.json
python agents/physcheck_agent.py post --result examples/md_result_PENDING.json
```

Expected: `MODIFY` for D230 (`NO_MATCHING_VERIFIED_CONVERSION`) and DDS composition (`NO_CURED_NETWORK_FOR_FORMULATION`). This is correct behavior on unverified inputs. `post` also outputs `MODIFY` because no real Tg calculation has occurred. Six tests include a **synthetic file/restart fixture**, proving the pipeline can generate a job in a fully populated test registry. Test fixtures are **not physically validated simulations**.

## Team handoff — bring a real network online

1. Add the original Fan data files into `raw/FAN_D230/`, or your cured DDS files to `raw/LIVRAGHI_DDS/`. Do not alter originals. A built/converted DDS network is a **new** `generated_cured_network` row, not one of the original uncured Packmol starts.
2. Resolve exact filename, SHA-256, conversion %, replica ID, formulation ID, atom count, source DOI and force-field family. Match the teammate's D230 network evidence to its actual file. A cured row is `record_kind=cured_network` or `generated_cured_network`; never set uncured inputs to cured merely to pass a check.
3. The MD scientist must verify topology and original pair/bond/angle/dihedral/improper/charge, 1–4 nonbonded and long-range conventions. Assign `force_field_id=FAN_CHARMM_REVIEWED` only for verified Fan-style networks or `DDS_AMBER_REVIEWED` for verified DDS-style networks. **These names are project gate IDs, not assertions of physical compatibility.** Set `topology_reviewed=true`, `force_field_reviewed=true` only after manual review.
4. Perform high-temperature NPT compression and equilibration **as a separate MD preparation step**, then save a LAMMPS `write_restart` file. Record it in `equilibrated_restart_path`, and set `npt_equilibrated=true` only after volume/energy convergence checks.
5. Review the temperature/cooling protocol and replacement LAMMPS header with simulation owner. Source-provided FF restart should hold the potential settings. Create a reviewed script without `REVIEW_REQUIRED` in its path/content, then update `template_file` and `reviewed_by_scientist=true` in the protocol. Check matching `units real`, pressure atmosphere, long-range Coulomb support and restart-version compatibility. **Do not mix D230 CHARMM with DDS AMBER**.
6. Run `Experiment Agent` with a BO proposal referencing that exact `structure_id`, or a unique verified conversion/replica. The job is written to `outputs/jobs/<request_id>/`, with `job.json`, `in.tg_scan.lmp` and `run_local.sh`. Job creation **does not launch MD**. Use `bash run_local.sh` only on a properly configured simulation machine after review.
7. Once the MD simulation and Tg fit complete, run `PhysCheck post` with a real result JSON and preserve the full trajectory, protocol, fitting uncertainty and trace path. Only `ACCEPT` results may become BO observations; `literature`, `synthetic`, and `UPET_singlepoint` results are never eligible.

### Example manifest row for a genuinely validated cured structure (illustrative format only)

```csv
structure_id,system_id,formulation_id,record_kind,relative_path,file_sha256,verified_cure_percent,replica_id,number_atoms,force_field_id,topology_reviewed,force_field_reviewed,npt_equilibrated,equilibrated_restart_path,cure_protocol_id
FAN_ACTUAL_ID,SYS_D230_FAN,DGEBA_D230_NOMINAL,cured_network,raw/FAN_D230/ACTUAL_FILE.data,ACTUAL_SHA256,80,1,ACTUAL_ATOMS,FAN_CHARMM_REVIEWED,true,true,true,raw/FAN_D230/ACTUAL_RELAXED.restart,
```

This illustrates required fields; **none of the placeholder values should be copied as a verified source record**.

## Agent response codes

| Code | Meaning | Action |
|---|---|---|
| `ACCEPT` | Source file exists, hash and header consistent, necessary review flags confirmed, applicable protocol approved and valid relaxed restart available | Prepare LAMMPS job (does not execute) |
| `MODIFY` | Incomplete file, metadata, equilibration, protocol, ambiguous replicate or absent cured DDS topology | Ask teammate for specific missing inputs, or invoke an external curing service |
| `REJECT` | Unsafe/invalid request, unsupported chemistry, uncured file masquerading as cured, path escape, hash mismatch or fundamentally broken header | Do not dispatch |

The PhysCheck Agent is **rule-based**, not an LLM claiming to have inspected bonds/charges. Human review flags are audit assertions, not validated automatically. No assumption is made that short Tg MD is converged.

## Interfaces for the integration developer and Nannan

`bo_proposal`:

```json
{"request_id":"EXP_0007","mode":"cure_conversion","system_id":"SYS_D230_FAN","structure_id":"VERIFIED_ARCHIVE_ID","target_tg_C":85}
```

For composition mode, include `mode=formulation`, stable `formulation_id`, `components: {"resin_id":"R001", "hardener_id":"H007"}`, and a cure conversion target. Only DGEBA/DDS composition is supported by this particular DDS backend; adding a new hardener requires a new system record and validation of its cure method. `replica_id` prevents ambiguous requests. The BO must not propose continuum conversions with no matching structure.

`experiment_decision`: `status`, `reason` or `checks`, `next_action`; when accepted `job_dir`.

`job.json`: `request_id`, `mode`, `system_id`, `structure_id`, `formulation_id`, `verified_conversion_pct`, `protocol_id`, `force_field_id`, `simulator`, `status`, `expected_data_pattern`, `output_dir`, `checks`.

`md_result`: `structure_id`, `source_type=real_classical_MD`, `property=Tg`, `tg_K`, `tg_uncertainty_K`, `fit_accepted`, `npt_converged`, `protocol_reviewed`, `cooling_trace_path`. Always retain raw thermal traces and fit report alongside the result. A real BO training pipeline should additionally confirm that `structure_id` is present in a verified manifest and that the simulation job uniquely matches it.

## What remains outside these two agents

- Automated DDS crosslinking and formulation-to-network chemistry (a separately validated Curing Agent).
- NPT network compression and equilibrium, classical LAMMPS execution, Tg fit/uncertainty scripts.
- BO Gaussian Process fitting, acquisition function and database storage; Nannan's notebook has a separate composition space.
- UPET PET-MAD-S calculations (separate physics diagnostic module, not a Tg label).
- Production operation: queue scheduling, timeouts, checkpointing, software versions, audit permissions and force-field benchmarking.

**Pragmatic scope:** bring the D230 branch online first, once its exact structure and relaxed restart are shared. Then use the same agents for DDS once genuinely cured networks and a reviewed DDS-specific LAMMPS template become available.

## NEW — DDS Curing Agent (2026-10-21 reported pilot)

`agents/curing_agent.py` is an **external-workflow adapter**, not a reaction chemistry implementation. The MD teammate reports a 6,350-atom DGEBA/DDS 100:50 system and an 8% curing pilot. HPC job `175498` was **reported running toward 90%** and `175499` was waiting; **90% achieved conversion and Tg have not been established**. See `data/dds_campaign_status_2026-10-21.json` (not a live scheduler query).

The agent supports both routes correctly:

* **DGEBA/D230:** the Experiment Agent resolves an independently verified *existing* cured network, and the Curing Agent is not called.
* **DGEBA/DDS conversion:** the Experiment Agent first seeks an already validated cured network; when absent, it can route a request to the Curing Agent. It prepares a **manifest-only job request** referencing the team's scientifically reviewed AMBER-style LAMMPS `bond/react` workflow. Actual curing inputs must be supplied by the MD teammate.
* **DGEBA/DDS composition:** supported in proposal schema, but altered DGEBA/DDS counts are deliberately blocked until a new starting box, reaction rules and force-field compatibility have been reviewed. A fixed 100:50 baseline alone does **not** support composition optimization.

### Example commands (run from the package root)

```bash
python agents/curing_agent.py prepare --proposal examples/bo_proposal_DDS_cure_request.json
python agents/curing_agent.py assess --event data/dds_campaign_status_2026-10-21.json
python -m unittest discover -s tests -v
```

The default curing config (`config/curing_DDS_EXTERNAL.json`) **must return MODIFY** until the MD team supplies `TEAM_REVIEWED_cure.lmp` and `TEAM_REVIEWED_start.data`, corrects their paths in the config, sets sha256 values, and confirms `scientist_reviewed: true`. It NEVER submits an HPC job. It never claims to validate LAMMPS reaction syntax or chemistry. Use the team's real Carbono SLURM procedure and job dependencies, not a generic fabricated submission script.

After curing, the external status checker requires an explicitly achieved conversion, reviewed bond graph, conserved charge, reviewed topology/force field and a structure hash. An `ACCEPT` from that check allows *consideration for importing the topology*; it **does not** mean the material is NPT equilibrated or Tg-BO eligible. The existing Experiment and PhysCheck agents still require an equilibrated restart, corresponding manifest record, trusted force field and Tg scan protocol.

To complete integration: implement a workflow dispatcher that calls `make_job` first and conditionally calls `prepare_curing_job` when its `next_action` equals `CALL_CURING_AGENT_OR_IMPORT_VALIDATED_NETWORK`. Separately implement import of science-reviewed curing outputs into the manifest; it is intentionally not automatic.

### One-command dual-path dispatcher

```bash
python agents/workflow.py --proposal examples/bo_proposal_D230.json
python agents/workflow.py --proposal examples/bo_proposal_DDS_cure_request.json
```

The dispatcher returns `PREBUILT_NETWORK_LOOKUP`, `DDS_CURING_REQUIRED`, or `EXISTING_CURED_NETWORK_TO_TG`. The DDS branch remains `MODIFY` without the actual reviewed scripts and molecular input data. Real files and evidence are never generated by these commands.
