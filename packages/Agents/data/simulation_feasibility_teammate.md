# FormulationLab — Simulation Feasibility Update

**Recommendation:** use classical LAMMPS as the main Tg simulation engine for the 48-hour hackathon, retaining UPET for short exploratory tests and density as the fallback.

## System and target protocol

We tested a **7,452-atom DGEBA/Jeffamine D230 epoxy network at approximately 80% cure conversion**.

The proposed protocol consists of **2 ns of high-temperature NPT equilibration plus 15 ns of cooling from 573 to 273 K at 20 K/ns**. At a 1 fs timestep, this requires **17 million MD steps**.

## What we tried

- **Model size:** PET-MAD-M exceeded the A40 GPU's available memory for the full structure. PET-MAD-S fitted and ran. Larger models, including PET-OAM-XL, were only smoke-tested on a tiny system; full-system feasibility and Tg accuracy were not established.
- **LAMMPS + UPET:** the initial estimate of approximately **290 days already used GPU force inference**, with CPU integration. A later warmed benchmark projected approximately 347 days.
- **Compiled UPET + NVIDIA ALCHEMI:** successfully ran conservative force evaluation and NPT integration on GPU. Observed throughput improved by approximately **3×** relative to the warmed LAMMPS + UPET test, but remained insufficient for a 17 ns trajectory within the hackathon.
- **Classical LAMMPS:** compiled a separate OpenMP installation and successfully ran a short NPT diagnostic using the full structure, bonded topology and supplied coefficients on **16 CPU cores**, without UPET or GPU acceleration.

## Measured throughput and runtime projections

- **LAMMPS + PET-MAD-S, NVIDIA A40:** approximately **1.5–1.8 s/step**, projecting **290–347 days for 17 ns**.
- **Compiled PET-MAD-S + ALCHEMI, NVIDIA A40:** approximately **0.58 s/step**, projecting **114 days for 17 ns**.
- **Classical LAMMPS, 16 CPU cores:** approximately **9.8 ns/day**, projecting **42 hours for 17 ns**, or **7.4 hours for 3 ns**.

**These are extrapolations from short benchmarks, not completed Tg calculations.** The classical timed segment was only 5 ps, and the cell was still compacting. Throughput must be reassessed after relaxation. GPU acceleration of the classical force field remains an additional option; it has not yet been benchmarked here.

## Proposed 48-hour workflow

1. **Validate the classical force-field setup**, including mixing rules, cutoffs and 1–4 interactions. The `charmm` style name alone does not establish the parameterization; the benchmark retained the user's input conventions without validating their physical provenance.
2. **Start from the existing cured networks**, avoiding implementation of a new curing workflow during the event.
3. **Relax at high temperature under NPT** until volume and energy show acceptable stability.
4. **Run a preliminary stepwise cooling scan:** an initial budget of 0.5–1 ns for relaxation, followed by 16 temperature plateaus from 573 to 273 K, separated by 20 K, with 100–200 ps per plateau. This totals approximately **2.1–4.2 ns per replica**, projecting **5–10 hours** at the measured classical throughput. Extend relaxation or the temperature range when required.
5. **Extract Tg from specific volume versus temperature**, discarding plateau transients and fitting the low- and high-temperature branches. Report Tg as unresolved if both branches are not adequately sampled.
6. **Use three independent network configurations** to assess variability. Running three replicas simultaneously requires **48 CPU cores**, subject to queue availability.

Reserve the final **12 hours** for integration, analysis, documentation and the demo. Keep **density at a single temperature** as the fallback, with its own relaxation and convergence checks.

## Scientific and project scope

- The faster screening scan corresponds to approximately **100–200 K/ns** between plateaus, rather than 20 K/ns. Its result must be reported as a **preliminary, cooling-rate-dependent Tg**. Shorter simulations do not guarantee convergence or experimental accuracy.
- The **27 available structures represent one chemistry at different cure conversions**, with three configurations per conversion. They support an initial cure-conversion study; new formulation compositions require additional network generation and consistent parameterization.
- Under this proposal, **the reported Tg would come from a classical force field**. UPET would remain an exploratory component, rather than the engine producing the full Tg curve.
- **No production Tg scan has been completed or left running.**

## References

- [DGEBA/Jeffamine D230 structure dataset](https://data.mendeley.com/datasets/wmsbdm9bh7/1)
- [Fan et al., Polymer 2020 — source study](https://doi.org/10.1016/j.polymer.2020.122477)
- [Orselly et al., ACS Omega 2022 — Tg methodology, Section 3.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9434774/)
- [LAMMPS force-field styles and available acceleration options](https://docs.lammps.org/pair_charmm.html)
