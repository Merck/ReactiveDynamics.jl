# Pharma R&D twin: guide for agents working in this folder

This folder is a runnable simulation of a synthetic pharma R&D portfolio. It turns a NocoBase snapshot of six programs at tick T-042 (January 2026) into a stochastic, resource-constrained model, built on ReactiveDynamics.jl 0.3 and AlgebraicAgents.jl 0.4. The model compares four portfolio policies over 24 months. Use it to produce new views: rerun the baseline, run what-ifs on an edited snapshot, and read the result files. All data are synthetic: a fictional portfolio with no real company, compound, patient or market data. Drug targets are shown as generic labels (Target A to F).

## The baseline answer (what a correct setup reproduces)

The closed loop (an automatic tox kill switch plus a monthly wet-lab re-plan) is worth **+$53.3M ± 1.7M** net over 24 months against the static annual plan. That is +3.9% on the static plan's net value; the prototype deck estimated +14.8%. The kill switch alone is worth +38.3 ± 1.7, the re-plan alone +14.1 ± 1.5, and their interaction +0.9. Split by source, +35.6 is spend avoided (independent of erosion) and +17.7 is value protected (scales with erosion). Net value is risk-adjusted eNPV minus spend, in $M. The figures are paired comparisons over 40 seeds per arm, root seed 2026, and ± is one standard error.

## Files

| Path | What it is |
|---|---|
| `talk.html` | The talk page (open in a browser). Its numbers are the baseline run in `out/`. |
| `talk.pdf` | The same page as a PDF, with every collapsed section open. |
| `data/nocobase_snapshot.json` | The input, in NocoBase table shape. Every value the deck does not give is listed under `_assumed`. |
| `twin.jl` | The bridge: snapshot → model, the policy arms (`ARMS`), the clinical-operations agent, the metrics. |
| `run_twin.jl` | The driver. Prints §1 to §6 and writes everything into the output folder. |
| `vendor/chartkit.jl` | Dependency-free SVG charts used by §6. |
| `LICENSE` | The MIT license of ReactiveDynamics.jl, which covers this folder. |
| `out/` | The baseline run's outputs. Treat as read-only reference. |

Output files written by every run:

| File | Contents |
|---|---|
| `ensemble_summary.json` | The headline numbers. `arms.<arm>` holds `riskadj_enpv`, `spend`, `net` (mean, se) and `delta_net_vs_static` (paired delta, se); `sensitivity` holds the closed-loop delta at erosion 0%, 0.75%, 1.5%, 3%, 5%. |
| `ensemble_per_seed.json` | One row per arm × seed (160 rows): `riskadj_enpv`, `spend`, `net`, `realized_enpv`, and per program `killed_at`, `slip_months`, `readout`, `readout_t`. Rows with the same `seed_index` are the same future under different policies, so pair on it. |
| `simulation_ticks.json` | Month-by-month portfolio eNPV, spend and pool utilization for one future (seed 1) per arm, as NocoBase `simulation_ticks` rows. |
| `simulation_events.json` | Kills, readouts, GMP batches, accrual lags and slips for that future, as `simulation_events` rows. |
| `resource_allocations_plan.json` | The wet-lab FTE plan at every re-plan, per program and arm. |
| `scenarios.json` | One row per arm (SCN-01 static … SCN-04 closed loop) with its policy parameters. |
| `model_<arm>.rdj.json` | The compiled model per arm, plain JSON with no executable code. |
| `chart_*.svg`, `net_*.svg`, `aa_wiring.svg` | Dashboard charts and model drawings. |

## Run it

Needs Julia 1.12 (`juliaup add 1.12`). The ensemble seeds come from Julia's `hash`, which can change between Julia versions, so the baseline numbers replay exactly only on 1.12. Graphviz (`dot` on the PATH) is optional; without it the four model drawings are skipped and everything else is produced.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'    # first time only
julia --project=. run_twin.jl --nseed=5 --out=scratch/smoke  # quick check, ~1 min
julia --project=. run_twin.jl                            # full baseline, 40 seeds per arm, ~2 min after the first compile
```

Flags: `--nseed=N` (seeds per arm, default 40), `--snapshot=path.json` (input, default `data/nocobase_snapshot.json`), `--out=dir` (output folder, default `out/`).

## Run a what-if

1. Copy the snapshot: `whatifs/<name>.json`. Edit only the fields the question needs.
2. Run: `julia --project=. run_twin.jl --snapshot=whatifs/<name>.json --out=runs/<name>`.
3. Compare `runs/<name>/ensemble_summary.json` with `out/ensemble_summary.json`, and say which fields changed.

Never write into `out/`; it is the reference that `talk.html` quotes. Use `--nseed=40` (the default) for any number you report. Smaller runs are for checking that an edit works.

## The inputs you can change

| Snapshot field | What it drives | Baseline |
|---|---|---|
| `valuation.erosion_per_month` | eNPV lost per month of slip. Applied when results are valued, so the simulation itself does not change. The biggest single assumption. | 0.015 |
| `kill_switch.herg_ic50_um_threshold` | The kill rule stops a study whose lead's hERG IC50 is below this (µM). | 1.0 |
| `compounds_candidates[].herg_ic50_um` | Each lead's hERG value. Only PRG-402 (0.8) is below the threshold. | per row |
| `compounds_candidates[].kits_on_hand`, `reorder_point_kits`, `batch_kits` | Kit stock, the level that triggers a GMP batch, and kits per batch. | per row |
| `gmp_batch.lead_time_months`, `gmp_batch.slots_per_month` | GMP batch duration, and Suite B slots a running batch holds. | 3, 3 |
| `therapeutics_programs[].stage_pos` | Probability of passing the readout that ends the current stage (assumed; the deck gives only cumulative PoS). | per row |
| `therapeutics_programs[].pos_estimate`, `enpv_usd_m` | Cumulative PoS and unadjusted eNPV ($M). | per row |
| `studies[].months_remaining_full_speed` | Months to readout at full speed. | per row |
| `clinical_trials[].site_count`, `accrual_per_site_month`, `target_enrollment`, `current_enrollment` | Enrollment: Poisson(site_count × accrual) patients a month against the need. | per row |
| `resources_capacity[].total_capacity_units`, `cost_usd_m_per_unit_month` | Pool sizes and unit costs (wet-lab FTE plus five shared pools). | per row |
| `resource_allocations[]` | Per program: wet-lab FTE held and needed at full speed, shared-pool draws, CRO $M a month. | per row |

Policy timing is code, not data: `ARMS` and `ANNUAL_REVIEW` in `twin.jl` set when the static plan acts on a tox flag (month 12) and how often each arm re-plans. To price a different wait, such as acting at month 3, copy `twin.jl` into the what-if folder, change `kill_at`, and say plainly that the policy, not the data, changed.

## Rules for reporting a number

- Every number you report must come from an output file of a run you made or from `out/`. If you edit an input, rerun before quoting anything.
- Report policy effects as paired deltas with their standard error, at 40 seeds. A delta smaller than about two standard errors is not a finding.
- One simulated future (the tick and event files, the talk page's matrix) shows how something happens, not how much it is worth. End-of-run eNPV varies by about $595M between futures, 11× the headline effect.
- Same snapshot, same seeds, same numbers. If a rerun of the unedited baseline does not reproduce +53.3 ± 1.7, first check `julia --version` is 1.12, then stop and report the mismatch with the package versions (`julia --project=. -e 'using Pkg; Pkg.status()'`).
- Do not change ReactiveDynamics or AlgebraicAgents code. Everything here runs on their released versions.
- Say which assumption a result depends on. Erosion, per-readout PoS and kit batch size and lead time are assumed, not taken from the deck.
- The data are synthetic. Do not present a result as a forecast for a real portfolio.

## Views worth producing

- **Erosion break-even.** At what erosion rate would the closed loop reach the deck's +14.8%? Edit `valuation.erosion_per_month` and rerun, or read it off `sensitivity` in `ensemble_summary.json`, and say which you did.
- **The cost of waiting on a tox flag.** Price acting at months 1, 3, 6 and 12 by changing `kill_at` in a what-if copy of `twin.jl`. This is the portfolio team's first question on the talk page.
- **PRG-518's kit supply.** Vary `batch_kits` and `gmp_batch.lead_time_months` for PRG-518 and show when manufacturing stops limiting it, using `slip_months` in `ensemble_per_seed.json`.
- **Per-readout PoS.** Replace the assumed `stage_pos` values with other benchmarks and show how much the headline moves.
- **Per-program view.** From `ensemble_per_seed.json`, the distribution of slip and readout timing per program under each arm.
- **A different portfolio.** Replace the snapshot's rows with another set of programs in the same shape, rerun, and compare.
