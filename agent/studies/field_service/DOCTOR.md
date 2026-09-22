# DOCTOR — plausibility gate for the field-service study

> The DIAGNOSE artifact (AGENT_WORKFLOW_CHARTER §3.3). Every line below is machine-produced evidence from a command recorded here, not an opinion. `validate` proved these four documents are well-formed; this phase asks whether they are about the operator's system.

## The command

Run once per arm, with the flags **pre-registered in `BRIEF.md` before any model existed**:

```bash
julia --project=. agent/scripts/rdx.jl diagnose agent/studies/field_service/<arm>.rdj.json --seed=1 \
  --expect-binds=laptop --expect-reusable=tech,laptop \
  --tally=breach_pharma,breach_acad,backlog_meter,fixed_pharma,fixed_acad
```

No `--host` flag: this study uses plain places only, no structured tokens, so there is no `host.jl` and no authored Julia in the loop at all.

## Verdict

| Arm | pass | warn | fail | exit |
|---|---|---|---|---|
| `model.rdj.json` (baseline: 7 techs, 4 laptops, strict pharma-first) | **8** | 1 | **0** | 0 |
| `arm_hire_2_techs.rdj.json` (9 techs) | **8** | 1 | **0** | 0 |
| `arm_buy_2_laptops.rdj.json` (6 laptops) | **8** | 1 | **0** | 0 |
| `arm_reprioritize.rdj.json` (pharma priority 10.0 → 1.5) | **8** | 1 | **0** | 0 |

**Exit 0 on all four arms — clear to report.** Every arm was diagnosed, not just the baseline, because a lever can move which resource binds.

## The checks, with evidence (baseline, seed 1)

| Check | Verdict | Evidence |
|---|---|---|
| `transitions_fire` | pass | all 11 transitions fired; lowest non-meter count is `pharma SLA breach accrual` at 51 |
| `completions_exist` | pass | 31,053 successful completions over the horizon |
| `contention_exists` | pass | `laptop` and `tech` both run to **0**; min in-flight fill fraction 1.0 |
| `intended_resource_binds` | pass | `laptop` runs to zero — **but see the caveat below, this check does not discriminate** |
| `modality_intent` | pass | "Every resource declared reusable is conserved on all its input arcs: **tech, laptop**" — both, on all six input arcs |
| `horizon_sufficient` | **warn** | 4 instances in flight at the horizon; `breach_*`, `backlog_meter`, `fixed_*` still trending |
| `bounded_pools` | pass | clean once the five cumulative counters were declared as tallies |
| `ledger_wired` | pass | reward −508,600, cost 0.0 — the penalty ledger moves |
| `stochasticity` | pass | "no explicit success probabilities" — see the caveat below |

## The caveat that matters most: `intended_resource_binds` does not discriminate

Re-running the baseline with the **counter-hypothesis** `--expect-binds=tech` also returns **pass**:

```
pass intended_resource_binds  {"message": "The intended resource binds: tech runs to zero.", ...}
```

The check asks only whether the named pool reaches zero, and under this load *both* pools do. So it confirms the brief named a resource that genuinely runs out, and it would have caught a model where laptops were never scarce — but **it is not evidence that laptops rather than technicians are the constraint.** That claim rests entirely on the arm comparison in `FINDINGS.md`, where the two capacity levers are priced against each other. Recording this explicitly because a check cited as proof of something it cannot establish is exactly the decoration this phase exists to prevent.

The discriminating signal visible *within* a single run is the queue asymmetry, and it is stark. Seed 1, baseline:

| Queue | init | max | final |
|---|---|---|---|
| `q_acad_rig` (academic, needs the diagnostic rig) | 10 | **40** | 1 |
| `q_acad_quick` (academic, no rig) | 2 | 6 | 0 |
| `q_pharma_rig` (pharma, needs the rig) | 2 | 5 | 0 |
| `q_pharma_quick` (pharma, no rig) | 0 | 5 | 0 |

The same seven technicians serve all four queues. Only the two rig-requiring queues build, and the unprotected one builds eightfold. Work that needs a laptop waits; work that does not, does not.

## Reading contention correctly here

`min_saturation` is 1.0, which looks like nothing is short. It is not a contradiction: a `conserved` resource rations by **refusing admission**, so a repair that cannot get a technician *and* a laptop never starts, and everything that did start was fully served. The evidence for contention is therefore the pool trough at 0 plus the building rig queue, not fractional saturation — exactly the case `rd-doctor` describes. Both pools trough at 0 on every arm.

## Warnings carried forward as scope, not cleared

**`horizon_sufficient` (warn, all four arms).** 3-5 repairs are mid-job at tick 2080 and the cumulative counters are still climbing. This is expected and is not a defect: the counters are cumulative by construction (that is what `--tally` declares), and 3-5 jobs in flight out of ~2,100 repairs a year is a rounding error, not truncated value. The consequence for the report: every figure is **"over 12 months," not a steady state.** The trend detail shows `q_acad_rig` at −0.29 per unit time in the final quarter — the academic rig queue is flat-to-shrinking at the calibrated operating point, so the baseline is not a system in runaway.

**`stochasticity` (pass, with a qualification worth stating).** The message is "no explicit success probabilities (all transitions complete)," because no transition uses `prob_of_success`. The model's randomness is real but comes from two other places: Poisson call arrivals, and the `choose`-valued arc multiplicities that route one job in four back for a revisit. So the ensemble error bars reflect **arrival timing and revisit incidence**, which are the right uncertainties for a queueing question — but they do not include any uncertainty in the operator's input estimates. Parameter uncertainty is handled separately, by the duration sensitivity recorded in `BRIEF.md` and `FINDINGS.md`.

## Checked by hand, beyond the nine

**Magnitudes are believable.** Baseline seed 1 completes 754 pharma + ~1,330 academic repairs over 12 months ≈ 2,090 jobs/year against 41 calls/week × 52 = 2,132 calls plus revisits. Throughput tracks demand, which is what a near-balanced service system should look like. Annual penalty on this seed is \$508,650 against the operator's stated status quo of \$427,800 — same order, one seed.

**The calibration anchor is met.** A 24-seed baseline ensemble (root seed 2026) gives mean waiting backlog **9.27 ± 0.70 jobs**, plus Little's-law in-flight work of 2,123 repairs/yr × 4 ticks ÷ 2,082 ticks = **4.08 jobs being worked**, for **13.35 total open tickets** against the operator's observed "about 14 or 15." The operator's count is of open *tickets*, which includes jobs a technician is currently on; the model's `backlog_meter` counts only *waiting* jobs, so the two are compared correctly only after adding work in progress. That interpretation is disclosed rather than absorbed.

**Two seeds agree qualitatively.** Seeds 1 and 2026-family members all show the same story: rig queues build, non-rig queues clear, both pools trough at zero. No seed reverses which queue accumulates.

**Single-seed arm comparison is NOT safe here, and this was nearly a trap.** The four arms at seed 1 give annual penalties of \$508,650 (baseline), \$685,600 (+2 techs), \$214,200 (+2 laptops), \$411,150 (reprioritize). Changing any parameter changes the RNG stream, so these are four different realizations of a near-critical system, not a controlled comparison — the apparent "hiring makes it \$177k worse" is not readable from one seed. All reported effects come from 48-member ensembles via `rdx effect` at a common root seed.

## Two engine/CLI defects found and worked around (not in `src/`)

1. **`validate` accepts a bare number inside a `choose` alt; `from_json_model` cannot construct it.** `rdx validate` returned exit 0 and zero diagnostics on a document whose `choose` alts held `"value": 1.0`; `rdx run` on the same file then died with `MethodError: no method matching node_from_dict(::Float64)` (`src/serialize.jl:84`). The bare-number shorthand documented for node positions is not honored inside `choose` alternatives. Worked around in this study by writing every alt value as an explicit `{"node":"const",...}`. This is the same class of validate/construct asymmetry as charter W14/W15 — a static pass that passes a document the constructor rejects. No `src/` change was made.

2. **Repeated CLI flags silently collapse to the last occurrence.** `--expect-reusable=tech --expect-reusable=laptop` checked only `laptop`, and five repeated `--tally` flags declared only the last, which left `bounded_pools` warning about four legitimately-cumulative counters. Comma lists (`--expect-reusable=tech,laptop`) work correctly. The first diagnose sweep reported 7 pass / 2 warn purely because of this; the corrected invocation gives 8 pass / 1 warn. `rdx.jl:53` documents the collapse for `--metric` but the `--expect-*`/`--tally` flags have the same behaviour without the same warning.

## Cycletime quantization — the grid was too coarse for the question, and this section is a correction

**Superseded figures, stated first so nobody quotes them.** An earlier version of this section reported that "the system is stable at 4 hours of resource occupancy and in runaway at 5," that the operator's "4 to 5 hours" therefore "admits exactly two distinct models," and that the 5-tick case gives **mean backlog 142 jobs and 1,023 academic breaches/month**. Those three claims are **retired.** They describe the authored `dt = 1.0` grid, not the operation. The surviving observation is narrower: job durations of 4.1, 4.2, 4.3 and 4.5 h produced **byte-identical** runs at `dt = 1.0` while 4.0 differed sharply — which is a property of the discretization and was mistaken for a property of the system.

**What is actually true.** `cycletime` is expressed in time units while progress accrues `dt` per tick (`firing.state += qs[i]*state.dt`, completion at `>= transCycleTime`, `src/solvers.jl:483,684`), so a firing occupies `ceil(cycletime/dt)` ticks. At `dt = 1.0` every duration in (4.0, 5.0] collapses onto 5.0 h. **The operator's own 4.3 h working figure therefore fell on the divergent side of the grid, while the calibrated baseline sat at 4.0 h** — the model was calibrated on a grid too coarse to represent the question it was built to answer, and the tie-break between the only two representable points was made by the observed backlog rather than by the duration.

**What refining it changed — and this paragraph has itself been corrected twice.** Re-authored at `dt = 0.25` (a one-line `meta.dt` change, since rates are per unit time and drawn as `rand(Poisson(dt*rate))`, `src/solvers.jl:316`, so neither rates nor `tspan` need rescaling), the same nominal 4.0 h moves **one tier and nothing else**. Matched 24-seed baseline ensembles at root seed 2026, both grids, both exit 0:

| Baseline, n = 24 both grids | `dt = 1.0` | `dt = 0.25` | difference | σ |
|---|---|---|---|---|
| **pharma breaches/yr** | 50.50 ± 1.61 | 37.62 ± 1.69 | **−12.88 ± 2.33** | **5.52 — real** |
| academic breaches/yr | 392.58 ± 32.80 | 396.00 ± 42.37 | +3.42 ± 53.58 | 0.06 |
| annual penalty | \$467,279 ± 24,502 | \$415,425 ± 30,965 | −\$51,854 ± 39,487 | 1.31 |
| backlog meter | 21,996 ± 1,699 | 21,603 ± 2,175 | −393 ± 2,760 | 0.14 |
| year-end acad rig queue | 12.50 ± 2.08 | 7.71 ± 1.93 | −4.79 ± 2.84 | 1.69 |

**Two claims from earlier drafts of this paragraph are retired.** It reported mean backlog **5.61** and **15.4** academic breaches/month at `dt = 0.25` against the `dt = 1.0` baseline's 10.13 and 31.2, and concluded that the coarse tick was charging up to a full tick of re-admission delay against every job so the baseline's *effective* occupancy exceeded its nominal 4.0 h. **Both are withdrawn.** The 5.61/15.4 pair was a **single seed**; the 24-seed means are 10.39 waiting and 33.00 academic breaches/month, against 10.57 and 32.71 on the coarse grid — the same operation to within noise. And the re-admission explanation was a capacity claim, which the integrated queue refutes at **0.14 σ** with an undetermined sign (−393 here, +507 when the two grids were compared at mismatched n). A later draft also cited \$4,200·Δpharma + \$650·Δacademic reproducing Δpenalty as evidence of a tier redistribution; that is an **identity** — this ledger defines `reward` as exactly that combination — and establishes nothing about mechanism.

**What is left is a real pharma improvement whose mechanism this design cannot resolve.** −12.88 breaches/year at 5.52 σ, with the academic tier flat to 0.06 σ. A redistribution needs a recipient and none is visible; nor could one be, because the academic tier's noise floor (±42.37/yr) is **3.3×** the size of the whole pharma gain, and resolving a transfer that small would need the difference's se down to 6.44/yr from 53.58 — about **1,700 seeds per arm**. The cheap discriminating test is the equal-priority arm run at both grids: if the pharma gap closes when the 10:1 weighting is flattened, the mechanism is the priority interaction.

**The finer-grid calibration probe is single-seed and is carried as an indication only.** At `dt = 0.125` with a 4.125 h hold (33 ticks exactly), seed 1 gives mean backlog 12.29, pharma breaches 4.00/month against a target of 4, and a penalty of \$553,250 against the `dt = 1.0` seed-1 \$508,650 — suggesting the observation-matching duration is near **4.10 h** of laptop checkout. One seed cannot carry that, and the earlier "mean backlog goes 5.61 → 55.1 between 4.00 and 4.25 h effective" range inherits the retired 5.61 endpoint; the defensible statement is that sensitivity to duration in this region is steep, with the slope unquantified.

**What refining it did not change.** The ranking of the three levers, and the sign and rough magnitude of every effect. The true cliff is at **4.366 h** by capacity arithmetic (36.64 laptop-needing jobs/week against 160 laptop-hours/week), the operator's 4.3 h sits at **98.5%** of laptop capacity with 0.57 jobs/week of headroom, and **two more technicians move that cliff by 0.000 h** while two more laptops move it to 5.463 h. On the effect side the finer grid reads **+\$288,644 ± \$31,502** for the two laptops against the published **+\$228,050 ± \$18,712** — a 1.65 σ difference, so it establishes no larger effect, but there is **no evidence the published headline is an overstatement**, and the lever clears by ≥9 σ on both grids. That one-sided statement replaces the earlier "±10% of tick-resolution uncertainty" caveat. Neither grid is demonstrably better calibrated: the operator's \$427,800 falls between them (+7.4% hot at n=48, +9.2% at n=24, −2.9% cold on the fine grid) and the matched-n penalties differ at 1.31 σ. The academic breach hazard, fitted at `dt = 1.0`, would need refitting before the finer-grid model is quoted for absolute breach levels.

## Run cost — do not parallelize `effect` against one depot

Recorded because it cost this study two and a half hours and would cost the next one the same. The three 48-seed `rdx effect` calls were first launched **concurrently**; they never completed. Three Julia processes against a single shared depot accumulated **116 seconds of CPU across 2 h 29 m of wall time — roughly 1% CPU** — i.e. blocked on depot/precompile contention rather than computing, leaving all three output JSONs at 0 bytes. Re-run **strictly sequentially** with identical commands, seeds and root seed, each call took **16-20 minutes at ~96% CPU** and completed normally. Budget ~1 hour of wall time for a three-arm, 48-seed sweep, and run the arms one at a time.
