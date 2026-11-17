# Pharma R&D digital twin: NocoBase → DyVE (ReactiveDynamics + AlgebraicAgents)

This folder turns a NocoBase prototype of a pharma R&D digital twin into a runnable DyVE model. The prototype has three parts, none of them included here: a written description of its collection schema and five workflows, the "Pharma R&D Digital Twin Prototype" deck at tick 42 (not public), and the "Corvane Twin" coordination mockup (for Corvane Therapeutics, a fictional company). The model asks the deck's own question: *what is a closed-loop portfolio policy (an automatic tox kill switch plus a monthly re-plan of wet-lab FTE) worth against a static annual budget, over 24 months?*

**All data are synthetic.** The six programs, their codes, targets, values and events come from the prototype's AI-generated synthetic data for a fictional company. Nothing here describes a real company, portfolio, compound, patient or market, and no result is a forecast for one. Drug targets are shown as generic labels (Target A to F). Every value the prototype does not state is an assumption, listed at the end of this file and in the snapshot under `_assumed`.

No ReactiveDynamics or AlgebraicAgents source is changed. Everything runs on the public API plus the documented host seam: one token kind and two registry functions.

## Run it

```bash
julia --project=demo/nocobase_pharma_twin -e 'using Pkg; Pkg.instantiate()'     # first time only
julia --project=demo/nocobase_pharma_twin demo/nocobase_pharma_twin/run_twin.jl  # ~2 min, 40 seeds per arm
julia --project=demo/nocobase_pharma_twin demo/nocobase_pharma_twin/run_twin.jl --nseed=10   # quick look
julia --project=demo/nocobase_pharma_twin demo/nocobase_pharma_twin/run_twin.jl --snapshot=whatif.json --out=runs/whatif   # a what-if
```

The script prints six sections and writes everything into `out/` (or `--out`), regenerated on every run. Same seed, same numbers, on Julia 1.12: ensemble seeds come from Julia's `hash`, which can change between Julia versions. The four model drawings need Graphviz `dot` on the PATH and are skipped without it.

`handoff/` holds a standalone copy for anyone who wants to run their own views with an agent: `handoff/build_handoff.sh` assembles `handoff/dist/pharma_twin/` (and a zip) with the talk page and its PDF (`handoff/make_pdf.sh`), the model code, the snapshot, the baseline `out/`, a `Project.toml` against the registered ReactiveDynamics 0.3 and AlgebraicAgents 0.4, the agent guide `AGENTS.md`, and a paste-ready `HANDOFF_PROMPT.md`.

| File | What it is |
|---|---|
| `data/nocobase_snapshot.json` | The input: the deck's collections at tick 42, in NocoBase shape. Every value the deck does not state is listed under `_assumed`. |
| `twin.jl` | The bridge. `compile_model` (snapshot → model document), the host seam (`ProgramToken`, `rebalance`, `sync_status`), the `ClinOpsAgent`, and the metrics. |
| `run_twin.jl` | The driver: §1 read and reconcile against the deck, §2 compile and validate, §3 one coupled run, §4 ensembles, §5 NocoBase write-back, §6 figures. |
| `out/model_<arm>.rdj.json` | The compiled, eval-free model document per scenario arm. It is plain JSON, so it can be validated, diffed, and stored back in NocoBase, and the agent layer's `rdx` can read it. |
| `out/scenarios.json`, `simulation_ticks.json`, `simulation_events.json`, `resource_allocations_plan.json` | The write-back, in the deck's table shapes, every row scoped by `scenario_id` (the Corvane convention). |
| `out/net_*.svg`, `out/aa_wiring.svg` | The model drawn: one program's Petri-net slice, the full net, the AlgebraicAgents hierarchy with its wires. |
| `out/ensemble_per_seed.json` | One row per arm × seed (160 rows): value, spend, net, realized eNPV, and per program the kill month, slip and readout. Rows sharing a `seed_index` are the same future under different policies. |
| `out/chart_*.svg` | Dashboard charts (lever effects, paired seeds, wet-lab and shared-pool utilization). |
| `talk.html` | The shareable talk page, in talk order: the answer and the decisions it informs, what limits each program, the model, what each lever is worth, how NocoBase records become the model, what the model leaves out, and next steps. Technical detail sits in collapsed sections. Its numbers are this run's 40-seed results; `handoff/build_handoff.sh` writes a copy that opens directly in a browser. |

## The result

Ensembles of 40 seeds per arm over 24 months. The arms share seeds, so each comparison is paired. Net value is risk-adjusted eNPV minus spend, both in $M.

| Arm | Risk-adj. eNPV | Spend | Net | Δ net vs. static |
|---|---|---|---|---|
| Static annual plan | 1824.1 | 466.7 | 1357.4 | — |
| Tox kill switch only | 1826.0 | 430.2 | 1395.7 | +38.3 ± 1.7 |
| Monthly re-plan only | 1837.1 | 465.7 | 1371.5 | +14.1 ± 1.5 |
| Closed loop (both) | 1841.8 | 431.1 | 1410.7 | **+53.3 ± 1.7** |

- **The closed loop is worth about +$53M net over 24 months, +3.9% on the static arm.** The prototype deck's first estimate was +14.8%. Two-thirds of the gain is spend avoided: PRG-402's Phase 1 is stopped at month 1 instead of running until the annual review at month 12. The remaining third is eNPV protected by less slip on the programs that get the freed FTE. (Split by lever instead, the kill switch is 72% of the gain: the two splits divide the same $53.3M differently.)
- **The two levers barely interact** (+0.9). Each can be justified on its own; the kill switch is the larger.
- **The number depends on one assumption above all:** how much eNPV a month of slip destroys (patent-life erosion). At 0%/month the gain is +$36M, which is the spend saving alone. At 1.5% (the default) it is +$53M; at 5% it is +$95M (+7%). Even at 5% erosion the model does not reach +14.8%. That figure would need a much larger erosion rate, or value sources the prototype does not yet model.
- **Freeing FTE is not enough for PRG-518.** The controller gives the cell-therapy program its full 16 FTE from month 0. Yet in the logged run it still stalls twice (months 6 and 15): kits run out while patients and FTE are available. With the assumed 15-kit batch and 3-month lead time, one Suite B batch line cannot keep a 5-kit/month study supplied. Its binding constraint is manufacturing, not people, which is the deck's own "GMP Slot Contention" warning. Whether that would hold for a real program of this kind depends on the batch size and lead time, which the deck does not give.
- **One simulated future shows nothing.** Realized eNPV of a single run (which readouts happened to pass) has a standard deviation of about $600M across seeds, 11× the effect. The deck's single-trajectory dashboard cannot show the policy effect; the paired ensemble can.

## How NocoBase maps onto DyVE, and how faithful it is

Verdicts: **faithful** means the NocoBase semantics are reproduced. **Approximated** means the structure is kept but the mechanism is simplified (stated). **Not mapped** means it is out of scope for a simulation engine, or not built here.

### Collections (the prototype's schema description and the deck's schema slide)

| NocoBase | DyVE construct | Verdict |
|---|---|---|
| `Programs` / `therapeutics_programs` row | A `ProgramToken`, a first-class AlgebraicAgents agent with code, modality, area, stage, `pos`, `stage_pos`, `enpv`, living in the structured place `Program` | Faithful |
| `Programs.current_stage` | A token attribute (the stage is the token's lifecycle state); a passed readout writes `readout_t` and the next stage is reported | Faithful for one stage transition. The 24-month horizon covers the readout of the current stage only; follow-on stages are not simulated (see gaps) |
| `Programs.pts_estimate` | The cumulative `pos`; the per-readout `stage_pos` is the transition's terminal success probability | Approximated: the deck gives only cumulative PoS, so `stage_pos` is assumed |
| `Programs.npv_projected` / `enpv_usd_m` | Token attribute; risk-adjusted value = `enpv × pos × (1 − ε·slip)` | Faithful (the erosion term is WF4's patent-life erosion, made explicit) |
| `Studies` row (the in-flight study) | A transition `study_<p>`, spawned once, bound to its program by `@select(Program, code == …)`, cycle time = remaining months × Gamma noise | Faithful |
| `Studies.status` (Planned/Active/Readout/Completed/Halted) | Planned = the enabler token `start_<p>`; Active = a firing in flight; Halted = the `go_<p>` gate at 0 (the firing freezes at 0% speed); Completed = finished; a failed readout soft-retires the token (= `Terminated`) | Faithful |
| `AssayReadouts.passed_gate` | The Binomial terminal draw of the study firing | Approximated: the readout is a single pass/fail, with no `observed_value` vs `threshold_value` |
| `Compounds` / `compounds_candidates` | `herg_ic50_um` becomes a param read by the kill rule; `kits_on_hand` becomes a place `kits_<p>` consumed per month | Faithful for the fields used; `smiles`, `ic50_nm`, `is_lead` are carried in the snapshot but drive nothing |
| `Departments`, `FTE_Pools` | Wet-lab FTE as per-program quota places `wet_<p>` refreshed from the plan (`q_<p>`); Clinical Ops FTE as a shared place `clin` | Faithful. `burn_rate_monthly` is the place's ledger cost |
| `StudyAllocations` / `resource_allocations` (M:N junction) | Arcs: study → pool with a per-month multiplicity. The controller rewrites the wet-lab column at every re-plan (`out/resource_allocations_plan.json`) | Faithful. The per-program split is assumed; it reproduces the deck's pool totals exactly |
| `resources_capacity` / `ComputeClusters` (mass spec, GMP Suite B, vivarium, GPU) | Shared places refreshed to monthly capacity and drawn `perstep`; contention is resolved by RD's priority-weighted allocator; utilization is read off the allocation log | Faithful for capacity and contention. Queue backlogs (14 runs, 6 lots) are not seeded |
| `BudgetLedgers` | Every place carries a unit cost; the engine's ledger gives spend per tick and per program (`program_ledger`), and the two reconcile exactly | Faithful, except that the ledger is monthly and not split into capex/opex/committed |
| `SimulationTicks` / `simulation_events` | Written back by the driver (§5) with `scenario_id`, `tick_index`, `sim_date`, `portfolio_enpv`, `snapshot_metrics` | Faithful as an output |
| Corvane: `trials`, `sites`, enrollment | `ClinOpsAgent` draws `Poisson(site_count × 1.8)` per month (the deck's own formula) and feeds a patient place per trial over an AlgebraicAgents wire | Faithful at trial level; site-level detail is not modeled |
| Corvane: `inventory_lots`, `purchase_orders`, `suppliers` | Kits per program plus a GMP batch transition (3-month lead time, occupies Suite B) | Approximated: no lots, expiry or FEFO |
| Corvane: `gl_transactions`, month-end close | Per-tick ledger rows | Approximated: there are no accrual or invoice-lag mechanics |
| Corvane: `scenarios` (self-FK, params jsonb) | Arms, written as scenario rows whose parent is the baseline | Faithful as an output |

### Workflows (the prototype's workflow description)

| NocoBase workflow | DyVE construct | Verdict |
|---|---|---|
| WF1 Material and capacity guardrail (block the status change) | Native: a study cannot run faster than its inputs allow; the allocator slows or stalls it, and the shortfall is logged as saturation < 1 | Faithful in effect. RD slows the study rather than refusing the status change, and no "Delayed" status row is written |
| WF2 PoTS and stage gate, pass branch | Terminal success of the study firing, then `readout_t` is written on the token and the value re-marked to `pos / stage_pos` | Approximated: a fixed re-mark, not a Bayesian update from `raw_payload` |
| WF2 tox branch with human approval (Terminate / Pivot) | The `kill_<p>` rule: hERG IC50 below 1.0 µM halts the study. The human approval is replaced by a policy: automatic (`kill_at = 0`) or acted on at the annual review (`kill_at = 12`) | Approximated: there is no human in the loop and no pivot to a backup compound |
| WF3 Inventory and budget burn, plus auto-reorder | Kits consumed per month at the study's speed; below the reorder point a GMP batch transition starts, occupies Suite B and restocks after 3 months | Faithful. The batch is a pre-declared transition rather than a new `Studies` row |
| WF4 Recruitment lag and shock propagation | Accrual below plan starves the study's patient input and slows it; slip feeds eNPV erosion; an "Accrual lag" event is written | Faithful |
| WF5 Portfolio rebalancing (AI / solver memo) | The `replan` rule invokes the host `rebalance`: wet-lab FTE goes to live studies by risk-adjusted eNPV per FTE, and allocator priorities are set the same way | Approximated: a greedy rule, not an MPC or LLM memo. It is optimal for the linear erosion objective used here |

### The Corvane coordination idea

The Corvane deck's architecture ("coordination, not universal simulation": several subsystem models talking through typed contracts) is exactly what AlgebraicAgents is for. Here the portfolio net is one RD agent and clinical operations is a separate AA agent with its own RNG, joined by four wires (`out/aa_wiring.svg`). The `enroll_to_supply` contract shape (range [0, 40] patients/month, reject and log on violation) is checked by the producer before it publishes. Supply, finance and workforce stay inside the RD net in this version; each could be split out the same way.

## What RD/AA do not cover (gaps)

These are the things the decks ask for that the engine does not provide. Each lists what this demo does instead.

1. **System of record, UI, roles, NL→SQL, LLM memos.** These are NocoBase's job, not a simulator's. This demo reads a NocoBase snapshot and writes NocoBase-shaped rows back. The "Ask the twin" Q&A and the AI trade-off memo would sit on those rows.
2. **Human-in-the-loop decisions.** RD runs a policy; it cannot pause a run for an approval. Model each human choice as a scenario arm (here: kill now vs. at the annual review), or step the AA root manually and stop at a decision.
3. **Bayesian PoTS update from assay payloads.** There is no inference engine in RD. A readout is a pass/fail draw with a fixed re-mark. A Bayesian or SciML updater would be another AA agent writing program parameters.
4. **Contracts on wires** (units, cadence, staleness, validity, ownership, the Corvane `check_contract`). AA wires carry values only. The range check here is host code in the producer agent; units and staleness are not checked.
5. **Learned residuals and algebraic constraint projection** (Corvane's ẋ = f_known + B·r_φ, then projection onto 0 = g(x,u)). RD has neither. Hard capacity constraints (ΣFTE ≤ headcount, stock ≥ 0) are native, because the allocator never over-allocates. A learned residual would be an AA agent (for example SciML) wired in.
6. **An MPC optimizer.** Not in RD. Re-planning is a registry function called by a rule, so a real optimizer (for example a JuMP model) can be dropped into the same seam unchanged.
7. **Queueing turnaround that grows nonlinearly with utilization** (the deck's Weibull lab queues). RD's contention slows work in proportion to the shortfall, with no queueing-theory delay. It could be approximated with a cycle time that reads utilization; that is not done here.
8. **New studies created mid-run.** RD has no live "append a transition" API yet (ADR 0004 growth half, [issue #27](https://github.com/Merck/ReactiveDynamics.jl/issues/27)). WF3's auto-created GMP study is therefore a pre-declared batch transition, and the next stage after a readout is not simulated within the horizon.
9. **Calendar time.** RD ticks are unitless months here. Dates, fiscal quarters and month-end closes are mapped by the driver.
10. **Topological producer-before-consumer order within a tick** (Corvane loop step ①). AA couples agents with a one-tick lag (Jacobi), so accrual drawn in month t reaches the net in month t+1.

## Engine findings from building this

Building the twin from JSON surfaced engine issues. None blocked the work: each has a workaround in `twin.jl`, and none touches `src/`. They are candidates for issues on this repository.

- **F1 `validate` accepts a bare number inside a `call`'s `args`, but the loader throws** (`node_from_dict(::Int64)`). Workaround: wrap every literal in a `const` node.
- **F2 `validate` accepts a string `const`, but the loader throws** (`Float64(::String)`). Workaround: use Symbol constants (`"symbol": true`) for codes.
- **F3 `transitions[].name` does not default to `id` on load** (`src/serialize.jl:166`, contrary to `agent/reference/model_schema.md`). `transName` stays `missing`, so `deactivate` by name and transition names in logs fail quietly. Workaround: always emit `name`.
- **F4 `RD.ensemble` fails on a member that is still entangled in a coupled AA hierarchy** (`KeyError: member_1`). Workaround: `disentangle!(rd)` after the coupled run, before returning it to `ensemble`.
- **F5 A rate above `capacity` is queued in a hidden spawn backlog** (`transToSpawn`) that outlives the condition that raised it. A conditional reorder ("kits < reorder point") kept firing at 395 kits. Workaround: gate each transition's rate on its own enabler place (`start_<p>`, `slot_<p>`).
- **F6 Genesis is not predicate-aware** (known). A `@select` transition spawns against the count of all tokens of the kind, and a spawn with no matching token still consumes resources and books unattributed cost. The enabler places above prevent it.
- **F7 A transition's `probability` cannot read a bound token's field**, so per-program PoS needs one transition per program. The table-driven compiler makes that cheap.
- **F8 Per-tick expressions are re-compiled every tick** (`wrap_fun`, about 0.5 ms each). `cache_compiles!` memoizes the compiler on the problem object, which makes a run about 10× faster and leaves the trajectory byte-identical; the driver checks this in §3.
- **F9 `PopulationEntry` attributes are evaluated as expressions**, so a bare Symbol is read as a variable name. Workaround: construct tokens directly and pass them to `population=`.

Already-known agent-layer gaps also apply: `from_json_model` ignores `d["population"]` (W14), and `pre_action`/`post_action` are dropped on import (W30); both are tracked in `spec/AGENT_WORKFLOW_CHARTER.md` §12. That is why all decisions here are `rules[]`.

## Assumptions (all in `data/nocobase_snapshot.json` → `_assumed`)

- **Per-readout success `stage_pos`:** from typical phase-transition rates, because the deck gives only cumulative PoS.
- **Months remaining per study:** not in the deck.
- **Per-program resource split:** chosen to reproduce the deck's pool totals exactly (48/50, 19/20, 35/40, 62/100, ClinOps 78%).
- **CRO cost and unit costs:** chosen so the tick-42 burn is $142.5M per quarter, against the deck's $142M QTD.
- **Enrollment and site counts:** chosen so PRG-104 accrues below 0.72 of plan (0.65 on average), as in event T-041.
- **Kit stocks:** chosen so PRG-518 starts below its reorder point, as in event T-039.
- **PRG-402's hERG IC50:** 0.8 µM, below the deck's 1.0 kill threshold (event T-042).
- **Patent-life erosion:** 1.5% of eNPV per month of slip (`valuation.erosion_per_month`), with the full sensitivity in §4. It is applied when results are valued, so it does not change the simulation.
