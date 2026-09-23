# BRIEF — What are three in-licensed assets worth to a development group that is already full?

> The FRAME artifact (AGENT_WORKFLOW_CHARTER §3.1). Written in the operator's register — headcount, phases, term sheets — but every slot maps 1:1 onto the object model, so AUTHOR can compile it without a second interview. Phases: this brief → the nine `*.rdj.json` documents → [`DOCTOR.md`](DOCTOR.md) → [`FINDINGS.md`](FINDINGS.md). Provenance for every number quoted downstream is [`README.md`](README.md); the workflow is [`../../README.md`](../../README.md).

## The decision

There is a term sheet for three clinical-stage assets — one entering Phase II, two behind it — and a banker's risk-adjusted NPV of roughly \$740M for the package. Development is funded at 40 FTE, and the budget is ample: money is not the binding constraint here, people are. The question is not whether the banker's arithmetic is right. It is **what that package is worth to us, given that our own pipeline already saturates the group that would have to develop it — and whether bidding obliges us to fund the integration capacity as well.**

Those are two cheques and this study prices them separately: the assets, and the people who would have to develop them. The deliverable is therefore two numbers with error bars — the eNPV the three assets add at today's capacity, and the eNPV twelve more FTE add *given* the deal — plus the gap between the first of those and the same three assets valued against unlimited capacity. That gap is the headline: it is the part of the seller's valuation that exists on their side of the table and not on ours.

**What would make us act differently.** (1) A bid ceiling: bid at or below what the assets deliver *to us*, which equals the standalone valuation only if capacity is free. (2) A go/no-go plus a dose on the second cheque: fund integration FTE alongside the bid, or not, and if so how many. (3) If the capacity gap turns out to be small, the banker's frame is roughly right, we bid against it and skip the headcount request — a null here is worth having.

## Why this is not a spreadsheet question

Three assets and unlimited people is arithmetic: multiply per-asset value by three and discount. Three assets and 40 shared people is a queue. Whether the second acquired asset ever reaches Phase III depends on what our own Phase III programmes are holding that month, on how long each holding lasts, and on who the allocator staffs when both want the same ten people — none of which is recoverable from the assets' average values. The quantity we want is the *difference* between a valuation that ignores contention and one that does not, and that difference only exists if the contention is real.

Qualifying test (charter §2.4): **something binds** — 40 development FTE; **something takes time** — a Phase III holds ten of them for 30 months before returning any; **something is decided mid-flight** — when three acquired assets land in a full group, a rule has to decide which programme waits. Three for three; this is a simulation question.

## Slots

| Slot | Value | Object-model target |
|---|---|---|
| **Scarce resource** | 40 development FTE, reusable — a programme *occupies* its people for the cycle time and returns them intact | plain place `fte`, `init = 40`, held by an `upfront`/`conserved`/`block` arc |
| **Metered resource** | development budget, drawn every month a programme is in flight, out of a deliberately ample pool | plain place `cash`, `init = 6000` (\$M), `cost = 1.0`, `perstep`/`consumed` arc |
| **Work in progress** | assets by stage; inside a stage they are interchangeable, so no per-asset attributes are needed | plain places `own_p1`/`own_p2`/`own_p3`/`own_filed`/`own_market` and `acq_p2`/`acq_p3`/`acq_filed`/`acq_market` — **no structured tokens, no `host.jl`, the whole study is data** |
| **Starting state** | 4 preclinical, 4 Phase II, 3 Phase III of our own; nothing acquired | `init` on the owned stage places; the acquired places start at 0 |
| **Intake** | our own discovery keeps feeding the group, ≈4 new preclinical candidates a year | transition `discovery`, Poisson `rate = 0.35`/month, no input arc |
| **The process** | per stage: take the asset, hold FTE for the whole cycle, meter cash monthly, emit into the next stage with the stage's probability of success | seven `adv_*` transitions — own P1 (3 FTE, 18 mo, 0.55), own P2 (6, 24, 0.38), own P3 (10, 30, 0.62), own filed (2, 12, 0.90), acq P2 (6, 14, 0.48), acq P3 (10, 30, 0.65), acq filed (2, 12, 0.90) |
| **Contention rule** | late-stage work outranks early-stage work; the acquired Phase II asset is slotted ahead of our own Phase II but behind our own Phase III | `priority` 1.0 / 2.0 / 3.0 / 4.0 on the owned ladder, 2.5 on `adv_acq_p2` |
| **Outcome that carries value** | a launch — an in-licensed launch is worth more per unit than a home-grown one | `reward = 900.0` on `own_market`, `reward = 1400.0` on `acq_market`, credited on production |
| **Horizon and tick** | ten years, monthly | `meta.tspan = 120.0`, `meta.dt = 1.0` |
| **The lever** | one decision point at month 6: `acq_p2 += 3` (the bid) and/or `fte += 12` (the integration hire), in every combination | one appended `rules` block per arm, `id: deal_closes`, guard `t > 6.0`, `fire_mode: "once"` — so `diff model.rdj.json <arm>.rdj.json` *is* the intervention |
| **Headline metric** | expected NPV in \$M over the decade, monthly flows discounted at 8%/yr | `enpv`, `--scale` left at 1.0 because `cost`/`reward` are already \$M |

**Units, pinned here and restated in the report:** one tick is one month, all money is \$M, all capacity is FTE. `cash` is a \$6,000M pool at `cost` 1.0 per \$M drawn; results are \$M of eNPV discounted at 8%/yr over 12 periods per year.

## The arms

Nine documents, because the decision has two levers and the lever that matters needs a control. Baseline `model` (40 FTE, no deal); `arm_acquire` (bid only) — what the assets are worth to us; `arm_capacity` (+12 FTE, no deal) — **the control that prices the people alone**, without which "capacity is valuable" is unfalsifiable; `arm_both` (bid + 12 FTE) — the programme; `arm_capacity_4` / `arm_capacity_24` / `arm_both_24` — dose response; and `uncapped_baseline` / `uncapped_acquire` at `fte.init = 4000`, which is how the seller's standalone valuation is reconstructed *inside the same model*, with the assets, the pipeline and the economics identical and only the constraint removed. Every comparison is paired at 120 ensemble members from root seed 2026, member *k* seeded from `hash((2026, k))`, so both sides of every contrast see the same seed stream.

## What could go wrong (pre-registered, for DIAGNOSE to check)

- **FTE doesn't actually bind.** If 40 people absorb three more assets without anything queueing, there is no capacity study and the answer is "bid the banker's number". → `contention_exists`, `--expect-binds=fte`.
- **FTE declared consumed rather than conserved.** The default arc modality spends a resource; headcount declared that way drains to zero and never comes back, which would understate capacity all run and *overstate* the value of the integration hire — the exact number we are trying to price. → `modality_intent`, `--expect-reusable=fte`.
- **Cash binds instead of people.** Then the headline is a financing number wearing a capacity number's label, and the memo would tell the VP to fund the wrong thing. → `intended_resource_binds`.
- **The three `adv_acq_*` transitions are dead in the baseline.** They must be: the counterfactual arm never acquires anything, so there is nothing for them to advance. We expect `transitions_fire` to **fail on `model` and pass on `arm_both`**, and we expect to *explain* it rather than fix it. Pre-registering it is the difference between a known counterfactual and a discovered embarrassment.
- **Launch counters look like runaway growth.** `own_market` and `acq_market` only ever accumulate. → `--tally=own_market`, plus `acq_market` on the arms that have one.
- **Horizon too short to realize the outcome.** A 30-month Phase III started in year eight cannot launch inside the window. We expect `horizon_sufficient` to warn; that is scope, not defect, and its consequence is that every eNPV here is a **lower bound** — to be stated in FINDINGS, not buried.
- **The effect is inside the ensemble noise.** Especially the two capacity-only arms and any interaction between the levers, which is a difference of differences and therefore the noisiest thing we will compute. → `significant_2se` on `rdx effect`, and an interaction standard error built as a root-sum-square of independent arm sems rather than by differencing two paired effect sems, which double-counts the variance of any arm the two effects share. Every arm entering a derived contrast therefore needs its own ensemble, and every such row is to be recomputed from those envelopes by script before it is quoted anywhere — a mistranscribed summary row is invisible to both `validate` and `diagnose`.

## What "done" looks like

Two priced cheques with standard errors, a stated gap between the seller's frame and ours, conservative lower bounds at Δ − 2·se for anything we would put in a board paper, and every null reported as a null. Then hand to `rd-author`.
