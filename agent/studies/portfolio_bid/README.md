# `portfolio_bid` — what three in-licensed assets are worth to an organisation that is already full

The provenance file for this study. It records the question, the model, the nine arms, **every command that produced a number**, and the numbers those commands produced, so that any figure quoted in [`FINDINGS.md`](FINDINGS.md) or in the session capture can be traced back to a command and re-run. Nothing in this directory was computed by hand; nothing quoted anywhere in the study is absent from this file.

## The question

A VP of corporate development has a term sheet for three clinical-stage assets — one entering Phase II, two behind it — and a banker's valuation that says the package is worth roughly \$740M in risk-adjusted NPV. Development is funded at 40 FTE and the budget is not the binding constraint. The decision is not "is the banker's number right"; it is **"what is that package worth *to us*, given that our own pipeline already saturates the group that would have to develop it — and if we are going to bid, do we also have to fund the integration capacity?"** Those are two different cheques, and the study prices them separately.

## Units and conventions

| Quantity | Unit | Note |
|---|---|---|
| Tick | one month | `meta.dt = 1.0`, `meta.tspan = 120.0` → a ten-year horizon |
| Money | \$M | `cash` is a \$6,000M budget; place `cost`/`reward` are \$M per unit |
| Capacity | FTE | `fte` is a pool of 40 development FTE |
| Result metric | `enpv` | ledger flows discounted at **8%/yr**, 12 periods per year, `--scale` left at 1.0 (flows are already \$M) |
| Ensemble | 120 members | root seed **2026**; member *k* is seeded from `hash((root_seed, k))`, so paired arms share the seed stream |

A firing occupies `ceil(cycletime/dt)` ticks, so with `dt = 1.0` every cycle time here is already a whole number of months and nothing is lost to quantization.

## The model

Eleven places, eight transitions, twenty-nine arcs, **no structured tokens and no `host.jl`** — the entire study is data. `rdx validate` reports zero diagnostics on all nine documents.

Places: `fte` (40, the funded development group), `cash` (6,000, `cost` 1.0 per \$M drawn), the owned pipeline `own_p1` 4 / `own_p2` 4 / `own_p3` 3 / `own_filed` 0 / `own_market` 0 (`reward` 900.0), and the acquired pipeline `acq_p2` / `acq_p3` / `acq_filed` / `acq_market` 0 (`reward` 1400.0 — the in-licensed assets are worth more per launch than the home-grown ones).

Transitions — each takes its stage place, holds FTE for the duration, meters cash per month, and emits into the next stage with probability `prob_of_success`:

| Transition | Stage | FTE held | \$M/month | Cycle | PoS | Priority |
|---|---|---|---|---|---|---|
| `adv_own_p1` | own preclinical → Ph II | 3 | 0.6 | 18 mo | 0.55 | 1.0 |
| `adv_own_p2` | own Ph II → Ph III | 6 | 1.2 | 24 mo | 0.38 | 2.0 |
| `adv_own_p3` | own Ph III → filed | 10 | 2.5 | 30 mo | 0.62 | 3.0 |
| `adv_own_filed` | own filed → market | 2 | 0.8 | 12 mo | 0.90 | 4.0 |
| `adv_acq_p2` | acquired Ph II → Ph III | 6 | 1.6 | 14 mo | 0.48 | 2.5 |
| `adv_acq_p3` | acquired Ph III → filed | 10 | 3.0 | 30 mo | 0.65 | 3.0 |
| `adv_acq_filed` | acquired filed → market | 2 | 0.9 | 12 mo | 0.90 | 4.0 |
| `discovery` | → own preclinical | — | — | — | — | — |

The seven `adv_*` transitions carry `rate = 3.0` with `rate_mode = "deterministic"`: intent is *start everything you can afford to start*, so the queue that forms is caused by FTE scarcity and not by a sampling rate. `discovery` is Poisson at 0.35/month (≈4 new preclinical candidates a year) and is the only transition without an input arc. Priorities put late-stage work ahead of early-stage work under contention, with the acquired Phase II asset slotted at 2.5 — ahead of our own Phase II, behind our own Phase III.

Three arc modalities carry the whole resource story:

| Modality | `allocation` | `return` | `blocking` | Used for |
|---|---|---|---|---|
| held | `upfront` | `conserved` | `block` | `fte` — claimed in full at start, returned at completion, and a firing that cannot claim its people does not start |
| metered | `perstep` | `consumed` | `block` | `cash` — drawn every month the work is in flight |
| taken | `upfront` | `consumed` | `block` | the stage place — the asset leaves the stage it is advancing out of |

`fte` being `conserved` rather than `consumed` is the modelling claim that makes the study a capacity study and not a budget study: people come back, money does not.

## The nine arms

Every intervention lands through **one rule** at month 6 (`fire_mode: "once"`, guard `t > 6.0`), so `diff model.rdj.json <arm>.rdj.json` shows the entire lever as one appended `rules` block. The rule is named `deal_closes` in every arm that has one; in the capacity-only arms it is the same decision point carrying only the headcount half.

| Document | `fte.init` | Rule at month 6 | Role |
|---|---|---|---|
| `model.rdj.json` | 40 | — | baseline: our pipeline, our capacity, no deal |
| `arm_acquire.rdj.json` | 40 | `acq_p2 += 3` | bid, capacity unchanged — **what the assets are worth to us** |
| `arm_capacity.rdj.json` | 40 | `fte += 12` | capacity with no deal — the control that prices the people alone |
| `arm_both.rdj.json` | 40 | `acq_p2 += 3`, `fte += 12` | the combined programme |
| `arm_capacity_4.rdj.json` | 40 | `fte += 4` | dose response, low |
| `arm_capacity_24.rdj.json` | 40 | `fte += 24` | dose response, high (no deal) |
| `arm_both_24.rdj.json` | 40 | `acq_p2 += 3`, `fte += 24` | dose response, high (with the deal) |
| `uncapped_baseline.rdj.json` | 4000 | — | capacity removed |
| `uncapped_acquire.rdj.json` | 4000 | `acq_p2 += 3` | capacity removed, with the deal — **the standalone valuation the seller is quoting** |

`fte.init = 4000` is how a standalone per-asset valuation is reproduced inside the same model: the assets are identical, the pipeline is identical, and the only thing removed is the constraint that two programmes cannot use the same person at the same time.

## Every command

Three batteries, run strictly sequentially — concurrent Julia processes contend on one package depot and block rather than compute.

```bash
# 1. the paired comparisons. 120 members per arm, same seed stream on both sides.
S=agent/studies/portfolio_bid; R=$S/results
M='enpv,net,reward,cost,final:own_market,final:acq_market,max:own_p1,successes:adv_own_p2,min:fte'
for pair in \
  "effect_integration    arm_acquire       arm_both" \
  "effect_acquire        model             arm_acquire" \
  "effect_standalone     uncapped_baseline uncapped_acquire" \
  "effect_capacity       model             arm_capacity" \
  "effect_both           model             arm_both" \
  "effect_integration_24 arm_acquire       arm_both_24" \
  "effect_capacity_4     model             arm_capacity_4"
do
  set -- $pair
  julia --project=. agent/scripts/rdx.jl effect $S/$2.rdj.json $S/$3.rdj.json \
      --seeds=120 --root-seed=2026 --metric=$M > $R/$1.json
done

# 2. arm levels and the per-member spread the error bars summarise. Every arm that appears in a
#    derived contrast needs one, because a contrast's se is built from arm sems, not from effect ses.
for a in model arm_acquire arm_capacity arm_both arm_both_24 uncapped_baseline uncapped_acquire; do
  julia --project=. agent/scripts/rdx.jl ensemble $S/$a.rdj.json --seeds=120 --root-seed=2026 \
      --per-seed --metric='enpv,net,final:acq_market,final:own_market,max:own_p1,min:fte' > $R/ens_$a.json
done

# 3. the two health reports, three single-seed traces, and the diagrams.
julia --project=. agent/scripts/rdx.jl diagnose $S/model.rdj.json --seed=1 \
    --expect-binds=fte --expect-reusable=fte --tally=own_market > $R/diagnose_model.json
julia --project=. agent/scripts/rdx.jl diagnose $S/arm_both.rdj.json --seed=1 \
    --expect-binds=fte --expect-reusable=fte --tally=own_market,acq_market > $R/diagnose_arm_both.json
T='fte,own_p1,own_p2,acq_p2,acq_p3,own_market,acq_market,cash'
for a in model arm_acquire arm_both; do
  julia --project=. agent/scripts/rdx.jl run $S/$a.rdj.json --seed=1 --trace=$T > $R/trace_$a.json
done
julia --project=. agent/scripts/rdx.jl draw $S/arm_acquire.rdj.json --exec-map --seed=1 --out=$S/figures/exec_map.svg
julia --project=. agent/scripts/rdx.jl draw $S/arm_acquire.rdj.json --out=$S/figures/net.svg
```

Then the audit and the four data figures, both read straight out of the envelopes above:

```bash
julia --project=. agent/studies/portfolio_bid/numbers.jl   # re-derives every table below
julia --project=. agent/studies/portfolio_bid/figures.jl   # redraws the four data figures
```

`numbers.jl` is the check on this file: it prints the arm table, the paired effects, the identity that justifies the derived contrasts, the contrasts themselves, the mechanism rows, the queue gauge and the traces, all recomputed from `results/`. **If its output and the tables below disagree, the tables are wrong.** That is not hypothetical — see the note under the derived contrasts.

## The numbers

### Measured paired effects — Δ eNPV over ten years, \$M

| Comparison | Arms | Δ eNPV | se | se from zero | Verdict |
|---|---|---|---|---|---|
| standalone | `uncapped_baseline` → `uncapped_acquire` | **+744** | 144 | 5.17 | clears |
| realized | `model` → `arm_acquire` | **+241** | 90 | 2.69 | clears |
| capacity alone, +12 FTE | `model` → `arm_capacity` | −8 | 91 | 0.09 | **null** |
| capacity alone, +4 FTE | `model` → `arm_capacity_4` | +117 | 90 | 1.31 | **null** |
| integration, +12 FTE given the deal | `arm_acquire` → `arm_both` | **+370** | 106 | 3.48 | clears |
| integration, +24 FTE given the deal | `arm_acquire` → `arm_both_24` | +306 | 116 | 2.64 | clears |
| programme | `model` → `arm_both` | **+611** | 101 | 6.07 | clears |

Arm levels — the eNPV mean of each arm and the standard error of that mean, over the same 120 members. This is the table every derived contrast is built from, because a paired `rdx effect` se **is** the root-sum-square of its two arms' sems (verified to three decimals on all four effects — see `numbers.jl`):

| Arm | eNPV mean | sem |
|---|---|---|
| `model` | 985.0 | 58.6 |
| `arm_acquire` | 1225.9 | 67.7 |
| `arm_capacity` | 977.0 | 69.0 |
| `arm_both` | 1595.9 | 81.9 |
| `arm_both_24` | 1532.3 | 94.2 |
| `uncapped_baseline` | 1269.9 | 84.4 |
| `uncapped_acquire` | 2013.8 | 116.3 |

Conservative lower bounds (Δ − 2·se): the deal alone **≥ \$62M**, the integration capacity given the deal **≥ \$158M**, the combined programme **≥ \$409M**.

### Derived contrasts

Each is a linear combination of arm means, so its standard error is the root-sum-square of the sems of the arms it combines — the same unpaired estimator `RD.treatment_effect` uses (`src/analysis.jl:401`: the two arms desync their RNG streams after the first divergent draw, so there is no pairing left to exploit). Every row below is recomputed from the arm table by `numbers.jl`; none is transcribed.

| Derived quantity | Value | se | se from zero | Verdict |
|---|---|---|---|---|
| **stranded value** = standalone Δ − realized Δ | +503 | 169 | 2.97 | clears |
| **stranded launches** = `uncapped_acquire` − `arm_acquire`, `final:acq_market` | +0.383 | 0.094 | 4.06 | clears |
| **interaction** = (`arm_both`−`arm_acquire`) − (`arm_capacity`−`model`) | **+378** | 140 | 2.71 | clears |
| **dose response** = (`arm_both`−`arm_acquire`) − (`arm_both_24`−`arm_acquire`) | +64 | 125 | 0.51 | **null** |

The interaction is the contrast that says capacity and the deal are **complements**: twelve FTE are worth +378 ± 140 more when you own the three assets than when you do not, which is the same statement as "−8 alone, +370 given the deal" with one error bar instead of two. **Build a contrast of contrasts from the arm sems, never from the sems of the effects it combines** — differencing two effect standard errors is valid only when the two effects share no arm. The interaction happens to be safe (four distinct arms, and the naive se agrees to 0.5%), but the dose-response row is not: +370 ± 106 and +306 ± 116 are both measured against `arm_acquire`, so the shared arm's variance has to come out once. Doing it wrong there gives ±157 instead of ±125 — conservative, and in this case harmless, since 0.41 se and 0.51 se are both nulls and the cheaper dose wins either way.

This row is also the study's caught error. Its first hand-written version carried a value from a different contrast and a standard error copied from the line above, and reported a clear result as a null. Nothing in `validate` or `diagnose` can see a mistake like that — it was in a transcribed summary, not in the model — and what caught it was re-deriving the table from the envelopes before publication. That is what `numbers.jl` is for, and it is the reason this file, not any document quoting it, is the place numbers live.

Launches of acquired assets over the decade (`final:acq_market`, 3 assets bought): `arm_acquire` 0.533 ± 0.061, `arm_both` 0.975 ± 0.072, `uncapped_acquire` 0.917 ± 0.072.

### The mechanism — crowd-out, and the control that proves it is the ceiling

| Comparison | `final:own_market` (own launches) | `successes:adv_own_p2` (own Ph II completions) |
|---|---|---|
| `model` → `arm_acquire` (deal, 40 FTE) | 2.842 → 2.542, Δ **−0.300 ± 0.144** (2.08 se, clears) | 4.108 → 3.392, Δ **−0.717 ± 0.171** (4.19 se, clears) |
| `uncapped_baseline` → `uncapped_acquire` (deal, no ceiling) | 4.575 → 4.708, Δ +0.133 ± 0.243 (null) | 8.383 → 8.558, Δ +0.175 ± 0.360 (null) |
| `model` → `arm_both` (deal + 12 FTE) | 2.842 → 2.808, Δ −0.033 ± 0.147 (null) | 4.108 → 4.825, Δ +0.717 ± 0.179 (4.01 se, gains) |

Taking on the three assets at today's capacity costs a third of an own-pipeline launch and three quarters of an own Phase II completion. Remove the ceiling and both effects vanish; fund twelve more FTE and the crowd-out is gone while the own pipeline finishes *more* Phase II work than before. Nothing changed between the first and second rows except `fte.init`.

`max:own_p1` is the queue gauge — the deepest the preclinical backlog gets: 25.7 baseline, 27.6 with the deal, 20.8 with the deal plus 12 FTE, 16.1 with 12 FTE and no deal, 14.9 with 24 FTE and the deal, and 4.0 with the ceiling removed. A preclinical backlog that peaks near 26 in a pipeline that opened with 4 is the shape of the constraint; with the ceiling removed it never exceeds 4.

### The traces — seed 1, the concrete mechanism

The three assets land in `acq_p2` at **month 8**. Then:

| Arm | Acquired assets start Phase II | FTE free ≤ 1 |
|---|---|---|
| `model` (no deal, 40 FTE) | — | 96 of 122 ticks (79%) |
| `arm_acquire` (deal, 40 FTE) | months **25, 31, 39** | 81 of 122 (66%) |
| `arm_both` (deal, +12 FTE) | months **9, 9, 25** | 52 of 122 (43%) |

The last of the three assets sits on the shelf for 31 months at today's capacity. Fund twelve more people and two of them start the month after the deal closes. Note that the *baseline* is the most saturated arm of the three: with no acquired work competing, the group is packed into its own pipeline, and a 6-FTE Phase II start cannot fit into the gaps a 40-FTE ceiling leaves. That is why the argument in the report is made from queue-wait and launch counts, not from a saturation percentage — the percentage moves the wrong way for a fragmentation reason that has nothing to do with the decision.

### `rdx diagnose`

`model` — 7 pass, 1 warn, **1 fail**. The failure is `transitions_fire`: `adv_acq_p2`, `adv_acq_p3` and `adv_acq_filed` never spawn, because in the *baseline* there are no acquired assets to advance. That is the counterfactual arm behaving correctly, and it is the one case where the check must be explained rather than fixed — the same check passes on `arm_both` (8 pass, 1 warn, 0 fail, all 8 transitions firing). A study that reports "diagnose is green" without saying which arm it ran on is hiding this.

Both arms warn `horizon_sufficient`: 7 (baseline) and 9 (`arm_both`) firings still in flight at month 120, with pools still trending in the final quarter. Every eNPV here is therefore **horizon-truncated and a lower bound** — work started before month 120 whose launch lands after it contributes nothing. This is also why `arm_both_24` measures +306 rather than more: 24 extra FTE start still more work, but the extra launches fall outside the window where the ledger can book them.

## What this study does not support

- **That 24 FTE beats 12.** +306 vs +370, a difference of 64 against a se of 125 (0.51 se). Indistinguishable, so the cheaper dose is the one to fund — and note that this is a *floor* argument, not a ceiling one: see the horizon warning below.
- **A dose-response curve.** Three points (+117 at 4 FTE, −8 at 12, +306/+370 with the deal) with sems near 90 do not determine a shape. The study supports "fund 12 with the deal", not "12 is optimal".
- **Anything about the horizon beyond month 120.** See the `horizon_sufficient` warning above; all figures are lower bounds.
- **A valuation of the assets themselves.** `reward` 1400/launch and the stage PoS values are illustrative inputs, not a diligence model. What the study measures is the *gap* between a valuation that ignores capacity and one that does not, and that gap is a property of the constraint, not of the price.
- **A number for a different organisation.** The 40-FTE ceiling, the stage FTE loads and the priority order are the inputs that produce the whole result. Change them and every figure moves.

## Figures

| File | Written by |
|---|---|
| `figures/exec_map.svg` | `rdx draw arm_acquire.rdj.json --exec-map --seed=1` — the network with contended places filled gold (9 of 11: `fte` and every stage place except `own_p1` and `cash`) |
| `figures/net.svg` | `rdx draw arm_acquire.rdj.json` — the same network, undecorated |
| `figures/capacity_trace.svg` | `figures.jl` from `results/trace_arm_acquire.json` + `trace_arm_both.json` |
| `figures/levers.svg` | `figures.jl` from the seven `results/effect_*.json` |
| `figures/spread.svg` | `figures.jl` from `results/ens_model.json` + `ens_arm_both.json` (`--per-seed`) |
| `figures/bridge.svg` | `figures.jl` from `effect_standalone`, `effect_acquire`, `effect_integration`, `effect_both` |

`figures.jl` computes no simulation result: every value it draws is read out of an envelope written by a command above, so a figure cannot drift away from the number it illustrates.

## The four phases

[`BRIEF.md`](BRIEF.md) frames the decision → the nine `*.rdj.json` documents are the model → [`DOCTOR.md`](DOCTOR.md) is the health report → [`FINDINGS.md`](FINDINGS.md) is the briefing. A rendered capture of the conversation that produces them is in [`../../sessions/portfolio_bid.html`](../../sessions/portfolio_bid.html).
