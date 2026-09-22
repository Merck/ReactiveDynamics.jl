---
name: rd-doctor
description: Phase 3 of an rd-study — prove a validate-clean ReactiveDynamics model is actually plausible before any number is quoted. Runs nine checks (transitions fire, contention exists, the intended resource binds, modality intent, horizon sufficiency, bounded pools, ledger wired, stochasticity) against the brief's pre-registered intent, and records the verdict in DOCTOR.md. Use after rd-author, before rd-report.
---

# rd-doctor — the plausibility gate

`validate` proved the document is well-formed. It proved nothing about whether the model is *about the user's system*. **A schema-clean model can be confident nonsense**, and the failure mode of this entire workflow is a beautifully formatted wrong answer. This phase is the gate that stands between the two.

Your output is `DOCTOR.md`: every line machine-produced evidence, not an opinion. `agent/studies/marginal_scientist/DOCTOR.md` is the exemplar.

## Run it

```bash
julia --project=. agent/scripts/rdx.jl diagnose <model> --host=<host.jl> --seed=1 \
  --expect-binds=scientist --expect-reusable=scientist --tally=launch
```

The `--expect-*` flags come from the brief's pre-registered failure modes — **do not invent them here**. A check whose expectation was written before the run existed is evidence; one you chose after seeing the output is decoration. If the brief didn't pre-register them, go back and fix the brief, then run.

Exit `0` = clear to report. Exit `1` = at least one failing check; **you may not quote a number**. Diagnose every arm, not just the baseline — an arm that changes the lever can change which resource binds.

## The nine checks

| Check | Severity | Asks |
|---|---|---|
| `transitions_fire` | fail | did every transition spawn at least once? |
| `completions_exist` | fail | did anything finish? |
| `contention_exists` | fail | **does anything actually run out?** |
| `intended_resource_binds` | fail | is the scarce thing the one the brief named? |
| `modality_intent` | fail | are the resources the brief called reusable declared `conserved`? |
| `horizon_sufficient` | warn | is anything still in flight at the end; are pools still trending? |
| `bounded_pools` | warn | is anything growing without limit that isn't a declared tally? |
| `ledger_wired` | warn | do reward and cost both move? |
| `stochasticity` | warn | is any outcome actually uncertain? |

### `contention_exists` is the one that matters

If nothing binds, **there is no study**. A model where every transition always gets what it asks for is an accounting identity with extra steps, and the marginal value it reports is just the accounting cost. When this check fails, do not tune the model until it passes — go back to the user and report the finding: *four scientists are enough; the fifth is worth their loaded cost and no more.* That is a real, valuable, publishable answer, and manufacturing scarcity to avoid delivering it is the worst thing you can do in this workflow.

### `transitions_fire` failures are usually one of three things

A dead transition means the model has a part that never ran, so any conclusion about it is vacuous. In order of likelihood: an lhs arc asks for a resource that is never available (a pool at 0, or drained because it was declared `consumed` when it should be `conserved`); a predicate never matches (check `"symbol": true` on every `:Symbol` constant — `"Onc"` will never equal `:Onc`); the rate is zero or the transition was never given one.

### `intended_resource_binds` protects the headline

If `budget` starves when the brief said `scientist` would, the number you are about to report is a *financing* number wearing a hiring number's label — and the memo would tell the user to do the wrong thing. This is not a modeling nicety; it is the check that keeps the answer about the question.

## Reading contention correctly

Saturation can be 1.0 while contention is completely real. A `conserved` resource rations by **refusing admission** — work that cannot get its people never starts, so everything that *did* start was fully served. The evidence there is the pool trough at 0 plus a growing queue, not fractional saturation. A `perstep` resource is the one that shows up as fractional saturation, because running work gets under-served. Both signals are checked and it only takes one; if you find yourself explaining why saturation looks fine, this is usually why. Say it in DOCTOR.md rather than leaving the reader to wonder.

## Warnings: carry them, don't clear them

A warning is usually a **scope statement**, and the right response is to state its consequence in the report, not to tune it away:

- **`horizon_sufficient`** — work in flight at the horizon means value not yet realized, so the headline is a **lower bound**. That is often exactly right: a hiring decision is made over the planning horizon, not over steady state. Extending `tspan` would raise the number; the decision does not need it to. Say that.
- **`bounded_pools`** — a cumulative counter (launches to date, total spend) legitimately only grows; declare it with `--tally=<place>`. Growth in a *queue* is different: that is the constraint binding, and it is evidence, not noise.
- **`ledger_wired`** — if reward never moves, no value is being realized and an NPV metric will be meaningless. Check the `reward` is on the **product** place and credited on production, not on the intake stock; rewarding intake pays for starting work and makes failure free.
- **`stochasticity`** — every `prob_of_success` at 1.0 means an error bar that reflects only arrival timing. Not wrong, but the ensemble spread means less than it looks like; say so.

Suppressing a warning to get a clean sheet is the one thing this phase cannot survive. Every warning either becomes a caveat in FINDINGS or gets fixed in AUTHOR — there is no third option.

## Look past the checks

The nine are a floor, not a ceiling. Also do this, by hand, once:

```bash
julia --project=. agent/scripts/rdx.jl run <model> --host=<host.jl> --seed=1
```

Read the spawn and completion counts and ask whether the *magnitudes* are believable in the user's world. Forty trials in four years for a four-person bench — plausible. Four hundred — something is off by an order of magnitude even though every check passed. Then run a second seed and confirm the story doesn't change qualitatively; if seed 1 and seed 2 disagree about which resource binds, the study needs more ensemble members before any of it means anything.

## Write DOCTOR.md

The command · the verdict line (`N pass · N warn · N fail`) · the table of checks with **exact evidence** (counts, troughs, totals — numbers, not adjectives) · a paragraph per warning explaining why it stands · anything you noticed by hand. Then hand to `rd-report`.
