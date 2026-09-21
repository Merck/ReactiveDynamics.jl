# BRIEF — What is the fifth scientist worth?

> The FRAME artifact (AGENT_WORKFLOW_CHARTER §3.1). Written in the operator's register — headcount, trials, queues — but every slot maps 1:1 onto the object model, so AUTHOR can compile it without a second interview. This is the reference example: a brief that is already a model in all but syntax.

## The decision

A discovery organization runs two therapeutic franchises — oncology and immunology — off one shared bench of four scientists. Both franchises have more candidate assets than the bench can staff. **Should we hire a fifth scientist?** The hiring case has to clear a number: the fully-loaded cost of the hire over the planning horizon. What it needs to be compared against is the *value the bench unlocks*, which nobody has computed.

**The answer is a single number with an error bar:** the expected-NPV difference between running the portfolio at four scientists and running it at five, over a 48-month horizon. Hire iff the loaded cost is below it.

## Why this is not a spreadsheet question

A spreadsheet can multiply a scientist's cost by four years. It cannot tell you what the *fifth* scientist is worth, because that value is not a property of the scientist — it is a property of the system. It depends on how often both franchises want the bench at once, which trial the allocator staffs when they collide, how long a trial holds its people, and which candidates therefore never get staffed at all. The quantity is a marginal response to relieving a constraint, so it exists only if the constraint binds and the contention has to be resolved.

Qualifying test (charter §2.4): **something binds** (the bench), **something takes time** (a trial holds its scientists for four months), **something is decided mid-flight** (which franchise gets the bench this tick). Three for three — this is a simulation question.

## Slots

| Slot | Value | Object-model target |
|---|---|---|
| **Scarce resource** | 4 scientists, reusable — a trial *occupies* them and returns them intact | plain place `scientist`, `init = 4`, held via a `conserved` input arc |
| **Metered resource** | trial budget, drawn per month while a trial runs, out of an ample pool | plain place `budget`, deep `init`, `perstep`/`consumed` input arc (a `@rate` draw), `cost = 1` |
| **Work in progress** | candidate assets, each carrying a franchise and a lifecycle phase | structured place `Project`, token kind `Project` with fields `program` / `phase` / `value` |
| **Starting portfolio** | 3 active oncology + 3 active immunology candidates | `population[]`: two entries, `count = 3` each |
| **Intake** | a steady stream of new candidates into each franchise | transitions `intake_onc` / `intake_imm`, rate 1.5/month, RHS `@structured` genesis |
| **The process that consumes capacity** | a trial: takes one active candidate of its franchise + 1 scientist + 20 budget/month, runs 4 months, succeeds 55% of the time | transitions `trial_onc` / `trial_imm`, `cycletime = 4`, `prob_of_success = 0.55` |
| **The contention rule** | oncology is protected — it gets the bench first when both want it | `priority` 3.0 vs 1.0, resolved by the priority-weighted allocator (ADR 0002) |
| **Outcome that carries value** | a launch, worth \$4M, realized only on trial success | plain place `launch`, `reward = 4000` (in \$k), produced on the RHS |
| **Horizon** | 48 months, monthly ticks | `meta.tspan = 48`, `meta.dt = 1.0` |
| **The lever** | `scientist.init`: 4 → 5 | one line of one JSON file |
| **Headline metric** | expected NPV in \$M, monthly ledger flows discounted at 8%/yr | `enpv` metric, `--scale=0.001` (pools are in \$k) |

## What could go wrong (pre-registered, for DIAGNOSE to check)

- **The bench doesn't actually bind.** If four scientists comfortably cover both franchises, the fifth is worth ≈ their accounting cost and there is no study. → `contention_exists`, `intended_resource_binds=scientist`.
- **Scientists get consumed instead of occupied.** The default arc modality is *consumed*; headcount declared that way would drain to zero and never come back, understating capacity for the whole run and overstating the value of a hire. → `modality_intent`, `expect_reusable=scientist`.
- **Budget binds instead of people.** Then the headline is a financing number, not a hiring number, and the memo would be wrong about what to do. → `intended_resource_binds`.
- **Reward credited at intake rather than at launch.** Putting the reward on the candidate stock instead of the `launch` product would pay for starting trials, not finishing them, and would make failure free.
- **Horizon too short to see launches.** A four-month trial inside a 48-month window is fine, but the portfolio queue is still growing at the end — the number is "value over 48 months", not a steady state. → `horizon_sufficient` (this one is expected to warn; it is a scope statement, not a defect).
- **Effect inside the noise.** A per-arm sample too small to separate a \$19M shift from ensemble spread. → the `effect_vs_noise` check on `rdx effect`.

## Acceptance

This study is the workflow's known-answer test. The published case study (`docs/literate/case_studies/marginal_scientist.jl`) reports a shadow price of **≈ +\$19M with a ±\$2.5M standard error** on 48-member ensembles from root seed 2026. The workflow reproduces it from the JSON documents in this directory; see `FINDINGS.md`.
