# BRIEF — Which lever fixes the field-service backlog: technicians, diagnostic laptops, or the priority rule?

> The FRAME artifact (AGENT_WORKFLOW_CHARTER §3.1). Written in the operator's register — technicians, jobs, laptops, SLAs, backlog — but every slot maps 1:1 onto the object model so AUTHOR can compile it without a second interview.

## The decision

A lab-equipment vendor runs field service for installed analyzers with **7 technicians** and **4 shared diagnostic laptops**. Two contract tiers compete for that capacity: **premium pharmaceutical** customers (24-hour fix SLA, \$4,200 penalty per breach) and **academic labs** (72-hour SLA, \$650 per breach). Dispatch currently runs **strict pharma-first**, which keeps the expensive breaches low and lets the academic queue starve. The backlog is the presenting complaint. A budget meeting in six weeks has to fund exactly one of three levers.

**The three candidate levers, funded against the same budget line:**

| Lever | Change | Assumed cost (operator to confirm) |
|---|---|---|
| **A — hire** | 7 → 9 technicians | ≈ \$190k/yr (2 × \$95k fully loaded) |
| **B — equip** | 4 → 6 diagnostic laptops | **\$38k one-time** (operator-confirmed) |
| **C — re-prioritize** | strict pharma-first → academic gets a protected share | ≈ \$0, a dispatch policy change |

**The answer is a ranking with error bars:** annual SLA-penalty dollars avoided per lever, each versus the status quo, over a 12-month horizon. The status quo costs **\$427,800/yr in breach penalties** (4 pharma + 29 academic per month, at the operator's stated penalties). Fund the lever whose avoided-penalty confidence interval clears its cost by the widest margin; if two arms overlap inside the noise, say so rather than ranking noise.

## Why this is not a spreadsheet question

A spreadsheet can divide 41 calls a week by 7 technicians and tell you the average technician is 82% busy. It cannot tell you why 15 jobs are sitting open right now, because the backlog is not made by the average — it is made by *coincidence*. Two resources (a technician AND a laptop) must be free at the same moment; when they are not, the job waits; while it waits its SLA clock runs; and a strict priority rule decides which tier eats the wait. The quantity the operator needs — the marginal breach reduction from one more unit of *which* resource — is a response to relieving a constraint under contention. It has no closed form, and critically, **it is not monotone in the obvious resource**: relieving the non-binding resource yields approximately nothing, which is precisely the \$190k error this study exists to prevent.

Qualifying test (charter §2.4): **something binds** (laptops at ≈95% utilization, technicians at ≈82%), **something takes time** (a job holds one technician and one laptop for 4.5 hours door-to-door), **something is decided mid-flight** (strict pharma-first when both tiers want the same free technician-and-laptop pair). Three for three — this is a simulation question.

## The envelope that reframed the study

Run from the operator's own figures before any model existed (`julia -e`, recorded in the transcript): 41 calls/week (15 pharma + 26 academic), ×1.25 for the one-in-four revisit rate, a 4.5-hour job, an 8-hour dispatch day.

Of the 10.25 revisits/week, 27.5% need the diagnostic rig again (the first visit did not nail the diagnosis) and 72.5% are parts installs that do not — so laptop-needing demand is 41×0.825 + 10.25×0.275 = **36.64 laptop-jobs/week** against a 4-laptop capacity of 35.6.

| job duration (operator's range is 4-5h) | technician utilization | laptop utilization |
|---|---|---|
| 4.0 h | 73.2% | 91.6% |
| **4.5 h (midpoint)** | **82.4%** | **103.1%** |
| 5.0 h | 91.5% | 114.5% |

**The laptops are the saturated resource, not the technicians — at every point in the operator's own stated duration range.** That robustness is the important part: the *level* of laptop load is uncertain by ±12 points because the job duration was given as a range, but the *ordering* (laptop load exceeds technician load by ~19 points throughout) is invariant. This is a hypothesis the model must confirm or refute, not a conclusion — but it is why the study is framed as a three-way comparison rather than a hire/don't-hire question. Note the two sub-arithmetics the arms will test: 9 technicians leaves laptop utilization *unchanged*, while 6 laptops drops it by a third.

At the 4.5h midpoint the envelope puts laptop demand slightly **past** capacity, which cannot be literally true of a surviving business — so job duration becomes the study's first calibration dial, fitted inside the operator-stated [4.0, 5.0] range against the observed backlog. See the calibration anchor below.

## Slots

| Slot | Value | Object-model target |
|---|---|---|
| **Scarce resource 1** | 7 technicians, reusable — a job *occupies* a tech and returns them intact, same day | plain place `tech`, `init = 7`, held via a `conserved` input arc |
| **Scarce resource 2** | 4 diagnostic laptops, reusable, held for the whole job | plain place `laptop`, `init = 4`, `conserved` input arc |
| **Reusable or spent?** | Both reusable. Neither is consumed by a repair. Nothing here is a consumable — parts are explicitly out of scope per the operator ("not really the blocker") | `conserved` on every `tech` and `laptop` input arc, never the default |
| **Metered resource** | None. No continuously-drawn budget or power in scope | no `perstep` arcs |
| **Work in progress** | Broken analyzers awaiting repair, in two interchangeable-within-tier pools. Jobs carry no per-item attributes that matter to the decision | plain places `queue_pharma`, `queue_academic` — **no structured tokens, so no `host.jl`** |
| **Starting state** | ≈14 open jobs today, predominantly academic (the starved tier) | `init`: `queue_pharma = 2`, `queue_academic = 12` |
| **Intake** | 15 pharma + 26 academic calls/week, Poisson | transitions `call_pharma` (`rate = 0.375`/tick), `call_academic` (`rate = 0.65`/tick) |
| **The process** | One job takes 1 waiting job + 1 technician + 1 laptop, runs 4.5 hours, and fixes the analyzer. 82.5% of jobs need a laptop; the other 17.5% (fuse, pump-module swap) need only a tech | transitions `repair_pharma`/`repair_academic` (`cycletime = 4.5`) for laptop jobs, plus `quick_pharma`/`quick_academic` with no `laptop` arc, split 0.825 / 0.175 by rate |
| **Revisits** | 1 in 4 jobs needs a second trip. 72.5% of those are parts installs that do not re-occupy a laptop; 27.5% need more diagnostics and do | RHS of a repair re-injects the job into a queue with probability 0.25, split 0.06875 back to the rig queue and 0.18125 to the no-rig queue, via `choose`-valued arc multiplicities. The re-injected job keeps accruing SLA exposure, which is the point |
| **Contention rule** | Today: strict pharma-first. Both tiers want the same free tech-and-laptop pair and pharma always wins | `priority` — pharma 10.0 vs academic 1.0 in the baseline (the ADR 0002 priority-weighted allocator) |
| **Outcome that carries value** | An SLA breach. Pharma \$4,200, academic \$650. Accrues when a job is still waiting past its clock | places `breach_pharma` (`reward = -4200`), `breach_academic` (`reward = -650`), produced by the breach transitions; the `reward`/`net` metric is therefore penalty dollars |
| **The SLA clock** | Pharma 24 clock-hours = 1 business day = **8 ticks**; academic 72 clock-hours = 3 business days = **24 ticks** | breach transitions with per-waiting-job hazard `queue/8` and `queue/24` |
| **Horizon and tick** | **1 tick = 1 dispatch hour.** 8 dispatch hours/day × 5 days × 52 weeks = **2080 ticks = 12 months**. No evening shifts | `meta.tspan = 2080`, `meta.dt = 1.0` |
| **Units** | Time in dispatch hours. Money in **whole dollars** (not \$k, not \$M) — penalties are \$4,200 and \$650 as stated | restated in FINDINGS to prevent a 1000× memo error |
| **The lever** | Arm A `tech.init` 7 → 9 · Arm B `laptop.init` 4 → 6 · Arm C `priority` 10.0/1.0 → 1.5/1.0 | one line of one JSON file per arm |
| **Headline metric** | Annual SLA-penalty dollars, and breaches/month per tier | `reward` (penalty dollars), `completions[:breach_pharma]`, `completions[:breach_academic]`, `sum:queue_*` for mean backlog |

## The calibration anchor

The operator supplied today's observed state, and the baseline has to reproduce it before any arm is quoted: **≈14-15 open jobs**, **4 pharma breaches/month**, **28-30 academic breaches/month**. This is the study's most valuable constraint — a model that cannot reproduce today has no standing to forecast next year.

The model has exactly **two fitted dials against those three observations**, and both are disclosed so the report cannot present a calibrated model as a predicted one:

1. **Job duration**, fitted inside the operator-stated [4.0, 5.0] hour range to reproduce the ≈14-job backlog. Not a free parameter so much as a range the operator gave rather than a point.
2. **The breach-hazard scale**, one per tier, fitted to reproduce 4 pharma and 29 academic breaches/month. Two hazards against two observations — exactly identified, no slack to overfit with.

Everything else — call volumes, tier split, laptop attach rate, revisit rate and its rig split, headcount, laptop count, penalties, SLA windows, the dispatch day — is an operator-supplied fact carried straight into the document.

## What could go wrong (pre-registered, for DIAGNOSE to check)

- **The laptops don't actually bind.** If 4 laptops comfortably cover demand, the envelope's 95% is wrong and the whole reframing collapses. → `--expect-binds=laptop`.
- **Technicians bind instead, and the hiring case was right all along.** The honest outcome if the model says so. Both resources get a binding check so the answer is not assumed. → `--expect-binds=tech` on the counter-hypothesis.
- **Technicians or laptops declared consumed rather than occupied.** The default arc modality is *consumed*; headcount declared that way drains to zero and never returns, understating capacity all run and inflating the value of any capacity lever. → `--expect-reusable=tech`, `--expect-reusable=laptop`.
- **The breach counters look like runaway growth.** They are legitimate cumulative counters over 2080 ticks and will grow monotonically. → `--tally=breach_pharma`, `--tally=breach_academic`.
- **The breach clock is a hazard, not a hard deadline.** A memoryless per-job hazard lets a job breach earlier than its SLA and survive later than it; only the *calibrated rate* is faithful. Mitigation: calibrate to the observed 4/29 per month, and lean on *differences between arms* rather than absolute breach levels. Recorded as a scope caveat, not a defect.
- **The job duration is a range, not a number.** Laptop load runs 91.6% → 114.5% across the operator's stated 4-5 hours, so absolute breach levels are soft. The ranking is not: laptop load exceeds technician load by ~19 points at every duration in the range. The report must lead with the ranking and treat absolute levels as calibrated, not predicted.
- **The revisit rig-split is load-bearing in arm B's favour.** If more revisits need the rig than the stated 27.5%, arm B gets *stronger*, not weaker — so the assumption is conservative with respect to the study's likely recommendation.
- **Horizon too short.** 2080 ticks against a 4.5-tick job is ample; the queue may still be off steady state at t=0 because of the 14 seeded jobs. → expect `horizon_sufficient` to pass; discard no warm-up without saying so.
- **The effect is inside the noise.** Especially plausible for arm C, a free policy change, and for pharma breaches, which are rare (4/month → ~48/yr) and therefore noisy. → the `effect_vs_noise` check on `rdx effect`, with enough seeds to separate the arms or an explicit "cannot separate" verdict.

## What "done" looks like

Three numbers the operator can put on one slide and defend: annual penalty dollars avoided by arm A, by arm B, and by arm C, each with a standard error, each divided by its cost. Plus the diagnostic sentence that makes the ranking legible — *which resource is actually the wall*. If the ranking is "the \$38k of laptops beats the \$190k of headcount," that claim must survive the binding check, the calibration anchor, and the noise check before it reaches FINDINGS.
