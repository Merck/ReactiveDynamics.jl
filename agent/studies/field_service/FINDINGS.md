# FINDINGS — Which lever fixes the field-service backlog?

> The REPORT artifact (AGENT_WORKFLOW_CHARTER §3.4). Every figure below traces to a command recorded in §7. A run is determined by `(model, seed)`: rerun the commands and you get these numbers, not numbers near them.

## The problem, in one paragraph

A lab-equipment vendor's field-service operation runs 7 technicians and 4 shared diagnostic laptops against 41 repair calls a week, under two contract tiers: premium pharmaceutical customers on a 24-hour fix SLA at \$4,200 a breach, and academic labs on 72 hours at \$650. Dispatch runs strict pharma-first. The backlog is chronic, roughly 14-15 open tickets at any time, and last month cost 4 pharma and 29 academic breaches. A budget meeting funds exactly one of three levers: **hire 2 more technicians** (≈\$190,000/yr), **buy 2 more diagnostic laptops** (\$38,000 one-time), or **change the dispatch priority rule** (free). Nobody knew which one the system was actually sensitive to, and the three cost figures differ by a factor of five, so guessing wrong is expensive in both directions.

## 1. Headline

## **+\$228,050 ± \$18,712 a year, for \$38,000 once — buy the laptops**

Two additional diagnostic laptops cut annual SLA penalties from **\$459,349** to **\$231,299** over a 12-month horizon: a saving of **\$228,050 with a standard error of \$18,712**, which is **12.2 standard errors from zero** on 48-member paired ensembles at root seed 2026. Every one of the five measured outcomes improves significantly. The purchase pays for itself in **2.0 months** and returns **\$6.00 per dollar spent** in year one alone, against a one-time cost.

The other two levers do not survive contact with the evidence:

| Lever | Cost | Annual penalty avoided | vs noise | Per \$1 spent | Verdict |
|---|---|---|---|---|---|
| **B — 2 diagnostic laptops** | **\$38,000 once** | **+\$228,050 ± \$18,712** | **12.2 σ** | **\$6.00** | **fund this** |
| A — 2 technicians | \$190,000/yr | +\$57,504 ± \$21,741 | 2.6 σ | \$0.30 | real but loses \$132,496/yr |
| C — drop strict pharma-first | \$0 | **+\$406 ± \$26,742** | **0.0 σ** | — | **cannot separate from doing nothing** |

**Arm C is a null result on all five metrics, and that is an answer, not a gap.** Penalties, pharma breaches, academic breaches, mean backlog and the year-end academic queue all land inside ±2 standard errors of zero at n = 48. `rdx effect` exited 1 with five `effect_vs_noise` findings. The correct statement is *"we cannot separate re-prioritizing from doing nothing"* — not "re-prioritizing saves \$406." The mechanism is plain once seen: with the laptop pool saturated, the total number of laptop-hours in a week is fixed, and a priority rule only decides the *order* in which work consumes them. Priority is a redistribution lever, not a capacity lever. Its theoretical best case is moving breaches from the expensive tier to the cheap one, and pharma breaches are already near their floor — the 51 a year that remain happen when *no* laptop is free at all, which no ordering rule can fix.

Two later results point the same way without relying on that mechanism. Refining the time grid measurably improved the priority tier — pharma breaches fall 12.88/year at 5.5 σ — while leaving the integrated queue unmeasurably changed at 0.14 σ, in a direction that is not even determined (§6). And adding two technicians at the stressed 4.75 h duration moved pharma breaches (−0.53/month, 3.3 σ) while leaving penalty, backlog, academic breaches and the year-end queue indistinguishable from zero (§4). **Three separate changes, each of which helped the premium tier or nothing at all; not one of them reduced the total backlog.** Only added laptop capacity ever did. Note what is *not* claimed here: §6 declines to identify the grid's mechanism, because the academic tier's noise floor is too wide to see a transfer that small. The pattern is three-for-three on the *outcome* — the total does not move — and agnostic about why.

**Arm A's \$57,504 is statistically real but economically wrong.** At 2.6 σ it clears the noise bar, and academic breaches and mean backlog both improve significantly — hiring is not useless. But 2 of its 5 metrics (pharma breaches, year-end academic queue) **cannot be separated from zero** at n = 48, and the saving recovers **\$0.30 of every dollar** of the \$190,000 annual salary cost. Funding it destroys about **\$132,500 a year** of value relative to doing nothing, and about \$170,500 a year relative to buying laptops instead.

## 2. The decision rule

**Buy the two diagnostic laptops iff they cost less than \$228,050 a year in avoided penalties. They cost \$38,000 once. Fund them.** The rule clears by a factor of six on the simulated baseline, and by a factor of 5.6 even after de-rating to the operator's own stated numbers (§6).

**Do not fund the hire on this evidence.** The rule for arm A is: hire iff the loaded cost of 2 technicians falls below \$57,504 a year. It is \$190,000. It fails by a factor of 3.3.

**Do not expect the priority change to move the needle, and do not spend political capital on it.** It is free, so there is no loss in adopting it for fairness reasons — academic customers currently starve, and that is a real service-quality problem. But it will not reduce the breach bill, and presenting it as a fix would be presenting noise as a result.

**Where laptop buying should stop — and this bounds the recommendation.** The rule is not "buy laptops without limit." Back-of-envelope arithmetic on the operator's own figures (not simulated — flagged as such): today the laptops run at 91.6% utilization and the technicians at 73.2%. At 6 laptops the laptops fall to **61.1%** while technicians stay at 73.2% — **the constraint moves off the laptops and onto the technicians.** Five laptops is the crossover point where the two are equally tight. So the 7th and 8th laptop should not be assumed to pay, and the hiring question genuinely reopens *after* the laptops land. That sequencing is the obvious follow-up study (§8) and it has not been simulated here.

**When to re-open the hiring question, grounded rather than guessed.** The system itself re-equilibrates in days, not quarters: at six laptops and 4.0 h the surplus is 93.4 laptop-hours/week, or 23.4 jobs/week of drain capacity, so a 10-job waiting queue clears in **0.43 weeks** of dispatch. The binding constraint on the re-run is therefore *measurement*, not settling. Distinguishing the true mean laptop-checkout time from the 4.366 h cliff (a 0.37 h gap) needs, at an assumed per-job spread of 1.0 h, about **100 jobs — two weeks at 51.25 jobs/week — for ±0.10 h**, and 400 jobs or about eight weeks for ±0.05 h; the assumed spread is the one input not supplied by the operator, and the checkout log will replace it. The new breach counts need a further month or two to accumulate enough events to be worth regressing (the model expects ~7 academic breaches/month post-purchase). **The grounded recommendation is to re-run 6-8 weeks after delivery**, with three measurements: the mean and spread of laptop checkout hours, whether jobs now wait on a technician rather than a laptop, and the new monthly breach counts by tier. An early trigger overrides the schedule: **if the mean waiting backlog has not fallen below ~5 jobs within three weeks of delivery, the model is wrong about something and should be re-run immediately.**

## 3. Why a spreadsheet cannot produce this

A spreadsheet can divide 41 calls a week by 7 technicians and report that the average technician is 73% busy. It cannot produce the number above, and the reason is specific to this system: all three qualifying conditions hold at once.

**Something genuinely runs out.** Both pools trough at zero during the run. Had the laptops been comfortably sufficient, the marginal laptop would be worth its purchase price and nothing more, and there would be no study.

**Work holds the resource for a while.** A repair occupies one technician *and* one laptop for four hours and returns both. That holding time is what converts a 92%-utilized resource into a 10-job queue, and it is why averages are useless here: the backlog is built by *coincidence*, not by the mean. Two things must be free at the same moment.

**Something has to decide who waits.** Pharma outranks academic 10-to-1, so when both tiers want the same free technician-and-laptop pair, the academic job eats the wait. That rule determines which tier's SLA clock runs out, and therefore which penalty gets charged.

Change any one of the three and the answer changes. The mechanism the arithmetic cannot reach is this: **a marginal laptop does not fix analyzers. It lets a repair *start* that would otherwise have queued** — and the value of that start depends on which tier's job was waiting, how long its clock had already run, whether it will bounce back as a revisit, and whether the *next* job to arrive would have found a free laptop anyway. That is a marginal response to relieving a constraint under contention, and it has no closed form. The proof that it needed simulating is arm A: the envelope alone would have said technicians were at 73% and left it there, whereas the ensemble puts an actual price — \$0.30 on the dollar — on the lever the organization was about to fund.

## 4. Text and visual result summary

Baseline, and each arm, on 48-member paired ensembles at root seed 2026. Bold means significant at 2·se; *italic* means the metric **cannot be separated from zero** and is therefore not a finding.

| Measure | Baseline | +2 laptops | +2 technicians | re-prioritize |
|---|---|---|---|---|
| Annual SLA penalty | \$459,349 | **\$231,299** | **\$401,845** | *\$458,943* |
| Pharma breaches / month | 4.28 | **3.50** | *4.21* | *4.21* |
| Academic breaches / month | 31.23 | **7.06** | **24.30** | *31.61* |
| Mean jobs waiting | 10.13 | **2.73** | **8.02** | *10.24* |
| Academic rig queue, year end | 11.98 | **1.67** | *8.94* | *9.29* |
| `rdx effect` exit code | — | **0 (clean)** | 1 (2 findings) | 1 (5 findings) |

The single most legible line is academic breaches: **31.2 a month today, 7.1 with two more laptops, 24.3 with two more technicians.** Two laptops remove 77% of them; two technicians remove 22% for five times the annual cost.

**What the queues say, and this is the evidence the binding claim actually rests on.** Baseline, seed 1:

| Queue | starts | peaks | ends |
|---|---|---|---|
| academic, **needs the rig** | 10 | **40** | 1 |
| academic, no rig | 2 | 6 | 0 |
| pharma, needs the rig | 2 | 5 | 0 |
| pharma, no rig | 0 | 5 | 0 |

The same 7 technicians serve all four queues. Only the two rig-requiring queues build, and the unprotected one builds eightfold. Work that needs a laptop waits; work that does not, does not. That asymmetry, plus the arm comparison above, is what identifies the laptops as the constraint.

### Duration robustness — the stress test at the operator's named worst case

n = 24 paired ensembles, root seed 2026, `dt = 0.25`, 4.75 h of laptop checkout: the duration the operator named as their worry ("4.7 or 4.8"), 8.8% above the 4.366 h laptop cliff, inside but not at the top of their stated 4-5 h range.

| At 4.75 h (n = 24) | Mean jobs waiting | Academic breaches/month | Pharma breaches/month | Year-end academic rig queue |
|---|---|---|---|---|
| 7 techs / 4 laptops | 114.48 | 377.85 | 4.33 | 218.38, **still climbing** |
| **7 techs / 6 laptops** | **8.10** | **25.17** | **2.74** | **2.67, settled** |
| 9 techs / 4 laptops | *108.04* | *359.09* | 3.80 | *206.38*, **still climbing** |

**The laptop arm is significant on all five metrics** — 28-30 σ on backlog, academic breaches and the year-end queue; `rdx effect` exit 0, zero findings. **The technician arm cannot be separated from zero on four of the five**: penalty Δ +\$173,052 ± \$155,526 (1.1 σ), academic breaches (0.9 σ), mean backlog (1.1 σ), year-end queue (1.0 σ); exit 1 with four `effect_vs_noise` findings. Only pharma breaches clear, at −0.53/month (3.3 σ) — the same signature seen elsewhere in this study, where a change helps the priority tier and leaves the total alone.

**No ratio is quoted between the two levers, because one of the two denominators cannot be separated from zero.** An earlier draft of this section carried "roughly 12:1 at the hot end," read off single runs; it is **retired**. Dividing one point estimate by another is not defensible when one of them is indistinguishable from nothing.

**The dollar figures at this duration are horizon-dependent, and they are not the finding.** The 4-laptop baseline never settles — its academic rig queue ends the 12-month horizon at 218 jobs and rising — so it has no well-defined annual penalty. Run 24 months instead and the \$3,165,258 baseline, the \$2,830,883 avoided, and the \$74.50-per-dollar that follows from them all grow, because the quantity they are measured against keeps growing. Quote the \$2.83M only with "over 12 months" attached, and prefer the horizon-free statement: **at 4.75 h, two laptops return the operation to a settled state (year-end queue 2.67, flat) while two technicians do not (206.38, still climbing).** That claim is robust to the horizon; the ratio is not.

One omission compounds in the same direction at this duration. The ledger prices SLA breaches and nothing else — the operator's \$2,800 and \$900 per completed repair were never modelled. In a settled system throughput equals arrival rate in every arm, so repair revenue cancels exactly in the difference and the omission is free. Here the baseline is *not* settled: its rig-needing academic queue alone ends 216 jobs above the 6-laptop arm's, so at least \$195k of earned value is missing from the comparison, and the true figure is larger because the other queues diverge too. The effect is to **understate** the laptop advantage at the hot duration, which is the safe direction, but it means the 4.75 h comparison is a lower bound on two independent counts.

<details>
<summary><b>Hero visual — the execution map (click to expand)</b></summary>

`figures/net.svg` — the network decorated with seed-1 run statistics:

```bash
julia --project=. agent/scripts/rdx.jl draw agent/studies/field_service/model.rdj.json \
  --out=agent/studies/field_service/figures/net.svg --exec-map --seed=1
```

What to look for: the `laptop` pool troughs at 0 and the two rig-requiring repair paths are the ones whose input queues carry stock. What *didn't* starve is equally informative — the two no-rig queues clear to 0 and stay there, despite being served by the same technician pool. That contrast is how you know the shortage is the rig and not the headcount.

</details>

## 5. The model

<details>
<summary><b>Formalism — places, transitions, modalities, units (click to expand)</b></summary>

**The model is data.** Four JSON documents, no authored Julia anywhere in this study — the model uses plain places only, no structured tokens, so there is no `host.jl` and no `--host` flag on any command.

| Operator's language | Formal role |
|---|---|
| 7 technicians | plain place `tech`, `init = 7`, held by a `conserved` input arc on all four repair transitions |
| 4 diagnostic laptops | plain place `laptop`, `init = 4`, `conserved` on both rig-requiring repair transitions |
| four job queues | plain places `q_pharma_rig`, `q_pharma_quick`, `q_acad_rig`, `q_acad_quick`; `init` 2 / 0 / 10 / 2 = 14 open tickets at t=0 |
| calls arriving | transitions `call_*`, Poisson, rates 0.309375 / 0.065625 / 0.53625 / 0.11375 per tick (= 15 + 26 per week, split 82.5% / 17.5% by whether the job needs the rig) |
| a repair | transitions `repair_*`, `cycletime = 4` ticks, consuming 1 queued job + 1 `tech` (+ 1 `laptop` on the rig paths) |
| revisits | `choose`-valued RHS arc multiplicities: 0.75 → fixed, 0.06875 → back to the rig queue, 0.18125 → back to the no-rig queue. The bounced job keeps accruing SLA exposure |
| strict pharma-first | `priority` 10.0 on pharma repairs vs 1.0 on academic, resolved by the priority-weighted allocator (ADR 0002) |
| SLA breaches | transitions `breach_*_tick`, rate = (tier queue) × hazard, accumulating into places `breach_pharma` (`reward = -4200`) and `breach_acad` (`reward = -650`), so the `reward` metric *is* penalty dollars |
| backlog measurement | transition `backlog_tick`, rate = sum of all four queues, accumulating `backlog_meter`; mean waiting backlog = `final:backlog_meter` ÷ 2082 ticks |

**Units, stated explicitly because this is the error that survives every check.** Time is in **dispatch hours**: 1 tick = 1 hour, 8 hours/day × 5 days × 52 weeks = **2,080 ticks = 12 months**. No evening shifts. Money is in **whole dollars** — \$4,200 and \$650 per breach, not \$k and not \$M. The `reward` metric is negative because a breach is a penalty; this report negates it and calls it "penalty."

**The arms are one-line diffs, and the diff is the lever.**

```diff
# arm_buy_2_laptops.rdj.json
      "name": "laptop",
-     "init": 4.0
+     "init": 6.0
```
```diff
# arm_hire_2_techs.rdj.json
      "name": "tech",
-     "init": 7.0
+     "init": 9.0
```
```diff
# arm_reprioritize.rdj.json  (both pharma repair transitions)
-     "priority": 10.0
+     "priority": 1.5
```

</details>

## 6. Caveats — scope, not apology

**The baseline reproduces today's operation, which is why the arms are believable.** The model was built from operational facts and then asked. It reports **10.13 jobs waiting** plus Little's-law work in progress of 4.08 jobs = **14.2 open tickets** against the operator's observed 14-15, **4.28 pharma breaches/month** against 4, and **31.2 academic/month** against 29. The model's annual penalty of \$459,349 runs **+7.4% hot** against the operator's own \$427,800.

**Consequence of running hot: the headline is mildly overstated, and it still clears easily.** De-rating the \$228,050 proportionally to the operator's own baseline gives **\$212,387**, a **5.6:1** return on \$38,000 instead of 6.0:1. Both versions support the same decision, which is the test that matters.

**Two fitted dials, disclosed.** (1) Job duration, fitted inside the operator-stated 4-5 hour range; (2) the breach-hazard scale, one per tier, fitted to the observed 4 and 29 breaches/month. Two hazards against two observations is exactly identified. Everything else — call volumes, tier split, laptop attach rate, revisit rate and its rig split, headcount, laptop count, penalties, SLA windows, the dispatch day — is an operator-supplied fact.

**The breach counter is proportional to queue-time, not a hard deadline.** It accrues breaches at a calibrated hazard per waiting job per tick rather than firing exactly when a 24- or 72-hour clock expires. So it is a sound instrument for *comparing* arms — which is what a lever ranking needs — while absolute breach levels are **calibrated, not predicted**. Lead with the ranking.

**`intended_resource_binds` does not discriminate between the two candidate resources.** The check passes with `--expect-binds=laptop`, but re-running it with the counter-hypothesis `--expect-binds=tech` **also passes**, because the check asks only whether the named pool reaches zero and under this load both do. It therefore confirms the brief named a resource that genuinely runs out, and would have caught a model where laptops were never scarce — but it is **not** evidence that laptops rather than technicians are the constraint. That claim rests entirely on the queue asymmetry in §4 and on the arm comparison in §1. Stated explicitly because a check cited as proof of something it cannot establish is exactly the decoration the DIAGNOSE phase exists to prevent.

**`horizon_sufficient` warns on all four arms — the figures are "over 12 months," not steady state.** 3-5 repairs are mid-job at tick 2,080 and the cumulative counters are still climbing. Out of ~2,100 repairs a year this is a rounding error rather than truncated value, but every number here is a 12-month figure, and the direction of the truncation makes the headline a **lower bound**.

**The engine's time grid was coarser than the question, and successive ensembles revised my reading of what that cost — three times.** `cycletime` is expressed in time units while progress accrues `dt` per tick (`firing.state += qs[i]*state.dt`, completion at `>= transCycleTime`, `src/solvers.jl:483,684`), so a firing occupies `ceil(cycletime/dt)` ticks; at the authored `dt = 1.0` every duration in (4.0, 5.0] collapses onto 5.0 h, which is why 4.1, 4.2, 4.3 and 4.5 h produced byte-identical runs — and it means **the operator's own 4.3 h working figure fell on the divergent side of the grid while the calibrated baseline sat at 4.0 h.** Re-authored at `dt = 0.25` (a one-line `meta.dt` change; rates are per unit time and drawn as `rand(Poisson(dt*rate))`, `src/solvers.jl:316`, so neither rates nor `tspan` rescale), what the refinement measurably moved is **one tier, and nothing else anywhere in the system**:

| Baseline, **matched n = 24** both grids, root seed 2026 | `dt = 1.0`, n=24 | `dt = 0.25`, n=24 | difference | σ | verdict |
|---|---|---|---|---|---|
| **pharma breaches / year** | 50.50 ± 1.61 | 37.62 ± 1.69 | **−12.88 ± 2.33** | **5.52** | **real** |
| academic breaches / year | 392.58 ± 32.80 | 396.00 ± 42.37 | +3.42 ± 53.58 | **0.06** | cannot be separated from zero |
| annual penalty | \$467,279 ± 24,502 | \$415,425 ± 30,965 | −\$51,854 ± 39,487 | 1.31 | cannot be separated from zero |
| backlog meter (integrated queue) | 21,996 ± 1,699 | 21,603 ± 2,175 | −393 ± 2,760 | **0.14** | cannot be separated from zero |
| year-end academic rig queue | 12.50 ± 2.08 | 7.71 ± 1.93 | −4.79 ± 2.84 | 1.69 | cannot be separated from zero |

**Read this table only against itself.** Both columns are `rdx ensemble`, n = 24, root seed 2026, both exit 0. The **published** study baseline is n = 48 (\$459,349, 51.38 pharma, 374.73 academic breaches/yr) and the n = 24 `dt = 1.0` subsample reads differently — \$467,279 and 392.58 academic. **Splicing the n = 48 headline figures into this comparison is what produced an earlier spurious "+21.27 academic breaches" term**, and that retracted artefact is the reason every row here carries its n.

**The integrated queue is unchanged at 0.14 σ, and its difference does not even have a determined sign** — it reads −393 here where the mismatched-n version read +507. That kills my first explanation. A single-seed reading of 5.61 waiting jobs led me to infer that the coarse tick charged a full tick of re-admission delay against every job, a capacity effect; **that explanation is withdrawn**, because a capacity effect could not leave the integrated queue unmoved to within a seventh of a standard error, in an undetermined direction.

**What survives is a real pharma improvement with no identified counterparty — and the honest position is that the mechanism is undetermined, not that it is free.** The pharma tier gains 12.88 breaches a year at 5.52 σ while the academic tier is flat to **0.06 σ**, the penalty to 1.31 σ and the queue to 0.14 σ. It is tempting to call this a redistribution — a finer tick gives the priority-weighted allocator four allocation opportunities an hour instead of one, so the 10:1 pharma weighting acts sooner — but **a redistribution needs a recipient, and none is visible.** Nor can one be, at this sample size: the academic tier's own noise floor is ±42.37 breaches/year, **3.3 times the size of the entire pharma gain**, and resolving a transfer that small against it at 2·se would need the difference's se down to 6.44/yr from 53.58 — roughly **1,700 seeds per arm** rather than 24. **The data therefore neither supports nor excludes the redistribution reading; the question is below the resolution of this design.** What is observed, and all that is claimed: the finer grid bought a real improvement in the tight-SLA tier and nothing else in the system measurably changed. The efficient test is not more seeds but the equal-priority arm run at both grids — if the pharma gap closes when the 10:1 weighting is flattened, the mechanism is the priority interaction, and that test needs no power in the academic tier at all.

**One methodological retraction to keep, because it applies to any future version of this table.** An earlier draft cited the fact that \$4,200·Δpharma + \$650·Δacademic reproduces Δpenalty as support for the redistribution reading. **That is an identity, not evidence.** This model's ledger defines `reward` as exactly −4200·pharma − 650·academic, so the sum closes for *any* two runs, signal or noise. It confirms the metrics were parsed correctly and establishes nothing about mechanism.

**Refining the grid did not reduce the laptops' measured value.** Re-run at `dt = 0.25` on the same 4.0 h calibrated duration, the same two laptops come out at **+\$288,644 ± \$31,502** (9.2 σ, all five metrics significant, exit 0) against the published **+\$228,050 ± \$18,712** (12.2 σ). The two estimates differ by \$60,594 ± \$36,640 — **1.65 σ, not separable** — so the finer grid does **not** establish a larger effect, and the earlier "±10% of tick-resolution uncertainty" caveat is retired in favour of the one-sided statement the data does support: **there is no evidence the published headline is an overstatement, and the lever clears by ≥9 σ at both resolutions.** Per dollar spent, \$6.00 at `dt = 1.0` and \$7.60 at `dt = 0.25`.

**Neither grid is demonstrably the better-calibrated one, and the spread across resolutions and sample sizes is about 10%.** The operator's stated \$427,800 falls *between* the two grids: `dt = 1.0` reads **+7.4% hot** at \$459,349 (n=48) and **+9.2% hot** at \$467,279 (n=24), while `dt = 0.25` reads **−2.9% cold** at \$415,425 (n=24). The matched-n penalties differ by \$51,854 ± \$39,487 — **1.31 σ, not separable** — so the apparent improvement in fit is not claimed. What is claimable: **every baseline reading lands within 10% of the operator's own figure, and they straddle it.** Scaling the finer-grid effect estimate to that figure gives \$297,242 (7.82:1); the fact that this de-rating moves the estimate *up* where the `dt = 1.0` de-rating moved it *down* reflects only the straddle, which is itself a difference inside the noise. The academic breach hazard was fitted at `dt = 1.0` and would need refitting before the finer-grid model is quoted for absolute breach levels — which is why the headline is not re-based on it.

**The cliff is arithmetic, it has an exact location, and hiring does not move it.** Capacity is division, not simulation: 36.64 laptop-needing jobs/week and 51.25 total jobs/week against 40 dispatch-hours per week per unit. The duration at which each resource hits 100% is therefore fixed, and the binding cliff is the lower of the two.

| Configuration | Laptop cliff | Technician cliff | Binding cliff |
|---|---|---|---|
| Today — 7 techs / 4 laptops | **4.366 h** | 5.463 h | **4.366 h — laptops** |
| +2 technicians — 9 / 4 | **4.366 h** | 7.024 h | **4.366 h — laptops, moved by zero** |
| **+2 laptops — 7 / 6** | 6.550 h | **5.463 h** | **5.463 h — technicians** |
| Both — 9 / 6 | **6.550 h** | 7.024 h | **6.550 h — laptops again** |

At 4.30 h the four-laptop operation runs at **98.5%** of laptop capacity, with **0.57 laptop-jobs/week** of headroom; at 4.75 h it is at **108.8%**, i.e. over capacity, where no allocation rule and no headcount can prevent divergence. **Six laptops stay sub-critical across the operator's entire stated 4.0-5.0 h range — 61.1% of capacity at 4.0 h rising to 76.3% at 5.0 h — while four laptops cross their own cliff at 4.366 h and are over capacity for the rest of that range.** And **two more technicians move the binding cliff by 0.000 h**, because headcount cannot relieve the resource that is actually short — the clearest single refutation of the lever the organization was about to fund.

**Three figures from the first draft of this study are retired, named here so they are not quoted from an earlier copy:** that the system is "stable at 4 hours of resource occupancy and in runaway at 5"; that the operator's stated 4-5 h range "admits exactly two distinct models"; and the 5-tick pair **mean backlog 142 jobs / 1,023 academic breaches per month**. All three describe the authored `dt = 1.0` grid rather than the operation. The surviving observation is narrower and is a discretization fact: at `dt = 1.0`, durations of 4.1, 4.2, 4.3 and 4.5 h produce byte-identical runs. `DOCTOR.md` §"Cycletime quantization" carries the same disclosure.

**Arm C's null is specific to the softening tested.** Pharma priority went 10.0 → 1.5, not to full parity or an inversion. A stronger rule was not tested. The mechanism bounds it, though: with laptop-hours fixed, no ordering rule creates capacity, so the best any priority change can do is move breaches between tiers.

**Uncertainty in the operator's own inputs is not in the error bars.** The ±\$18,712 reflects call-arrival timing and revisit incidence — the right uncertainties for a queueing question — not uncertainty in the 41 calls/week or the 82.5% attach rate. The duration sensitivity above is the substitute for a full parameter sweep.

**Not simulated:** the combined arm (laptops *and* technicians), any laptop count above 6, and a stronger priority rule. The \$0.30-on-the-dollar verdict for hiring holds *at today's 4 laptops*; it should be re-tested after the laptops land, because that is precisely when technicians become the tighter resource.

**One engine defect was found and worked around in this study's own files, with no change to `src/`.** `rdx validate` returned exit 0 and zero diagnostics on a document whose `choose` alternatives held bare numbers, and `rdx run` on the same file then died with `MethodError: no method matching node_from_dict(::Float64)` (`src/serialize.jl:84`) — the bare-number shorthand is not honored inside `choose` alts. Every alt value in these documents is written as an explicit `const` node. Separately, repeated CLI flags (`--expect-reusable`, `--tally`) silently collapse to the last occurrence; comma lists work. Both are recorded in `DOCTOR.md`.

## 7. How we got it, and provenance

1. **FRAME** → [`BRIEF.md`](BRIEF.md) — the decision, the three levers, 17 slots, 8 pre-registered failure modes with their DIAGNOSE flags, the calibration anchor.
2. **AUTHOR** → [`model.rdj.json`](model.rdj.json), [`arm_buy_2_laptops.rdj.json`](arm_buy_2_laptops.rdj.json), [`arm_hire_2_techs.rdj.json`](arm_hire_2_techs.rdj.json), [`arm_reprioritize.rdj.json`](arm_reprioritize.rdj.json) — all four `rdx validate` clean, zero diagnostics.
3. **DIAGNOSE** → [`DOCTOR.md`](DOCTOR.md) — all four arms **8 pass · 1 warn · 0 fail**, exit 0.
4. **REPORT** → this document.

The command that produced the headline, verbatim and runnable from the repository root:

```bash
julia --project=. agent/scripts/rdx.jl effect \
  agent/studies/field_service/model.rdj.json \
  agent/studies/field_service/arm_buy_2_laptops.rdj.json \
  --seeds=48 --root-seed=2026 \
  --metric=reward,final:breach_pharma,final:breach_acad,final:backlog_meter,final:q_acad_rig
```

It prints `delta: 228050.0, se: 18711.52` on `reward` and exits **0**. The other two arms, same invocation with the arm file swapped, exit **1** and carry their `effect_vs_noise` findings in `findings[]`. Member *k* of each arm is seeded from `hash((root_seed, k))`, so both arms see the same seed stream — the paired comparison is what makes a \$57k effect detectable at all. **A run is determined by `(model, seed)`**: rerun the above and you get 228050.0, not a number near it.

**The second round — the tick-grid rebuild — added eight more documents, all `rdx validate` clean.** They are listed here because §4 and §6 quote them and a reader should be able to find the input behind each number.

| Document | `dt` | Duration | Config | Quoted in |
|---|---|---|---|---|
| [`dt25_baseline`](dt25_baseline.rdj.json) | 0.25 | 4.00 h | 7 techs / 4 laptops | §6 grid table, §6 headline re-run |
| [`dt25_laptops`](dt25_laptops.rdj.json) | 0.25 | 4.00 h | 7 / **6** | §6 headline re-run (+\$288,644) |
| [`dt25_hot_baseline`](dt25_hot_baseline.rdj.json) | 0.25 | **4.75 h** | 7 / 4 | §4 stress test |
| [`dt25_hot_laptops`](dt25_hot_laptops.rdj.json) | 0.25 | 4.75 h | 7 / **6** | §4 stress test |
| [`dt25_hot_techs`](dt25_hot_techs.rdj.json) | 0.25 | 4.75 h | **9** / 4 | §4 stress test (the four-way null) |
| [`probe_cal_4.10h`](probe_cal_4.10h.rdj.json) · [`probe_cal_4.30h`](probe_cal_4.30h.rdj.json) | 1.0 | 4.10 / 4.30 h | 7 / 4 | §6 — the pair that proved the (4.0, 5.0] collapse by running byte-identically |
| [`probe_anchor_4.125h_dt0.125`](probe_anchor_4.125h_dt0.125.rdj.json) | **0.125** | 4.125 h | 7 / 4 | `DOCTOR.md` — the single-seed calibration probe, carried as an indication only |

Reproducing §4 and §6 is the same `effect` and `ensemble` invocations with these files substituted, at `--seeds=24 --root-seed=2026`. The raw JSON envelopes are **not** committed: they are outputs, fully determined by `(model, seed)`, and the commands that regenerate them are in the tool blocks of [`../../sessions/field_service.json`](../../sessions/field_service.json) and in `DOCTOR.md`.

Run cost, for anyone reproducing this: each 48-seed two-arm `effect` call takes **16-20 minutes at ~96% CPU** when run **sequentially**. Do not run the three concurrently against one Julia depot — three parallel invocations spent 116 seconds of CPU across 2 h 29 m at 1% CPU, blocked rather than computing, and produced nothing. See `DOCTOR.md` §"Run cost."

## 8. What to do next, if you want it

The infrastructure is built, so each of these is one file and one command:

- **The sequencing question.** Simulate 6 laptops + 9 technicians. The envelope says technicians become the tighter resource once the laptops land, which would make the hire that fails today a reasonable second purchase. This is the highest-value follow-up.
- **Where laptop buying stops.** Arms at 5, 7 and 8 laptops to trace the diminishing-returns curve and find the last laptop that pays.
- **The knife edge, priced.** Sweep call volume at 4 and at 6 laptops to find how many extra calls a week each configuration absorbs before the backlog diverges. This converts "we're running close to the edge" into a headroom number, which is a better capacity-planning instrument than the breach count.
- **A stronger priority rule.** Full parity, or protecting a fixed academic share, to confirm the null is about priority as a class rather than about the particular weights tested.
- **~~Error bars on the duration stress test~~ — done, and it changed three claims.** The 24-seed ensembles at 4.75 h landed (§4): the laptop arm clears at 28-30 σ, the technician arm cannot be separated from zero on four of five metrics, and the single-seed "12:1" ratio was retired rather than refined. Left here as a record that the follow-up was run, not just proposed.
- **Re-calibrate at a resolved tick, and refit the academic hazard while you are there.** The observation-matching duration (~4.10 h) is not representable at `dt = 0.25`; `dt = 0.125` reaches it, but the academic breach hazard was fitted at `dt = 1.0` and must be refitted before the finer-grid model is quoted for absolute breach levels. This is now the one open thread that could move a number rather than add one: the cross-grid comparison in §6 leaves a real −12.88/year pharma improvement with no identified counterparty, and a refit is how you find out whether that is a calibration artefact or a genuine allocator effect.
