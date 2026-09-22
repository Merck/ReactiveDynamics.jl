# `field_service` — the operator's pre-commitment sheet, and what the interview recovered

This is the control for the [`field_service`](field_service.json) capture. The operator persona wrote the sheet below **before the conversation started** and did not modify it afterwards; the assistant agent never saw it, and could not — it had no access to the directory the sheet was written in. Publishing it next to the transcript is what turns the session from a demo into a measurement: you can check, line by line, which facts the FRAME interview actually pulled out and which it left on the table.

The audit is deliberately unflattering. Three things were missed, one of them material, and they are named first.

## What the interview did not recover

| Fact on the sheet | What happened | Does it matter? |
|---|---|---|
| **Job duration is 4.3 h** | The interview elicited a *range* — "4-5 hours" — not the point value, and fitted 4.0 h inside it against the observed backlog. The operator's actual working figure surfaced only at turn 11, when they pushed back on the knife-edge claim. | **Yes, and this is the session's most instructive moment.** 4.3 h sits inside the `dt = 1.0` quantization dead zone (every duration in (4.0, 5.0] collapses onto 5.0 h), so the operator's own number was unrepresentable on the grid the assistant had authored. Asking "what is the average job length?" instead of "give me a range" would have exposed the grid defect in phase 2 instead of phase 4. |
| **A completed repair is worth \$2,800 (pharma) / \$900 (academic)** | Never modeled. The ledger carries only the two breach penalties (`breach_pharma` reward `-4200`, `breach_acad` reward `-650`); there is no positive value term anywhere in `model.rdj.json`. | **Not at the calibrated operating point; yes at the hot one.** In a settled system throughput equals arrival rate in every arm, so repair revenue is identical across arms and cancels exactly in the difference — omitting it costs nothing. At 4.75 h the baseline *diverges*: its rig-needing academic queue alone ends the horizon 216 jobs higher than the 6-laptop arm's, so **at least** \$195k of earned value is missing from that comparison, and the true figure is larger because the other queues diverge too. The omission biases the hot-case laptop benefit **downward**, which is the safe direction, but it should have been caught. |
| **Three quantities the operator flagged as gut estimates**: ~40% of academic jobs could tolerate an extra day; hiring would cut the backlog ~30%; the true cost of a pharma miss is ~\$8,400, about double the contractual penalty | None were carried into the study as sensitivities. | **Mixed.** The ~30% hiring estimate is the one the study went on to refute, so leaving it implicit cost nothing — refuting it *was* the finding. The \$8,400 figure is the real gap: the whole study is denominated in contractual penalties, and doubling the pharma weight is exactly the kind of re-weighting that could change a ranking. It was flagged by the operator as uncertain and never tested. The 40% academic-tolerance estimate would have changed the breach *definition*, not just a coefficient. |

The pattern is consistent: the interview was strong on **structure** and on the numbers that were **contractual** (headcount, unit counts, SLA windows, penalties, costs), and weak on the numbers the operator themselves had labelled soft. A FRAME skill that asked "which of these did you look up, and which did you estimate?" — and then pre-registered the estimated ones as sensitivities — would have closed all three gaps.

## What it recovered correctly

Every structural fact, and every figure the decision actually turned on.

| Fact on the sheet | In the study |
|---|---|
| 7 field technicians | `tech.init = 7` |
| 4 diagnostic laptops, required for most jobs, job cannot start without one | `laptop.init = 4`, `conserved` input arc; the 82.5% laptop-needing share was elicited in conversation, not on the sheet |
| Pharma tier: 24-hour SLA | breach hazard `queue/8`, 8 dispatch hours = 1 business day |
| Academic tier: 72-hour SLA | breach hazard `queue/24` |
| Pharma breach penalty −\$4,200 | `breach_pharma` reward `-4200.0` |
| Academic breach penalty −\$650 | `breach_acad` reward `-650.0` |
| Pharma 2.1 calls/day, academic 3.7 calls/day, Poisson | `call_pharma` rate 0.375/tick, `call_academic` 0.65/tick — 15 + 26 per week, i.e. arrivals on a 7-day week served on a 5-day dispatch week |
| 27% of jobs need a second visit | revisit probability 0.25, split 27.5% rig-needing / 72.5% parts-install |
| Strict pharma-first, FIFO within tier | ADR 0002 priority-weighted allocator, pharma 10.0 vs academic 1.0 |
| Hire 2 technicians at \$95k loaded | arm A, ≈\$190k/yr |
| 2 laptops at \$19k each | arm B, \$38k one-time — the operator corrected the assistant's first figure to this at turn 7 |
| Change the priority rule | arm C, priority 10.0/1.0 → 1.5/1.0 |
| Budget review in 6 weeks, 12-month planning horizon | `meta.tspan = 2080` = 8 h × 5 d × 52 w |
| Spare parts: \$18,500/week, \$340-\$520 per repair | Correctly scoped **out** — the operator judged parts "not really the blocker" and the brief records that as an explicit exclusion rather than an oversight |

Two figures the model was calibrated against — ≈14 open jobs today, 4 pharma and 28-30 academic breaches per month — are **not on the sheet**. The persona supplied them under questioning as elaborations consistent with it. They are the study's calibration anchor, so the reader should know they were generated during the interview rather than pre-committed: they cannot be used to score the interview, only to check the model reproduces what it was told.

## The sheet, verbatim

Reproduced exactly as written before turn 1, including its own uncertainty markers.

````markdown
# Field Service Operations Facts
## Pre-commitment sheet — do not modify after initial write

**Company**: Mid-sized industrial equipment company
**Equipment**: Laboratory analytical instruments (chromatography systems, mass spectrometers, spectrophotometers)
**My role**: Head of Field Service

## Customer Tiers

**Tier 1 (Premium)**: Pharmaceutical manufacturing companies
- Contractual SLA: 24-hour response time
- Higher service contract value

**Tier 2 (Standard)**: Academic institutions, research labs, quality-control labs
- Contractual SLA: 72-hour response time
- Standard service contract

## Current Resources

**Field technicians**: 7 technicians on staff

**Diagnostic laptops**: 4 units available
- Specialized equipment, expensive ($19,000 each)
- Required for most repair jobs
- Technician cannot start a job without one

**Spare parts budget**: $18,500 per week
- Average repair consumes $340 in parts (Standard)
- Average repair consumes $520 in parts (Premium - more expensive instruments)

## Job Characteristics

**Job duration**: Average 4.3 hours per job (includes drive time + on-site work)
- After completion, technician and equipment are immediately available

**First-visit success rate**: 73%
- 73% of jobs are resolved on first visit
- 27% require a second visit (for a part not on truck, additional diagnostics, etc.)

**Call arrival rate**:
- Tier 1 (Premium): 2.1 calls per day on average (Poisson distribution)
- Tier 2 (Standard): 3.7 calls per day on average (Poisson distribution)

## Financial Impact

**Value of completed repair**:
- Tier 1: $2,800 (maintains contract value + customer relationship)
- Tier 2: $900

**Cost of missed SLA**:
- Tier 1: -$4,200 (contractual penalty)
- Tier 2: -$650 (contractual penalty)

## Current Operating Rules

**Priority/contention rule**: When multiple jobs compete for available technician + diagnostic laptop:
1. All Tier 1 jobs have absolute priority over Tier 2
2. Within same tier, first-in-first-out (FIFO)

## The Decision

**Context**: Annual budget review coming up in 6 weeks. Need to justify any major spending requests.

**Current problem**: Backlog is consistently too long. Missing too many SLAs, especially on Tier 2, and occasionally on Tier 1 when we get a cluster of premium calls.

**Options under consideration**:
1. **Hire 2 additional field technicians** — $95,000/year salary each (loaded cost)
2. **Purchase 2 additional diagnostic laptops** — $38,000 total (2 × $19,000)
3. **Change the priority rule** — implement some kind of value-weighted or time-sensitive prioritization instead of strict Tier 1 first

**Planning horizon**: Next fiscal year (12 months)

## Uncertain Quantities (marked with ~ for gut estimate)

1. **~ Tier 2 SLA flexibility**: Exactly what percentage of Tier 2 calls could actually wait an extra day beyond the 72hr SLA without real damage to the relationship?
   - **Gut estimate: ~40%** (some academic labs are pretty flexible, but not all)

2. **~ Hiring impact**: If we hire 2 more technicians, how much would that actually reduce our backlog?
   - **Gut estimate: ~30% reduction** (but not sure if that's optimistic)

3. **~ True cost of Tier 1 SLA miss**: What's the real total cost of missing a Tier 1 SLA, including reputation damage and risk of contract non-renewal, beyond just the contractual penalty?
   - **Gut estimate: ~$8,400 total** (roughly double the penalty when you factor in relationship damage)

---

*Sheet finalized before conversation start. Extensions logged below if facts need elaboration.*
````

## What this sheet does and does not license you to conclude

It supports one narrow claim: **an operator's messy brief can be turned into a validate-clean, calibrated model without the modeller having privileged access to the answer.** The facts in the right-hand column of the recovery table were genuinely extracted through questioning, not read off a spec.

It does not show that the workflow handles a *resistant* informant. This persona answered cooperatively and in good faith, supplied elaborations on request, and volunteered its own uncertainty markers. A real operator forgets, rounds, contradicts themselves between meetings, and does not know which of their numbers is load-bearing. The three misses above are the closest this capture comes to that condition, and two of them happened because the persona had flagged something as soft and the interview did not press. That is the failure mode to design against, and one cooperative session is not enough evidence about it.
