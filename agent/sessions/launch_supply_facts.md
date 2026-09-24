# `launch_supply` — the operator's pre-commitment sheet, and what the interview recovered

This is the control for the [`launch_supply`](launch_supply.json) capture. The operator persona wrote the sheet below **before the conversation started** and did not modify it afterwards. The assistant agent never saw it, and could not: it was forbidden from reading the directory the sheet was written in, and it held to that even when the coordinator mistakenly pointed it at a compute file in the relay directory (turn 6, logged in `agent/studies/launch_supply/capture_log.md`). Publishing the sheet next to the transcript turns the session from a demo into a measurement, because you can check line by line which facts the FRAME interview pulled out and which it left on the table.

The audit leads with what was missed.

## What the interview did not recover, or did not test

| Fact on the sheet | What happened | Does it matter? |
|---|---|---|
| **Demand forecasts: Kestrel trusted about ±30%, Tern "lower confidence than Kestrel"** | Both were elicited and the ±30% was written into `BRIEF.md`, but demand was never varied. The one value-side sensitivity, year-1 against peak value per batch, changes what a batch is *worth*, not how many batches are *needed*. | **Yes, most of all.** The reservation's whole value sits in the both-succeed world, where Tern's shortfall is a race between demand and spare line capacity. A 30% higher Kestrel ramp takes line time from Tern and makes the reservation worth more. A 30% lower one frees the line and could make waiting enough. This is the same pattern as `field_service`: the numbers the informant labelled soft were recorded and then not tested. |
| **Marketed volume growth ~5%/yr, labelled ESTIMATE** | Used as a point value. The study lists it in the caveats as an input "without paperwork" but did not vary it. | **Yes, late in the horizon.** The turn-3 headroom table shows the in-house line running out of spare capacity by 2032 from marketed growth alone, so the build options and the late Tern ramp depend on it directly. The assistant flagged it as the input that matters most to that table, then did not re-run it. |
| **Capex Class 4 range $100–140M** | Recorded in the brief. Only the *phasing* was tested. | **No.** Even when both succeed, building now is $73–105M worse than waiting and building after the readouts $82–87M worse (8 paths each), so a $20M swing in capex cannot change its sign. Not testing it was correct, but the report does not say why. |
| **Effective capacity ~13/month includes two 2-week shutdowns a year** | Modelled as a flat 12.8 batches/month average. The shutdowns are not placed in the calendar. | **Possibly, at launch.** A shutdown that lands inside Kestrel's six-month stock build or Tern's PPQ window takes two weeks of capacity at exactly the contended moment. Averaging it away is the kind of simplification the "Kestrel never short at launch" finding rests on, and it was not examined. |
| **The two readouts "treated as independent", labelled ESTIMATE** | Taken as given, and every outcome weight (0.27 / 0.33 / 0.18 / 0.22) is built on it. | **Minor, but unexamined.** A positive correlation, such as a shared platform or a shared mechanism of failure, raises the both-succeed and both-fail weights together. The first is where the reservation earns and the second is where it costs a near-fixed $9–11M, so the net effect is ambiguous in sign and was never checked. |

The interview was strong on the regulatory mechanics and the contract terms, and on getting the *alternatives* right. It did not convert the operator's uncertainty labels into sensitivities, except for the Kestrel PoS that the operator had challenged out loud. This is the same gap the `field_service` audit found, and the same fix applies: the FRAME skill should ask which inputs were looked up and which were estimated, and pre-register the estimated ones as sensitivities.

## What it recovered correctly

Every structural fact, and every figure the recommendation turned on.

| Fact on the sheet | In the study |
|---|---|
| Kestrel readout ~March 2027, could slip 1–2 months, PoS 60% ("a bit generous") | readout month 5 + `k_slip` drawn 0–2; PoS 0.60 with a 0.50 re-test, proposed by the assistant in turn 2 after the operator called 0.60 generous |
| Tern readout ~Dec 2027 ±3 months (event-driven), PoS 45% | `t_shift` drawn −3…+3; PoS 0.45 |
| Readout → BLA ~3 months, 12-month review, launch ~1 month after approval | `K_LAUNCH = 21 + k_slip` (July 2028), Tern ~month 30 |
| **Filed site needs 3 PPQ batches before submission, or a ~6-month PAS during which that site cannot ship** — "people forget this" | Recovered in turn 2 and made central. It is the mechanism of both headline findings: waiting cannot use a late CMO for Tern without a PAS, and a reservation that makes the CMO Tern's only filed site creates the late-filing risk the review rule then removes |
| Kestrel PPQ on the own line in Q1 2027; Tern has no fill site locked | own line is Kestrel's filed site; Tern's filed site is the CMO when reserved, the own line otherwise |
| Nominal 16, effective ~13 (trailing 12.8) batches/month; marketed ~9/month +5%/yr; ~3% rejects | `line.init = 128` tenths, marketed demand grossed up for rejects, `REJECT_POS = 0.97` |
| Kestrel 1 → 3 at m12 → 4.5 peak; Tern 1 → 2.5 → 4; 6 months of forward stock at launch | the same piecewise ramps; launch-stock tier in the planner |
| Kestrel peak $1.2B, first 12 months ~$250M; Tern peak $0.9B, first 12 months ~$150M; 85% margin | value per tenth-batch derived from both, with year-1 as base and peak as a sensitivity. The assistant found that the sales and batch forecasts disagree by about 2× per batch (turn 3) and asked which to trust. Calibration reproduced $250.0M, 54 batches/yr and $151.1M |
| Supply-constrained first 6 months costs Kestrel ~15% of peak share permanently | `k_share = 0.85` latched by a rule. It never triggered, and the report says so |
| In-house $150k/batch; 9% WACC; horizon end-2032 | `varcost 15` per tenth; discounting 1.09^(t/12); `tspan` through month 75 |
| CMO: $5M/yr in advance, up to 4 batches/month, 3-year minimum, 12 months' notice forfeiting that year's fee (walk-away $5–10M), $400k/batch, 12–15 months tech transfer, offer holds to 30 Nov 2026, next slot not before mid-2028 | all encoded. `cmo_lead` drawn 12–15, walk-away $10M base and $5M sensitivity, and the "wait" arm takes the mid-2028 opening |
| Second line: ~$120M capex, 30–36 months, ~12 batches/month effective, ~$10M/yr fixed | `build_lead` drawn 30–36, `LINE2_TENTHS = 120`, capex phased 20/50/30, `l2_fixed` $10M/yr |
| **Operator's prior: reserve the CMO "sized for Tern"** | **Tested and overturned.** The study showed a Tern-only reservation is the worst CMO option without a second contract (Kestrel runs short on every both-succeed path), and the operator took it off the recommendation |
| **Operator's worry: is Kestrel alone fine on the in-house line once launch stock is counted?** | Answered directly: Kestrel is never short at launch on any path, under any option |
| Operator's opinion: the second line is too late to help either launch | Confirmed by the numbers, with the correct caveat that the 2032 horizon understates it |

Several facts the recommendation depends on are **not on the sheet**. The persona supplied them under questioning as elaborations consistent with it: the patient-first allocation rule that replaced the literal "marketed → Kestrel → Tern" order (turn 5), the possibility of a second CMO contract if both succeed, that the reserved slot can carry Kestrel if Tern fails, the $0.5M/batch cost of a backup qualification, and, most consequentially, **the review rule at Kestrel's readout** (turn 7). The review rule came from the operator, not the assistant. The assistant had proposed a backup qualification in the both-succeed world. The operator pointed out that Kestrel reads out nine months before Tern, so the reservation can be dropped cheaply when Kestrel fails. That one rule moved the expected value from roughly break-even (−$1.1M ± $3.4M when it was proposed, −$0.2M ± $1.5M once the Tern-only outcome had 10 paths) to +$3.8M ± $1.1M. This is the session's most instructive moment: the model made the *cost* of the rigid structure visible (a $72M tail on one Tern-only path), and the domain owner supplied the *option* that removed it. None of these elaborations can be used to score the interview. They show that the model was built to take a business rule as an input.

## Numerical provenance

No option-level ensembles were run (run budget: each 96-seed ensemble took 53 minutes, tracked as W24 in `spec/AGENT_WORKFLOW_CHARTER.md`). Every comparison is a set of paired single-seed paths with readout outcomes forced and weighted by the stated odds. There are 37 paths when both succeed, 30 for Tern only, and 3 each for Kestrel only and both fail. One 96-seed ensemble of the wait option was used as an independent cross-check. Its level figures agree, and three of its members, replayed on their exact seeds, reproduce it exactly. The coordinator re-derived the headline figures after each revision by running the study's own scripts (`analyze.py`, `paired.py`, `paired_both.py`, `exit_review.py`, `ens_check.py`) over the committed `runs/`, and every published number matched.

## The sheet, verbatim

## PRE-COMMITMENT FACT SHEET — Launch supply decision (fill-finish capacity)

**THIS SHEET MUST NOT BE MODIFIED AFTER IT IS WRITTEN (written 2026-09-24, before the first message to the assistant).** During the conversation I answer only from this sheet or from elaborations clearly consistent with it. Anything not here that I am asked, I give as a labelled estimate.

Role: VP Launch Readiness & Technical Operations, mid-sized commercial-stage biopharma (~$3B revenue, three marketed biologics). Executive committee (ExCo) decision meeting: ~5 weeks out (week of 2026-10-29).

### The two Phase-3 assets

- **Asset A, codename "Kestrel" (KST-301)**, subcutaneous mAb, immunology. Phase-3 topline readout expected **~March 2027** (about 6 months from now). [LOOKED-UP — clinical operations plan; could slip 1–2 months, ESTIMATE]
- **Kestrel probability of Phase-3 success: 60%.** [ESTIMATE — our portfolio committee's number; I think it is a bit generous]
- **Asset B, codename "Tern" (TRN-118)**, IV/SC biologic, oncology-adjacent. Phase-3 topline readout expected **~December 2027** (about 15 months from now). [LOOKED-UP — event-driven trial, so timing is soft; ±3 months ESTIMATE]
- **Tern probability of Phase-3 success: 45%.** [ESTIMATE — portfolio committee]
- The two readouts are treated as independent. [ESTIMATE]
- Readout to BLA submission: **~3 months.** [ESTIMATE, based on our last filing]
- BLA submission to approval: **~12 months standard review** (10 months + 60-day filing period). Priority review unlikely for either. [LOOKED-UP — FDA PDUFA timelines; priority review unlikely is ESTIMATE]
- So approvals if successful: **Kestrel ~June 2028**, **Tern ~March 2029**. [derived]
- Commercial launch ~1 month after approval. [ESTIMATE]

### Regulatory constraint that matters (people forget this)

- The commercial fill-finish site must be in the BLA with **3 process performance qualification (PPQ) batches** run at that site before submission (or added later via a Prior Approval Supplement, **~6 months FDA review**, during which that site cannot ship). [LOOKED-UP on the 6 months PAS review; that we must run PPQ before filing is our regulatory team's standard view]
- PPQ at our own line for Kestrel is already planned for Q1 2027 (in-house line is the filed site for Kestrel). [LOOKED-UP]
- Tern: no fill site locked in yet. [LOOKED-UP]

### In-house fill-finish line (single sterile line, vials)

- Nominal capacity: **16 batches/month.** [LOOKED-UP — site master plan]
- Effective capacity after changeovers, cleaning, environmental monitoring, planned shutdowns: **~13 batches/month** (two 2-week shutdowns per year are in that average). [LOOKED-UP — last 12 months actuals averaged 12.8]
- Committed to marketed portfolio today: **~9 batches/month.** [LOOKED-UP — S&OP plan]
- Marketed portfolio volume growth: **~5% per year.** [ESTIMATE — S&OP long-range]
- Spare today therefore **~4 batches/month**, shrinking to ~3 by 2029. [derived]
- Batch size: **~40,000 vials per batch** (both assets would use the same vial format). [LOOKED-UP for current products; ESTIMATE for the new assets]
- Batch failure / reject rate: **~3% of batches.** [LOOKED-UP]

### Demand after launch (batches per month on our line)

- **Kestrel**: ramps roughly linearly from ~1 batch/month at launch to ~3 batches/month at month 12, ~4.5 batches/month at peak (~year 3). [ESTIMATE — commercial forecast, which I trust about ±30%]
- **Tern**: ~1 batch/month at launch to ~2.5/month at month 12, ~4 batches/month at peak (~year 3). [ESTIMATE — commercial forecast, lower confidence than Kestrel]
- **Launch stock**: we need ~**6 months of forward demand in inventory at launch** (channel fill + safety stock). For Kestrel that is ~12 batches in the ~6 months before launch. [LOOKED-UP — our launch-readiness policy; batch count ESTIMATE]
- Drug substance is not the bottleneck; fill-finish is. [ESTIMATE, but I'm fairly confident]

### Money

- **Kestrel peak net sales ~$1.2B/year**; first-12-months revenue ~$250M. [ESTIMATE — commercial forecast]
- **Tern peak net sales ~$0.9B/year**; first-12-months revenue ~$150M. [ESTIMATE]
- Gross margin on both: **~85%.** [ESTIMATE — in line with our marketed biologics, LOOKED-UP for those]
- **Stock-out / supply-constrained launch cost**: lost sales for the months short, PLUS permanent share loss — a competitor (a large pharma) launches in Kestrel's class ~6–9 months after our expected approval; commercial's view is that a supply-constrained first 6 months costs **~15% of peak share permanently**. [ESTIMATE — commercial team's view, I think it's directionally right]
- Fill-finish cost per batch in-house: ~$150k variable. [LOOKED-UP]
- Discount rate: **9% WACC.** [LOOKED-UP — finance]
- Planning horizon: through **end of 2032**. [ESTIMATE — what ExCo usually looks at]

### Option 1 — Reserve a CMO fill-finish slot now

- CMO: a large established CDMO in Europe (term sheet in hand). [LOOKED-UP]
- **Reservation fee: $5M per year**, for dedicated capacity of **up to 4 batches/month**, **3-year minimum term**, payable annually in advance → **$15M committed**. [LOOKED-UP — term sheet]
- Early termination: allowed with 12 months' notice, but the fee for the notice year is forfeited; so if we walk away after a failed readout we lose roughly **$5–10M**. [LOOKED-UP on the terms; $ range derived]
- Per-batch fee on top: **~$400k per batch** filled. [LOOKED-UP]
- **Tech transfer + PPQ at the CMO: 12–15 months** from signing. [ESTIMATE — CMO says 12, our experience says 15]
- **Offer holds the slot until 30 November 2026**; after that they release it to another client, and the next available slot they quote is not before mid-2028. [LOOKED-UP — CMO BD contact, by email]

### Option 2 — Build a second in-house fill-finish line

- **Capex: ~$120M** (engineering pre-study range $100–140M). [ESTIMATE — engineering pre-study, Class 4]
- **Build + qualification: 30–36 months** to commercial-ready (including media fills and PPQ). If approved at ExCo, earliest ready ~Q2 2029. [ESTIMATE — engineering]
- Effective capacity: **~12 batches/month.** [ESTIMATE]
- Fixed operating cost once running: **~$10M/year.** [ESTIMATE]
- Board has been lukewarm on large capex this year. [LOOKED-UP — CFO told me]

### Option 3 — Wait for readouts

- No spend now. [fact]
- After Kestrel readout (~March 2027), CMO slot likely gone; next slot mid-2028 + 12–15 months tech transfer → CMO ready ~mid/late 2029 at earliest. [derived from CMO email + tech transfer estimate]
- If both succeed, the in-house line alone would have to carry marketed (~9–10/month) + Kestrel ramp + Tern ramp + both launch stock builds. [derived]

### My priors / opinions

- I am **leaning toward Option 1 (reserve the CMO slot now)**, sized for Tern, because the $5–10M downside if things fail looks small against a constrained launch, and waiting means lead times kill us. [opinion]
- I think Option 2 is too late to help either launch and is really a long-term portfolio question, not a launch question. [opinion]
- The CFO is pushing Option 3 ("don't pay for capacity on a 45% asset"). I need an answer that survives the CFO. [fact about the room]
- I'm not sure Kestrel alone is actually fine on the in-house line once you count launch stock build and the shutdowns — that's what I want checked. [opinion / worry]
- What I want for ExCo: a clear recommendation, what it costs, what it protects, and under what conditions it flips. [goal]

## What this sheet does and does not license you to conclude

It supports one narrow claim: **a decision owner with a regulatory constraint, a term sheet and a set of soft forecasts can be taken, by interview alone, to a checked model that reproduces their own calibration figures, overturns one of their priors (the Tern-only reservation) and confirms another (don't build), and whose recommendation survived three rounds of pushback without a number being retracted silently.** Every retraction in the transcript (the "0.18" late-filing claim in turn 6, the "worst case is $9–10M" framing in turn 6, and the "$1M understated" validation caveat in the final turn) is made out loud, with the evidence that forced it.

It does not show that the workflow handles a *resistant* informant. Like the `field_service` persona, this one answered in good faith and volunteered its uncertainty labels. Nor does it show the conclusion is robust to the untested soft inputs above: demand ±30%, marketed growth and readout correlation. A reader taking the recommendation into a room should treat those as the open questions.
