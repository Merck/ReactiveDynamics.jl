# Fill-finish capacity for Kestrel and Tern: what to take to ExCo

> The REPORT artifact (phase 4 of 4), final revision (turn 9). The recommendation comes first, then the numbers, then how they were computed. **Read this before any number below.** Each figure comes from single simulated paths, one per seed, with the readout outcomes forced. Paired paths per outcome: both succeed 37 (seeds 1–37); Tern only 30 (seeds 1–30, for the recommended rule); Kestrel only and both fail 3 each; the two build options 8 in the both-succeed outcome. The same seed gives every option the same readout dates, slips and lead times. The four outcomes are combined using the stated odds (Kestrel 0.60, Tern 0.45, independent: 0.27 / 0.33 / 0.18 / 0.22). The ± figures are standard errors across paths, and the probability ranges are 95% Wilson intervals on the path counts. The relayed 96-seed ensemble of the wait option serves as an independent check.

## **Reserve the CMO slot now, usable by either asset, and review at Kestrel's readout: +$3.8M ± $1.1M vs waiting, with a bounded cost in every outcome except both-succeed**

**Recommendation: sign the reservation before 30 November 2026, usable for either asset. At Kestrel's readout (March 2027), if Kestrel has failed, move Tern's PPQ in-house, file Tern on our own line, give the CMO notice and pay the forfeit. If Kestrel has succeeded, keep the reservation. Do not build a second line. Do not use the Tern-only structure unless the second CMO contract is in writing.**

**The wording the user will take to ExCo:**

> We recommend reserving the CMO fill-finish slot now, before 30 November 2026, usable by either Kestrel or Tern, and reviewing it at Kestrel's readout in March 2027. If Kestrel fails, Tern's validation batches move to our own line, Tern files there, and we give the CMO notice and pay the forfeit. If Kestrel succeeds, the reservation stands. Modelled over paired simulated paths, this is worth about +$3.8M (±$1.1M) against waiting at the stated odds, and about +$25M if a batch earns its peak-year value. At a 50% chance of Kestrel succeeding it is roughly break-even, around +$1.6M. Its purpose is supply protection. If both products succeed (27% chance at the stated odds), it is worth about +$39M and prevents a Tern launch shortfall that otherwise occurs on 26 of 37 paths, costing on average about 11 batches and $51M. If either product fails, the cost is a bounded $9–11M. We do not recommend building a second line. Even when both products succeed, it is $73–105M worse than waiting if built now and $82–87M worse if built after the readouts (8 paths each). We do not recommend a Tern-only reservation unless the second CMO contract is in writing. The one remaining risk is the CMO qualifying late: if both succeed, Tern's filing can be delayed, and 3 of 37 paths lose money, the worst about $14M.

- **Money is positive and passes the brief's test.** Against waiting: **+$3.8M ± $1.1M**, plausible range +$1.6M to +$6.0M, so Δ − 2·se > 0. The value sits in the both-succeed world. Every other world is a near-fixed cost, measured at $8.8M–9.7M. The ExCo wording says $9–11M, a conservative bound; see Caveats.

  | Outcome | Chance | Reserve now + review vs wait | Paths |
  |---|---|---|---|
  | Both succeed | 0.27 | +$38.7M ± $4.1M | 37 |
  | Kestrel only | 0.33 | −$8.8M ± $0.07M | 3 |
  | Tern only | 0.18 | −$9.4M ± $0.02M (range −$9.7M to −$9.3M) | 30 |
  | Both fail | 0.22 | −$9.4M ± $0.00M | 3 |

- **The review rule removes the Tern-only tail.**
  - Without the rule, reserving lost $31.7M ± $6.0M in the Tern-only world over 10 paths. Four of the 10 fell below −$20M, the worst lost $72.1M, and 4 of 10 filed late.
  - With the rule, the loss is all fee and forfeit. Tern's revenue matches waiting on all 30 paths, 0 of 30 file late, and 0 of 30 lose demand.
  - The rule never fires when Kestrel succeeds. Checked: it gives the same result as plain reserve-now in the both-succeed and Kestrel-only worlds.
- **What it buys.** If both succeed and we wait, Tern is short at launch on **26 of 37 paths**, which is **0.19 overall (0.15–0.22)**. With the reservation it is short on **0 of 37**.
- **What it does not buy: an on-time Tern filing.** When both succeed, Tern files late on 14 of 37 paths with the reservation (1–3 months) and on 13 of 37 under wait (2–3 months). Overall, the chance of a late Tern filing is 0.102 with reservation and review, and 0.095 without.
- **Kestrel is never short at launch** under either option (0 of 37 paths).
- **The remaining downside is a late CMO qualification.** 3 of 37 both-succeed paths lose money (−$14.3M, −$3.7M, −$5.5M). In each, Tern reads out before the CMO has qualified, and the CMO is Tern's only filed site.
- **Two confirmations the user holds.** Legal confirms that the slot can be reassigned; regulatory confirms that Tern can be filed on our line if Kestrel fails. Both are assumed by the model. If either reverses, the recommendation must be re-run.

![Options against waiting](figures/options_vs_wait.svg)

## What each option costs and protects against

eNPV is in $M, at 9% WACC through end-2032, at the year-1 value per batch.

| Option | vs wait, overall | vs wait, if both succeed | Tern short at launch, if both succeed | Tern files late, if both succeed | Kestrel short at launch |
|---|---|---|---|---|---|
| Wait | — | — | 26/37 | 13/37 (2–3 months) | 0/37 |
| **Reserve now + review at Kestrel readout** | **+3.8 ± 1.1** | **+38.7 ± 4.1** | **0/37** | 14/37 (1–3 months) | 0/37 |
| Reserve now, no review | −0.2 ± 1.5 (Tern only on 10 paths) | +38.7 ± 4.1 | 0/37 | 14/37 | 0/37 |
| Build after readouts | −23.0 (3 paths, turn 6) | −84.3 ± 0.6; range −86.8 to −82.3 (8 paths) | 6/8 | — | 0 |
| Build now | −105.8 (3 paths, turn 6) | −86.9 ± 3.1; range −104.7 to −72.8 (8 paths) | 6/8 | — | 0 |

**What a Tern shortfall costs.** When both succeed and we wait, Tern loses demand on 37 of 37 paths: on average 10.7 batches, worth $50.9M PV (range 5.9–20.3 batches, $27.8M–95.9M), over 5 to 10 months just after launch. With the reservation Tern loses no demand on any of the 37.

**Kestrel.** Under wait, Kestrel runs slightly short later in its life on 4 of 37 both-succeed paths (0.5 to 3.7 batches, never in its first six months). With the reservation, 0 of 37.

## Sensitivities

On the same paired paths, with the review rule:

| Case | Reserve now + review vs wait, overall | If both succeed |
|---|---|---|
| Base | +$3.8M ± $1.1M (+$1.6M to +$6.0M) | +$38.7M ± $4.1M |
| Value per batch at the peak rate | +$25.0M ± $2.8M (+$19.3M to +$30.7M) | +$117.3M ± $10.5M, negative on 1 of 37 |
| Kestrel PoS 0.50 (the user's own estimate) | +$1.6M ± $0.9M (−$0.2M to +$3.4M): roughly break-even, narrowly fails the test | +$38.7M ± $4.1M |

At the peak rate, the Tern-only loss without review would be −$51.5M on average (worst −$157.0M over 10 paths). With review it stays at −$9.4M. When both succeed, the build options stay negative even at the peak rate: −$88.7M and −$85.0M over 8 paths.

**Still on 3 paths from turn 6, indicative only, without the review rule:**
- No second CMO contract: reserve-for-both +$16.3M, reserve-for-Tern-only −$72.5M. Deferred at the user's direction, because the Tern-only structure is off the table without a written contract.
- Literal allocation rule: +$8.8M.
- A $5M walk-away: about +$1.0M.
- Even capex phasing: $0.95M worse.
- The 2-month Kestrel floor: no change.

## Cross-check against the 96-seed wait ensemble: resolved

- **Replays match exactly.** Ensemble members 20, 96 and 6 were replayed with their exact seeds, both unforced and forced to both-succeed. The two versions match each other byte for byte and reproduce the ensemble's late-filing months (2, 2, 0) and lost Tern demand (92, 59, 91 tenths).
- **The eNPV offset is harmless.** eNPV differs by a constant factor of 1.0072 = 1.09^(1/12), a one-month discounting-offset convention that cancels in every comparison.
- **The late-filing gap is sampling chance.** Among the ensemble's 20 both-succeed members, 2 filed late (2 months each); our paths show 13 of 37. Fisher exact two-sided p = 0.059: chance, not a model or forcing error.
- **Pooled late filing:** 15 of 57 (0.26, range 0.17–0.39).
- **Pooled Tern short at launch:** 40 of 57 (0.70), which is 0.19 overall (0.155–0.217), in line with our 0.19.
- **Scope:** the other turn-7 checks agreed (eNPV, P(Tern short), P(Kestrel short), Tern demand lost). The late-filing figures remain the least precise in this report. The ensemble covers only the wait option.

## What would change the answer

- **Kestrel's odds.** At PoS 0.50, which is the user's own estimate, the result is roughly break-even: +$1.6M ± $0.9M.
- **The value per batch.** At the peak rate the case is strong: +$25.0M.
- **Legal or regulatory reversing.** Either reassignment of the slot or filing Tern on our own line after a Kestrel failure. The review rule depends on the latter.
- **Whether a second CMO contract really exists.** It separates reserve-for-both from reserve-for-Tern-only.
- **A backup Tern qualification on our own line when both succeed.** This would address the remaining downside (3 of 37 paths). Not modelled, at the user's direction.

## Caveats

- **Horizon.** The model ends in 2032, so value earned after 2032 is not counted; the build options are understated as a result.
- **Validation cost, and a correction.** Revision 3 said the Tern-only cost was understated by about $1M because own-line PPQ is charged at $0.15M/batch rather than regulatory's $0.5M/batch. That overstated the gap. The Tern-only difference against waiting is all CMO fee plus forfeit (cost −$9.40M, Tern revenue +$0.00M). Under both wait and review, the own line is Tern's only filed site. So a higher own-line PPQ cost falls on both options and largely cancels. The $0.5M/batch applies to qualifying a backup site as well, which is not modelled. The ExCo wording's "$9–11M" is a conservative bound on the measured $8.8M–9.7M.
- **Small samples.** The Kestrel-only and both-fail worlds rest on 3 paths each. Their spread is under $0.1M, because the cost there is fee and forfeit.
- **Unverified inputs.** Two remain without paperwork: the second contract and +5%/yr marketed growth. Reassignment of the slot is now confirmed by legal, per the user.
- **Kestrel's share-loss penalty.** The 15% penalty never triggered.

## Why a spreadsheet cannot produce this

The reservation's value exists only in the joint outcome where both succeed. There it depends on a race between Tern's readout date and two other timelines: the CMO's qualification lead time, and Kestrel's six-month stock build on the same own line. The review rule works because Kestrel reads out nine months before Tern. It turns a commitment whose worst case depended on that race into one whose worst case is a known fee plus forfeit. Both effects come from month-by-month allocation under drawn dates, not from any per-asset expected value.

<details>
<summary>The model (click to expand)</summary>

The model is monthly, from 1 October 2026 to end-2032, with capacity counted in tenths of a batch. It includes the own line, the reserved CMO slot, a second CMO contract if both succeed, a possible second line, the readout draws, the regulatory clock including PAS, and the patient-first planner. It consists of 23 JSON documents from `gen_models.py` plus forced-outcome copies in `scenarios/`. Each option differs from waiting in one or two parameters. The recommended option is `arm_cmo_exit` (`cmo_now` = 1, `exit_at_k` = 1). Its rule `cmo_exit_k` fires once at Kestrel's readout if Kestrel failed. It releases the CMO slot and CMO PPQ, returns Tern to own-line PPQ and filing, and charges the forfeit.

</details>

## How the numbers were made

1. **Frame**: `BRIEF.md`.
2. **Author**: `gen_models.py`. All documents pass validation.
3. **Diagnose**: `DOCTOR.md`. `arm_cmo_exit` diagnoses 7 · 1 · 1, identical to `arm_cmo`.
4. **Report**: single-seed `rdx run` paths, every one `ok`, plus 6 ensemble-seed replays. The command, output excerpt and exit code for each are in `capture_log.md`, turns 6–9. The figure is `figures/options_vs_wait.json`, rendered by `agent/scripts/chartkit.jl`.

The headline numbers come from:

```bash
(cd agent/studies/launch_supply && python3 exit_review.py)            # recommended rule, base case
(cd agent/studies/launch_supply && python3 exit_review.py --peak)     # peak value per batch
(cd agent/studies/launch_supply && python3 exit_review.py --pk=0.5)   # Kestrel PoS 0.50
(cd agent/studies/launch_supply && python3 paired_both.py)            # both-succeed detail, build options
(cd agent/studies/launch_supply && python3 ens_check.py)              # ensemble cross-check
```

Each path is determined by (document, seed), and the verbatim outputs are in `runs/` and `runs_ens/`.
