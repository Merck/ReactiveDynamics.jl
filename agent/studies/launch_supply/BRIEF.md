# BRIEF — Fill-finish capacity for the Kestrel and Tern launches

> The FRAME artifact. Written in the operator's register — batches, readouts, PPQ, term sheets — with every slot mapped onto the object model so AUTHOR can compile it without a second interview. Every fact below came from the user in turns 1–5 of the session (see `capture_log.md`); nothing is taken from another study. Status: signed off in turn 4 ("with those changes, go ahead and build it"), with the turn-4 and turn-5 changes folded in below.

## The decision

ExCo meets in about five weeks (early November 2026). Two Phase-3 biologics are coming: **Kestrel** (reads out ~March 2027, launch ~mid-2028) and **Tern** (reads out ~December 2027 ± 3 months, launch ~early-to-mid 2029). There is one in-house sterile fill-finish line, mostly committed to marketed products. The three options: (1) pay a CMO a reservation fee now to hold a slot that expires 30 November 2026; (2) build a second in-house line; (3) wait for the readouts and decide then. The binding problem is lead time: CMO tech transfer + PPQ (12–15 months), a new line (30–36 months) and a Prior Approval Supplement for any site not in the BLA (~6 months) are all longer than the readout-to-launch window, and launch requires six months of forward stock on the shelf.

**Deliverables.** (a) Expected NPV difference of each option against waiting, $M at 9% WACC through end-2032, with ± se; (b) the probability of a supply-constrained Kestrel or Tern launch under each option — ExCo's first question; (c) what would change the answer.

**What would make us act differently.** Take the CMO slot (or approve the build) if the eNPV it protects exceeds its cost and the gap holds at the low end of the uncertainty band (Δ − 2·se > 0). Re-test at Kestrel PoS 0.50 (the user's own belief is 0.50–0.55). If the gap does not hold, waiting is the answer and the CFO is right — that is a reportable finding.

## Why this is not a spreadsheet question

A per-asset expected value — the CFO's frame — prices each asset against capacity on its own. The exposure is in the joint outcomes: both readouts succeed (≈ 0.60 × 0.45 under independence) and their stock builds, Tern's PPQ batches and the first-year ramps land on the same line in the same months, with a regulatory rule (site must be in the BLA or wait ~6 months for a PAS) that decides which capacity is even usable when. When the line is short, the allocation rule (marketed → Kestrel → Tern) decides who goes without, and a Kestrel shortfall in its first six months costs 15% of peak share permanently. The value of reserved capacity is a marginal response to relieving a timed, contended constraint under two independent binary readouts with uncertain dates. No closed form.

**Qualifying test: three of three.** Something binds — the in-house line (≈12.8 effective batches/month, ~9 of them marketed and growing). Something takes time — tech transfer, PPQ, build, PAS, a six-month stock build. Something is decided mid-flight — who gets the line when Kestrel, Tern and marketed products all want it, and whether to sign/terminate/approve after each readout.

## A finding from the user's own numbers, before any model (turn 3)

Arithmetic on the user's figures only (logged in `capture_log.md`, turn 3): good-batch headroom on the in-house line after marketed demand (12.8 × 0.97 − marketed) is ≈2.6/month at Kestrel's launch and falls to ≈2.1 by mid-2029 and ≈1.1 by mid-2031, while Kestrel alone needs 3/month at month 12 and 4.5 at peak. **So the in-house line carries Kestrel alone only for roughly its first year; "it can carry one launch" does not hold past ~mid-2029 on these numbers.** Confirmed by the user in turn 4 (no debottleneck planned, no marketed volume moving), so the question is not "one launch or two" but "where does Kestrel's growth and all of Tern come from".

## Slots

Time origin t = 0 is 1 October 2026; one tick = one month; horizon t = 75 (end 2032). All money $M, all capacity batches (the model counts in tenths of a batch internally so fractional monthly rates are exact to 0.1 batch).

| Slot | Value (source: user) | Object-model target |
|---|---|---|
| **Scarce resource** | In-house line, effective 12.8 batches/month (nominal 16; ~13 after changeovers and shutdowns; trailing 12 months 12.8). Reusable: a batch occupies the line for its month and gives it back | plain place `line`, `init = 128` (tenths), held `conserved`, 10 units per batch, `cycletime = 1` |
| **Second resource (option 1)** | CMO slot, up to 4 batches/month, usable 12–15 months after signing (CMO says 12, experience says 15) | place `cmo` (40 tenths) created at `sign + U(12,15)` by a rule; `conserved` like the line |
| **Third resource (option 2)** | Second in-house line, commercial-ready 30–36 months after approval (earliest ~Q2 2029). 12 batches/month effective (turn 4) | place `line2` (120 tenths) created at `approve + U(30,36)` |
| **Background load** | Marketed products ~9 batches/month today, +5%/yr (S&OP estimate — soft), never bumped | transition `mkt`, integer monthly count from the cumulative curve `9·1.05^(t/12)/0.97` (grossed up for rejects), taken off the line first |
| **Rejects** | ~3% of batches rejected | `prob_of_success = 0.97` on every batch transition |
| **Readouts** | Kestrel PoS 0.60 (portfolio committee; user thinks generous → test 0.50), readout ~March 2027 slipping 0–2 months; Tern PoS 0.45, readout ~Dec 2027 ± 3 months (event-driven). Independent | params drawn once at t = 0 by rules (`Bernoulli`, `Uniform`) so every arm sees the same readout outcomes and dates per seed |
| **Regulatory clock** | Readout → +3 months BLA filing → +12 months standard review (no priority review) → +1 month launch. Filed fill site must have 3 PPQ batches before filing; otherwise the site is added by PAS (~6 months) after approval and cannot ship until approved | rules keyed on the drawn readout time; per-asset per-site "approved to ship" flags |
| **Kestrel's filed site** | Own line, PPQ planned Q1 2027 | Kestrel may ship from `line` from launch; from `cmo`/`line2` only after a PAS |
| **Tern's filed site** | None locked in. Without the CMO, Tern's 3 PPQ batches run on the own line between readout and filing, competing for the same batches as Kestrel's stock build | transition `tern_ppq_own` (3 batches, priority below Kestrel); filing waits until PPQ is done |
| **Launch stock** | Six months of forward demand on the shelf at launch: Kestrel ≈12 batches, Tern ≈9–10, built in the six months before launch | stock target = 6 × demand(month + 6), production starts 6 months before launch |
| **Demand (batches/month)** | Kestrel 1 at launch → 3 at month 12 (linear) → 4.5 at peak ~year 3 (±30%). Tern 1 → 2.5 at month 12 → 4 at peak ~year 3 (less certain) | `min(1 + a/6, 3 + (a−12)/16, 4.5)` and `min(1 + a/8, 2.5 + (a−12)/16, 4)`, a = months since launch |
| **Contention rule** | Marketed first, then Kestrel, then Tern, on every line; Tern PPQ on the own line competes with Kestrel (turn 3) | a `planner` rule at the end of each month allocates next month's free capacity in tiers — Kestrel, Tern PPQ, Tern stock, (Kestrel buffer top-up when `k_first = 0`) — each tier filled own line → line 2 → CMO (cheapest first). Turn 5 correction, base case `k_first = 0`: "the priority was about supplying patients" — Kestrel's launch build and each month's demand first, then Tern's PPQ and stock, then Kestrel's buffer top-up; the literal rule (`k_first = 1`, Kestrel's full six-month buffer before any Tern batch) is shown as a sensitivity. Guardrail in both: Kestrel never below 2 months of forward stock (`k_floor = 2`; checked against `k_floor = 0`) |
| **Shortfall cost** | Kestrel: short in first six months ⇒ 15% of peak share lost for good (competitor launches 6–9 months after our approval). Tern: shortfall is lost sales only. Unmet monthly demand is lost, not backlogged (turn 4: no delayed variant) | rule sets `k_share = 0.85` once if any Kestrel demand goes unmet with a < 6; unmet orders drain to a `*_lost` tally |
| **Value of a batch sold** | Kestrel peak $1.2B/yr, first 12 months $250M; Tern peak $0.9B/yr, first 12 months $150M; margin 85%. Base case: year-1 rate ($10.9M/batch Kestrel, $7.4M Tern, before margin); upside sensitivity: peak rate ($22.2M, $18.8M) — turn 4: "if the recommendation only holds at the peak rate, I need to know that" | `reward` per tenth-batch shipped on `k_sold`/`t_sold` |
| **Costs** | CMO: $5M/yr paid annually in advance, 3-year minimum, +$0.4M per batch filled; terminate on 12 months' notice forfeiting that year's fee (walk-away after a failed readout: $10M base = two annual fees, $5M sensitivity = the number quoted in the room, turn 4). In-house: $0.15M/batch. Second line: $120M capex (Class 4, $100–140M), phased 20/50/30 over three years from approval (evenly spread as a sensitivity, turn 4), $10M/yr fixed once operating | `cost` on a deep `cash` pool drawn by fee/capex/batch transitions |
| **Reassigning the CMO slot** | Term sheet reserves capacity, not a product — the user's reading, not confirmed. If Tern fails, Kestrel can move in via PAS (CMO won't be ready before Kestrel files) | tested both ways: `arm_cmo` (reassign allowed) and `arm_cmo_tern_only` |
| **Second CMO contract** | Turn 5: if both assets succeed, sign a second contract the same week on the same terms ($5M/yr, 3-year minimum, $0.4M/batch, 12–15 months to ready) — the user's estimate, no term sheet. Assumed not to open before mid-2028 like the first | place `cmo2` (40 tenths), signed at max(Tern readout, t = 20) to start the next month, ready after its own `U(12,15)` draw; both assets reach it by PAS. `second_cmo = 0` sensitivity |
| **Horizon and tick** | Through end-2032, monthly | `meta.tspan = 75.0`, `meta.dt = 1.0` |
| **Headline metrics** | eNPV at 9% WACC ($M); P(supply-constrained launch) per asset | `enpv --annual=0.09 --periods-per-year=12`; `final:k_short`, `final:t_short` (0/1 per run → ensemble mean is the probability) |

## The arms

All arms are the same network and differ only in `params` (diff any arm against `model.rdj.json` to see the lever). All business uncertainty is drawn first, at t = 0, identically in every arm, so every comparison is paired on the readouts and dates; effects are computed as paired per-seed differences. Base case in every arm (turn 5): patient-first allocation with the 2-month Kestrel floor, and a second CMO contract if both assets succeed.

- `model` — **wait**: nothing signed now; if either readout succeeds, sign the CMO's mid-2028 opening (t = 21), qualify it (12–15 months) and add it to each approved asset by PAS. Tern PPQs on the own line.
- `arm_cmo` — **reserve the CMO now for both** (sign by 30 Nov 2026, t = 2); Tern files with the CMO as its site; Kestrel may use the slot via PAS; if both fail the slot is terminated (walk-away $10M).
- `arm_cmo_tern_only` — **reserve now for Tern, add a slot later**: the reserved slot carries Tern only; Kestrel gets the second contract if both succeed, or the mid-2028 opening if Tern fails and Kestrel succeeds.
- `arm_build_now` — approve the second line at ExCo (t = 2); both assets reach it by PAS.
- `arm_build_after` — wait, and approve the build after Tern's readout if both succeed.
- Sensitivities: literal allocation rule (`*_literal`), no second contract (`*_no2`), Kestrel PoS 0.50 (`*_pk50`), peak-rate value per batch (`*_peak`); cost-only variants `arm_cmo_walk5`, `arm_cmo_tern_only_walk5`, `arm_build_now_even` (identical trajectories, different cash timing — their eNPV shift is computed exactly from the per-seed outcomes rather than re-simulated); `*_nofloor` checks whether the 2-month floor matters.

## Analysis plan

Ensembles with `--per-seed` on every simulated arm with one shared `--root-seed`, 96 seeds; the effect of each option is the mean of the per-seed difference against `model` (paired), ± its standard error, computed by a logged script. Metrics: eNPV at 9% (`enpv --annual=0.09 --periods-per-year=12`); `final:k_short` and `final:t_short` (P of a supply-constrained launch, first six months); `final:k_short_any` (Kestrel short at any time, turn 5); `final:t_file_delay` (months Tern's filing slips past readout + 3 because its filed site's PPQ is not done — P(> 0) and mean, turn 5); `final:k_lost`, `final:t_lost`; `final:k_success`/`final:t_success` (pairing check and state breakdown).

## What could go wrong (pre-registered, for DIAGNOSE)

- **The line doesn't bind** → no study. `--expect-binds=line`. Known blind spot: a batch holds the line for one month and returns it inside the same month, so the end-of-month snapshot always shows the line full and the check cannot see it bind; the model carries `own_used`/`cmo_used`/`l2_used` tallies as the direct evidence (capacity used per month vs available).
- **The line declared consumed instead of conserved** → capacity drains to zero and every option looks priceless. `--expect-reusable=line` (and `cmo`, `line2` on the arms that have them).
- **The wrong thing binds** — e.g. the `cash` pool or the order-routing plumbing rather than line capacity. `--expect-binds=line`; cash must be deep enough never to bind.
- **Dead transitions by design.** In any seed where an asset fails, its batch/ship transitions never fire; on `model` the CMO/line2 transitions never fire unless both succeed. Expect `transitions_fire` findings to be explained, not fixed.
- **Tallies that only grow** (`k_sold`, `t_sold`, `k_lost`, `t_lost`, launch flags). `--tally=` each.
- **Horizon truncation.** A second line approved after the readouts comes online ~2031 and runs ~18 months inside the window; every eNPV is a lower bound on long-run value, and the build arms especially. Expect `horizon_sufficient` to warn: scope, not defect.
- **Effect inside the noise.** The option value lives in the ≈27% of seeds where both succeed; paired draws help, but the CMO-vs-wait contrast may still be noisy. `significant_2se` on `rdx effect`; seeds ≥ 96.
- **Unit error.** Tenths of a batch internally; a factor-of-10 slip in `reward` would scale every NPV. Check `rdx run` ledger totals against the user's first-year and peak sales before any effect is quoted.

## What "done" looks like

For each option: Δ eNPV vs waiting ± se, its conservative bound Δ − 2·se, and P(short launch) for Kestrel and for Tern; the Kestrel-PoS-0.50 re-run of the headline contrast; and a stated list of what would flip the answer (CMO lead time, reassignment right, Tern timing, demand forecast, marketed growth).
