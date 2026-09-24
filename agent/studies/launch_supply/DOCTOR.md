# DOCTOR — launch_supply

Every line below is machine-produced evidence from `rdx diagnose` / `rdx run` (full excerpts in `capture_log.md`, turn 4). Seed 1 throughout. Besides the 11 study documents, forced-state copies (same document, only `p_k`/`p_t` pinned to 0 or 1, and `k_first` to 0 where marked) were diagnosed so every mechanism is exercised in at least one run; they live in `/tmp/lstest/` and are never used for a quoted number.

## Command

```
julia --project=. agent/scripts/rdx.jl diagnose <doc> --seed=1 --expect-binds=own_free --expect-reusable=line,cmo,line2 \
  --tally=k_sold,t_sold,k_lost,t_lost,own_used,cmo_used,l2_used,mkt_done,k_age,t_age,t_file_delay,k_orders,t_orders
```

## Verdict

Every study arm: **7 pass · 1 warn · 1 fail** (`transitions_fire` fail, `horizon_sufficient` warn), except the four build arms which add a `bounded_pools` warn (6 · 2 · 1). The one fail is design-dead transitions within a single seed; the union over the forced-state runs fires all 27 transitions (table below), so no mechanism is silently missing. Cleared for REPORT with the caveats below.

| Document | Verdict | Ledger (seed 1, undiscounted $M) |
|---|---|---|
| model (wait) | 7 · 1 · 1 | reward 1704.0, cost 87.87 |
| arm_cmo | 7 · 1 · 1 | reward 1704.0, cost 98.07 (re-diagnosed turn 8; 97.98 was the turn-4 value, before the turn-5 fee timing fix) |
| arm_cmo_exit (turn 8) | 7 · 1 · 1 | reward 1704.0, cost 98.07 (identical to arm_cmo on seed 1: Kestrel succeeds, so the exit rule does not fire) |
| arm_cmo_tern_only | 7 · 1 · 1 | reward 1704.0, cost 98.08 |
| arm_cmo_walk5 | 7 · 1 · 1 | reward 1704.0, cost 97.98 |
| arm_build_now | 6 · 2 · 1 | reward 1704.0, cost 212.7 |
| arm_build_now_even | 6 · 2 · 1 | reward 1704.0, cost 212.7 |
| arm_build_after | 7 · 1 · 1 | reward 1704.0, cost 87.87 |
| model_pk50 / arm_cmo_pk50 | 7 · 1 · 1 | reward 1704.0, cost 87.87 / 97.98 |
| model_peak / arm_cmo_peak | 7 · 1 · 1 | reward 3483.0, cost 87.87 / 97.98 |
| model_both (forced) | 7 · 1 · 1 | reward 1806.0, cost 95.54 |
| model_both_kf0 (forced) | 7 · 1 · 1 | reward 1882.0, cost 95.54 |
| arm_cmo_both (forced) | 7 · 1 · 1 | reward 1990.0, cost 116.8 |
| arm_cmo_tern_only_both (forced) | 7 · 1 · 1 | reward 1488.0, cost 104.3 |
| arm_build_now_both (forced) | 6 · 2 · 1 | reward 2322.0, cost 232.0 |
| arm_build_after_both (forced) | 6 · 2 · 1 | reward 2145.0, cost 231.2 |
| model_tonly (forced, Tern only) | 7 · 1 · 1 | reward 741.7, cost 67.64 |
| arm_cmo_tonly (forced, Tern only) | 7 · 1 · 1 | reward 741.7, cost 91.52 |

Seed 1 draws Kestrel success and Tern failure, so the study arms differ only in cost on that seed; that is expected (the CMO/build is insurance whose payoff is in other states), not a defect.

## Checks

| Check | Result | Evidence |
|---|---|---|
| transitions_fire | fail (explained) | per-seed dead sets are the mechanisms of the state not drawn: Tern's production/ship/lose chain when Tern fails, line-2 capex/fixed/production outside the build arms, `cmo_fee_y2` outside the reserve-now arms, `k_buf_*` whenever `k_first = 1`. First run firing each: `t_*` chain `model_tonly` (t_make_own 566, t_ppq_own 31, t_ship 1178); `t_lose` and `t_delay_clock` `model_both` (935, 2); `k_lose` `arm_cmo_tern_only_both` (1036); `k_base_l2`, `l2_capex_y1/y2/y3` (12/12/12), `l2_fixed` (38) `arm_build_now`; `l2_capex_even` `arm_build_now_even` (36); `t_make_l2` `arm_build_now_both` (1170); `k_buf_own`/`k_buf_cmo` `model_both_kf0` (27/25); `k_buf_l2` `arm_build_now_both_kf0` (125); `cmo_fee_y2` `arm_cmo` (1). Never fired in any run: none |
| completions_exist | pass | model 14091, arm_cmo 14092, model_both 16570, arm_cmo_both 16991 successful completions |
| contention_exists | pass | pools run to zero incl. `own_free`, `cmo_free` |
| intended_resource_binds | pass | `own_free` runs to zero on every arm (pre-registered in BRIEF, see blind spot below) |
| modality_intent | pass | `line`, `cmo`, `line2` conserved on all input arcs |
| horizon_sufficient | warn (stands) | sales tallies and capacity-used tallies still rising in the final quarter |
| bounded_pools | pass / warn on build arms | build arms: `cmo_free` (and `line2` in `arm_build_after_both`) end at their maximum — see below |
| ledger_wired | pass | see table |
| stochasticity | pass | 10 of 27 transitions carry a success probability below 1 (the 3% reject rate on every batch) |

## Why the bind is checked on `own_free`, not on `line`

A batch holds its site for one month and gives it back inside the same month (conserved, cycletime 1), so the end-of-month snapshot always shows `line` at 128 tenths and the check can never see the line run dry. The planner rule therefore mirrors the capacity it left unallocated for next month into `own_free` / `cmo_free` / `l2_free`; `own_free = 0` means the line is fully booked. The direct evidence from actual firings agrees: on `model_both`, successful own-line tenths per month run 92–100 (marketed only) until the Kestrel stock build starts, then 127/128/128/127/126 and 125–128 for every month from t ≈ 27 to the horizon (128 available, the shortfall is the 3% rejects), and the CMO runs 37–40 of 40 from the month it qualifies.

## Warnings

- **horizon_sufficient** — the horizon is end-2032 by the user's decision rule. Both assets are still ramping or at peak at the horizon and a line approved after the readouts has run under two years by then, so every eNPV is a within-horizon number and the build arms' is a lower bound on their long-run value. Scope, not defect; carried into FINDINGS.
- **bounded_pools on the build arms** — `cmo_free` ends at its maximum because once line 2 is up, the cheaper line 2 takes the load and the mid-2028 CMO contract sits idle (the planner fills own → line 2 → CMO). That is real: in those arms the wait-path CMO is paid for and under-used. `line2` in `arm_build_after_both` is a capacity set once at qualification, not a queue. `k_orders`/`t_orders` (declared tallies here) are the current month's demand, drained the next month — they track the demand curve (Kestrel orders 12 → 45 tenths, Tern 14 → 40 on `model_both`), not a growing queue.

## Noticed by hand

- **Calibration** (forced both-succeed `arm_cmo`): Kestrel's first 12 shipping months 230 tenths = $250.0M revenue, last 12 months 540 tenths = 54 batches/yr; Tern's first 12 months 204 tenths = $151.1M. User: $250M, 54/yr at peak, $150M.
- **Kestrel is almost never short.** Kestrel lost demand in only one of the forced runs (`arm_cmo_tern_only_both`, 1036 tenths lost, all after its first six months). The priority rule (marketed → Kestrel → Tern) puts the whole shortage on Tern: `model_both` loses 935 tenths of Tern demand and files Tern 2 months late because Tern's PPQ waits behind Kestrel's stock build. So P(supply-constrained Kestrel launch) is likely near zero in every arm and the 15% share-loss mechanism may never trigger — to be read off the ensemble, not assumed.
- **Tern alone is not short.** With Kestrel failed, Tern fits on the own line under wait (`model_tonly` reward 741.7 = `arm_cmo_tonly` 741.7, `t_lose` dead in both). The option's value sits in the both-succeed state.
- Fixed during DIAGNOSE (AUTHOR re-entered): the first CMO fee and first capex month were not charged (rules act at the end of a tick, so capex ran 11 months in year 1); and a Tern-only reservation terminated after Tern fails had no fallback to the mid-2028 opening, which made `arm_cmo_tern_only` lose 953 of Kestrel's reward on seed 1 (751.1 vs 1704.0). Both fixed; post-fix capex runs 12/12/12 months and `arm_cmo_tern_only` matches `model` on reward (1704.0).

## Turn 5 re-check (patient-first allocation, 2-month floor, second CMO contract)

AUTHOR re-entered for the user's turn-5 corrections; 22 documents regenerated, all validate with `error_count 0`. Re-diagnosed `model`, `arm_cmo`, `arm_cmo_tern_only`, `arm_build_now` (seed 1) and forced both-succeed `model_both`, `arm_cmo_tern_only_both`, `arm_cmo_tern_only_no2_both`: every one 7·1·1 or 6·2·1 with the same explained `transitions_fire` fail and `horizon_sufficient` warn; `intended_resource_binds` passes on `own_free` everywhere. New mechanisms fire: second contract `cmo2_fee`, `k_base_cmo2`, `t_make_cmo2` on `model_both` (cmo2 used 1072 tenths). The new `bounded_pools` warn is `cmo2_used`, a capacity-used tally like `own_used`. `t_make_l2` and `l2_capex_even` did not fire in this subset (no forced build arm diagnosed this round); both fired before the turn-5 edit, which did not touch them, and `arm_build_now_both` shows line-2 use of 2228 tenths on `rdx run`.

Forced-state `rdx run`, seed 1 (tenths of a batch; full lines in `capture_log.md` turn 5):

- **The 2-month floor changes nothing on these paths**: `model_nofloor_both` and `arm_cmo_nofloor_both` reproduce the base runs exactly (reward 2356.0 and 2445.4). Kestrel's post-launch stock never falls below 91 tenths once launched with capacity available, since its launch build is never drawn down below the floor. The floor could only bind when Kestrel itself lacks capacity, and there it cannot be met anyway (`arm_cmo_tern_only_no2_both`: stock falls to 2).
- **Kestrel short: one setup.** Only `arm_cmo_tern_only_no2_both` — slot reserved for Tern, no second contract, both succeed — loses Kestrel demand: 1037 tenths, none in its first six months (`k_short 0`, `k_short_any 1`).
- **Tern late filing** only in own-line-PPQ arms (wait, build): 2 months on `model_both`, `arm_build_now_both`, `arm_build_after_both`; 0 on the reserve-now arms.
- The `*_nofloor` arms are therefore NOT sent to the ensemble; the floor's effect is reported from this evidence.

## Turn 6: evidence without ensembles

Ensembles were cancelled by the operator (see BRIEF, Analysis plan). Wherever the sections above say "to be read off the ensemble", read "odds-weighted forced-world paths" instead. The inputs behind every quoted number are the 44 forced-world copies in `scenarios/` (regenerated byte-for-byte by `python3 gen_models.py`, so not committed) and the 51 verbatim `rdx run` outputs in `runs/`. Each run returned `ok: true` and no findings, and `analyze.py` asserts that its recomputed reward and cost match every run's ledger. The structural checks above apply unchanged to the copies, because a copy differs from its arm only in the forced readout outcomes.

Seeds 2 and 3 in the both-succeed world refine two of the seed-1 statements in the turn-5 section:

- **Kestrel short.** The one setup where Kestrel loses material demand is still `arm_cmo_tern_only_no2`. It does so on all 3 seeds, losing 1037, 1016 and 1069 tenths, and never in Kestrel's first six months. Two smaller cases also appear. `arm_cmo_tern_only` (with a second contract) loses 8 tenths on seed 3 at the hand-over to the second contract. `model_no2` (wait, no second contract) loses 32 and 13 tenths on seeds 2 and 3.
- **Tern late filing is not only own-line PPQ.** On seed 2, all three reserve-now arms file Tern 3 months late (`t_file_delay 3`) even though Tern's PPQ is at the CMO. That path is consistent with the CMO's long lead-time draw (12–15 months) landing after Tern's readout + 3. The wait and build arms file 2 months late on seeds 1 and 3, when the PPQ waits behind Kestrel on the own line, and on time on seed 2.
- These are 3 paths per arm in that world and 1 in the others, so they show what can happen and roughly how often, not a measured probability with a spread.
