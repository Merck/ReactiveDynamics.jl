# capture log — launch_supply

Every command run during the live session, with excerpts of its output. Compute requests relayed back by the operator are recorded verbatim in their key lines.

### turn 1

No commands run. Read the operating procedure (`agent/skills/rd-study/SKILL.md`), the phase-1 skill (`agent/skills/rd-frame/SKILL.md`) and the exemplar brief `agent/studies/portfolio_bid/BRIEF.md`. Ran the qualifying test from the user's opening message (verdict: three of three, see reply) and opened phase 1 (FRAME) with the first round of questions.

### turn 2

User supplied odds (Kestrel 0.60, Tern 0.45, independent), sales/margin/WACC/horizon, lead times, CMO/build terms, and line numbers (nominal 16, effective ~13, trailing 12.8 batches/mo; marketed ~9/mo growing ~5%/yr; 40,000 vials/batch; ~3% reject). Back-of-envelope headroom check on the in-house line, from the user's numbers only:

```
python3 -c "for y,label in [(1.75,'mid-2028'),(2.5,'early-2029'),(6.25,'end-2032')]: m=9*1.05**y; print(label,'marketed %.1f/mo'%m,'headroom vs 12.8 = %.1f, vs 13 = %.1f'%(12.8-m,13-m))"
mid-2028 marketed 9.8/mo headroom vs 12.8 = 3.0, vs 13 = 3.2
early-2029 marketed 10.2/mo headroom vs 12.8 = 2.6, vs 13 = 2.8
end-2032 marketed 12.2/mo headroom vs 12.8 = 0.6, vs 13 = 0.8
exit 0
```

Opened FRAME round 2: batch demand per asset, allocation rule, Tern's filed site, the wait-path timing, and the decision threshold.

### turn 3

User supplied per-asset batch demand (Kestrel 1 → 3 at m12 → 4.5 peak; Tern 1 → 2.5 → 4), launch stock policy, share-loss scope (Kestrel only), allocation order (marketed → Kestrel → Tern), Tern-PPQ competition, decision rule, Kestrel PoS test at 0.50, both build timings, readout slip ranges; correction: take ~3% off for rejects.

Read engine semantics before mapping the brief (no model yet): `agent/README.md`, `agent/reference/model_schema.md`, `agent/skills/rd-author/SKILL.md`, `spec/CONTRACT_DRAFT.md` (deterministic rate floors; unmet spawn proposals are dropped, not queued — `src/solvers.jl` `evolve!`), `src/serialize.jl` (set_marking/set_params take value nodes).

Headroom vs demand check, from the user's numbers only (t = months from Oct 2026; Kestrel launch t=21, Tern launch t=30; good-batch headroom = 12.8×0.97 − marketed):

```
python3 -c "def dk(a): ...; def dt(a): ...; for t,lab in [...]: m=9*1.05**(t/12); good=12.8*0.97-m; ..."
Jul-2028 K launch    marketed  9.80  good-batch headroom 2.61  Kestrel need 1.00  Tern need 0.00  gap(K only) +1.61  gap(K+T) +1.61
Apr-2029 T launch    marketed 10.17  good-batch headroom 2.25  Kestrel need 2.50  Tern need 1.00  gap(K only) -0.25  gap(K+T) -1.25
Jul-2029 K+12        marketed 10.29  good-batch headroom 2.12  Kestrel need 3.00  Tern need 1.38  gap(K only) -0.88  gap(K+T) -2.25
Apr-2030 T+12        marketed 10.68  good-batch headroom 1.74  Kestrel need 3.56  Tern need 2.50  gap(K only) -1.82  gap(K+T) -4.32
Jul-2031 K+36 peak   marketed 11.35  good-batch headroom 1.07  Kestrel need 4.50  Tern need 3.44  gap(K only) -3.43  gap(K+T) -6.87
Dec-2032             marketed 12.21  good-batch headroom 0.21  Kestrel need 4.50  Tern need 4.00  gap(K only) -4.29  gap(K+T) -8.29
exit 0
```

Implied value per batch, sales forecast vs batch forecast:

```
python3 -c "k1=sum(dk(a) for a in range(12)); t1=sum(dt(a) for a in range(12)); ..."
Kestrel yr1 batches 23.0 -> 10.9 $M/batch; peak 54/yr -> 22.2 $M/batch
Tern    yr1 batches 20.2 -> 7.4 $M/batch; peak 48/yr -> 18.8 $M/batch
exit 0
```

Wrote `BRIEF.md` (draft, pending sign-off; three OPEN items: second-line capacity, $/batch reconciliation, capex phasing).

### turn 4

User confirmed the headroom table (no debottleneck planned, no marketed volume moving); changes: wait takes the mid-2028 CMO opening if EITHER readout succeeds; second line 12 batches/month effective; capex 20/50/30 (even spread as a sensitivity); year-1 value per batch as base, peak as upside; lost not delayed; walk-away $10M base and $5M sensitivity. "With those changes, go ahead and build it."

AUTHOR: `gen_models.py` writes 11 documents (same network, arms differ only in `params`). Production is allocated by an end-of-month `planner` rule in tiers (Kestrel → Tern PPQ on own line → Tern stock → Kestrel buffer top-up if `k_first = 0`), each tier filled own line → line 2 → CMO. Fixes made during authoring, each found by a run below: bare numbers in call args must be const nodes; transitions need `name`; `t_file_delay` off-by-one (`ge` → `gt`); whole-batch flooring wasted headroom (moved to tenths); first CMO fee and first capex month were never charged because rules act at the end of a tick (sign/approve rules now fire one month early with the start date set to the intended month); a Tern-only reservation terminated after Tern fails now falls back to the mid-2028 opening for Kestrel (otherwise `arm_cmo_tern_only` was a straw man, reward 751.1 vs 1704.0 on the same seed).

```
python3 agent/studies/launch_supply/gen_models.py
wrote 11 documents; value per tenth-batch: {'year1': (0.923913043478261, 0.6296296296296295), 'peak': (1.8888888888888888, 1.59375)}
exit 0
```

```
for f in agent/studies/launch_supply/*.rdj.json; do julia --project=. agent/scripts/rdx.jl validate $f; done
arm_build_after / arm_build_now / arm_build_now_even / arm_cmo / arm_cmo_peak / arm_cmo_pk50 / arm_cmo_tern_only / arm_cmo_walk5 / model / model_peak / model_pk50: "error_count": 0
exit 0 (each)
```

Calibration, forced both-succeed copy of `arm_cmo` (p_k = p_t = 1, /tmp/lstest/arm_cmo_both.rdj.json):

```
julia --project=. agent/scripts/rdx.jl run /tmp/lstest/arm_cmo_both.rdj.json --seed=1 --trace=k_sold,t_sold,k_age,t_age,own_used,own_free,cmo_used
K sold by month from first age>0: [0.0, 10.0, 12.0, 13.0, 15.0, 17.0, 18.0, 20.0, 22.0, 23.0, 25.0, 27.0, 28.0, 30.0]
T sold by month: [0.0, 10.0, 11.0, 13.0, 14.0, 15.0, 16.0, 18.0, 19.0, 20.0, 21.0, 23.0, 24.0, 25.0]
K 12 shipping months 230.0 revenue $M 250.0 | T 204.0 revenue $M 151.1
K last-12-month sold tenths 540.0
ledger {'reward': 1989.5475040257652, 'cost': 116.81999999999996, 'net': 1872.7275040257653}
exit 0
```

(User: Kestrel ~$250M first 12 months, peak 54 batches/yr ↔ 540 tenths; Tern ~$150M.)

Line binds, forced both-succeed wait (`model_both`), capacity used per month (tenths, successful batches; 128 available):

```
julia --project=. agent/scripts/rdx.jl run /tmp/lstest/model_both.rdj.json --seed=1 --trace=...,own_used,cmo_used,...
model_both own/month [92.0, 94.0, 94.0, 94.0, 94.0, 95.0, 95.0, 96.0, 96.0, 96.0, 97.0, 98.0, 97.0, 98.0, 99.0, 98.0, 100.0, 127.0, 128.0, 128.0, 127.0, 126.0, 117.0, 103.0, 113.0, 121.0, 126.0, 128.0, 127.0, ... 128.0]
  cmo/month [0.0 ×40, 39.0, 40.0, 40.0, 38.0, ...]
exit 0
```

Forced-state runs, seed 1 (spawns in tenths; lost = unmet demand in tenths):

```
model_both            ledger reward 1806.3 cost 90.6   K y1 230  T y1 49   klost 0  tlost 935   t_file_delay 2  k_short 0 t_short 1
model_both_kf0        ledger reward 1881.9 cost 90.6   K y1 230  T y1 142  klost 0  tlost 815   t_file_delay 2  k_short 0 t_short 1
arm_cmo_both          ledger reward 1983.9 cost 111.7  K y1 230  T y1 204  klost 0  tlost 733   t_file_delay 0  k_short 0 t_short 0
arm_build_now_both    ledger reward 2322.0 cost 225.1  K y1 230  T y1 88   klost 0  tlost 116   t_file_delay 2  k_short 0 t_short 1
arm_build_after_both  ledger reward 2148.8 cost 224.1  K y1 230  T y1 54   klost 0  tlost 391   t_file_delay 2  k_short 0 t_short 1
arm_cmo_tern_only_both ledger reward 1490.1 cost 99.0  K y1 230  T y1 204  klost 1034 tlost 0   t_file_delay 0  k_short 0 t_short 0
exit 0 (each)
```

(These predate the fee/capex timing fix and the Tern-only fallback; rewards are unchanged by those fixes, costs rose — see the diagnose ledgers below.)

DIAGNOSE, every arm plus forced-state copies, seed 1:

```
bash /tmp/lstest/diag_all.sh <arm> ...   # rdx diagnose <arm> --seed=1 --expect-binds=own_free --expect-reusable=line,cmo,line2 --tally=k_sold,t_sold,k_lost,t_lost,own_used,cmo_used,l2_used,mkt_done,k_age,t_age,t_file_delay,k_orders,t_orders
arm_build_after | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 87.87 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_build_after_both | 6 pass · 2 warn · 1 fail | reward 2145.0, cost 231.2 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in line2, cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or a
arm_build_now | 6 pass · 2 warn · 1 fail | reward 1704.0, cost 212.7 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an intak
arm_build_now_both | 6 pass · 2 warn · 1 fail | reward 2322.0, cost 232.0 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an intak
arm_build_now_both_kf0 | 6 pass · 2 warn · 1 fail | reward 2356.0, cost 233.2 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an intak
arm_build_now_even | 6 pass · 2 warn · 1 fail | reward 1704.0, cost 212.7 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an intak
arm_cmo | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 97.98 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_both | 7 pass · 1 warn · 1 fail | reward 1990.0, cost 116.8 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_konly | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 97.98 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_peak | 7 pass · 1 warn · 1 fail | reward 3483.0, cost 97.98 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_pk50 | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 97.98 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 98.08 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only_both | 7 pass · 1 warn · 1 fail | reward 1488.0, cost 104.3 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only_konly | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 98.08 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only_tonly | 7 pass · 1 warn · 1 fail | reward 741.7, cost 91.52 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tonly | 7 pass · 1 warn · 1 fail | reward 741.7, cost 91.52 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_walk5 | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 97.98 | nonpass=['transitions_fire', 'horizon_sufficient']
model | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 87.87 | nonpass=['transitions_fire', 'horizon_sufficient']
model_both | 7 pass · 1 warn · 1 fail | reward 1806.0, cost 95.54 | nonpass=['transitions_fire', 'horizon_sufficient']
model_both_kf0 | 7 pass · 1 warn · 1 fail | reward 1882.0, cost 95.54 | nonpass=['transitions_fire', 'horizon_sufficient']
model_konly | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 87.87 | nonpass=['transitions_fire', 'horizon_sufficient']
model_peak | 7 pass · 1 warn · 1 fail | reward 3483.0, cost 87.87 | nonpass=['transitions_fire', 'horizon_sufficient']
model_pk50 | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 87.87 | nonpass=['transitions_fire', 'horizon_sufficient']
model_tonly | 7 pass · 1 warn · 1 fail | reward 741.7, cost 67.64 | nonpass=['transitions_fire', 'horizon_sufficient']
never fired in any run: []
exit 1 (each: transitions_fire fails on design-dead transitions; see DOCTOR.md)
```

### turn 5

User: base case = patient-first allocation (Kestrel launch build + monthly demand, then Tern, then Kestrel buffer), literal rule as sensitivity; guardrail Kestrel never below ~2 months forward stock in both (drop if complicated, but report whether it matters); allow a second CMO contract on the same terms if both succeed (user's estimate), no-second-contract as sensitivity; flag the setup where Kestrel runs short at any time; report odds Tern files late because its PPQ is stuck behind Kestrel, per option. "Run it."

AUTHOR: `k_first` default 0, new `k_floor` (2), `second_cmo` (1), places `cmo2`, `k_short_any`, draw `cmo2_lead` appended last in the draws rule (earlier draws unchanged). 22 documents.

```
python3 agent/studies/launch_supply/gen_models.py
wrote 22 documents; value per tenth-batch: {'year1': (0.923913043478261, 0.6296296296296295), 'peak': (1.8888888888888888, 1.59375)}
exit 0
for f in agent/studies/launch_supply/*.rdj.json; do rdx validate $f; done
all 22: "error_count": 0, exit 0
```

Forced-state runs (/tmp/lstest copies with p_k/p_t pinned), seed 1:

```
bash /tmp/lstest/run_all.sh ...   # rdx run <doc> --seed=1 --trace=k_sold,t_sold,k_lost,t_lost,k_age,t_age,k_stock,t_stock,own_used,cmo_used,cmo2_used,l2_used,t_file_delay
arm_build_after_both         rew  2356.0 cost  272.0 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 142.0 Tlost 62 Tshort 1 Tdelay 2 | cmo 785 cmo2 475 l2 1228 fees 5+5 | Kstock min post-launch 92.0
arm_build_now_both           rew  2356.0 cost  259.9 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 142.0 Tlost 62 Tshort 1 Tdelay 2 | cmo 193 cmo2 67 l2 2228 fees 5+5 | Kstock min post-launch 91.0
arm_cmo_both                 rew  2445.4 cost  179.9 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1678 cmo2 924 l2 0 fees 7+5 | Kstock min post-launch 121.0
arm_cmo_konly                rew  1703.7 cost   98.1 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@None Ty1 None Tlost 0 Tshort 0 Tdelay 0 | cmo 1234 cmo2 0 l2 0 fees 7+0 | Kstock min post-launch 121.0
arm_cmo_literal_both         rew  2445.4 cost  180.0 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1675 cmo2 929 l2 0 fees 7+5 | Kstock min post-launch 121.0
arm_cmo_nofloor_both         rew  2445.4 cost  179.9 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1678 cmo2 924 l2 0 fees 7+5 | Kstock min post-launch 121.0
arm_cmo_tern_only_both       rew  2445.4 cost  179.8 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1368 cmo2 1241 l2 0 fees 7+5 | Kstock min post-launch 46.0
arm_cmo_tern_only_konly      rew  1703.7 cost   98.4 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@None Ty1 None Tlost 0 Tshort 0 Tdelay 0 | cmo 1237 cmo2 0 l2 0 fees 7+0 | Kstock min post-launch 92.0
arm_cmo_tern_only_no2_both   rew  1487.3 cost  104.0 K@24 Ky1 230.0 Klost 1037 Kshort6 0 Kany 1 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1370 cmo2 0 l2 0 fees 7+0 | Kstock min post-launch 2.0
arm_cmo_tern_only_tonly      rew   741.7 cost   91.5 K@None Ky1 None Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1375 cmo2 0 l2 0 fees 7+0 | Kstock min post-launch None
arm_cmo_tonly                rew   741.7 cost   91.5 K@None Ky1 None Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 1375 cmo2 0 l2 0 fees 7+0 | Kstock min post-launch None
model_both                   rew  2356.0 cost  164.5 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 142.0 Tlost 62 Tshort 1 Tdelay 2 | cmo 1395 cmo2 1072 l2 0 fees 5+5 | Kstock min post-launch 92.0
model_konly                  rew  1703.7 cost   88.1 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@None Ty1 None Tlost 0 Tshort 0 Tdelay 0 | cmo 1235 cmo2 0 l2 0 fees 5+0 | Kstock min post-launch 93.0
model_literal_both           rew  2334.0 cost  163.1 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 107.0 Tlost 97 Tshort 1 Tdelay 2 | cmo 1399 cmo2 1033 l2 0 fees 5+5 | Kstock min post-launch 94.0
model_no2_both               rew  1881.9 cost   95.5 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 142.0 Tlost 815 Tshort 1 Tdelay 2 | cmo 1408 cmo2 0 l2 0 fees 5+0 | Kstock min post-launch 91.0
model_nofloor_both           rew  2356.0 cost  164.5 K@24 Ky1 230.0 Klost 0 Kshort6 0 Kany 0 | T@37 Ty1 142.0 Tlost 62 Tshort 1 Tdelay 2 | cmo 1395 cmo2 1072 l2 0 fees 5+5 | Kstock min post-launch 92.0
model_tonly                  rew   741.7 cost   67.6 K@None Ky1 None Klost 0 Kshort6 0 Kany 0 | T@35 Ty1 204.0 Tlost 0 Tshort 0 Tdelay 0 | cmo 825 cmo2 0 l2 0 fees 5+0 | Kstock min post-launch None
exit 0 (each)
```

```
bash /tmp/lstest/diag_all.sh model arm_cmo arm_cmo_tern_only arm_build_now model_both arm_cmo_tern_only_both arm_cmo_tern_only_no2_both
arm_build_now | 6 pass · 2 warn · 1 fail | reward 1704.0, cost 212.7 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo_free: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an intak
arm_cmo | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 98.07 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 98.39 | nonpass=['transitions_fire', 'horizon_sufficient']
arm_cmo_tern_only_both | 6 pass · 2 warn · 1 fail | reward 2445.0, cost 179.8 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo2_used: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an inta
arm_cmo_tern_only_no2_both | 7 pass · 1 warn · 1 fail | reward 1487.0, cost 104.0 | nonpass=['transitions_fire', 'horizon_sufficient']
model | 7 pass · 1 warn · 1 fail | reward 1704.0, cost 88.1 | nonpass=['transitions_fire', 'horizon_sufficient']
model_both | 6 pass · 2 warn · 1 fail | reward 2356.0, cost 164.5 | nonpass=['transitions_fire', 'horizon_sufficient', 'bounded_pools']
      bounded_pools Unbounded growth in cmo2_used: each ends at its maximum, well above its initial level, and is still rising. Check for a missing downstream transition or an inta
never fired in any run: ['l2_capex_even', 't_make_l2']
exit 1 (each; transitions_fire design-dead, see DOCTOR.md)
```

Cost-only variants, exact PV shift at 9% using the same discounting as `rdx enpv` (flow / 1.09^(t/12)):

```
python3 -c "d=lambda t:1.09**(-t/12); ... "
walk-away 5M saving PV at t=14: 4.522
capex PV phased 20/50/30: 103.64  even: 104.60  diff (even - phased cost): 0.95
exit 0
```

### turn 6

Compute relay from the operator, not from the user. The ensemble batch was cancelled after the first of 15 ensembles, because each 96-seed ensemble took 53 minutes on this machine: `rdx ensemble` rebuilds and recompiles the model for every seed. The operator's decision was "no ensembles. Keep it simple." Instead, use single-seed `rdx run` calls in the foreground with the same seeds across options and forced readout worlds, combined with the odds (Kestrel 0.60, Tern 0.45) by arithmetic, stated plainly, with no spread that was not measured. The relay pointed to the one completed ensemble (baseline `model.rdj.json`, 96 seeds, root seed 2609) as a file under `~/rd_relay/`. **It was not read.** Blind rule 1 forbids reading anything under that path, so no number from it appears here or in FINDINGS.md.

Forced-world copies of the 11 decision arms (the 5 options, plus the `_literal` and `_no2` variants of `model`/`arm_cmo`/`arm_cmo_tern_only`) were generated. Each world sets the readout outcomes: both = Kestrel ok and Tern ok; konly = Kestrel ok and Tern fails; tonly = Kestrel fails and Tern ok; neither = both fail.

```
python3 agent/studies/launch_supply/gen_models.py
wrote 22 documents + 44 scenario copies; value per tenth-batch: {'year1': (0.923913043478261, 0.6296296296296295), 'peak': (1.8888888888888888, 1.59375)}
exit 0
```

Each path was run one at a time in the foreground by `run_scenarios.sh ARM WORLD SEED ...`. That script calls `rdx run scenarios/<arm>__<world>.rdj.json --seed=S --trace=cash,k_sold,t_sold,k_lost,t_lost,k_age,t_age,k_stock,t_stock,k_short,k_short_any,t_short,t_file_delay,own_used,cmo_used,cmo2_used,l2_used` and saves the verbatim JSON output to `runs/<arm>__<world>__s<S>.json`. Batches were kept to 12 or fewer runs. Every batch line reported `exit 0`. The coverage is shown below, verified over the saved outputs:

```
python3 - (parse every runs/*.json, check "ok")
51 run files, unparseable: []
{('both', 's1'): 11, ('both', 's2'): 11, ('both', 's3'): 11, ('konly', 's1'): 8, ('neither', 's1'): 5, ('tonly', 's1'): 5}
each: {'error': None, 'ok': True, 'command': 'run', 'findings': []}
exit 0
```

Per-path summary of the both-succeed world. The eNPV is recomputed from the traced flows by `analyze.load`, and launch months are the first month with k_age/t_age > 0:

```
arm_build_after__both__s1  npv 1380.1 Klaunch 24 Tlaunch 37 Tdelay 2 Kshort_any 0 Klost 0 Tlost 62 Tshort 1
arm_build_after__both__s2  npv 1447.3 Klaunch 25 Tlaunch 30 Tdelay 0 Kshort_any 0 Klost 0 Tlost 114 Tshort 0
arm_build_after__both__s3  npv 1430.9 Klaunch 23 Tlaunch 35 Tdelay 2 Kshort_any 0 Klost 0 Tlost 92 Tshort 1
arm_build_now__both__s1    npv 1376.7 Klaunch 24 Tlaunch 37 Tdelay 2 Kshort_any 0 Klost 0 Tlost 62 Tshort 1
arm_build_now__both__s2    npv 1446.8 Klaunch 25 Tlaunch 30 Tdelay 0 Kshort_any 0 Klost 0 Tlost 112 Tshort 0
arm_build_now__both__s3    npv 1429.4 Klaunch 23 Tlaunch 35 Tdelay 2 Kshort_any 0 Klost 0 Tlost 92 Tshort 1
arm_cmo__both__s1          npv 1516.2 Klaunch 24 Tlaunch 35 Tdelay 0 Kshort_any 0 Klost 0 Tlost 0 Tshort 0
arm_cmo__both__s2          npv 1518.7 Klaunch 25 Tlaunch 33 Tdelay 3 Kshort_any 0 Klost 0 Tlost 0 Tshort 0
arm_cmo__both__s3          npv 1581.2 Klaunch 23 Tlaunch 33 Tdelay 0 Kshort_any 0 Klost 0 Tlost 0 Tshort 0
arm_cmo_tern_only__both__s1 npv 1516.6 Klaunch 24 Tlaunch 35 Tdelay 0 Kshort_any 0 Klost 0 Tlost 0 Tshort 0
arm_cmo_tern_only__both__s2 npv 1519.0 Klaunch 25 Tlaunch 33 Tdelay 3 Kshort_any 0 Klost 0 Tlost 0 Tshort 0
arm_cmo_tern_only__both__s3 npv 1576.2 Klaunch 23 Tlaunch 33 Tdelay 0 Kshort_any 1 Klost 8 Tlost 0 Tshort 0
arm_cmo_tern_only_no2__both__s1 npv 956.0 Klaunch 24 Tlaunch 35 Tdelay 0 Kshort_any 1 Klost 1037 Tlost 0 Tshort 0
arm_cmo_tern_only_no2__both__s2 npv 971.6 Klaunch 25 Tlaunch 33 Tdelay 3 Kshort_any 1 Klost 1016 Tlost 0 Tshort 0
arm_cmo_tern_only_no2__both__s3 npv 1000.9 Klaunch 23 Tlaunch 33 Tdelay 0 Kshort_any 1 Klost 1069 Tlost 0 Tshort 0
arm_cmo_no2__both__s1      npv 1291.0 ... Tlost 688 ; s2 1283.1 Tdelay 3 Tlost 715 ; s3 1338.6 Tlost 733 (Kshort_any 0 all)
model__both__s1            npv 1464.3 Klaunch 24 Tlaunch 37 Tdelay 2 Kshort_any 0 Klost 0 Tlost 62 Tshort 1
model__both__s2            npv 1533.1 Klaunch 25 Tlaunch 30 Tdelay 0 Kshort_any 0 Klost 0 Tlost 110 Tshort 0
model__both__s3            npv 1516.0 Klaunch 23 Tlaunch 35 Tdelay 2 Kshort_any 0 Klost 0 Tlost 92 Tshort 1
model_no2__both__s1        npv 1214.7 Tdelay 2 Kshort_any 0 Klost 0 Tlost 815 Tshort 1
model_no2__both__s2        npv 1186.3 Tdelay 0 Kshort_any 1 Klost 32 Tlost 1060 Tshort 0
model_no2__both__s3        npv 1237.8 Tdelay 2 Kshort_any 1 Klost 13 Tlost 903 Tshort 1
model_literal__both__s1..3 npv 1449.3 / 1488.8 / 1487.6, Tlost 97 / 210 / 157
exit 0
```

Scenario weighting by `analyze.py`. The ledger check asserts that the recomputed reward and cost match each run's ledger. The seed mean is taken within each world, and the worlds are weighted by the odds. Worlds that are missing for a sensitivity arm are borrowed from its base arm and marked `(= base)`. For `_no2` this holds by construction, because the second contract can only be signed when both assets succeed. For `_literal`, the tonly and neither worlds are borrowed; with one asset or none there is no contention, so the literal rule and the patient-first rule should coincide, and konly bears that out (1096.3 vs 1096.2).

```
python3 agent/studies/launch_supply/analyze.py
weights (p_k=0.6, p_t=0.45): {'both': 0.27, 'konly': 0.33, 'tonly': 0.18, 'neither': 0.22} | value per tenth-batch (0.923913043478261, 0.6296296296296295)
per world, mean over seeds shown (single paths when one seed): eNPV $M | K short launch / K short ever / T short launch / T filing delay (months)
  both     model                        seeds s1,s2,s3 eNPV   1504.5 | 0.00 0.00 0.67 1.33
  both     arm_cmo                      seeds s1,s2,s3 eNPV   1538.7 | 0.00 0.00 0.00 1.00
  both     arm_cmo_tern_only            seeds s1,s2,s3 eNPV   1537.2 | 0.00 0.33 0.00 1.00
  both     arm_build_now                seeds s1,s2,s3 eNPV   1417.6 | 0.00 0.00 0.67 1.33
  both     arm_build_after              seeds s1,s2,s3 eNPV   1419.4 | 0.00 0.00 0.67 1.33
  both     model_literal                seeds s1,s2,s3 eNPV   1475.2 | 0.00 0.00 0.67 1.33
  both     arm_cmo_literal              seeds s1,s2,s3 eNPV   1538.8 | 0.00 0.00 0.00 1.00
  both     arm_cmo_tern_only_literal    seeds s1,s2,s3 eNPV   1538.8 | 0.00 0.33 0.00 1.00
  both     model_no2                    seeds s1,s2,s3 eNPV   1212.9 | 0.00 0.67 0.67 1.33
  both     arm_cmo_no2                  seeds s1,s2,s3 eNPV   1304.2 | 0.00 0.00 0.00 1.00
  both     arm_cmo_tern_only_no2        seeds s1,s2,s3 eNPV    976.2 | 0.00 1.00 0.00 1.00
  konly    model                        seeds s1       eNPV   1096.2 | 0.00 0.00 0.00 0.00
  konly    arm_cmo                      seeds s1       eNPV   1087.4 | 0.00 0.00 0.00 0.00
  konly    arm_cmo_tern_only            seeds s1       eNPV   1086.6 | 0.00 0.00 0.00 0.00
  konly    arm_build_now                seeds s1       eNPV    989.7 | 0.00 0.00 0.00 0.00
  konly    arm_build_after              seeds s1       eNPV   1096.2 | 0.00 0.00 0.00 0.00
  konly    model_literal                seeds s1       eNPV   1096.3 | 0.00 0.00 0.00 0.00
  konly    arm_cmo_literal              seeds s1       eNPV   1087.5 | 0.00 0.00 0.00 0.00
  konly    arm_cmo_tern_only_literal    seeds s1       eNPV   1086.8 | 0.00 0.00 0.00 0.00
  konly    model_no2                    seeds s1       eNPV   1096.2 | 0.00 0.00 0.00 0.00  (= model)
  konly    arm_cmo_no2                  seeds s1       eNPV   1087.4 | 0.00 0.00 0.00 0.00  (= arm_cmo)
  konly    arm_cmo_tern_only_no2        seeds s1       eNPV   1086.6 | 0.00 0.00 0.00 0.00  (= arm_cmo_tern_only)
  tonly    model                        seeds s1       eNPV    438.9 | 0.00 0.00 0.00 0.00
  tonly    arm_cmo                      seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00
  tonly    arm_cmo_tern_only            seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00
  tonly    arm_build_now                seeds s1       eNPV    328.0 | 0.00 0.00 0.00 0.00
  tonly    arm_build_after              seeds s1       eNPV    438.9 | 0.00 0.00 0.00 0.00
  tonly    model_literal                seeds s1       eNPV    438.9 | 0.00 0.00 0.00 0.00  (= model)
  tonly    arm_cmo_literal              seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00  (= arm_cmo)
  tonly    arm_cmo_tern_only_literal    seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00  (= arm_cmo_tern_only)
  tonly    model_no2                    seeds s1       eNPV    438.9 | 0.00 0.00 0.00 0.00  (= model)
  tonly    arm_cmo_no2                  seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00  (= arm_cmo)
  tonly    arm_cmo_tern_only_no2        seeds s1       eNPV    420.3 | 0.00 0.00 0.00 0.00  (= arm_cmo_tern_only)
  neither  model                        seeds s1       eNPV      0.0 | 0.00 0.00 0.00 0.00
  neither  arm_cmo                      seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00
  neither  arm_cmo_tern_only            seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00
  neither  arm_build_now                seeds s1       eNPV   -123.9 | 0.00 0.00 0.00 0.00
  neither  arm_build_after              seeds s1       eNPV      0.0 | 0.00 0.00 0.00 0.00
  neither  model_literal                seeds s1       eNPV      0.0 | 0.00 0.00 0.00 0.00  (= model)
  neither  arm_cmo_literal              seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00  (= arm_cmo)
  neither  arm_cmo_tern_only_literal    seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00  (= arm_cmo_tern_only)
  neither  model_no2                    seeds s1       eNPV      0.0 | 0.00 0.00 0.00 0.00  (= model)
  neither  arm_cmo_no2                  seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00  (= arm_cmo)
  neither  arm_cmo_tern_only_no2        seeds s1       eNPV     -9.4 | 0.00 0.00 0.00 0.00  (= arm_cmo_tern_only)

weighted by the odds (arithmetic over the paths above; not an ensemble):
  model                        eNPV    847.0  vs model:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo                      eNPV    847.9  vs model:    +0.9 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only            eNPV    847.2  vs model:    +0.3 | P(K short launch) 0.00  P(K short ever) 0.09  P(T short launch) 0.00  P(T files late) 0.09
  arm_build_now                eNPV    741.1  vs model:  -105.8 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_build_after              eNPV    824.0  vs model:   -23.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  model_literal                eNPV    839.1  vs model_literal:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo_literal              eNPV    847.9  vs model_literal:    +8.8 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only_literal    eNPV    847.7  vs model_literal:    +8.6 | P(K short launch) 0.00  P(K short ever) 0.09  P(T short launch) 0.00  P(T files late) 0.09
  model_no2                    eNPV    768.2  vs model_no2:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.18  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo_no2                  eNPV    784.6  vs model_no2:   +16.3 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only_no2        eNPV    695.7  vs model_no2:   -72.5 | P(K short launch) 0.00  P(K short ever) 0.27  P(T short launch) 0.00  P(T files late) 0.09
exit 0
```

```
python3 agent/studies/launch_supply/analyze.py --peak     (same paths revalued at peak-year value; the dynamics do not depend on value)
weighted by the odds (arithmetic over the paths above; not an ensemble):
  model                        eNPV   1896.6  vs model:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo                      eNPV   1917.0  vs model:   +20.4 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only            eNPV   1915.9  vs model:   +19.3 | P(K short launch) 0.00  P(K short ever) 0.09  P(T short launch) 0.00  P(T files late) 0.09
  arm_build_now                eNPV   1790.7  vs model:  -105.9 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_build_after              eNPV   1873.4  vs model:   -23.2 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  model_literal                eNPV   1876.0  vs model_literal:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo_literal              eNPV   1917.1  vs model_literal:   +41.1 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only_literal    eNPV   1916.8  vs model_literal:   +40.8 | P(K short launch) 0.00  P(K short ever) 0.09  P(T short launch) 0.00  P(T files late) 0.09
  model_no2                    eNPV   1677.9  vs model_no2:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.18  P(T short launch) 0.18  P(T files late) 0.18
  arm_cmo_no2                  eNPV   1739.0  vs model_no2:   +61.1 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.09
  arm_cmo_tern_only_no2        eNPV   1591.7  vs model_no2:   -86.2 | P(K short launch) 0.00  P(K short ever) 0.27  P(T short launch) 0.00  P(T files late) 0.09
exit 0
```

```
python3 agent/studies/launch_supply/analyze.py --pk=0.5   (same paths, worlds reweighted for Kestrel at 50%)
weighted by the odds (arithmetic over the paths above; not an ensemble):
  model                        eNPV    738.7  vs model:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.15  P(T files late) 0.15
  arm_cmo                      eNPV    737.2  vs model:    -1.5 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.07
  arm_cmo_tern_only            eNPV    736.7  vs model:    -2.0 | P(K short launch) 0.00  P(K short ever) 0.07  P(T short launch) 0.00  P(T files late) 0.07
  arm_build_now                eNPV    630.9  vs model:  -107.9 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.15  P(T files late) 0.15
  arm_build_after              eNPV    719.6  vs model:   -19.1 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.15  P(T files late) 0.15
  model_literal                eNPV    732.2  vs model_literal:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.15  P(T files late) 0.15
  arm_cmo_literal              eNPV    737.3  vs model_literal:    +5.1 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.07
  arm_cmo_tern_only_literal    eNPV    737.1  vs model_literal:    +4.9 | P(K short launch) 0.00  P(K short ever) 0.07  P(T short launch) 0.00  P(T files late) 0.07
  model_no2                    eNPV    673.1  vs model_no2:    +0.0 | P(K short launch) 0.00  P(K short ever) 0.15  P(T short launch) 0.15  P(T files late) 0.15
  arm_cmo_no2                  eNPV    684.5  vs model_no2:   +11.3 | P(K short launch) 0.00  P(K short ever) 0.00  P(T short launch) 0.00  P(T files late) 0.07
  arm_cmo_tern_only_no2        eNPV    610.4  vs model_no2:   -62.7 | P(K short launch) 0.00  P(K short ever) 0.23  P(T short launch) 0.00  P(T files late) 0.07
exit 0
```

Figure, from the weighted table above:

```
julia --project=. agent/scripts/chartkit.jl agent/studies/launch_supply/figures/options_vs_wait.json --out=agent/studies/launch_supply/figures/options_vs_wait.svg
chartkit: wrote agent/studies/launch_supply/figures/options_vs_wait.svg — bars_error, 6.2 KiB, self-contained
exit 0
```

Path-by-path check of the weighted deltas (same run files; `paired.py` imports `analyze.runs`/`analyze.W`). This is the variation across the 3 both-world paths; it is not a standard error:

```
python3 agent/studies/launch_supply/paired.py   (run from the study directory)
both-world per-seed difference vs its wait reference, and the odds-weighted delta if that path were the both-world path (other worlds: seed 1)
  arm_cmo                  vs model      both-world d s1..s3   +51.9   -14.3   +65.2 | weighted per path   +5.7  -12.2   +9.3
  arm_cmo_tern_only        vs model      both-world d s1..s3   +52.2   -14.1   +60.2 | weighted per path   +5.5  -12.4   +7.7
  arm_build_now            vs model      both-world d s1..s3   -87.7   -86.3   -86.6 | weighted per path -106.0 -105.7 -105.8
  arm_build_after          vs model      both-world d s1..s3   -84.3   -85.7   -85.1 | weighted per path  -22.8  -23.1  -23.0
  arm_cmo_no2              vs model_no2  both-world d s1..s3   +76.3   +96.8  +100.7 | weighted per path  +12.3  +17.8  +18.9
  arm_cmo_tern_only_no2    vs model_no2  both-world d s1..s3  -258.7  -214.7  -236.9 | weighted per path  -78.4  -66.6  -72.5
walk5: PV of $5M at t=14 4.522; arm_cmo +0.99 (P doomed = neither 0.22); arm_cmo_tern_only +2.49 (P doomed = Tern fails 0.55)
exit 0
```

```
julia --project=. agent/scripts/rdx.jl draw agent/studies/launch_supply/scenarios/model__both.rdj.json --out=agent/studies/launch_supply/figures/exec_map_wait_both_s1.svg --exec-map --seed=1
"ok": true, "layer": "exec_map", "artifact": ".../figures/exec_map_wait_both_s1.svg"
exit 0
```

### turn 7

The user asked three things:

1. Put a range on the probabilities by running more both-succeed paths, plus the 12 single-asset runs offered last turn: "I'd rather wait a day than present two runs out of three."
2. Say what "Tern short" costs in batches, $M and duration, and explain the −$12.2M path.
3. Check this framing: "Sign the CMO reservation before 30 November, usable for either product. Expected value is neutral, around break-even. Worst case, if both trials fail, we lose about $9–10M. In return, if both succeed, Tern launches on time and in full instead of short and two months late. Building a second line loses money in every case we ran."

He will have legal confirm reassignment and will ask the CMO for indicative terms on a second slot. He agrees to keep Tern-only off the recommendation unless the second contract is in writing.

Compute relay, not from the user. The 96-seed wait ensemble summary was pasted into the conversation, and relayed compute output is usable, so it is recorded here verbatim:

```
rdx ensemble agent/studies/launch_supply/model.rdj.json --seeds=96 --root-seed=2609 --per-seed --annual=0.09 --periods-per-year=12 (exit 0, 3171 s)
| metric | mean | sem | median | q25 | q75 | n |
| enpv | 765.9561 | 61.6450 | 1072.5145 | 0.0000 | 1135.8647 | 96 |
| final:k_success | 0.5417 | 0.0511 | 1.0000 | 0.0000 | 1.0000 | 96 |
| final:t_success | 0.3542 | 0.0491 | 0.0000 | 0.0000 | 1.0000 | 96 |
| final:k_short | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 96 |
| final:k_short_any | 0.0208 | 0.0147 | 0.0000 | 0.0000 | 0.0000 | 96 |
| final:t_short | 0.1458 | 0.0362 | 0.0000 | 0.0000 | 0.0000 | 96 |
| final:t_file_delay | 0.0417 | 0.0293 | 0.0000 | 0.0000 | 0.0000 | 96 |
| final:k_lost | 0.1771 | 0.1574 | 0.0000 | 0.0000 | 0.0000 | 96 |
| final:t_lost | 24.0417 | 4.9599 | 0.0000 | 0.0000 | 0.0000 | 96 |
```

The compute rule still stands: no ensemble or effect calls, and about 100 single-seed runs this turn. The per-seed samples are available on request.

This turn made 90 runs, all with `run_scenarios.sh` and all `exit 0`, taking about 56 s each:

- wait (`model`) and reserve-now (`arm_cmo`) with both readouts succeeding, seeds 4–37 (68 runs);
- the same two options in the Kestrel-only, Tern-only and neither worlds, seeds 2–3 (12 runs);
- `arm_build_now` and `arm_build_after` with both succeeding, seeds 4–8 (10 runs).

```
./run_scenarios.sh model konly 2 arm_cmo konly 2 model tonly 2 arm_cmo tonly 2 model neither 2 arm_cmo neither 2 model konly 3 arm_cmo konly 3 model tonly 3 arm_cmo tonly 3   -> 10 x exit 0 (9:18.92 total)
./run_scenarios.sh model neither 3 arm_cmo neither 3 model both 4 arm_cmo both 4 ... model both 7                -> 9 x exit 0
./run_scenarios.sh arm_cmo both 7 {model,arm_cmo} both 8..11                                                       -> 9 x exit 0
./run_scenarios.sh {model,arm_cmo} both 12..15 model both 16                                                        -> 9 x exit 0 (see note)
./run_scenarios.sh arm_cmo both 16 {model,arm_cmo} both 17..20                                                      -> 9 x exit 0 (see note)
./run_scenarios.sh {model,arm_cmo} both 21..24 model both 25                                                        -> 9 x exit 0
./run_scenarios.sh arm_cmo both 25 {model,arm_cmo} both 26..29                                                      -> 9 x exit 0
./run_scenarios.sh {model,arm_cmo} both 30..32 model both 33                                                        -> 7 x exit 0 (4:31.79)
./run_scenarios.sh arm_cmo both 33 {model,arm_cmo} both 34..36                                                      -> 7 x exit 0 (4:50.13)
./run_scenarios.sh model both 37 arm_cmo both 37 {arm_build_now,arm_build_after} both 4..6                          -> 8 x exit 0 (5:06.15)
./run_scenarios.sh {arm_build_now,arm_build_after} both 7..8                                                        -> 4 x exit 0 (2:32.29)
ls runs | wc -l   -> 141
```

**Protocol note.** The seed 12–16 and seed 16–20 batches exceeded the 600 s tool timeout. The harness moved them to the background, and the next two batches started while they were still running, so up to three runs overlapped. I did not choose to background anything. Afterwards, batches were cut to 8 or fewer runs. Both background tasks finished with every line `exit 0`, and each wrote different files. A run is determined by (document, seed), so the overlap cannot change any result.

`paired_both.py` computes the paired comparison: the same seed gives the same readout dates, slips and lead times in both options. Its standard errors are across paths.

```
python3 agent/studies/launch_supply/paired_both.py
BOTH SUCCEED, wait vs reserve-now, 37 paired paths (s1..s37)
 seed  d_eNPV |  dK_rev  dT_rev  d_cost | wait: Tlaunch delay Tshort Tlost(batches) TlostPV months_short | cmo: Tlaunch delay Tshort Tlost
 s1     +51.9 |    +0.0   +65.0   -13.2 | 37 2 1   6.2   29.2 38-42 (5)      | 35 0 0 0.0
 s2     -14.3 |    +0.0    -3.6   -10.8 | 30 0 0  11.0   51.9 38-44 (7)      | 33 3 0 0.0
 s3     +65.2 |    +0.0   +79.9   -14.7 | 35 2 1   9.2   43.6 36-42 (7)      | 33 0 0 0.0
 s4     +47.1 |    +0.0   +61.2   -14.1 | 31 0 1  12.9   61.2 36-43 (8)      | 31 0 0 0.0
 s5     +57.9 |    +0.0   +71.5   -13.6 | 36 2 1   7.5   35.5 37-42 (6)      | 34 0 0 0.0
 s6     +27.9 |    +0.0   +40.9   -13.0 | 29 0 0  12.5   59.6 36-42 (7)      | 30 1 0 0.0
 s7     +51.8 |    +0.0   +65.0   -13.3 | 37 2 1   6.2   29.2 38-42 (5)      | 35 0 0 0.0
 s8     +42.2 |    +0.0   +56.0   -13.9 | 30 0 1  15.7   74.6 35-44 (10)     | 31 1 0 0.0
 s9     +20.9 |    +0.0   +33.0   -12.1 | 31 0 1  14.7   69.8 35-44 (10)     | 33 2 0 0.0
 s10    +91.0 |    +0.0  +107.0   -16.0 | 35 3 1  11.1   52.4 36-43 (8)      | 32 0 0 0.0
 s11    +51.0 |    +0.0   +63.6   -12.5 | 37 2 1   5.9   27.8 38-42 (5)      | 35 0 0 0.0
 s12    +27.8 |    +0.0   +40.8   -13.0 | 31 0 0   8.6   40.8 37-41 (5)      | 31 0 0 0.0
 s13    +85.6 |    +3.4   +97.9   -15.7 | 35 2 1  13.1   61.6 36-44 (9)      | 33 0 0 0.0
 s14     +9.6 |    +0.0   +21.4   -11.8 | 31 0 0   8.4   39.9 37-41 (5)      | 32 1 0 0.0
 s15     +4.4 |    +0.0   +16.3   -11.9 | 29 0 0  15.2   72.1 36-43 (8)      | 32 3 0 0.0
 s16    +51.6 |   +25.0   +40.4   -13.8 | 30 0 1  20.3   95.9 35-44 (10)     | 33 3 0 0.0
 s17    +68.8 |    +0.0   +83.0   -14.2 | 38 3 1   6.3   29.5 39-43 (5)      | 35 0 0 0.0
 s18    +28.1 |    +0.0   +40.9   -12.8 | 29 0 0  12.5   59.6 36-41 (6)      | 30 1 0 0.0
 s19    +50.4 |    +0.0   +63.6   -13.2 | 37 2 1   5.9   27.8 38-42 (5)      | 35 0 0 0.0
 s20    +14.6 |    +0.0   +26.2   -11.7 | 31 0 0  13.4   63.1 37-44 (8)      | 33 2 0 0.0
 s21    +50.5 |    +0.0   +64.0   -13.5 | 36 2 1   5.9   28.0 37-41 (5)      | 34 0 0 0.0
 s22    +35.7 |    +0.0   +48.3   -12.6 | 33 0 1  10.2   48.3 37-42 (6)      | 33 0 0 0.0
 s23    +72.8 |   +11.5   +76.8   -15.5 | 33 0 1  16.3   76.8 35-44 (10)     | 33 0 0 0.0
 s24     +7.8 |    +0.0   +20.0   -12.2 | 30 0 0  16.0   75.5 36-44 (9)      | 33 3 0 0.0
 s25    +30.8 |    +0.0   +44.4   -13.6 | 30 0 0  13.3   63.0 36-43 (8)      | 31 1 0 0.0
 s26    +40.2 |    +0.0   +52.7   -12.4 | 33 0 1  11.1   52.7 36-42 (7)      | 33 0 0 0.0
 s27     -3.7 |    +0.0    +6.6   -10.3 | 30 0 1  13.0   62.0 35-44 (9)      | 33 3 0 0.0
 s28    +51.8 |    +0.0   +64.5   -12.8 | 37 2 1   6.1   28.8 38-42 (5)      | 35 0 0 0.0
 s29    +48.9 |    +0.0   +61.8   -13.0 | 32 0 1  13.0   61.8 35-42 (8)      | 32 0 0 0.0
 s30    +57.8 |    +0.0   +71.5   -13.7 | 36 2 1   7.5   35.5 37-42 (6)      | 34 0 0 0.0
 s31    +41.1 |    +0.0   +53.2   -12.1 | 33 0 1  11.2   53.2 36-42 (7)      | 33 0 0 0.0
 s32     -5.5 |    +0.0    +4.9   -10.4 | 31 0 0   8.8   41.8 37-44 (6)      | 33 2 0 0.0
 s33    +27.4 |    +4.0   +36.3   -13.0 | 31 0 0  15.6   73.1 37-44 (8)      | 33 2 0 0.0
 s34    +26.7 |    +0.0   +38.5   -11.8 | 33 0 1   8.1   38.5 37-41 (5)      | 33 0 0 0.0
 s35    +26.6 |    +0.0   +38.5   -11.9 | 33 0 1   8.1   38.5 37-41 (5)      | 33 0 0 0.0
 s36    +65.5 |    +0.0   +79.8   -14.3 | 36 2 1   9.3   43.8 37-43 (7)      | 34 0 0 0.0
 s37    +24.6 |    +0.0   +35.6   -11.1 | 34 0 1   7.5   35.6 37-41 (5)      | 34 0 0 0.0

 d_eNPV (reserve - wait | both): mean +38.7 se 4.1 n 37; min -14.3 max +91.0; negative on 3 of 37
 wait     Tern short at launch 26/37 = 0.70 (95% Wilson 0.54-0.83) | files late 13/37 = 0.35 (0.22-0.51), months [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0] | any Tern demand lost 37/37; lost batches mean 10.7 se 0.6; lost value PV mean 50.9 se 2.9 $M | Kestrel short ever 4/37, at launch 0/37
 reserve  Tern short at launch 0/37 = 0.00 (95% Wilson 0.00-0.09) | files late 14/37 = 0.38 (0.24-0.54), months [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0] | any Tern demand lost 0/37; lost batches mean 0.0 se 0.0; lost value PV mean 0.0 se 0.0 $M | Kestrel short ever 0/37, at launch 0/37

 reserve-now late filings (seed: months late, Tern launch reserve vs wait): [('s2', 3.0, 33, 30), ('s6', 1.0, 30, 29), ('s8', 1.0, 31, 30), ('s9', 2.0, 33, 31), ('s14', 1.0, 32, 31), ('s15', 3.0, 32, 29), ('s16', 3.0, 33, 30), ('s18', 1.0, 30, 29), ('s20', 2.0, 33, 31), ('s24', 3.0, 33, 30), ('s25', 1.0, 31, 30), ('s27', 3.0, 33, 30), ('s32', 2.0, 33, 31), ('s33', 2.0, 33, 31)]

SINGLE-ASSET AND NEITHER WORLDS, paired, seeds available in both options
 konly    s1: d_eNPV   -8.9 = dK   +0.0 dT   +0.0 dCost   -8.9 | Tlaunch wait None reserve None, delay wait 0 reserve 0
 konly    s2: d_eNPV   -8.8 = dK   +0.0 dT   +0.0 dCost   -8.8 | Tlaunch wait None reserve None, delay wait 0 reserve 0
 konly    s3: d_eNPV   -8.6 = dK   +0.0 dT   +0.0 dCost   -8.6 | Tlaunch wait None reserve None, delay wait 0 reserve 0
 tonly    s1: d_eNPV  -18.6 = dK   +0.0 dT   +0.0 dCost  -18.6 | Tlaunch wait 35 reserve 35, delay wait 0 reserve 0
 tonly    s2: d_eNPV  -72.1 = dK   +0.0 dT  -55.4 dCost  -16.7 | Tlaunch wait 30 reserve 33, delay wait 0 reserve 3
 tonly    s3: d_eNPV  -19.5 = dK   +0.0 dT   +0.0 dCost  -19.5 | Tlaunch wait 33 reserve 33, delay wait 0 reserve 0
 neither  s1: d_eNPV   -9.4 = dK   +0.0 dT   +0.0 dCost   -9.4 | Tlaunch wait None reserve None, delay wait 0 reserve 0
 neither  s2: d_eNPV   -9.4 = dK   +0.0 dT   +0.0 dCost   -9.4 | Tlaunch wait None reserve None, delay wait 0 reserve 0
 neither  s3: d_eNPV   -9.4 = dK   +0.0 dT   +0.0 dCost   -9.4 | Tlaunch wait None reserve None, delay wait 0 reserve 0

WEIGHTED reserve-now vs wait (odds weights {'both': 0.27, 'konly': 0.33, 'tonly': 0.18, 'neither': 0.22} )
 both +38.7 (se 4.1, n 37) | konly -8.8 (se 0.07, n 3) | tonly -36.7 (se 17.68, n 3) | neither -9.4 (se 0.00, n 3)
 weighted delta -1.1, se 3.4; delta - 2se -7.8, delta + 2se +5.6
 P(Tern short at launch) wait = 0.27 x 26/37 = 0.190 (0.27 x Wilson: 0.146-0.223); reserve = 0.27 x 0/37
 arm_build_now: both-world d vs wait mean -86.9 se 3.1 n 8 range -104.7..-72.8; Tern short 6/8
 arm_build_after: both-world d vs wait mean -84.3 se 0.6 n 8 range -86.8..-82.3; Tern short 6/8

 wait eNPV | both: mean 1511.0 se 7.6 n 37
 wait eNPV | konly: mean 1096.1 n 3
 wait eNPV | tonly: mean 479.6 n 3
 wait eNPV | neither: mean 0.0 n 3

WEIGHTED P(Tern files late) and mean months late, by world (paths):
 model    P(late) 0.095  mean months late 0.204  [both 13/37, konly 0/3, tonly 0/3, neither 0/3]
 arm_cmo  P(late) 0.162  mean months late 0.384  [both 14/37, konly 0/3, tonly 1/3, neither 0/3]

ENSEMBLE CROSS-CHECK (wait): forced-world eNPV at realized rates 751.2 vs ensemble 765.9561 (sem 61.6450); at design odds 856.0
  t_short: 0.1919 x 0.703 = 0.135 vs ensemble 0.1458 (sem 0.0362)
  t_file_delay months: 0.1919 x 0.757 = 0.145 vs ensemble 0.0417 (sem 0.0293)
  k_short_any: 0.1919 x 0.108 = 0.021 vs ensemble 0.0208 (sem 0.0147); t_lost tenths: 0.1919 x 107.5 = 20.6 vs ensemble 24.0417 (sem 4.9599)
exit 0
```

Kestrel's lost demand under wait on the four both-succeed paths where it ran short:

```
python3 -c "... from paired_both import load, seeds ..."
s13 wait Kestrel lost 0.5 batches; first-6-month short 0.0
s16 wait Kestrel lost 3.7 batches; first-6-month short 0.0
s23 wait Kestrel lost 1.7 batches; first-6-month short 0.0
s33 wait Kestrel lost 0.6 batches; first-6-month short 0.0
exit 0
```

The same paired paths, revalued at the peak value per batch and reweighted for Kestrel PoS 0.50. `paired_both.py` now passes its flags through to `analyze`. After the label fix, the P(Tern short) line prints the world weight rather than a hard-coded 0.27; the base-case numbers are unchanged.

```
python3 agent/studies/launch_supply/paired_both.py --peak   (excerpt)
 d_eNPV (reserve - wait | both): mean +117.3 se 10.5 n 37; min -19.8 max +255.0; negative on 1 of 37
 both +117.3 (se 10.5, n 37) | konly -8.8 (se 0.07, n 3) | tonly -65.0 (se 45.98, n 3) | neither -9.4 (se 0.00, n 3)
 weighted delta +15.0, se 8.8; delta - 2se -2.5, delta + 2se +32.5
 P(Tern short at launch) wait = 0.270 x 26/37 = 0.190 (0.270 x Wilson: 0.146-0.223); reserve = 0.270 x 0/37
 arm_build_now: both-world d vs wait mean -88.7 se 7.3 n 8 range -132.2..-55.7; Tern short 6/8
 arm_build_after: both-world d vs wait mean -85.0 se 0.9 n 8 range -89.1..-81.6; Tern short 6/8
 model    P(late) 0.095  mean months late 0.204  [both 13/37, konly 0/3, tonly 0/3, neither 0/3]
 arm_cmo  P(late) 0.162  mean months late 0.384  [both 14/37, konly 0/3, tonly 1/3, neither 0/3]
exit 0
```

```
python3 agent/studies/launch_supply/paired_both.py --pk=0.5   (excerpt)
 d_eNPV (reserve - wait | both): mean +38.7 se 4.1 n 37; min -14.3 max +91.0; negative on 3 of 37
 both +38.7 (se 4.1, n 37) | konly -8.8 (se 0.07, n 3) | tonly -36.7 (se 17.68, n 3) | neither -9.4 (se 0.00, n 3)
 weighted delta -4.5, se 4.1; delta - 2se -12.7, delta + 2se +3.6
 P(Tern short at launch) wait = 0.225 x 26/37 = 0.158 (0.225 x Wilson: 0.122-0.186); reserve = 0.225 x 0/37
 arm_build_now: both-world d vs wait mean -86.9 se 3.1 n 8 range -104.7..-72.8; Tern short 6/8
 arm_build_after: both-world d vs wait mean -84.3 se 0.6 n 8 range -86.8..-82.3; Tern short 6/8
 model    P(late) 0.079  mean months late 0.170  [both 13/37, konly 0/3, tonly 0/3, neither 0/3]
 arm_cmo  P(late) 0.160  mean months late 0.395  [both 14/37, konly 0/3, tonly 1/3, neither 0/3]
exit 0
```

### turn 8

User request (summary, not numbers): Kestrel reads out in March 2027, nine months before Tern. Model the rule "reserve the CMO now; if Kestrel fails, move Tern's PPQ in-house, file Tern on our own line, give the CMO notice and pay the forfeit; if Kestrel succeeds, Tern stays at the CMO". His guess: this cuts the Tern-only loss to roughly the both-fail cost and removes the $72M tail. Backup qualification (from regulatory, via the user): both sites can be listed if both have PPQ data, about $0.5M/batch all-in, about $1.5M for 3; only relevant when Kestrel has failed; not to be modelled separately unless cheap. Budget: spend it on the Tern-only world under the exit rule; the no-second-contract sensitivity waits. He will not quote late-filing numbers until the mismatch with the 96-seed baseline is resolved.

Relayed by the coordinator (from the 96-member wait ensemble, root seed 2609): 20 both-succeed members; among them t_file_delay > 0 on 2 of 20 (2 months each) and t_short > 0 on 14 of 20; members 20 / 96 / 6 had enpv 1526.7 / 1444.8 / 1529.5, t_file_delay 2 / 2 / 0, t_lost 92 / 59 / 91. Compute rules: no `ensemble`/`effect`, single-seed runs in batches of 8 or fewer, about 80 runs.

Protocol note: `ens_check.py` was rewritten this turn as a self-contained script (an earlier, uncommitted draft of the same check was overwritten); its output below is the one of record.

Late-filing mismatch. Ensemble member k is seeded `hash((2609, k))`; members 20, 96 and 6 were replayed on the unforced model and on the forced both-succeed copy (6 runs, one at a time, foreground):

```
julia --project=. agent/scripts/rdx.jl run agent/studies/launch_supply/model.rdj.json --seed=3689622802375696111 --trace=... > runs_ens/model__natural__m20.json      exit 0
julia --project=. agent/scripts/rdx.jl run agent/studies/launch_supply/scenarios/model__both.rdj.json --seed=3689622802375696111 --trace=... > runs_ens/model__both__m20.json   exit 0
(same for member 96, seed 8272771591343840987, and member 6, seed 4787811836286466114; all exit 0)
```

```
python3 agent/studies/launch_supply/ens_check.py
member 20: forced == unforced True; our run delay 2 t_lost 92 (ensemble 2 / 92); our eNPV 1515.8 vs ensemble 1526.7 -> ratio 1.00721 (1.09^(1/12) = 1.00721)
member 96: forced == unforced True; our run delay 2 t_lost 59 (ensemble 2 / 59); our eNPV 1434.4 vs ensemble 1444.8 -> ratio 1.00724 (1.09^(1/12) = 1.00721)
member 6: forced == unforced True; our run delay 0 t_lost 91 (ensemble 0 / 91); our eNPV 1518.5 vs ensemble 1529.5 -> ratio 1.00723 (1.09^(1/12) = 1.00721)
wait, both succeed: late filing ours 13/37, ensemble 2/20; Fisher exact two-sided p = 0.059
pooled 15/57 = 0.263 (95% Wilson 0.166-0.390)
Tern short at launch, wait, both succeed: ours 26/37, ensemble 14/20, pooled 40/57 = 0.702 (95% Wilson 0.573-0.805); x 0.27 = 0.189 (0.155-0.217)
reserve now, both succeed: late filing 14/37 = 0.378 (95% Wilson 0.241-0.539); no ensemble exists for this option
exit 0
```

Reading: forcing reproduces the ensemble member exactly (same late-filing months, same lost Tern demand); the eNPV differs by the constant factor 1.09^(1/12), a one-month discounting-offset convention that cancels in every comparison. The 2/20 vs 13/37 late-filing gap is at the edge of chance (Fisher two-sided p = 0.059), not a model error; pooled 15/57.

Exit rule. `gen_models.py` edits: params `exit_at_k` and `cmo_exited` (both 0 by default); the year-2 CMO fee is charged as the forfeit when `cmo_exited`; new rule `cmo_exit_k` (fires once at Kestrel readout if `cmo_now` and `exit_at_k` and Kestrel failed: CMO capacity, reservation and CMO PPQ set to 0, `cmo_now` 0, `cmo_exited` 1, which returns the arm to the wait posture: Tern PPQ in-house, Tern filed on own-line PPQ); new arm `arm_cmo_exit.rdj.json` (`cmo_now` 1, `exit_at_k` 1).

```
python3 agent/studies/launch_supply/gen_models.py
wrote 23 documents + 48 scenario copies
exit 0
julia --project=. agent/scripts/rdx.jl validate agent/studies/launch_supply/arm_cmo_exit.rdj.json
ok, 0 findings
exit 0
```

Regeneration check (the edit must not move any stored run): `arm_cmo` both s1 and `model` tonly s2 were re-run from the regenerated documents and compared byte for byte against the stored `runs/arm_cmo__both__s1.json` and `runs/model__tonly__s2.json`: identical, both exit 0.

Runs this turn via `run_scenarios.sh` (foreground, batches of 8 or fewer, every line `exit 0`): `arm_cmo_exit` both s1, konly s1, neither s1, tonly s1-30; `model` tonly s4-30; `arm_cmo` tonly s4-10. With the 6 ensemble replays and the 2 regeneration checks, 75 runs this turn. All 79 turn-relevant run files parse.

```
python3 agent/studies/launch_supply/exit_review.py [ | --peak | --pk=0.5]
$ exit_review.py 
inertness check, arm_cmo_exit vs arm_cmo on the same seed (must be 0): [('both', 's1', 0.0), ('konly', 's1', 0.0), ('neither', 's1', 0.0)]

TERN-ONLY WORLD (Kestrel fails, Tern succeeds), paired against wait
 arm_cmo_exit  vs wait: mean -9.40 se 0.02 n 30; min -9.66 max -9.25; paths below -20: 0/30 | Tern files late 0/30 (wait 0/30); Tern demand lost on 0/30
 arm_cmo       vs wait: mean -31.67 se 5.95 n 10; min -72.09 max -18.32; paths below -20: 4/10 | Tern files late 4/10 (wait 0/10); Tern demand lost on 0/10
  arm_cmo_exit decomposition: cost diff mean -9.40 (min -9.66, max -9.25); Tern revenue diff mean +0.00 (min +0.00 max +0.00)

WEIGHTED 'reserve now, review at Kestrel readout' vs wait, odds {'both': 0.27, 'konly': 0.33, 'tonly': 0.18, 'neither': 0.22}
  both +38.71 (se 4.06, n 37, from arm_cmo) | konly -8.76 (se 0.07, n 3, from arm_cmo) | tonly -9.40 (se 0.02, n 30, from arm_cmo_exit) | neither -9.38 (se 0.00, n 3, from arm_cmo)
  weighted delta +3.81, se 1.10; -2se +1.62, +2se +6.00
  loss given each non-both world (mean): {'konly': -8.76, 'tonly': -9.4, 'neither': -9.38}
  P(Tern short at launch): wait 0.270 x 26/37 = 0.190 (Wilson 0.146-0.223); reserve 0.270 x 0/37; Tern-only world short: wait 0/30, review 0/30
exit 0
$ exit_review.py --peak
inertness check, arm_cmo_exit vs arm_cmo on the same seed (must be 0): [('both', 's1', 0.0), ('konly', 's1', 0.0), ('neither', 's1', 0.0)]

TERN-ONLY WORLD (Kestrel fails, Tern succeeds), paired against wait
 arm_cmo_exit  vs wait: mean -9.40 se 0.02 n 30; min -9.66 max -9.25; paths below -20: 0/30 | Tern files late 0/30 (wait 0/30); Tern demand lost on 0/30
 arm_cmo       vs wait: mean -51.52 se 15.43 n 10; min -156.99 max -18.32; paths below -20: 4/10 | Tern files late 4/10 (wait 0/10); Tern demand lost on 0/10
  arm_cmo_exit decomposition: cost diff mean -9.40 (min -9.66, max -9.25); Tern revenue diff mean +0.00 (min +0.00 max +0.00)

WEIGHTED 'reserve now, review at Kestrel readout' vs wait, odds {'both': 0.27, 'konly': 0.33, 'tonly': 0.18, 'neither': 0.22}
  both +117.25 (se 10.51, n 37, from arm_cmo) | konly -8.76 (se 0.07, n 3, from arm_cmo) | tonly -9.40 (se 0.02, n 30, from arm_cmo_exit) | neither -9.38 (se 0.00, n 3, from arm_cmo)
  weighted delta +25.01, se 2.84; -2se +19.34, +2se +30.69
  loss given each non-both world (mean): {'konly': -8.76, 'tonly': -9.4, 'neither': -9.38}
  P(Tern short at launch): wait 0.270 x 26/37 = 0.190 (Wilson 0.146-0.223); reserve 0.270 x 0/37; Tern-only world short: wait 0/30, review 0/30
exit 0
$ exit_review.py --pk=0.5
inertness check, arm_cmo_exit vs arm_cmo on the same seed (must be 0): [('both', 's1', 0.0), ('konly', 's1', 0.0), ('neither', 's1', 0.0)]

TERN-ONLY WORLD (Kestrel fails, Tern succeeds), paired against wait
 arm_cmo_exit  vs wait: mean -9.40 se 0.02 n 30; min -9.66 max -9.25; paths below -20: 0/30 | Tern files late 0/30 (wait 0/30); Tern demand lost on 0/30
 arm_cmo       vs wait: mean -31.67 se 5.95 n 10; min -72.09 max -18.32; paths below -20: 4/10 | Tern files late 4/10 (wait 0/10); Tern demand lost on 0/10
  arm_cmo_exit decomposition: cost diff mean -9.40 (min -9.66, max -9.25); Tern revenue diff mean +0.00 (min +0.00 max +0.00)

WEIGHTED 'reserve now, review at Kestrel readout' vs wait, odds {'both': 0.225, 'konly': 0.275, 'tonly': 0.225, 'neither': 0.275}
  both +38.71 (se 4.06, n 37, from arm_cmo) | konly -8.76 (se 0.07, n 3, from arm_cmo) | tonly -9.40 (se 0.02, n 30, from arm_cmo_exit) | neither -9.38 (se 0.00, n 3, from arm_cmo)
  weighted delta +1.61, se 0.91; -2se -0.22, +2se +3.44
  loss given each non-both world (mean): {'konly': -8.76, 'tonly': -9.4, 'neither': -9.38}
  P(Tern short at launch): wait 0.225 x 26/37 = 0.158 (Wilson 0.122-0.186); reserve 0.225 x 0/37; Tern-only world short: wait 0/30, review 0/30
exit 0
```

Diagnose for the new arm, with the flags used for every other arm (a first unflagged `rdx diagnose arm_cmo_exit.rdj.json` returned exit 1, 4 pass · 4 warn · 1 fail, the extra warns being the missing `--expect-binds`/`--expect-reusable` and a `bounded_pools` warn on the default tally; superseded by the flagged run below). `arm_cmo` was re-diagnosed alongside because DOCTOR's table still carried the turn-4 cost:

```
rdx diagnose <arm> --seed=1 --expect-binds=own_free --expect-reusable=line,cmo,line2 --tally=k_sold,t_sold,k_lost,t_lost,own_used,cmo_used,l2_used,mkt_done,k_age,t_age,t_file_delay,k_orders,t_orders
arm_cmo_exit exit 1
7 1 1 ['transitions_fire', 'horizon_sufficient'] {'reward': 1703.7, 'cost': 98.07, 'net': 1605.63}
arm_cmo exit 1
7 1 1 ['transitions_fire', 'horizon_sufficient'] {'reward': 1703.7, 'cost': 98.07, 'net': 1605.63}
```

Runs this turn, final count: 75 simulation runs plus 3 diagnose runs = 78.

Overall chance of a late Tern filing under the review rule (konly and neither worlds are 0/3 under both options, and tonly is 0/30 under wait):

```
python3 -c "... from paired_both import load, seeds, W ..."
review: both late 14 / 37 tonly late 0 / 30 -> overall P(late) 0.102
wait: both late 13 / 37 overall 0.095
exit 0
```

### turn 9

User request (final turn, summary): three edits to the ExCo paragraph. (1) Add that at Kestrel PoS 50% the result is roughly break-even, around +$1.6M. (2) Replace "about $84–87M worse" for building with the range, which he read as −$73M to −$105M across 8 paths. (3) Say "$9–11M" for the failure cost, to cover the roughly $1M of understated in-house validation cost. On his side: legal confirms the slot can be reassigned; regulatory confirms Tern can be filed on our line if Kestrel fails. No new simulation runs this turn (coordinator).

Check of each edit against stored runs (no simulation):

```
python3 agent/studies/launch_supply/paired_both.py   (excerpt)
 both +38.7 (se 4.1, n 37) | konly -8.8 (se 0.07, n 3) | tonly -31.7 (se 5.95, n 10) | neither -9.4 (se 0.00, n 3)
 weighted delta -0.2, se 1.5; delta - 2se -3.3, delta + 2se +2.9
 arm_build_now: both-world d vs wait mean -86.9 se 3.1 n 8 range -104.7..-72.8; Tern short 6/8
 arm_build_after: both-world d vs wait mean -84.3 se 0.6 n 8 range -86.8..-82.3; Tern short 6/8
exit 0
```

- Edit 1 holds: `exit_review.py --pk=0.5` (turn 8) gives weighted delta +1.61, se 0.91, -2se -0.22, +2se +3.44.
- Edit 2 is right for build-now only: build now -104.7 to -72.8, build after -86.8 to -82.3 (8 paths each, both succeed). Both ranges are quoted.
- Edit 3 is accepted as a conservative bound, with a correction to my turn-8 caveat. The measured non-both costs are 8.76 / 9.40 (range 9.25–9.66) / 9.38. The Tern-only decomposition (turn 8) is all CMO fee plus forfeit: cost difference -9.40, Tern revenue difference +0.00. So own-line PPQ is charged identically under wait and under review. The user's $0.5M/batch was for a backup-site qualification; under the review rule the own line is Tern's only filed site, as it is under wait, so a higher own-line PPQ cost would fall on both options and largely cancel. My "understated by about $1M" overstated the gap in the comparison.
- The no-review row in FINDINGS is updated from turn 7's -1.1 ± 3.4 (Tern only on 3 paths) to paired_both's current -0.2 ± 1.5 (Tern only on 10 paths).

Figure regenerated from the turn-8 numbers:

```
julia --project=. agent/scripts/chartkit.jl agent/studies/launch_supply/figures/options_vs_wait.json --out=agent/studies/launch_supply/figures/options_vs_wait.svg
chartkit: wrote agent/studies/launch_supply/figures/options_vs_wait.svg — bars_error, 12.9 KiB, self-contained
exit 0
```

(A first render with a subtitle reading "error bars 2 standard errors" was replaced after reading chartkit.jl, which draws ±1 se dark and ±2 se grey whiskers.)
