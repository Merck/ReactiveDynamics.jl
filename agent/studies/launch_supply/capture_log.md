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
