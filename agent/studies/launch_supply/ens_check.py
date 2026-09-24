"""Turn 8: the late-filing mismatch against the relayed 96-seed wait ensemble (root seed 2609).
Ensemble member k is seeded hash((2609, k)); three both-succeed members were re-run with that seed on
the unforced model and on the forced both-succeed copy (runs_ens/). Relayed per-sample facts used here
(verbatim from the coordinator's relay, turn 8): both-succeed members 20; t_file_delay among them
2 of 20 (2 months each); t_short 14 of 20; members 20/96/6 had enpv 1526.7/1444.8/1529.5,
t_file_delay 2/2/0, t_lost 92/59/91."""
import os, sys, json, math
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv = sys.argv[:1]
from analyze import load as aload
from paired_both import wilson, seeds, load
REL = {20: (1526.7, 2.0, 92.0), 96: (1444.8, 2.0, 59.0), 6: (1529.5, 0.0, 91.0)}
for k, (e, dl, tl) in REL.items():
    nat, frc = (aload(os.path.join(HERE, "runs_ens", f"model__{w}__m{k}.json")) for w in ("natural", "both"))
    same = all(abs(nat[x] - frc[x]) < 1e-9 for x in ("npv", "t_file_delay", "t_lost", "t_short"))
    print(f"member {k}: forced == unforced {same}; our run delay {nat['t_file_delay']:.0f} t_lost {nat['t_lost']:.0f} (ensemble {dl:.0f} / {tl:.0f}); "
          f"our eNPV {nat['npv']:.1f} vs ensemble {e} -> ratio {e/nat['npv']:.5f} (1.09^(1/12) = {1.09**(1/12):.5f})")
S = seeds("model", "both")
ours = sum(load("model", "both", s)["t_file_delay"] > 0 for s in S); n = len(S)
N, K, m = n + 20, ours + 2, 20
p = [comb(K, x) * comb(N - K, m - x) / comb(N, m) for x in range(m + 1)]
print(f"wait, both succeed: late filing ours {ours}/{n}, ensemble 2/20; Fisher exact two-sided p = {sum(q for q in p if q <= p[2] + 1e-12):.3f}")
lo, hi = wilson(K, N); print(f"pooled {K}/{N} = {K/N:.3f} (95% Wilson {lo:.3f}-{hi:.3f})")
ts = sum(load("model", "both", s)["t_short"] > 0 for s in S); lo2, hi2 = wilson(ts + 14, n + 20)
print(f"Tern short at launch, wait, both succeed: ours {ts}/{n}, ensemble 14/20, pooled {ts+14}/{n+20} = {(ts+14)/(n+20):.3f} (95% Wilson {lo2:.3f}-{hi2:.3f}); x 0.27 = {0.27*(ts+14)/(n+20):.3f} ({0.27*lo2:.3f}-{0.27*hi2:.3f})")
r = sum(load("arm_cmo", "both", s)["t_file_delay"] > 0 for s in S); lo3, hi3 = wilson(r, n)
print(f"reserve now, both succeed: late filing {r}/{n} = {r/n:.3f} (95% Wilson {lo3:.3f}-{hi3:.3f}); no ensemble exists for this option")
