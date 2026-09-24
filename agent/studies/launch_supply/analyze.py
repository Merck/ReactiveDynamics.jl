"""Scenario-weighted analysis of the launch_supply forced-outcome paths (turn 6: no ensembles).

Every number here is a SINGLE simulated path per (option, outcome world, seed), read from
runs/<arm>__<world>__s<seed>.json (the verbatim `rdx run` output). eNPV is recomputed from the
traced monthly flows with the same discounting as `rdx enpv` (flow / 1.09^(t/12)):
  reward_t = (k_sold_t - k_sold_{t-1}) * k_val + (t_sold_t - t_sold_{t-1}) * t_val
  cost_t   = (cash_{t-1} - cash_t) * 0.001
and checked against the run's undiscounted ledger totals. The expected value of an option is the
average over the available seeds within each world, weighted by the stated readout odds (Kestrel p_k,
Tern p_t, independent). It is arithmetic on a handful of paths, not an ensemble mean: no standard
error is claimed. Usage: python3 analyze.py [--peak] [--pk=0.5]
"""
import json, glob, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gen_models import VAL

PEAK = "--peak" in sys.argv
PK = next((float(a.split("=")[1]) for a in sys.argv if a.startswith("--pk=")), 0.60)
PT = 0.45
KV, TV = VAL["peak" if PEAK else "year1"]
W = {"both": PK * PT, "konly": PK * (1 - PT), "tonly": (1 - PK) * PT, "neither": (1 - PK) * (1 - PT)}


def series(tr, n):
    x = tr[n]
    return x if isinstance(x, list) else x["values"]


def load(path):
    r = json.load(open(path))["result"]
    tr, pl = r["trace"], r["places"]
    t, cash, ks, ts = series(tr, "t"), series(tr, "cash"), series(tr, "k_sold"), series(tr, "t_sold")
    npv = rew = cost = 0.0
    for i in range(1, len(t)):
        rw = (ks[i] - ks[i - 1]) * KV + (ts[i] - ts[i - 1]) * TV
        cs = (cash[i - 1] - cash[i]) * 0.001
        rew += rw; cost += cs
        npv += (rw - cs) / 1.09 ** (t[i] / 12)
    if not PEAK:  # the ledger check only holds at the values baked into the document
        led = r["ledger"]
        assert abs(led["reward"] - rew) < 1e-6 * max(1, rew) + 1e-6, (path, led["reward"], rew)
        assert abs(led["cost"] - cost) < 1e-6 * max(1, cost) + 1e-3, (path, led["cost"], cost)
    ka = series(tr, "k_age")
    f = {k: pl[k]["final"] for k in ("k_short", "k_short_any", "t_short", "t_file_delay", "k_lost", "t_lost")}
    f["npv"] = npv; f["cost"] = cost
    return f


runs = defaultdict(lambda: defaultdict(dict))   # arm -> world -> seed -> metrics
for p in sorted(glob.glob(os.path.join(HERE, "runs", "*.json"))):
    arm, world, seed = os.path.basename(p)[:-5].split("__")
    runs[arm][world][seed] = load(p)


def world_mean(arm, world, key, seeds=None):
    d = runs[arm][world]
    ss = [s for s in d if seeds is None or s in seeds]
    return sum(d[s][key] for s in ss) / len(ss), ss


def fallback(arm, world):
    # sensitivities differ from their base arm only when both succeed (no2: by construction, the
    # second contract is only signed then) — other worlds are borrowed from the base arm, stated
    if runs[arm][world]:
        return arm
    base = arm.replace("_literal", "").replace("_no2", "")
    return base


if __name__ == "__main__":
    print(f"weights (p_k={PK}, p_t={PT}):", {k: round(v, 4) for k, v in W.items()}, "| value per tenth-batch", (KV, TV))
    arms = [a for a in ["model", "arm_cmo", "arm_cmo_tern_only", "arm_build_now", "arm_build_after",
                        "model_literal", "arm_cmo_literal", "arm_cmo_tern_only_literal",
                        "model_no2", "arm_cmo_no2", "arm_cmo_tern_only_no2"] if runs.get(a)]
    print("\nper world, mean over seeds shown (single paths when one seed): eNPV $M | K short launch / K short ever / T short launch / T filing delay (months)")
    for w in W:
        for a in arms:
            src = fallback(a, w)
            if not runs[src][w]:
                continue
            npv, ss = world_mean(src, w, "npv")
            vals = [world_mean(src, w, k)[0] for k in ("k_short", "k_short_any", "t_short", "t_file_delay")]
            tag = "" if src == a else f"  (= {src})"
            print(f"  {w:8s} {a:28s} seeds {','.join(ss):8s} eNPV {npv:8.1f} | {vals[0]:.2f} {vals[1]:.2f} {vals[2]:.2f} {vals[3]:.2f}{tag}")
    print("\nweighted by the odds (arithmetic over the paths above; not an ensemble):")
    base = {}
    for a in arms:
        tot = 0.0; ok = True
        pk_short = pt_short = pk_any = p_late = 0.0
        for w, wt in W.items():
            src = fallback(a, w)
            if not runs[src][w]:
                ok = False; break
            tot += wt * world_mean(src, w, "npv")[0]
            pk_short += wt * world_mean(src, w, "k_short")[0]
            pk_any += wt * world_mean(src, w, "k_short_any")[0]
            pt_short += wt * world_mean(src, w, "t_short")[0]
            p_late += wt * (sum(1.0 for s in runs[src][w].values() if s["t_file_delay"] > 0) / len(runs[src][w]))
        if not ok:
            continue
        ref = "model_literal" if a.endswith("_literal") else "model_no2" if a.endswith("_no2") else "model"
        base[a] = tot
        d = f"  vs {ref}: {tot - base[ref]:+7.1f}" if ref in base else ""
        print(f"  {a:28s} eNPV {tot:8.1f}{d} | P(K short launch) {pk_short:.2f}  P(K short ever) {pk_any:.2f}  P(T short launch) {pt_short:.2f}  P(T files late) {p_late:.2f}")
