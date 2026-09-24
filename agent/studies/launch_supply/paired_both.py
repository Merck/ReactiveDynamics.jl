"""Turn 7: the both-succeed world on many paired paths, plus decompositions for the user's questions.
Every number is arithmetic over runs/*.json (verbatim `rdx run` output). Paired = same seed, so same
readout dates, slips and lead-time draws in both options. Standard errors are across paths (sd/sqrt(n))."""
import json, glob, os, sys, math
# flags pass through to analyze: --peak (peak value per batch), --pk=0.5 (Kestrel odds)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analyze import series, W, KV, TV

def load(a, w, s):
    p = os.path.join(HERE, "runs", f"{a}__{w}__{s}.json")
    if not os.path.exists(p): return None
    r = json.load(open(p))["result"]; tr = r["trace"]
    t = series(tr, "t"); g = lambda n: series(tr, n)
    ks, ts, cash, tl, kl, ta, ka = g("k_sold"), g("t_sold"), g("cash"), g("t_lost"), g("k_lost"), g("t_age"), g("k_age")
    d = lambda i: 1.09 ** (-t[i] / 12)
    o = dict(kpv=0.0, tpv=0.0, cpv=0.0, tlost_pv=0.0, t_short_months=[], )
    for i in range(1, len(t)):
        o["kpv"] += (ks[i] - ks[i-1]) * KV * d(i); o["tpv"] += (ts[i] - ts[i-1]) * TV * d(i)
        o["cpv"] += (cash[i-1] - cash[i]) * 0.001 * d(i)
        if tl[i] > tl[i-1]:
            o["tlost_pv"] += (tl[i] - tl[i-1]) * TV * d(i); o["t_short_months"].append(int(t[i]))
    o["npv"] = o["kpv"] + o["tpv"] - o["cpv"]
    o["t_launch"] = next((int(t[i]) for i, v in enumerate(ta) if v > 0), None)
    o["k_launch"] = next((int(t[i]) for i, v in enumerate(ka) if v > 0), None)
    fin = r["places"]
    for k in ("t_short", "k_short", "k_short_any", "t_file_delay", "t_lost", "k_lost"): o[k] = fin[k]["final"]
    return o

def mse(xs):
    n = len(xs); m = sum(xs) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1)) if n > 1 else float("nan")
    return m, sd / math.sqrt(n), n

def wilson(k, n, z=1.96):
    p = k / n; den = 1 + z*z/n; c = (p + z*z/(2*n)) / den; h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return c - h, c + h

seeds = lambda a, w: sorted({os.path.basename(p).split("__")[2][:-5] for p in glob.glob(os.path.join(HERE, "runs", f"{a}__{w}__s*.json"))}, key=lambda s: int(s[1:]))

if __name__ == "__main__":
    S = [s for s in seeds("model", "both") if s in seeds("arm_cmo", "both")]
    M = {s: load("model", "both", s) for s in S}; C = {s: load("arm_cmo", "both", s) for s in S}
    print(f"BOTH SUCCEED, wait vs reserve-now, {len(S)} paired paths ({S[0]}..{S[-1]})")
    print(" seed  d_eNPV |  dK_rev  dT_rev  d_cost | wait: Tlaunch delay Tshort Tlost(batches) TlostPV months_short | cmo: Tlaunch delay Tshort Tlost")
    for s in S:
        m, c = M[s], C[s]
        ms = m["t_short_months"]; rng = f"{ms[0]}-{ms[-1]} ({len(ms)})" if ms else "-"
        print(f" {s:4s} {c['npv']-m['npv']:+7.1f} | {c['kpv']-m['kpv']:+7.1f} {c['tpv']-m['tpv']:+7.1f} {-(c['cpv']-m['cpv']):+7.1f} | "
              f"{m['t_launch']} {m['t_file_delay']:.0f} {m['t_short']:.0f} {m['t_lost']/10:5.1f} {m['tlost_pv']:6.1f} {rng:14s} | {c['t_launch']} {c['t_file_delay']:.0f} {c['t_short']:.0f} {c['t_lost']/10:.1f}")
    D = [C[s]["npv"] - M[s]["npv"] for s in S]
    m, se, n = mse(D); print(f"\n d_eNPV (reserve - wait | both): mean {m:+.1f} se {se:.1f} n {n}; min {min(D):+.1f} max {max(D):+.1f}; negative on {sum(x<0 for x in D)} of {n}")
    for lab, X in (("wait", M), ("reserve", C)):
        ks = sum(X[s]["t_short"] > 0 for s in S); kl = sum(X[s]["t_file_delay"] > 0 for s in S); ka = sum(X[s]["t_lost"] > 0 for s in S)
        dl = [X[s]["t_file_delay"] for s in S if X[s]["t_file_delay"] > 0]
        lo, hi = wilson(ks, n); lo2, hi2 = wilson(kl, n)
        tl = mse([X[s]["t_lost"] / 10 for s in S]); tpv = mse([X[s]["tlost_pv"] for s in S])
        print(f" {lab:8s} Tern short at launch {ks}/{n} = {ks/n:.2f} (95% Wilson {lo:.2f}-{hi:.2f}) | files late {kl}/{n} = {kl/n:.2f} ({lo2:.2f}-{hi2:.2f}), months {sorted(dl)} | "
              f"any Tern demand lost {ka}/{n}; lost batches mean {tl[0]:.1f} se {tl[1]:.1f}; lost value PV mean {tpv[0]:.1f} se {tpv[1]:.1f} $M | Kestrel short ever {sum(X[s]['k_short_any']>0 for s in S)}/{n}, at launch {sum(X[s]['k_short']>0 for s in S)}/{n}")
    # the late-filing mechanism in the reserve arm: Tern's filing waits for the CMO PPQ
    print("\n reserve-now late filings (seed: months late, Tern launch reserve vs wait):", [(s, C[s]["t_file_delay"], C[s]["t_launch"], M[s]["t_launch"]) for s in S if C[s]["t_file_delay"] > 0])
    # other worlds
    print("\nSINGLE-ASSET AND NEITHER WORLDS, paired, seeds available in both options")
    other = {}
    for w in ("konly", "tonly", "neither"):
        ss = [s for s in seeds("model", w) if s in seeds("arm_cmo", w)]
        dd = []
        for s in ss:
            a, b = load("model", w, s), load("arm_cmo", w, s); dd.append(b["npv"] - a["npv"])
            print(f" {w:8s} {s}: d_eNPV {b['npv']-a['npv']:+6.1f} = dK {b['kpv']-a['kpv']:+6.1f} dT {b['tpv']-a['tpv']:+6.1f} dCost {-(b['cpv']-a['cpv']):+6.1f} | Tlaunch wait {a['t_launch']} reserve {b['t_launch']}, delay wait {a['t_file_delay']:.0f} reserve {b['t_file_delay']:.0f}")
        other[w] = mse(dd)
    print("\nWEIGHTED reserve-now vs wait (odds weights", {k: round(v, 2) for k, v in W.items()}, ")")
    mb = mse(D)
    tot = W["both"] * mb[0] + sum(W[w] * other[w][0] for w in other)
    var = (W["both"] * mb[1]) ** 2 + sum((W[w] * other[w][1]) ** 2 for w in other if other[w][2] > 1)
    print(f" both {mb[0]:+.1f} (se {mb[1]:.1f}, n {mb[2]}) | " + " | ".join(f"{w} {other[w][0]:+.1f} (se {other[w][1]:.2f}, n {other[w][2]})" for w in other))
    print(f" weighted delta {tot:+.1f}, se {math.sqrt(var):.1f}; delta - 2se {tot-2*math.sqrt(var):+.1f}, delta + 2se {tot+2*math.sqrt(var):+.1f}")
    ks = sum(M[s]["t_short"] > 0 for s in S); lo, hi = wilson(ks, n)
    print(f" P(Tern short at launch) wait = {W['both']:.3f} x {ks}/{n} = {W['both']*ks/n:.3f} ({W['both']:.3f} x Wilson: {W['both']*lo:.3f}-{W['both']*hi:.3f}); reserve = {W['both']:.3f} x {sum(C[s]['t_short']>0 for s in S)}/{n}")
    # build options, both world
    for a in ("arm_build_now", "arm_build_after"):
        ss = [s for s in seeds(a, "both") if s in S]
        dd = [load(a, "both", s)["npv"] - M[s]["npv"] for s in ss]; mm = mse(dd)
        print(f" {a}: both-world d vs wait mean {mm[0]:+.1f} se {mm[1]:.1f} n {mm[2]} range {min(dd):+.1f}..{max(dd):+.1f}; Tern short {sum(load(a,'both',s)['t_short']>0 for s in ss)}/{len(ss)}")
    # model level in the both world, for the ensemble cross-check
    mm = mse([M[s]["npv"] for s in S]); print(f"\n wait eNPV | both: mean {mm[0]:.1f} se {mm[1]:.1f} n {mm[2]}")
    for w in ("konly", "tonly", "neither"):
        v = mse([load("model", w, s)["npv"] for s in seeds("model", w)]); print(f" wait eNPV | {w}: mean {v[0]:.1f} n {v[2]}")
    # weighted late filing, with the path count behind each world
    print("\nWEIGHTED P(Tern files late) and mean months late, by world (paths):")
    for a in ("model", "arm_cmo"):
        tot = mon = 0.0; parts = []
        for w in ("both", "konly", "tonly", "neither"):
            ss = S if w == "both" else [s for s in seeds("model", w) if s in seeds("arm_cmo", w)]
            L = [load(a, w, s)["t_file_delay"] for s in ss]
            k = sum(x > 0 for x in L); tot += W[w] * k / len(L); mon += W[w] * sum(L) / len(L); parts.append(f"{w} {k}/{len(L)}")
        print(f" {a:8s} P(late) {tot:.3f}  mean months late {mon:.3f}  [{', '.join(parts)}]")
    # cross-check against the relayed 96-seed wait ensemble (root seed 2609): the forced-world means
    # reweighted at the success rates that ensemble actually realized (independence assumed)
    ek, et = 0.5417, 0.3542
    Wr = {"both": ek * et, "konly": ek * (1 - et), "tonly": (1 - ek) * et, "neither": (1 - ek) * (1 - et)}
    lv = {"both": mm[0]}; lv.update({w: mse([load("model", w, s)["npv"] for s in seeds("model", w)])[0] for w in ("konly", "tonly", "neither")})
    tsh = sum(M[s]["t_short"] > 0 for s in S) / n; dly = sum(M[s]["t_file_delay"] for s in S) / n
    dly_t = sum(W["tonly"] * 0 for _ in [0])
    print(f"\nENSEMBLE CROSS-CHECK (wait): forced-world eNPV at realized rates {sum(Wr[w]*lv[w] for w in Wr):.1f} vs ensemble 765.9561 (sem 61.6450); at design odds {sum(W[w]*lv[w] for w in W):.1f}")
    print(f"  t_short: {Wr['both']:.4f} x {tsh:.3f} = {Wr['both']*tsh:.3f} vs ensemble 0.1458 (sem 0.0362)")
    print(f"  t_file_delay months: {Wr['both']:.4f} x {dly:.3f} = {Wr['both']*dly:.3f} vs ensemble 0.0417 (sem 0.0293)")
    ka = sum(M[s]["k_short_any"] > 0 for s in S) / n
    print(f"  k_short_any: {Wr['both']:.4f} x {ka:.3f} = {Wr['both']*ka:.3f} vs ensemble 0.0208 (sem 0.0147); t_lost tenths: {Wr['both']:.4f} x {sum(M[s]['t_lost'] for s in S)/n:.1f} = {Wr['both']*sum(M[s]['t_lost'] for s in S)/n:.1f} vs ensemble 24.0417 (sem 4.9599)")
