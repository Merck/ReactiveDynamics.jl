"""Turn 8: 'reserve now, review at Kestrel's readout' (arm_cmo_exit) against wait, and the late-filing
cross-check against the relayed 96-seed wait ensemble. Arithmetic over runs/*.json and runs_ens/*.json
(verbatim `rdx run` output). The exit rule only acts when Kestrel fails, so in the both-succeed and
Kestrel-only worlds arm_cmo_exit IS arm_cmo (checked byte-for-byte on seed 1 of each, and on the
neither world); those worlds reuse the arm_cmo paths. Flags as analyze.py: --peak, --pk=0.5."""
import os, sys, math, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from analyze import W, load as aload
from paired_both import load, mse, wilson, seeds

def paired(a, b, w):
    ss = [s for s in seeds(a, w) if s in seeds(b, w)]
    return ss, [load(b, w, s)["npv"] - load(a, w, s)["npv"] for s in ss]

if __name__ == "__main__":
    print("inertness check, arm_cmo_exit vs arm_cmo on the same seed (must be 0):",
          [(w, s, round(load("arm_cmo_exit", w, s)["npv"] - load("arm_cmo", w, s)["npv"], 6)) for w, s in (("both", "s1"), ("konly", "s1"), ("neither", "s1"))])
    print("\nTERN-ONLY WORLD (Kestrel fails, Tern succeeds), paired against wait")
    for b in ("arm_cmo_exit", "arm_cmo"):
        ss, d = paired("model", b, "tonly"); m, se, n = mse(d)
        late = [load(b, "tonly", s)["t_file_delay"] for s in ss]; tl = [load(b, "tonly", s)["t_lost"] for s in ss]
        mlate = [load("model", "tonly", s)["t_file_delay"] for s in ss]
        print(f" {b:13s} vs wait: mean {m:+.2f} se {se:.2f} n {n}; min {min(d):+.2f} max {max(d):+.2f}; paths below -20: {sum(x < -20 for x in d)}/{n}"
              f" | Tern files late {sum(x > 0 for x in late)}/{n} (wait {sum(x > 0 for x in mlate)}/{n}); Tern demand lost on {sum(x > 0 for x in tl)}/{n}")
    ss, d = paired("model", "arm_cmo_exit", "tonly")
    dc = [-(load("arm_cmo_exit", "tonly", s)["cpv"] - load("model", "tonly", s)["cpv"]) for s in ss]
    dt = [load("arm_cmo_exit", "tonly", s)["tpv"] - load("model", "tonly", s)["tpv"] for s in ss]
    print(f"  arm_cmo_exit decomposition: cost diff mean {sum(dc)/len(dc):+.2f} (min {min(dc):+.2f}, max {max(dc):+.2f}); Tern revenue diff mean {sum(dt)/len(dt):+.2f} (min {min(dt):+.2f} max {max(dt):+.2f})")
    print("\nWEIGHTED 'reserve now, review at Kestrel readout' vs wait, odds", {k: round(v, 3) for k, v in W.items()})
    src = {"both": "arm_cmo", "konly": "arm_cmo", "tonly": "arm_cmo_exit", "neither": "arm_cmo"}
    tot = var = 0.0; parts = []
    for w, a in src.items():
        ss, d = paired("model", a, w); m, se, n = mse(d)
        tot += W[w] * m; var += (W[w] * se) ** 2 if n > 1 else 0.0
        parts.append(f"{w} {m:+.2f} (se {se:.2f}, n {n}, from {a})")
    print("  " + " | ".join(parts))
    print(f"  weighted delta {tot:+.2f}, se {math.sqrt(var):.2f}; -2se {tot-2*math.sqrt(var):+.2f}, +2se {tot+2*math.sqrt(var):+.2f}")
    worst = min((W[w], w) for w in ("konly", "tonly", "neither"))
    print("  loss given each non-both world (mean):", {w: round(mse(paired('model', src[w], w)[1])[0], 2) for w in ("konly", "tonly", "neither")})
    # probabilities with path counts
    S = seeds("model", "both")
    ks = sum(load("model", "both", s)["t_short"] > 0 for s in S); lo, hi = wilson(ks, len(S))
    print(f"  P(Tern short at launch): wait {W['both']:.3f} x {ks}/{len(S)} = {W['both']*ks/len(S):.3f} (Wilson {W['both']*lo:.3f}-{W['both']*hi:.3f}); "
          f"reserve {W['both']:.3f} x {sum(load('arm_cmo','both',s)['t_short']>0 for s in S)}/{len(S)}; Tern-only world short: wait {sum(load('model','tonly',s)['t_short']>0 for s in seeds('model','tonly'))}/{len(seeds('model','tonly'))}, review {sum(load('arm_cmo_exit','tonly',s)['t_short']>0 for s in seeds('arm_cmo_exit','tonly'))}/{len(seeds('arm_cmo_exit','tonly'))}")
