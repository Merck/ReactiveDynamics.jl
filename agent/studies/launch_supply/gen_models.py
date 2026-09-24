"""Generate the launch_supply model documents (pure JSON, no host.jl).

Every arm is the SAME network; arms differ only in `params` values (policy switches and
sensitivities), so `diff model.rdj.json arm_X.rdj.json` is the lever. Run from anywhere:
    python3 agent/studies/launch_supply/gen_models.py
Units: one tick = one month, t = 0 is 1 Oct 2026, horizon t = 75 (end 2032). Capacity and
stock in TENTHS of a batch (every batch transition moves one tenth, so line capacity is continuous
and no fraction of a month is wasted). Cash in $0.001M units (place cost 0.001 $M per unit).
Rewards in $M.
Every number traces to the user's turns 1-4 (see BRIEF.md / capture_log.md).
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))

# ── node helpers ──────────────────────────────────────────────────────────────────────────────
def C(v): return {"node": "const", "value": float(v)}
def P(n): return {"node": "ref", "kind": "param", "name": n}
def S(n): return {"node": "ref", "kind": "place", "name": n}
T = {"node": "timeref"}
def op(o, *a): return {"node": "call", "op": o, "args": [x if isinstance(x, dict) else C(x) for x in a]}
def add(*a):
    out = a[0]
    for x in a[1:]: out = op("+", out, x)
    return out
def mul(*a):
    out = a[0]
    for x in a[1:]: out = op("*", out, x)
    return out
def sub(a, b): return op("-", a, b)
def div(a, b): return op("/", a, b)
def mx(a, b): return op("max", a, b)
def mn(a, b): return op("min", a, b)
def fl(a): return op("floor", a)
def ce(a): return op("ceil", a)
def ge(a, b): return op(">=", a, b)
def gt(a, b): return op(">", a, b)
def lt(a, b): return op("<", a, b)
def eq(a, b): return op("==", a, b)
def AND(*a):
    out = a[0]
    for x in a[1:]: out = op("&&", out, x)
    return out
def OR(a, b): return op("||", a, b)
def NOT(a): return op("!", a)
def on(x): return gt(x, 0.5)          # boolean view of a 0/1 param or flag
def ind(b): return mul(1.0, b)        # numeric 0/1 from a boolean (Bool * Float)
def smp(dist, *a): return {"node": "sample", "dist": dist, "args": [x if isinstance(x, dict) else C(x) for x in a]}

# ── the user's numbers ────────────────────────────────────────────────────────────────────────
LINE_TENTHS = 128.0        # 12.8 effective batches/month (trailing 12 months)
CMO_TENTHS = 40.0          # up to 4 batches/month
LINE2_TENTHS = 120.0       # 12 batches/month effective (engineering pre-study)
MKT0, MKT_G = 9.0, 1.05    # marketed batches/month today, +5%/yr
REJECT_POS = 0.97          # ~3% of batches rejected
MARGIN = 0.85
# value per batch shipped: base = year-1 rate, upside = peak rate (turn 4)
K_Y1_BATCHES = sum(min(1 + a / 6, 3 + (a - 12) / 16, 4.5) for a in range(12))   # 23.0
T_Y1_BATCHES = sum(min(1 + a / 8, 2.5 + (a - 12) / 16, 4.0) for a in range(12))  # 20.25
VAL = {
    "year1": (250.0 / K_Y1_BATCHES * MARGIN / 10, 150.0 / T_Y1_BATCHES * MARGIN / 10),
    "peak": (1200.0 / 54 * MARGIN / 10, 900.0 / 48 * MARGIN / 10),
}

# cumulative marketed demand in tenths, grossed up for rejects, so floor(C(t+1)) - floor(C(t)) is an
# exact integer monthly count with the right long-run mean
LN = math.log(MKT_G)
def cum_mkt(tnode):
    k = 10.0 * MKT0 * 12.0 / LN / REJECT_POS
    return mul(k, sub(op("exp", mul(LN / 12.0, tnode)), 1.0))
MKT_NOW = sub(fl(cum_mkt(add(T, 1.0))), fl(cum_mkt(T)))

def dK(a): return mn(mn(add(1.0, div(a, 6.0)), add(3.0, div(sub(a, 12.0), 16.0))), 4.5)
def dT(a): return mn(mn(add(1.0, div(a, 8.0)), add(2.5, div(sub(a, 12.0), 16.0))), 4.0)

# ── derived timeline (months from Oct 2026) ───────────────────────────────────────────────────
K_READ = add(5.0, P("k_slip"))            # Mar 2027 + 0..2
K_PRE = add(15.0, P("k_slip"))            # launch - 6
K_LAUNCH = add(21.0, P("k_slip"))         # +3 file +12 review +1 launch = Jul 2028
K_APPROVE = add(20.0, P("k_slip"))
T_READ = add(14.0, P("t_shift"))          # Dec 2027 +/- 3
T_APPROVE = add(P("t_file_t"), 12.0)

def base_params():
    return {
        # business uncertainty (drawn at t = 0)
        "p_k": 0.60, "p_t": 0.45,
        "k_ok": 0.0, "t_ok": 0.0, "k_slip": 0.0, "t_shift": 0.0, "cmo_lead": 12.0, "build_lead": 30.0, "cmo2_lead": 12.0,
        "cmo2_start": 999.0, "cmo2_ready_t": 999.0,
        # state written by rules
        "t_file_t": 999.0, "cmo_start": 999.0, "cmo_resigned": 0.0, "cmo_ready_t": 999.0, "l2_start": 999.0, "l2_ready_t": 999.0,
        "k_share": 1.0,
        # value per tenth-batch shipped, $M margin
        "k_val": VAL["year1"][0], "t_val": VAL["year1"][1],
        # policy switches (the levers)
        "cmo_now": 0.0,       # 1 = reserve the CMO slot now (sign t = 2)
        "cmo_reassign": 1.0,  # 1 = the reserved slot can carry Kestrel (term sheet reserves capacity)
        "wait_cmo": 1.0,      # 1 = if either readout succeeds, take the CMO's mid-2028 opening
        "build_now": 0.0,     # 1 = approve the second line at ExCo (t = 2)
        "build_after": 0.0,   # 1 = approve the second line after Tern's readout if both succeed
        "k_first": 0.0,       # 0 = base (turn 5): Kestrel's launch build + this month's demand + 2-month floor,
                              #     then Tern, then Kestrel's buffer; 1 = literal rule: Kestrel's full
                              #     six-month buffer before any Tern batch
        "k_floor": 2.0,       # months of forward stock Kestrel never drops below (turn 5 guardrail)
        "second_cmo": 1.0,    # 1 = if both succeed, sign a second CMO contract on the same terms (turn 5)
        "walk2": 1.0,         # 1 = walking away after failure costs two annual fees ($10M), 0 = one ($5M)
        "capex_even": 0.0,    # 1 = spread the $120M evenly over 36 months instead of 20/50/30
    }

def build(overrides, name):
    prm = base_params(); prm.update(overrides)

    places = [
        {"name": "cash", "init": 1.0e12, "cost": 0.001},         # $0.001M units; deep, must never bind
        {"name": "line", "init": LINE_TENTHS},
        {"name": "cmo"}, {"name": "line2"}, {"name": "cmo2"},
        {"name": "cmo2_active"}, {"name": "cmo2_ready"}, {"name": "k_cmo2_ok"}, {"name": "t_cmo2_ok"},
        {"name": "k_short_any"},                                   # Kestrel lost demand at any time (turn 5)
        {"name": "mkt_done"},
        {"name": "k_stock"}, {"name": "t_stock"}, {"name": "k_orders"}, {"name": "t_orders"},
        {"name": "k_sold", "reward": 0.0}, {"name": "t_sold", "reward": 0.0},  # reward set below
        {"name": "k_lost"}, {"name": "t_lost"},
        {"name": "k_pre"}, {"name": "k_launched"}, {"name": "k_age"},
        {"name": "t_pre"}, {"name": "t_launched"}, {"name": "t_age"},
        {"name": "t_ppq_done"}, {"name": "t_filed"}, {"name": "t_file_delay"},
        {"name": "k_cmo_ok"}, {"name": "t_cmo_ok"}, {"name": "k_l2_ok"}, {"name": "t_l2_ok"},
        {"name": "cmo_active"}, {"name": "cmo_ppq"}, {"name": "l2_active"}, {"name": "l2_ready"},
        {"name": "k_short"}, {"name": "t_short"},
        {"name": "k_success"}, {"name": "t_success"},             # mirrors of the draws, for pairing checks
        {"name": "own_used"}, {"name": "cmo_used"}, {"name": "l2_used"},  # tallies of capacity-months used
        {"name": "own_free"}, {"name": "cmo_free"}, {"name": "l2_free"}, {"name": "cmo2_free"}, {"name": "cmo2_used"},  # planner mirrors: unallocated capacity
    ]
    # reward is a static per-place column: bake the arm's value per tenth-batch in
    for p in places:
        if p["name"] == "k_sold": p["reward"] = prm["k_val"]
        if p["name"] == "t_sold": p["reward"] = prm["t_val"]

    CONS = {"allocation": "upfront", "return": "conserved", "blocking": "block"}
    trs, arcs = [], []
    def tr(id_, rate, cycletime=0.0, pos=1.0, priority=1.0):
        trs.append({"id": id_, "name": id_, "rate": rate, "rate_mode": "deterministic",
                    "cycletime": cycletime, "prob_of_success": pos, "priority": priority})
    def arc(t, place, side, m=1.0, conserved=False):
        a = {"transition": t, "place": place, "side": side, "multiplicity": float(m)}
        if conserved: a["modality"] = CONS
        arcs.append(a)

    # ── the production planner ──────────────────────────────────────────────────────────────
    # A rule at the END of every tick computes next tick's production proposals into params, step
    # by step (each set_params sees the ones before it), so shared sub-expressions stay small. The
    # markings it reads are exactly those next tick's genesis sees (nothing moves in between), and
    # the line is full again then: a batch holds its site for the month and returns it in-tick.
    T1 = add(T, 1.0)
    MKT_NEXT = sub(fl(cum_mkt(add(T1, 1.0))), fl(cum_mkt(T1)))
    plan = []
    def let(name, node):
        plan.append({"verb": "set_params", "assigns": [{"name": name, "value": node}]})
        prm.setdefault(name, 0.0)
        return P(name)
    free = {"own": let("pl_free_own", mx(0.0, sub(S("line"), MKT_NEXT))), "l2": S("line2"), "cmo": S("cmo"), "cmo2": S("cmo2")}
    k_gap = let("pl_k_gap", mx(0.0, ce(sub(mul(60.0, P("k_share"), dK(add(S("k_age"), 6.0))), S("k_stock")))))
    t_gap = let("pl_t_gap", mx(0.0, ce(sub(mul(60.0, dT(add(S("t_age"), 6.0))), S("t_stock")))))
    # Kestrel's need in two tiers. `k_first = 1` (the user's literal rule): Kestrel's whole need,
    # buffer included, before any Tern batch. `k_first = 0`: Kestrel's launch build and this month's
    # demand first, then Tern, then Kestrel's buffer top-up.
    # Guardrail (both versions): never below `k_floor` months of forward stock. Under the literal
    # rule the six-month gap already covers it.
    k_dem = ce(mul(10.0, P("k_share"), dK(S("k_age"))))
    k_floor_gap = let("pl_k_floor", mx(0.0, ce(sub(mul(10.0, P("k_floor"), P("k_share"), dK(add(S("k_age"), P("k_floor")))), S("k_stock")))))
    k_base_need = let("pl_k_base", add(mul(P("k_first"), k_gap), mul(sub(1.0, P("k_first")),
                      add(mul(sub(1.0, S("k_launched")), k_gap), mul(S("k_launched"), mn(k_gap, mx(k_dem, k_floor_gap)))))))
    k_buf_need = let("pl_k_buf", mx(0.0, sub(k_gap, k_base_need)))
    ppq_need = let("pl_ppq_need", mul(sub(1.0, P("cmo_now")), P("t_ok"), ind(ge(T1, T_READ)), sub(1.0, S("t_filed")),
                   mx(0.0, sub(30.0, S("t_ppq_done")))))
    t_own_ok = mul(S("t_pre"), sub(1.0, P("cmo_now")))
    ok = {
        "k": {"own": S("k_pre"), "l2": S("k_l2_ok"), "cmo": S("k_cmo_ok"), "cmo2": S("k_cmo2_ok")},
        "t": {"own": t_own_ok, "l2": S("t_l2_ok"), "cmo": S("t_cmo_ok"), "cmo2": S("t_cmo2_ok")},
    }
    rates = {}
    def fill(tag, asset, need):
        rem = need
        for site in ("own", "l2", "cmo", "cmo2"):
            use = let(f"r_{tag}_{site}", mul(ok[asset][site], mn(rem, free[site])))
            rates[(tag, site)] = use
            free[site] = let(f"pl_free_{site}", mx(0.0, sub(free[site], use)))
            rem = let(f"pl_rem_{tag}", mx(0.0, sub(rem, use)))
    # tiers in priority order, each filled site by site (own line, line2, CMO: cheapest first)
    fill("k_base", "k", k_base_need)
    ppq_use = let("r_ppq", mn(ppq_need, free["own"]))           # Tern's PPQ, ahead of Tern's stock
    free["own"] = let("pl_free_own", mx(0.0, sub(free["own"], ppq_use)))
    fill("t_make", "t", t_gap)
    fill("k_buf", "k", k_buf_need)
    # mirrors of the capacity left unallocated for next month: 0 = the site binds. A batch holds its
    # site for one month and returns it inside that month, so the site itself never reads empty in
    # the end-of-month snapshot; these are what `diagnose --expect-binds` can see.
    for site in ("own", "l2", "cmo", "cmo2"):
        plan.append({"verb": "set_marking", "name": f"{site}_free", "value": free[site], "mode": "set"})

    PM = 10000.0
    PRIO = {"k_base": 100.0, "t_make": 1.0, "k_buf": 0.5}

    # marketed products — never bumped
    tr("mkt", MKT_NOW, 1.0, 1.0, PM)
    arc("mkt", "line", "lhs", 1, True); arc("mkt", "mkt_done", "rhs", 1); arc("mkt", "own_used", "rhs", 1)

    # batch transitions (each firing = one tenth of a batch on a site for one month; 3% rejected)
    for tag in ("k_base", "t_make", "k_buf"):
        asset = tag[0]
        for site, place, varcost in (("own", "line", 15), ("l2", "line2", 15), ("cmo", "cmo", 40), ("cmo2", "cmo2", 40)):  # $0.15M / $0.4M per batch
            tid = f"{tag}_{site}"
            tr(tid, rates[(tag, site)], 1.0, REJECT_POS, PRIO[tag])
            arc(tid, place, "lhs", 1, True)
            arc(tid, "cash", "lhs", varcost)
            arc(tid, f"{asset}_stock", "rhs", 1)
            arc(tid, {"line": "own_used", "cmo": "cmo_used", "line2": "l2_used", "cmo2": "cmo2_used"}[place], "rhs", 1)

    # Tern PPQ on the own line (only when Tern is NOT filed on the reserved CMO)
    tr("t_ppq_own", ppq_use, 1.0, REJECT_POS, PRIO["t_make"])
    arc("t_ppq_own", "line", "lhs", 1, True); arc("t_ppq_own", "cash", "lhs", 15)
    arc("t_ppq_own", "t_ppq_done", "rhs", 1); arc("t_ppq_own", "own_used", "rhs", 1)

    # demand, shipping, lost sales (unsupplied patients are lost — turn 4)
    for asset, dem, share in (("k", dK, P("k_share")), ("t", dT, C(1.0))):
        tr(f"{asset}_clock", S(f"{asset}_launched"))
        arc(f"{asset}_clock", f"{asset}_age", "rhs", 1)
        tr(f"{asset}_demand", mul(S(f"{asset}_launched"), fl(add(mul(10.0, share, dem(S(f"{asset}_age"))), 0.5))))
        arc(f"{asset}_demand", f"{asset}_orders", "rhs", 1)
        tr(f"{asset}_ship", 1000.0, 0.0, 1.0, 10.0)
        arc(f"{asset}_ship", f"{asset}_orders", "lhs", 1); arc(f"{asset}_ship", f"{asset}_stock", "lhs", 1)
        arc(f"{asset}_ship", f"{asset}_sold", "rhs", 1)
        tr(f"{asset}_lose", 1000.0, 0.0, 1.0, 0.0)
        arc(f"{asset}_lose", f"{asset}_orders", "lhs", 1); arc(f"{asset}_lose", f"{asset}_lost", "rhs", 1)

    # Tern filing delay tally: months past readout+3 that Tern is not yet filed
    tr("t_delay_clock", mul(P("t_ok"), ind(gt(T, add(T_READ, 3.0))), sub(1.0, S("t_filed"))))
    arc("t_delay_clock", "t_file_delay", "rhs", 1)

    # CMO fees ($5M = 5000 units, annually in advance while the contract is live)
    anniv = AND(ge(T, P("cmo_start")), eq(mul(12.0, fl(div(sub(T, P("cmo_start")), 12.0))), sub(T, P("cmo_start"))))
    doomed = AND(NOT(on(P("t_ok"))), OR(NOT(on(P("cmo_reassign"))), NOT(on(P("k_ok")))))
    y2_now = AND(on(P("cmo_now")), eq(T, 14.0))
    tr("cmo_fee", ind(AND(on(S("cmo_active")), anniv, NOT(y2_now))))
    arc("cmo_fee", "cash", "lhs", 5000)
    # year-2 fee on the reserved slot: paid unless the slot is doomed AND walk-away costs one fee
    tr("cmo_fee_y2", ind(AND(y2_now, OR(NOT(doomed), on(P("walk2"))))))
    arc("cmo_fee_y2", "cash", "lhs", 5000)
    # second contract, same terms; never terminated (it is only signed when both assets succeed)
    anniv2 = AND(ge(T, P("cmo2_start")), eq(mul(12.0, fl(div(sub(T, P("cmo2_start")), 12.0))), sub(T, P("cmo2_start"))))
    tr("cmo2_fee", ind(AND(on(S("cmo2_active")), anniv2)))
    arc("cmo2_fee", "cash", "lhs", 5000)

    # second line: capex 20/50/30 over three years from approval (or even), then $10M/yr fixed
    m = sub(T, P("l2_start"))
    live = mul(S("l2_active"), ind(ge(m, 0.0)))
    back = sub(1.0, P("capex_even"))
    for tid, lo, hi, amt, w in (("l2_capex_y1", 0, 12, 2000, back), ("l2_capex_y2", 12, 24, 5000, back),
                                ("l2_capex_y3", 24, 36, 3000, back), ("l2_capex_even", 0, 36, 3333, P("capex_even"))):
        tr(tid, mul(live, w, ind(AND(ge(m, float(lo)), lt(m, float(hi))))))
        arc(tid, "cash", "lhs", amt)   # $M/month x 1000
    tr("l2_fixed", S("l2_ready"))
    arc("l2_fixed", "cash", "lhs", 833)  # $10M/yr

    # ── rules ────────────────────────────────────────────────────────────────────────────────
    def rule(id_, guard, stmts, mode="once"):
        return {"id": id_, "guard": guard, "action": {"verb": "seq", "stmts": stmts}, "fire_mode": mode}
    def setm(n, v, mode="set"): return {"verb": "set_marking", "name": n, "value": v if isinstance(v, dict) else C(v), "mode": mode}
    def setp(**kw): return {"verb": "set_params", "assigns": [{"name": k, "value": v if isinstance(v, dict) else C(v)} for k, v in kw.items()]}

    rules = [
        # 0. all business uncertainty drawn FIRST, identically in every arm → paired arms per seed
        rule("draws", ge(T, 0.0), [
            setp(k_ok=mul(1.0, smp("Bernoulli", P("p_k")))),
            setp(t_ok=mul(1.0, smp("Bernoulli", P("p_t")))),
            setp(k_slip=fl(smp("Uniform", 0.0, 3.0))),        # 0..2 months
            setp(t_shift=fl(smp("Uniform", -3.0, 4.0))),      # -3..+3 months
            setp(cmo_lead=fl(smp("Uniform", 12.0, 16.0))),    # 12..15 months
            setp(build_lead=fl(smp("Uniform", 30.0, 37.0))),  # 30..36 months
            setp(cmo2_lead=fl(smp("Uniform", 12.0, 16.0))),   # second contract, 12..15 months (drawn last: earlier draws unchanged)
            setm("k_success", P("k_ok")), setm("t_success", P("t_ok")),
        ]),
        # Kestrel: own line is the filed site (PPQ Q1 2027)
        rule("k_stock_build", AND(on(P("k_ok")), ge(T, K_PRE)), [setm("k_pre", 1.0)]),
        rule("k_launch", AND(on(P("k_ok")), ge(T, K_LAUNCH)), [setm("k_launched", 1.0)]),
        # Tern filing: 3 months after readout, once its filed site has PPQ done
        rule("t_file", AND(on(P("t_ok")), ge(T, add(T_READ, 3.0)), NOT(on(S("t_filed"))),
                           OR(AND(on(P("cmo_now")), on(S("cmo_ppq"))),
                              AND(NOT(on(P("cmo_now"))), ge(S("t_ppq_done"), 30.0)))),
             [setm("t_filed", 1.0), setp(t_file_t=T)]),
        rule("t_stock_build", AND(on(S("t_filed")), ge(T, add(P("t_file_t"), 7.0))),
             [setm("t_pre", 1.0)] + [setm("t_cmo_ok", P("cmo_now"))]),
        rule("t_launch", AND(on(S("t_filed")), ge(T, add(P("t_file_t"), 13.0))), [setm("t_launched", 1.0)]),
        # CMO: reserve now (sign t = 2) ...
        # (rules act at the END of a tick, so a contract or build starting in month m is decided in
        # month m-1: its first fee / capex month is then charged in month m itself)
        rule("cmo_sign_now", AND(on(P("cmo_now")), ge(T, 1.0)), [setm("cmo_active", 1.0), setp(cmo_start=2.0)]),
        # ... or wait and take the mid-2028 opening if either readout succeeded (both readouts are
        # in by t = 17). A Tern-only reservation terminated after Tern fails takes it for Kestrel.
        rule("cmo_sign_wait", AND(on(P("wait_cmo")), ge(T, 20.0),
                                  OR(AND(NOT(on(P("cmo_now"))), OR(on(P("k_ok")), on(P("t_ok")))),
                                     AND(on(P("cmo_now")), NOT(on(S("cmo_active"))), on(P("k_ok"))))),
             [setm("cmo_active", 1.0), setp(cmo_start=21.0), setp(cmo_resigned=1.0)]),
        rule("cmo_ready", AND(on(S("cmo_active")), NOT(on(P("cmo_resigned"))), ge(T, add(P("cmo_start"), P("cmo_lead")))),
             [setm("cmo", CMO_TENTHS), setm("cmo_ppq", 1.0), setp(cmo_ready_t=T)]),
        rule("cmo_ready_resigned", AND(on(S("cmo_active")), on(P("cmo_resigned")), ge(T, add(P("cmo_start"), P("cmo_lead")))),
             [setm("cmo", CMO_TENTHS), setm("cmo_ppq", 1.0), setp(cmo_ready_t=T)]),
        # walk away from the reserved slot once it has no product
        rule("cmo_terminate", AND(on(P("cmo_now")), on(S("cmo_active")), ge(T, T_READ), doomed),
             [setm("cmo_active", 0.0), setm("cmo", 0.0), setm("cmo_ppq", 0.0)]),
        # Kestrel onto the CMO by PAS: its own tech transfer from readout (reserved slot) or with the
        # new contract (wait), then 6 months after both PPQ and approval
        rule("k_cmo_pas", AND(on(P("k_ok")), on(S("cmo_ppq")),
                              OR(OR(NOT(on(P("cmo_now"))), on(P("cmo_reassign"))), on(P("cmo_resigned"))),
                              ge(T, add(mx(mx(P("cmo_ready_t"), add(K_READ, P("cmo_lead"))), K_APPROVE), 6.0))),
             [setm("k_cmo_ok", 1.0)]),
        # Tern onto a CMO it was not filed with (wait arm) by PAS
        rule("t_cmo_pas", AND(NOT(on(P("cmo_now"))), on(S("cmo_ppq")), on(S("t_filed")),
                              ge(T, add(mx(P("cmo_ready_t"), T_APPROVE), 6.0))),
             [setm("t_cmo_ok", 1.0)]),
        # second line
        rule("l2_approve_now", AND(on(P("build_now")), ge(T, 1.0)), [setm("l2_active", 1.0), setp(l2_start=2.0)]),
        rule("l2_approve_after", AND(on(P("build_after")), on(P("k_ok")), on(P("t_ok")), ge(T, T_READ)),
             [setm("l2_active", 1.0), setp(l2_start=add(T, 1.0))]),
        rule("l2_ready_r", AND(on(S("l2_active")), ge(T, add(P("l2_start"), P("build_lead")))),
             [setm("line2", LINE2_TENTHS), setm("l2_ready", 1.0), setp(l2_ready_t=T)]),
        rule("k_l2_pas", AND(on(P("k_ok")), on(S("l2_ready")), ge(T, add(mx(P("l2_ready_t"), K_APPROVE), 6.0))),
             [setm("k_l2_ok", 1.0)]),
        rule("t_l2_pas", AND(on(S("t_filed")), on(S("l2_ready")), ge(T, add(mx(P("l2_ready_t"), T_APPROVE), 6.0))),
             [setm("t_l2_ok", 1.0)]),
        # a supply-constrained launch: any unmet demand in the first six months on the market
        # (lost orders surface one tick after they are placed, hence age <= 6.5)
        rule("k_short_launch", AND(gt(S("k_lost"), 0.5), lt(S("k_age"), 6.5)),
             [setm("k_short", 1.0), setp(k_share=0.85)]),
        rule("t_short_launch", AND(gt(S("t_lost"), 0.5), lt(S("t_age"), 6.5)), [setm("t_short", 1.0)]),
        rule("k_short_ever", gt(S("k_lost"), 0.5), [setm("k_short_any", 1.0)]),
        # second CMO contract (turn 5): once both readouts are in and both succeeded, sign the next
        # opening (not before mid-2028, t = 21); both assets reach it by PAS
        rule("cmo2_sign", AND(on(P("second_cmo")), on(P("k_ok")), on(P("t_ok")), ge(T, mx(T_READ, 20.0))),
             [setm("cmo2_active", 1.0), setp(cmo2_start=add(T, 1.0))]),
        rule("cmo2_ready_r", AND(on(S("cmo2_active")), ge(T, add(P("cmo2_start"), P("cmo2_lead")))),
             [setm("cmo2", CMO_TENTHS), setm("cmo2_ready", 1.0), setp(cmo2_ready_t=T)]),
        rule("k_cmo2_pas", AND(on(S("cmo2_ready")), ge(T, add(mx(P("cmo2_ready_t"), K_APPROVE), 6.0))), [setm("k_cmo2_ok", 1.0)]),
        rule("t_cmo2_pas", AND(on(S("cmo2_ready")), on(S("t_filed")), ge(T, add(mx(P("cmo2_ready_t"), T_APPROVE), 6.0))),
             [setm("t_cmo2_ok", 1.0)]),
    ]

    rules.append(rule("planner", ge(T, 0.0), plan, mode="every_tick"))   # LAST: sees this tick's final state
    params = [{"name": k, "value": float(v)} for k, v in prm.items()]
    doc = {
        "rd_format": "reactive-dynamics-model", "version": 1.0,
        "meta": {"tspan": 75.0, "dt": 1.0},
        "params": params, "places": places, "transitions": trs, "arcs": arcs, "rules": rules,
    }
    with open(os.path.join(HERE, name), "w") as f:
        json.dump(doc, f, indent=1)

ARMS = {
    # the options (base: turn-5 allocation, 2-month floor, second contract if both succeed)
    "model.rdj.json": {},                                                   # wait
    "arm_cmo.rdj.json": {"cmo_now": 1.0},                                   # reserve now, usable by both
    "arm_cmo_tern_only.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0},    # reserve now for Tern
    "arm_build_now.rdj.json": {"build_now": 1.0},
    "arm_build_after.rdj.json": {"build_after": 1.0},
    # sensitivities
    "model_literal.rdj.json": {"k_first": 1.0},
    "arm_cmo_literal.rdj.json": {"cmo_now": 1.0, "k_first": 1.0},
    "arm_cmo_tern_only_literal.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0, "k_first": 1.0},
    "model_no2.rdj.json": {"second_cmo": 0.0},
    "arm_cmo_no2.rdj.json": {"cmo_now": 1.0, "second_cmo": 0.0},
    "arm_cmo_tern_only_no2.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0, "second_cmo": 0.0},
    "model_pk50.rdj.json": {"p_k": 0.50},
    "arm_cmo_pk50.rdj.json": {"cmo_now": 1.0, "p_k": 0.50},
    "arm_cmo_tern_only_pk50.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0, "p_k": 0.50},
    "model_peak.rdj.json": {"k_val": VAL["peak"][0], "t_val": VAL["peak"][1]},
    "arm_cmo_peak.rdj.json": {"cmo_now": 1.0, "k_val": VAL["peak"][0], "t_val": VAL["peak"][1]},
    "arm_cmo_tern_only_peak.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0, "k_val": VAL["peak"][0], "t_val": VAL["peak"][1]},
    "model_nofloor.rdj.json": {"k_floor": 0.0},
    "arm_cmo_nofloor.rdj.json": {"cmo_now": 1.0, "k_floor": 0.0},
    # cost-only variants: same trajectories, different cash timing
    "arm_cmo_walk5.rdj.json": {"cmo_now": 1.0, "walk2": 0.0},
    "arm_cmo_tern_only_walk5.rdj.json": {"cmo_now": 1.0, "cmo_reassign": 0.0, "walk2": 0.0},
    "arm_build_now_even.rdj.json": {"build_now": 1.0, "capex_even": 1.0},
}

if __name__ == "__main__":
    for name, ov in ARMS.items():
        build(ov, name)
    print("wrote", len(ARMS), "documents; value per tenth-batch:", VAL)
