#!/usr/bin/env julia
#
# numbers.jl — re-derive every number `README.md` quotes, straight from the envelopes in `results/`.
#
#   julia --project=. agent/studies/portfolio_bid/numbers.jl
#
# Nothing here simulates: it reads the JSON the CLI wrote and prints the rows, so a disagreement
# between this output and README.md is a bug in README.md. That is not hypothetical — the first
# hand-written version of the interaction row carried a value from a different contrast and a standard
# error copied from the line above, and so reported a clear result as a null. No validator can see a
# mistake like that, because it is in a transcribed summary and not in the model. Running this is what
# caught it.
#
# The one piece of statistics here: a derived contrast is a signed sum of arm means, so its standard
# error is the root-sum-square of the sems of the arms it combines. `rdx effect` uses exactly that for
# a pair (`RD.treatment_effect`, `src/analysis.jl:401` — the arms desync their RNG streams after the
# first divergent draw, so there is no pairing left to exploit), and the identity check below confirms
# the envelope's se against the two arm sems before any of it is trusted. The consequence worth
# stating: differencing two EFFECT ses is valid only when those effects share no arm.

using JSON, Printf

const RESULTS = joinpath(@__DIR__, "results")

env(name) = JSON.parsefile(joinpath(RESULTS, name * ".json"))["result"]
metrics(name) = env(name)["metrics"]
has(name) = isfile(joinpath(RESULTS, name * ".json"))

armmean(a, metric = "enpv") = metrics("ens_" * a)[metric]["mean"]
armsem(a, metric = "enpv") = metrics("ens_" * a)[metric]["sem"]

rss(xs) = sqrt(sum(abs2, xs))
verdict(d, s) = abs(d) > 2s ? "clears" : "NULL"

function row(label, d, s; digits = 3)
    num(x) = string(round(x; digits = digits))
    @printf("  %-38s %10s  ± %7s   %5.2f se   %s\n", label, (d > 0 ? "+" : "") * num(d), num(s), abs(d) / s, verdict(d, s))
    return nothing
end

# A contrast is a signed list of arms. `dose response` is the case that matters: written out as
# (both − acquire) − (both_24 − acquire) the shared arm cancels algebraically, which is the same thing
# as removing its variance once instead of counting it twice.
function contrast(terms, metric = "enpv")
    d = sum(c * armmean(a, metric) for (c, a) in terms)
    s = rss([armsem(a, metric) for (_, a) in terms])
    return (d, s)
end

const ARMS = [
    "model", "arm_acquire", "arm_capacity", "arm_both", "arm_both_24",
    "uncapped_baseline", "uncapped_acquire",
]

# label, effect envelope, baseline arm, variant arm
const PAIRS = [
    ("standalone", "standalone", "uncapped_baseline", "uncapped_acquire"),
    ("realized", "acquire", "model", "arm_acquire"),
    ("capacity alone, +12 FTE", "capacity", "model", "arm_capacity"),
    ("capacity alone, +4 FTE", "capacity_4", "model", "arm_capacity_4"),
    ("integration, +12 FTE given the deal", "integration", "arm_acquire", "arm_both"),
    ("integration, +24 FTE given the deal", "integration_24", "arm_acquire", "arm_both_24"),
    ("programme", "both", "model", "arm_both"),
]

println("── arm levels, eNPV, 120 members each ────────────────────────────────────────────")
for a in ARMS
    has("ens_" * a) || continue
    @printf("  %-22s mean %9.1f   sem %6.1f\n", a, armmean(a), armsem(a))
end

println("\n── paired effects, Δ eNPV ────────────────────────────────────────────────────────")
for (label, key, _, _) in PAIRS
    m = metrics("effect_" * key)["enpv"]
    row(label, m["delta"], m["se"]; digits = 1)
end

println("\n── identity: does the envelope se equal rss(arm sems)? ───────────────────────────")
for (label, key, b, v) in PAIRS
    (has("ens_" * b) && has("ens_" * v)) || continue
    envse, armse = metrics("effect_" * key)["enpv"]["se"], rss([armsem(b), armsem(v)])
    ok = abs(envse - armse) < 1.0e-3 ? "ok" : "MISMATCH"
    @printf("  %-38s %8.3f vs %8.3f   %s\n", label, envse, armse, ok)
end

println("\n── derived contrasts ─────────────────────────────────────────────────────────────")
row(
    "stranded value = standalone − realized",
    contrast([(1, "uncapped_acquire"), (-1, "uncapped_baseline"), (-1, "arm_acquire"), (1, "model")])...;
    digits = 1
)
row(
    "stranded launches, final:acq_market",
    contrast([(1, "uncapped_acquire"), (-1, "arm_acquire")], "final:acq_market")...
)
row(
    "interaction = (both−acq) − (cap−model)",
    contrast([(1, "arm_both"), (-1, "arm_acquire"), (-1, "arm_capacity"), (1, "model")])...;
    digits = 1
)
row("dose response = both − both_24", contrast([(1, "arm_both"), (-1, "arm_both_24")])...; digits = 1)

println("\n  the same dose response done WRONG, by differencing the two effect ses:")
let d = metrics("effect_integration")["enpv"], e = metrics("effect_integration_24")["enpv"]
    row("  (shares arm_acquire — counted twice)", d["delta"] - e["delta"], rss([d["se"], e["se"]]); digits = 1)
end

println("\n── mechanism: crowd-out, and the control ─────────────────────────────────────────")
for metric in ("final:own_market", "successes:adv_own_p2")
    println("  $metric")
    for (label, key) in (
            ("deal at 40 FTE (model → arm_acquire)", "acquire"),
            ("deal with no ceiling (uncapped pair)", "standalone"),
            ("deal + 12 FTE (model → arm_both)", "both"),
        )
        m = metrics("effect_" * key)[metric]
        @printf(
            "    %-42s %+7.3f ± %6.3f   %5.2f se   %s\n",
            label, m["delta"], m["se"], abs(m["delta"]) / m["se"], verdict(m["delta"], m["se"])
        )
    end
end

println("\n── queue gauge: max:own_p1, the deepest the preclinical backlog gets ──────────────")
let seen = Dict{String, Float64}()
    for (_, key, b, v) in PAIRS
        m = metrics("effect_" * key)["max:own_p1"]
        seen[b] = m["baseline"]
        seen[v] = m["variant"]
    end
    for a in ARMS
        haskey(seen, a) && @printf("  %-22s %6.1f\n", a, seen[a])
    end
end

println("\n── traces, seed 1: when the acquired assets start Phase II ────────────────────────")
for a in ("model", "arm_acquire", "arm_both")
    tr = env("trace_" * a)["trace"]
    cap = a == "arm_both" ? 52 : 40
    starts = Int[]
    prev = 0.0
    for (i, v) in enumerate(tr["acq_p2"])                  # a start consumes one from acq_p2 …
        i > 1 && v < prev - 1.0e-9 &&
            append!(starts, fill(round(Int, tr["t"][i]), round(Int, prev - v)))   # … two in one tick, twice
        prev = v
    end
    tight = count(<=(1.0), tr["fte"])
    @printf(
        "  %-14s starts at %-14s   FTE free ≤1 in %3d of %3d ticks (%2.0f%%, ceiling %d)\n",
        a, isempty(starts) ? "—" : join(starts, ", "), tight, length(tr["fte"]),
        100 * tight / length(tr["fte"]), cap
    )
end
