#!/usr/bin/env julia
#
# figures.jl — the study's five figures, rebuilt from the rdx envelopes in `results/`.
#
#   julia --project=. agent/studies/portfolio_bid/figures.jl
#
# Nothing here computes a simulation result. Every number is read out of an envelope written by a
# command recorded in `README.md`, so a figure cannot drift away from the number it illustrates: move
# the study and the figures move with it, or the read fails loudly. The fifth figure (`exec_map.svg`)
# is Graphviz output from `rdx draw` and is not rebuilt here.
#
# One derived quantity appears in the bridge — the stranded value — and it is a linear combination of
# four INDEPENDENT arm means, so its se is the root-sum-square of the arm sems. That is the same
# unpaired estimator `RD.treatment_effect` uses (`src/analysis.jl:401`: the arms desync their RNG
# streams after the first divergent draw, so there is no pairing left to exploit).

using JSON

const HERE = @__DIR__
const RESULTS = joinpath(HERE, "results")
const FIGDIR = joinpath(HERE, "figures")

include(joinpath(HERE, "..", "..", "scripts", "chartkit.jl"))

mkpath(FIGDIR)

envelope(name) = JSON.parsefile(joinpath(RESULTS, name * ".json"))["result"]

# `effect` → (delta, se) for one metric; `ensemble` → (mean, sem); `run --trace` → the series.
eff(name, metric) = let m = envelope("effect_" * name)["metrics"][metric]
    (m["delta"], m["se"])
end
arm(name, metric) = let m = envelope("ens_" * name)["metrics"][metric]
    (m["mean"], m["sem"])
end
samples(name, metric) = envelope("ens_" * name)["metrics"][metric]["samples"]
trace(name) = envelope("trace_" * name)["trace"]

rss(xs) = sqrt(sum(abs2, xs))

write_fig(name, svg) = begin
    path = joinpath(FIGDIR, name * ".svg")
    write(path, svg)
    println("figures: wrote $name.svg — $(round(sizeof(svg) / 1024; digits = 1)) KiB")
end

# ── 1. the ceiling, and what moving it does ───────────────────────────────────────────────────────
# Committed FTE is the ceiling minus the free pool. Two arms share one axis so the second line can be
# seen crossing the line the board drew — which is what funding twelve more people looks like.

const FTE_BASE = 40
const FTE_BOTH = 52

let acq = trace("arm_acquire"), both = trace("arm_both")
    months = acq["t"]
    committed(tr, cap) = [cap - v for v in tr["fte"]]
    write_fig(
        "capacity_trace",
        ChartKit.trace_ceiling(
            months,
            [
                Dict(
                    "label" => "Bid only — 40 FTE", "emphasis" => true,
                    "y" => committed(acq, FTE_BASE),
                ),
                Dict(
                    "label" => "Bid + 12 FTE funded", "emphasis" => false,
                    "y" => committed(both, FTE_BOTH),
                ),
            ];
            ceiling = FTE_BASE,
            ceiling_label = "40 FTE — today's ceiling",
            shade_at = 0.975, unit = "FTE", digits = 0,
            x_label = "Month", y_label = "FTE committed",
            title = "What the ceiling does to the work you just bought",
            subtitle = "Development FTE committed to in-flight programmes, seed 1. Shaded where the " *
                "bid-only arm sits within one FTE of the dashed ceiling; the grey arm has 52.",
        )
    )
end

# ── 2. every lever, with its error bar, nulls drawn as nulls ───────────────────────────────────────

let
    levers = [
        ("Bid + fund 12 FTE", "both", "vs. do nothing"),
        ("Fund 12 FTE, given the bid", "integration", "the integration budget"),
        ("Fund 24 FTE, given the bid", "integration_24", "no better than 12"),
        ("Bid, capacity unchanged", "acquire", "what the assets are worth to you"),
        ("Fund 4 FTE, no bid", "capacity_4", "inside the noise"),
        ("Fund 12 FTE, no bid", "capacity", "inside the noise"),
    ]
    series = [
        let (d, se) = eff(key, "enpv")
            Dict("label" => label, "value" => d, "se" => se, "note" => note)
        end for (label, key, note) in levers
    ]
    write_fig(
        "levers",
        ChartKit.bars_error(
            series; unit = "\$M", digits = 0, sortby = "desc", null_k = 2.0,
            x_label = "Δ eNPV over ten years, \$M",
            title = "Six ways to spend the money, ranked",
            subtitle = "120 seeds per arm, root seed 2026, discounted 8%/yr on monthly ticks. " *
                "Bars inside ±2 se are drawn hollow: those are nulls, not small effects.",
        )
    )
end

# ── 3. the spread the error bars summarise ────────────────────────────────────────────────────────
# The figure that earns the rest of them: one simulated decade tells you almost nothing, because a
# launch is a coin flip worth several hundred million. The mean of 120 is what carries the decision.

let
    b, bse = arm("model", "enpv")
    v, vse = arm("arm_both", "enpv")
    write_fig(
        "spread",
        ChartKit.paired_dots(
            [
                Dict("label" => "Do nothing", "samples" => samples("model", "enpv")),
                Dict("label" => "Bid + 12 FTE", "samples" => samples("arm_both", "enpv")),
            ];
            delta = v - b, delta_se = rss((bse, vse)),
            unit = "\$M", digits = 0, x_label = "eNPV of the decade, \$M",
            title = "Why one scenario is not an answer",
            subtitle = "Every dot is one simulated decade. The decades overlap almost completely; " *
                "the difference of their means does not.",
        )
    )
end

# ── 4. the bridge from the seller's number to yours ────────────────────────────────────────────────
# A reconciliation, not an inference: each quantity on it was measured separately, and the two steps
# are the two measured effects that connect the three measured levels.

let
    stand, sstand = eff("standalone", "enpv")
    real, sreal = eff("acquire", "enpv")
    integ, sinteg = eff("integration", "enpv")
    prog, sprog = eff("both", "enpv")
    stranded = real - stand
    write_fig(
        "bridge",
        ChartKit.waterfall(
            Dict("label" => "Standalone value\nof the three assets", "value" => stand, "se" => sstand),
            [
                Dict(
                    "label" => "Stranded by the\n40-FTE ceiling", "value" => stranded,
                    "se" => rss((sstand, sreal))
                ),
                Dict(
                    "label" => "Recovered by funding\n12 more FTE", "value" => integ,
                    "se" => sinteg
                ),
            ],
            Dict("label" => "Value to you,\ncapacity funded", "value" => prog, "se" => sprog);
            unit = "\$M", digits = 0, y_label = "Δ eNPV over ten years, \$M",
            title = "What the seller's valuation assumes, and what you would actually get",
            subtitle = "Standalone is the same three assets against unlimited capacity. Each level is " *
                "its own 120-seed measurement and the steps are their differences; thick bar ±1 se, " *
                "thin ±2 se.",
        )
    )
end

println("figures: done — 4 rebuilt here, exec_map.svg comes from `rdx draw --exec-map`")
