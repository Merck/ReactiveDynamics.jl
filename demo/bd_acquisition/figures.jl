# BD acquisition-impact demo — presentation figures (the "technical-executive" deck layer).
#
# Generates three annotated PNGs from the SAME engine run as run_demo.jl (no separate model):
#   1. delta_waterfall.png   — Δ-rNPV bridge S0 → pipeline → +each synergy → Full (the answer)
#   2. synergy_decomp.png    — marginal synergy contributions with ±1 SE bars (the insight)
#   3. rnpv_distribution.png  — per-seed rNPV spread S0 vs S5 (honest stochastic view)
#   4. capital_tradeoff.png   — the binding-cash view: financing dip + pool low-water marks
#
# Run:  julia --project=. demo/bd_acquisition/figures.jl   (writes to demo/bd_acquisition/figures/)

using Printf, Statistics
using Plots
gr()

const HERE = @__DIR__
include(joinpath(HERE, "host.jl"))
include(joinpath(HERE, "analysis.jl"))
include(joinpath(HERE, "run_demo.jl"))

const OUTDIR = joinpath(HERE, "figures")
isdir(OUTDIR) || mkdir(OUTDIR)

# Brand-ish palette (muted, deck-friendly).
const NAVY = RGB(0.13, 0.20, 0.36)
const TEAL = RGB(0.10, 0.52, 0.55)
const GOLD = RGB(0.83, 0.62, 0.18)
const GREY = RGB(0.62, 0.65, 0.70)
const RUST = RGB(0.72, 0.34, 0.22)

# ── Run the full grid once and cache per-member metric vectors ──────────────────────────
function run_grid(; root_seed = 2026, nseed = 160)
    scenarios = [:S0, :S1, :S2, :S3, :S4, :S5]
    Dict(s => ensemble(seed -> run_scenario(s, seed); root_seed = root_seed, nseed = nseed,
                       acq_price = (s == :S0 ? 0.0 : ACQ_PRICE)) for s in scenarios)
end

function make_figures(; nseed = 160)
    @info "Running S0–S5 grid for figures ($nseed seeds)…"
    R = run_grid(; nseed = nseed)
    base = R[:S0]
    Δ(s)   = treatment_effect(base, R[s]).delta_rnpv
    SEΔ(s) = treatment_effect(base, R[s]).se_delta_rnpv

    # ── Figure 1: Δ-rNPV waterfall (bridge) ──────────────────────────────────────────────
    # Steps: baseline level → + pipeline → + resource → + capability → + op-eff → +interaction → Full.
    # Drawn with explicit rectangles (Shape) so each bar genuinely FLOATS from its running total —
    # a true waterfall — rather than all starting at zero.
    s1 = Δ(:S1)
    marg_res = Δ(:S2) - s1
    marg_pos = Δ(:S3) - s1
    marg_eff = Δ(:S4) - s1
    interaction = Δ(:S5) - (s1 + marg_res + marg_pos + marg_eff)  # super/sub-additivity remainder
    steps   = ["Baseline\n(S0)", "+ Pipeline\nprograms", "+ Resource\nsynergy",
               "+ Capability\n/PoS synergy", "+ Op-eff.\nsynergy", "+ Inter-\naction", "Full deal\n(S5)"]
    incr    = [s1, marg_res, marg_pos, marg_eff, interaction]      # the five additive blocks
    barcols = [TEAL, GOLD, RUST, GREY, RGB(0.55,0.58,0.63)]
    base0   = mean_rnpv(base)
    total   = base0 + sum(incr)
    n       = length(steps)
    hw      = 0.38                                                  # half-width of each bar

    rect(i, lo, hi) = Shape([i-hw, i+hw, i+hw, i-hw], [lo, lo, hi, hi])
    p1 = plot(; legend = false, title = "Where the deal's value comes from — Δ-rNPV bridge",
        ylabel = "risk-adjusted portfolio NPV", size = (980, 560), titlefontsize = 13,
        guidefontsize = 10, tickfontsize = 8, grid = :y, framestyle = :box,
        xticks = (1:n, steps), xlims = (0.4, n+0.6), bottom_margin = 7Plots.mm, left_margin = 7Plots.mm)
    # endpoint totals (bar 1 = baseline level, bar n = full-deal level), both float from 0
    plot!(p1, rect(1, 0, base0); fillcolor = NAVY, linecolor = :white, fillalpha = 0.92)
    plot!(p1, rect(n, 0, total); fillcolor = NAVY, linecolor = :white, fillalpha = 0.92)
    annotate!(p1, 1, base0 + 110, text(@sprintf("%.0f", base0), 8, :black))
    annotate!(p1, n, total + 110, text(@sprintf("%.0f", total), 9, :black))
    # floating marginal blocks + connector lines
    runtot = base0
    for (k, dv) in enumerate(incr)
        i = k + 1
        lo, hi = runtot, runtot + dv
        plot!(p1, rect(i, min(lo,hi), max(lo,hi)); fillcolor = barcols[k], linecolor = :white, fillalpha = 0.92)
        plot!(p1, [i-1+hw, i-hw], [runtot, runtot]; color = GREY, lw = 1, linestyle = :dot)  # connector
        annotate!(p1, i, max(lo,hi) + 110, text(@sprintf("%+.0f", dv), 8, :black))
        runtot = hi
    end
    plot!(p1, [n-1+hw, n-hw], [runtot, runtot]; color = GREY, lw = 1, linestyle = :dot)
    savefig(p1, joinpath(OUTDIR, "delta_waterfall.png"))

    # ── Figure 2: synergy decomposition with ±1 SE bars ─────────────────────────────────
    # Marginal Δ of each synergy over the pipeline-only deal (S1), plus the full deal Δ over S0.
    labels2 = ["Pipeline\nonly (S1)", "Resource\n(capital+HC)", "Capability\n/PoS", "Op-efficiency", "Full deal\n(S5)"]
    # values: S1 over S0; marginals over S1; S5 over S0
    vals2 = [s1, marg_res, marg_pos, marg_eff, Δ(:S5)]
    # SEs: S1 and S5 are vs-S0 SEs; marginals use the paired SE over S1
    paired_se(s) = (d = [m.rnpv for m in R[s]] .- [m.rnpv for m in R[:S1]]; std(d)/sqrt(length(d)))
    ses2 = [SEΔ(:S1), paired_se(:S2), paired_se(:S3), paired_se(:S4), SEΔ(:S5)]
    cols2 = [TEAL, GOLD, RUST, GREY, NAVY]
    n2 = length(labels2)
    p2 = plot(; legend = false, title = "Synergy decomposition (Δ-rNPV, ±1 SE)",
        ylabel = "Δ risk-adjusted NPV (S1 & S5 vs baseline; synergies vs S1)", size = (960, 560),
        titlefontsize = 13, guidefontsize = 9, tickfontsize = 9, grid = :y, framestyle = :box,
        xticks = (1:n2, labels2), xlims = (0.4, n2 + 0.6), ylims = (0, maximum(vals2) * 1.18),
        bottom_margin = 7Plots.mm, left_margin = 7Plots.mm)
    rect2(i, h) = Shape([i-0.4, i+0.4, i+0.4, i-0.4], [0, 0, h, h])
    for (i, v) in enumerate(vals2)
        plot!(p2, rect2(i, v); fillcolor = cols2[i], linecolor = :white, fillalpha = 0.92)
        plot!(p2, [i, i], [v - ses2[i], v + ses2[i]]; color = :black, lw = 1.5)             # error bar
        plot!(p2, [i-0.08, i+0.08], [v + ses2[i], v + ses2[i]]; color = :black, lw = 1.5)    # top cap
        plot!(p2, [i-0.08, i+0.08], [v - ses2[i], v - ses2[i]]; color = :black, lw = 1.5)    # bottom cap
        annotate!(p2, i, v + ses2[i] + 55, text(@sprintf("%+.0f", v), 9, :black))
    end
    annotate!(p2, 3, vals2[3] + ses2[3] + 200, text("← the BD insight:\nvalue lifts the pipeline\nyou ALREADY own", 8, RUST, :center))
    savefig(p2, joinpath(OUTDIR, "synergy_decomp.png"))

    # ── Figure 3: rNPV distribution S0 vs S5 (honest stochastic spread) ──────────────────
    # Overlaid histograms (base Plots; StatsPlots/violin not in the env) — the point is that the
    # deal shifts the WHOLE distribution right, not just the mean, and the arms overlap (so the
    # headline Δ is a mean shift over wide per-seed spread, consistent with the ±SE on the table).
    r0 = [m.rnpv for m in R[:S0]]
    r5 = [m.rnpv for m in R[:S5]]
    edges = range(min(minimum(r0), minimum(r5)), max(maximum(r0), maximum(r5)); length = 28)
    # headroom: tallest bin (the baseline mode) must not clip under the legend
    ymax = 1.18 * maximum(vcat(
        [count(v -> edges[b] <= v < edges[b+1], r0) for b in 1:length(edges)-1],
        [count(v -> edges[b] <= v < edges[b+1], r5) for b in 1:length(edges)-1]))
    p3 = histogram(r0; bins = edges, color = GREY, alpha = 0.55, linecolor = :white,
        label = @sprintf("Baseline S0 (mean %.0f)", mean(r0)),
        title = "Outcomes are distributions, not points (per-seed rNPV)",
        xlabel = "portfolio rNPV", ylabel = "ensemble members", size = (880, 540),
        titlefontsize = 13, guidefontsize = 10, tickfontsize = 9, framestyle = :box, ylims = (0, ymax),
        legend = :topright, left_margin = 6Plots.mm, bottom_margin = 5Plots.mm)
    histogram!(p3, r5; bins = edges, color = TEAL, alpha = 0.55, linecolor = :white,
        label = @sprintf("Full deal S5 (mean %.0f)", mean(r5)))
    vline!(p3, [mean(r0)]; color = GREY, lw = 2, linestyle = :dash, label = "")
    vline!(p3, [mean(r5)]; color = TEAL, lw = 2, linestyle = :dash, label = "")
    savefig(p3, joinpath(OUTDIR, "rnpv_distribution.png"))

    # ── Figure 4: which constraint binds, and when the deal relieves it ──────────────────
    # The two pools' LOW-WATER MARKS (independent signals — unlike the financing dip, which is just
    # 150 − cash trough). Both run near the floor through S4: the pipeline is cash- AND
    # capacity-constrained, and resource synergy alone (S2) is absorbed into more parallel programs.
    # Only the FULL deal (S5) lifts both troughs off the floor — injected capital is run, not banked,
    # until capability synergy also pushes programs out to launch and frees the pools.
    scen = [:S0, :S1, :S2, :S3, :S4, :S5]
    names4 = ["S0", "S1", "S2", "S3", "S4", "S5"]
    cashlo = [mean_cash_trough(R[s]) for s in scen]
    scilo  = [mean_sci_trough(R[s]) for s in scen]
    p4 = plot(names4, cashlo; color = GOLD, lw = 3, marker = :circle, markersize = 6,
        label = "cash trough (budget pool low-water mark)", legend = :topleft,
        ylabel = "pool units remaining at the trough",
        title = "Both constraints bind until the FULL deal lifts them off the floor",
        size = (920, 520), titlefontsize = 12, guidefontsize = 10, tickfontsize = 9,
        ylims = (0, 45), framestyle = :box, left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
    plot!(p4, names4, scilo; color = TEAL, lw = 3, marker = :diamond, markersize = 6,
        label = "scientist trough (headcount pool low-water mark)")
    hline!(p4, [0]; color = RUST, lw = 1, linestyle = :dash, label = "")
    annotate!(p4, 3.5, 26, text("pinned near the floor through S4 — resource synergy alone (S2)", 7, NAVY, :center))
    annotate!(p4, 3.5, 23.5, text("is absorbed into more parallel programs, not banked", 7, NAVY, :center))
    annotate!(p4, 5.55, scilo[6] + 3.0, text("S5 frees\nboth pools", 7, TEAL, :center))
    savefig(p4, joinpath(OUTDIR, "capital_tradeoff.png"))

    @info "Wrote 4 figures to $OUTDIR"
    for f in ("delta_waterfall.png", "synergy_decomp.png", "rnpv_distribution.png", "capital_tradeoff.png")
        println("  ", joinpath(OUTDIR, f))
    end
    return R
end

if abspath(PROGRAM_FILE) == @__FILE__
    make_figures()
end
