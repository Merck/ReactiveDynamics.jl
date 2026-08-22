# RDPlotsExt — result-plot recipes + the generic `_draw` reduction (ADR 0014 §A / CONTRACT §15.1).
#
# Loaded automatically when `Plots` is available (declared in Project.toml `[weakdeps]`/`[extensions]`).
# Holds the model-agnostic Plots.jl `@recipe` set, each keyed off a thin wrapper type defined in the
# core (src/analysis.jl) over a raw result artifact (`prob.sol`, `prob.log`, `program_ledger(prob)`,
# the ADR 0013 `token_trajectory`/`ensemble`) so the recipes work on ANY model (§15.4 Invariant 5),
# plus the live `_draw` species-trajectory reduction relocated from src/interface/plots.jl.
module RDPlotsExt

using ReactiveDynamics
using ReactiveDynamics: MarkingPlot, SaturationPlot, ValuationPlot, LedgerPlot,
    TokenTrajectoryPlot, EnsembleBar, TreatmentEffectPlot, ThroughputPlot,
    log_scalar_series, log_count_series, token_trajectory, trajectory_envelope,
    summarize, ReactionNetworkProblem
using ReactiveDynamics: DataFrames
using Plots
using Statistics

# ── Recipe 1: marking trajectory (species/token counts over time) ───────────────────────
@recipe function f(m::MarkingPlot)
    xguide --> "time"
    yguide --> "quantity"
    sol = m.prob.sol
    for var in m.vars
        @series begin
            label --> var
            sol[!, "t"], sol[!, var]
        end
    end
end

# ── Recipe 2: resource utilization / saturation (pool levels; troughs = where it binds) ──
@recipe function f(s::SaturationPlot)
    xguide --> "time"
    yguide --> "pool level"
    title --> "resource utilization (troughs = binding constraint)"
    sol = s.prob.sol
    for var in s.vars
        @series begin
            label --> var
            sol[!, "t"], sol[!, var]
        end
    end
end

# ── Recipe 3: valuation / rNPV curve (portfolio value, cost, reward over time) ───────────
@recipe function f(v::ValuationPlot)
    xguide --> "time"
    yguide --> "value"
    title --> "valuation"
    for (tag, lbl) in (
            (:valuation, "portfolio valuation"),
            (:valuation_cost, "cost"),
            (:valuation_reward, "reward"),
        )
        t, y = log_scalar_series(v.prob, tag)
        isempty(t) && continue
        @series begin
            label --> lbl
            t, y
        end
    end
end

# ── Recipe 4: per-program ledger (cost/reward bars per program) ──────────────────────────
@recipe function f(l::LedgerPlot)
    led = ReactiveDynamics.program_ledger(l.prob)
    seriestype --> :bar
    xguide --> "program"
    yguide --> "value"
    title --> "per-program ledger"
    xs = led.program
    @series begin
        label --> "cost_incurred"
        xs, led.cost_incurred
    end
    @series begin
        label --> "reward_realized"
        xs, led.reward_realized
    end
end

# ── Recipe 5: token trajectory + the typical envelope band ───────────────────────────────
@recipe function f(tt::TokenTrajectoryPlot)
    df = tt.pred === nothing ? token_trajectory(tt.prob) : token_trajectory(tt.prob, tt.pred)
    field = tt.field
    xguide --> "time"
    yguide --> string(field)
    title --> "token trajectory — $(field)"
    # individual token paths (thin, no legend spam)
    if string(field) in names(df)
        for g in DataFrames.groupby(df, :program)
            @series begin
                label := ""
                seriesalpha --> 0.35
                linewidth --> 1
                g[!, "t"], g[!, string(field)]
            end
        end
    end
    # the median + IQR envelope band over the cohort
    env = tt.pred === nothing ? trajectory_envelope(tt.prob) : trajectory_envelope(tt.prob, tt.pred)
    fenv = env[env.field .== field, :]
    if DataFrames.nrow(fenv) > 0
        @series begin
            label --> "median ± IQR"
            linewidth --> 3
            ribbon --> (fenv.median .- fenv.q25, fenv.q75 .- fenv.median)
            fenv.align, fenv.median
        end
    end
end

# ── Recipe 6a: ensemble distribution (per-run metric histogram across members) ───────────
@recipe function f(eb::EnsembleBar)
    xs = Float64[Float64(eb.metric(m)) for m in eb.ens.members]
    seriestype --> :histogram
    xguide --> "per-run metric"
    yguide --> "members"
    title --> "ensemble distribution (n=$(length(xs)))"
    label --> "metric"
    xs
end

# ── Recipe 6b: treatment effect (baseline vs deal distributions, Δ annotated) ────────────
@recipe function f(te::TreatmentEffectPlot)
    b = Float64[Float64(te.metric(m)) for m in te.baseline.members]
    d = Float64[Float64(te.metric(m)) for m in te.deal.members]
    Δ = mean(d) - mean(b)
    seriestype --> :histogram
    fillalpha --> 0.5
    xguide --> "per-run metric"
    yguide --> "members"
    title --> "treatment effect  Δ = $(round(Δ, digits = 1))"
    @series begin
        label --> "baseline (mean $(round(mean(b), digits = 1)))"
        b
    end
    @series begin
        label --> "deal (mean $(round(mean(d), digits = 1)))"
        d
    end
end

# ── Recipe 7: throughput (firings / terminations per tick) ───────────────────────────────
@recipe function f(tp::ThroughputPlot)
    xguide --> "time"
    yguide --> "instances / tick"
    title --> "throughput"
    for (tag, lbl) in (
            (:new_transitions, "firings"),
            (:terminated_all, "terminations"),
            (:terminated_success, "successful terminations"),
        )
        t, y = log_count_series(tp.prob, tag)
        isempty(t) && continue
        @series begin
            label --> lbl
            t, y
        end
    end
end

# ── The live generic `_draw` reduction, relocated from src/interface/plots.jl ────────────
# Species trajectories from `prob.sol` (generalized by recipe 1, kept as the AA `draw(prob)` entry).
function ReactiveDynamics.AlgebraicAgents._draw(
        prob::ReactionNetworkProblem,
        vars = string.(prob.acs[:, :placeName]);
        kwargs...,
    )
    p = Plots.plot()
    for var in vars
        p = Plots.plot!(
            p, prob.sol[!, "t"], prob.sol[!, var];
            label = "$var", xlabel = "time", ylabel = "quantity", kwargs...
        )
    end
    return p
end

end # module RDPlotsExt
