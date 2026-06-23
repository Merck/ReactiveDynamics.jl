# BD acquisition-impact demo — analysis / post-processing (MVP §5).
#
# Outcome metrics computed from the final token population + ledger streams: portfolio rNPV
# (the headline), # launches, P(launch), peak capital, scientist utilization, and the
# treatment effect Δ-rNPV = rNPV(Sₖ) − rNPV(S0), ensemble-averaged (MVP §4.1 — single RNG
# stream kept, per-entity substreams are §4.6 future work). Discounting/roll-ups are pure
# post-processing; the engine does no discounting (MVP §5 / finding D).

using ReactiveDynamics
using ReactiveDynamics: inners, getagent, get_species, find_index
using Statistics

# A program is RETIRED when its species (kind) has been flipped to :removed on failure/timeout
# (ADR 0006 soft-retire). Its `phase` ATTRIBUTE records how far it got (phase-as-attribute).
is_active(t) = get_species(t) != :removed
reached_market(t) = t.phase == :Market

tokens(prob) = collect(values(inners(getagent(prob, "structured"))))

# Portfolio rNPV (MVP §5 headline): expected value of the live pipeline net of the acquisition
# price. A launched program (reached :Market) has realized its peak value; an in-flight program
# is risk-adjusted by its remaining probability of success and discounted for the expected time
# still to launch (a coarse phase→years-to-market map). The engine does no discounting (MVP §5 /
# finding D) — this is the post-processing roll-up over the final token population.
const YEARS_TO_MARKET = Dict(
    :Discovery => 9.0, :Phase1 => 7.0, :Phase2 => 5.0, :Phase3 => 3.0, :Filed => 1.0, :Market => 0.0,
)
function portfolio_rnpv(prob; discount = 0.10, acq_price = 0.0)
    gross = 0.0
    for t in tokens(prob)
        is_active(t) || continue
        if reached_market(t)
            gross += t.npv_peak                                    # realized at launch
        else
            ttm = get(YEARS_TO_MARKET, t.phase, 8.0)
            gross += t.pos_remaining * t.npv_peak / (1 + discount)^ttm  # risked + time-discounted
        end
    end
    return gross - acq_price
end

# Headline scalar metrics for one finished run.
function metrics(prob; acq_price = 0.0)
    toks = tokens(prob)
    launches = count(reached_market, toks)
    active = count(is_active, toks)
    retired = count(t -> !is_active(t), toks)
    # Peak capital requirement = the largest cumulative budget burn (the financing ask). Budget is
    # a pool that financing refills; the peak draw-down below the initial level proxies the ask.
    budget = prob.sol[!, "budget"]
    peak_capital = maximum(budget) - minimum(budget)
    sci = prob.sol[!, "scientist"]
    utilization = 1 .- sci ./ maximum(sci)   # fraction of scientists in use over time
    return (
        rnpv = portfolio_rnpv(prob; acq_price = acq_price),
        launches = launches,
        active = active,
        retired = retired,
        peak_capital = peak_capital,
        mean_utilization = mean(utilization),
    )
end

# ── Ensemble + counterfactual (MVP §4.1) ────────────────────────────────────────────────
# Run a scenario over an ensemble of seeds derived from a root seed (§4 D8), returning the
# per-member metric vector. `build` is a closure (root_seed-independent) returning a finished prob.
function ensemble(build_and_run; root_seed = 2026, nseed = 20, acq_price = 0.0)
    return [metrics(build_and_run(hash((root_seed, k))); acq_price = acq_price) for k = 1:nseed]
end

mean_rnpv(ms) = mean(m.rnpv for m in ms)
mean_launches(ms) = mean(m.launches for m in ms)
p_launch(ms) = mean(m.launches >= 1 for m in ms)

# Δ-rNPV treatment effect: ensemble-averaged difference of means (NOT per-seed paired, MVP §4.1
# finding A — the two scenarios desync the shared RNG after the deal). Returns the deal's
# attributable change in portfolio value, launches, and launch probability.
function treatment_effect(baseline_ms, deal_ms)
    return (
        delta_rnpv = mean_rnpv(deal_ms) - mean_rnpv(baseline_ms),
        delta_launches = mean_launches(deal_ms) - mean_launches(baseline_ms),
        delta_p_launch = p_launch(deal_ms) - p_launch(baseline_ms),
        baseline_rnpv = mean_rnpv(baseline_ms),
        deal_rnpv = mean_rnpv(deal_ms),
    )
end
