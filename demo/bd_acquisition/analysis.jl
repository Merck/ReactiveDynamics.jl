# BD acquisition-impact demo — analysis / post-processing (MVP §5).
#
# Outcome metrics computed from the final token population + ledger streams: portfolio rNPV
# (the headline), # launches, P(launch), the injection-robust financing dip + pool troughs, and the
# treatment effect Δ-rNPV = rNPV(Sₖ) − rNPV(S0), ensemble-averaged (MVP §4.1 — single RNG
# stream kept, per-entity substreams are §4.6 future work). Discounting/roll-ups are pure
# post-processing; the engine does no discounting (MVP §5 / finding D).
#
# COLLAPSED onto the engine APIs (Phase-0.6, ADR 0013 §14.2). The ensemble fan-out, the cross-run
# `summarize`, the unpaired `treatment_effect`, and the per-member seeding `hash((root_seed,k))`
# now live in `src/analysis.jl` (`RD.ensemble`/`RD.summarize`/`RD.treatment_effect`, with
# `EnsembleProblem` an AA-readable result node). This file keeps ONLY what is genuinely demo
# modeling: the rNPV discounting roll-up (`portfolio_rnpv`), the per-member scalar metrics
# (`run_metrics`), the thin scenario adapters over the engine ensemble (the `acq_price` netting the
# engine's symmetric `treatment_effect` does not model), and the per-program ledger views.

using ReactiveDynamics
using ReactiveDynamics: inners, getagent, get_species, find_index
using ReactiveDynamics: program_ledger, program_ledger_entries
using Statistics
using DataFrames

const RD = ReactiveDynamics

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
function portfolio_rnpv(prob; discount = 0.1, acq_price = 0.0)
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

# Per-run scalar metrics, as engine-ensemble METRIC closures (each `prob -> Real`, the shape
# `RD.summarize`/`RD.treatment_effect` consume). The headline rNPV nets the acquisition price (the
# deal scenarios carry it, S0 does not); the rest read the pools / final population.
#
# Peak capital requirement (the financing ask) = the deepest draw-down of the budget reserve BELOW
# ITS STARTING LEVEL: budget[1] - min(budget). This is injection-robust — a resource-synergy deal
# that injects +capital (SetSpecies(:budget,+Δ)) raises max(budget), so the naive max-min would
# conflate the injected capital with the acquirer's own financing dip and report a spuriously larger
# ask for exactly the scenarios that ease it. Measuring the dip below start answers the real
# question: "how much of my own reserve did I have to burn through?" The pool troughs are reported
# instead of a utilization ratio because injecting +scientists changes the pool size, so a
# 1 - sci/max(sci) ratio is not comparable across scenarios either.
m_rnpv(acq_price) = prob -> portfolio_rnpv(prob; acq_price = acq_price)
m_launches(prob) = count(reached_market, tokens(prob))
m_cash_trough(prob) = minimum(prob.sol[!, "budget"])      # lowest the cash reserve got (0 ⇒ starved)
m_sci_trough(prob) = minimum(prob.sol[!, "scientist"])    # lowest the scientist pool got (0 ⇒ booked)
m_peak_capital(prob) = (b = prob.sol[!, "budget"]; b[1] - minimum(b))

# All headline scalars for one finished run, as a NamedTuple (used by the per-program-ledger demo
# block and any caller wanting the whole row at once). `acq_price` is netted into rNPV.
function run_metrics(prob; acq_price = 0.0)
    toks = tokens(prob)
    return (
        rnpv = portfolio_rnpv(prob; acq_price = acq_price),
        launches = count(reached_market, toks),
        active = count(is_active, toks),
        retired = count(t -> !is_active(t), toks),
        peak_capital = m_peak_capital(prob),
        cash_trough = m_cash_trough(prob),
        sci_trough = m_sci_trough(prob),
    )
end

# ── Ensemble + counterfactual (MVP §4.1), now over the engine runner ─────────────────────
# `scenario_ensemble(build_and_run; …) -> EnsembleProblem` is a thin wrapper over `RD.ensemble`:
# member k is built+run by `build_and_run(hash((root_seed,k)))` — the SAME per-member seeding (§4
# D8), now owned by the engine. `build_and_run` returns a finished `ReactionNetworkProblem`, so
# `max_t` is left unset (the closure already simulates, as `run_scenario` does). The returned
# `EnsembleProblem` holds the members + realized seeds and is itself an AA-readable node.
function scenario_ensemble(build_and_run; root_seed = 2026, nseed = 20)
    return RD.ensemble(build_and_run; nseed = nseed, root_seed = root_seed)
end

# Cross-run means/SEs read straight off the engine `summarize` (mean + sem = std/√n). One metric
# closure per headline scalar; `acq_price` flows into the rNPV metric.
mean_rnpv(ens; acq_price = 0.0) = RD.summarize(ens, m_rnpv(acq_price)).mean
sem_rnpv(ens; acq_price = 0.0) = RD.summarize(ens, m_rnpv(acq_price)).sem
mean_launches(ens) = RD.summarize(ens, m_launches).mean
sem_launches(ens) = RD.summarize(ens, m_launches).sem
p_launch(ens) = RD.summarize(ens, prob -> m_launches(prob) >= 1 ? 1.0 : 0.0).mean
mean_peak_capital(ens) = RD.summarize(ens, m_peak_capital).mean
mean_cash_trough(ens) = RD.summarize(ens, m_cash_trough).mean
mean_sci_trough(ens) = RD.summarize(ens, m_sci_trough).mean

# Per-member metric vectors (the raw distributions the figures/export need) — read off the engine
# ensemble's `members` in seed order, applying the same metric closures.
rnpv_samples(ens; acq_price = 0.0) = Float64[m_rnpv(acq_price)(m) for m in ens.members]
launch_samples(ens) = Float64[m_launches(m) for m in ens.members]
cash_trough_samples(ens) = Float64[m_cash_trough(m) for m in ens.members]
sci_trough_samples(ens) = Float64[m_sci_trough(m) for m in ens.members]

# Δ-rNPV treatment effect: the engine's unpaired `treatment_effect` (difference of means with
# se = √(var_b/n_b + var_d/n_d), MVP §4.1 finding A — the two scenarios desync the shared RNG after
# the deal). The engine applies ONE metric to both arms; the deal nets the acquisition price while
# the baseline does not, so we net it here: the deal price is a deterministic constant offset, so
# Δ_net = Δ_gross − deal_price and the SE is unchanged. Returns the same fields the demo reported.
function treatment_effect(baseline_ens, deal_ens; deal_price = 0.0)
    te = RD.treatment_effect(baseline_ens, deal_ens, m_rnpv(0.0))    # gross (price-free) Δ + SE
    tl = RD.treatment_effect(baseline_ens, deal_ens, m_launches)
    pl = RD.treatment_effect(
        baseline_ens, deal_ens,
        prob -> m_launches(prob) >= 1 ? 1.0 : 0.0
    )
    return (
        delta_rnpv = te.delta - deal_price,                          # net of the acquisition price
        se_delta_rnpv = te.se,                                       # SE invariant to a constant offset
        delta_launches = tl.delta,
        delta_p_launch = pl.delta,
        baseline_rnpv = te.baseline,                                 # gross baseline (S0 price is 0)
        deal_rnpv = te.deal - deal_price,                            # net deal mean
    )
end

# ── Engine-level per-program ledger (MVP finding D — src/ledger.jl) ──────────────────────
# Historically this demo RECONSTRUCTED per-program economics in post (portfolio_rnpv walks the
# final token population and reads host fields). Finding D asked the ENGINE to attribute the ledger
# per program DURING the run. It now does: with `placeCost` priced on `budget` (host.jl), each
# program's capital burn is accrued onto it at every advance it sits in, and `program_ledger(prob)`
# returns the per-program cost/reward/valuation summary in deterministic (species, creation_index)
# order. These functions surface that engine ledger and CROSS-CHECK it against the aggregate row and
# the post-hoc rNPV roll-up.

# The total capital burned over the run, from the AGGREGATE ledger (the sum of the per-tick
# :valuation_cost rows the engine pushes in evolve!). This is the pool-level spend.
total_engine_cost(prob) = sum(r[3] for r in prob.log if r[1] == :valuation_cost; init = 0.0)

# The per-program ledger as a DataFrame, joined with each live token's current phase / npv_peak /
# acquired flag so the engine-attributed `cost_incurred` sits next to the program's modeling
# descriptors. Programs that never bound a costed transition show cost_incurred 0.
function program_economics(prob)
    led = program_ledger(prob)                       # engine ledger: program, species, cost, ...
    toks = Dict(ReactiveDynamics.AlgebraicAgents.getname(t) => t for t in tokens(prob))
    led.phase = [haskey(toks, n) ? toks[n].phase : :removed for n in led.program]
    led.npv_peak = [haskey(toks, n) ? toks[n].npv_peak : NaN for n in led.program]
    led.acquired = [haskey(toks, n) ? toks[n].acquired : false for n in led.program]
    led.reached_market = led.phase .== :Market
    return led
end

# Reconciliation: the per-program `cost_incurred` summed over ALL programs, PLUS the network
# UNATTRIBUTED bucket (capital burned by advance instances that bound no program — the documented
# finding-D boundary), equals the aggregate :valuation_cost total. Returns the three numbers + the
# residual, so a caller (or test) can assert the engine ledger is internally consistent.
function ledger_reconciliation(prob)
    led = program_ledger(prob)
    per_program = sum(led.cost_incurred)
    unattributed = prob.unattributed_cost
    aggregate = total_engine_cost(prob)
    return (
        per_program = per_program,
        unattributed = unattributed,
        attributed_plus_unattributed = per_program + unattributed,
        aggregate = aggregate,
        residual = (per_program + unattributed) - aggregate,    # ≈ 0 (the SUM invariant, §8.5)
    )
end

# A compact engine-ledger view for one finished run: total capital, how much landed on programs
# that REACHED MARKET (the value-creating spend) vs in-flight vs retired, and the reconciliation
# residual. This is the engine-level analogue of the demo's headline, now sourced from the ledger
# the engine built during the run rather than reconstructed in post.
function program_ledger_summary(prob)
    led = program_economics(prob)
    rec = ledger_reconciliation(prob)
    launched_cost = sum(led.cost_incurred[led.reached_market]; init = 0.0)
    active_cost = sum(led.cost_incurred[led.species .!= :removed .&& .!led.reached_market]; init = 0.0)
    retired_cost = sum(led.cost_incurred[led.species .== :removed]; init = 0.0)
    return (
        n_programs = nrow(led),
        total_program_cost = rec.per_program,
        launched_cost = launched_cost,           # capital that reached a launched program
        active_cost = active_cost,               # capital sunk into still-in-flight programs
        retired_cost = retired_cost,             # capital sunk into programs that failed/timed out
        unattributed_cost = rec.unattributed,    # pool-level burn with no bound program (boundary)
        aggregate_cost = rec.aggregate,
        reconciliation_residual = rec.residual,  # ≈ 0 — the per-program rows reconcile to the aggregate
    )
end
