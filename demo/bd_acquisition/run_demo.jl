# BD acquisition-impact demo — driver (MVP §7 Milestone-1).
#
# Wires the host model + the endogenous acquisition lever + the rNPV post-processor into a
# reproducible counterfactual: run the same pipeline once WITHOUT a deal (S0) and once WITH it
# (S1–S5, toggling synergies), over an ensemble of seeds, and read off the deal's attributable
# Δ-rNPV / Δ-launches. The lever lives IN the model (an ADR-0010 Rule); the driver only sets
# scenario params + seed, so each cell is one reproducible (model, scenario, seed) triple.
#
# Run:  julia --project=. demo/bd_acquisition/run_demo.jl

using Printf

const HERE = @__DIR__
include(joinpath(HERE, "host.jl"))
include(joinpath(HERE, "analysis.jl"))

# Build + run one scenario to completion under a given seed. A scenario is a choice of which
# synergies the acquisition rule arms (MVP §2.1); S0 arms no rule at all (no deal).
function run_scenario(scenario::Symbol, seed; tspan = 40.0, T_acq = 8.0)
    acs = build_pipeline_model()
    # The starting portfolio is the declarative initial marking (ADR 0007 §B) — passed to the
    # constructor, instantiated before t=0, reproducible as part of (model, population, seed).
    prob = ReactionNetworkProblem(
        acs;
        tspan = tspan,
        dt = 1.0,
        seed = seed,
        registry = PROJECT_REGISTRY,
        population = initial_population(),
    )
    if scenario != :S0
        # synergy toggles per scenario (MVP §4 grid)
        res = scenario in (:S2, :S5)           # resource synergy: +scientists
        pos = scenario in (:S3, :S5)           # capability/PoS synergy
        eff = scenario in (:S4, :S5)           # operational-efficiency synergy
        push!(
            prob.rules,
            acquisition_rule(;
                T_acq = T_acq,
                n_programs = 3,
                extra_scientists = res ? 15 : 0,
                synergy_pos = pos,
                synergy_eff = eff,
            ),
        )
    end
    simulate(prob)
    return prob
end

# The acquisition price (what the deal costs) — netted out of the deal scenarios' rNPV.
# (A price sweep finds breakeven, MVP §5; this default is below breakeven so the full-synergy
# deal is value-accretive and the synergy decomposition is visible.)
const ACQ_PRICE = 400.0

function main(; root_seed = 2026, nseed = 24)
    println("=" ^ 78)
    println("BD ACQUISITION-IMPACT DEMO — acquisition effect on a pharma pipeline portfolio")
    println("  ensemble: $nseed seeds from root $root_seed  |  horizon 40 ticks  |  acq price $(ACQ_PRICE)")
    println("=" ^ 78)

    scenarios = [:S0, :S1, :S2, :S3, :S4, :S5]
    labels = Dict(
        :S0 => "Baseline (no deal)",
        :S1 => "Deal, pipeline-only",
        :S2 => "+ resource synergy",
        :S3 => "+ capability/PoS synergy",
        :S4 => "+ op-efficiency synergy",
        :S5 => "Full (all synergies)",
    )

    results = Dict{Symbol,Any}()
    for s in scenarios
        price = s == :S0 ? 0.0 : ACQ_PRICE
        ms = ensemble(seed -> run_scenario(s, seed); root_seed = root_seed, nseed = nseed, acq_price = price)
        results[s] = ms
    end

    base = results[:S0]
    println("\nScenario                      mean rNPV    mean launches   P(≥1 launch)   Δ-rNPV vs S0")
    println("-" ^ 92)
    for s in scenarios
        ms = results[s]
        te = treatment_effect(base, ms)
        @printf(
            "%-28s  %9.1f    %9.2f       %8.2f      %+10.1f\n",
            labels[s],
            mean_rnpv(ms),
            mean_launches(ms),
            p_launch(ms),
            s == :S0 ? 0.0 : te.delta_rnpv,
        )
    end

    println("\nHeadline (Full deal S5 vs Baseline S0):")
    te = treatment_effect(base, results[:S5])
    @printf("  Δ-rNPV (deal value, net of price)  : %+.1f\n", te.delta_rnpv)
    @printf("  Δ-launches (extra programs to mkt) : %+.2f\n", te.delta_launches)
    @printf("  Δ-P(launch)                        : %+.2f\n", te.delta_p_launch)
    println("\nReproducible: each cell is a (model, scenario, seed) triple; the lever is in-model")
    println("(an ADR-0010 Rule), so re-running with the same root seed gives identical numbers.")
    return results
end

# Run when executed as a script (not when included by a test).
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
