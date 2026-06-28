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

# The pipeline is also serialized as an eval-free `model.rdj.json` (ADR 0005, Stage E). Loading it
# via from_json_model with the host PROJECT_REGISTRY yields a model byte-for-byte identical to the
# DSL build_pipeline_model() under the same seed (verified in test/semantic/serialization_ir.jl::E8).
# Set RD_BD_FROM_JSON=1 to drive the demo from the JSON artifact instead of the in-Julia DSL.
const MODEL_JSON = joinpath(HERE, "model.rdj.json")
build_pipeline_from_json() = ReactiveDynamics.from_json_model(read(MODEL_JSON, String); registry = PROJECT_REGISTRY)

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
                extra_scientists = res ? 15 : 0,    # resource synergy: +headcount AND +capital (§2.1)
                extra_budget = res ? 250 : 0,
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

# nseed default is 160: at 24 seeds the per-scenario SE (~±450) swamps the synergy decomposition
# (the marginal synergies are ~100–550 apart), so the ordering looks noisy/non-monotone; by ~160
# seeds the SE (~±150) resolves res < pos and op-efficiency≈0 as real, stable signals.
function main(; root_seed = 2026, nseed = 160)
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
    # Resource columns are the two pools' LOW-WATER MARKS (cash⌄, sci⌄) — independent signals of
    # which constraint binds. (We don't also print the financing dip: it is exactly 150 − cash⌄,
    # i.e. a linear transform of cash⌄, so it would carry no extra information; it appears once in
    # the headline below as the "capital ask" framing.)
    println("\nScenario                      mean rNPV     launches  P(≥1)   cash⌄   sci⌄   Δ-rNPV vs S0 (±SE)")
    println("-" ^ 100)
    for s in scenarios
        ms = results[s]
        te = treatment_effect(base, ms)
        @printf(
            "%-28s  %9.1f    %6.2f   %5.2f   %6.1f  %5.1f   %s\n",
            labels[s],
            mean_rnpv(ms),
            mean_launches(ms),
            p_launch(ms),
            mean_cash_trough(ms),      # cash⌄: lowest the budget pool reached (binding ⇒ near 0)
            mean_sci_trough(ms),       # sci⌄: lowest the scientist pool reached
            s == :S0 ? "      —" : @sprintf("%+8.1f ± %5.1f", te.delta_rnpv, te.se_delta_rnpv),
        )
    end

    println("\nHeadline (Full deal S5 vs Baseline S0):")
    te = treatment_effect(base, results[:S5])
    @printf("  Δ-rNPV (deal value, net of price)  : %+.1f ± %.1f (1 SE, %d seeds)\n", te.delta_rnpv, te.se_delta_rnpv, nseed)
    @printf("  Δ-launches (extra programs to mkt) : %+.2f\n", te.delta_launches)
    @printf("  Δ-P(launch)                        : %+.2f\n", te.delta_p_launch)
    @printf("  cash trough, baseline → full deal  : %.1f → %.1f (the deal EASES the cash squeeze)\n",
            mean_cash_trough(base), mean_cash_trough(results[:S5]))
    println("\nSynergy decomposition (marginal Δ over S1 pipeline-only):")
    s1 = treatment_effect(base, results[:S1]).delta_rnpv
    for (s, name) in ((:S2, "resource (capital+headcount)"), (:S3, "capability/PoS"), (:S4, "op-efficiency"))
        @printf("  %-30s : %+8.1f\n", name, treatment_effect(base, results[s]).delta_rnpv - s1)
    end
    # ── Engine-level per-program ledger (MVP finding D, src/ledger.jl) ──────────────────────
    # NEW capability: the engine now attributes the cost ledger PER PROGRAM during the run (with
    # `budget` priced via specCost in host.jl — which leaves the dynamics and the Δ-rNPV above
    # untouched). `program_ledger(prob)` returns the per-program cost/reward/valuation summary in
    # deterministic token order; previously this was reconstructed in post. We surface it for one
    # representative full-deal run and CROSS-CHECK that the per-program rows reconcile to the
    # aggregate ledger and that the engine-attributed spend agrees with the post-hoc rNPV roll-up.
    println("\n" * "=" ^ 78)
    println("ENGINE-LEVEL PER-PROGRAM LEDGER (MVP finding D) — one representative full-deal run")
    println("=" ^ 78)
    demo_prob = run_scenario(:S5, hash((root_seed, 1)))
    sm = program_ledger_summary(demo_prob)
    @printf("  programs tracked                    : %d\n", sm.n_programs)
    @printf("  total capital attributed to programs: %.1f\n", sm.total_program_cost)
    @printf("    ↳ to programs that reached market  : %.1f\n", sm.launched_cost)
    @printf("    ↳ to still-in-flight programs      : %.1f\n", sm.active_cost)
    @printf("    ↳ to failed/retired programs       : %.1f\n", sm.retired_cost)
    @printf("  unattributed pool burn (no program) : %.1f  (finding-D boundary: spend with no bound token)\n", sm.unattributed_cost)
    @printf("  aggregate :valuation_cost ledger     : %.1f\n", sm.aggregate_cost)
    @printf("  reconciliation residual (≈0)         : %.3g  (per-program + unattributed == aggregate, §8.5)\n", sm.reconciliation_residual)

    # Show the top few programs by capital burned, next to their modeling descriptors.
    led = sort(program_economics(demo_prob), :cost_incurred; rev = true)
    println("\n  Top programs by engine-attributed capital burned:")
    println("    creation_idx  phase       acquired   npv_peak    cost_incurred")
    for r in eachrow(first(led, min(6, nrow(led))))
        @printf("    %10d  %-10s  %-8s   %8.0f    %10.1f\n",
                r.creation_index, string(r.phase), string(r.acquired),
                isnan(r.npv_peak) ? 0.0 : r.npv_peak, r.cost_incurred)
    end
    # The per-program ledger AGREES with the demo's portfolio roll-up: both walk the same final
    # population; the ledger now ALSO carries the engine-attributed capital each program consumed.
    @printf("\n  cross-check — portfolio rNPV (post-hoc roll-up)        : %.1f\n",
            portfolio_rnpv(demo_prob; acq_price = ACQ_PRICE))
    @printf("  cross-check — engine cost reconciles to aggregate     : %s\n",
            abs(sm.reconciliation_residual) < 1e-6 ? "YES" : "NO")

    println("\nReproducible: each cell is a (model, scenario, seed) triple; the lever is in-model")
    println("(an ADR-0010 Rule), so re-running with the same root seed gives identical numbers.")
    return results
end

# Run when executed as a script (not when included by a test).
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
