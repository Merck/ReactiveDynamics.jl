# Phase-0.6 semantic tests — Analysis & observability (ADR 0013 / CONTRACT §14).
#
# Covers the three additive read-only mechanisms:
#   §14.1 the per-token trajectory log — `log_token_fields` opt-in, `push_token_trajectory_row!` at
#         the per-program-ledger seam, the `token_trajectory` read API (+ predicate scope), the
#         "typical" helpers `representative_token`/`trajectory_envelope`, determinism (§4 D4) and the
#         `_reinit!` reset (§4 D7).
#   §14.2 the ensemble runner — `ensemble`/`summarize`/`treatment_effect`, per-member hash seeding
#         (§4 D8), reproducibility (§4 D9), and `EnsembleProblem` as an AA-readable node (Invariant 4).
#   §14.3 the export bundle — `export_run`/`export_ensemble` JSON+CSV core, the JSON round-trip
#         (§8 / Invariant 5), and the Arrow weakdep sibling (when RDArrowExt is loaded).
#
# These run against the built engine — real @tests, not pins. Default runs are entropy-seeded, so
# every problem is built with an explicit seed=.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames, JSON
using Arrow   # trigger RDArrowExt so the §14.3 Arrow-sibling assertions run (not skip)

RD = ReactiveDynamics

# A structured token kind with a `phase` AND a numeric `value` field, so the trajectory log records
# both a Symbol path (phase) and a Real path (value, for the envelope/medoid numeric metric).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct TrajProjectToken
        phase::Symbol
        value::Float64
    end
    function TrajProjectToken(phase, value)
        return TrajProjectToken(
            "TP" * string(rand(1:10^9)),
            :Project,
            nothing,
            Tuple{Symbol,Float64,ReactiveDynamics.Transition}[],
            phase,
            value,
        )
    end
end

# The trajectory-log opt-in for this kind (ADR 0013 §A1): record phase + value each tick. The kind
# is defined into the ReactiveDynamics module by @register, so it is referenced as `RD.TrajProjectToken`.
RD.log_token_fields(t::RD.TrajProjectToken) = (; phase = t.phase, value = t.value)

# A coarse advance model identical in spirit to the program-ledger test: @select a Phase1 Project,
# burn budget at rate, @advance to Phase2. Used to exercise the trajectory log over real ticks.
function traj_model(; budget0 = 100, cost = 1.0)
    acs = @ReactionNetworkSchema begin
        @deterministic(1.0),
        @select(Project, phase == :Phase1) + 2 * @rate(budget) --> @advance(phase, :Phase2),
        name => adv, cycletime => 1.0, probability => 1.0
    end
    RD.register_structured_species!(acs, :Project)
    bi = findfirst(==(:budget), acs[:, :specName])
    acs[bi, :specInitVal] = Float64(budget0)
    acs[bi, :specCost] = cost
    @prob_meta acs tspan = 5 dt = 1.0
    return acs
end

build_traj_prob(seed; pop = [RD.TrajProjectToken(:Phase1, 1.0), RD.TrajProjectToken(:Phase1, 2.0)]) =
    ReactionNetworkProblem(traj_model(); seed = seed, population = pop)

# Strip the random program names so two seeded runs are compared on field-value CONTENT (the demo's
# token names are randstring-based and intentionally non-deterministic; §4 D4 is about field values).
content(df) = select(df, Not(:program))

@testset "Analysis & observability (ADR 0013 / CONTRACT §14)" begin

    # ── §14.1 per-token trajectory log ───────────────────────────────────────────────────
    @testset "§14.1 trajectory log records opted-in fields each tick, in token order" begin
        p = build_traj_prob(1)
        simulate(p)
        df = token_trajectory(p)
        @test nrow(df) > 0
        # the long form carries t, program, species, then the opted-in fields
        @test issubset(["t", "program", "species", "phase", "value"], names(df))
        # one token's life is a subset filtered by name
        nm = df[1, :program]
        one = token_trajectory(p, nm)
        @test all(==(nm), one.program)
        @test nrow(one) >= 1
    end

    @testset "§14.1 opt-in is bounded: a non-logging kind contributes no rows (Invariant 2)" begin
        # A bare problem with NO opted-in tokens logs nothing.
        acs = traj_model()
        p = ReactionNetworkProblem(acs; seed = 2, population = [])
        simulate(p)
        @test nrow(token_trajectory(p)) == 0
    end

    @testset "§14.1 trajectory determinism under (model, seed) (§4 D4)" begin
        a = build_traj_prob(7); simulate(a)
        b = build_traj_prob(7); simulate(b)
        @test content(token_trajectory(a)) == content(token_trajectory(b))
    end

    @testset "§14.1 reinit! clears the trajectory log and a re-run reproduces it (§4 D7)" begin
        p = build_traj_prob(11); simulate(p)
        first_run = content(token_trajectory(p))
        AlgebraicAgents.reinit!(p)
        @test nrow(token_trajectory(p)) == 0        # cleared on reinit
        simulate(p)
        @test content(token_trajectory(p)) == first_run
    end

    @testset "§14.1 typical helpers: representative_token + trajectory_envelope" begin
        p = build_traj_prob(13); simulate(p)
        rep = representative_token(p)
        @test rep isa AbstractString          # a token name (medoid)
        env = trajectory_envelope(p)
        @test issubset(["align", "field", "median", "q25", "q75", "n"], names(env))
        @test :value in env.field             # the numeric field has an envelope band
        @test all(env.q25 .<= env.median .<= env.q75)
    end

    @testset "§14.1 predicate scope reuses @select machinery (§9.5)" begin
        p = build_traj_prob(17); simulate(p)
        # Phase2 tokens are those that advanced; the scoped trajectory has only their rows.
        adv = token_trajectory(p, RD.TokenPredicate(:Project, [RD.Clause(:phase, :(==), QuoteNode(:Phase2))]))
        names_adv = Set(adv.program)
        all_df = token_trajectory(p)
        @test issubset(names_adv, Set(all_df.program))
    end

    # ── §14.2 ensemble runner ────────────────────────────────────────────────────────────
    @testset "§14.2 ensemble seeds members hash((root_seed,k)) and is reproducible (§4 D8/D9)" begin
        e1 = ensemble(s -> (p = build_traj_prob(s); simulate(p); p); nseed = 5, root_seed = 2026)
        @test length(e1.members) == 5
        @test e1.seeds == UInt64[hash((2026, k)) for k in 1:5]
        @test e1.mode === :rebuild
        # reproducible: same (root_seed, nseed, build) → same per-run metric
        metric(p) = nrow(token_trajectory(p))
        e2 = ensemble(s -> (p = build_traj_prob(s); simulate(p); p); nseed = 5, root_seed = 2026)
        @test summarize(e1, metric).mean == summarize(e2, metric).mean
    end

    @testset "§14.2 summarize reports mean/sem/quantiles" begin
        e = ensemble(s -> (p = build_traj_prob(s); simulate(p); p); nseed = 6, root_seed = 5)
        s = summarize(e, p -> last(p.sol.budget))
        @test s.n == 6
        @test s.q25 <= s.median <= s.q75
        @test s.sem >= 0.0
    end

    @testset "§14.2 treatment_effect is the unpaired Δ with se=√(var_b/n_b+var_d/n_d)" begin
        # baseline: cost 1; deal: cost 0 (no burn). Final budget differs → a real treatment effect.
        base = ensemble(s -> (p = ReactionNetworkProblem(traj_model(; cost = 1.0); seed = s,
            population = [RD.TrajProjectToken(:Phase1, 1.0)]); simulate(p); p); nseed = 6, root_seed = 9)
        deal = ensemble(s -> (p = ReactionNetworkProblem(traj_model(; cost = 0.0); seed = s,
            population = [RD.TrajProjectToken(:Phase1, 1.0)]); simulate(p); p); nseed = 6, root_seed = 9)
        te = treatment_effect(base, deal, p -> last(p.sol.budget))
        @test te.n_baseline == 6 && te.n_deal == 6
        @test te.deal >= te.baseline           # zero-cost deal leaves more budget
        @test te.se >= 0.0
        @test isapprox(te.delta, te.deal - te.baseline; atol = 1e-9)
    end

    @testset "§14.2 EnsembleProblem is a readable AA node (Invariant 4)" begin
        e = ensemble(s -> (p = build_traj_prob(s); simulate(p); p); nseed = 4, root_seed = 3)
        obs = AlgebraicAgents.observables(e)
        @test :budget in obs
        @test AlgebraicAgents.getobservable(e, :budget) isa Real
        # unknown name is a hard error, not AA's silent @error fall-through
        @test_throws Exception AlgebraicAgents.getobservable(e, :no_such_obs)
        # the hierarchy must hold ALL members as distinct children (members share the default name,
        # so they are renamed member_<k> before entangle! — else inners collapses to one).
        @test length(AlgebraicAgents.inners(e)) == 4
    end

    # ── §14.3 export bundle ──────────────────────────────────────────────────────────────
    @testset "§14.3 export_run writes the JSON+CSV core and round-trips the JSON streams (§8)" begin
        p = build_traj_prob(21); simulate(p)
        dir = mktempdir()
        export_run(p, dir)
        @test isfile(joinpath(dir, "trajectory.csv"))
        @test isfile(joinpath(dir, "ledger.csv"))
        @test isfile(joinpath(dir, "events.json"))
        @test isfile(joinpath(dir, "tokens.json"))
        @test isfile(joinpath(dir, "run.json"))
        # manifest pins (model_hash, seed) — traceable to a replayable run (Invariant 5, §4 D6)
        man = JSON.parsefile(joinpath(dir, "run.json"))
        @test man["model_hash"] == string(hash(p.acs))
        @test man["seed"] == string(p.seed)
        # the JSON streams round-trip structurally
        ev = JSON.parsefile(joinpath(dir, "events.json"))
        @test ev isa AbstractVector && !isempty(ev)
        tk = JSON.parsefile(joinpath(dir, "tokens.json"))
        @test tk isa AbstractDict
    end

    @testset "§14.3 export_ensemble writes per-member bundles + ensemble.json with summary" begin
        e = ensemble(s -> (p = build_traj_prob(s); simulate(p); p); nseed = 3, root_seed = 4)
        dir = mktempdir()
        export_ensemble(e, dir; metric = p -> last(p.sol.budget))
        @test isdir(joinpath(dir, "member_1"))
        @test isfile(joinpath(dir, "member_3", "run.json"))
        ej = JSON.parsefile(joinpath(dir, "ensemble.json"))
        @test ej["nseed"] == 3
        @test ej["mode"] == "rebuild"
        @test haskey(ej, "summary")
        @test length(ej["seeds"]) == 3
    end

    # ── §14.3 Arrow weakdep sibling (only when RDArrowExt is loaded) ──────────────────────
    @testset "§14.3 Arrow sibling is written when Arrow is available (RDArrowExt)" begin
        if RD._arrow_available()
            p = build_traj_prob(23); simulate(p)
            dir = mktempdir()
            export_run(p, dir)
            @test isfile(joinpath(dir, "trajectory.arrow"))
            @test isfile(joinpath(dir, "ledger.arrow"))
        else
            @test_skip "Arrow not loaded — .arrow siblings skipped (core writes CSV/JSON only)"
        end
    end
end
