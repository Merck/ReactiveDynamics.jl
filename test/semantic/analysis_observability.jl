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
            "TP" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Firing}[],
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
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) + 2 * @rate(budget) --> @advance(phase, :Phase2),
            name => adv, cycletime => 1.0, probability => 1.0
    end
    RD.register_token_kind!(net, :Project)
    bi = findfirst(==(:budget), net[:, :placeName])
    net[bi, :placeInitVal] = Float64(budget0)
    net[bi, :placeCost] = cost
    @prob_meta net tspan = 5 dt = 1.0
    return net
end

build_traj_prob(seed; pop = [RD.TrajProjectToken(:Phase1, 1.0), RD.TrajProjectToken(:Phase1, 2.0)]) =
    ReactionNetworkProblem(traj_model(); seed = seed, population = pop)

# A registry so the PopulationEntry authoring form can resolve `:Project` to a host constructor
# (ADR 0006 §C). Used by the mode-(a)≡mode-(b) equivalence test to exercise the SEED-DEPENDENT
# initial-attribute path: `value` is a seeded `rand` draw, so a member's t=0 attributes depend on
# its seed — the reseed ordering (new stream installed BEFORE instantiate_population!) is what makes
# a reinit-reseeded member's sampled attributes match a fresh build(seed) (§14.2 equivalence).
const TRAJ_REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RD.TrajProjectToken(get(f, :phase, :Phase1), get(f, :value, 1.0)),
)
build_traj_prob_pe(seed) = ReactionNetworkProblem(
    traj_model(); seed = seed, registry = TRAJ_REGISTRY,
    population = [
        RD.PopulationEntry(
            :Project, :Project; count = 3,
            attributes = Dict(:phase => QuoteNode(:Phase1), :value => :(rand(state.rng, Normal(5.0, 2.0))))
        ),
    ],
)

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
        # the long form carries t, program, place, then the opted-in fields
        @test issubset(["t", "program", "place", "phase", "value"], names(df))
        # one token's life is a subset filtered by name
        nm = df[1, :program]
        one = token_trajectory(p, nm)
        @test all(==(nm), one.program)
        @test nrow(one) >= 1
    end

    @testset "§14.1 opt-in is bounded: a non-logging kind contributes no rows (Invariant 2)" begin
        # A bare problem with NO opted-in tokens logs nothing.
        net = traj_model()
        p = ReactionNetworkProblem(net; seed = 2, population = [])
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
        base = ensemble(
            s -> (
                p = ReactionNetworkProblem(
                    traj_model(; cost = 1.0); seed = s,
                    population = [RD.TrajProjectToken(:Phase1, 1.0)]
                ); simulate(p); p
            ); nseed = 6, root_seed = 9
        )
        deal = ensemble(
            s -> (
                p = ReactionNetworkProblem(
                    traj_model(; cost = 0.0); seed = s,
                    population = [RD.TrajProjectToken(:Phase1, 1.0)]
                ); simulate(p); p
            ); nseed = 6, root_seed = 9
        )
        te = treatment_effect(base, deal, p -> last(p.sol.budget))
        @test te.n_baseline == 6 && te.n_deal == 6
        @test te.deal >= te.baseline           # zero-cost deal leaves more budget
        @test te.se >= 0.0
        @test isapprox(te.delta, te.deal - te.baseline; atol = 1.0e-9)
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

    # ── §14.2 ensemble mode (b): reinit-reseed member reuse (ADR 0013 §14.2, WS-B1) ────────
    @testset "§14.2 reinit-reseed reproduces a fresh build(seed) — single-member A/B (§4 D7/D8)" begin
        # Reinit-reseeding member 1 to member 2's seed must reproduce EXACTLY what a fresh
        # build(seed₂) produces: same final marking AND same seed-sampled initial attributes (the
        # reseed installs the new stream BEFORE the t=0 marking is re-sampled). Uses the SEED-DEPENDENT
        # PopulationEntry form so the initial `value` draws are part of what must match.
        s1, s2 = hash((99, 1)), hash((99, 2))
        reused = build_traj_prob_pe(s1); simulate(reused)
        fresh = build_traj_prob_pe(s2); simulate(fresh)
        AlgebraicAgents.reinit!(reused; seed = s2)
        @test nrow(token_trajectory(reused)) == 0            # reinit clears the log
        simulate(reused)
        @test reused.seed == s2                              # the realized seed is updated (D6/D8)
        @test content(token_trajectory(reused)) == content(token_trajectory(fresh))
        @test last(reused.sol.budget) == last(fresh.sol.budget)
        # the seed-sampled initial `value`s match too (the reseed-ordering guarantee)
        v_reused = sort(collect(skipmissing(token_trajectory(reused).value)))
        v_fresh = sort(collect(skipmissing(token_trajectory(fresh).value)))
        @test v_reused == v_fresh
    end

    @testset "§14.2 the no-seed reinit! path is unchanged (byte-identical re-run, §4 D7)" begin
        # WS-B1 must not perturb the plain reinit! contract (this backs the existing green reinit
        # test): reinit!(p) with NO seed still restores the construction stream and reproduces.
        p = build_traj_prob_pe(hash((99, 5))); simulate(p)
        first_run = content(token_trajectory(p)); first_sol = copy(p.sol)
        seed0 = p.seed
        AlgebraicAgents.reinit!(p)
        @test p.seed == seed0                                # seed unchanged on the no-seed path
        simulate(p)
        @test content(token_trajectory(p)) == first_run
        @test p.sol == first_sol
    end

    @testset "§14.2 mode-(a) ≡ mode-(b) equivalence — member-for-member (THE headline gate)" begin
        metric(m) = last(m.sol.budget)
        # Verify the equivalence for BOTH population forms (handoff Risk note): the SEED-DEPENDENT
        # PopulationEntry form (initial attrs re-sampled through the reseeded stream) AND the
        # explicit host-token form (seed-INDEPENDENT initial marking via restore_token_snapshot!,
        # only the simulation draws differ).
        for build in (build_traj_prob_pe, build_traj_prob)
            ea = ensemble(s -> (p = build(s); simulate(p); p); nseed = 6, root_seed = 7, mode = :rebuild)
            eb = ensemble(s -> (p = build(s); simulate(p); p); nseed = 6, root_seed = 7, mode = :reinit)
            @test ea.mode === :rebuild && eb.mode === :reinit
            # drawable-node invariant holds for BOTH modes (snapshots entangle as child nodes)
            @test length(AlgebraicAgents.inners(ea)) == 6
            @test length(AlgebraicAgents.inners(eb)) == 6
            @test length(eb.members) == 6
            @test eb.seeds == ea.seeds
            # member-for-member: same per-member metric, same summarize, same trajectory content
            @test [metric(m) for m in ea.members] == [metric(m) for m in eb.members]
            @test summarize(ea, metric) == summarize(eb, metric)
            @test [content(token_trajectory(m)) for m in ea.members] ==
                [content(token_trajectory(m)) for m in eb.members]
        end
    end

    @testset "§14.2 treatment_effect Δ is identical under mode (a) vs mode (b)" begin
        # The BD A/B lever comparison: baseline (cost 1) vs deal (cost 0). The unpaired Δ must be
        # the SAME whether the two arms are built by rebuild or by reinit-reseed.
        base_pop() = [RD.TrajProjectToken(:Phase1, 1.0)]
        base_build(cost) = s -> (
            p = ReactionNetworkProblem(
                traj_model(; cost = cost);
                seed = s, population = base_pop()
            ); simulate(p); p
        )
        rnpv(m) = last(m.sol.budget)
        te = Dict(
            mode => treatment_effect(
                    ensemble(base_build(1.0); nseed = 6, root_seed = 9, mode = mode),
                    ensemble(base_build(0.0); nseed = 6, root_seed = 9, mode = mode),
                    rnpv,
                ) for mode in (:rebuild, :reinit)
        )
        @test te[:rebuild].delta == te[:reinit].delta
        @test te[:rebuild].se == te[:reinit].se
        @test te[:rebuild].baseline == te[:reinit].baseline
        @test te[:rebuild].deal == te[:reinit].deal
    end

    @testset "§14.2 ensemble rejects an unknown mode" begin
        @test_throws Exception ensemble(
            s -> (p = build_traj_prob(s); simulate(p); p);
            nseed = 2, root_seed = 1, mode = :bogus
        )
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
        @test man["model_hash"] == string(hash(p.network))
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
