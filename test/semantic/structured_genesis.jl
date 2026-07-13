# Phase-1 semantic tests — structured-token genesis as a transition PRODUCT (`@structured` RHS).
#
# The `@structured(ctor(...))` RHS op is the "birth" leg of the structured-token lifecycle, the
# counterpart to token_filtration.jl's @select (bind) and @advance (mutate) legs: a transition
# whose PRODUCT is a freshly-constructed host token, entangled live into the structured pool by
# `structured_rhs` (src/solvers.jl:551-562 — context_eval the ctor expr, then entangle!). It is
# genesis-as-first-class-transition-product — structurally parallel to `∅ --> plain_species` — as
# distinct from the imperative `AddToken`-in-a-Rule decision-channel path (ADR 0010, rules_decisions.jl).
#
# Parsing: recognized as a RHS macrocall alongside @move/@advance (reaction_parser.jl:102,
# create.jl:302). Serialization boundary: a `@structured` body is host Julia (an Expr), so it does
# NOT round-trip through the eval-free JSON IR — to_json_model raises the documented error
# (serialize.jl:901-905). All assertions run against the built engine (real @tests, not pins).

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

# A ProjectToken structured kind (RD-scoped via the @register/@aagent idiom — the engine's genesis
# path context_evals the ctor inside RD's namespace, so a bare `GenProjectToken` resolves there).
# `born` records the clock at construction, to prove the ctor sees live state; `npv` may be drawn.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct GenProjectToken
        phase::Symbol
        npv::Float64
        born::Float64
    end
    function GenProjectToken(phase, npv, born)
        return GenProjectToken(
            "GP" * string(rand(1:10^9)),
            :Project,
            nothing,
            Tuple{Symbol,Float64,ReactiveDynamics.Transition}[],
            phase,
            npv,
            born,
        )
    end
end

livetokens(p) = collect(values(RDX.inners(RDX.getagent(p, "structured"))))
nphase(p, ph) = count(t -> t.phase == ph, livetokens(p))

@testset "Structured-token genesis as a transition product (@structured RHS)" begin

    # ── an empty-LHS source mints one token per tick; identity is fresh each time ─────────
    @testset "∅ --> @structured(ctor(...)) mints a fresh token per firing, tracked in state.u" begin
        acs = @ReactionNetworkSchema begin
            @deterministic(1.0),
            ∅ --> @structured(GenProjectToken(:Phase1, 100.0, @t())),
            name => genesis
        end
        RDX.register_structured_species!(acs, :Project)
        @prob_meta acs tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(acs; seed = 1)
        @test isempty(livetokens(p))                       # nothing at t=0
        simulate(p)
        toks = livetokens(p)
        # one @deterministic(1.0) genesis per tick over the run (6 tick boundaries t=0..5)
        @test length(toks) == 6
        @test all(t -> RDX.get_species(t) == :Project, toks)
        @test all(t -> t.phase == :Phase1 && t.npv == 100.0, toks)
        # each minted token has its OWN identity — the entangle! pool is keyed by unique name
        @test length(unique(AlgebraicAgents.getname.(toks))) == length(toks)
        # the plain state.u count tracks the live structured population (update_u_structured!)
        ci = RDX.find_index(:Project, p)
        @test p.u[ci] == 6.0
    end

    # ── the constructor sees LIVE state (@t()) and may DRAW from the seeded RNG ────────────
    @testset "genesis ctor reads @t() and draws from state.rng; reproducible under seed" begin
        function genesis_dynamic(seed)
            acs = @ReactionNetworkSchema begin
                @deterministic(1.0),
                ∅ --> @structured(GenProjectToken(:Phase1, rand(state.rng, Normal(100.0, 10.0)), @t())),
                name => genesis
            end
            RDX.register_structured_species!(acs, :Project)
            @prob_meta acs tspan = 4 dt = 1.0
            p = ReactionNetworkProblem(acs; seed = seed)
            simulate(p)
            p
        end
        p = genesis_dynamic(1)
        toks = livetokens(p)
        # @t() is captured at construction: one token born at each of t = 0,1,2,3,4
        @test sort([t.born for t in toks]) == [0.0, 1.0, 2.0, 3.0, 4.0]
        # the sampled npvs are genuinely varied (not the constant literal case above)
        @test length(unique(round.([t.npv for t in toks]; digits = 6))) > 1
        # same seed ⇒ identical draws (genesis draws from the state's seeded rng)
        npvs(pp) = sort([t.npv for t in livetokens(pp)])
        @test npvs(p) ≈ npvs(genesis_dynamic(1))
        @test npvs(p) != npvs(genesis_dynamic(2))          # a different seed diverges
    end

    # ── a minted token then FLOWS through a downstream @select/@advance pipeline ───────────
    @testset "genesis feeds a downstream @select/@advance leg (birth → select → advance)" begin
        acs = @ReactionNetworkSchema begin
            @deterministic(1.0),
            ∅ --> @structured(GenProjectToken(:Phase1, 100.0, @t())),
            name => genesis
            @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 1.0
        end
        RDX.register_structured_species!(acs, :Project)
        @prob_meta acs tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(acs; seed = 1)
        simulate(p)
        # every minted Phase1 token is bindable by the downstream leg and advances to Phase2;
        # the pipeline keeps only the most-recently-born token in Phase1 (steady one-per-tick flow).
        @test nphase(p, :Phase2) >= 1
        @test length(livetokens(p)) == nphase(p, :Phase1) + nphase(p, :Phase2)
        @test all(t -> t.phase in (:Phase1, :Phase2), livetokens(p))
    end

    # ── serialization boundary: a @structured RHS body is host Julia, NOT eval-free JSON ───
    @testset "to_json_model raises on a @structured RHS (host Expr body, not serializable)" begin
        acs = @ReactionNetworkSchema begin
            @deterministic(1.0),
            ∅ --> @structured(GenProjectToken(:Phase1, 100.0, @t())),
            name => genesis
        end
        RDX.register_structured_species!(acs, :Project)
        @prob_meta acs tspan = 2 dt = 1.0
        p = ReactionNetworkProblem(acs; seed = 1)
        # the eval-free JSON IR covers @advance / a typed AddToken rule; a host-built @structured
        # ctor cannot round-trip, and the exporter surfaces that explicitly rather than silently.
        @test_throws Exception RDX.to_json_model(p)
    end
end
