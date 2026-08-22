# Phase-1 Stage D semantic tests — declarative initial state + checkpoint (ADR 0007, §10).
#
# Covers: the declarative `population[]` initial marking (both authoring forms — PopulationEntry
# count+attribute-exprs, and explicit host token list); reflected-count consistency at t=0;
# marking determinism under (model, seed); the completed reinit! that rebuilds the token
# population so structured runs replay (§4 D7); the Live-phase guard refusing reindexers
# (equalize!); and dump_state/restore round-trip + resume. These make the BD demo's starting
# portfolio reproducible serializable input (MVP finding H).

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

# A ProjectToken kind with a phase attribute + npv, in RD scope (the @register/@aagent idiom).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct InitProjectToken
        phase::Symbol
        npv::Float64
    end
    function InitProjectToken(phase, npv)
        return InitProjectToken(
            "IP" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            npv,
        )
    end
end

# registry constructor the declarative population[] resolves `kind` against (ADR 0006 §C).
const INIT_REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RDX.InitProjectToken(get(f, :phase, :Phase2), get(f, :npv, 100.0)),
)

# A phase-advance model used across the tests.
function init_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => adv, cycletime => 1.0, probability => 1.0
    end
    RDX.register_token_kind!(net, :Project)
    return net
end

phases(p) = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(p, "structured")))]))
ntok(p) = length(collect(values(RDX.inners(RDX.getagent(p, "structured")))))

@testset "Declarative initial state + checkpoint (ADR 0007, §10)" begin

    # ── (B) declarative population[]: count + attribute exprs form ──────────────────────
    @testset "population[] PopulationEntry form instantiates count tokens with seeded attributes" begin
        pop = [
            RDX.PopulationEntry(
                :Project, :Project; count = 5,
                attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 200.0)
            ),
        ]
        p = ReactionNetworkProblem(
            init_model(); tspan = 4, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = pop
        )
        @test ntok(p) == 5                                   # 5 instances created at construction
        @test all(==(:Phase2), [t.phase for t in values(RDX.inners(RDX.getagent(p, "structured")))])
        @test sort(collect(values(p.creation_index))) == [1, 2, 3, 4, 5]   # per-species creation indices
    end

    # ── (B) explicit host-token list form ──────────────────────────────────────────────
    @testset "population[] explicit host-token list form entangles the given agents" begin
        pop = [RDX.InitProjectToken(:Phase2, 111.0), RDX.InitProjectToken(:Phase1, 222.0)]
        p = ReactionNetworkProblem(
            init_model(); tspan = 4, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = pop
        )
        @test ntok(p) == 2
        @test phases(p) == ["Phase1", "Phase2"]
    end

    # ── invariant 2: reflected-count consistency at t=0 (u == active population) ─────────
    @testset "structured u column at t=0 equals the active initial population" begin
        pop = [
            RDX.PopulationEntry(
                :Project, :Project; count = 3,
                attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 100.0)
            ),
        ]
        p = ReactionNetworkProblem(
            init_model(); tspan = 4, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = pop
        )
        @test p.u[RDX.find_index(:Project, p)] == 3.0
    end

    # ── invariant 1: marking determinism under (model, seed) ────────────────────────────
    @testset "the initial marking (count, sampled attrs, creation indices) is reproducible" begin
        # a sampled attribute draws from the seeded stream — two same-seed constructions agree
        mk(seed) = begin
            pop = [
                RDX.PopulationEntry(
                    :Project, :Project; count = 4,
                    attributes = Dict(
                        :phase => QuoteNode(:Phase2),
                        :npv => :(rand(state.rng, Normal(100.0, 10.0)))
                    )
                ),
            ]
            p = ReactionNetworkProblem(
                init_model(); tspan = 2, dt = 1.0, seed = seed,
                registry = INIT_REGISTRY, population = pop
            )
            sort([t.npv for t in values(RDX.inners(RDX.getagent(p, "structured")))])
        end
        @test mk(123) == mk(123)        # same seed ⇒ identical sampled NPVs
        @test mk(123) != mk(456)        # different seed ⇒ (almost surely) different
    end

    # ── invariant 6: reinit! rebuilds the marking ⇒ structured runs replay (§4 D7) ──────
    @testset "reinit! rebuilds the initial token population so a structured run reproduces" begin
        pop = [
            RDX.PopulationEntry(
                :Project, :Project; count = 4,
                attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 100.0)
            ),
        ]
        p = ReactionNetworkProblem(
            init_model(); tspan = 5, dt = 1.0, seed = 7,
            registry = INIT_REGISTRY, population = pop
        )
        simulate(p); sol1 = copy(p.sol); ph1 = phases(p)
        AlgebraicAgents._reinit!(p)
        @test ntok(p) == 4              # population rebuilt to the t=0 marking (not end-state)
        @test phases(p) == ["Phase2", "Phase2", "Phase2", "Phase2"]   # all back at Phase2
        simulate(p)
        @test p.sol == sol1             # trajectory reproduced
        @test phases(p) == ph1          # final phase distribution reproduced
    end

    # ── invariant 6 (explicit-token form): reinit! restores the SAME objects' t=0 attributes ──
    @testset "reinit! restores explicit host-token attributes (advanced/retired tokens reset to t=0)" begin
        # explicit-host-token population: the SAME objects are re-entangled on reinit, so their
        # mutated fields (phase advanced to :Phase3, or species soft-retired to :removed) must be
        # restored to the captured t=0 snapshot — else the second run does not reproduce.
        toks = [RDX.InitProjectToken(:Phase2, 100.0), RDX.InitProjectToken(:Phase2, 200.0)]
        p = ReactionNetworkProblem(
            init_model(); tspan = 5, dt = 1.0, seed = 3,
            registry = INIT_REGISTRY, population = toks
        )
        simulate(p); ph1 = phases(p)
        @test any(!=(:Phase2), [t.phase for t in toks])   # at least one advanced/changed during the run
        AlgebraicAgents._reinit!(p)
        # the same token objects are back at their t=0 phase/species
        @test all(==(:Phase2), [t.phase for t in toks])
        @test all(==(:Project), [RDX.get_species(t) for t in toks])
        simulate(p)
        @test phases(p) == ph1                             # second run reproduces the first
    end

    # ── (A) Live-phase guard: equalize! refuses on a constructed model ──────────────────
    @testset "equalize! (a rem_parts! reindexer) refuses on a live constructed model (ADR 0004 INV-2)" begin
        p = ReactionNetworkProblem(
            init_model(); tspan = 3, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = []
        )
        @test p.live == true
        @test_throws Exception RDX.equalize!(p, [])
        # the authoring-phase equalize! on a schema stays legal (not tested here; it's the §7 path)
    end

    # ── (C) dump_state / restore round-trip + resume (the checkpoint, §10.5) ────────────
    @testset "dump_state at a tick boundary restores to an identical state and resumes equally" begin
        spec = @reaction_network begin
            @deterministic(1.0),
                @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
                name => adv, cycletime => 0.0, probability => 1.0
        end
        RDX.register_token_kind!(spec, :Project)
        pop = [
            RDX.PopulationEntry(
                :Project, :Project; count = 4,
                attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 150.0)
            ),
        ]
        p = ReactionNetworkProblem(
            spec; tspan = 10, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = pop
        )
        simulate(p, 2)                          # ct=0 ⇒ no in-flight at the boundary
        @test isempty(p.ongoing_transitions)
        d = RDX.dump_state(p)
        @test d.t == 2.0 && length(d.tokens) == 4
        p2 = RDX.restore(spec, d; registry = INIT_REGISTRY)
        @test p2.t == 2.0 && ntok(p2) == 4 && p2.u == p.u
        simulate(p); simulate(p2)
        @test phases(p) == phases(p2)           # resume reproduces the continuation
    end

    # ── dump_state refuses mid-cycle (Milestone-1 scope: clean tick boundary only) ──────
    @testset "dump_state refuses when an in-flight transition is mid-cycle (documented M1 scope)" begin
        spec = @reaction_network begin
            @deterministic(1.0),
                @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
                name => adv, cycletime => 5.0, probability => 1.0    # long cycle ⇒ in-flight mid-run
        end
        RDX.register_token_kind!(spec, :Project)
        pop = [
            RDX.PopulationEntry(
                :Project, :Project; count = 2,
                attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 100.0)
            ),
        ]
        p = ReactionNetworkProblem(
            spec; tspan = 10, dt = 1.0, seed = 1,
            registry = INIT_REGISTRY, population = pop
        )
        simulate(p, 1)                          # an instance is now mid-cycle
        @test !isempty(p.ongoing_transitions)
        @test_throws Exception RDX.dump_state(p)
    end
end
