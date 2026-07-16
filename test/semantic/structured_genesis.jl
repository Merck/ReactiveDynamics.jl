# Phase-1 semantic tests — structured-token genesis as a transition PRODUCT (`@structured` RHS).
#
# The `@structured` RHS op is the "birth" leg of the structured-token lifecycle, the counterpart to
# token_filtration.jl's @select (bind) and @advance (mutate) legs: a transition whose PRODUCT is a
# freshly-constructed host token, entangled live into the structured pool by `structured_rhs`
# (src/solvers.jl). It is genesis-as-first-class-transition-product — structurally parallel to
# `∅ --> plain_species` — distinct from the imperative `AddToken`-in-a-Rule decision-channel path
# (ADR 0010, rules_decisions.jl).
#
# There is ONE authoring form (ADR 0005 §39 — the structured escape hatch promoted to a typed node):
#   `@structured(:Kind, field = value, …)` — registry-resolved (the RHS-product twin of AddToken,
#   sharing its `(state, fields::Dict) -> token` contract). It carries the kind NAME + field-value
#   ExprNodes, never the constructor, so it is eval-free-SERIALIZABLE and round-trips through
#   to_json_model / from_json_model / validate. The registry resolves the host constructor by name
#   at firing time (ADR 0006 §C).
#
# The RAW `@structured(Ctor(…))` constructor form was REMOVED: it was the only reactant construct
# that could not round-trip eval-free, so forbidding it makes eval-free serialization a TOTAL
# invariant (every genesis product is data). It is now rejected at CONSTRUCTION (create.jl), tested
# below. All assertions run against the built engine (real @tests, not pins).

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
import JSON

const RDX = ReactiveDynamics

# A ProjectToken structured kind (RD-scoped via the @register/@aagent idiom). The registry below
# resolves `:Project` to this constructor BY NAME at genesis time. `born` records the clock at
# construction, to prove the field exprs see live state; `npv` may be a seeded draw.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct GenProjectToken
        phase::Symbol
        npv::Float64
        born::Float64
    end
    function GenProjectToken(phase, npv, born)
        return GenProjectToken(
            "GP" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            npv,
            born,
        )
    end
end

# Registry the NAMED @structured form (and AddToken) resolve the kind through, BY NAME (ADR 0006
# §C): key => a `(state, fields::Dict) -> token` host constructor — the SAME convention AddToken
# uses (actions.jl:203-211), so the named genesis product and the rule action share one contract.
const GEN_REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RDX.GenProjectToken(
        get(f, :phase, :Phase1), get(f, :npv, 0.0), get(f, :born, -1.0)
    ),
)

livetokens(p) = collect(values(RDX.inners(RDX.getagent(p, "structured"))))
nphase(p, ph) = count(t -> t.phase == ph, livetokens(p))
phases_of(p) = sort(string.([t.phase for t in livetokens(p)]))

@testset "Structured-token genesis as a transition product (@structured RHS)" begin

    # ── an empty-LHS source mints one token per tick; identity is fresh each time ─────────
    @testset "∅ --> @structured(:Kind, …) mints a fresh token per firing, tracked in state.u" begin
        net = @reaction_network begin
            @deterministic(1.0),
                ∅ --> @structured(:Project, phase = :Phase1, npv = 100.0, born = @t()),
                name => genesis
        end
        RDX.register_structured_species!(net, :Project)
        @prob_meta net tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1, registry = GEN_REGISTRY)
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

    # ── field exprs see LIVE state (@t()) and may DRAW from the seeded RNG ─────────────────
    @testset "genesis field exprs read @t() and draw from state.rng; reproducible under seed" begin
        function genesis_dynamic(seed)
            net = @reaction_network begin
                @deterministic(1.0),
                    ∅ --> @structured(
                        :Project, phase = :Phase1,
                        npv = rand(state.rng, Normal(100.0, 10.0)), born = @t()
                    ),
                    name => genesis
            end
            RDX.register_structured_species!(net, :Project)
            @prob_meta net tspan = 4 dt = 1.0
            p = ReactionNetworkProblem(net; seed = seed, registry = GEN_REGISTRY)
            simulate(p)
            p
        end
        p = genesis_dynamic(1)
        toks = livetokens(p)
        # @t() is captured at genesis: one token born at each of t = 0,1,2,3,4
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
        net = @reaction_network begin
            @deterministic(1.0),
                ∅ --> @structured(:Project, phase = :Phase1, npv = 100.0, born = @t()),
                name => genesis
            @deterministic(1.0),
                @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
                name => adv12, cycletime => 1.0, probability => 1.0
        end
        RDX.register_structured_species!(net, :Project)
        @prob_meta net tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1, registry = GEN_REGISTRY)
        simulate(p)
        # every minted Phase1 token is bindable by the downstream leg and advances to Phase2;
        # the pipeline keeps only the most-recently-born token in Phase1 (steady one-per-tick flow).
        @test nphase(p, :Phase2) >= 1
        @test length(livetokens(p)) == nphase(p, :Phase1) + nphase(p, :Phase2)
        @test all(t -> t.phase in (:Phase1, :Phase2), livetokens(p))
    end

    # ── genesis round-trips through the eval-free JSON IR (the ADR 0005 §39 payoff) ────────
    @testset "@structured round-trips: export → validate → reload is trajectory-identical" begin
        function det_model()   # deterministic npv so DSL and JSON runs are bit-identical
            net = @reaction_network begin
                @deterministic(1.0),
                    ∅ --> @structured(:Project, phase = :Phase1, npv = 100.0, born = @t()),
                    name => genesis
                @deterministic(1.0),
                    @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
                    name => adv12, cycletime => 1.0, probability => 1.0
            end
            RDX.register_structured_species!(net, :Project)
            @prob_meta net tspan = 5 dt = 1.0
            net
        end
        p_dsl = ReactionNetworkProblem(det_model(); seed = 7, registry = GEN_REGISTRY)
        simulate(p_dsl)

        exported = RDX.to_json_model(p_dsl)                      # live model → eval-free JSON
        doc = JSON.parse(exported)
        # the genesis reactant emitted a typed structured{kind, fields} row (no host Expr)
        srow = only(r for r in doc["reactants"] if haskey(r, "structured"))
        @test srow["structured"]["kind"] == "Project"
        @test Set(f["name"] for f in srow["structured"]["fields"]) == Set(["phase", "npv", "born"])

        @test isempty(RDX.validate(doc; registry = GEN_REGISTRY))   # validates clean

        p_json = RDX.from_json_model(exported; seed = 7, registry = GEN_REGISTRY)
        simulate(p_json)
        @test phases_of(p_json) == phases_of(p_dsl)              # identical trajectory
        @test sort([t.born for t in livetokens(p_json)]) == sort([t.born for t in livetokens(p_dsl)])
        # idempotent: re-exporting the reload reproduces the same document (round-trip fixed point)
        @test JSON.parse(RDX.to_json_model(p_json)) == doc
    end

    # ── validate flags a genesis whose kind is not in the registry / not structured ───────
    @testset "validate rejects a named @structured with an unknown kind" begin
        net = @reaction_network begin
            @deterministic(1.0),
                ∅ --> @structured(:Project, phase = :Phase1, npv = 100.0, born = @t()),
                name => genesis
        end
        RDX.register_structured_species!(net, :Project)
        @prob_meta net tspan = 2 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1, registry = GEN_REGISTRY)
        doc = JSON.parse(RDX.to_json_model(p))
        for r in doc["reactants"]
            haskey(r, "structured") && (r["structured"]["kind"] = "Ghost")   # dangling kind
        end
        diags = RDX.validate(doc; registry = GEN_REGISTRY)
        @test !isempty(diags)
        @test any(d -> occursin("Ghost", string(d)), diags)
    end

    # ── the RAW `@structured(Ctor(…))` form is REJECTED at construction (removed) ──────────
    @testset "raw @structured(Ctor(…)) is rejected at construction — only the named form exists" begin
        # The raw constructor form was removed: it was the sole reactant construct that could not
        # round-trip eval-free, so forbidding it makes eval-free serialization TOTAL. The rejection
        # is at CONSTRUCTION (recursively_find_reactants!, create.jl) — fail fast at model build,
        # not deep in a simulation — and the message points at the named replacement.
        err = try
            @reaction_network begin
                @deterministic(1.0),
                    ∅ --> @structured(GenProjectToken(:Phase1, 100.0, @t())),   # inline host ctor
                    name => genesis
            end
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing                       # construction threw
        @test occursin("named form", err)           # message names the supported form
        @test occursin("@structured(:Kind", err)     # …and shows its shape
    end

    # ── the two-arg raw form `@structured(token, species)` is likewise rejected ────────────
    @testset "raw @structured(token, species) two-arg form is also rejected at construction" begin
        err = try
            @reaction_network begin
                @deterministic(1.0),
                    ∅ --> @structured(GenProjectToken(:Phase1, 100.0, @t()), :Project),
                    name => genesis
            end
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("named form", err)
    end
end
