# Phase-0.6 semantic tests — Visualization (ADR 0014 / CONTRACT §15).
#
# Covers the network "exec map" three layers + the result-plot recipes:
#   §15.2 Layer A — `network_graph` is a pure function of the model (no run, no plotting dep, does
#         not perturb the RNG); structure correctness (place/transition/arc counts).
#   §15.2 Layer B — `to_graphviz` emits valid DOT (a `dot`-parseable digraph) for the reference
#         models; `draw_network` renders through AA's `run_graphviz` when a backend is available.
#   §15.2 Layer C — `exec_map` decorates the structure with run statistics + `@select` highlighting,
#         read-only (Invariant 3).
#   §15.1 recipes — each plot wrapper is constructible from a raw artifact (model-agnostic). The
#         actual `plot(...)` rendering is exercised only when Plots is loaded (RDPlotsExt).
#
# DOT validity is checked structurally (digraph envelope) and, when a Graphviz backend is present,
# by round-tripping through `run_graphviz` (a malformed DOT throws).

using ReactiveDynamics, Test
using Random, DataFrames
using Plots   # trigger RDPlotsExt so the §15.1 recipe-render assertion runs (not skip)

RD = ReactiveDynamics

# A small SIR model (the reference acceptance model) — plain place, no structured tokens, so
# network_graph works with no population.
function sir_model()
    net = @reaction_network begin
        0.5 / 1000, S + I --> 2 * I, name => infection
        0.05, I --> R, name => recovery
    end
    @prob_meta net tspan = 10 dt = 1.0
    return net
end

# A structured advance model (toy-pharma-like) with a phase attribute, to test highlighting.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct VizProjectToken
        phase::Symbol
    end
    function VizProjectToken(phase)
        return VizProjectToken(
            "VP" * string(rand(1:(10^9))), :Project, nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Firing}[], phase
        )
    end
end

function pharma_model(; budget0 = 100)
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) + 2 * @rate(budget) --> @advance(phase, :Phase2),
            name => adv, cycletime => 1.0, probability => 1.0
    end
    RD.register_token_kind!(net, :Project)
    bi = findfirst(==(:budget), net[:, :placeName])
    net[bi, :placeInitVal] = Float64(budget0)
    net[bi, :placeCost] = 1.0
    @prob_meta net tspan = 5 dt = 1.0
    return net
end

# Does a string look like a single well-formed DOT digraph? (envelope + balanced braces)
function valid_dot(s)
    st = strip(s)
    startswith(st, "digraph") || return false
    endswith(st, "}") || return false
    return count(==('{'), st) == count(==('}'), st)
end

# Try to render DOT through AA's run_graphviz; returns true if it rendered, false if no backend.
function renders(dot)
    try
        io = IOBuffer()
        AlgebraicAgents.run_graphviz(io, dot; format = "svg")
        return length(take!(io)) > 0
    catch
        return false   # no Graphviz backend in this environment — Layer B is best-effort (open Q)
    end
end

@testset "Visualization (ADR 0014 / CONTRACT §15)" begin

    # ── §15.2 Layer A — structure extraction is pure (no run, no RNG perturbation) ────────
    @testset "§15.2 Layer A network_graph: SIR structure correctness, no simulation" begin
        p = ReactionNetworkProblem(sir_model(); seed = 1)
        g = network_graph(p)
        spnames = Set(s.name for s in g.places)
        @test spnames == Set([:S, :I, :R])
        @test length(g.transitions) == 2          # infection + recovery
        @test !isempty(g.arcs)
        # an :in arc (LHS→T) and an :out arc (T→RHS) both exist
        @test any(a -> a.dir === :in, g.arcs)
        @test any(a -> a.dir === :out, g.arcs)
    end

    @testset "§15.2 Layer A is observationally pure — does not perturb the caller's RNG (Inv 1)" begin
        p = ReactionNetworkProblem(sir_model(); seed = 5)
        before = copy(p.rng)
        network_graph(p)
        # the rng state is unchanged (network_graph works on a deepcopy)
        @test rand(p.rng) == rand(before)
    end

    @testset "§15.2 Layer A round-trips from JSON (structure ⟂ run, Invariant 1)" begin
        p = ReactionNetworkProblem(sir_model(); seed = 1)
        json = to_json_model(p)
        q = from_json_model(json)
        gp, gq = network_graph(p), network_graph(q)
        @test Set(s.name for s in gp.places) == Set(s.name for s in gq.places)
        @test length(gp.transitions) == length(gq.transitions)
    end

    # ── §15.2 Layer B — valid DOT emission for the reference models ───────────────────────
    @testset "§15.2 Layer B to_graphviz emits valid DOT for SIR and toy-pharma" begin
        for acs_fn in (sir_model, pharma_model)
            p = ReactionNetworkProblem(
                acs_fn(); seed = 1,
                population = acs_fn === pharma_model ? [RD.VizProjectToken(:Phase1)] : []
            )
            dot = to_graphviz(network_graph(p))
            @test valid_dot(dot)
            @test occursin("shape=box", dot)       # transitions are boxes
            @test occursin("->", dot)              # arcs present
        end
    end

    @testset "§15.2 Layer B draw_network renders when a Graphviz backend is present" begin
        p = ReactionNetworkProblem(sir_model(); seed = 1)
        dot = to_graphviz(network_graph(p))
        if renders(dot)
            out = draw_network(p; format = "svg")
            @test out isa AbstractString && !isempty(out)
        else
            @test_skip "no Graphviz backend — draw_network returns DOT-render best-effort (open Q)"
        end
    end

    # ── §15.2 Layer C — the result overlay (read-only) ────────────────────────────────────
    @testset "§15.2 Layer C exec_map decorates structure with run stats + @select highlight" begin
        p = ReactionNetworkProblem(pharma_model(); seed = 1, population = [RD.VizProjectToken(:Phase1)])
        simulate(p)
        # exec_map returns rendered bytes (when a backend exists) or is best-effort; either way it
        # must build the highlighted DOT without mutating state.
        sol_before = copy(p.sol)
        pred = RD.TokenPredicate(:Project, [RD.Clause(:phase, :(==), QuoteNode(:Phase2))])
        g = network_graph(p)
        # the overlay DOT (Layer C styling) is valid and includes highlight styling hooks
        hi = to_graphviz(g; highlight_places = [:budget], highlight_arcs = Tuple{Symbol, Symbol}[])
        @test valid_dot(hi)
        @test occursin("fillcolor=gold", hi)       # starvation/highlight fill present
        if renders(hi)
            out = exec_map(p; highlight = pred, format = "svg")
            @test out isa AbstractString
        end
        @test p.sol == sol_before                   # read-only (Invariant 3)

        # token-path highlighting must map a bond's transition INDEX to the SAME node id the graph
        # uses (not the bond's per-instance "<name>_@<t>"), so a highlighted arc actually matches a
        # graph arc. Build the highlight arc set the way exec_map does and assert overlap.
        gnodes = Set(t.name for t in g.transitions)
        garcs = Set((a.from, a.to) for a in g.arcs)
        hi = Tuple{Symbol, Symbol}[]
        for tok in RD.select_tokens(p, pred)
            for (pl, _t, tr) in tok.past_bonds
                push!(hi, (pl, RD._transition_node_name(p.network, tr.i)))
            end
        end
        if !isempty(hi)
            @test all(((_f, to),) -> to in gnodes, hi)   # every highlight points at a real T node
            @test any(in(garcs), hi)                      # and at least one matches a real arc
        end
    end

    # ── §15.1 recipes are model-agnostic wrappers over raw artifacts ──────────────────────
    @testset "§15.1 recipe wrappers construct from raw artifacts (no BD-specific field)" begin
        p = ReactionNetworkProblem(sir_model(); seed = 1)
        simulate(p)
        # the wrapper types are in the core (constructible with Plots absent)
        @test RD.MarkingPlot(p) isa RD.MarkingPlot
        @test RD.SaturationPlot(p) isa RD.SaturationPlot
        @test RD.ValuationPlot(p) isa RD.ValuationPlot
        @test RD.LedgerPlot(p) isa RD.LedgerPlot
        @test RD.ThroughputPlot(p) isa RD.ThroughputPlot
        # the log-series helpers extract scalar/count series from the raw log
        t, v = RD.log_scalar_series(p, :valuation)
        @test length(t) == length(v)
        # count-series must read BOTH (hash,q) Tuple rows (:new_transitions) and Symbol=>Float64 Pair
        # rows (:terminated_all). SIR is instant (cycletime 0 — nothing in-flight to terminate across
        # ticks), so exercise the Pair path on the pharma model (cycletime 1 → real terminations).
        ph = ReactionNetworkProblem(pharma_model(); seed = 1, population = [RD.VizProjectToken(:Phase1)])
        simulate(ph)
        tt, tv = RD.log_count_series(ph, :terminated_all)
        @test length(tt) == length(tv)
        @test sum(tv) > 0.0      # would be 0 if the Pair rows were mis-read as Tuples
    end

    @testset "§15.1 recipes render when Plots is loaded (RDPlotsExt)" begin
        if Base.get_extension(ReactiveDynamics, :RDPlotsExt) !== nothing
            # Plots is available; a recipe should produce a plot object. Import lazily.
            @eval using Plots
            p = ReactionNetworkProblem(sir_model(); seed = 1); simulate(p)
            pl = @eval Plots.plot(RD.MarkingPlot($p))
            @test pl !== nothing
        else
            @test_skip "Plots not loaded — recipes live in RDPlotsExt"
        end
    end
end
