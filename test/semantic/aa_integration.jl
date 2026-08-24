# Phase-1 semantic tests — AlgebraicAgents integration & external coupling (ADR 0012, CONTRACT §13).
#
# Covers, in BOTH coupling directions:
#   (A) OUTBOUND — RD as a readable hierarchy node: observables(rd) lists place + named
#       observables; getobservable(rd, name|i) resolves a place count and a named observable's
#       sampled value; _getparameters/_setparameters! round-trip state.p (param-only, index-safe).
#   (B) INBOUND  — the ExternalRef leaf round-trips (node_to_dict/from_dict and to_expr/from_expr);
#       validate rule 8 flags an undeclared inputs[] port and passes a declared one; the _prestep!
#       latch makes a coupled read DETERMINISTIC and ONE-TICK-LAGGED (the latched value equals the
#       source's previous-tick-boundary projection and is identical across read sites in a tick,
#       independent of AA sibling step order — the §4-D4 hazard the latch closes); _reinit!
#       restores external_inputs to the declared defaults.
#
# Two tiers (per the suite convention): all tests here run against the BUILT engine — they are
# real @tests of the ADR-0012 surface, not @test_broken pins.
#
# NB Julia 1.12: `const` is a syntax error in a testset body, so the module alias and the helper
# agent type are defined at FILE scope; testsets use the alias `RD`.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
import MacroTools
import JSON

const RD = ReactiveDynamics

# A minimal sibling SOURCE agent (a stand-in macro/market signal), defined at file scope. It
# exports one observable `signal` and advances it deterministically each step, so a wire from it
# into an RD port lets us assert the one-tick coupling lag exactly. It implements only the AA
# surface the test needs: observables/getobservable (outbound), _step!/_projected_to (stepping).
@aagent struct AATestSource
    val::Float64
    dt::Float64
    t::Float64
end
# The custom constructor MUST type its first arg `::AbstractString` so dispatch picks it over the
# @aagent-generated positional `AATestSource(name, args...)` (else the macro ctor intercepts it).
AATestSource(name::AbstractString, v0::Real, dt::Real) =
    AATestSource(name, Float64(v0), Float64(dt), 0.0)
AlgebraicAgents.observables(s::AATestSource) = [:signal]
AlgebraicAgents.getobservable(s::AATestSource, ::Union{Symbol, AbstractString}) = s.val
AlgebraicAgents.getobservable(s::AATestSource, ::Int) = s.val
AlgebraicAgents._step!(s::AATestSource) = (s.val += 1.0; s.t += s.dt; s.t)
AlgebraicAgents._projected_to(s::AATestSource) = s.t > 5.0 ? true : s.t

# A small classical RD net: a `cash` pool that a deterministic source fills at a rate read from
# the external port `ext_rate`. Used by the latch tests. `cycletime`=0 ⇒ clean tick boundaries.
const _COUPLED_JSON = """
{ "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
  "params":[],
  "inputs":[{"port":"ext_rate","default":{"node":"const","value":0.0}}],
  "places":[{"name":"A","init":0}],
  "transitions":[{"id":"t1","name":"t1","rate":{"node":"externalref","port":"ext_rate"},"rate_mode":"deterministic"}],
  "arcs":[{"transition":"t1","side":"rhs","place":"A","multiplicity":1}] }
"""

@testset "AlgebraicAgents integration & external coupling (ADR 0012)" begin

    # ── (A) OUTBOUND read surface — RD as a readable hierarchy node ─────────────────────
    @testset "A: observables/getobservable resolve place + named observables" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[{"name":"beta","value":0.4}],
          "places":[{"name":"A","init":100},{"name":"B","init":7}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic",
                          "cycletime":0.0,"prob_of_success":1.0}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        p = RD.from_json_model(json; seed = 1)

        # observables() lists every place name (named observables would follow; none here).
        obs = AlgebraicAgents.observables(p)
        @test :A in obs && :B in obs
        @test obs[1:2] == [:A, :B]                       # place first, in declared order

        # getobservable by name = the live place count; identical via Symbol and String.
        @test AlgebraicAgents.getobservable(p, :A) == 100.0
        @test AlgebraicAgents.getobservable(p, "A") == 100.0
        @test AlgebraicAgents.getobservable(p, :B) == 7.0
        # getobservable by Int indexes the SAME ordered list (Invariant 1 / §A).
        @test AlgebraicAgents.getobservable(p, 1) == AlgebraicAgents.getobservable(p, :A)

        # an unknown name is a hard error — a diagnostic, never AA's silent @error fall-through.
        @test_throws Exception AlgebraicAgents.getobservable(p, :nonexistent)

        # outbound reads track the live trajectory after a step.
        simulate(p)
        @test AlgebraicAgents.getobservable(p, :B) == p.u[RD.find_index(:B, p)]
    end

    @testset "A: a NAMED observable surfaces in observables() and reads its sampled value" begin
        # a lean-explicit token/derived aggregate is exported as a named observable (§A open Q).
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":4.0,"dt":1.0},
          "params":[],
          "places":[{"name":"A","init":50}],
          "transitions":[{"id":"t1","name":"t1","rate":0.0,"rate_mode":"deterministic"}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1}],
          "observables":[{"name":"a_level","every":1.0,
                          "range":[{"weight":1.0,"value":{"node":"ref","kind":"place","name":"A"}}]}] }
        """
        p = RD.from_json_model(json; seed = 1)
        @test :a_level in AlgebraicAgents.observables(p)   # named observable is exported
        simulate(p)                                        # resample fills .sampled
        # reads the observable's last-sampled value (== A's level, since the obs is just A).
        @test AlgebraicAgents.getobservable(p, :a_level) == p.observables[:a_level].sampled
    end

    @testset "A: _getparameters / _setparameters! round-trip state.p (param-only, index-safe)" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":3.0,"dt":1.0},
          "params":[{"name":"beta","value":0.4},{"name":"gamma","value":2.0}],
          "places":[{"name":"A","init":1}],
          "transitions":[{"id":"t1","name":"t1","rate":0.0,"rate_mode":"deterministic"}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1}] }
        """
        p = RD.from_json_model(json; seed = 1)
        ps = AlgebraicAgents._getparameters(p)
        @test ps === p.p                                   # exposes the live param dict
        @test ps[:beta] == 0.4 && ps[:gamma] == 2.0

        # a param-only patch merges into state.p; structure (place count) is untouched.
        nS_before = RD.nrows(p, :S)
        AlgebraicAgents._setparameters!(p, Dict(:beta => 0.9, :delta => 3.0))
        @test p.p[:beta] == 0.9                            # existing param overwritten
        @test p.p[:delta] == 3.0                           # new param added
        @test p.p[:gamma] == 2.0                           # untouched param preserved
        @test RD.nrows(p, :S) == nS_before                # no structural change (Invariant 5)
    end

    # ── (B) INBOUND — ExternalRef leaf round-trips ──────────────────────────────────────
    @testset "B: ExternalRef round-trips through node dict and to_expr/from_expr" begin
        n = RD.ExternalRef(:ext_rate)
        # JSON node dict round-trip (the eval-free serialization boundary).
        d = RD.node_to_dict(n)
        @test d["node"] == "externalref" && d["port"] == "ext_rate"
        @test RD.node_from_dict(d) == n

        # to_expr lowers to the per-tick buffer read; from_expr is the structural inverse.
        ex = RD.to_expr(n)
        @test ex == :(state.external_inputs[:ext_rate])
        # MacroTools.striplines: the lowered Expr carries no LineNumberNodes here, but compare
        # robustly anyway (the suite convention for Expr equality).
        @test MacroTools.striplines(RD.from_expr(ex) |> RD.to_expr) ==
            MacroTools.striplines(ex)
        @test RD.from_expr(ex) == n                        # exact leaf recovery

        # structural ==/hash work (so dedup / round-trip tests compare by value).
        @test RD.ExternalRef(:p) == RD.ExternalRef(:p)
        @test RD.ExternalRef(:p) != RD.ExternalRef(:q)
        @test hash(RD.ExternalRef(:p)) == hash(RD.ExternalRef(:p))

        # ExternalRef nested in a Call (a rate `spend * ExternalRef(ext_rate)`) round-trips.
        call = RD.Call(:*, [RD.NodeRef(:param, :spend), RD.ExternalRef(:ext_rate)])
        @test RD.node_from_dict(RD.node_to_dict(call)) == call
    end

    # ── (B) validate rule 8 — every ExternalRef port must be a declared inputs[] port ───
    @testset "B: validate rule 8 flags an undeclared port and passes a declared one" begin
        base = JSON.parse(
            """
            { "meta":{"tspan":5.0,"dt":1.0},"params":[],
              "inputs":[{"port":"ext_rate","default":{"node":"const","value":0.0}}],
              "places":[{"name":"A","init":0}],
              "transitions":[{"id":"t1","rate":{"node":"externalref","port":"ext_rate"},"rate_mode":"deterministic"}],
              "arcs":[{"transition":"t1","place":"A","side":"rhs","multiplicity":1}] }
            """
        )
        # the declared port validates clean.
        @test isempty(RD.validate(base))

        # an UNDECLARED port in the rate is a rule-8 diagnostic.
        bad = deepcopy(base); bad["transitions"][1]["rate"]["port"] = "undeclared_port"
        diags = RD.validate(bad)
        @test any(d -> occursin("not a declared inputs[] port", d.msg), diags)
        # from_json_model gates on validation.
        @test_throws Exception RD.from_json_model(JSON.json(bad))

        # an ExternalRef in a RULE GUARD against a declared port also validates clean,
        # and an undeclared one in the guard is flagged (the "guard value" context, §B2).
        guarded = deepcopy(base)
        guarded["inputs"] = [Dict("port" => "sentiment", "default" => Dict("node" => "const", "value" => 0.0))]
        guarded["transitions"][1]["rate"] = 0.0
        guarded["rules"] = [
            Dict(
                "id" => "lever", "fire_mode" => "once",
                "guard" => Dict(
                    "node" => "call", "op" => ">",
                    "args" => [
                        Dict("node" => "externalref", "port" => "sentiment"),
                        Dict("node" => "const", "value" => 0.5),
                    ]
                ),
                "action" => Dict(
                    "verb" => "set_marking", "name" => "A", "mode" => "inc",
                    "value" => Dict("node" => "const", "value" => 1)
                )
            ),
        ]
        @test isempty(RD.validate(guarded))
        bad_guard = deepcopy(guarded); bad_guard["rules"][1]["guard"]["args"][1]["port"] = "ghost"
        @test any(d -> occursin("not a declared inputs[] port", d.msg), RD.validate(bad_guard))
    end

    # ── (B) the _prestep! latch — deterministic, one-tick-lagged, read-site-identical ──
    @testset "B: _prestep! latches the source's previous-boundary value (one-tick Jacobi lag)" begin
        rd = RD.from_json_model(_COUPLED_JSON; seed = 1)
        src = AATestSource("src", 10.0, 1.0)               # signal starts at 10, +1 each step
        root = FreeAgent("root")
        entangle!(root, src); entangle!(root, rd)
        add_wire!(root; from = src, to = rd, from_var_name = "signal", to_var_name = "ext_rate")

        # before any step the buffer holds the declared default (pre-wire fallback, §B3).
        @test rd.external_inputs[:ext_rate] == 0.0

        # one root step: _prestep! (whole-hierarchy first phase) latches src.val = 10 BEFORE src
        # steps to 11 — the previous-tick-boundary projection (Invariant 2), a one-tick lag.
        step!(root)
        @test rd.external_inputs[:ext_rate] == 10.0
        @test src.val == 11.0                              # source advanced this same tick
        # the RD rate read that latched 10.0 ⇒ A advanced by 10 over the tick.
        @test rd.sol.A[end] == 10.0

        step!(root)
        @test rd.external_inputs[:ext_rate] == 11.0         # next tick latches the new boundary
        @test src.val == 12.0
    end

    @testset "B: a coupled run is reproducible under (hierarchy, seed)" begin
        function build()
            rd = RD.from_json_model(_COUPLED_JSON; seed = 7)
            src = AATestSource("src", 5.0, 1.0)
            root = FreeAgent("root")
            entangle!(root, src); entangle!(root, rd)
            add_wire!(root; from = src, to = rd, from_var_name = "signal", to_var_name = "ext_rate")
            return root, rd
        end
        r1, rd1 = build(); simulate(r1)
        r2, rd2 = build(); simulate(r2)
        @test rd1.sol == rd2.sol                            # identical coupled trajectory (Invariant 3)
    end

    @testset "B: every ExternalRef read site in a tick agrees (identical buffered value)" begin
        # two transitions both reading the SAME port in one tick must see the SAME latched value
        # (Invariant 2: read-site-identical within a tick, independent of sibling step order).
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":3.0,"dt":1.0},
          "params":[],
          "inputs":[{"port":"ext_rate","default":{"node":"const","value":0.0}}],
          "places":[{"name":"A","init":0},{"name":"B","init":0}],
          "transitions":[
            {"id":"ta","name":"ta","rate":{"node":"externalref","port":"ext_rate"},"rate_mode":"deterministic"},
            {"id":"tb","name":"tb","rate":{"node":"externalref","port":"ext_rate"},"rate_mode":"deterministic"}],
          "arcs":[{"transition":"ta","side":"rhs","place":"A","multiplicity":1},
                 {"transition":"tb","side":"rhs","place":"B","multiplicity":1}] }
        """
        rd = RD.from_json_model(json; seed = 1)
        src = AATestSource("src", 4.0, 1.0)
        root = FreeAgent("root")
        entangle!(root, src); entangle!(root, rd)
        add_wire!(root; from = src, to = rd, from_var_name = "signal", to_var_name = "ext_rate")
        step!(root)
        # both transitions read the same latched 4.0 ⇒ A and B advance identically.
        @test rd.sol.A[end] == rd.sol.B[end] == 4.0
    end

    # ── (B) reinit! restores external_inputs to the declared defaults ──────────────────
    @testset "B: _reinit! restores external_inputs to the declared defaults" begin
        rd = RD.from_json_model(_COUPLED_JSON; seed = 1)
        src = AATestSource("src", 10.0, 1.0)
        root = FreeAgent("root")
        entangle!(root, src); entangle!(root, rd)
        add_wire!(root; from = src, to = rd, from_var_name = "signal", to_var_name = "ext_rate")
        step!(root); step!(root)
        @test rd.external_inputs[:ext_rate] != 0.0          # a stale latched value is present

        AlgebraicAgents._reinit!(rd)
        # reinit drops the latched wire value and restores the pre-wire default (§B3 / §4 D7).
        @test rd.external_inputs[:ext_rate] == 0.0
        @test rd.external_inputs == rd.external_input_defaults
    end

end
