# Phase-1 Stage E semantic tests — typed ExprNode IR + JSON serialization (ADR 0005).
#
# Grows step by step (E1..E9). The IR is eval-free: a model authored as a typed tree (or as JSON)
# lowers via to_expr to EXACTLY the Expr the DSL produces, then through the unchanged compiler.
# Round-trip + validate + the BD model as JSON are the acceptance bars.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

@testset "Typed ExprNode IR + JSON serialization (ADR 0005)" begin

    # ── E1: to_expr lowers to the exact DSL Expr; from_expr is the structural inverse ──
    @testset "E1: ExprNode to_expr / from_expr (Const / Ref / Call arithmetic core)" begin
        # scalar leaves
        @test RDX.to_expr(RDX.Const(0.3)) === 0.3
        @test RDX.to_expr(RDX.Const(5)) === 5
        @test RDX.to_expr(RDX.Const(true)) === true
        @test RDX.to_expr(RDX.Const(:Phase2)) == QuoteNode(:Phase2)   # literal symbol ⇒ QuoteNode
        @test RDX.to_expr(RDX.NodeRef(:species, :X)) === :X                # bare symbol
        @test RDX.to_expr(RDX.NodeRef(:param, :beta)) === :beta

        # arithmetic: 0.3 * beta
        n = RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        @test RDX.to_expr(n) == :(0.3 * beta)

        # comparison nested under short-circuit && uses the :&& HEAD, not a :call
        g = RDX.Call(:&&, [RDX.Call(:>, [RDX.NodeRef(:species, :X), RDX.Const(0)]), RDX.Const(true)])
        gx = RDX.to_expr(g)
        @test gx.head == :(&&)
        @test gx.args[1] == :(X > 0)

        # unary !
        @test RDX.to_expr(RDX.Call(:!, [RDX.Const(true)])) == :(!true)

        # from_expr round-trip on arithmetic (species/param sets classify bare symbols)
        rt = RDX.from_expr(:(0.3 * beta); species = Set{Symbol}(), params = Set([:beta]))
        @test rt == RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        # to_expr ∘ from_expr == identity on a compound arithmetic/comparison expr
        ex = :((X + 2) * beta > 0)
        @test RDX.to_expr(RDX.from_expr(ex; species = Set([:X]), params = Set([:beta]))) == ex

        # op outside the whitelist is rejected at to_expr
        @test_throws Exception RDX.to_expr(RDX.Call(:foncall, [RDX.Const(1)]))
    end

    # ── E2: node_to_dict/from_dict + model envelope + from_json_model round-trip ────────
    @testset "E2: ExprNode JSON dict round-trip" begin
        for n in (
            RDX.Const(0.3), RDX.Const(5), RDX.Const(true), RDX.Const(:Phase2),
            RDX.NodeRef(:param, :beta),
            RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)]),
            RDX.Call(:&&, [RDX.Call(:>, [RDX.NodeRef(:species, :X), RDX.Const(0)]), RDX.Const(true)]),
        )
            @test RDX.node_from_dict(RDX.node_to_dict(n)) == n
        end
    end

    @testset "E2: from_json_model builds a runnable model; params are numbers, not eval'd" begin
        json = """
        {
          "rd_format": "reactive-dynamics-model", "version": "1.0",
          "meta": { "tspan": 5.0, "dt": 1.0 },
          "params": [ { "name": "beta", "value": 0.4 } ],
          "species": [ { "name": "A", "init": 100 }, { "name": "B", "init": 0 } ],
          "transitions": [
            { "id": "t1", "name": "t1", "rate": 1.0, "rate_mode": "deterministic",
              "cycletime": 0.0, "prob_of_success": 1.0 } ],
          "reactants": [
            { "transition": "t1", "species": "A", "side": "lhs", "stoich": 1 },
            { "transition": "t1", "species": "B", "side": "rhs", "stoich": 1 } ]
        }
        """
        p = RDX.from_json_model(json; seed = 1)
        @test p.p[:beta] === 0.4                       # JSON number, not an eval'd string
        @test p.u[RDX.find_index(:A, p)] == 100.0      # species init
        simulate(p)
        @test p.sol.B[end] > 0.0                       # the A --> B transition fired
        # the assembled reaction line is the expected Expr
        @test p.acs[1, :trans] == :(A → B)
    end

    @testset "E2: model_to_dict ∘ build_acs round-trips on the parsed Dict" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{},
          "params":[{"name":"k","value":0.5}],
          "species":[{"name":"A","init":10},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic"}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","stoich":1},
                       {"transition":"t1","species":"B","side":"rhs","stoich":1}] }
        """
        import JSON
        d = JSON.parse(json)
        acs = RDX.build_acs_from_dict(d)
        # params + species survive the round-trip through the acset
        back = RDX.model_to_dict(acs)
        @test any(pr -> pr["name"] == "k" && pr["value"] == 0.5, back["params"])
        @test Set(sp["name"] for sp in back["species"]) == Set(["A", "B"])
    end

end
