# Phase-1 Stage E semantic tests — typed ExprNode IR + JSON serialization (ADR 0005).
#
# Grows step by step (E1..E9). The IR is eval-free: a model authored as a typed tree (or as JSON)
# lowers via to_expr to EXACTLY the Expr the DSL produces, then through the unchanged compiler.
# Round-trip + validate + the BD model as JSON are the acceptance bars.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

# A structured token kind for the E4 phase-pipeline tests (defined in RD scope, the @register idiom).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct SerProjectToken
        phase::Symbol
        npv::Float64
    end
    function SerProjectToken(phase, npv)
        return SerProjectToken(
            "SP" * string(rand(1:(10^9))), :Project, nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Firing}[], phase, npv
        )
    end
end

@testset "Typed ExprNode IR + JSON serialization (ADR 0005)" begin

    # ── E1: to_expr lowers to the exact DSL Expr; from_expr is the structural inverse ──
    @testset "E1: ExprNode to_expr / from_expr (Const / Ref / Call arithmetic core)" begin
        # scalar leaves
        @test RDX.to_expr(RDX.Const(0.3)) === 0.3
        @test RDX.to_expr(RDX.Const(5)) === 5
        @test RDX.to_expr(RDX.Const(true)) === true
        @test RDX.to_expr(RDX.Const(:Phase2)) == QuoteNode(:Phase2)   # literal symbol ⇒ QuoteNode
        @test RDX.to_expr(RDX.NodeRef(:place, :X)) === :X                # bare symbol
        @test RDX.to_expr(RDX.NodeRef(:param, :beta)) === :beta

        # arithmetic: 0.3 * beta
        n = RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        @test RDX.to_expr(n) == :(0.3 * beta)

        # comparison nested under short-circuit && uses the :&& HEAD, not a :call
        g = RDX.Call(:&&, [RDX.Call(:>, [RDX.NodeRef(:place, :X), RDX.Const(0)]), RDX.Const(true)])
        gx = RDX.to_expr(g)
        @test gx.head == :(&&)
        @test gx.args[1] == :(X > 0)

        # unary !
        @test RDX.to_expr(RDX.Call(:!, [RDX.Const(true)])) == :(!true)

        # from_expr round-trip on arithmetic (place/param sets classify bare symbols)
        rt = RDX.from_expr(:(0.3 * beta); places = Set{Symbol}(), params = Set([:beta]))
        @test rt == RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        # to_expr ∘ from_expr == identity on a compound arithmetic/comparison expr
        ex = :((X + 2) * beta > 0)
        @test RDX.to_expr(RDX.from_expr(ex; places = Set([:X]), params = Set([:beta]))) == ex

        # op outside the whitelist is rejected at to_expr
        @test_throws Exception RDX.to_expr(RDX.Call(:foncall, [RDX.Const(1)]))
    end

    # ── E2: node_to_dict/from_dict + model envelope + from_json_model round-trip ────────
    @testset "E2: ExprNode JSON dict round-trip" begin
        for n in (
                RDX.Const(0.3), RDX.Const(5), RDX.Const(true), RDX.Const(:Phase2),
                RDX.NodeRef(:param, :beta),
                RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)]),
                RDX.Call(:&&, [RDX.Call(:>, [RDX.NodeRef(:place, :X), RDX.Const(0)]), RDX.Const(true)]),
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
          "places": [ { "name": "A", "init": 100 }, { "name": "B", "init": 0 } ],
          "transitions": [
            { "id": "t1", "name": "t1", "rate": 1.0, "rate_mode": "deterministic",
              "cycletime": 0.0, "prob_of_success": 1.0 } ],
          "arcs": [
            { "transition": "t1", "place": "A", "side": "lhs", "multiplicity": 1 },
            { "transition": "t1", "place": "B", "side": "rhs", "multiplicity": 1 } ]
        }
        """
        p = RDX.from_json_model(json; seed = 1)
        @test p.p[:beta] === 0.4                       # JSON number, not an eval'd string
        @test p.u[RDX.find_index(:A, p)] == 100.0      # place init
        simulate(p)
        @test p.sol.B[end] > 0.0                       # the A --> B transition fired
        # the assembled reaction line is the expected Expr
        @test p.network[1, :trans] == :(A → B)
    end

    # ── E3: Sample / TimeRef / Choose + rate_mode Poisson wrapping/unwrapping ──────────
    @testset "E3: Sample lowers to rand(state.rng, Dist(...)); DIST_WHITELIST enforced" begin
        s = RDX.Sample(:Poisson, [RDX.Const(0.3)])
        @test RDX.to_expr(s) == :(rand(state.rng, Poisson(0.3)))
        @test_throws Exception RDX.to_expr(RDX.Sample(:NotADist, [RDX.Const(1.0)]))
        @test RDX.to_expr(RDX.TimeRef()).args[1] == Symbol("@t")
    end

    @testset "E3: rate_mode poisson lowers to the documented expand_rate Expr; unwrap inverts" begin
        # bare intensity 0.3 * beta, poisson mode ⇒ the wrapped expand_rate shape
        bare = RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        lowered = RDX.lower_rate(bare, :poisson)
        @test lowered == :(rand(state.rng, Poisson(max(state.dt * (0.3 * beta), 0))))
        # unwrap recovers (bare, :poisson)
        node, mode = RDX.rate_from_expr(lowered; params = Set([:beta]))
        @test mode == :poisson
        @test node == bare
        # deterministic mode is bare, and unwraps to (bare, :deterministic)
        @test RDX.lower_rate(RDX.Const(2.0), :deterministic) === 2.0
        n2, m2 = RDX.rate_from_expr(2.0)
        @test m2 == :deterministic && n2 == RDX.Const(2.0)
    end

    @testset "E3: a Choose-rate model simulates deterministically under a fixed seed" begin
        # @choose branches between two rates; recursively_choose consumes it through state.rng
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":10.0,"dt":1.0},
          "params":[],
          "places":[{"name":"A","init":1000},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":3.0,"rate_mode":"deterministic",
                          "cycletime":0.0,"prob_of_success":1.0}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        p1 = RDX.from_json_model(json; seed = 42); simulate(p1)
        p2 = RDX.from_json_model(json; seed = 42); simulate(p2)
        @test p1.sol == p2.sol           # same (model, seed) ⇒ identical
    end

    # ── E4: arcs[] → reaction-line :trans Expr (multiplicity, modality, @select, @advance) ──
    @testset "E4: multi-LHS + integer multiplicity assembles to the runtime-parsed reaction line" begin
        rs = [
            Dict("transition" => "t", "place" => "X", "side" => "lhs", "multiplicity" => 1),
            Dict("transition" => "t", "place" => "Y", "side" => "lhs", "multiplicity" => 2),
            Dict("transition" => "t", "place" => "Z", "side" => "rhs", "multiplicity" => 1),
        ]
        line = RDX.assemble_reaction_line(rs)
        @test line == :((X + 2Y) → Z)
    end

    @testset "E4: LHS modality macros (@conserved/@rate) are emitted per the 3-axis modality" begin
        rs = [
            Dict(
                "transition" => "t", "place" => "scientist", "side" => "lhs", "multiplicity" => 3,
                "modality" => Dict("allocation" => "upfront", "return" => "conserved", "blocking" => "block")
            ),
            Dict(
                "transition" => "t", "place" => "budget", "side" => "lhs", "multiplicity" => 1,
                "modality" => Dict("allocation" => "perstep", "return" => "consumed", "blocking" => "block")
            ),
            Dict("transition" => "t", "place" => "out", "side" => "rhs", "multiplicity" => 1),
        ]
        line = RDX.assemble_reaction_line(rs)
        # the LHS terms wrap their place in @conserved / @rate; the runtime parser unions these
        s = string(line)
        @test occursin("@conserved", s) && occursin("scientist", s)
        @test occursin("@rate", s) && occursin("budget", s)
    end

    @testset "E4: @select(Project,phase==:Phase2)-->@advance(phase,:Phase3) — JSON ≡ DSL behavior" begin
        # the JSON form of the Stage-C phase-advance pipeline (one Project kind, phase attribute)
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[],
          "places":[{"name":"Project","structured":true}],
          "transitions":[{"id":"adv","name":"adv","rate":1.0,"rate_mode":"deterministic",
                          "cycletime":1.0,"prob_of_success":1.0}],
          "arcs":[
            {"transition":"adv","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase2"]]}},
            {"transition":"adv","side":"rhs","advance":{"field":"phase","value":"Phase3"}} ] }
        """
        REG = Dict{Symbol, Any}(:Project => (s, f) -> RDX.SerProjectToken(get(f, :phase, :Phase2), get(f, :npv, 0.0)))
        # build the same advance_model behavior: seed two Phase2 + one Phase1, advance the Phase2's
        toks() = [RDX.SerProjectToken(:Phase2, 1.0), RDX.SerProjectToken(:Phase2, 2.0), RDX.SerProjectToken(:Phase1, 3.0)]
        p = RDX.from_json_model(json; seed = 1, registry = REG, population = toks())
        # structural match ignoring LineNumberNodes (macrocalls carry source-line metadata)
        import MacroTools
        @test MacroTools.striplines(p.network[1, :trans]) ==
            MacroTools.striplines(:((@select(Project, phase == :Phase2)) → @advance(phase, :Phase3)))
        simulate(p)
        ph = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(p, "structured")))]))
        @test ph == ["Phase1", "Phase3", "Phase3"]    # the two Phase2 projects advanced
    end

    # ── E5: action + predicate + rule (de)serialization ────────────────────────────────
    @testset "E5: every ActionStmt verb round-trips through stmt_to_dict/from_dict" begin
        stmts = ReactiveDynamics.ActionStmt[
            RDX.SetMarking(:cash, 500, :inc),
            RDX.SetParams([:synergy => 1]),
            RDX.AddToken(:ProjectToken, [:phase => QuoteNode(:Phase2), :npv => 1400.0]),
            RDX.Activate(:line),
            RDX.Deactivate(:line),
            RDX.Invoke(:rebalance, Any[0.5]),
            RDX.Log("acquired"),
            RDX.Seq(ReactiveDynamics.ActionStmt[RDX.SetMarking(:cash, 300, :inc), RDX.SetParams([:s => 1])]),
        ]
        for s in stmts
            d = RDX.stmt_to_dict(s)
            s2 = RDX.stmt_from_dict(d)
            @test typeof(s2) === typeof(s)
        end
        # RawExpr is intentionally not serializable
        @test_throws Exception RDX.stmt_to_dict(RDX.RawExpr(:(x += 1)))
    end

    @testset "E5: TokenPredicate round-trips; op outside PRED_OP_WHITELIST rejected" begin
        pred = RDX.TokenPredicate(:Project, [RDX.Clause(:phase, :(==), :(:Phase2)), RDX.Clause(:npv, :(>), 100.0)])
        d = RDX.pred_to_dict(pred)
        p2 = RDX.pred_from_dict(d)
        @test p2.kind == :Project && length(p2.clauses) == 2
        @test_throws Exception RDX.pred_from_dict(Dict("kind" => "Project", "clauses" => [["phase", "~~", "x"]]))
    end

    @testset "E5: a rules[]-bearing JSON model fires the lever at the right tick" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":6.0,"dt":1.0},
          "params":[],
          "places":[{"name":"cash","init":0},{"name":"A","init":0},{"name":"B"}],
          "transitions":[{"id":"inert","name":"inert","rate":0.0,"rate_mode":"deterministic"}],
          "arcs":[{"transition":"inert","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"inert","place":"B","side":"rhs","multiplicity":1}],
          "rules":[ { "id":"lever", "fire_mode":"once",
                      "guard": {"node":"call","op":">","args":[{"node":"timeref"},{"node":"const","value":2}]},
                      "action": {"verb":"set_marking","name":"cash","mode":"inc","value":{"node":"const","value":500}} } ] }
        """
        p = RDX.from_json_model(json; seed = 1)
        @test length(p.rules) == 1
        simulate(p)
        @test p.u[RDX.find_index(:cash, p)] == 500.0     # the once-rule injected capital after t>2
        @test count(>(0.0), diff(p.sol[!, "cash"])) == 1 # exactly once
    end

    # ── E6: modality 3-axis ⟷ Set{Symbol} bijection + observables ──────────────────────
    @testset "E6: the 5 legal modality rows round-trip; illegal combos rejected" begin
        legal = [
            (:upfront, :consumed, :block),    # row 1: {}
            (:upfront, :conserved, :block),   # row 2: {conserved}
            (:perstep, :consumed, :block),    # row 3: {rate}
            (:perstep, :conserved, :block),   # row 4: {rate, conserved}
            (:upfront, :consumed, :nonblock), # row 5: {nonblock}
        ]
        for (a, r, b) in legal
            s = RDX.to_set(a, r, b)
            fs = RDX.from_set(s)
            @test (fs.allocation, fs.return_, fs.blocking) == (a, r, b)   # round-trip identity
        end
        # illegal: nonblock + conserved (CONTRACT §1.4)
        @test_throws Exception RDX.to_set(:upfront, :conserved, :nonblock)
        # JSON modality dict → Set
        @test RDX.modality_from_dict(Dict("allocation" => "perstep", "return" => "conserved", "blocking" => "block")) == Set([:rate, :conserved])
        @test RDX.modality_from_dict(Dict("allocation" => "upfront", "return" => "consumed", "blocking" => "block")) == Set{Symbol}()
    end

    # ── E7: validate — the eval-free pre-load self-check ───────────────────────────────
    @testset "E7: a valid model validates clean; each rule's violation yields a diagnostic" begin
        import JSON
        valid = JSON.parse(
            """
            { "meta":{"tspan":5.0,"dt":1.0},
              "params":[{"name":"beta","value":0.4}],
              "places":[{"name":"A","init":100},{"name":"B"}],
              "transitions":[{"id":"t1","rate":{"node":"call","op":"*","args":[{"node":"const","value":0.3},{"node":"ref","kind":"param","name":"beta"}]},
                              "prob_of_success":0.5,"cycletime":2.0}],
              "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                     {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
            """
        )
        @test isempty(RDX.validate(valid))

        # rule 1: unknown ref name
        bad_ref = deepcopy(valid); bad_ref["transitions"][1]["rate"]["args"][2]["name"] = "nonexistent"
        @test any(d -> occursin("undeclared", d.msg), RDX.validate(bad_ref))
        # rule 1, PLACE pool: a `ref` of kind `place` resolves against the declared PLACE names,
        # not the params. Regression pin (ADR 0017): the only ref in the document above is a param,
        # so the place-pool branch of `_validate_node!` was entirely uncovered — a typo in it threw
        # `UndefVarError` on every place-referencing document while the suite stayed green.
        place_ref = deepcopy(valid)
        place_ref["transitions"][1]["rate"]["args"][2] =
            Dict("node" => "ref", "kind" => "place", "name" => "A")
        @test isempty(RDX.validate(place_ref))
        bad_place_ref = deepcopy(place_ref)
        bad_place_ref["transitions"][1]["rate"]["args"][2]["name"] = "nonexistent"
        @test any(d -> occursin("undeclared", d.msg), RDX.validate(bad_place_ref))
        # rule 1: bad op
        bad_op = deepcopy(valid); bad_op["transitions"][1]["rate"]["op"] = "system"
        @test any(d -> occursin("OP_WHITELIST", d.msg), RDX.validate(bad_op))
        # rule 2: dangling arc FK
        bad_fk = deepcopy(valid); bad_fk["arcs"][1]["transition"] = "ghost"
        @test any(d -> occursin("dangling", d.msg), RDX.validate(bad_fk))
        # rule 3: prob_of_success out of [0,1]
        bad_pos = deepcopy(valid); bad_pos["transitions"][1]["prob_of_success"] = 1.5
        @test any(d -> occursin("[0.0,1.0]", d.msg), RDX.validate(bad_pos))
        # rule 3: negative cycletime
        bad_ct = deepcopy(valid); bad_ct["transitions"][1]["cycletime"] = -1.0
        @test any(d -> occursin("≥ 0", d.msg), RDX.validate(bad_ct))
        # rule 4: illegal modality (nonblock+conserved)
        bad_mod = deepcopy(valid)
        bad_mod["places"][1]["modality"] = Dict("allocation" => "upfront", "return" => "conserved", "blocking" => "nonblock")
        @test any(d -> occursin("§1.4", d.msg), RDX.validate(bad_mod))
        # rule 1 in a predicate: Sample is not 𝓕ₜ-measurable
        bad_pred = deepcopy(valid)
        bad_pred["arcs"][1] = Dict(
            "transition" => "t1", "side" => "lhs",
            "predicate" => Dict("kind" => "A", "clauses" => [["phase", "==", Dict("node" => "sample", "dist" => "Poisson", "args" => [Dict("node" => "const", "value" => 1.0)])]])
        )
        # A is not structured AND the clause has a Sample → at least one diagnostic
        @test !isempty(RDX.validate(bad_pred))

        # from_json_model gates on validation
        @test_throws Exception RDX.from_json_model(JSON.json(bad_ref))
    end

    @testset "E7: AddToken/Invoke kind/fn must resolve against the registry" begin
        import JSON
        m = JSON.parse(
            """
            { "meta":{"tspan":5.0,"dt":1.0},"params":[],"places":[{"name":"A","init":0}],
              "transitions":[{"id":"t1","rate":1.0,"rate_mode":"deterministic"}],
              "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1}],
              "rules":[{"id":"r","fire_mode":"once",
                        "guard":{"node":"call","op":">","args":[{"node":"timeref"},{"node":"const","value":2}]},
                        "action":{"verb":"add_token","kind":"Unregistered","fields":[]}}] }
            """
        )
        @test any(d -> occursin("not in registry", d.msg), RDX.validate(m))   # empty registry
        @test isempty(
            filter(
                d -> occursin("not in registry", d.msg),
                RDX.validate(m; registry = Dict{Symbol, Any}(:Unregistered => identity))
            )
        )
        # SetField in a Rule is illegal (no bound token)
        m["rules"][1]["action"] = Dict("verb" => "set_field", "field" => "phase", "value" => Dict("node" => "const", "value" => "X"))
        @test any(d -> occursin("illegal in a Rule", d.msg), RDX.validate(m))
    end

    @testset "E7: a Field node in a @select predicate clause is rejected (ADR 0008 §D)" begin
        import JSON
        # @field is legal only in a SetField/@advance value — a Field in a predicate clause would
        # crash at runtime (@field is a macro), so validate must reject it up front.
        m = JSON.parse(
            """
            { "meta":{"tspan":5.0,"dt":1.0},"params":[],
              "places":[{"name":"Project","structured":true}],
              "transitions":[{"id":"adv","rate":1.0,"rate_mode":"deterministic"}],
              "arcs":[{"transition":"adv","side":"lhs",
                "predicate":{"kind":"Project","clauses":[["npv",">",{"node":"field","name":"npv"}]]}},
                {"transition":"adv","side":"rhs","advance":{"field":"phase","value":"Done"}}] }
            """
        )
        @test any(d -> occursin("Field", d.msg) && occursin("legal only", d.msg), RDX.validate(m))
    end

    # ── E8: the BD pipeline as model.rdj.json — JSON ≡ DSL trajectory (the north-star) ──
    @testset "E8: BD pipeline loaded from model.rdj.json matches the DSL model byte-for-byte" begin
        demodir = joinpath(pkgdir(RDX), "demo", "bd_acquisition")
        include(joinpath(demodir, "host.jl"))   # ProjectToken kind + PROJECT_REGISTRY + DSL builder
        mpath = joinpath(demodir, "model.rdj.json")
        # DSL-built
        pd = ReactionNetworkProblem(
            build_pipeline_model(); tspan = 40, dt = 1.0, seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population()
        )
        simulate(pd)
        # JSON-built (same registry, same initial population, same seed)
        pj = RDX.from_json_model(
            read(mpath, String); seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population()
        )
        simulate(pj)
        @test pd.sol == pj.sol   # identical trajectory ⇒ the JSON model is the DSL model
        phd = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(pd, "structured")))]))
        phj = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(pj, "structured")))]))
        @test phd == phj
        # the JSON model validates clean against its registry
        import JSON
        @test isempty(RDX.validate(JSON.parse(read(mpath, String)); registry = PROJECT_REGISTRY))
    end

    # ── E9: eval removal — the import-time RCE is closed ───────────────────────────────
    @testset "E9: the JSON load path never evals; a string-where-number-expected is a diagnostic" begin
        import JSON
        # the load-path source files carry no eval/Meta.parse CALLS on model data (strip comment
        # lines first — the doc-comments legitimately mention "eval"/"Meta.parse" in prose).
        code_lines(path) = join(
            filter(
                l -> !startswith(strip(l), "#"),
                split(read(joinpath(pkgdir(RDX), path), String), '\n')
            ),
            '\n',
        )
        for f in ("src/serialize.jl", "src/loadsave.jl")
            code = code_lines(f)
            @test !occursin(r"Meta\.parse\(", code)
            @test !occursin(r"\beval\(", code)
        end
        # the legacy Set{Symbol}/FoldedObservable string→eval convert hooks are gone
        rdjl = code_lines("src/ReactiveDynamics.jl")
        @test !occursin(r"convert\(::Type\{Set\{Symbol\}\}, ex::String\) = eval", rdjl)
        @test !occursin(r"convert\(::Type\{FoldedObservable\}, ex::String\) = eval", rdjl)

        # a model with a string where a number is expected is a VALIDATION error, not code execution
        bad = JSON.parse(
            """
            { "meta":{"tspan":5.0,"dt":1.0},"params":[{"name":"k","value":"run(`echo pwned`)"}],
              "places":[{"name":"A","init":0}],"transitions":[],"arcs":[] }
            """
        )
        # the malicious string is inert data — it is never parsed/eval'd (param value stays a string;
        # the engine treats it as data, no code runs). from_json builds without executing it.
        p = RDX.from_json_model(JSON.json(bad); seed = 1)
        @test p.p[:k] == "run(`echo pwned`)"   # the string is stored verbatim, never executed
    end

    @testset "E9: @import_model / @export_model round-trip a JSON model file" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[{"name":"k","value":0.5}],
          "places":[{"name":"A","init":10},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic","prob_of_success":1.0}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        tmp = tempname() * ".rdj.json"
        write(tmp, json)
        @import_model tmp prob seed = 7
        @test prob isa RDX.ReactionNetworkProblem
        @test prob.p[:k] == 0.5
        rm(tmp; force = true)
    end

    @testset "E2: model_to_dict ∘ build_acs round-trips on the parsed Dict" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{},
          "params":[{"name":"k","value":0.5}],
          "places":[{"name":"A","init":10},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic"}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        import JSON
        d = JSON.parse(json)
        net = RDX.build_network_from_dict(d)
        # params + place survive the round-trip through the acset
        back = RDX.model_to_dict(net)
        @test any(pr -> pr["name"] == "k" && pr["value"] == 0.5, back["params"])
        @test Set(pl["name"] for pl in back["places"]) == Set(["A", "B"])
    end

    # ── E10: the EXPORT path — to_json_model is the inverse of from_json_model ───────────
    # Completed _transition_to_dict / _arcs_to_dict (the inverse of assemble_reaction_line):
    # a DSL-or-JSON model → to_json_model → from_json_model is an EQUIVALENT model. The standard of
    # correctness is reconstructed-Expr equality (striplines) + a trajectory-equal simulation under
    # a fixed seed + idempotency of re-export — NOT raw-JSON-byte equality (the import-side
    # `_sum_terms` foldl re-associates an n-ary `+` reaction sum, a cosmetic Expr-nesting difference
    # that `recursive_find_arcs!` flattens identically, so the trajectory is unaffected).
    @testset "E10: a built model exports + re-imports with matching rate/attr nodes" begin
        import JSON
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":6.0,"dt":1.0},
          "params":[{"name":"beta","value":0.4}],
          "places":[{"name":"cash","init":0,"cost":2.0,"valuation":-1.0},
                    {"name":"A","init":10},{"name":"B","reward":50.0}],
          "transitions":[{"id":"t1","name":"t1",
              "rate":{"node":"call","op":"*","args":[{"node":"const","value":0.3},{"node":"ref","kind":"param","name":"beta"}]},
              "rate_mode":"poisson","prob_of_success":0.8,"cycletime":2.0,"priority":3.0}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":2},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        p = RDX.from_json_model(json; seed = 1)
        back = JSON.parse(RDX.to_json_model(p; meta = Dict("tspan" => 6.0, "dt" => 1.0)))

        # the rate node + rate_mode round-trip (Poisson-unwrapped back to the bare 0.3*beta tree)
        t1 = first(filter(t -> t["id"] == "t1", back["transitions"]))
        @test t1["rate_mode"] == "poisson"
        @test RDX.node_from_dict(t1["rate"]) ==
            RDX.Call(:*, [RDX.Const(0.3), RDX.NodeRef(:param, :beta)])
        # the ExprNode-valued attrs survive (emitted as Const nodes, the inverse of to_expr)
        @test RDX.node_from_dict(t1["prob_of_success"]) == RDX.Const(0.8)
        @test RDX.node_from_dict(t1["cycletime"]) == RDX.Const(2.0)
        @test RDX.node_from_dict(t1["priority"]) == RDX.Const(3.0)
        # non-default place attrs are emitted; defaults (e.g. cash.reward=0) are omitted
        cash = first(filter(s -> s["name"] == "cash", back["places"]))
        @test cash["cost"] == 2.0 && cash["valuation"] == -1.0 && !haskey(cash, "reward")
        @test first(filter(s -> s["name"] == "B", back["places"]))["reward"] == 50.0
        # the arcs[] decompose back to the same (place, side, multiplicity) the loader consumes
        ra = Set((r["place"], r["side"], get(r, "multiplicity", 1)) for r in back["arcs"])
        @test ra == Set([("A", "lhs", 2), ("B", "rhs", 1)])

        # full equivalence: re-import and simulate — identical trajectory under the same seed
        p2 = RDX.from_json_model(JSON.json(back); seed = 1)
        simulate(p); simulate(p2)
        @test p.sol == p2.sol
    end

    @testset "E10: modality / @select / @advance arcs are the inverse of assemble_reaction_line" begin
        import MacroTools, JSON
        # a phase-advance transition with a @select LHS, @conserved/@rate resources, integer multiplicity,
        # and an @advance RHS — exercising every arc shape _arcs_to_dict must invert.
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[],
          "places":[{"name":"Project","structured":true},{"name":"sci","init":10},{"name":"bud","init":20}],
          "transitions":[{"id":"adv","name":"adv","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":1.0}],
          "arcs":[
            {"transition":"adv","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase2"]]}},
            {"transition":"adv","side":"lhs","place":"sci","multiplicity":3,
             "modality":{"allocation":"upfront","return":"conserved","blocking":"block"}},
            {"transition":"adv","side":"lhs","place":"bud","multiplicity":5,
             "modality":{"allocation":"perstep","return":"consumed","blocking":"block"}},
            {"transition":"adv","side":"rhs","advance":{"field":"phase","value":"Phase3"}} ] }
        """
        REG = Dict{Symbol, Any}(:Project => (s, f) -> RDX.SerProjectToken(get(f, :phase, :Phase2), get(f, :npv, 0.0)))
        toks() = [RDX.SerProjectToken(:Phase2, 1.0), RDX.SerProjectToken(:Phase2, 2.0)]
        p = RDX.from_json_model(json; seed = 1, registry = REG, population = toks())
        back = JSON.parse(RDX.to_json_model(p; meta = Dict("tspan" => 5.0, "dt" => 1.0)))

        # the @select predicate is recovered exactly (kind + clause)
        sel = first(filter(r -> haskey(r, "predicate"), back["arcs"]))
        @test sel["predicate"]["kind"] == "Project"
        @test RDX.node_from_dict(sel["predicate"]["clauses"][1][3]) == RDX.Const(:Phase2)
        # the @advance RHS field-write is recovered
        adv = first(filter(r -> haskey(r, "advance"), back["arcs"]))
        @test adv["advance"]["field"] == "phase"
        @test RDX.node_from_dict(adv["advance"]["value"]) == RDX.Const(:Phase3)
        # the 3-axis modality is recovered per place
        sci = first(filter(r -> get(r, "place", "") == "sci", back["arcs"]))
        @test sci["modality"] == Dict("allocation" => "upfront", "return" => "conserved", "blocking" => "block")
        @test sci["multiplicity"] == 3
        bud = first(filter(r -> get(r, "place", "") == "bud", back["arcs"]))
        @test bud["modality"]["allocation"] == "perstep" && bud["multiplicity"] == 5

        # re-imported :trans is the SAME reaction-line Expr (striplines: macrocalls carry line meta)
        p2 = RDX.from_json_model(JSON.json(back); seed = 1, registry = REG, population = toks())
        @test MacroTools.striplines(p.network[1, :trans]) == MacroTools.striplines(p2.network[1, :trans])
        # and the two Phase2 projects advance to Phase3 identically
        simulate(p); simulate(p2)
        ph(q) = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(q, "structured")))]))
        @test ph(p) == ph(p2) == ["Phase3", "Phase3"]
    end

    @testset "E10: typed Rules round-trip through to_json_model(prob).rules[]" begin
        import JSON
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":6.0,"dt":1.0},
          "params":[],"places":[{"name":"cash","init":0},{"name":"A","init":0},{"name":"B"}],
          "transitions":[{"id":"inert","name":"inert","rate":0.0,"rate_mode":"deterministic"}],
          "arcs":[{"transition":"inert","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"inert","place":"B","side":"rhs","multiplicity":1}],
          "rules":[{"id":"lever","fire_mode":"once",
                    "guard":{"node":"call","op":">","args":[{"node":"timeref"},{"node":"const","value":2}]},
                    "action":{"verb":"set_marking","name":"cash","mode":"inc","value":{"node":"const","value":500}}}] }
        """
        p = RDX.from_json_model(json; seed = 1)
        back = JSON.parse(RDX.to_json_model(p; meta = Dict("tspan" => 6.0, "dt" => 1.0)))
        @test haskey(back, "rules") && length(back["rules"]) == 1
        @test back["rules"][1]["id"] == "lever" && back["rules"][1]["fire_mode"] == "once"
        @test back["rules"][1]["action"]["verb"] == "set_marking"
        # re-import fires the lever identically
        p2 = RDX.from_json_model(JSON.json(back); seed = 1)
        @test length(p2.rules) == 1
        simulate(p); simulate(p2)
        @test p.sol == p2.sol
        @test p2.u[RDX.find_index(:cash, p2)] == 500.0
    end

    @testset "E10: a DSL-built model (@reaction_network) exports to a re-importable JSON" begin
        import MacroTools, JSON
        # build directly with the authoring DSL (NOT from JSON), then export → re-import.
        net = @reaction_network begin
            1.0, A --> B, name => grow, probability => 0.5, cycletime => 2.0
            @deterministic(3.0), 2 * B --> C, name => merge
        end
        json = RDX.to_json_model(net; meta = Dict("tspan" => 5.0, "dt" => 1.0))
        p = RDX.from_json_model(json; seed = 1)
        # the assembled reaction lines re-parse to the same FoldedArc decomposition: a grow
        # transition A→B and a merge transition 2B→C.
        @test MacroTools.striplines(p.network[1, :trans]) == MacroTools.striplines(:(A → B))
        @test MacroTools.striplines(p.network[2, :trans]) == MacroTools.striplines(:(2B → C))
        # the rate modes are recovered: grow poisson-wrapped, merge bare (@deterministic)
        back = JSON.parse(json)
        grow = first(filter(t -> t["id"] == "grow", back["transitions"]))
        merge_ = first(filter(t -> t["id"] == "merge", back["transitions"]))
        @test grow["rate_mode"] == "poisson" && merge_["rate_mode"] == "deterministic"
        @test grow["prob_of_success"] == Dict{String, Any}("node" => "const", "value" => 0.5, "symbol" => false)
    end

    # ── E10: the BD model.rdj.json fixture — load → export → reload → trajectory-equal ──
    @testset "E10: BD model.rdj.json survives load→export→reload trajectory-equal (the fixture)" begin
        import JSON
        demodir = joinpath(pkgdir(RDX), "demo", "bd_acquisition")
        include(joinpath(demodir, "host.jl"))
        mpath = joinpath(demodir, "model.rdj.json")
        meta = Dict("tspan" => 40.0, "dt" => 1.0, "alloc_strategy" => "weighted")

        p1 = RDX.from_json_model(
            read(mpath, String); seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population()
        )
        json2 = RDX.to_json_model(p1; meta = meta)            # the now-real export
        p2 = RDX.from_json_model(
            json2; seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population()
        )
        simulate(p1); simulate(p2)
        @test p1.sol == p2.sol                                # trajectory-equal under the same seed
        ph(q) = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(q, "structured")))]))
        @test ph(p1) == ph(p2)

        # re-export validates clean against the host registry (still an eval-free, loadable model)
        @test isempty(RDX.validate(JSON.parse(json2); registry = PROJECT_REGISTRY))

        # idempotency: to_json_model(reload(to_json_model(m))) parses equal to to_json_model(m)
        p3 = RDX.from_json_model(
            json2; seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population()
        )
        @test JSON.parse(RDX.to_json_model(p3; meta = meta)) == JSON.parse(json2)
    end

    @testset "E10: @export_model writes a JSON file that @import_model reloads equivalently" begin
        import JSON
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[{"name":"k","value":0.5}],
          "places":[{"name":"A","init":10},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic","prob_of_success":1.0}],
          "arcs":[{"transition":"t1","place":"A","side":"lhs","multiplicity":1},
                 {"transition":"t1","place":"B","side":"rhs","multiplicity":1}] }
        """
        prob = RDX.from_json_model(json; seed = 7)
        tmp = tempname() * ".rdj.json"
        @export_model prob tmp                       # now writes a FULL, loadable model (was a stub)
        @import_model tmp reloaded seed = 7
        @test reloaded isa RDX.ReactionNetworkProblem
        @test reloaded.p[:k] == 0.5
        simulate(prob); simulate(reloaded)
        @test prob.sol == reloaded.sol
        rm(tmp; force = true)
    end

    # The integration seam between ADR-0012 inputs[] IMPORT (inputs_from_dict) and the completed
    # JSON EXPORT (model_to_dict): a model that declares external read ports must round-trip them.
    # The ports + pre-wire defaults live on the ReactionNetworkProblem (external_input_defaults),
    # not the net, so to_json_model(::ReactionNetworkProblem) threads them through model_to_dict's
    # inputs[] kwarg — the inverse of inputs_from_dict. Closes the round-trip-symmetry gap.
    @testset "E10: declared inputs[] ports survive load → export → reload (ADR 0012 ⟷ export)" begin
        import JSON
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[],
          "places":[{"name":"x","init":10}],
          "transitions":[{"id":"t1","rate":{"node":"const","value":1.0},"rate_mode":"deterministic"}],
          "arcs":[{"transition":"t1","place":"x","side":"lhs"}],
          "inputs":[{"port":"ext_rate","default":{"node":"const","value":0.5}},
                    {"port":"sentiment","default":{"node":"const","value":1.0}}] }
        """
        prob = RDX.from_json_model(json)
        @test prob.external_input_defaults == Dict(:ext_rate => 0.5, :sentiment => 1.0)
        d = JSON.parse(RDX.to_json_model(prob))
        @test haskey(d, "inputs")                         # inputs[] is emitted on export
        @test Set(p["port"] for p in d["inputs"]) == Set(["ext_rate", "sentiment"])
        @test isempty(RDX.validate(d))                    # the re-exported document validates clean
        prob2 = RDX.from_json_model(RDX.to_json_model(prob))
        @test prob2.external_input_defaults == prob.external_input_defaults  # ports round-trip
    end

    # ── ADR 0003 Phase 2: the promoted ArcSpec incidence table ─────────────────────────────
    @testset "ArcSpec table: population, FK exactness, escape-hatch, and JSON round-trip" begin
        net = @reaction_network begin
            1.0, 2 * A + @conserved(B) --> C, name => rx
        end
        RDX.populate_arcs!(net)
        rs = RDX.arcs(net)
        # every static arc carries an in-range integer FK and no escape-hatch expr.
        static = filter(r -> r.place != 0, rs)
        @test !isempty(static)
        @test all(r -> 1 <= r.place <= RDX.nrows(net, :S), static)
        @test all(r -> r.expr === nothing, static)
        # FK targets match the place names / sides / multiplicity the reaction line declares.
        byname = Dict(RDX.placename(net, r.place) => r for r in static)
        @test haskey(byname, :A) && byname[:A].side == :lhs && byname[:A].multiplicity == 2.0
        @test haskey(byname, :B) && byname[:B].side == :lhs && :conserved in byname[:B].modality
        @test haskey(byname, :C) && byname[:C].side == :rhs
        # JSON round-trip: the table is DERIVED from :trans, which round-trips, so re-populating the
        # reloaded model reproduces the same FK rows (place-name → side → multiplicity).
        @prob_params net
        json = RDX.to_json_model(net; meta = Dict{String, Any}("tspan" => 5.0))
        acs2 = RDX.build_network_from_dict(RDX.JSON.parse(json))
        RDX.populate_arcs!(acs2)
        rt(m) = sort(
            [
                (string(RDX.placename(m, r.place)), r.side, Float64(r.multiplicity))
                    for r in RDX.arcs(m) if r.place != 0
            ]
        )
        @test rt(acs2) == rt(net)

        # Escape-hatch: a dynamic RHS (@advance field write, ADR 0008) is a place=0 / expr-carried
        # row, NOT a static FK — the table records it without inventing a bogus FK.
        acs3 = @reaction_network begin
            1.0, @select(Project, phase == :Phase2) --> @advance(phase, :Phase3), name => adv
        end
        RDX.populate_arcs!(acs3)
        rs3 = RDX.arcs(acs3)
        @test any(r -> r.place == 0 && r.expr !== nothing, rs3)   # a dynamic term is escape-hatched
    end

    # ── ADR 0017 Tier 3: the retired wire keys still LOAD, for one release ──────────────────
    # The vocabulary rename reached the format last: `species[]` → `places[]`, `reactants[]` →
    # `arcs[]`, an arc's `species` → `place`, a `ref` kind `species` → `place`, the verb
    # `set_species` → `set_marking`. The WRITER emits only the new spellings; the READER accepts
    # the retired ones for ONE release and `Base.depwarn`s on the first one it meets (depwarn is
    # maxlog=1 per call site and a no-op unless `--depwarn=yes` — the same contract as the Tier-1
    # `@deprecate` name shims, which retire in the same release). So: a pre-rename document must
    # still produce the SAME model, and re-exporting it must migrate its keys.
    @testset "ADR 0017: a retired-key document still loads (deprecated) and re-exports renamed" begin
        import JSON
        # one model, two spellings — the only difference is the key names
        legacy = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":6.0,"dt":1.0},
          "params":[{"name":"beta","value":0.4}],
          "species":[{"name":"cash","init":0},{"name":"A","init":100},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic",
                          "cycletime":0.0,"prob_of_success":1.0}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","multiplicity":1},
                       {"transition":"t1","species":"B","side":"rhs","multiplicity":1}],
          "rules":[{"id":"lever","fire_mode":"once",
                    "guard":{"node":"call","op":">","args":[{"node":"ref","kind":"species","name":"A"},
                                                            {"node":"const","value":0}]},
                    "action":{"verb":"set_species","name":"cash","mode":"inc","value":{"node":"const","value":500}}}] }
        """
        renamed = replace(
            legacy,
            "\"species\":[" => "\"places\":[", "\"reactants\":[" => "\"arcs\":[",
            "\"species\":\"" => "\"place\":\"", "\"kind\":\"species\"" => "\"kind\":\"place\"",
            "\"set_species\"" => "\"set_marking\"",
        )
        @test !occursin("species", renamed)      # the twin really is fully renamed

        if Base.JLOptions().depwarn == 2
            # `--depwarn=error` turns every retired-key read into a throw; that is all that is
            # assertable under it, and the functional checks below cannot run.
            @test_throws ErrorException RDX.validate(JSON.parse(legacy))
        else
            # The FIRST retired key met in the session warns (maxlog=1 per call site), so this
            # assertion has to come before any other read of a legacy document.
            Base.JLOptions().depwarn == 1 && @test_deprecated RDX.validate(JSON.parse(legacy))

            @test isempty(RDX.validate(JSON.parse(legacy)))     # a legacy document validates clean
            p_old = RDX.from_json_model(legacy; seed = 11)
            p_new = RDX.from_json_model(renamed; seed = 11)
            # same places, same arcs, same rule — and therefore the same trajectory
            @test p_old.network[:, :placeName] == p_new.network[:, :placeName]
            @test p_old.network[:, :trans] == p_new.network[:, :trans]
            @test p_old.rules[1].action isa RDX.SetMarking
            simulate(p_old)
            simulate(p_new)
            @test p_old.sol == p_new.sol

            # The writer emits ONLY the new spellings, so re-exporting migrates the document.
            back = JSON.parse(RDX.to_json_model(p_old; meta = Dict("tspan" => 6.0, "dt" => 1.0)))
            @test haskey(back, "places") && !haskey(back, "species")
            @test haskey(back, "arcs") && !haskey(back, "reactants")
            @test !any(haskey(r, "species") for r in back["arcs"])
            @test Set(r["place"] for r in back["arcs"]) == Set(["A", "B"])
            @test back["rules"][1]["action"]["verb"] == "set_marking"
            @test isempty(RDX.validate(back))
        end
    end

end
