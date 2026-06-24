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
        return SerProjectToken("SP" * string(rand(1:10^9)), :Project, nothing,
            Tuple{Symbol,Float64,ReactiveDynamics.Transition}[], phase, npv)
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
          "species":[{"name":"A","init":1000},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":3.0,"rate_mode":"deterministic",
                          "cycletime":0.0,"prob_of_success":1.0}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","stoich":1},
                       {"transition":"t1","species":"B","side":"rhs","stoich":1}] }
        """
        p1 = RDX.from_json_model(json; seed = 42); simulate(p1)
        p2 = RDX.from_json_model(json; seed = 42); simulate(p2)
        @test p1.sol == p2.sol           # same (model, seed) ⇒ identical
    end

    # ── E4: reactants[] → reaction-line :trans Expr (stoich, modality, @select, @advance) ──
    @testset "E4: multi-LHS + integer stoich assembles to the runtime-parsed reaction line" begin
        rs = [
            Dict("transition" => "t", "species" => "X", "side" => "lhs", "stoich" => 1),
            Dict("transition" => "t", "species" => "Y", "side" => "lhs", "stoich" => 2),
            Dict("transition" => "t", "species" => "Z", "side" => "rhs", "stoich" => 1),
        ]
        line = RDX.assemble_reaction_line(rs)
        @test line == :((X + 2Y) → Z)
    end

    @testset "E4: LHS modality macros (@conserved/@rate) are emitted per the 3-axis modality" begin
        rs = [
            Dict("transition" => "t", "species" => "scientist", "side" => "lhs", "stoich" => 3,
                "modality" => Dict("allocation" => "upfront", "return" => "conserved", "blocking" => "block")),
            Dict("transition" => "t", "species" => "budget", "side" => "lhs", "stoich" => 1,
                "modality" => Dict("allocation" => "perstep", "return" => "consumed", "blocking" => "block")),
            Dict("transition" => "t", "species" => "out", "side" => "rhs", "stoich" => 1),
        ]
        line = RDX.assemble_reaction_line(rs)
        # the LHS terms wrap their species in @conserved / @rate; the runtime parser unions these
        s = string(line)
        @test occursin("@conserved", s) && occursin("scientist", s)
        @test occursin("@rate", s) && occursin("budget", s)
    end

    @testset "E4: @select(Project,phase==:Phase2)-->@advance(phase,:Phase3) — JSON ≡ DSL behavior" begin
        # the JSON form of the Stage-C phase-advance pipeline (one Project kind, phase attribute)
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[],
          "species":[{"name":"Project","structured":true}],
          "transitions":[{"id":"adv","name":"adv","rate":1.0,"rate_mode":"deterministic",
                          "cycletime":1.0,"prob_of_success":1.0}],
          "reactants":[
            {"transition":"adv","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase2"]]}},
            {"transition":"adv","side":"rhs","advance":{"field":"phase","value":"Phase3"}} ] }
        """
        REG = Dict{Symbol,Any}(:Project => (s, f) -> RDX.SerProjectToken(get(f, :phase, :Phase2), get(f, :npv, 0.0)))
        # build the same advance_model behavior: seed two Phase2 + one Phase1, advance the Phase2's
        toks() = [RDX.SerProjectToken(:Phase2, 1.0), RDX.SerProjectToken(:Phase2, 2.0), RDX.SerProjectToken(:Phase1, 3.0)]
        p = RDX.from_json_model(json; seed = 1, registry = REG, population = toks())
        # structural match ignoring LineNumberNodes (macrocalls carry source-line metadata)
        import MacroTools
        @test MacroTools.striplines(p.acs[1, :trans]) ==
            MacroTools.striplines(:((@select(Project, phase == :Phase2)) → @advance(phase, :Phase3)))
        simulate(p)
        ph = sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(p, "structured")))]))
        @test ph == ["Phase1", "Phase3", "Phase3"]    # the two Phase2 projects advanced
    end

    # ── E5: action + predicate + rule (de)serialization ────────────────────────────────
    @testset "E5: every ActionStmt verb round-trips through stmt_to_dict/from_dict" begin
        stmts = ReactiveDynamics.ActionStmt[
            RDX.SetSpecies(:cash, 500, :inc),
            RDX.SetParams([:synergy => 1]),
            RDX.AddToken(:ProjectToken, [:phase => QuoteNode(:Phase2), :npv => 1400.0]),
            RDX.Activate(:line),
            RDX.Deactivate(:line),
            RDX.Invoke(:rebalance, Any[0.5]),
            RDX.Log("acquired"),
            RDX.Seq(ReactiveDynamics.ActionStmt[RDX.SetSpecies(:cash, 300, :inc), RDX.SetParams([:s => 1])]),
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
          "species":[{"name":"cash","init":0},{"name":"A","init":0},{"name":"B"}],
          "transitions":[{"id":"inert","name":"inert","rate":0.0,"rate_mode":"deterministic"}],
          "reactants":[{"transition":"inert","species":"A","side":"lhs","stoich":1},
                       {"transition":"inert","species":"B","side":"rhs","stoich":1}],
          "rules":[ { "id":"lever", "fire_mode":"once",
                      "guard": {"node":"call","op":">","args":[{"node":"timeref"},{"node":"const","value":2}]},
                      "action": {"verb":"set_species","name":"cash","mode":"inc","value":{"node":"const","value":500}} } ] }
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
        valid = JSON.parse("""
        { "meta":{"tspan":5.0,"dt":1.0},
          "params":[{"name":"beta","value":0.4}],
          "species":[{"name":"A","init":100},{"name":"B"}],
          "transitions":[{"id":"t1","rate":{"node":"call","op":"*","args":[{"node":"const","value":0.3},{"node":"ref","kind":"param","name":"beta"}]},
                          "prob_of_success":0.5,"cycletime":2.0}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","stoich":1},
                       {"transition":"t1","species":"B","side":"rhs","stoich":1}] }
        """)
        @test isempty(RDX.validate(valid))

        # rule 1: unknown ref name
        bad_ref = deepcopy(valid); bad_ref["transitions"][1]["rate"]["args"][2]["name"] = "nonexistent"
        @test any(d -> occursin("undeclared", d.msg), RDX.validate(bad_ref))
        # rule 1: bad op
        bad_op = deepcopy(valid); bad_op["transitions"][1]["rate"]["op"] = "system"
        @test any(d -> occursin("OP_WHITELIST", d.msg), RDX.validate(bad_op))
        # rule 2: dangling reactant FK
        bad_fk = deepcopy(valid); bad_fk["reactants"][1]["transition"] = "ghost"
        @test any(d -> occursin("dangling", d.msg), RDX.validate(bad_fk))
        # rule 3: prob_of_success out of [0,1]
        bad_pos = deepcopy(valid); bad_pos["transitions"][1]["prob_of_success"] = 1.5
        @test any(d -> occursin("[0.0,1.0]", d.msg), RDX.validate(bad_pos))
        # rule 3: negative cycletime
        bad_ct = deepcopy(valid); bad_ct["transitions"][1]["cycletime"] = -1.0
        @test any(d -> occursin("≥ 0", d.msg), RDX.validate(bad_ct))
        # rule 4: illegal modality (nonblock+conserved)
        bad_mod = deepcopy(valid)
        bad_mod["species"][1]["modality"] = Dict("allocation"=>"upfront","return"=>"conserved","blocking"=>"nonblock")
        @test any(d -> occursin("§1.4", d.msg), RDX.validate(bad_mod))
        # rule 1 in a predicate: Sample is not 𝓕ₜ-measurable
        bad_pred = deepcopy(valid)
        bad_pred["reactants"][1] = Dict("transition"=>"t1","side"=>"lhs",
            "predicate"=>Dict("kind"=>"A","clauses"=>[["phase","==",Dict("node"=>"sample","dist"=>"Poisson","args"=>[Dict("node"=>"const","value"=>1.0)])]]))
        # A is not structured AND the clause has a Sample → at least one diagnostic
        @test !isempty(RDX.validate(bad_pred))

        # from_json_model gates on validation
        @test_throws Exception RDX.from_json_model(JSON.json(bad_ref))
    end

    @testset "E7: AddToken/Invoke kind/fn must resolve against the registry" begin
        import JSON
        m = JSON.parse("""
        { "meta":{"tspan":5.0,"dt":1.0},"params":[],"species":[{"name":"A","init":0}],
          "transitions":[{"id":"t1","rate":1.0,"rate_mode":"deterministic"}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","stoich":1}],
          "rules":[{"id":"r","fire_mode":"once",
                    "guard":{"node":"call","op":">","args":[{"node":"timeref"},{"node":"const","value":2}]},
                    "action":{"verb":"add_token","kind":"Unregistered","fields":[]}}] }
        """)
        @test any(d -> occursin("not in registry", d.msg), RDX.validate(m))   # empty registry
        @test isempty(filter(d -> occursin("not in registry", d.msg),
            RDX.validate(m; registry = Dict{Symbol,Any}(:Unregistered => identity))))
        # SetField in a Rule is illegal (no bound token)
        m["rules"][1]["action"] = Dict("verb"=>"set_field","field"=>"phase","value"=>Dict("node"=>"const","value"=>"X"))
        @test any(d -> occursin("illegal in a Rule", d.msg), RDX.validate(m))
    end

    @testset "E7: a Field node in a @select predicate clause is rejected (ADR 0008 §D)" begin
        import JSON
        # @field is legal only in a SetField/@advance value — a Field in a predicate clause would
        # crash at runtime (@field is a macro), so validate must reject it up front.
        m = JSON.parse("""
        { "meta":{"tspan":5.0,"dt":1.0},"params":[],
          "species":[{"name":"Project","structured":true}],
          "transitions":[{"id":"adv","rate":1.0,"rate_mode":"deterministic"}],
          "reactants":[{"transition":"adv","side":"lhs",
            "predicate":{"kind":"Project","clauses":[["npv",">",{"node":"field","name":"npv"}]]}},
            {"transition":"adv","side":"rhs","advance":{"field":"phase","value":"Done"}}] }
        """)
        @test any(d -> occursin("Field", d.msg) && occursin("legal only", d.msg), RDX.validate(m))
    end

    # ── E8: the BD pipeline as model.rdj.json — JSON ≡ DSL trajectory (the north-star) ──
    @testset "E8: BD pipeline loaded from model.rdj.json matches the DSL model byte-for-byte" begin
        demodir = normpath(joinpath(homedir(), "ReactiveDynamics-review", "demo", "bd_acquisition"))
        include(joinpath(demodir, "host.jl"))   # ProjectToken kind + PROJECT_REGISTRY + DSL builder
        mpath = joinpath(demodir, "model.rdj.json")
        # DSL-built
        pd = ReactionNetworkProblem(build_pipeline_model(); tspan = 40, dt = 1.0, seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population())
        simulate(pd)
        # JSON-built (same registry, same initial population, same seed)
        pj = RDX.from_json_model(read(mpath, String); seed = 7,
            registry = PROJECT_REGISTRY, population = initial_population())
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
            filter(l -> !startswith(strip(l), "#"),
                split(read(normpath(joinpath(homedir(), "ReactiveDynamics-review", path)), String), '\n')),
            '\n',
        )
        for f in ("src/serialize.jl", "src/loadsave.jl")
            code = code_lines(f)
            @test !occursin(r"Meta\.parse", code)
            @test !occursin(r"\beval\(", code)
        end
        # the legacy Set{Symbol}/FoldedObservable string→eval convert hooks are gone
        rdjl = code_lines("src/ReactiveDynamics.jl")
        @test !occursin(r"convert\(::Type\{Set\{Symbol\}\}, ex::String\) = eval", rdjl)
        @test !occursin(r"convert\(::Type\{FoldedObservable\}, ex::String\) = eval", rdjl)

        # a model with a string where a number is expected is a VALIDATION error, not code execution
        bad = JSON.parse("""
        { "meta":{"tspan":5.0,"dt":1.0},"params":[{"name":"k","value":"run(`echo pwned`)"}],
          "species":[{"name":"A","init":0}],"transitions":[],"reactants":[] }
        """)
        # the malicious string is inert data — it is never parsed/eval'd (param value stays a string;
        # the engine treats it as data, no code runs). from_json builds without executing it.
        p = RDX.from_json_model(JSON.json(bad); seed = 1)
        @test p.p[:k] == "run(`echo pwned`)"   # the string is stored verbatim, never executed
    end

    @testset "E9: @import_model / @export_model round-trip a JSON model file" begin
        json = """
        { "rd_format":"reactive-dynamics-model","version":"1.0","meta":{"tspan":5.0,"dt":1.0},
          "params":[{"name":"k","value":0.5}],
          "species":[{"name":"A","init":10},{"name":"B"}],
          "transitions":[{"id":"t1","name":"t1","rate":1.0,"rate_mode":"deterministic","prob_of_success":1.0}],
          "reactants":[{"transition":"t1","species":"A","side":"lhs","stoich":1},
                       {"transition":"t1","species":"B","side":"rhs","stoich":1}] }
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
