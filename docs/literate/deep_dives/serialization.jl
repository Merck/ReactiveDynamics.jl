# # Deep dive: a model is data — the eval-free JSON document
#
# **What this deep-dive is about.** In ReactiveDynamics a model is not opaque Julia code that happens to run — it is a single JSON *document*, and every claim the engine makes about reproducibility and safe exchange follows from that. Here we take a structured-token pipeline, export it to JSON, reload it, and show the reload reproduces the *same trajectory* under the same seed; we run the eval-free `validate` pass on a clean model and on a deliberately broken one; we show that a hostile string riding in a model field loads as inert data, never executed; and we open up the two mechanisms that make all of this possible — the closed, typed `ExprNode` IR and the host-function registry. We close on `@structured`, the eval-free RHS genesis twin of a rule's `AddToken`.
#
# **This is a topic-driven deep-dive, not a cumulative tutorial.** Each section is self-contained and seed-pinned. For the end-to-end structured-token modeling workflow these constructs come from, see the [advanced tutorial](../tutorials/advanced.md); for the exact document shape, the [JSON model schema](../reference/json_schema.md) reference; for the full API, the [serialization reference](../reference/serialization.md); and for the *why* behind the trust boundary, [ADR 0005](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0005-serialization-json-ir.md) and contract [§8](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/CONTRACT_DRAFT.md).

using ReactiveDynamics
using ReactiveDynamics: validate, inners, getagent
import JSON
using Distributions            # Normal — evaluated at firing time in a genesis field expression
using DataFrames               # the solution trajectory is a DataFrame

# A short qualified alias — the structured-token TYPE and the engine's selection/advancement
# machinery live in the ReactiveDynamics module, so a token kind must be defined and referenced there.
const RD = ReactiveDynamics

# ## The structured-token kind and its registry
#
# A model that carries projects needs a structured-token *kind*. We define it in the engine's own scope with the `@register`/`@aagent` idiom (the selection/advancement machinery must see the type), giving each `:Project` a lifecycle `phase` and a net-present value `npv`. The four leading constructor arguments are the `@aagent` protocol fields — name, kind tag, `bound_firing`, `past_bonds` — followed by our modeling attributes.

@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
        phase::Symbol
        npv::Float64
    end
    function ProjectToken(phase, npv)
        return ProjectToken(
            "Proj" * string(rand(1:(10^9))),                       # name
            :Project,                                              # kind tag
            nothing,                                               # bound_firing
            Tuple{Symbol, Float64, ReactiveDynamics.Firing}[], # past_bonds
            phase,
            npv,
        )
    end
end

# The `REGISTRY` maps a kind symbol to a constructor `(state, fields::Dict) -> token`. This is the boundary the whole design turns on: a serialized model references host token kinds BY NAME, and the registry is how those names resolve to real Julia constructors WITHOUT the document carrying any code. The same registry serves the JSON loader, the declarative population, and a restored checkpoint.

const REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RD.ProjectToken(get(f, :phase, :Phase1), get(f, :npv, 100.0)),
)

# A couple of small helpers to read the live token pool by phase.
livetokens(p) = collect(values(inners(getagent(p, "structured"))))
phases_of(p) = sort(string.([t.phase for t in livetokens(p)]))

# The model itself: a phase-as-attribute pipeline where a project is advanced in place by predicate-selected transitions, one risky gate (`adv23` at 60%) among them.
function pipeline_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 1.0
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => adv23, cycletime => 1.0, probability => 0.6
        @deterministic(1.0),
            @select(Project, phase == :Phase3) --> @advance(phase, :Launched),
            name => adv3L, cycletime => 1.0, probability => 0.9
    end
    register_token_kind!(net, :Project)
    return net
end

# A fixed opening portfolio, shared by every run below so the comparisons hold everything but the *authoring form* constant.
shared_pop() = [
    RD.ProjectToken(:Phase1, 120.0),
    RD.ProjectToken(:Phase2, 200.0),
    RD.ProjectToken(:Phase2, 150.0),
    RD.ProjectToken(:Phase3, 300.0),
]

# ## 1. The round-trip: `to_json_model` → `from_json_model`
#
# `to_json_model(prob)` walks a constructed model's stored columns — the rate (Poisson-unwrapped to its bare intensity plus a `rate_mode`), the `ExprNode`-valued attributes, and each reaction line decomposed back into the flat `arcs[]` list — and emits a single JSON string. `from_json_model(json; seed, registry, population)` is its inverse: parse → `validate` → build → construct. Because a run is fully determined by `(model, population, seed)`, the reload must reproduce the original run *exactly*.
#
# We build the model in the DSL, simulate it under a fixed seed, export it, reload the emitted JSON under the *same* seed and population, and compare the two solution trajectories.

p_dsl = ReactionNetworkProblem(
    pipeline_model(); tspan = 6, dt = 1.0, seed = 7,
    registry = REGISTRY, population = shared_pop(),
)
simulate(p_dsl)

exported = to_json_model(p_dsl)                       # live model → eval-free JSON document
p_reload = from_json_model(exported; seed = 7, registry = REGISTRY, population = shared_pop())
simulate(p_reload)

println("DSL final phases   : ", phases_of(p_dsl))
println("reload final phases: ", phases_of(p_reload))
sol_identical = p_dsl.sol == p_reload.sol
println("trajectories byte-identical (sol == sol): ", sol_identical)

# The reload is not merely *equivalent* — its trajectory is identical to the DSL-built run's, so `to_json_model ∘ from_json_model` is a faithful reconstruction, not an approximation. The round-trip is also a fixed point: exporting the reload reproduces the same document.

idempotent = JSON.parse(to_json_model(p_reload)) == JSON.parse(exported)
println("export idempotent (re-export == export): ", idempotent)

# ## 2. `validate`: a clean model returns no diagnostics; a broken one returns a `Diagnostic`
#
# `validate(dict; registry)` is a PURE walk over the parsed document — no `node_from_dict`, no `to_expr`, no `eval`. It returns a `Vector{Diagnostic}`; an empty vector means clean, and `from_json_model` gates construction on exactly that. This is the pass an agent (or a colleague) runs to self-check a model *before* loading it.
#
# We reuse the pipeline as a hand-authored JSON document — structured places carry `"structured": true`, and a pipeline step's arcs are an LHS `predicate` plus an RHS `advance`.

const PIPELINE_JSON = """
{ "rd_format":"reactive-dynamics-model", "version":"1.0",
  "meta":{ "tspan":6.0, "dt":1.0 },
  "params":[],
  "places":[ {"name":"Project","structured":true} ],
  "transitions":[
    {"id":"adv12","name":"adv12","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":1.0},
    {"id":"adv23","name":"adv23","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":0.6},
    {"id":"adv3L","name":"adv3L","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":0.9} ],
  "arcs":[
    {"transition":"adv12","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase1"]]}},
    {"transition":"adv12","side":"rhs","advance":{"field":"phase","value":"Phase2"}},
    {"transition":"adv23","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase2"]]}},
    {"transition":"adv23","side":"rhs","advance":{"field":"phase","value":"Phase3"}},
    {"transition":"adv3L","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase3"]]}},
    {"transition":"adv3L","side":"rhs","advance":{"field":"phase","value":"Launched"}} ] }
"""

# A clean document validates to the empty diagnostic vector:

clean_diags = validate(JSON.parse(PIPELINE_JSON); registry = REGISTRY)
println("validate(clean model)  -> ", isempty(clean_diags) ? "OK (no diagnostics)" : clean_diags)

# Now break it deliberately: point a arc's foreign key at a transition that does not exist. `validate` reports it as a diagnostic — it does not throw, and it certainly does not evaluate anything.

broken = JSON.parse(PIPELINE_JSON)
broken["arcs"][1]["transition"] = "ghost"  # no transition with id "ghost"
broken_diags = validate(broken; registry = REGISTRY)
println("validate(broken model) -> ", length(broken_diags), " diagnostic(s):")
for d in broken_diags
    println("    ", string(d))
end

# The broken document is *inert*: the dangling foreign key surfaces as a `Diagnostic` a caller can read and act on, and `from_json_model` would refuse to construct from it. Nothing in the document ran.

# ## 3. The registry and the RCE boundary — a hostile string loads as inert data
#
# The loader NEVER `Meta.parse`s or `eval`s a model field. Every value is either a JSON scalar taken verbatim or a node-tagged tree drawn from a closed whitelist; every host construct a model needs — a token kind, a value helper — is referenced by NAME through the registry the host program supplies. So the worst a malicious field can do is fail validation. We prove it directly: a parameter whose value is the *string* `"run(`echo pwned`)"` loads, and stays a string.

const MALICIOUS_JSON = """
{ "rd_format":"reactive-dynamics-model", "version":"1.0", "meta":{"tspan":3.0,"dt":1.0},
  "params":[ {"name":"k","value":"run(`echo pwned`)"} ],
  "places":[ {"name":"A","init":0} ],
  "transitions":[], "arcs":[] }
"""

p_mal = from_json_model(MALICIOUS_JSON; seed = 1)
println("param k loaded as: ", repr(p_mal.p[:k]))
println("still a String (never evaluated): ", p_mal.p[:k] isa String)

# The string is stored verbatim as the parameter's value — `from_json_model` treated it as data, not as a command. There is no code path from a model field to `eval`, which is what lets an untrusted party (or a language model) *author* a model as JSON and hand it over: the document is not a program.

# ## 4. The typed `ExprNode` IR — the closed node algebra
#
# The mechanism under the round-trip is the `ExprNode` IR: every rate/attribute expression is one of a closed set of typed nodes — `Const`, `NodeRef`, `Call` (op ∈ `OP_WHITELIST`), `Sample` (dist ∈ `DIST_WHITELIST`), `TimeRef`, `Choose`, `Field`, `ExternalRef`. `to_expr` lowers a node tree to exactly the Julia `Expr` the authoring macros already produce; `from_expr` recovers the tree from a stored `Expr`. Neither ever parses a source string.
#
# Take the canonical example: the rate `Poisson(0.3 * beta) * Preclinical`, built as a tree by hand.

rate_tree = Call(
    :*, [
        Sample(:Poisson, [Call(:*, [Const(0.3), NodeRef(:param, :beta)])]),
        NodeRef(:place, :Preclinical),
    ]
)
lowered = to_expr(rate_tree)
println("to_expr(rate_tree) = ", lowered)

# The tree serializes to a JSON dict — a nested `"node"`-tagged object, no Julia anywhere — and parses back losslessly. We check the reconstruction by lowering both to their `Expr` (structural `Expr` equality, since two separately-built node trees are not `===`).

node_dict = node_to_dict(rate_tree)
recovered = node_from_dict(node_dict)
println("node dict round-trips: ", to_expr(recovered) == to_expr(rate_tree))

# `from_expr` is the inverse direction — a stored attribute `Expr` classified back into the typed tree, with the place/param name sets telling a bare symbol which kind of `NodeRef` it is:

back = from_expr(lowered; places = Set([:Preclinical]), params = Set([:beta]))
println("from_expr ∘ to_expr is identity (on Expr): ", to_expr(back) == to_expr(rate_tree))

# The algebra is *closed*: an operator outside `OP_WHITELIST` cannot enter the IR. `from_expr` rejects it rather than admitting arbitrary calls — the whitelist is the gate, not a convention.

println("OP_WHITELIST = ", OP_WHITELIST)
rejected = try
    from_expr(:(sin(beta)); params = Set([:beta]))
    "UNEXPECTEDLY accepted"
catch e
    "rejected: " * sprint(showerror, e)
end
println("from_expr(:(sin(beta))): ", rejected)

# ## 5. `@structured` — the eval-free RHS genesis twin of `AddToken`
#
# A token can be BORN as a transition product, written `∅ --> @structured(:Kind, field = …, …)`. This is the RHS-product twin of a rule's `AddToken`: the reaction line carries the registry KIND plus field-value *expressions*, and the host constructor is resolved by NAME through the same `(state, fields::Dict) -> token` registry at firing time. The line never carries the constructor — so, like everything else, it round-trips through the eval-free JSON IR. (The named form is the only `@structured` form the engine accepts; a raw inline constructor is rejected at construction, which is what makes eval-free serialization a total invariant.)
#
# The field expressions are evaluated inside the run's context, so a drawn `npv` uses the run's seeded RNG and is reproducible under the seed.

function genesis_model()
    net = @reaction_network begin
        @deterministic(1.0),
            ∅ --> @structured(:Project, phase = :Phase1, npv = rand(state.rng, Normal(120.0, 20.0))),
            name => genesis
        @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 1.0
    end
    register_token_kind!(net, :Project)
    return net
end

pg = ReactionNetworkProblem(genesis_model(); tspan = 5, dt = 1.0, seed = 1, registry = REGISTRY)
simulate(pg)
gnpvs = round.(sort([t.npv for t in livetokens(pg)]); digits = 1)
println("genesis minted ", length(livetokens(pg)), " tokens; npvs (seeded, reproducible): ", gnpvs)

# Because `@structured` carries only the kind name plus typed field nodes, the genesis model exports and reloads loss-free — the same round-trip as §1, here exercised on an RHS genesis product. The reload reproduces the same births under the same seed.

gjson = to_json_model(pg)
pg_reload = from_json_model(gjson; seed = 1, registry = REGISTRY)
simulate(pg_reload)
genesis_reproduced = sort([t.npv for t in livetokens(pg)]) ≈ sort([t.npv for t in livetokens(pg_reload)])
println("genesis model validates clean: ", isempty(validate(JSON.parse(gjson); registry = REGISTRY)))
println("genesis reload reproduces births: ", genesis_reproduced)

# ## The demonstrated invariant
#
# Everything above rests on one computed fact, worth stating plainly as the thing a reader can rely on:

println()
println("The DSL-built model and the JSON-reloaded model produce byte-identical solutions")
println("under the same seed (sol == sol): ", sol_identical)
println("The deliberately broken model produced ", length(broken_diags), " validation diagnostic(s),")
println("returned as data — never an eval.")

# A model in ReactiveDynamics is therefore data in both directions: a document can be authored as JSON and loaded, or built in the DSL and exported, and either path lands on the *same* run under the same seed. The typed `ExprNode` IR plus the by-name registry is the trust boundary that makes loading inert — the property that lets models be exchanged, validated, and agent-authored without turning a data file into a program. For the exact document shape see the [JSON model schema](../reference/json_schema.md); for the full round-trip API, the [serialization reference](../reference/serialization.md); and for the rationale, [ADR 0005](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0005-serialization-json-ir.md) and contract [§8](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/CONTRACT_DRAFT.md).
