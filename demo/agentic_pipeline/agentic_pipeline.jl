# ════════════════════════════════════════════════════════════════════════════════════════
# ReactiveDynamics.jl — a literate tour of the PRODUCTION & AGENTIC capabilities
# ════════════════════════════════════════════════════════════════════════════════════════
#
# Run it:   julia --project=. demo/agentic_pipeline/agentic_pipeline.jl
#
# This script is an essay-with-code. Each section opens with a block comment that explains a
# concept and WHY the engine models it the way it does, then runs a small, self-contained
# example, then narrates the result with `println`. The numbers are deliberately tiny and
# didactic — the point is the machinery, not a heavy case study (for the heavy one, see the
# sibling `demo/bd_acquisition/`, which exercises the SAME machinery on a real pharma deal).
#
# The running theme: a small R&D PROJECT PORTFOLIO. Each project is a *structured token* with
# a lifecycle `phase` (Phase1 → Phase2 → Phase3 → Launched) and a net-present value `npv`.
# Projects are advanced through the pipeline by predicate-selected transitions; an in-model
# DECISION RULE acts as a management lever; the whole model is then serialized to an eval-free
# JSON document, validated, checkpointed, and replayed.
#
# Why ONE script for the whole tour? Because the strongest claim ReactiveDynamics makes is that
# a model — its species, its pipeline transitions, its management levers, AND its initial
# portfolio — is reproducible DATA, fully determined by `(model, population, seed)`. A single
# command that builds, runs, serializes, validates, checkpoints, and replays the very same
# model is that claim, executable.
#
# Everything below uses only constructs that appear in the engine's PASSING semantic test suite
# (test/semantic/{token_filtration,rules_decisions,initial_state,serialization_ir}.jl). No
# invented API.

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, register_structured_species!, add_structured_token!,
    Rule, Seq, SetSpecies, SetParams, SetTokens, AddToken, Activate, Deactivate, Log,
    get_species, inners, getagent, find_index, TokenPredicate, Clause, PopulationEntry,
    from_json_model, to_json_model, validate, dump_state, restore, apply_action!, set_guard!
using Random, Distributions, DataFrames
import JSON

# A short qualified alias — the structured-token TYPE and several helpers live in the
# ReactiveDynamics module (that is where the bind/advance machinery can see them).
const RDX = ReactiveDynamics

banner(title) = (println(); println("="^78); println(title); println("="^78))

# ════════════════════════════════════════════════════════════════════════════════════════
# §0. The structured-token KIND — a project is a first-class entity, not a count
# ════════════════════════════════════════════════════════════════════════════════════════
#
# A *classical* reaction-network species is a scalar: a single Float64 saying "how many A are
# there." That is perfect for indistinguishable molecules, but a project is not a molecule. We
# want each project to carry ATTRIBUTES (its current phase, its value) and a stable IDENTITY
# (the SAME object as it advances Phase1 → Phase2 → …, so a downstream report can follow it).
#
# A *structured token* gives us exactly that: a first-class agent carrying host-Julia fields,
# whose identity (uuid / kind / creation_index) is preserved as the engine mutates its fields.
# We define the kind in ReactiveDynamics' own scope via the @register/@aagent idiom, because the
# engine's selection / advancement machinery lives there and must see the type. (This exact
# shape is required; the bare-name `@structured_token` documentation form fails outside RD.)
#
# The four leading constructor arguments are the @aagent protocol fields, in order:
#   name::String            — a unique token name
#   species::Symbol         — the KIND tag (here :Project; every project shares one kind)
#   bound_transition        — nothing (the engine sets this when a transition binds the token)
#   past_bonds              — an empty Tuple{Symbol,Float64,Transition}[] history vector
# …followed by our modeling attributes: `phase` and `npv`.

@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
        phase::Symbol     # lifecycle stage — the canonical "phase-as-attribute" (ADR 0008 §D)
        npv::Float64      # the program's net-present value (a plain descriptor field)
    end
    function ProjectToken(phase, npv)
        return ProjectToken(
            "Proj" * string(rand(1:(10^9))),                       # name
            :Project,                                            # kind (one kind for all phases)
            nothing,                                             # bound_transition
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[], # past_bonds
            phase,
            npv,
        )
    end
end

# The REGISTRY maps a kind symbol to a constructor `(state, fields::Dict) -> token`. This is the
# ADR 0006 §C boundary: a serialized model (JSON, or a declarative population, or a restored
# checkpoint) references host token kinds BY NAME, and the registry is how those names resolve
# to real Julia constructors WITHOUT the data file carrying any code. The same registry serves
# the declarative population, the JSON loader, and `restore`.
const REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RDX.ProjectToken(get(f, :phase, :Phase1), get(f, :npv, 100.0)),
)

# A small helper: collect the live token pool as a Vector (the structured container is a Dict
# keyed by token name; we usually want the values).
livetokens(p) = collect(values(inners(getagent(p, "structured"))))
# Count live (non-retired) projects in a given phase.
nphase(p, ph) = count(t -> get_species(t) == :Project && t.phase == ph, livetokens(p))
# How many soft-retired (failed a stage gate; species flipped to :removed)?
nretired(p) = count(t -> get_species(t) == :removed, livetokens(p))

banner("§0. The structured-token KIND")
println("Defined kind :Project as a ProjectToken{phase::Symbol, npv::Float64}.")
println(
    "One demo token: ", let t = RDX.ProjectToken(:Phase1, 120.0)
        "phase=$(t.phase), npv=$(t.npv), kind=$(get_species(t))"
    end
)
println("Why a structured token and not a Float64 count? A count is anonymous and stateless;")
println("a token carries attributes (phase, npv) AND a stable identity preserved across phase")
println("advances — so we can select projects by attribute and follow each one to launch.")

# ════════════════════════════════════════════════════════════════════════════════════════
# §1. The phase-as-attribute pipeline + the declarative initial portfolio
# ════════════════════════════════════════════════════════════════════════════════════════
#
# A naive model would make a SPECIES per phase (Phase1, Phase2, …) and "advance" by destroying
# a Phase1 token and creating a Phase2 token. That breaks identity (the new token is a different
# object) and multiplies the species count. The canonical ReactiveDynamics design (ADR 0008 §D)
# is PHASE-AS-ATTRIBUTE: there is ONE :Project kind, and `phase` is a field. A pipeline step is
#
#     @select(Project, <clause>) --> @advance(phase, :NextPhase)
#
# `@select(Project, clauses)` binds only tokens of kind Project whose attributes satisfy the
# conjunctive `&&` clause (ops: == != < <= > >= in). `@advance(phase, :Phase2)` then writes the
# bound token's `phase` field IN PLACE — same object, identity preserved.
#
# A stage gate can also FAIL: `probability => q` makes each advance a Binomial(·, q) trial. On
# failure the bound token SOFT-RETIRES — its species flips to `:removed`, and its `phase` field
# records how far it got (a killed Phase2 program stays at phase==:Phase2 but species==:removed).
#
# We seed the starting portfolio with the DECLARATIVE `population[]` initial marking (ADR 0007
# §B), passed to the constructor and instantiated before t=0. This is preferred over an
# imperative `add_structured_token!` loop because the run is then reproducible from
# `(model, population, seed)` — the portfolio is reproducible INPUT, not post-construction host
# code. Two authoring forms exist; we use both here.

function pipeline_model()
    net = @reaction_network begin
        # Phase1 -> Phase2 : a sure step (probability 1) for didactic clarity.
        @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 1.0
        # Phase2 -> Phase3 : a RISKY gate — only 60% of attempts succeed; the rest soft-retire.
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => adv23, cycletime => 1.0, probability => 0.6
        # Phase3 -> Launched : a fairly safe late stage.
        @deterministic(1.0),
            @select(Project, phase == :Phase3) --> @advance(phase, :Launched),
            name => adv3L, cycletime => 1.0, probability => 0.9
    end
    register_structured_species!(net, :Project)
    return net
end

# Form A — an explicit host-token list (the maintainer's "instantiate as a list of structures"):
# a fixed, hand-authored opening portfolio.
explicit_portfolio() = [
    RDX.ProjectToken(:Phase1, 120.0),
    RDX.ProjectToken(:Phase1, 90.0),
    RDX.ProjectToken(:Phase2, 200.0),
    RDX.ProjectToken(:Phase2, 150.0),
    RDX.ProjectToken(:Phase3, 300.0),
]

banner("§1. Phase-as-attribute pipeline + declarative initial portfolio")
p1 = ReactionNetworkProblem(
    pipeline_model(); tspan = 6, dt = 1.0, seed = 1,
    registry = REGISTRY, population = explicit_portfolio()
)
println("Seeded an explicit 5-project portfolio (population[] Form A — a list of structs).")
println(
    "t=0 by phase:  Phase1=$(nphase(p1, :Phase1))  Phase2=$(nphase(p1, :Phase2))  ",
    "Phase3=$(nphase(p1, :Phase3))  Launched=$(nphase(p1, :Launched))"
)
simulate(p1)
println("After 6 ticks (Phase2->Phase3 gate is only 60% — some programs fail and soft-retire):")
println("  Launched = $(nphase(p1, :Launched))   still-in-flight Phase3 = $(nphase(p1, :Phase3))")
println("  soft-retired (:removed) = $(nretired(p1))   (their phase records how far they got)")
println("Identity note: a project that reached Launched is the SAME object that started in")
println("Phase1 — @advance writes the field in place, it does not create a new token.")

# Form B — the PopulationEntry "count + attribute exprs" form: "N programs with these attrs",
# registry-built and seeded. Symbol literals use QuoteNode; a value can be a SAMPLED expression
# drawn from the run's seeded RNG (`state.rng`), so a same-seed construction is reproducible.
function sampled_portfolio()
    return [
        # 4 Phase1 programs with NPV sampled ~ Normal(100, 15), reproducible under the seed:
        PopulationEntry(
            :Project, :Project; count = 4,
            attributes = Dict(
                :phase => QuoteNode(:Phase1),
                :npv => :(rand(state.rng, Normal(100.0, 15.0)))
            )
        ),
        # 3 Phase2 programs at a fixed NPV:
        PopulationEntry(
            :Project, :Project; count = 3,
            attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 200.0)
        ),
    ]
end

p1b = ReactionNetworkProblem(
    pipeline_model(); tspan = 4, dt = 1.0, seed = 42,
    registry = REGISTRY, population = sampled_portfolio()
)
println()
println("population[] Form B (PopulationEntry: count + seeded attribute exprs):")
println("  built ", length(livetokens(p1b)), " projects = 4 Phase1 (NPV sampled) + 3 Phase2.")
println(
    "  sampled Phase1 NPVs (seeded, reproducible): ",
    round.(sort([t.npv for t in livetokens(p1b) if t.phase == :Phase1]); digits = 1)
)

# ════════════════════════════════════════════════════════════════════════════════════════
# §2. Predicate selection — advancing only the qualifying subset
# ════════════════════════════════════════════════════════════════════════════════════════
#
# The power of `@select` is that the LHS predicate is a 𝓕ₜ-measurable FILTER over token
# attributes: only the matching subset is bindable. A continuous clause like `npv > θ` lets a
# transition act on a value threshold — e.g. "only fast-track high-value Phase2 programs."
#
# When several tokens match but the transition can only fire on a few per tick, WHICH bind first
# is deterministic: equal-priority ties break by (species, creation_index) — the earlier-added
# token wins — NOT by the agent dictionary's hash order or the tokens' random names. So a
# predicate-selected pipeline reproduces exactly under `(model, seed)`.

function fasttrack_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(phase, :Phase3),
            name => fasttrack, cycletime => 1.0, probability => 1.0
    end
    register_structured_species!(net, :Project)
    return net
end

banner("§2. Predicate selection (a continuous npv > θ clause)")
p2 = ReactionNetworkProblem(
    fasttrack_model(); tspan = 3, dt = 1.0, seed = 1,
    registry = REGISTRY,
    population = [
        RDX.ProjectToken(:Phase2, 100.0),   # below θ=150 — stays in Phase2
        RDX.ProjectToken(:Phase2, 220.0),   # above θ — fast-tracked to Phase3
        RDX.ProjectToken(:Phase2, 180.0),   # above θ — fast-tracked to Phase3
        RDX.ProjectToken(:Phase1, 999.0),   # wrong phase — never selected (clause is conjunctive)
    ]
)
println("Predicate: @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(:Phase3)")
println("Before:  Phase1=$(nphase(p2, :Phase1))  Phase2=$(nphase(p2, :Phase2))  Phase3=$(nphase(p2, :Phase3))")
simulate(p2)
println("After:   Phase1=$(nphase(p2, :Phase1))  Phase2=$(nphase(p2, :Phase2))  Phase3=$(nphase(p2, :Phase3))")
println("Only the two high-NPV Phase2 programs advanced; the npv=100 program stayed (below θ),")
println("and the Phase1 program was never eligible (the && clause gates BOTH phase and npv).")
println("Bind order is deterministic: equal-priority ties break by creation_index (first added")
println("wins), so this selection reproduces exactly under the same (model, seed).")

# ════════════════════════════════════════════════════════════════════════════════════════
# §3. The in-model DECISION RULE — a management lever that lives in the model
# ════════════════════════════════════════════════════════════════════════════════════════
#
# This is the agentic heart of the engine (ADR 0010). A management decision — "if conditions
# hold, take an action" — is itself part of the model, as a typed `Rule`, NOT host patch code
# that reaches in and mutates the run from outside. A Rule has:
#
#     Rule(id, guard::Expr, action; fire_mode = :once | :every_tick)
#
# The guard is evaluated against the live state (`@t()` is the clock; species/params are in
# scope). `fire_mode = :once` fires the action exactly once, the first tick its guard holds, then
# latches OFF (`p.rules[i].enabled == false`); `_reinit!` re-arms it. The ACTION family — all
# verified — composes via `Seq`:
#   SetSpecies(:cash, 500, :inc)            — inject into a resource pool (:inc or :set)
#   SetParams([:synergy => 1])              — flip a model parameter the transitions read
#   AddToken(:Project, [...])               — inject a brand-new token via the registry (BY KIND;
#                                             the registry key, here :Project — see REGISTRY in §0)
#   Activate(:line) / Deactivate(:line)     — soft-gate a transition on/off
#   Seq([...]) / Log("msg")                 — compose / annotate
#
# Our lever models a "Series-B raise + portfolio expansion": once t > 2, inject capital into a
# `cash` pool, flip a `synergy` parameter, and add a fresh Phase2 project — all in one Seq, all
# IN the model. We also gate a `fund` transition on `cash >= 50` via `set_guard!`, so that line
# only comes alive once the raise lands.
#
# Because `cash` must be a real species column, we build a small model with a `cash` pool and a
# `fund` line that converts cash into a `report` (a stand-in for "spend the raise"), plus the
# Phase2->Phase3 pipeline step so AddToken has somewhere to land.

function lever_model()
    net = @reaction_network begin
        # a spendable line, gated below on cash >= 50 (genesis withheld until funded)
        @deterministic(1.0), cash --> report, name => fund
        # the pipeline step the injected project will flow through
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => adv23, cycletime => 1.0, probability => 1.0
    end
    @prob_init net cash = 0 report = 0
    @prob_params net synergy = 0
    register_structured_species!(net, :Project)
    return net
end

# The lever: a once-rule firing at t > 2 that raises capital, flips synergy, and adds a project.
raise_lever() = Rule(
    :series_b, :(@t() > 2.0),
    Seq(
        [
            SetSpecies(:cash, 500, :inc),                                    # +500 capital
            SetParams([:synergy => 1]),                                      # flip the synergy flag
            AddToken(:Project, [:phase => QuoteNode(:Phase2), :npv => 175.0]), # add a Phase2 program
            Log("Series-B raised: +500 cash, synergy on, +1 Phase2 program"),
        ]
    );
    fire_mode = :once
)

banner("§3. The in-model decision rule (a management lever, ADR 0010)")
p3 = ReactionNetworkProblem(
    lever_model(); tspan = 6, dt = 1.0, seed = 1,
    registry = REGISTRY,
    population = [RDX.ProjectToken(:Phase2, 200.0)],   # one organic Phase2 program at t=0
    rules = [raise_lever()]
)
set_guard!(p3, :fund, :(cash >= 50))   # the `fund` line only fires once the raise lands
ci_cash = find_index(:cash, p3)
println("t=0:  cash = $(p3.u[ci_cash])   synergy = $(p3.p[:synergy])   Phase2 projects = $(nphase(p3, :Phase2))")
println("Rule: once @t() > 2, Seq[ +500 cash, synergy:=1, AddToken(Phase2 npv=175) ]")
println("Guard: set_guard!(:fund, cash >= 50) — the spend line is withheld until funded.")
simulate(p3)
println("After the run:")
cash_series = p3.sol[!, "cash"]
report_series = p3.sol[!, "report"]
println("  cash injected once   : peak cash reached $(maximum(cash_series)) (0 -> 500 in a single jump)")
println("  positive cash jumps  : $(count(>(0.0), diff(cash_series))) (a :once rule fires exactly once)")
println("  synergy param flipped: $(p3.p[:synergy])")
println("  rule latched off     : enabled = $(p3.rules[1].enabled) (re-armed by _reinit!)")
println("  `report` produced    : $(report_series[end]) (>0 ⇒ the funded line did fire once cash≥50)")
println("The decision lives in the MODEL — the driver only set the seed and armed the rule; no")
println("host code reached in mid-run to mutate the state. That is what makes the scenario a")
println("reproducible (model, rules, seed) triple rather than an imperative script.")

# ════════════════════════════════════════════════════════════════════════════════════════
# §4. Genesis as a transition PRODUCT — the agentic constructor on the RHS
# ════════════════════════════════════════════════════════════════════════════════════════
#
# §3 added a token via an `AddToken` action buried inside a Rule's `Seq` — an imperative
# side-effect on the decision channel. There is a second, more primitive way to BIRTH a token,
# and it is first-class in the dynamics: a transition whose RHS PRODUCT is a token, written
#
#   ∅ --> @structured(:Project, phase = :Phase1, npv = rand(state.rng, …), born = @t())
#
# This is structurally parallel to a plain source reaction `∅ --> budget` (which mints a scalar) —
# except the product is a full agentic token, entangled live into the structured pool. Genesis is
# thus a transition PRODUCT, not a rule action.
#
# `@structured` is the RHS-product TWIN of §3's `AddToken`: the reaction line carries the registry
# KIND (`:Project`) plus field-value expressions, and the host constructor is resolved BY NAME
# through the registry at firing time — the reaction line NEVER carries the constructor. It shares
# AddToken's exact `(state, fields::Dict) -> token` contract, so BOTH genesis paths use one
# registry. Because it carries only a name + typed field nodes, it ROUND-TRIPS through the
# eval-free JSON IR (we prove that below, and again in §6) — genesis-as-product and
# genesis-as-rule-action are the same operation in two positions. This is the ONLY @structured
# form: an inline-constructor form `@structured(Ctor(…))` — the only reactant construct that could
# not serialize eval-free — was removed, so eval-free serialization is now a TOTAL invariant (every
# genesis product is data). A raw constructor on the RHS is rejected at construction (see the note
# at the end of this section).
#
# The field values are evaluated at firing time INSIDE the run's context, so they can read live
# state — here `@t()` (the clock, stamped into `born`) and `state.rng` (the seeded RNG, so a drawn
# `npv` is reproducible under the seed). We give the kind a `born` field to SEE the clock captured
# at each birth.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct GenesisProjectToken
        phase::Symbol
        npv::Float64
        born::Float64      # the clock value @t() captured at construction — proves genesis is live
    end
    function GenesisProjectToken(phase, npv, born)
        return GenesisProjectToken(
            "Gen" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            npv,
            born,
        )
    end
end

# The registry the NAMED @structured form resolves `:Project` through — the SAME `(state, fields)`
# convention §3's AddToken uses, so the genesis product and the rule action share one host contract.
const GENESIS_REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RDX.GenesisProjectToken(
        get(f, :phase, :Phase1), get(f, :npv, 0.0), get(f, :born, state.t)
    ),
)

# The pipeline: a genesis source that BIRTHS one Phase1 project per tick (npv drawn from the seeded
# RNG, birth-time stamped from the clock) via the named, registry-resolved form, and a downstream
# @select/@advance leg the newborns flow through. Birth → select → advance, all in the dynamics.
function genesis_model()
    net = @reaction_network begin
        @deterministic(1.0),
            ∅ --> @structured(
                :Project, phase = :Phase1,
                npv = rand(state.rng, Normal(120.0, 20.0)), born = @t()
            ),
            name => genesis
        @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 1.0
    end
    register_structured_species!(net, :Project)
    return net
end

banner("§4. Genesis as a transition product (@structured RHS — the agentic constructor)")
pg = ReactionNetworkProblem(genesis_model(); tspan = 5, dt = 1.0, seed = 1, registry = GENESIS_REGISTRY)
println("t=0 live tokens: ", length(livetokens(pg)), " (the source has not fired yet).")
simulate(pg)
gtoks = livetokens(pg)
println("RHS: ∅ --> @structured(:Project, phase=:Phase1, npv=rand(state.rng,·), born=@t())")
println("After the run the source minted ", length(gtoks), " tokens, each a distinct agent:")
println("  born times (from @t() at construction): ", sort([t.born for t in gtoks]))
println("  npv values (drawn from the seeded RNG): ", round.(sort([t.npv for t in gtoks]); digits = 1))
println(
    "  every name unique (independent identity): ",
    length(unique(RDX.getname.(gtoks))) == length(gtoks)
)
println(
    "  phases now (newborns flowed through @select/@advance to Phase2): ",
    "Phase1=$(nphase(pg, :Phase1)) Phase2=$(nphase(pg, :Phase2))"
)
# same seed ⇒ identical births (the field exprs draw from the run's seeded rng)
pg2 = ReactionNetworkProblem(genesis_model(); tspan = 5, dt = 1.0, seed = 1, registry = GENESIS_REGISTRY)
simulate(pg2)
println(
    "Reproducible: same-seed npvs identical? ",
    sort([t.npv for t in gtoks]) ≈ sort([t.npv for t in livetokens(pg2)])
)

# Genesis is DATA: because @structured carries only the kind name + typed field nodes (not the
# constructor), the genesis model exports to the eval-free JSON IR and reloads loss-free — the
# same round-trip §6 makes for the whole pipeline, here exercised on a @structured RHS. Because the
# raw inline-constructor form was removed, THIS ALWAYS HOLDS: any model the engine accepts exports.
gjson = to_json_model(pg)
pg_rt = from_json_model(gjson; seed = 1, registry = GENESIS_REGISTRY); simulate(pg_rt)
println(
    "Genesis round-trips: exported model validates clean? ",
    isempty(validate(JSON.parse(gjson); registry = GENESIS_REGISTRY)),
    "; reload reproduces births? ",
    sort([t.born for t in gtoks]) == sort([t.born for t in livetokens(pg_rt)])
)
println("Contrast §3: there a token was ADDED by a rule ACTION (AddToken, decision channel); here")
println("it is BORN as a transition PRODUCT (@structured), the agentic analogue of ∅ --> species —")
println("sharing AddToken's registry, so it serializes as eval-free data too.")

# A raw inline constructor on the RHS is REJECTED at construction — @structured is named-only, so
# eval-free serialization is a total invariant. Show the rejection (caught, for the demo).
raw_line = :(@structured(GenesisProjectToken(:Phase1, 100.0, 0.0)))
try
    @eval @reaction_network begin
        @deterministic(1.0), ∅ --> $raw_line, name => genesis
    end
    println("Raw @structured(Ctor(…)): UNEXPECTEDLY accepted")
catch e
    msg = sprint(showerror, e)
    println(
        "Raw @structured(Ctor(…)) rejected at construction: ",
        occursin("named form", msg) ? "✓ (points to the named form)" : msg
    )
end

# ════════════════════════════════════════════════════════════════════════════════════════
# §5. Population queries & writes — SetTokens with @field
# ════════════════════════════════════════════════════════════════════════════════════════
#
# Sometimes a decision must rewrite an ATTRIBUTE across a selected sub-population — e.g. "on a
# bad competitor readout, write down every Phase-2 valuation by 10%." `SetTokens` is the action
# that does this: it takes a TokenPredicate (the same selection logic as `@select`) and a list
# of field-update expressions. The right-hand value may reference `@field(name)`, which reads
# the selected token's OWN current field value.
#
# Two important rules, both enforced by the engine and its validator:
#   • `@field(npv)` is legal ONLY in a SetField/@advance/SetTokens VALUE — it is a syntactic
#     marker substituted to a literal before evaluation. It is NOT legal in a @select predicate
#     clause (a Field there would crash, since @field is a macro). The validator rejects that.
#   • a predicate clause matching a SYMBOL literal writes it as `:(:Phase2)` (a QuoteNode-bearing
#     expr); a numeric clause uses the bare number, e.g. `Clause(:npv, :(>), 150.0)`.

banner("§5. Population write — SetTokens(@field) writes down a selected sub-population")
p4 = ReactionNetworkProblem(
    pipeline_model(); tspan = 3, dt = 1.0, seed = 1,
    registry = REGISTRY,
    population = [
        RDX.ProjectToken(:Phase2, 100.0),
        RDX.ProjectToken(:Phase2, 200.0),
        RDX.ProjectToken(:Phase1, 50.0),   # not Phase2 — must be left untouched
    ]
)
writedown = SetTokens(
    TokenPredicate(:Project, [Clause(:phase, :(==), :(:Phase2))]),
    [:npv => :(@field(npv) * 0.9)],   # each selected token's own npv, ×0.9
)
println("Before:  ", sort([(string(t.phase), t.npv) for t in livetokens(p4)]))
apply_action!(p4, nothing, writedown)   # apply the population write directly
println("After :  ", sort([(string(t.phase), t.npv) for t in livetokens(p4)]))
println("Every Phase2 valuation written down 10% (100->90, 200->180); the Phase1 program (50)")
println("was not selected and is untouched. @field read each token's OWN npv before scaling it.")

# ════════════════════════════════════════════════════════════════════════════════════════
# §6. Model-as-DATA — the eval-free JSON model (ADR 0005)
# ════════════════════════════════════════════════════════════════════════════════════════
#
# A model in ReactiveDynamics is not just Julia code — it is DATA. The same pipeline can be
# written as a JSON document and loaded with `from_json_model`. The loader NEVER `eval`s or
# `Meta.parse`s a model field: it walks a typed ExprNode tree and lowers it through the unchanged
# compiler. This is the security / agentic-authoring story: an untrusted party (or an LLM) can
# AUTHOR a model as JSON, and the worst a malicious field can do is fail validation — it cannot
# execute code at load time. Host token kinds are referenced BY NAME and resolved through the
# registry; the JSON carries no Julia.
#
# We (a) write the SAME pipeline as JSON, (b) `validate` it (clean), (c) show a deliberately
# broken model producing a diagnostic, (d) load it and confirm the JSON-built run matches the
# DSL-built run under the same seed, (e) make the security point with an inert malicious string,
# (f) load the document from a file with `@import_model`, and (g) the INVERSE — export a LIVE
# model back to a full JSON document with `to_json_model`/`@export_model` and reload it loss-free
# (a model is DATA in BOTH directions: author-as-JSON → load, and build/load → export → reload).

# The pipeline as a JSON model. Structured species carry "structured": true; a pipeline step's
# reactants are an LHS predicate + an RHS advance. (Bare string "Phase2" in clause arrays.)
const PIPELINE_JSON = """
{ "rd_format":"reactive-dynamics-model", "version":"1.0",
  "meta":{ "tspan":6.0, "dt":1.0 },
  "params":[],
  "species":[ {"name":"Project","structured":true} ],
  "transitions":[
    {"id":"adv12","name":"adv12","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":1.0},
    {"id":"adv23","name":"adv23","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":0.6},
    {"id":"adv3L","name":"adv3L","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":0.9} ],
  "reactants":[
    {"transition":"adv12","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase1"]]}},
    {"transition":"adv12","side":"rhs","advance":{"field":"phase","value":"Phase2"}},
    {"transition":"adv23","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase2"]]}},
    {"transition":"adv23","side":"rhs","advance":{"field":"phase","value":"Phase3"}},
    {"transition":"adv3L","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase3"]]}},
    {"transition":"adv3L","side":"rhs","advance":{"field":"phase","value":"Launched"}} ] }
"""

banner("§6. Model-as-data — the eval-free JSON model (ADR 0005)")

# (b) validate the JSON model — `validate` returns a Vector of diagnostics; [] means clean.
diags_ok = validate(JSON.parse(PIPELINE_JSON); registry = REGISTRY)
println("(b) validate(clean model) -> ", isempty(diags_ok) ? "OK (no diagnostics)" : diags_ok)

# (c) a deliberately broken model: a dangling reactant foreign-key (transition that doesn't exist).
broken = JSON.parse(PIPELINE_JSON)
broken["reactants"][1]["transition"] = "ghost"   # no transition with id "ghost"
diags_bad = validate(broken; registry = REGISTRY)
println("(c) validate(broken: dangling FK) -> ", length(diags_bad), " diagnostic(s):")
for d in diags_bad
    println("      ", string(d))
end

# (d) JSON-built run ≡ DSL-built run under the same seed and population.
shared_pop() = [
    RDX.ProjectToken(:Phase1, 120.0),
    RDX.ProjectToken(:Phase2, 200.0),
    RDX.ProjectToken(:Phase2, 150.0),
    RDX.ProjectToken(:Phase3, 300.0),
]
p_dsl = ReactionNetworkProblem(
    pipeline_model(); tspan = 6, dt = 1.0, seed = 7,
    registry = REGISTRY, population = shared_pop()
)
p_json = from_json_model(PIPELINE_JSON; seed = 7, registry = REGISTRY, population = shared_pop())
simulate(p_dsl)
simulate(p_json)
phases_of(p) = sort(string.([t.phase for t in livetokens(p)]))
println("(d) DSL final phases : ", phases_of(p_dsl))
println("    JSON final phases: ", phases_of(p_json))
println(
    "    trajectories identical (sol == sol): ", p_dsl.sol == p_json.sol,
    "  ⇒ the JSON model IS the DSL model under the same seed."
)

# (e) the security point: a malicious string param value is stored as INERT DATA, never executed.
malicious = """
{ "rd_format":"reactive-dynamics-model", "version":"1.0", "meta":{"tspan":3.0,"dt":1.0},
  "params":[ {"name":"k","value":"run(`echo pwned`)"} ],
  "species":[ {"name":"A","init":0} ],
  "transitions":[], "reactants":[] }
"""
p_mal = from_json_model(malicious; seed = 1)
println("(e) security: a param whose value is the string \"run(`echo pwned`)\" loads as inert data:")
println("      p.p[:k] = ", repr(p_mal.p[:k]), "  (stored verbatim — from_json_model never eval'd it)")

# (f) the file path: a JSON model document lives ON DISK and is loaded with `@import_model`. This
# is the canonical model-as-data direction (an LLM or a colleague AUTHORS the document; the engine
# LOADS it) — exactly how demo/bd_acquisition/model.rdj.json is consumed. `@import_model path name
# kw=val…` reads the file and calls `from_json_model`; we re-supply the registry + population (host
# Julia, referenced by name) to make it runnable, then simulate.
tmp = tempname() * ".rdj.json"
write(tmp, PIPELINE_JSON)                                      # the authored document on disk
@import_model tmp p_fromfile seed = 7 registry = REGISTRY population = shared_pop()
simulate(p_fromfile)
println(
    "(f) @import_model from a file: rebuilt a ", typeof(p_fromfile).name.name,
    "; final phases match the in-memory build? ", phases_of(p_fromfile) == phases_of(p_dsl)
)
rm(tmp; force = true)

# (g) the INVERSE direction — emit a LIVE model back to a full JSON document. `to_json_model`
# (and the `@export_model` macro) is the structural inverse of `from_json_model`: it walks the
# constructed model's stored columns — the rate (Poisson-unwrapped to its bare intensity + a
# rate_mode), the ExprNode-valued attrs, and the reaction line decomposed back into reactants[]
# (the inverse of the import-time reaction-line assembly) — and re-emits the eval-free document.
# A model is therefore DATA in both directions: author-as-JSON → load, AND build/load → export.
# We export the live DSL-built model, reload the emitted JSON, and confirm the reload is the SAME
# model — identical trajectory under the same seed (and re-validating clean).
exported = to_json_model(p_dsl)                            # live model → eval-free JSON document
p_roundtrip = from_json_model(exported; seed = 7, registry = REGISTRY, population = shared_pop())
simulate(p_roundtrip)
println("(g) to_json_model(live model) -> reload -> simulate:")
println(
    "    re-exported model validates clean?  ",
    isempty(validate(JSON.parse(exported); registry = REGISTRY))
)
println("    reloaded final phases match DSL?     ", phases_of(p_roundtrip) == phases_of(p_dsl))
println(
    "    trajectories identical (sol == sol): ", p_dsl.sol == p_roundtrip.sol,
    "  ⇒ build/load → export → reload is loss-free; a model is DATA in both directions."
)
# Idempotency: exporting the reload reproduces the same document (round-trip is a fixed point).
println(
    "    export idempotent (re-export == export)? ",
    JSON.parse(to_json_model(p_roundtrip)) == JSON.parse(exported)
)

# (h) the run's OUTPUT (the solution trajectory) is a SEPARATE artifact from the model document —
# the model is the reproducible input, the trajectory is its result (ADR 0005 §76).
soltable = @export_solution_as_table p_fromfile           # the trajectory as a DataFrame
println(
    "(h) solution trajectory exported as a ", size(soltable, 1), "×", size(soltable, 2),
    " DataFrame via @export_solution_as_table (the run OUTPUT, distinct from the model)."
)

# ════════════════════════════════════════════════════════════════════════════════════════
# §7. Checkpoint & replay — dump_state / restore + reinit determinism
# ════════════════════════════════════════════════════════════════════════════════════════
#
# `dump_state(p)` serializes a LIVE run at a clean tick boundary into an eval-free, JSON-able
# artifact: the clock, the RNG state, creation counters, the plain `u`, the full token population
# with current field values, and the once-rule latches. `restore(spec, dump; registry)` rebuilds
# an identical problem from it. Resuming the restored problem reproduces the original's
# continuation exactly.
#
# A Milestone-1 requirement: the dump only supports a CLEAN tick boundary (no in-flight
# instance mid-cycle). We design this section's model with `cycletime => 0.0`, so every advance
# completes within its tick and `ongoing_transitions` is empty at the boundary.
#
# Finally we show full reproducibility: `_reinit!` resets the state, rebuilds the t=0 population,
# and re-arms once-rules, so re-running from the same seed reproduces the first trajectory.

function instant_pipeline()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 0.0, probability => 1.0   # ct=0 ⇒ clean boundary every tick
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => adv23, cycletime => 0.0, probability => 1.0
    end
    register_structured_species!(net, :Project)
    return net
end

banner("§7. Checkpoint & replay (dump_state / restore + reinit determinism)")
spec = instant_pipeline()
cp_pop() = [
    PopulationEntry(
        :Project, :Project; count = 4,
        attributes = Dict(:phase => QuoteNode(:Phase1), :npv => 150.0)
    ),
]
pc = ReactionNetworkProblem(
    spec; tspan = 10, dt = 1.0, seed = 1,
    registry = REGISTRY, population = cp_pop()
)
simulate(pc, 2)                          # step to a clean boundary (ct=0 ⇒ no in-flight)
println(
    "Simulated 2 ticks. ongoing_transitions empty? ", isempty(pc.ongoing_transitions),
    "   t = ", pc.t
)
d = dump_state(pc)
println("dump_state: t=$(d.t), $(length(d.tokens)) tokens captured (eval-free, JSON-able).")
pc2 = restore(spec, d; registry = REGISTRY)
println("restore  : t=$(pc2.t), tokens=$(length(livetokens(pc2))), u matches? ", pc2.u == pc.u)
simulate(pc)                             # resume the original
simulate(pc2)                            # resume the restored copy
println(
    "Resume both to t=10 — continuation matches? ", phases_of(pc) == phases_of(pc2),
    "   (", phases_of(pc), ")"
)

# reinit determinism: same model + seed ⇒ identical trajectory after a reset.
pr = ReactionNetworkProblem(
    pipeline_model(); tspan = 6, dt = 1.0, seed = 3,
    registry = REGISTRY, population = explicit_portfolio()
)
simulate(pr); sol1 = copy(pr.sol); ph1 = phases_of(pr)
AlgebraicAgents._reinit!(pr)             # reset state + rebuild t=0 population + re-arm once-rules
simulate(pr)
println(
    "reinit replay: trajectory reproduced (sol == sol)? ", pr.sol == sol1,
    "   final phases reproduced? ", phases_of(pr) == ph1
)

# ════════════════════════════════════════════════════════════════════════════════════════
# §8. Recap — production & agentic capabilities, and the ADRs they realize
# ════════════════════════════════════════════════════════════════════════════════════════
banner("§8. Recap")
println(
    """
    We toured, on ONE small R&D portfolio, the engine's production & agentic machinery:

      §0  Structured tokens          first-class project entities (attributes + identity)   ADR 0006/0008
      §1  Phase-as-attribute + pop[]  one :Project kind, phase is a field; declarative input  ADR 0008/0007
      §2  Predicate selection         @select(npv > θ) binds a subset; deterministic ties     ADR 0008
      §3  In-model decision rule      a once-Rule lever: SetSpecies+SetParams+AddToken in Seq  ADR 0010
      §4  Genesis as a product        ∅ --> @structured(:Kind, …): a token BORN on the RHS      ADR 0006/0008
      §5  Population write            SetTokens(@field) revalues a selected sub-population      ADR 0011
      §6  Model-as-data (JSON)        eval-free load + export round-trip; JSON ≡ DSL ≡ reload    ADR 0005
      §7  Checkpoint & replay         dump_state/restore at a clean boundary; reinit determinism ADR 0007

    The through-line: a model — species, pipeline, levers, and starting portfolio — is reproducible
    DATA, fully determined by (model, population, rules, seed). A token can enter the run three ways —
    declaratively at t=0 (population[], §1), imperatively via a rule action (AddToken, §3), or as a
    first-class transition product (@structured, §4) — the decision logic lives IN the model as typed
    Rules, the model serializes to an eval-free document an untrusted author can safely produce, and
    any point of a run can be checkpointed and replayed bit-for-bit.
    """
)
println("Done.")
