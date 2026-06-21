# ADR 0005: Single JSON serialization + typed ExprNode IR (agentic artifact)

Status: Accepted — pending Phase-0/Phase-1 implementation

Date: 2026-06-18

Relates to: ADR 0001 (keep the native discrete-event engine `ReactionNetworkProblem`), ADR 0002 (priority-weighted water-filling allocation), ADR 0003 (drop ACSets for a dependency-free typed-struct-of-columns IR; this ADR is the serialization counterpart of that IR and depends on its Phase 2 promotion of the reactant relation). The maintainer has chosen a single JSON format, dropping the TOML/CSV/JLD2 zoo, with NO backward on-disk compatibility. This ADR decides only the serialization/authoring boundary; the runtime hot path (compiled closures, allocator, clock) is unaffected.

## Context

ReactiveDynamics is a timed, stochastic, resource-constrained GSPN / discrete-event system used to model business/R&D processes (pipelines, budgeting, ledgers, rNPV, what-if). The authoring layer must now serve two roles at once: an on-disk model format, and an agentic-authoring artifact that an LLM can emit under structured output and self-validate before it is loaded. The current serialization layer cannot serve either role safely.

Ground-truth facts, verified against source on branch ref-agents:

- The current loader is a multi-format zoo in `src/loadsave.jl`: TOML/CSV import builds the ACSet by reflection over `propertynames(acs.subparts)` with string-prefix routing, and solutions are persisted as JLD2 binary blobs (`loadsave.jl:205,221`) plus an optional `CSV.write` (`loadsave.jl:254`). There is no schema, no validation, and three incompatible file formats to maintain.
- The loader executes arbitrary code at import time. Parameter values are `eval`d if they arrive as strings — `attrval = attrval isa String ? eval(Meta.parseall(attrval)) : attrval` (`loadsave.jl:65`) — and an entire `registered` source block is `eval`d wholesale — `eval(Meta.parseall(row["body"]))` (`loadsave.jl:72`). This is remote-code-execution on load: a model file is a Julia program.
- The same eval-on-assignment behavior is wired into the type system via `Base.convert` hooks: `Base.convert(::Type{SampleableValues}, ex::String) = MacroTools.striplines(Meta.parse(ex))` (`ReactiveDynamics.jl:99`), `Base.convert(::Type{Set{Symbol}}, ex::String) = eval(Meta.parse(ex))` (`ReactiveDynamics.jl:101`), and the same for `FoldedObservable` (`ReactiveDynamics.jl:102`). So even attribute assignment parses/evals strings.
- Expression-valued attributes (rate, stoich, priority, cost, cycletime, prob-of-success, observable ranges, event trigger/action) are today Julia `Expr` in a catch-all union: `const SampleableValues = Union{Expr,Symbol,AbstractString,Float64,Int,Function}` (`ReactiveDynamics.jl:10`). An `Expr` or a `Function` is not JSON-representable, so a naive "dump the Expr as a string" approach would re-introduce the very `Meta.parse`/`eval` round-trip we are trying to remove.
- Those Exprs are compiled to closures by positionally substituting species/param names: `compile_attrs` builds a construction-time `varmap` of `name => :(state.u[$i])` and `name => :(state.p[$(QuoteNode(name))])` (`compilers.jl:148-153`), then maps every non-skip attribute through `wrap_fun` (`compilers.jl:160-180`). The substitution is the hard part to preserve — but it operates on an `Expr`, so any new front-end only needs to produce that same `Expr`.
- The transition↔reactant relation is NOT stored as structure today; it is re-parsed from the `trans` Expr at runtime by `extract_reactants` (`reaction_parser.jl:32`), called per-tick for the LHS and for RHS emission. ADR 0003 Phase 2 promotes this to a first-class typed `ReactantSpec` incidence table (`docs/adr/0003-data-store.md:37`). Serialization must mirror that promotion.
- The dynamic escape-hatches are `@structured`/`@move` (recognized as macrocalls at `reaction_parser.jl:67` and passed through at `create.jl:300`) and `@choose` (a separate branch, recognized at `reaction_parser.jl:13` inside `recursively_choose` and at `create.jl:292`). These are the only legitimately dynamic reactant idioms and must become explicit typed nodes rather than arbitrary code.
- Genesis is Poisson per tick: `expand_rate` wraps a rate as `:(rand(Poisson(max(state.dt * $rate, 0))))` unless the `@deterministic` macro is present, in which case the bare rate is used (`create.jl:150-155`). The verified runtime distribution draws are: the Poisson spawn draw (`create.jl:151`), the Binomial probability-of-success draw (`solvers.jl:413` — the engine's ONLY Binomial), and the event-firing Poisson draw (`solvers.jl:321`).
- The solution is an output, not part of the model: `sol = DataFrame("t" => Float64[], (specName => Float64[])...)` (`solvers.jl:588-591`), and a `log` ledger accumulates `(:allocation, …)` and `(:valuation_cost, …)` rows (`solvers.jl:304-312`) and a `(:valuation, …)` row each tick (`solvers.jl:659-666`). The ledger is exactly what the rNPV/BD demo consumes.
- `meta` should absorb both the `:M` keyword rows and the solver kwargs because both already merge into one `keywords` bag at construction (`solvers.jl:544-552`).
- JSON tooling is already in the stack: `JSON` (the older JSON.jl) is a declared dependency (`Project.toml:16,48`) and `Distributions` is too (`Project.toml:12,53`). NOTE: the proposed mechanism below uses `JSON3`/`StructTypes` (and `Arrow` for solutions), which are NET-NEW dependencies absent from `Project.toml`; the "low-cost" framing is scoped to "JSON is already an accepted format in this stack," not "zero new deps." The validation rules are pre-specified in CONTRACT_DRAFT.md §1.4 (illegal modality combos) and §5 (ranges, integrality, time-varying-expression policy).

## Decision

Adopt ONE canonical serialization format, JSON (file extension `model.rdj.json`), serving as both the on-disk model file and the agentic-authoring artifact, and remove the TOML/CSV/JLD2 model zoo and every eval-on-load path. The format is EVAL-FREE by construction: no field is ever `Meta.parse`d or `eval`d. This replaces `loadsave.jl:65,72`, the `Base.convert` eval hooks at `ReactiveDynamics.jl:99,101,102`, and the JLD2 solution path at `loadsave.jl:205,221`.

### Document shape (mirrors the ADR-0003 typed IR one-to-one)

A model is one JSON object with a `meta` object and the top-level arrays `params[]`, `species[]`, `transitions[]`, `reactants[]`, `observables[]`, `events[]`. The mapping to the verified schema objects:

- `meta` absorbs the `:M` keyword rows and the solver kwargs (`tspan`, `dt` or `tstops`, `tunit`, `seed`, `alloc_strategy`), all of which already merge into the single `keywords` bag at `solvers.jl:544-552`. This eliminates the `eval`d-string `metaVal` path.
- `params[]` maps to `:P` (`prmName`/`prmVal`). Values are JSON numbers, not strings — directly replacing the `eval(Meta.parseall(prmVal))` at `loadsave.jl:65`.
- `species[]` maps to `:S`. `species.modality` uses the orthogonal 3-axis form from CONTRACT §1.1 (`allocation∈{upfront,perstep}`, `return∈{consumed,conserved}`, `blocking∈{block,nonblock}`), not the unvalidated `Set{Symbol}`. It may be omitted (defaults to `{upfront,consumed,block}`). A translation layer maps the 3-axis form back to the legacy `Set{Symbol}` the engine still consumes until the §1 re-model lands.
- `transitions[]` maps to `:T`, each with `id` (the FK target for reactants), `name`, a `rate` ExprNode, a `rate_mode∈{poisson,deterministic}` flag (the `@deterministic` discriminator at `create.jl:150-155`), and ExprNode-valued `cycletime`/`prob_of_success`/`capacity`/`priority`/`max_lifetime`, plus optional `pre_action`/`post_action` statement trees.
- `reactants[]` is the PROMOTED first-class ReactantSpec incidence table (ADR 0003 Phase 2): each row carries `side∈{lhs,rhs}`, FK `transition` (id), FK `species` (name), an ExprNode `stoich`, and (lhs only) a `modality`. This is the structural promotion of the relation that is today re-parsed from the `trans` Expr at runtime (`extract_reactants`, `reaction_parser.jl:32`). The `@choose`/`@structured`/`@move` escape hatches (`reaction_parser.jl:67`, `create.jl:300`) become explicit ExprNode variants so this table stays the single source of truth.
- `observables[]` maps to `:obs` (`FoldedObservable`), `events[]` to `:E` (`eventTrigger`/`eventAction`). The BD acquisition lever — a mid-simulation injection — is an event whose `action` is a `SetSpecies`/`SetParams` statement node.

### ExprNode IR (closed, JSON-representable sum type)

Every "TVE? = yes" attribute in CONTRACT §5 is encoded NOT as a Julia source string but as a closed tagged-union tree, keyed on a `node` tag:

```julia
abstract type ExprNode end
struct Const   <: ExprNode; value::Union{Float64,Int,Bool} end
struct Ref     <: ExprNode; kind::Symbol; name::Symbol end           # kind ∈ REF_KINDS
struct Call    <: ExprNode; op::Symbol;   args::Vector{ExprNode} end # op ∈ OP_WHITELIST
struct Sample  <: ExprNode; dist::Symbol; args::Vector{ExprNode} end # dist ∈ DIST_WHITELIST
struct TimeRef <: ExprNode end                                       # @t() -> time(state)
struct Choose  <: ExprNode; alts::Vector{Tuple{Float64,ExprNode}} end

const OP_WHITELIST   = (:+, :-, :*, :/, :^, :>, :<, :>=, :<=, :(==), :&&, :||, :!,
                        :min, :max, :exp, :log, :floor, :ceil, :abs)
const DIST_WHITELIST = (:Poisson, :Binomial, :Normal, :Uniform, :Exponential,
                        :Bernoulli, :LogNormal, :Gamma, :Beta)
const REF_KINDS      = (:species, :param, :obs)
```

The whitelists are closed: `OP_WHITELIST` is arithmetic/comparison only — no `apply`, no `eval`, no `Expr`-head smuggling. `DIST_WHITELIST` is a subset of `Distributions.jl` (`Project.toml:12`); `Poisson` (spawn, `create.jl:151`; event firing, `solvers.jl:321`) and `Binomial` (probability-of-success, `solvers.jl:413`) are the verified runtime draws. Action statements (`pre_action`/`post_action` and event `action`) use a small separate statement set — `SetSpecies{name, ExprNode}`, `SetParams{[(name, ExprNode)]}`, `Log{msg}`, `Seq{[stmt]}` — covering the reserved query verbs `set_params`/`log` (`compilers.jl:67`) without arbitrary code.

### from_json / to_json / validate

`to_expr(::ExprNode)::Expr` lowers a type-checked, whitelisted tree to exactly the species-name/param-name `Expr` that today's authoring macros already produce, then hands it to the unchanged `wrap_fun`/`compile_attrs` (`compilers.jl:148-180`), which positionally substitutes `X -> state.u[i]`, `β -> state.p[:β]` and evals to a closure ONCE at construction. The LLM never emits Julia source; `to_expr` emits Julia from a closed allow-list of interned symbols, so there is no RCE. The runtime hot path is therefore untouched — the format change is confined to the authoring/IR boundary, exactly as ADR 0003 scopes it.

- `from_json(io; seed)::ReactionNetworkProblem` does parse → `validate` → `build_store` → construct. `build_store` fills the ADR-0003 typed-columnar store append-only (preserving position-indexed closures per the runtime facts) and lowers each ExprNode via `to_expr` into the same `SampleableValues` `Expr` the store already holds.
- `to_json(m)::String` round-trips the document via JSON3.
- `validate(spec)::Vector{Diagnostic}` is a pure, eval-free pass enforcing: (1) ExprNode walk — unknown `Ref` names, `op∉OP_WHITELIST`, `dist∉DIST_WHITELIST`, arity; (2) dangling reactant FKs (`r.transition∉ids || r.species∉names`); (3) CONTRACT §5 ranges/sign/integrality (`rate≥0`, `prob_of_success∈[0,1]`, `cycletime≥0`, `capacity≥0`, integer `capacity`/`init`/`stoich` for unstructured species); (4) CONTRACT §1.4 illegal modality combos (`nonblock ⇒ consumed`; `perstep ⇒ !structured && cycletime>0`); (5) CONTRACT §5 A3 TVE policy — a "no" attribute (`init`, `init_uncertainty`, `structured`, `modality`) must be a literal `Const`, not a non-trivial tree.

### Julia mechanism

JSON3.jl + StructTypes.jl for the flat record envelope (`StructTypes.StructType(::Type{ModelSpec}) = StructTypes.Struct()`); hand-rolled custom (de)serialization for the recursive `node`-tagged ExprNode/action union, because StructTypes' built-in abstract-type machinery maps awkwardly onto a recursive `Vector{ExprNode}`. The JSON Schema published for LLM structured output is AUTO-DERIVED from the `ModelSpec`/ExprNode types (the ADR-0003 `const SCHEMA` single-source-of-truth) and round-trip-tested, so schema and loader cannot drift.

### Solutions are a separate concern

The model JSON is the reproducible input; `(model.json, seed)` fully determines a run (CONTRACT §4 D1/D6). The `sol` DataFrame (`solvers.jl:588-591`) and the `log` ledger (`solvers.jl:304-312,659-666`) are OUTPUTS and do NOT live in the model file. They are persisted in a columnar store — Apache Arrow/Parquet (or CSV for human inspection) — replacing the Julia-version-fragile JLD2 blob (`loadsave.jl:205,221`). Recommended layout: `runs/<model_content_hash>/<seed>/{trajectory.arrow, ledger.arrow, run.json}`, where `run.json` is a tiny header `{model_hash, seed, rd_version, tspan, dt}` so a run self-describes and can be matched back to its model. The ledger (cost/reward/valuation per tick) is its own Arrow table — it is what the rNPV/BD demo reads.

## Concrete JSON example

A two-phase pharma pipeline. The first transition carries an expression-valued rate, `Poisson(0.3 * beta) * Preclinical`, serialized as an ExprNode tree (`Call{*}(Sample{Poisson}(Call{*}(Const 0.3, Ref{param} beta)), Ref{species} Preclinical)`):

```json
{
  "rd_format": "reactive-dynamics-model",
  "version": "1.0",
  "meta": {
    "name": "pharma_pipeline_2phase",
    "tspan": 1000.0,
    "dt": 1.0,
    "tunit": 1.0,
    "seed": 1234,
    "alloc_strategy": "weighted"
  },
  "params": [
    { "name": "beta", "value": 0.4 },
    { "name": "discount_rate", "value": 0.1 }
  ],
  "species": [
    { "name": "Preclinical", "init": 10, "cost": 2.0, "reward": 0.0, "valuation": -1.0,
      "modality": { "allocation": "upfront", "return": "consumed", "blocking": "block" } },
    { "name": "Phase1",   "init": 0, "cost": 5.0, "valuation": -2.0 },
    { "name": "Approved", "init": 0, "reward": 100.0, "valuation": 50.0 }
  ],
  "transitions": [
    { "id": "t_pre_to_p1", "name": "preclinical->phase1",
      "rate": { "node": "call", "op": "*", "args": [
        { "node": "sample", "dist": "Poisson", "args": [
          { "node": "call", "op": "*", "args": [
            { "node": "const", "value": 0.3 },
            { "node": "ref", "kind": "param", "name": "beta" } ] } ] },
        { "node": "ref", "kind": "species", "name": "Preclinical" } ] },
      "rate_mode": "poisson",
      "cycletime":       { "node": "const", "value": 180.0 },
      "prob_of_success": { "node": "const", "value": 0.6 },
      "capacity":        { "node": "const", "value": 5 },
      "priority":        { "node": "const", "value": 1 },
      "max_lifetime":    { "node": "const", "value": 730.0 },
      "pre_action": null, "post_action": null },
    { "id": "t_p1_to_appr", "name": "phase1->approved",
      "rate":            { "node": "const", "value": 0.05 },
      "rate_mode": "poisson",
      "cycletime":       { "node": "const", "value": 365.0 },
      "prob_of_success": { "node": "const", "value": 0.3 } }
  ],
  "reactants": [
    { "transition": "t_pre_to_p1", "species": "Preclinical", "side": "lhs",
      "stoich": { "node": "const", "value": 1 },
      "modality": { "allocation": "upfront", "return": "consumed", "blocking": "block" } },
    { "transition": "t_pre_to_p1", "species": "Phase1",      "side": "rhs",
      "stoich": { "node": "const", "value": 1 } },
    { "transition": "t_p1_to_appr", "species": "Phase1",     "side": "lhs",
      "stoich": { "node": "const", "value": 1 } },
    { "transition": "t_p1_to_appr", "species": "Approved",   "side": "rhs",
      "stoich": { "node": "const", "value": 1 } }
  ],
  "observables": [
    { "name": "active_p1", "every": 30.0,
      "on":    [ { "node": "call", "op": ">", "args": [ { "node": "timeref" }, { "node": "const", "value": 0 } ] } ],
      "range": [ { "weight": 1.0, "value": { "node": "ref", "kind": "species", "name": "Phase1" } } ] }
  ],
  "events": [
    { "id": "acquisition_lever",
      "trigger": { "node": "call", "op": "&&", "args": [
        { "node": "call", "op": ">", "args": [ { "node": "timeref" }, { "node": "const", "value": 500 } ] },
        { "node": "call", "op": "<", "args": [ { "node": "timeref" }, { "node": "const", "value": 501 } ] } ] },
      "action": { "node": "set_species", "name": "Preclinical", "value": { "node": "const", "value": 20 } } }
  ]
}
```

The rate field above is the canonical demonstration that an expression-valued attribute serializes as a typed tree, never as a Julia source string: the loader walks it, type-checks `op:"*"` against `OP_WHITELIST` and `dist:"Poisson"` against `DIST_WHITELIST`, lowers it via `to_expr` to `:(rand(Poisson(max(state.dt * (0.3 * beta), 0))) * Preclinical)` (the same `Expr` that `expand_rate` produces at `create.jl:150-155`), and hands that to the untouched `wrap_fun`.

## Consequences

- Closes the import-time RCE. The three `eval`/`Meta.parse` sites (`loadsave.jl:65,72`; `ReactiveDynamics.jl:99,101,102`) are deleted; a model file is now inert data, not a Julia program. A malicious or malformed model can no longer execute code on load.
- Enables LLM emission and self-validation. The JSON Schema (auto-derived from `ModelSpec`) is a valid structured-output / tool-call `input_schema`, so an LLM is constrained to emit conformant models; `validate(spec)` is the eval-free self-check before `from_json`. This directly serves the agentic north-star (BD acquisition-lever demo as a mid-sim event).
- One format, two roles, no zoo. TOML, batched CSV, and JLD2 model paths are removed (`loadsave.jl`), cutting three serializers to one. Human-diffable, round-trippable, schema-validatable.
- Runtime untouched. ExprNodes lower to the existing `SampleableValues` `Expr`; `compile_attrs`/`wrap_fun` (`compilers.jl:148-180`) and the whole hot path are unchanged, so this carries no engine-performance regression risk.
- Costs: a closed whitelist REJECTS some currently expressible models (arbitrary `@register`d user functions, `create.jl`-style raw call exprs). This is intended and is the eval-free trade; the cutover must enumerate dropped idioms (`@choose`/`@move`/`@structured`/`set_params`/`log` are covered by explicit nodes; arbitrary registered bodies are not). A temporary translation layer maps the 3-axis JSON modality back to the legacy `Set{Symbol}` the engine still consumes (`ReactiveDynamics.jl:144`, `state.jl:204`) until the CONTRACT §1 re-model lands — a coupling that is a regression risk if the two mappings disagree.

## Open questions

- Action statement coverage. Is `{SetSpecies, SetParams, Log, Seq}` sufficient for real pre/post and event actions, or do existing tutorials (`interface/update.jl` `@register` bodies) need richer statements? Under-scoping breaks the BD acquisition lever if it needs more than `set_species`/`set_params`. Audit tutorials before freezing the statement set. **ANSWERED:** the action family was extended — `SetField`/`@advance` (ADR 0008), `AddToken`/`Activate`/`Deactivate` (ADR 0010 §C), and `SetTokens` + the general-code `Invoke` escape hatch (ADR 0011) — to the canonical `{SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}`. `Invoke` (and registered value-helpers, §C) cover arbitrary host logic by name without reopening eval (§S4 holds).
- Integer vs float in JSON. JSON has one number type; `stoich`/`capacity`/`init` that must be integer (`solvers.jl:413` `Int` cast) need `validate()` to enforce integrality, since emitters (LLMs especially) often write `1.0` for an integer. `validate` must coerce/check, not trust the parsed number subtype.
- Mutation-during-simulation patch shape (maintainer ruling 2). `from_json` yields a static spec, but add-reaction-mid-run must re-derive `wrap_fun` against a refreshed varmap when a new transition references a NEW species/param (the construction-time varmap is frozen, `compilers.jl:148-153`). A JSON patch document (add-only species/transitions/reactants, never reorder/delete, to preserve position-indexed `state.u[i]`) should be specified alongside the full-document shape, with an `apply_patch(spec_delta)` API.
- `@register` user functions. ADR 0003 already flags the `registered` path (`loadsave.jl:72`) for removal. Decide between outright removal and a named-function registry keyed by a closed whitelist — the latter preserves some legacy expressivity without re-opening eval.
- Genesis primitive (maintainer open question 5). The rate ExprNode currently lowers to `Poisson(dt·rate)` per `expand_rate` (`create.jl:150-155`), with a `deterministic` mode. For business processes other genesis modes (scheduled/cohort arrivals, capacity-gated pull, fixed batches) may be needed; these would be additional `rate_mode` values, orthogonal to the serialization format decided here.
- Schema/loader drift mitigation. Confirm the JSON Schema is auto-derived from `ModelSpec` and round-trip-tested at build time rather than hand-maintained as a separate `.schema.json`.
