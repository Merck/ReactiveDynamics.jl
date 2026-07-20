# The JSON model document (`model.rdj.json`)

A ReactiveDynamics model is **data**, not code. Its canonical serialized form is a single JSON document — conventionally `model.rdj.json` — that is at once the on-disk model file and the artifact an agent can emit under structured output and self-check before loading. Nothing in the document is ever `Meta.parse`d or `eval`d on load: every expression is a node-tagged tree drawn from a closed whitelist, and every host construct a model needs — a structured-token type, a custom value helper, a callback body — is referenced **by name** against a per-network *registry* that the host program populates, never as source text carried in the file. Loading is therefore inert: a malformed or hostile document can produce diagnostics or a construction error, but it cannot execute arbitrary code. The full rationale for this boundary lives in the [Serialization deep-dive](../deep_dives/serialization.md) and [ADR 0005](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0005-serialization-json-ir.md); this reference documents the **shape** of the document and the exact keys the loader reads and the exporter emits.

The two entry points are [`from_json_model`](@ref) (JSON string → constructed `ReactionNetworkProblem`, via parse → `validate` → build → construct) and [`to_json_model`](@ref) (a `ReactionNetwork` or a constructed `ReactionNetworkProblem` → JSON string). `validate` is a pure, eval-free static pass returning a `Vector{Diagnostic}`; `from_json_model` gates construction on it being empty. The ground truth for the schema is the round-trip code in `src/serialize.jl` (`build_network_from_dict` reads, `model_to_dict` writes) — where this page and the prose in [ADR 0005](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0005-serialization-json-ir.md) diverge, the code is authoritative.

## Top-level envelope

The document is one JSON object. Two string keys tag the format; the rest are the model's arrays and its `meta` bag.

| Key | Type | Role |
|---|---|---|
| `rd_format` | string | Always `"reactive-dynamics-model"` (format tag; emitted, informational). |
| `version` | string | Schema version, currently `"1.0"`. |
| `meta` | object | Solver settings + free-form keywords (see below). `meta.tspan` is **required** by `from_json_model`. |
| `params` | array | Named scalar parameters (`:P` rows). |
| `species` | array | Resource pools and structured-token kinds (`:S` rows). |
| `transitions` | array | The recipes (`:T` rows): rate, cycle time, success probability, capacity, actions. |
| `reactants` | array | The promoted incidence table: one row per LHS/RHS term, FK'd to a transition. |
| `observables` | array | Folded observables (`:obs` rows). Optional. |
| `rules` | array | The endogenous decision channel — guarded actions (ADR 0010). Optional. |
| `inputs` | array | Declared external read ports + their pre-wire defaults (ADR 0012). Optional. |
| `population` | array | Declarative initial structured-token instances (ADR 0007). Optional; validated. |

Two notes where the code differs from older prose. First, the endogenous channel the loader consumes is `rules[]`, not `events[]`: a legacy `:E` event is lifted to a `RawExpr`-action `Rule` at construction, and because `RawExpr` is the non-typed bridge (deliberately not JSON-serializable), the exporter emits only typed rules and never an `events[]` array. Second, `meta.alloc_strategy` (and its aliases `strategy`/`schedule`) is **accepted-and-ignored**: ADR 0002 made priority-weighted progressive filling the single allocation policy, so there is no longer a strategy switch to honor — the key round-trips harmlessly.

### `meta`

`meta` merges solver settings and free-form model keywords into one bag (they already merge into a single `keywords` bag at construction). Recognized keys:

- `tspan` — the simulation horizon (a scalar `tend`; **required** — `from_json_model` errors without it).
- `dt` — the fixed time step.
- `tunit`, `seed` — the time unit and the RNG seed (a `seed` in `meta` is used when the caller passes none to `from_json_model`).
- `alloc_strategy` / `strategy` / `schedule` — legacy no-op keys, symbolized and ignored (see above).

`to_json_model(::ReactionNetworkProblem)` reconstructs `meta` from the constructed model's solver fields (`tspan[2]`, `dt`, and `seed` when set) so the exported document is complete and re-importable on its own; a caller-supplied `meta` always takes precedence.

### `params[]`

Each entry is `{ "name": <string>, "value": <number> }`. Values are plain JSON numbers, never expression strings — this is the direct replacement for the old `eval`-on-load parameter path. `value` becomes the `prmVal` of a `:P` row.

### `species[]`

Each entry carries a `name` and optional scalar attributes; all of these are **time-invariant literals** (validate rule 5 rejects a non-trivial expression tree here):

- `init` — the initial count / marking at t = 0 (`specInitVal`; omitted when 0).
- `cost`, `reward`, `valuation` — the per-unit ledger coefficients (omitted when 0.0, their default).
- `structured` — `true` marks the species a **structured-token kind** (its instances are live agents, not a counted `Float64` stock); omitted otherwise.
- `modality` — the resource-claim modality as a 3-axis object (see [Modality](#Modality-3-axis)); omitted when the modality set is empty.

### `transitions[]`

Each transition is the recipe for one reaction. Keys:

- `id` — the identifier the `reactants[]` rows FK against (in an exported document this is the transition's `name`, e.g. `"adv_discovery"`, falling back to a positional `"t<i>"` for an unnamed transition).
- `name` — the display name (optional).
- `rate` — the **bare firing intensity** as an ExprNode (or a bare literal). It is *not* pre-wrapped: the loader applies the genesis wrapping itself per `rate_mode`.
- `rate_mode` — `"poisson"` (default) wraps the bare rate as a per-tick Poisson draw `rand(rng, Poisson(max(dt * rate, 0)))`; `"deterministic"` uses the bare rate as-is (the `@deterministic` path).
- `cycletime`, `prob_of_success`, `capacity`, `priority`, `max_lifetime`, `multiplier` — ExprNode-valued (or bare-literal) attributes; each is emitted only when it differs from its construction default, so a minimal document omits them.
- `pre_action`, `post_action` — optional typed action-statement trees run before/after firing (see [Action statements](#Action-statements)). Only typed `ActionStmt`s serialize; a raw-`Expr` action is dropped on export.

### `reactants[]`

This is the promoted first-class incidence table — the transition↔reactant relation made explicit on disk (rather than re-parsed from the reaction line at runtime). Each row is one LHS or RHS term:

- `transition` — FK to a `transitions[].id` (a dangling FK is validate rule 2).
- `side` — `"lhs"` (a consumed/held reactant) or `"rhs"` (an emitted product).
- Exactly one *atom* form:
  - `species` — a plain species name (the ordinary case).
  - `predicate` — an LHS token filter `@select(Kind, …)`: `{ "kind": <structured species>, "clauses": [[field, op, value], …] }`, ops drawn from `PRED_OP_WHITELIST` (`== != < <= > >= in`), AND-joined. See [Predicates](#Predicates).
  - `advance` — an RHS lifecycle field-write `@advance`: `{ "field": <string>, "value": <node|literal> }`.
  - `structured` — an RHS **named genesis product** `@structured(:Kind, field = …)`: `{ "kind": <registered structured species>, "fields": [{ "name", "value" }, …] }` (the raw-constructor form was removed so serialization is total).
  - `move` — a species relabel `@move`: `{ "from": <species>, "to": <species> }`.
- `stoich` — the integer stoichiometric coefficient (or an ExprNode for the rare expression-valued case); omitted when 1.
- `modality` — the 3-axis modality (LHS terms only); omitted when empty.

## The ExprNode IR

Every expression-valued attribute — a rate, cycle time, success probability, observable range, guard, predicate-clause value, action value — is serialized as a **node-tagged tree**, never as a Julia source string. A scalar attribute may also be authored as a bare JSON number, bool, or string, which the loader normalizes to a constant. The tree is a closed tagged union keyed on the `"node"` tag:

| `node` | Fields | Meaning |
|---|---|---|
| `const` | `value`, `symbol` (bool) | A literal. `symbol: true` marks a `value` that is a Symbol (e.g. a phase name `"Phase2"`) rather than a numeric/string scalar. |
| `ref` | `kind`, `name` | A named reference; `kind ∈ REF_KINDS` = `species`, `param`, `obs`. |
| `call` | `op`, `args[]` | An operator application; `op ∈ OP_WHITELIST` (arithmetic/comparison/boolean only). |
| `sample` | `dist`, `args[]` | A distribution draw; `dist ∈ DIST_WHITELIST` (`Poisson`, `Binomial`, `Normal`, …). |
| `timeref` | — | The current simulation time `@t()`. |
| `choose` | `alts[]` | A weighted choice; `alts` is a list of `{ "weight", "value" }`. |
| `field` | `name` | A read of the firing instance's own bound-token field (legal only in a `SetField`/`@advance` value). |
| `externalref` | `port` | A read of a declared `inputs[]` port's latched value (ADR 0012). |

The whitelists are closed and are the trust boundary: `OP_WHITELIST` is `+ - * / ^ > < >= <= == != && || ! min max exp log floor ceil abs` — arithmetic and comparison only, with no `apply`, no `eval`, and no `Expr`-head smuggling — and `DIST_WHITELIST` is `Poisson Binomial Normal Uniform Exponential Bernoulli LogNormal Gamma Beta`. On load, each node is lowered by `to_expr` to exactly the species-name/param-name `Expr` the authoring macros produce, then handed to the unchanged compile step; an out-of-whitelist `op`/`dist` or an undeclared `ref` name is a `validate` diagnostic.

## Modality (3-axis)

An LHS reactant's (or a species') resource-claim modality serializes as an orthogonal 3-axis object rather than an unvalidated symbol set:

```json
{ "allocation": "upfront" | "perstep",
  "return":     "consumed" | "conserved",
  "blocking":   "block" | "nonblock" }
```

The loader translates this to the internal `Set{Symbol}` the engine consumes (and back on export). Only the five legal combinations are accepted: `blocking = nonblock` requires `return = consumed` (validate rule 4 rejects the illegal `nonblock + conserved` combination). All three axes default (`upfront` / `consumed` / `block`), so a `block`-`upfront`-`consumed` claim needs no `modality` key at all.

## Observables

Each `observables[]` entry is a folded observable: `{ "name", "every": <number>, "on": [<node>, …], "range": [{ "weight", "value": <node> }, …] }`. `every` is the sampling period, `on` a list of gating conditions, and `range` the weighted value expressions folded at each sample. All expression slots are ExprNodes.

## Rules and action statements

`rules[]` is the endogenous decision channel (ADR 0010): guarded actions that fire on state- or time-contingent conditions inside the model. Each rule is `{ "id", "guard": <node>, "action": <stmt>, "fire_mode": <string> }` — a guard ExprNode (a boolean fires once; a numeric fires `Poisson(v)` times), a typed action statement, and a firing mode (`"every_tick"` by default).

### Action statements

An action statement is tagged by a `"verb"` from the closed `ACTION_VERBS` set. Its value slots are ExprNodes (or bare literals):

| `verb` | Fields | Effect |
|---|---|---|
| `set_species` | `name`, `value`, `mode` | Set (`mode: "set"`) or increment (`"inc"`) a plain-species count. |
| `set_params` | `assigns: [{name, value}, …]` | Reassign named parameters. |
| `set_field` | `field`, `value` | Write a field of the firing instance's bound token (legal only in a transition post-action, not a `Rule`). |
| `set_tokens` | `predicate`, `assigns` | Write `assigns` over every token matching a predicate (the population generalization of `set_field`). |
| `add_token` | `kind`, `fields: [{name, value}, …]` | Instantiate a new structured token of a registered `kind` (the in-model "acquisition"). |
| `activate` / `deactivate` | `transition` | Flip a transition's activation (soft; in-flight instances still finish). |
| `invoke` | `fn`, `args[]` | The general-code escape hatch — call a registered host function by NAME (its body is host Julia, never file bytes). |
| `log` | `msg` | Emit a log entry. |
| `seq` | `stmts: [<stmt>, …]` | Run a sequence of statements. |

An unregistered `add_token.kind` or `invoke.fn`, an out-of-set verb, a `set_field` in a `Rule`, or a `ref` to an undeclared name is a `validate` diagnostic — never an eval. `RawExpr` (the legacy non-typed bridge) has no verb and is intentionally not serializable.

### Predicates

A `@select` predicate — used in a `reactants[]` LHS `predicate` and in a `set_tokens` action — is `{ "kind": <structured species>, "clauses": [[field, op, value], …] }`. Each clause is a `[field, op, value]` triple with `op ∈ PRED_OP_WHITELIST` (`== != < <= > >= in`); clause `value`s must be time-measurable (a `sample`/RNG node is rejected in a predicate — the filtration must not depend on draw order).

## Other arrays

- `inputs[]` — declared external read ports (ADR 0012): `{ "port": <string>, "default": <Const|literal> }`. The `default` seeds the external-input buffer before any host `add_wire!` has delivered; it must be a literal value, and an `externalref` node to an undeclared port is a diagnostic. The wiring topology itself is host-side, never in the document.
- `population[]` — the declarative initial structured-token marking (ADR 0007): entries naming a structured `species` and a registered `kind`, so a structured model's initial state is in the document (making `(model, seed)` reproducible). The loader validates these; the instances are built by the initial-marking machinery.

## A worked example

The excerpt below is trimmed from the BD acquisition demo model at `demo/bd_acquisition/model.rdj.json`. It shows: an envelope with a minimal `meta`; two plain params; a `structured` species (`Project`) alongside two plain resource pools (`scientist`, `budget`); a mix of literal and ExprNode transition attributes; and reactant rows exercising the `predicate`, `species`+`modality`, and `advance` atom forms.

```json
{
  "rd_format": "reactive-dynamics-model",
  "version": "1.0",
  "meta": { "tspan": 40.0, "dt": 1.0, "alloc_strategy": "weighted" },
  "params": [
    { "name": "synergy_pos", "value": 0 },
    { "name": "synergy_eff", "value": 0 }
  ],
  "species": [
    { "name": "Project", "structured": true },
    { "name": "scientist", "init": 40 },
    { "name": "budget", "init": 150 }
  ],
  "transitions": [
    { "id": "adv_discovery", "name": "adv_discovery", "rate": 2.0, "rate_mode": "deterministic",
      "cycletime": 1.0, "prob_of_success": 0.45, "priority": 1.0 },
    { "id": "adv_phase2", "name": "adv_phase2", "rate": 2.0, "rate_mode": "deterministic",
      "cycletime": { "node": "call", "op": "-", "args": [ { "node": "const", "value": 2.0 },
                     { "node": "call", "op": "*", "args": [ { "node": "const", "value": 0.5 },
                       { "node": "ref", "kind": "param", "name": "synergy_eff" } ] } ] },
      "prob_of_success": { "node": "call", "op": "+", "args": [ { "node": "const", "value": 0.4 },
                     { "node": "call", "op": "*", "args": [ { "node": "const", "value": 0.2 },
                       { "node": "ref", "kind": "param", "name": "synergy_pos" } ] } ] },
      "priority": 2.0 },
    { "id": "financing", "name": "financing", "rate": 16.0, "rate_mode": "deterministic" }
  ],
  "reactants": [
    { "transition": "adv_discovery", "side": "lhs",
      "predicate": { "kind": "Project", "clauses": [ ["phase", "==", "Discovery"] ] } },
    { "transition": "adv_discovery", "side": "lhs", "species": "scientist", "stoich": 2,
      "modality": { "allocation": "upfront", "return": "conserved", "blocking": "block" } },
    { "transition": "adv_discovery", "side": "lhs", "species": "budget", "stoich": 2,
      "modality": { "allocation": "perstep", "return": "consumed", "blocking": "block" } },
    { "transition": "adv_discovery", "side": "rhs",
      "advance": { "field": "phase", "value": "Phase1" } },

    { "transition": "financing", "side": "rhs", "species": "budget", "stoich": 1 }
  ]
}
```

Reading the `adv_phase2` transition: its `prob_of_success` is not a fixed number but the tree `0.4 + 0.2 * synergy_pos` — a `call{+}` of a `const 0.4` and a `call{*}` of `const 0.2` and `ref{param} synergy_pos`. On load this lowers to the same `Expr` the DSL would produce, so the acquisition lever (setting `synergy_pos`/`synergy_eff` via a param) reprices the pipeline with no code carried in the file. The `adv_discovery` reactant rows show the three atom shapes: an LHS `predicate` picking `Project` tokens in the `Discovery` phase, two plain-species claims with explicit modalities (`scientist` held-and-returned `conserved`; `budget` spent `perstep`/`consumed`), and an RHS `advance` writing the token's `phase` field to `Phase1`.

## Why eval-free

No field of a `model.rdj.json` document is ever `Meta.parse`d or `eval`d on load — the import-time remote-code-execution surface of the old TOML/CSV loaders is closed, and a model file is inert data. The only place Julia source is produced is `to_expr`, which emits from the closed whitelist of interned symbols and hands the result to a single construction-time compile; there is no `Expr`-head smuggling and no `apply`/`eval` op. The two ways a model reaches arbitrary host logic — a registered value helper and an `invoke` callback — do not weaken this: the file carries only a NAME plus ExprNode args, lowering to `registry[name](…)`, so the trust boundary is the host program that populated the registry, not the document. An `invoke` body is the one trusted-but-unverified tier (prefer the declarative verbs, which `validate` can check statically). The API entry points — `from_json_model`/`to_json_model`, `validate`, the `*_to_dict`/`*_from_dict` helpers, and the `ExprNode` IR types — are documented on the [Serialization reference](serialization.md) page; the full argument and its mapping to the runtime compile path are in the [Serialization deep-dive](../deep_dives/serialization.md) and [ADR 0005](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0005-serialization-json-ir.md).
