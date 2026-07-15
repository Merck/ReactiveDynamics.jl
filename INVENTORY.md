# ReactiveDynamics.jl — Ground-Truth Inventory (branch `ref-agents`)

## Orientation

ReactiveDynamics.jl (RD) is NOT a chemical reaction network despite its DSL surface and Catalyst-derived parser. It is a timed, stochastic, resource-constrained Petri net / discrete-event system: transitions are stateful recipes that spawn instances at a Poisson rate, occupy resources (species) over a cycle time, and complete with a terminal Binomial probability-of-success that emits RHS products; resources carry a modality (`:nonblock`/`:conserved`/`:rate`) governing allocation, and a cost/reward/valuation ledger accumulates into a per-step `log`. On the `ref-agents` branch the package has fully pivoted from the old SciML embedding: the simulation output is now the NATIVE engine type `ReactionNetworkProblem` (an AlgebraicAgents `@aagent`, not a `DiscreteProblem`); AlgebraicAgents (AA) is reexported from `state.jl` and stepping is driven by AA's `simulate`/`_step!`; the ACSet schema/macro were renamed to `ReactionNetworkSchema`/`@ReactionNetworkSchema` and migrated from a Catlab `@present FreeSchema` to an ACSets `BasicSchema`; a new structured-token agent layer (`interface/agents.jl`, `plots.jl`) was added and `optim.jl` was removed. This document is the canonical map of current source and SUPERSEDES the architecture map in `REVIEW.md`, which was written against an earlier (SciML/`DiscreteProblem`) snapshot — every `file:line` here is verified against the current `ref-agents` tree, and REVIEW.md's line numbers must not be trusted.

Section index: [Module Map](#module-map) · [Public-API Audit](#public-api-audit) · [ACSet Schema](#acset-schema) · [Stepping Trace](#stepping-trace) · [AlgebraicAgents Touchpoints](#algebraicagents-touchpoints) · [What Changed Since REVIEW.md](#what-changed-since-reviewmd) · [REVIEW.md Defect Reconciliation](#reviewmd-defect-reconciliation)

---

## Module Map

The module root `src/ReactiveDynamics.jl` declares the package deps (`ACSets`, `Reexport`, `MacroTools`, `ComponentArrays`, and `@reexport using GeneratedExpressions`), defines the core type aliases, the ACSet schema, the construction/merge helpers, and orchestrates the include order. It has NO `export` of its own. Include order (`ReactiveDynamics.jl:222-247`, grown since the original snapshot): `state.jl`, `exprnode.jl`, `compilers.jl` are explicit, THEN `interface/`, `utils/`, `operators/` are bulk-included via `include.(readdir(...; join=true))` (filesystem-alphabetical), then `solvers.jl`, `predicates.jl`, `actions.jl`, `ledger.jl`, `serialize.jl`, `#include("optim.jl")` (commented — file removed), `loadsave.jl`, and finally the Phase-0.6 analysis/viz trio `analysis.jl`/`export.jl`/`visualize.jl` (placed last so they see the run-state, ledger, predicate/sortkey, and serializer). The readdir bulk-includes have a real-but-currently-benign load-order dependency: `interface/agents.jl` (first alphabetically) subtypes `AbstractAlgebraicAgent` and names `ReactiveDynamics.Transition`, both of which exist only because `state.jl` was included first at line 222. Adding a new alphabetically-earlier interface/operator file with top-level executable code that calls a later sibling would break the build.

| Load # | File | Lines | Purpose | ACSet | AA | Stepping |
|---|---|---|---|---|---|---|
| (root) | `src/ReactiveDynamics.jl` | 249 | Module root: deps/reexports, type aliases, the `BasicSchema` + concrete acset type, `Base.convert` coercions, schema metadata tables, `FoldedObservable`, and `assign_defaults!`/`add_obs!`/`merge_acs!`/constructors. Orchestrates includes (`:222-247`). | Schema home (defines it) | reexports `GeneratedExpressions` only; NOT AA | none (static construction) |
| 222 | `src/state.jl` | 366 | Live simulation state and sub-structs: `ReactionNetworkProblem` (top-level `@aagent`), in-flight `Transition`, `Observable` sampler, `UnfoldedReactant`. Query/indexing API, ACSet pass-throughs, observable resampling, `sample_transitions!`. (Grown well past the 258 of the original snapshot — RNG field, structured-token sync, `sol`/`log` reshape.) | heavy reader via `state.acs` | `@reexport using AlgebraicAgents` (line 1); 3 `@aagent` structs | provides per-step machinery (`update_observables`, `sample_transitions!`, `save!`) |
| 223 | `src/exprnode.jl` | 181 | **NEW (Phase-1, ADR 0005):** the typed, eval-free expression IR. A closed tagged-union `ExprNode` tree is the authoring/serialization form of every time-varying attribute (rate/stoich/cycletime/PoS/priority/action value); `to_expr` lowers a tree to EXACTLY the `Expr` the DSL produces (then compiled ONCE by `compilers.jl`), `from_expr` is the structural inverse. The closed whitelists here are the only Julia ever produced from an inert model document — the runtime hot path is untouched. | none (operates on `Expr`, not the acset) | none | none (authoring/serialization boundary only) |
| 224 | `src/compilers.jl` | 202 | Expression-compilation layer: rewrites DSL attribute exprs (rates/stoich/actions/observables) into runnable `(state, transition)` closures; species→`u[i]`, params→`p[:name]`, query macros→state calls; bulk `compile_attrs`. | reads schema attrs | none | setup (builds `transition_recipes`); `wrap_expr` also fires per-step via `state.wrap_fun` |
| 225 | `src/interface/agents.jl` | 169 | NEW structured-token agent layer: `AbstractStructuredToken`/`BaseStructuredToken`, `@structured_token`, `register_structured_species!`/`add_structured_token!`, and the (unexported) binding/priority/species accessor protocol. (Grown past the 72 of the original snapshot.) | `register_structured_species!` mutates `:S`/`:specStructured` | heavy: `@aagent FreeAgent`, no-op `_step!`/`_projected_to` | supplies the protocol the loop calls; tokens are clock-inert |
| 225 | `src/interface/create.jl` | 333 | CREATE half of the DSL: `@ReactionNetworkSchema` + `get_data`/`get_data!` parser pipeline (arrow normalization, rate→Poisson, keyword-attr extraction, event expansion, `@register`/`@take` rewriting); `@append_transitions`. | emits `ReactionNetworkSchema(get_data(...)...)`; uses `prettynames`/`defargs` | none | parse-time; `expand_rate` AUTHORS the spawn-rate Poisson expr |
| 225 | `src/interface/parsing_utils.jl` | 41 | Low-level AST helpers (from Catalyst): tuple arity/index, multiplication flattening (`multiplex`), and `recursively_expand_dots` (alias of `underscorize`). | none | none | none (parse-time) |
| 225 | `src/interface/plots.jl` | 11 | **Phase-0.6:** reduced to a comment-only stub — `plot_df` deleted and `_draw` RELOCATED to `ext/RDPlotsExt.jl` (Plots demoted to a weakdep). | none | — | none |
| 225 | `src/interface/reaction_parser.jl` | 117 | Reactant extraction (from Catalyst): resolves `@choose`, escapes refs, folds a reaction side into `FoldedReactant`s carrying species/stoich/modality. | reads `state[:, :specName]` | indirect (consumes the `@aagent` state's getindex) | sampling/parse-support; `recursively_choose` re-randomizes per call |
| 225 | `src/interface/aa_coupling.jl` | 138 | **NEW (Phase-1, ADR 0012 / CONTRACT §13):** the AlgebraicAgents READ/COUPLING surface a bare `entangle!` did not provide — RD-level `getobservable`/wire bridge for bidirectional observables between a `ReactionNetworkProblem` and sibling AA agents (explicit one-tick-lag Jacobi coupling). | reads run-state observables | heavy: implements AA observable/wire hooks | read surface; coupling interleaves at the AA time gate |
| 225 | `src/interface/checkpoint.jl` | 116 | **NEW (Phase-1, ADR 0007 §C / CONTRACT §10.5):** state dump/restore — serialize a live run at a TICK BOUNDARY into an eval-free JSON-representable artifact (declarative initial marking + dynamic run state: clock, RNG, creation counters, `u`, the token population with current field values, rule latches) and reconstruct it. A zero-tick dump with no in-flight transitions IS an initial marking (the "initial marking = checkpoint" identity). | reads/reconstructs the acset + run-state | dumps/restores the AA token population | boundary-only (not in the hot loop) |
| 225 | `src/interface/solve.jl` | ~55 | **Phase-0.6 + WS-4 + WS-B2:** the dead SciML `@plot`/`plot_summary`/`plot_ensemble_sol`/`first_sol` were DELETED (ADR 0014 §A1); WS-4 then deleted the dangling `export @agentize`. **WS-B2 re-added `@agentize` as thin constructor sugar** — it expands to exactly one `ReactionNetworkProblem(acs[, u0, p]; …)` call (no second construction path), with auto-naming from a bare-binding acs (`@agentize net` ⇒ `name = "net"`). The `export @agentize` is co-located here (solve.jl:1). | none | expands to the constructor call | none |
| 225 | `src/interface/update.jl` | 518 | UPDATE half of the DSL: mutation/config macros (`@push`, `@add_species`, `@mode`, `@prob_init`/`@prob_params`/`@prob_meta`, `@cost`/`@reward`/`@valuation`, `@periodic`/`@jump`, `@aka`, `@register`, ...) that reshape syntax into ACSet mutations. **WS-4:** removed the dangling `@prob_role`/`@list_by_role`/`@list_roles` exports and the vestigial `@prob_check_verbose` macro (it called an undefined `check_params`). | heavy reader/mutator | none | seeds step inputs (`@jump`/`@periodic`/`@prob_meta`) |
| 226 | `src/utils/utils.jl` | 102 | AST/kwarg grab-bag (macroname, blockize, underscorize, args_kwargs, ...) + one dead engine helper `get_bound_agent`. NO exports. | none | only via dead `get_bound_agent` | none |
| 227 | `src/operators/equalize.jl` | 84 | Species identification/collapse: `equalize!` merges matched species rows + rewrites refs; `@equalize` parses `A = B = C` blocks. | heavy mutator | none | none (model transform) |
| 227 | `src/operators/joins.jl` | 238 | Model composition: `union_acs!` merges one network into another in place; `@join` builds + folds models. | heavy reader/mutator | none | none (model transform) |
| 228 | `src/solvers.jl` | 1067 | CORE discrete-event engine: `ReactionNetworkProblem` constructor + the per-tick step loop (spawn/allocate/advance/finish), structured-token RHS, cost/reward/valuation ledger, tspan/dt resolution, and the AA `_step!`/`_reinit!`/`_projected_to` contract. (Grown far past the 685 of the original snapshot — progressive-fill allocator, RNG threading, rule/predicate hooks.) | heavy reader/writer | implements ReactionNetworkProblem's AA stepping contract | THIS IS THE STEP LOOP |
| 229 | `src/predicates.jl` | 108 | **NEW (Phase-1, ADR 0008 / CONTRACT §9.5):** token filtration — a `TokenPredicate{kind, clauses}` narrows the candidate set a structured LHS reactant binds (in front of the priority/creation-index sort + integer take) and powers the population-level `SetTokens` write. Phase is an ATTRIBUTE (`@select(Project, phase==:Phase2)` selects, `@advance` advances) — no per-phase kinds. | reads token/species attrs | operates over structured-token agents | binding-time filter in the step loop |
| 230 | `src/actions.jl` | 253 | **NEW (Phase-1, ADR 0010/0011 / CONTRACT §12):** the endogenous decision channel — the closed action family `{SetSpecies, SetParams, AddToken, Activate, Deactivate, Log, Seq, ...}` + the `Rule` record (guard, action, fire_mode) evaluated once/tick at step 10, plus the stateless transition `guard` AND-ed with the latching `transActivated` gate. Action values evaluate through the RNG-threaded `context_eval` hot path. | reads/writes species/params/tokens | none directly | fires at a fixed point in `_step!` (the repaired event channel) |
| 234 | `src/ledger.jl` | 288 | **NEW (Phase-2, MVP finding D/D-bis):** the per-program (per-structured-token) ledger — attributes cost at the bind site (`evolve!`) and reward at the finishing transition (`finish!`), split across a transition's bound programs, with a network UNATTRIBUTED bucket for plain spend/reward; enforces the SUM-CONSISTENCY invariant vs the aggregate §8.5 rows. Query API `program_ledger(state)::DataFrame` / `program_ledger_entries(state, name)`. Placed after `actions.jl` so it sees the token types/accessors; hooks in `solvers.jl` call into it. | reads run-state | reads structured-token agents | accrues economics onto in-flight programs during stepping |
| 235 | `src/serialize.jl` | 919 | **NEW (Phase-1, ADR 0005):** single-JSON, eval-free model serialization. `from_json_model` parses one JSON object (`meta` + `params`/`species`/`transitions`/`reactants`/`observables`/`events` arrays of node-tagged `ExprNode` dicts) → validates (`_validate_node!`/`_validate_action!`) → lowers each node via `to_expr` into the same `Expr` columns the DSL fills → constructs a `ReactionNetworkProblem`. NEVER `Meta.parse`s/`eval`s a model field (closes the import-time RCE). Hand-rolled node-tagged (de)serialization on JSON.jl. | writes the acset from parsed IR | constructs the AA problem | none (authoring/serialization boundary) |
| 237 | `src/loadsave.jl` | 72 | Model + solution persistence, now a THIN eval-free wrapper (ADR 0005): `@import_model`/`@export_model` delegate to the `src/serialize.jl` typed-IR loader (no `Meta.parse`/`eval` of model fields — the old TOML/CSV/JLD2 loader + its import-time RCE are REMOVED); `@export_solution_as_table`/`@export_solution_as_csv` write the solution DataFrame for human inspection. | reader/writer via serialize.jl | none | none (I/O) |
| 245 | `src/analysis.jl` | 492 | **NEW (Phase-0.6, ADR 0013 §14):** `push_token_trajectory_row!` + the `token_trajectory`/`representative_token`/`trajectory_envelope` read API; `ensemble`/`summarize`/`treatment_effect` + `EnsembleProblem <: AbstractAlgebraicAgent` (with `observables`/`getobservable`); the §15.1 recipe wrapper types + `log_scalar_series`/`log_count_series`. | reads `program_ledger`/`prob.sol` | `@aagent EnsembleProblem`, AA read surface | reads finished state only |
| 246 | `src/export.jl` | 172 | **NEW (Phase-0.6, ADR 0013 §14.3):** `export_run`/`export_ensemble` directory bundles (CSV+JSON core); the `_arrow_write` empty-generic weakdep seam + `_arrow_available`. | reads `prob.acs` for the model hash | none | reads finished state only |
| 247 | `src/visualize.jl` | 285 | **NEW (Phase-0.6, ADR 0014 §15.2):** Layer A `network_graph`/`NetworkGraph`; Layer B `to_graphviz`/`draw_network`; Layer C `exec_map`. | walks `transLHS`/`transRHS` incidence | composes with AA's `run_graphviz`/`wiring_diagram` | Layer A runs `sample_transitions!` on a `deepcopy` (RNG-pure); else read-only |
| ext | `ext/RDPlotsExt.jl` | 175 | **NEW (Phase-0.6, ADR 0014 §15.1):** the seven Plots.jl `@recipe`s + the relocated `AlgebraicAgents._draw`. Loads only when `Plots` is present (weakdep). | reads raw artifacts | defines `_draw` | none (post-hoc viz) |
| ext | `ext/RDArrowExt.jl` | 21 | **NEW (Phase-0.6, ADR 0013 §14.3):** the SOLE `ReactiveDynamics._arrow_write` method (calls `Arrow.write`). Loads only when `Arrow` is present (weakdep). | none | none | none (I/O) |

> Line-number convention: the `Lines` column and all `file:line` citations in this document are 1-indexed last-line numbers (matching the editor/Read view). For files lacking a trailing newline this reads +1 vs `wc -l`. The WS-4 refresh (below) rows the previously-missing files and re-syncs the grown counts to current `wc -l` (so those rows now match `wc` rather than adding the +1); older untouched citations may still carry the +1 convention. Load-# values were re-synced to the current include order (`ReactiveDynamics.jl:222-247`).
>
> **Drift note (2026-06-30, RESOLVED by WS-4 2026-07-09).** This table was authored against an early `ref-agents` snapshot (pre-Phase-1/2). The 2026-06-30 note flagged five un-rowed files (`src/exprnode.jl`, `src/predicates.jl`, `src/actions.jl`, `src/ledger.jl`, `src/serialize.jl`) and several undercounts. WS-4 (housekeeping) has now ADDED rows for all five plus the two previously-missing interface files (`aa_coupling.jl`, `checkpoint.jl`), re-synced the grown counts (`state.jl` 366, `solvers.jl` 1067, `agents.jl` 169, `create.jl` 333, `reaction_parser.jl` 117 by `wc`), refreshed the Load-# column to the current include order, and rewrote the `solve.jl`/`update.jl`/`loadsave.jl` rows to reflect the WS-4 dangling-export cleanup and the ADR-0005 eval-free loader. The Phase-0.6 rows (`analysis.jl`/`export.jl`/`visualize.jl` + `ext/`) were already current. Line numbers in citations elsewhere in the doc that predate this refresh should still be re-verified before use.

---

## Public-API Audit

The module root declares **no `export`**; every public symbol originates in an included file. Symbols the root defines (`ReactionNetworkSchema`, `assign_defaults!`, `merge_acs!`, `add_obs!`, `TheoryReactionNetwork`, the type aliases) are public-by-use but UNEXPORTED. The RD-owned export surface has grown to **98 symbols** across the Phase-1/2 files (`predicates.jl`, `actions.jl`, `ledger.jl`, `serialize.jl`, `analysis.jl`, `export.jl`, `visualize.jl`, `checkpoint.jl`, `aa_coupling.jl` all add exports beyond the original snapshot's 41).

**WS-4 status (2026-07-09): every RD-owned export now resolves — the historical dangling exports are eliminated.** The exhaustive check is now an automated invariant, not a hand-maintained list: `test/semantic/exports_resolve.jl` iterates `names(ReactiveDynamics)` and asserts every RD-owned export `isdefined` and (for macros) resolves to a `Function`. That test SCOPES OUT symbols merely reexported from upstream (`@reexport using AlgebraicAgents`/`GeneratedExpressions`), because a few of those are dangling in the UPSTREAM package (e.g. AA's `@derived`, `@integration`, `AgentCall`, `add_integration_to_path`) — not RD's exports to police.

The curated table below maps the CORE public API + the historically-flagged exports (it is not the full 98-symbol enumeration — the test is the exhaustive guarantee). Definition-line citations predating the WS-4 refresh may have drifted; the `Resolves?` column is current. `@cost`/`@reward`/`@valuation` share one row (all defined via an include-time `eval` loop).

| Symbol | Kind | Export site | Definition | Resolves? |
|---|---|---|---|---|
| `ReactionNetworkProblem` | function/type | `solvers.jl:4` | struct `state.jl:40`; constructor `solvers.jl:536` | yes |
| `@ReactionNetworkSchema` | macro | `create.jl:3` | `create.jl:61` → `make_ReactionNetwork` `:78` | yes |
| `@append_transitions` | macro | `create.jl:4` | `create.jl:314` | yes |
| `AbstractStructuredToken` | type | `agents.jl:1` | `agents.jl:6` | yes |
| `BaseStructuredToken` | type | `agents.jl:1` | `agents.jl:11` | yes |
| `@structured_token` | macro | `agents.jl:2` | `agents.jl:30` | yes |
| `register_structured_species!` | function | `agents.jl:3` | `agents.jl:18` | yes |
| `add_structured_token!` | function | `agents.jl:3` | `agents.jl:42` | yes |
| `@push` | macro | `update.jl:3` | `update.jl:46` | yes |
| `@name_transition` | macro | `update.jl:3` | `update.jl:61` | yes |
| `@mode` | macro | `update.jl:3` | `update.jl:131` | yes |
| `@add_species` | macro | `update.jl:3` | `update.jl:213` | yes |
| `@periodic` | macro | `update.jl:4` | `update.jl:483` | yes |
| `@jump` | macro | `update.jl:4` | `update.jl:496` | yes |
| `@prob_init` | macro | `update.jl:5` | `update.jl:240` | yes |
| `@prob_uncertainty` | macro | `update.jl:5` | `update.jl:303` | yes |
| `@prob_params` | macro | `update.jl:5` | `update.jl:373` | yes |
| `@prob_meta` | macro | `update.jl:5` | `update.jl:411` | yes |
| ~~`@prob_role`~~ | macro | — | **REMOVED (WS-4)** — was a dangling `export` at `update.jl:6`; no definition, no `specRole` attr | n/a (export deleted) |
| ~~`@list_by_role`~~ | macro | — | **REMOVED (WS-4)** — dangling `export`, no definition | n/a (export deleted) |
| ~~`@list_roles`~~ | macro | — | **REMOVED (WS-4)** — dangling `export`, no definition | n/a (export deleted) |
| ~~`@prob_check_verbose`~~ | macro | — | **REMOVED (WS-4)** — macro + export deleted; its body called an undefined `check_params` (threw at call time). Vestigial: no code/test/demo/doc reference; CONTRACT specifies no such check | n/a (removed) |
| `@aka` | macro | `update.jl:6` | `update.jl:442` | yes |
| `@register` | macro | `update.jl:7` | `update.jl:516` | yes |
| `@cost` / `@reward` / `@valuation` | macro | `update.jl:168` (+ redundant `:172`) | generated by `eval` loop `update.jl:169-202` (macro body `:181`) | yes |
| `union_acs!` | function | `joins.jl:2` | `joins.jl:10` | yes |
| `@join` | macro | `joins.jl:2` | `joins.jl:201` | yes |
| `equalize!` | function | `equalize.jl:1` | `equalize.jl:24` | yes |
| `@equalize` | macro | `equalize.jl:1` | `equalize.jl:77` | yes |
| `@agentize` | macro | `solve.jl:13` | **implemented — thin constructor sugar (WS-B2)** — expands to exactly one `ReactionNetworkProblem(acs[, u0, p]; …)` call (no second construction path); auto-names from a bare-binding acs. Was REMOVED in WS-4 (dangling export); ADR 0001's "thin sugar OR delete" alternative, re-added per ADR 0012 | yes |
| `@import_model` | macro | `loadsave.jl:8` | `loadsave.jl:30` (delegates to `serialize.jl` typed-IR loader) | yes |
| `@export_model` | macro | `loadsave.jl:8` | `loadsave.jl:48` | yes |
| `@export_solution_as_table` | macro | `loadsave.jl:9` | `loadsave.jl:59` | yes |
| `@export_solution_as_csv` | macro | `loadsave.jl:9` | `loadsave.jl:70` | yes |

> **INVENTORY was STALE here (corrected WS-4).** The prior audit listed `@import_network`/`@export_network`/`@load_models`/`@import_solution`/`@export_solution` and flagged `@export`/`@import` as dangling at `loadsave.jl:5`. All of that predates the ADR-0005 pivot: the current `loadsave.jl` (72 lines) exports ONLY `@import_model`, `@export_model`, `@export_solution_as_table`, `@export_solution_as_csv` (all DEFINED). There is no `export @export, @import` line and no `@..._network`/`@..._solution` macro on `ref-agents` — the TOML/CSV/JLD2 model zoo (and its import-time RCE) was removed. So the historical `@export`/`@import` dangling-export finding is MOOT (nothing to delete).

### Historically exported-but-UNDEFINED — now all RESOLVED
The six danglers the prior audit flagged are eliminated as of WS-4 (2026-07-09):
- **`@agentize`** (`solve.jl:1`) — the WS-4 dangling `export` was DELETED, then WS-B2 RE-ADDED the macro as thin constructor sugar (expands to `ReactionNetworkProblem(acs[, u0, p]; …)`, auto-names from a bare binding). It is once again exported AND now defined; the exports-resolve guard pins it.
- **`@prob_role`, `@list_by_role`, `@list_roles`** (`update.jl:6`) — dangling `export` DELETED; no definitions, no `specRole` attribute, the "role" concept was dropped. ADR 0009's `PortRole` is an unrelated per-Species field authored in `@ReactionNetworkSchema`, NOT a repurpose of these config-macro exports (coordinated with WS-2).
- **`@export`, `@import`** (`loadsave.jl:5`) — MOOT: the export line no longer exists on `ref-agents` (see stale-note above); INVENTORY had drifted.
- **`@prob_check_verbose`** (`update.jl:469`) — the macro + its export were REMOVED (it called an undefined `check_params`; vestigial with zero references).

This is now enforced by `test/semantic/exports_resolve.jl` (see intro), so the class of bug cannot silently reappear.

---

## Schema / static store (post-WS-1)

> **Refreshed 2026-07-09 for the ADR-0003 store swap (WS-1).** ACSets is GONE (dropped from `[deps]`). The static authoring/IR store is now a dependency-free typed-struct-of-columns, defined in the module root `src/ReactiveDynamics.jl`. Line numbers below are approximate — re-grep before use.

- **`const SCHEMA`** — a `NamedTuple` of the six objects (`:S` species, `:T` transitions, `:E` events, `:obs` observables, `:P` params, `:M` meta — zero homs), each mapping `attr ⇒ element-type`. The single source of truth replacing the old ACSets `BasicSchema`; declaration order is load-bearing (it fixes `propertynames(acs.subparts)` order for the eight `"spec"`/`"trans"`-substring reflection loops). Derived: `ATTR2OBJ`, `ALLATTRS` (27 attrs), `columns(SCHEMA[,obj])`.
- **`struct ReactionNetworkSchema{parts::Dict{Symbol,Int}, subparts::NamedTuple, reactants::Vector{ReactantSpec}}`** — the static network container (unchanged NAME, so every `state.acs::ReactionNetworkSchema` annotation still compiles). `subparts` is a `NamedTuple` of `AttrColumn{T}(v::Vector{T}, def::Vector{Bool})` columns in `ALLATTRS` order; an undefined cell reads `nothing`. An explicit typed inner constructor prevents the auto-`(Any,Any,Any)` field ctor from colliding with the legacy 3-arg `ReactionNetworkSchema(transitions,reactants,obs)`. Content-based `Base.hash` (the derived `reactants` table is excluded so `_model_hash` is stable across promotion).
- **Shim verbs** (RD-owned, EXPORTED — no longer ACSets': `nparts`/`parts`/`dom_parts`/`incident`/`subpart`/`set_subpart!`/`add_part!`/`add_parts!`/`rem_parts!`). `getindex`/`setindex!` in 4 forms; `[:,attr]`/`subpart` return COPIES; `incident` = `findall(isequal)`. **`rem_parts!` is SWAP-AND-POP** (last row fills each freed slot, reverse-iterating sorted victims — cloning ACSets 0.2.29 exactly; equalize.jl relies on the surviving-row order).

Object/attribute model (element types):

| Object | Attributes (col → type) |
|---|---|
| `:S` species | `specName`→Symbol, `specModality`→`Set{Symbol}`, `specInitVal`/`specInitUncertainty`/`specCost`/`specReward`/`specValuation`→SampleableValues, `specStructured`→Bool, **`specRole`→Symbol (WS-2, ADR 0009 §A port role, default `:private`)** |
| `:T` transitions | `trans`/`transPriority`/`transRate`/`transCycleTime`/`transProbOfSuccess`/`transCapacity`/`transMaxLifeTime`/`transPreAction`/`transPostAction`/`transMultiplier`→SampleableValues, `transName`→`Union{String,Symbol,Missing}` |
| `:E` events | `eventTrigger`/`eventAction`→SampleableValues |
| `:obs` observables | `obsName`→Symbol, `obsOpts`→`FoldedObservable` |
| `:P` params | `prmName`→Symbol, `prmVal`→Any |
| `:M` meta | `metaKeyword`→Symbol, `metaVal`→SampleableValues |

**Promoted reactant relation (WS-1 Phase 2, ADR 0003):** `struct ReactantSpec{trans::Int (FK→:T), species::Int (FK→:S, 0=dynamic), stoich, side::Symbol, modality::Set{Symbol}, expr}` — a `Vector{ReactantSpec}` FIELD on the store (NOT a 7th SCHEMA object, so it stays out of the reflection loops). `populate_reactant_specs!` (serialize.jl) derives it from the `:trans` column via the eval-free `_split_reaction_line`/`_static_reactants` decomposition; the `expr` field is the escape-hatch for dynamic reactants. Accessors: `reactant_specs`/`specname`/`find_index(::Symbol,::ReactionNetworkSchema)`. `equalize!` re-derives it post-merge for a structural FK-repoint.

Supporting root tables/routines: `prettynames` (user alias → canonical attr), `defargs` (per-object defaults, incl. `:specStructured => false`, `:specRole => :private`), `species_modalities` (`[:nonblock,:conserved,:rate]`), `assign_defaults!`, `add_obs!`, `merge_acs!`, the 3/4-arg constructors. `compilable_attrs` was DELETED (dead). String→attribute coercion is `convert(Symbol,String)` on assignment (the only live hook; the old eval-based Set/FoldedObservable convert hooks were already removed per ADR 0005).

---

## Stepping Trace

Stepping is owned entirely by the native engine in `src/solvers.jl`, driven through the AlgebraicAgents interface. There is no separate standalone driver: both "standalone" and "AA-wrapped" usage call AA's reexported `simulate(prob, max_t=Inf)`, which loops `step!` while the projected time is continuable and below `max_t`. RD defines no `simulate` of its own.

### Driver chain (AA → engine)
> AA line numbers below are for **AlgebraicAgents v0.3.27** (the version resolved at audit time). The repo pins only `AlgebraicAgents = "0.3"` (compat, no `Manifest.toml`), so these cites are version-dependent — e.g. on v0.3.25 `simulate`/`step!` sit at `interface.jl:158`/`181`. RD's own `file:line` cites are exact.
1. `simulate(prob[, n])` (AA `interface.jl:160-167`) loops `step!` while `iscontinuable(ret) && ret < max_t` (`iscontinuable` = not-`nothing` and not-`Bool`). Tutorials call `simulate(prob)` / `simulate(prob, N)`.
2. `step!` (AA `interface.jl:183-206`) descends into INNER agents first (`foreach inners step!`), then runs the LOCAL `_step!(state)` only when `_projected_to(state) == t` (the least projected time across the hierarchy). So co-stepped siblings (a user `Controller`, structured tokens) step before the network body for a given tick. Structured tokens are clock-inert (`_step!`/`_projected_to` no-ops, `agents.jl:49-50`).
3. **Entry point: `AlgebraicAgents._step!(state::ReactionNetworkProblem)`** at `solvers.jl:642`.

### One tick — `_step!` body (`solvers.jl:642-673`), in order
1. `update_u_structured!` (`:643`) — sync structured-species counts from the agent layer.
2. `save!(state)` only if `sol` empty (`:644-646`) — record the initial `(t, u)` row.
3. `free_blocked_species!` (`:648`) — release `:nonblock` resources. **BROKEN**: undefined `q` at `solvers.jl:512` → `UndefVarError` on any ongoing `:nonblock` LHS token (so the `:nonblock` release path is effectively dead).
4. `update_u_structured!` (`:649`).
5. `update_observables` (`:650` → `state.jl:144`) — resample observables with `(t - last) >= every`.
6. `sample_transitions!` (`:651` → `state.jl:174`) — clear `state.transitions`, eval activated recipes, unfold LHS into `UnfoldedReactant`. (Resets `transToSpawn .= 0` INSIDE the per-recipe loop, `state.jl:215`.)
7. `evolve!` (`:652` → `solvers.jl:133`) — SPAWN + ALLOCATE + advance ongoing (detail below).
8. `update_u_structured!` (`:653`).
9. `finish!` (`:654` → `solvers.jl:400`) — TERMINATE + PoS + emit RHS (detail below).
10. `update_u_structured!` (`:655`).
11. `event_action!` (`:657`) — events run once/tick but `eventAction` is fetched and never evaluated (`solvers.jl:323`), so events are a no-op.
12. push `:valuation` ledger row (`:659-666`).
13. **Advance clock**: `state.t += state.dt` (`:668`) — the SINGLE clock advance, after all tick work.
14. `save!(state)` (`:670`) — records the post-increment `(t, u)` row.

Termination/progress: `_projected_to` (`solvers.jl:675`) returns `true` once `state.t > state.tspan[2]`, else `state.t`.

### Spawning (Poisson) — `evolve!` (`solvers.jl:133-229`)
- The rate is authored at parse time as `rand(Poisson(max(state.dt*rate, 0)))` (`create.jl:151`, `expand_rate`); `@ct`/`@cycletime` rewrite to `1/arg`.
- Per transition: `qs[i] = transRate * transMultiplier`, `ceil` to Int (`:140-144`); apply `transCapacity` (overflow deferred via `add_to_spawn!` at `:152` — but `add_to_spawn!` is BROKEN, see below); compute init requirements (`get_reqs_init!`, `:156`, skips `:rate`/`:nonblock` tokens); allocate (`get_allocs!`, `:158`); floor to whole instances (`get_init_satisfied`, `:160`); construct `Transition` heap entries (`:178-189`) with `t=state.t`, `q=qs[i]`, `state=0.0`; bind structured tokens by descending `priority` (`:202-223`); run `transPreAction` (`:227`). Ongoing transitions then advance: `transition.state += qs[i]*state.dt` (`:261`).

### Completion (cycle time + PoS) — `finish!` (`solvers.jl:400-508`)
- Terminate when lifetime exceeded `(t - trans.t) >= transMaxLifeTime` OR cycle complete `trans.state >= transCycleTime` (`:408-410`); draw success count `q = rand(Binomial(Int(trans.q), transProbOfSuccess))` (`:412-416`); emit RHS products into `state.u` (`:418-435`); return `:conserved`/`:nonblock` resources (`:437-483`); run `transPostAction` (`:485`); prune (`:501`). NOTE: `filter!` at `:501` keeps `s.state < s[:transCycleTime]`, so transitions terminated solely by MAX-LIFETIME (state still < cycleTime) are NOT removed and get re-evaluated/re-emitted every tick — a latent bug.

### Allocation under contention — `get_allocs!` (`solvers.jl:53`)
Dispatch on `state.p[:strategy]` (default `:weighted`, set at `:550`/`:599`): `:weighted` → `alloc_weighted!` (`:61`, demand×priority proportional to supply); anything else → `alloc_greedy!` (`:74`, descending-priority fill until supply exhausted).

### Clock / dt / tspan notes
- Single authoritative clock is `state.t`, advanced only at `solvers.jl:668`.
- **dt/tstep naming split**: runtime field is `dt` (`state.jl:53`, used at `:668`/`:37`/`:261`), but the constructor fills it positionally from `get(keywords,:tstep,1)` (`solvers.jl:603`) while `get_tcontrol` returns the derived value into `keywords[:tstep]` (`:552`). Reconciled by position, not name.
- **tspan** is not hardcoded (from `@prob_meta tspan=...` or the `tspan` kwarg; kwarg overrides meta), but `keywords[:tspan]` is indexed with no default at `solvers.jl:552` → `KeyError` if neither is supplied.
- **`add_to_spawn!`** (`state.jl:251-256`, called `solvers.jl:152`) is doubly broken: `findfirst` is given a scalar `length(...)` Int instead of a range, and on match it increments `:transHash` (`+= n`) rather than `:transToSpawn` — capacity-overflow deferral is non-functional.
- **Intra-tick token snapshot staleness**: `evolve!` snapshots the structured-token set once (`solvers.jl:173`) and reuses it for both the spawn and ongoing-advance loops, so tokens entangled mid-tick are not visible to later allocation within the same tick.

---

## AlgebraicAgents Touchpoints

AA is integrated at three layers; all references verified against current source.

### Import / reexport
- `@reexport using AlgebraicAgents` at **`state.jl:1`** is the single AA entry point, reexporting the full AA public surface to RD users (`@aagent`, `entangle!`, `getagent`, `inners`, `FreeAgent`, `AbstractAlgebraicAgent`, and the driver/viz wrappers `simulate`, `step!`, `projected_to`, `draw`).
- A redundant `import AlgebraicAgents` at `agents.jl:46` makes the qualified `AlgebraicAgents._step!/._projected_to` extensions and `AlgebraicAgents.aagent(...)` call explicit. Harmless.

### `@aagent` structs (4) + abstract type
| Type | Location | Notes |
|---|---|---|
| `Transition` | `state.jl:14` | in-flight instance; holds `Vector{AbstractAlgebraicAgent}` resource fields |
| `Observable` | `state.jl:31` | per-observable sampler agent |
| `ReactionNetworkProblem` | `state.jl:40` | the top-level agent = the engine state |
| `BaseStructuredToken` (`@aagent FreeAgent`) | `agents.jl:11` | base structured-token struct |
| `AbstractStructuredToken <: AbstractAlgebraicAgent` | `agents.jl:6` | root of the token hierarchy (abstract type) |

`@aagent` injects `(uuid, name, parent, inners, relpathrefs, opera)` and a `(name::AbstractString, args...)` constructor; this is why every RD agent is constructed name-first (`Transition(...)` at `solvers.jl:178`, `Observable(...)` at `state.jl:112`, `ReactionNetworkProblem(...)` at `solvers.jl:593`). No `@aagent_call`/`@call` is used anywhere; user token subtypes are generated via the programmatic `AlgebraicAgents.aagent(...)` (`agents.jl:32`).

### AA interface methods implemented
| Method | Type | Location |
|---|---|---|
| `_step!` | `ReactionNetworkProblem` | `solvers.jl:642` |
| `_reinit!` | `ReactionNetworkProblem` | `solvers.jl:619` |
| `_projected_to` | `ReactionNetworkProblem` | `solvers.jl:675` |
| `_draw` | `ReactionNetworkProblem` | `plots.jl:13` (reached via exported `draw`) |
| `_step!` (no-op) | `AbstractStructuredToken` | `agents.jl:50` |
| `_projected_to` (no-op) | `AbstractStructuredToken` | `agents.jl:49` |

NOT implemented anywhere in RD: `getobservable`, `gettimeobservable`, `getparameters`, `setparameters!`, `_getparameters`, `_setparameters!`. There is NO AA wire/bidirectional-observable mechanism — RD's `observables::Dict{Symbol,Observable}` + `compile_observables`/`update_observables` (`state.jl:94`/`:144`) are an internal sampler, unrelated to AA's observable/wire API.

### Agentization & token hierarchy
- **`@agentize` is exported (`solve.jl:1`) but undefined.** Agentization is implicit: the `ReactionNetworkProblem` outer constructor (`solvers.jl:536`) builds the `@aagent` object and entangles a `FreeAgent("structured")` child (`solvers.jl:612`). The real public contract is `prob = ReactionNetworkProblem(acs[, u0, p]; name, kwargs...)` then `simulate(prob[, n])`.
- Structured tokens are first-class AA agents under the `"structured"` `FreeAgent` container; attached via `entangle!(getagent(problem,"structured"), agent)` (`add_structured_token!`, `agents.jl:42-44`); the solver reaches them via `inners(getagent(state,"structured"))` (`solvers.jl:173,631`) and entangles new/moved tokens during stepping (`solvers.jl:342,356,376,393`).
- No dedicated network back-reference on the state; the network↔token link is the AA `parent`/`inners` hierarchy. A token's `bound_transition` (`agents.jl:13`) points to its in-flight `Transition`, not the network.

> KNOWN BUGS in the token-binding paths (introduced by the pivot): `set_bound_transition!` is called with a `.bound_transition` field (a `Transition`/`Nothing`) instead of the token agent at `solvers.jl:288,452` (the method at `agents.jl:62` expects an `AbstractStructuredToken`; correct calls at `:215,489`); `delete!` is called on a `Vector` with an Int index at `solvers.jl:455` (should be `deleteat!`).

---

## What Changed Since REVIEW.md

REVIEW.md describes the OLD design (SciML `DiscreteProblem` front-end, `optim.jl` fitting, Catlab schema). The `ref-agents` branch is a full pivot. All `file:line` below are current unless prefixed `main`.

### Native engine replaces DiscreteProblem
- On `main`: `export DiscreteProblem`, `using DiffEqBase, DifferentialEquations`, a `DiffEqBase.DiscreteProblem` transform (`main solvers.jl:3,5,310,371`), and an `EnsembleProblem` solve path.
- On `ref-agents`: `solvers.jl` imports only `Distributions`/`Random` and `export ReactionNetworkProblem`. `ReactionNetworkProblem` is an `@aagent` (`state.jl:40`) with constructor `solvers.jl:536` and the AA stepping contract (`_step!:642`, `_reinit!:619`, `_projected_to:675`).
- `Project.toml` dropped `DiffEqBase`/`DifferentialEquations`/`OrdinaryDiffEq`/`NLopt`/`Catlab` and added `AlgebraicAgents`/`ACSets`.

### Schema/macro renames + ACSets migration
- `const ReactionNetwork` → `const ReactionNetworkSchema` (`main ReactiveDynamics.jl:78` → `:80`).
- `@ReactionNetwork` → `@ReactionNetworkSchema` (`main create.jl:60` → current `create.jl:61`).
- `@present TheoryReactionNetwork(FreeSchema)` → `BasicSchema(...)` (`main ReactiveDynamics.jl:31` → `:30`).

### New structured-token layer
- `interface/agents.jl` (added commit `8652f3d`): `AbstractStructuredToken`/`BaseStructuredToken`, `@structured_token`, `register_structured_species!`/`add_structured_token!`, binding/priority accessors.
- New `specStructured` Bool schema attr (`ReactiveDynamics.jl:51`, default `false` at `:135`). Tokens entangled under `FreeAgent("structured")` (`solvers.jl:612`); `_step!`/`_projected_to` no-ops.

### Removals + new viz
- `optim.jl` REMOVED (deleted commit `4e93e19`; had NO exports — all symbols module-internal). Dropped `build_parametrized_solver`, `optim!`, `build_loss_objective[_datapoints]`, `prep_params!`/`prep_u0!`, `get_free_vars`/`get_vars`. Include left commented at `ReactiveDynamics.jl:225`. The `@optimize`/`@fit`/`@fit_and_plot` macros survive only in the stale `src/macros_overview/macros.mmd` diagram.
- The SciML front-end macros are gone: `main solve.jl` exported `@problematize/@solve/@optimize/@fit/@fit_and_plot/@build_solver`; `ref-agents` keeps only `@plot` + helpers and the dead `export @agentize`.
- New `plots.jl` (added `4e93e19`) supplies the AA `_draw` overload (`plots.jl:13`).

### State struct reshape
- `mutable struct ReactiveDynamicsState` with `solverargs::Any` (`main state.jl:32,47`) → `@aagent struct ReactionNetworkProblem` with explicit `tspan::Tuple`/`dt::Float64`/`sol::DataFrame` (`state.jl:40,52,53,62`); `solverargs[:tstep]` replaced by `state.dt`.

### Brief corrections (claims that DO NOT hold)
- **`@agentize` is exported but NEVER defined** (`solve.jl:1`) — a dangling export, not a working macro. The contract is `ReactionNetworkProblem(...) |> simulate`.
- **`transMultiplier` was NOT added on this branch.** It existed on `main` (`main ReactiveDynamics.jl:57`, default `:transMultiplier => 1` at `:119`); on `ref-agents` it is the same column re-expressed as a `BasicSchema` tuple (`:62`, default `:124`). What changed nearby is the spawn computation (`ceil` of `qs` at `solvers.jl:144`; dt-scaled Poisson moved into `expand_rate`, `create.jl:151`).

### Stale-against-native cautions
- `solve.jl`'s `@plot` is dead-on-arrival: reads `sol.prob.p[:__state__]` (`:88,91,94`) and references unimported SciML `EnsembleSummary`/`EnsembleSolution` (`:9,38,45`). Native viz flows through `plots.jl` `_draw` (`draw(prob)`).
- **Log record-type drift (narrower than first thought).** The engine DOES emit `:new_transitions` (`solvers.jl:165`), `:saturation` (`:247`), `:allocation` (`:304`), `:valuation_cost` (`:308`), `:terminated_all`/`:terminated_success`/`:valuation_reward` (`:503-505`), and `:valuation` (`:662`). `complete_log_valuations` (`solve.jl:191`) reads `:valuation`/`:valuation_cost`/`:valuation_reward` — all emitted. The ONLY genuine orphans are `solve.jl:180`'s `:terminated_transitions` selector (engine emits `:terminated_all`/`:terminated_success` instead) and `solve.jl:189`'s `:valuations` (a plot-type selector, not a log key). `complete_log_valuations` also allocates 4 columns but never fills column 4 ("total balance"), so that series is always zeros.
- The `log` field is declared `Vector{Tuple}` (`state.jl:57`) but allocated `NamedTuple[]` (`solvers.jl:579`) while all pushes are plain Tuples — a contract inconsistency to pin down (NamedTuple is not `<: Tuple`).

---

## REVIEW.md Defect Reconciliation

REVIEW.md's 30-item catalog reconciled against current `ref-agents` source. Tally: **3 FIXED/MOOT-as-fixed (#1, #19), 4 MOOT (#8, #22; plus the optim-coupled parts of #19/#20), ~5 PARTIALLY-applicable, and ~18 STILL-APPLICABLE** (several now worse, several with new sibling bugs). All REVIEW.md line numbers are stale; current lines cited.

| # | Topic | Verdict | Current evidence |
|---|---|---|---|
| 1 | join mid-loop return / printlns | **FIXED** | `joins.jl:14-28` clean loop, `:58` returns acs1. Residual: obs (`:obs`) merge still absent (`prepend_obs` dead `:87-97`), `:E` not merged |
| 2 | eval-on-import RCE | **STILL APPLIES** | `loadsave.jl:65,72`; `ReactiveDynamics.jl:101-102`; export drops `registered` section; `@export`/`@import` still undefined exports (`loadsave.jl:5`) |
| 3 | per-step eval of fresh closure | **STILL APPLIES (worse)** | `compilers.jl:122-126` eval'd per-step via `state.wrap_fun` at `solvers.jl:227,340,420,430,485`, `state.jl:201` |
| 4 | bare `specModality` (not `:specModality`) | **STILL APPLIES** | `update.jl:108` (bare) vs `:109` (`:specModality`); masked by `assign_defaults!` seeding |
| 5 | undefined `q` in `free_blocked_species!` | **STILL APPLIES** | `solvers.jl:512`; called every step `:648` |
| 6 | `resample!` assigns nonexistent `o.val` | **STILL APPLIES** | `state.jl:137`; field is `sampled` (`:37`) |
| 7 | `add_to_spawn!` findfirst/`transHash`; reset-in-loop | **STILL APPLIES** | `state.jl:251-256`, `:215`; caller `solvers.jl:152` |
| 8 | optim warning-path crash | **MOOT** | optim.jl removed (`ReactiveDynamics.jl:225`); no `prep_params!`/`prep_u0!` |
| 9 | silently-dropped attr names | **PARTIALLY** | mechanism at `create.jl:179-182`; `transMultiplier` has no alias (`ReactiveDynamics.jl:104-115`); exact tutorial offenders not located |
| 10 | eval-coupled, untyped construction | **STILL APPLIES** | `ReactiveDynamics.jl:10`; `create.jl:37,61`; `compilers.jl:67`; no `validate(acs)` |
| 11 | no RNG seeding / reproducibility | **STILL APPLIES** | global rand: `solvers.jl:141,321,413`, `state.jl:123`; no AbstractRNG/seed threaded |
| 12 | type-unstable state / deepcopy | **PARTIALLY** (SciML/optim deepcopies moot) | `state.jl:43-47,55,80,219-221`; `p::Any`, Dict caches, string-routing getindex |
| 13 | non-orthogonal, unvalidated modalities | **STILL APPLIES** | `ReactiveDynamics.jl:144,45`; inconsistent solver branching `solvers.jl:20,35-43,442`; new `:conserved+:nonblock` guard `:461-465` |
| 14 | ledger "total balance" never computed; unit mix | **STILL APPLIES** | `solve.jl:146-163,190` (col 4 unfilled); flow vs stock mix `solvers.jl:308-311,505,659-666` |
| 15 | eval in module/global scope | **STILL APPLIES** | `update.jl:514,413,360,163` |
| 16 | assertion-free smoke tests | **STILL APPLIES** | `test/` is tutorial-include smoke tests with no assertions (per `safeinclude.jl`); underlying defect real |
| 17 | exported-but-undefined macros + phantom docstrings | **STILL APPLIES (expanded)** | `loadsave.jl:5`; `update.jl:6,470`; PLUS new `@agentize` (`solve.jl:1`); stale docstrings `create.jl:40,58` |
| 18 | missing-vs-nothing desync; `metaVal` latent crash | **PARTIALLY** | `ReactiveDynamics.jl:137-138,74,149`; guarded by `isnothing` at `solvers.jl:545-546,573` → crash latent/conditional, not guaranteed |
| 19 | hardcoded (0.0,2.0) + remake; dual clock | **MOOT / FIXED** | single clock `state.t` (`solvers.jl:668`), `get_tcontrol :526-534`, `_projected_to :675`; no remake/DiscreteProblem |
| 20 | kwarg leak / hygiene | **PARTIALLY** (inner ctors moot) | `utils.jl:18-31`; `solve.jl:89,92,99`; blast radius subsumed by #21 |
| 21 | `plot_from_log` no-else + STALE SciML | **STILL APPLIES (worse)** | `solve.jl:179-201` no else; `:9,38,45,88,91,94,96` reference removed SciML; whole `@plot` path dead-on-arrival |
| 22 | two near-duplicate solver builders; biased u0 sampling | **MOOT** | single ctor `solvers.jl:536`; no sign-clamp/`sample_u0!` |
| 23 | `compilable_attrs` broken/unused; compile-by-substring | **STILL APPLIES** | `ReactiveDynamics.jl:141-142` (dead); `compilers.jl:140-143,158`; `state.jl:80` |
| 24 | macro boilerplate dup; dead helpers/types | **STILL APPLIES** | `update.jl:136-148,182-198,246-261,309-323,376-391`; dead `create.jl:15,18,20,26,37,99,200-217,222-230` |
| 25 | string-matched pseudo-macros, no allowlist | **STILL APPLIES** | `create.jl:150,157,292-304`; `reaction_parser.jl` `@error`-continue |
| 26 | `@join` → undefined `include_model`; equalize precedence | **STILL APPLIES** | `joins.jl:226,228` (undefined); `equalize.jl:28` implicit first-term alias |
| 27 | non-const mutable global tables; redundant defaults | **STILL APPLIES** | `ReactiveDynamics.jl:104,117,141,144,157-161`; `compilers.jl:67`; `update.jl:450` |
| 28 | heavy hard deps; core type/builders unexported | **PARTIALLY** (type now exported) | `ReactionNetworkProblem` exported (`solvers.jl:4`); Plots/Distributions/CSV/JLD2 etc. still hard `[deps]`; `merge_acs!`/`add_obs!`/function `import_network`/`export_network` unexported |
| 29 | documentation drift | **PARTIALLY** (@per_step moot in src) | `macros.mmd:21+` still shows `@problematize`→DiscreteProblem/`@ReactionNetwork`; `create.jl:40,58` stale |
| 30 | assorted latent | **MIXED** | `Int(trans_.q)` Binomial InexactError `solvers.jl:413` (applies); `prune_r_line` bare `state` `state.jl:165` (applies); `convert` try/catch→`string(ex)` `ReactiveDynamics.jl:92-97` (applies); `history_u` & `optim!` (moot) |

### NEW defects introduced by the pivot (not in REVIEW.md's 30)
- `solvers.jl:288,452` — `set_bound_transition!` passed a `.bound_transition` (`Transition`/`Nothing`) instead of the token agent (method expects `AbstractStructuredToken`, `agents.jl:62`; correct calls `:215,489`).
- `solvers.jl:455` — `delete!` on a `Vector` with Int index (should be `deleteat!`, cf. `:383,479`).
- `solvers.jl:40` — `get_reqs_ongoing!` error string interpolates undefined `trans`.
- `solvers.jl:323` — `event_action!` fetches `:eventAction` but never evaluates it; events are a no-op.
- `solvers.jl:501` — lifetime-terminated transitions (state < cycleTime) are never pruned and re-emit RHS every subsequent tick.
