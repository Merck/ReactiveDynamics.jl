# ReactiveDynamics.jl — Ground-Truth Inventory (branch `rework`)

## Orientation

ReactiveDynamics.jl (RD) is NOT a chemical reaction network despite its DSL surface and Catalyst-derived parser. It is a timed, stochastic, resource-constrained Petri net / discrete-event system: transitions are stateful recipes that spawn instances at a Poisson (or deterministic) rate, occupy resources (**places**, in Petri-net terms) over a cycle time, and complete with a terminal Binomial probability-of-success that emits RHS products; resources carry a modality (`:nonblock`/`:conserved`/`:rate`) governing allocation, and a cost/reward/valuation ledger accumulates into a per-step `log` and a per-program ledger. The simulation output is the NATIVE engine type `ReactionNetworkProblem`, an AlgebraicAgents (AA) `@aagent` stepped by AA's reexported `simulate`/`_step!` — NOT a SciML `DiscreteProblem`. The static authoring/IR store is a dependency-free typed struct-of-columns: ACSets and Catlab were dropped entirely (ADR 0003) and the ACSets-lineage vocabulary was retired (ADR 0015), while ADR 0017 retired the chemical-reaction-network vocabulary (a resource pool is a **place**, its quantity a **marking**, a transition↔place participation an **arc**), so the store type is **`ReactionNetwork`** (its type-level object model is `const SCHEMA`), the authoring macro is **`@reaction_network`**, the conventional variable is **`net`**, the `ReactionNetworkProblem` field is **`.network`**, the store fields are **`.counts`**/**`.columns`**, and the store verbs are `nrows`/`row_ids`/`column`/`cell`/`find_rows`/`add_row!`/… — all UNEXPORTED (internal, `RD.`-qualified). A structured-token agent layer (`interface/agents.jl`), a typed eval-free `ExprNode` IR (`exprnode.jl`), single-JSON serialization (`serialize.jl`), token filtration (`predicates.jl`), the endogenous decision channel (`actions.jl`), the per-program ledger (`ledger.jl`), refinement/composition (`operators/refine.jl`), and the analysis/visualization layer (`analysis.jl`/`export.jl`/`visualize.jl` + two extensions) were all added. The old `optim.jl` was removed and the `GeneratedExpressions` dependency was dropped.

This document maps the CURRENT source on branch `rework`. Everything specified in the contract (§1–§15) and ADRs 0001–0015 is implemented and green — the semantic suite runs **801 pass / 0 broken** (see `spec/STATUS.md` for the feature→commit table). Re-verify any `file:line` before acting — line numbers drift; grep to confirm.

Section index: [Module Map](#module-map) · [Public-API Audit](#public-api-audit) · [Static Store (`ReactionNetwork` / `SCHEMA`)](#static-store-reactionnetwork--schema) · [Stepping Trace](#stepping-trace) · [AlgebraicAgents Touchpoints](#algebraicagents-touchpoints)

---

## Module Map

The module root `src/ReactiveDynamics.jl` declares the package deps (`Reexport`, `MacroTools`, `ComponentArrays`; AA and the rest are pulled in by included files), the core type aliases (`SampleableValues`/`ActionableValues`/`SampleableRange`), the `const SCHEMA` object model + the `ReactionNetwork` store type and its unexported store verbs, and the construction/merge helpers. Its ONLY `export` is `ArcSpec, arcs, placename` (`:15`); every other public symbol originates in an included file.

Include order (`src/ReactiveDynamics.jl:443-468`): `state.jl`, `exprnode.jl`, `compilers.jl` are explicit, THEN `interface/`, `utils/`, `operators/` are bulk-included via `include.(readdir(...; join=true))` (filesystem-alphabetical), then `solvers.jl`, `predicates.jl`, `actions.jl`, `ledger.jl`, `serialize.jl`, `#include("optim.jl")` (commented — file removed), `loadsave.jl`, and finally the analysis/viz trio `analysis.jl`/`export.jl`/`visualize.jl` (placed last so they see the run-state, ledger, predicate/sortkey, and serializer). The `readdir` bulk-includes carry a real-but-benign load-order dependency: `interface/agents.jl` (first alphabetically) subtypes `AbstractAlgebraicAgent` and names `ReactiveDynamics.Transition`, both of which exist only because `state.jl` was included first. Adding an alphabetically-earlier `interface/`/`operators/` file with top-level code that calls a later sibling would break the build.

| File | Lines | Purpose | Store | AA | Stepping |
|---|---|---|---|---|---|
| `src/ReactiveDynamics.jl` | 470 | Module root: deps, the type aliases, `const SCHEMA` (+ `ATTR2OBJ`/`ALLATTRS`), the `AttrColumn{T}`/`ArcSpec`/`ReactionNetwork` types, the unexported store verbs (`nrows`/`row_ids`/`column`/`cell`/`set_cell!`/`add_row!`/`add_rows!`/`rem_rows!`/`find_rows`), `PORT_ROLES`/`port_role`, `assign_defaults!`/`add_obs!`/`merge_networks!`, the deprecation shims. Orchestrates includes (`:443-468`). | Store home (defines it) | none | none (static construction) |
| `src/state.jl` | 371 | Live simulation state and sub-structs: `ReactionNetworkProblem` (top-level `@aagent`, ~30 fields incl. `network`, `rng`/`seed`/`initial_rng`, `rules`, `registry`, `population`, `program_ledgers`, `external_inputs`, `token_trajectory`), in-flight `Transition`, `Observable` sampler, `UnfoldedArc`, `ProgramLedger`. `getindex` pass-through to the store, `context_eval`, observable resampling, `sample_transitions!`, `add_to_spawn!`. | heavy reader via `state.network` | `@reexport using AlgebraicAgents` (`:1`); 3 `@aagent` structs | per-step machinery (`update_observables`, `sample_transitions!`, `save!`) |
| `src/exprnode.jl` | 181 | The typed, eval-free expression IR (ADR 0005). A closed tagged-union `ExprNode` tree (`Const`/`NodeRef`/`Call`/`Sample`/`TimeRef`/`Choose`/`Field`/`ExternalRef`) is the authoring/serialization form of every time-varying attribute; `to_expr` lowers a tree to EXACTLY the `Expr` the DSL produces (then compiled ONCE by `compilers.jl`), `from_expr` is the structural inverse. The closed `OP_WHITELIST`/`DIST_WHITELIST`/`REF_KINDS` are the only Julia ever produced from an inert model document. | none (operates on `Expr`) | none | none (authoring/serialization boundary) |
| `src/compilers.jl` | 205 | Expression-compilation layer: rewrites DSL attribute exprs (rates/stoich/actions/observables) into runnable `(state, transition)` closures; place→`u[i]`, params→`p[:name]`, query markers→state calls; bulk `compile_attrs`. | reads store attrs | none | setup (builds `transition_recipes`); `wrap_expr` also fires per-step via `state.wrap_fun` |
| `src/interface/agents.jl` | 173 | Structured-token agent layer (ADR 0006): `AbstractStructuredToken`/`BaseStructuredToken`, `@structured_token`, `register_token_kind!`/`add_structured_token!`, the `PopulationEntry` + `instantiate_population!` declarative initial marking (ADR 0007 §B), the `log_token_fields` hook (ADR 0013 §14.1), and the unexported binding/priority/place accessor protocol. Tokens are clock-inert (`_step!`/`_projected_to` no-ops, `:138-139`). | `register_token_kind!` mutates `:S`/`:placeStructured` | `@aagent FreeAgent`; no-op step | supplies the protocol the loop calls |
| `src/interface/create.jl` | 365 | CREATE half of the DSL: `@reaction_network` (`:66-78`) + the `get_data` parser pipeline (arrow normalization, rate→Poisson via `expand_rate`, keyword-attr extraction, event expansion, `@register`/`@take`/`@select`/`@advance`/`@structured` rewriting); `@append_transitions`; the deprecated `@ReactionNetworkSchema` alias (`:83`, depwarns → `@reaction_network`). | emits the `ReactionNetwork(get_data(...)...)` build | none | parse-time; `expand_rate` AUTHORS the spawn-rate Poisson expr |
| `src/interface/parsing_utils.jl` | 42 | Low-level AST helpers (from Catalyst): tuple arity/index, multiplication flattening (`multiplex`), dot expansion. | none | none | none (parse-time) |
| `src/interface/plots.jl` | 11 | Comment-only stub — `plot_df` deleted and `_draw` RELOCATED to `ext/RDPlotsExt.jl` (Plots demoted to a weakdep, ADR 0014). | none | — | none |
| `src/interface/reaction_parser.jl` | 117 | Arc extraction (from Catalyst): folds a reaction side into `FoldedArc`s carrying place/stoich/modality; parses the `@select(Kind, clauses)` LHS into a `TokenPredicate` (`:45-69`, `PRED_OP_WHITELIST`); resolves `@choose`, escapes refs. | reads `net[:, :placeName]` | indirect (consumes state `getindex`) | sampling/parse-support |
| `src/interface/aa_coupling.jl` | 138 | The AA READ/COUPLING surface (ADR 0012 / CONTRACT §13): OUTBOUND `observables`/`getobservable` (RD as a wire source) + `_getparameters`/`_setparameters!` (param read/patch); INBOUND `_prestep!` latch feeding `ExternalRef` leaves from incoming AA wires (explicit one-tick-lag Jacobi coupling). | reads run-state observables | heavy: implements AA observable/wire/param hooks | read surface; coupling at the AA time gate |
| `src/interface/checkpoint.jl` | 118 | State dump/restore (ADR 0007 §C / CONTRACT §10.5): `StateDump`/`dump_state`/`restore` serialize a live run at a TICK BOUNDARY into an eval-free JSON-representable artifact (initial marking + dynamic state: clock, RNG, creation counters, `u`, token population, rule latches). `dump_state` refuses mid-cycle (errors unless `ongoing_transitions` empty) — a scoped Milestone-1 constraint. | reads/reconstructs the store + run-state | dumps/restores the AA token population | boundary-only (not in the hot loop) |
| `src/interface/solve.jl` | 64 | `@agentize` thin constructor sugar (`:45`): expands to exactly one `ReactionNetworkProblem(net[, u0, p]; …)` call (no second construction path), auto-naming from a bare-binding net. The dead SciML `@plot`/`plot_summary`/`plot_ensemble_sol`/`first_sol`/`plot_from_log` were deleted (ADR 0014 §A1). `export @agentize` co-located at `:13`. | none | expands to the constructor call | none |
| `src/interface/update.jl` | 524 | UPDATE half of the DSL: mutation/config macros (`@push`, `@add_place`, `@mode`, `@name_transition`, `@prob_init`/`@prob_uncertainty`/`@prob_params`/`@prob_meta`, `@cost`/`@reward`/`@valuation`, `@periodic`/`@jump`, `@aka`, `@register`) that reshape syntax into store mutations. | heavy reader/mutator | none | seeds step inputs (`@jump`/`@periodic`/`@prob_meta`) |
| `src/utils/utils.jl` | 104 | AST/kwarg grab-bag (`macroname`, `blockize`, `underscorize`, `args_kwargs`, …). NO exports. | none | none | none |
| `src/operators/equalize.jl` | 94 | Place identification/collapse (authoring-time): `equalize!` merges matched place rows + rewrites refs and re-derives the `ArcSpec` table; `@equalize` parses `A = B = C` blocks. | heavy mutator | none | none (model transform) |
| `src/operators/joins.jl` | 282 | Model composition (authoring-time): `merge_networks!` merges one network into another in place (S/P/M/E/obs by name, T appended); `@join` builds + folds models. | heavy reader/mutator | none | none (model transform) |
| `src/operators/refine.jl` | 423 | Hierarchical refinement & open-port composition (ADR 0009 / CONTRACT §11): `refine`/`refine!`/`abstract`/`abstract_transitions` (boundary-matched FK-splice, substitutable granularity), `set_port_role!`/`@port` (place port roles), `@compose`/`compose`, `@pipeline`/`@process` compact authoring, `refinement_diagnostics`. | reads/rewrites the store + `ArcSpec` FKs | none | none (authoring-time transform) |
| `src/solvers.jl` | 1230 | CORE discrete-event engine: the `ReactionNetworkProblem` constructor (`:912`) + the per-tick step loop `_step!` (`:1170`), the ADR-0002 weighted-progressive-fill allocator (`build_requirements!` `:100`, `progressive_fill!` `:162`, `spawn_integer!` `:268`), `evolve!` (`:305`) / `finish!` (`:668`) / `free_blocked_places!` (`:809`), structured-token RHS (`structured_rhs` `:561`), construction-time modality validation (`validate_modalities` `:856`), `get_tcontrol` (`:826`), and the AA contract `_reinit!`/`reinit!`/`_projected_to` (`:1080`/`:1150`/`:1219`). | heavy reader/writer | implements the AA stepping contract | THIS IS THE STEP LOOP |
| `src/predicates.jl` | 108 | Token filtration (ADR 0008 / CONTRACT §9.5): `TokenPredicate{kind, clauses}` (+ `Clause`, `matches`) narrows the candidate set a structured LHS binds; `select_tokens` (`:71`) filters then sorts by `token_sortkey` (`:105`, the `(place, creation_index)` total order, §4 D4). Phase is an ATTRIBUTE (`@select` selects, `@advance` advances) — no per-phase kinds. | reads token/place attrs | operates over structured-token agents | binding-time filter in the step loop |
| `src/actions.jl` | 253 | The endogenous decision channel (ADR 0010/0011 / CONTRACT §12): the closed `ActionStmt` family `{SetMarking, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}` + the `Rule{guard, action, fire_mode}` record; `fire_rules!` runs once/tick at step 10, `apply_action!`, `activate!`/`deactivate!`/`set_guard!`. `Invoke` is the eval-free general-code escape hatch (registry-by-name). | reads/writes place markings/params/tokens | none directly | fires at a fixed point in `_step!` |
| `src/ledger.jl` | 288 | The per-program (per-structured-token) ledger (MVP finding D): `ProgramLedger`, `program_ledger(state)::DataFrame`, `program_ledger_entries`; attributes cost at the bind site and reward at `finish!`, split across a transition's bound programs, with a network UNATTRIBUTED bucket, reconciling to the aggregate §8.5 rows. Included after `actions.jl` so it sees the token types/accessors. | reads run-state | reads structured-token agents | accrues economics onto in-flight programs during stepping |
| `src/serialize.jl` | 1049 | Single-JSON, eval-free model serialization (ADR 0005): `from_json_model` (`:167`) parses one JSON object → validates (`_validate_node!` `:593` / `_validate_action!` `:757`) → lowers each `ExprNode` dict via `to_expr` into the same `Expr` columns the DSL fills → constructs a `ReactionNetworkProblem`; `to_json_model` (`:259`) is the inverse; `build_network_from_dict` (`:78`), `node_to_dict`/`node_from_dict`/`model_to_dict`. NEVER `Meta.parse`s/`eval`s a model field (closes the import-time RCE). | writes the store from parsed IR | constructs the AA problem | none (authoring/serialization boundary) |
| `src/loadsave.jl` | 72 | Thin eval-free I/O wrapper: `@import_model`/`@export_model` (`:30`/`:48`) delegate to the `serialize.jl` typed-IR loader; `@export_solution_as_table`/`@export_solution_as_csv` (`:59`/`:70`) write the solution DataFrame. The old TOML/CSV/JLD2 model zoo + its import-time RCE are REMOVED. | reader/writer via serialize.jl | none | none (I/O) |
| `src/analysis.jl` | 557 | Analysis & observability (ADR 0013 / CONTRACT §14): `token_trajectory`/`representative_token`/`trajectory_envelope` read API + `push_token_trajectory_row!`; `ensemble`/`summarize`/`treatment_effect` + `EnsembleProblem <: AbstractAlgebraicAgent` (with `observables`/`getobservable`, `:431`/`:448`); the §15.1 recipe wrapper types (`MarkingPlot`/…) + `log_scalar_series`/`log_count_series`. | reads `program_ledger`/`prob.sol` | `@aagent EnsembleProblem`, AA read surface | reads finished state only |
| `src/export.jl` | 174 | Results export bundle (ADR 0013 §14.3): `export_run`/`export_ensemble` directory bundles (CSV+JSON core); the `_arrow_write` empty-generic weakdep seam + `_arrow_available`. | reads `prob.network` for the model hash | none | reads finished state only |
| `src/visualize.jl` | 295 | The three-layer network exec map (ADR 0014 §15.2): Layer A `network_graph`/`NetworkGraph`; Layer B `to_graphviz`/`draw_network`; Layer C `exec_map`. | walks `transLHS`/`transRHS` incidence | composes with AA's `run_graphviz`/`wiring_diagram` | Layer A runs `sample_transitions!` on a `deepcopy` (RNG-pure); else read-only |
| `ext/RDPlotsExt.jl` | 181 | The seven Plots.jl `@recipe`s + the relocated `AlgebraicAgents._draw` (ADR 0014 §15.1). Loads only when `Plots` is present (weakdep). | reads raw artifacts | defines `_draw` | none (post-hoc viz) |
| `ext/RDArrowExt.jl` | 21 | The SOLE `ReactiveDynamics._arrow_write` method (calls `Arrow.write`, ADR 0013 §14.3). Loads only when `Arrow` is present (weakdep). | none | none | none (I/O) |

> Line counts are `wc -l` as of this refresh; `file:line` citations are the current 1-indexed positions. `optim.jl` was removed (include commented at `src/ReactiveDynamics.jl:457`); its old `@optimize`/`@fit`/`@fit_and_plot`/`build_parametrized_solver`/`optim!` symbols are gone.

---

## Public-API Audit

The module root's only `export` is `ArcSpec, arcs, placename`; every other public symbol comes from an included file's `export` line. Symbols the root defines but does NOT export (`ReactionNetwork`, `SCHEMA`, `assign_defaults!`, `merge_networks!`, `add_obs!`, the store verbs `nrows`/`row_ids`/`column`/`cell`/`find_rows`/…, the type aliases) are public-by-use but UNEXPORTED — the ADR-0015 store-shim is deliberately internal (`RD.`-qualified).

**Every RD-owned export resolves — enforced as an automated invariant, not a hand-maintained list.** `test/semantic/exports_resolve.jl` iterates `names(ReactiveDynamics)` and asserts every RD-owned export is `isdefined` and (for macros) resolves. It SCOPES OUT symbols merely reexported from `@reexport using AlgebraicAgents` (a few of which are dangling upstream — AA's `@derived`/`@integration`/`AgentCall` — not RD's to police). The historical dangling exports (`@agentize`, `@prob_role`/`@list_by_role`/`@list_roles`, `@prob_check_verbose`, the old `@export`/`@import`) are all eliminated: the role/`@prob_check_verbose` exports were deleted, and `@agentize` was re-added as thin constructor sugar (defined at `solve.jl:45`).

The full export surface, by file (verbatim from the `export` lines):

| Export site | Exported symbols |
|---|---|
| `ReactiveDynamics.jl:15` | `ArcSpec`, `arcs`, `placename` |
| `state.jl:1` (reexport) | the full AlgebraicAgents surface — `@aagent`, `entangle!`, `getagent`, `inners`, `FreeAgent`, `AbstractAlgebraicAgent`, `simulate`, `step!`, `getobservable`, `draw`, … |
| `exprnode.jl:68-69` | `ExprNode`, `Const`, `NodeRef`, `Call`, `Sample`, `TimeRef`, `Choose`, `Field`, `ExternalRef`; `OP_WHITELIST`, `DIST_WHITELIST`, `REF_KINDS`, `to_expr`, `from_expr` |
| `agents.jl:1-3,57,158` | `AbstractStructuredToken`, `BaseStructuredToken`, `@structured_token`, `register_token_kind!`, `add_structured_token!`, `PopulationEntry`, `log_token_fields` |
| `create.jl:3,6,7` | `@reaction_network`, `@ReactionNetworkSchema` (deprecated alias), `@append_transitions` |
| `update.jl:3-7,176` | `@push`, `@name_transition`, `@mode`, `@add_place`, `@periodic`, `@jump`, `@prob_init`, `@prob_uncertainty`, `@prob_params`, `@prob_meta`, `@aka`, `@register`, `@cost`, `@reward`, `@valuation` |
| `solve.jl:13` | `@agentize` |
| `checkpoint.jl:15` | `StateDump`, `dump_state`, `restore` |
| `solvers.jl:4` | `ReactionNetworkProblem` |
| `predicates.jl` | (`TokenPredicate`/`select_tokens`/`matches`/`token_sortkey` are internal, `RD.`-qualified — no `export`) |
| `actions.jl:14-16` | `Rule`, `ActionStmt`, `SetMarking`, `SetParams`, `SetField`, `SetTokens`, `AddToken`, `Activate`, `Deactivate`, `Invoke`, `Log`, `Seq`, `apply_action!`, `fire_rules!`, `activate!`, `deactivate!`, `set_guard!` |
| `ledger.jl:57` | `ProgramLedger`, `program_ledger`, `program_ledger_entries` |
| `serialize.jl:10-11` | `node_to_dict`, `node_from_dict`, `model_to_dict`, `build_network_from_dict`, `from_json_model`, `to_json_model` |
| `loadsave.jl:8-9` | `@import_model`, `@export_model`, `@export_solution_as_table`, `@export_solution_as_csv` |
| `operators/equalize.jl:1` | `equalize!`, `@equalize` |
| `operators/joins.jl:2` | `merge_networks!`, `@join` |
| `operators/refine.jl:12` | `refine!`, `refine`, `abstract!`, `abstract_transitions`, `set_port_role!`, `@port`, `@compose`, `@pipeline`, `@process`, `compose`, `refinement_diagnostics`, `port_role` |
| `analysis.jl:17-18,482` | `token_trajectory`, `representative_token`, `trajectory_envelope`, `ensemble`, `summarize`, `treatment_effect`, `EnsembleProblem`, `MarkingPlot`, `SaturationPlot`, `ValuationPlot`, `LedgerPlot`, … (recipe wrapper types) |
| `export.jl:24` | `export_run`, `export_ensemble` |
| `visualize.jl:16` | `NetworkGraph`, `network_graph`, `to_graphviz`, `draw_network`, `exec_map` |

**Pseudo-macro markers (recognized by string, not defined as macros).** As with the modality markers (`@conserved`/`@rate`/`@nonblock`) and rate markers (`@ct`/`@deterministic`), the token-DSL markers `@select`, `@advance`/`@move`, `@structured`, `@field`, `@take` are NOT defined macros — they are recognized by `macroname(ex)` string comparison during parsing (`create.jl:320,334-337`, `reaction_parser.jl:45-69,97`) and lowered to `TokenPredicate`/`SetField`/genesis IR. They resolve only inside a `@reaction_network`/`@push` body, never as standalone macro calls.

---

## Static Store (`ReactionNetwork` / `SCHEMA`)

> ACSets/Catlab are GONE (ADR 0003); the ACSets-lineage names are retired (ADR 0015). The store is a dependency-free typed-struct-of-columns defined in the module root `src/ReactiveDynamics.jl`. `GeneratedExpressions` and its brace-comprehension pass were dropped.

- **`const SCHEMA`** (`ReactiveDynamics.jl:45`) — a `NamedTuple` of the six objects (`:S` places, `:T` transitions, `:E` events, `:obs` observables, `:P` params, `:M` meta — zero homs), each mapping `attr ⇒ element-type`. The single source of truth for the object model; declaration order is load-bearing (it fixes column order for the `"place"`/`"trans"`-substring reflection loops). Derived: `ATTR2OBJ` (`:81`), `ALLATTRS` (`:84`), `columns(SCHEMA[, obj])` (`:89-90`).
- **`struct AttrColumn{T}`** (`:95`) — a typed column `v::Vector{T}` + a `def::Vector{Bool}` presence mask; an undefined cell reads `nothing`.
- **`struct ReactionNetwork`** (`:145`) — the static network container: `counts::Dict{Symbol,Int}` (rows per object), `columns::NamedTuple` (typed `AttrColumn`s in `ALLATTRS` order), `arcs::Vector{ArcSpec}` (the promoted incidence table). A typed inner constructor (`:153`) prevents the auto field-ctor from colliding with the legacy 3-arg semantic constructor. Content-based `Base.hash` (`:206`) EXCLUDES the derived `arcs` table so the model hash is stable across promotion. `Base.@deprecate_binding ReactionNetworkSchema ReactionNetwork` (`:168`) keeps old annotations compiling for one release.
- **Store verbs (UNEXPORTED, `RD.`-qualified):** `nrows` (`:220`), `row_ids` (`:221`), `col_row_ids` (`:222`), `getindex`/`setindex!` in 4 forms (`:225-231`), `column`/`cell`/`set_cell!` (`:234-236`), `find_rows` = `findall(isequal)` (`:240`, a `findall`, NOT an FK follow), `add_row!` (`:243`)/`add_rows!` (`:254`)/`rem_rows!` (`:266`). `[:,attr]`/`column` return COPIES. `rem_rows!` is SWAP-AND-POP (last row fills each freed slot; `equalize!` relies on the surviving-row order).
- **Port roles (ADR 0009 §A):** `const PORT_ROLES = (:private, :input, :output, :shared)` (`:187`); `port_role(net, i|name)` (`:192-196`) reads the `placeRole` attr (default `:private`).

Object/attribute model (element types, from `SCHEMA`):

| Object | Attributes (col → type) |
|---|---|
| `:S` places | `placeName`→Symbol, `placeModality`→`Set{Symbol}`, `placeInitVal`/`placeInitUncertainty`/`placeCost`/`placeReward`/`placeValuation`→SampleableValues, `placeStructured`→Bool, `placeRole`→Symbol (ADR 0009 §A port role, default `:private`) |
| `:T` transitions | `trans`/`transPriority`/`transRate`/`transCycleTime`/`transProbOfSuccess`/`transCapacity`/`transMaxLifeTime`/`transPreAction`/`transPostAction`/`transMultiplier`→SampleableValues, `transName`→`Union{String,Symbol,Missing}` |
| `:E` events | `eventTrigger`/`eventAction`→SampleableValues |
| `:obs` observables | `obsName`→Symbol, `obsOpts`→`FoldedObservable` |
| `:P` params | `prmName`→Symbol, `prmVal`→Any |
| `:M` meta | `metaKeyword`→Symbol, `metaVal`→SampleableValues |

**Promoted arc relation (ADR 0003 Phase 2):** `struct ArcSpec` (`:135`) — `trans::Int` (FK→`:T`), `place::Int` (FK→`:S`, 0=dynamic), `stoich` (the arc weight), `side::Symbol`, `modality::Set{Symbol}`, `expr` (escape-hatch for a dynamic arc) — a `Vector{ArcSpec}` FIELD on the store (NOT a 7th SCHEMA object, so it stays out of the reflection loops). Derived by `populate_arcs!` (serialize.jl) from the `:trans` column via the eval-free static arc decomposition; `equalize!`/`refine` re-derive it post-transform for a structural FK-repoint. Accessors: `arcs`, `placename`, `find_index` (exported from the root).

Supporting root tables/routines: `prettynames` (user alias → canonical attr), `defargs` (per-object defaults, incl. `:placeStructured => false`, `:placeRole => :private`), `place_modalities` (`[:nonblock, :conserved, :rate]`), `assign_defaults!`, `add_obs!`, `merge_networks!`, the constructors. String→attribute coercion is `convert(Symbol, String)` on assignment (the eval-based Set/FoldedObservable convert hooks were removed with ADR 0005).

---

## Stepping Trace

Stepping is owned entirely by the native engine in `src/solvers.jl`, driven through the AlgebraicAgents interface. There is no separate driver: both "standalone" and "AA-wrapped" usage call AA's reexported `simulate(prob, max_t=Inf)`, which loops `step!` while the projected time is continuable and below `max_t`. RD defines no `simulate` of its own. A run is fully determined by `(model, seed)` (CONTRACT §4) — every draw routes through the state-owned `state.rng` (seeded from the `seed=` kwarg, snapshotted as `initial_rng` for `_reinit!`), never the global RNG.

### Driver chain (AA → engine)
1. `simulate(prob[, max_t])` (AA) loops `step!` while `iscontinuable(ret) && ret < max_t`.
2. `step!` descends into INNER agents first, then runs the LOCAL `_step!(state)` only when `_projected_to(state)` equals the least projected time across the hierarchy. So co-stepped siblings (a user controller, coupled AA agents) step before the network body for a given tick. Structured tokens are clock-inert (`_step!`/`_projected_to` no-ops, `agents.jl:138-139`).
3. **Entry point: `AlgebraicAgents._step!(state::ReactionNetworkProblem)`** at `solvers.jl:1170`.

### One tick — `_step!` body (`solvers.jl:1170-1218`), in order
1. `update_u_structured!` — sync structured-place counts from the agent layer.
2. `save!(state)` only if `sol` empty — record the initial `(t, u)` row.
3. `free_blocked_places!` (`:809`) — release `:nonblock` resources.
4. `update_u_structured!`.
5. `update_observables` — resample observables whose interval elapsed.
6. `sample_transitions!` (`state.jl:269`) — clear `state.transitions`, eval activated recipes (guard AND-ed with the `transActivated` latch), unfold LHS into `UnfoldedArc`.
7. `evolve!` (`:305`) — SPAWN + ALLOCATE (progressive-fill) + advance ongoing (detail below).
8. `update_u_structured!`.
9. `finish!` (`:668`) — TERMINATE + PoS + emit RHS + return `:conserved`/`:nonblock` + prune (detail below).
10. `update_u_structured!`.
11. **`fire_rules!` (`actions.jl`)** — the endogenous decision channel (ADR 0010), evaluated once/tick; replaces the old no-op `event_action!` slot. Rules see this tick's post-finish state; their writes land on this tick's ledger row and the next tick's genesis/guards.
12. `update_u_structured!`.
13. push the aggregate `:valuation` ledger row.
14. `attribute_valuation!` + `push_program_ledger_row!` (`ledger.jl`) — mark each live program to market and push the per-program ledger row in deterministic token order.
15. `push_token_trajectory_row!` (`analysis.jl`) — snapshot each opted-in token's declared fields at the same seam/order (ADR 0013 §14.1).
16. **Advance clock**: `state.t += state.dt` — the SINGLE clock advance, after all tick work.
17. `save!(state)` — record the post-increment `(t, u)` row.

Termination/progress: `_projected_to` (`solvers.jl:1219`) returns `true` once `state.t > state.tspan[2]`, else `state.t`.

### Allocation under contention — the ADR-0002 weighted progressive-fill allocator
`build_requirements!` (`:100`) assembles the demand matrix; `progressive_fill!` (`:162`, with `_fill_tier!` `:191`) does priority-weighted water-filling — work-conserving, deterministic, conjunctive-consistent; `spawn_integer!` (`:268`) is the integer spawn wrapper (progressive-fill with `fmax = q_desired`). This replaced the old per-resource-split allocation tangle.

### Spawning + completion
- `evolve!` (`:305`): per transition, compute the (Poisson/deterministic) desired count, apply `transCapacity`, allocate init requirements via the progressive-fill allocator, floor to whole instances, construct `Transition` heap entries (`t=state.t`, `q`, `state=0.0`), bind structured tokens in `select_tokens`/`token_sortkey` order, run `transPreAction`; ongoing instances then advance by `q*dt`.
- `finish!` (`:668`): terminate on lifetime-exceeded OR cycle-complete; draw `Binomial(q, PoS)` success; emit RHS products (incl. `structured_rhs` `:561` for `@structured` genesis); return `:conserved`/`:nonblock` resources; run `transPostAction`; prune terminated instances.

### Construction-time modality validation
`validate_modalities(net)` (`solvers.jl:856`) runs in the constructor BEFORE any closure compiles or tick runs, rejecting the three illegal §1.4 configs (`{:nonblock, :conserved}`; `:rate` with a concrete `cycletime == 0`; `:rate` on a structured place) with a clear `ArgumentError` (CONTRACT §1.4).

### Clock / dt notes
- Single authoritative clock `state.t`, advanced only at the end of `_step!`.
- The runtime field is `dt` (`state.jl`), and the accepted meta keyword is `dt`; the legacy `tstep` meta keyword is an accepted-but-deprecated alias mapped to `dt` (depwarn, one release). `get_tcontrol` (`:826`) resolves `tspan`/`dt` from `@prob_meta` or the constructor kwargs.

---

## AlgebraicAgents Touchpoints

AA is the sole simulation substrate; RD is a node in the AA hierarchy. AlgebraicAgents is the published registry release **0.4** (`[compat] AlgebraicAgents = "0.4"`); the earlier `Merck/AlgebraicAgents.jl@main` `[sources]` pin was dropped.

### Import / reexport
- `@reexport using AlgebraicAgents` at **`state.jl:1`** is the single AA entry point, reexporting the full AA public surface to RD users (`@aagent`, `entangle!`, `getagent`, `inners`, `FreeAgent`, `AbstractAlgebraicAgent`, and the driver/viz wrappers `simulate`, `step!`, `getobservable`, `draw`). `GeneratedExpressions` is NO LONGER reexported (dependency dropped).

### `@aagent` structs (4) + abstract types
| Type | Location | Notes |
|---|---|---|
| `Transition` | `state.jl:37` | in-flight instance; holds resource-agent fields |
| `Observable` | `state.jl:54` | per-observable sampler agent |
| `ReactionNetworkProblem` | `state.jl:63` | the top-level agent = the engine state |
| `EnsembleProblem` | `analysis.jl:262` | ensemble runner as a drawable/walkable AA node |
| `BaseStructuredToken` (`@aagent FreeAgent`) | `agents.jl:11` | base structured-token struct |
| `AbstractStructuredToken <: AbstractAlgebraicAgent` | `agents.jl:6` | root of the token hierarchy (abstract type) |

`@aagent` injects `(uuid, name, parent, inners, relpathrefs, opera)` and a name-first constructor; user token subtypes are generated programmatically via `AlgebraicAgents.aagent(...)` (`agents.jl:30`).

### AA interface methods implemented
| Method | Type | Location |
|---|---|---|
| `_step!` | `ReactionNetworkProblem` | `solvers.jl:1170` |
| `_reinit!` / `reinit!` | `ReactionNetworkProblem` | `solvers.jl:1080` / `:1150` (`reinit!` takes a `seed=` kwarg for ensemble mode (b) reseed) |
| `_projected_to` | `ReactionNetworkProblem` | `solvers.jl:1219` |
| `observables` / `getobservable` | `ReactionNetworkProblem` | `aa_coupling.jl:40` / `:59,71,74` (by name/string/index — ADR 0012 §A) |
| `_getparameters` / `_setparameters!` | `ReactionNetworkProblem` | `aa_coupling.jl:87` / `:100` |
| `_prestep!` | `ReactionNetworkProblem` | `aa_coupling.jl:129` (latches incoming AA wires into `ExternalRef` leaves; one-tick Jacobi lag) |
| `_step!` / `_projected_to` (no-op) | `AbstractStructuredToken` | `agents.jl:139` / `:138` |
| `_step!` / `_projected_to` (no-op) | `EnsembleProblem` | `analysis.jl:421` / `:422` |
| `observables` / `getobservable` | `EnsembleProblem` | `analysis.jl:431` / `:448,461,464` (cross-run aggregate reductions) |
| `_draw` | `ReactionNetworkProblem` | `ext/RDPlotsExt.jl` (weakdep; reached via exported `draw`) |

The ADR-0012 AA coupling surface (`getobservable`/`observables`/params outbound, `ExternalRef` + `_prestep!` inbound) IS implemented in `aa_coupling.jl` — a reactive network can be a wire SOURCE and can read sibling/parent state through declared wires. What remains future work (ADR 0012): AA `Opera`-level implicit/fixed-point (within-tick algebraic-loop) coupling — the current coupling is explicit one-tick-lag Jacobi.

### Agentization & token hierarchy
- **`@agentize` is exported (`solve.jl:13`) AND defined (`:45`)** as thin constructor sugar — it expands to exactly one `ReactionNetworkProblem(net[, u0, p]; …)` call (no second construction path), auto-naming from a bare-binding net. The core public contract is `prob = ReactionNetworkProblem(net[, u0, p]; seed, registry, population, name, …)` then `simulate(prob[, max_t])`.
- Structured tokens are first-class AA agents under a `"structured"` `FreeAgent` container; attached via `entangle!` in `add_structured_token!` (`agents.jl:46`); the solver reaches them via `inners(getagent(state, "structured"))`. The declarative `population[]` initial marking is instantiated at construction by `instantiate_population!` (`agents.jl:86`) in DECLARED order, assigning per-place creation indices through the seeded RNG (ADR 0007 §B). A token's `bound_transition` points to its in-flight `Transition`; the network↔token link is the AA `parent`/`inners` hierarchy.
