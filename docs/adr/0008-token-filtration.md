# ADR 0008 — Agentic species under a filtration: predicate-based token selection (`TokenPredicate`, `@select`)

- Status: **Implemented** (proposed 2026-06-21; revised 2026-06-21 — maintainer ruling: PHASE IS AN ATTRIBUTE, canonically). Generalizes a structured-token LHS reactant from "select by KIND only" to "select the 𝓕ₜ-measurable sub-population a predicate carves out," via one closed eval-free IR node (`TokenPredicate`) and one DSL form (`@select`). LANDED: `TokenPredicate{kind, clauses}` (`src/predicates.jl:17`) as the binding filter, `SetField`/`@advance` the single-token field write and `SetTokens` its population generalization (`src/actions.jl:35`), phase-as-attribute canonical. Commit `896bc45` (Phase-1 Stage C); exercised by `test/semantic/token_filtration.jl`. Answers the maintainer's Q4 directly: "a transition may want to take a project in a given state / subject to a filtration expression — more complex than a previous species." The maintainer has ruled that a token's lifecycle state (pipeline phase, status, etc.) is modeled as a TOKEN ATTRIBUTE, not as a distinct species KIND — so this ADR adopts phase-as-attribute as THE canonical structured-token model. The KIND axis (`get_species(a)==kind`) remains the degenerate predicate clause, so this is a unification (one bind path), not a second mechanism; advancing a token's state is a field write (`@advance`/`SetField`, §D), generalizing the existing `@move` (which writes the `species` field) to write any field.
- Verification note: `file:line` read against `ref-agents` (Julia 1.12.5). The binding filter that this ADR generalizes is literally `filter(a -> get_species(a) == type && !isblocked(a), structured_token)` at the genesis pass (`src/solvers.jl:202-205`) and again at the ongoing pass (`src/solvers.jl:274-277`); both then `sort!` by `priority(a, transName)` descending (`:207-211`, `:279-283`) and take `allocs[j,i]` instances. The selector is KIND + not-blocked; there is no attribute filter anywhere on the bind path. The `priority` hook (`src/interface/agents.jl:67`) lets you ORDER by an attribute but never SELECT by one.
- Relates to: ADR 0006 (structured tokens: the token protocol `get_species`/`set_species!`/`isblocked`/`priority`/`past_bonds`, `agents.jl:49-72`; the `(species, creation_index)` total order; the `TokenAgg{reducer,species,field,active_only}` READ node — this ADR adds the BIND analogue `TokenPredicate` and the field-WRITE primitive `SetField`/`@advance` that generalizes `@move`'s `set_species!`), ADR 0005 (closed `ExprNode` IR, `OP_WHITELIST`/`DIST_WHITELIST`, eval-free `validate`; the action family `{SetSpecies,SetParams,Log,Seq}` which this ADR extends with `SetField`; `TokenPredicate`/`Field` are new closed nodes like `TokenAgg`), ADR 0002 (priority-weighted allocation — the predicate filters the candidate set the allocator then orders), ADR 0007 (the initial marking that creates the tokens a predicate selects over; the `population[]` instances carry the `phase` field this ADR makes canonical). CONTRACT §1 (modalities — `perstep` still illegal for structured species), §4 (D1/D4/D5 — the predicate is 𝓕ₜ-measurable, never reads the future or a later-in-tick RNG draw; a `SetField` write may draw, consuming the seeded stream deterministically), §9 (the structured-token subsystem this extends). Becomes a CONTRACT §9.5 + an ADR 0005 amendment (action family + two nodes).

## Context

A structured species today is a flat KIND. An LHS reactant naming it consumes/holds tokens chosen purely by `get_species(a) == kind && !isblocked(a)` (`src/solvers.jl:202-205`), ordered by the host-overridable `priority(a, transition)` (`agents.jl:67`), taking the granted integer count (ADR 0006). You can therefore answer "give this transition any unbound Project, preferring high-priority ones" but NOT "give it a Project whose `phase == :Phase2` and `npv_estimate > θ`." The maintainer's Q4 is exactly this: a transition wants a token "in a given state / subject to a filtration expression."

Three ways the current engine forces you to fake it, all bad:
1. **Phase-as-species explosion.** Encode every selectable state as its own KIND (`:Project_Phase2_HighNPV`). Combinatorial blow-up; `npv > θ` with continuous `θ` is impossible; `@move` (`src/solvers.jl:360-389`) must shuttle tokens through a lattice of kinds.
2. **Smuggle the test into `priority`.** Return `priority = -Inf` for non-matching tokens so they sort last — but they are still SELECTED if the granted count exceeds the matching count, so a transition can bind a token it should have rejected. `priority` orders; it cannot exclude.
3. **Host-side pre-filter.** Not expressible in the eval-free JSON model at all (ADR 0005); it would live in host code, breaking the §8.2 S2 "(model.json, seed) determines the run" guarantee, exactly the failure ADR 0007 §B fixes for initial populations.

The maintainer's own word — *filtration* — is the right formalism and the ADR honors it literally. In a stochastic process the natural filtration 𝓕ₜ is "everything knowable by time t." A token-selection predicate MUST be 𝓕ₜ-measurable: it may read token fields, params, observables, the clock, and sibling-pool aggregates AS OF the allocation point, but never the future and never an RNG draw made later in the same tick — otherwise determinism (CONTRACT §4 D1) is ill-defined because the bind set would depend on draw order.

There is a latent half of this already in ADR 0006: `TokenAgg{reducer, species, field, active_only}` lets a RATE READ the token pool ("number of active Phase-2 projects") non-consumingly. What is missing is the BIND analogue — a way to SELECT which tokens a transition consumes. This ADR adds it, deliberately mirroring `TokenAgg` so the two share one predicate vocabulary.

## Decision

### (A) One closed, eval-free IR node: `TokenPredicate`

Add to the ADR 0005 `ExprNode` family one node, structurally identical in spirit to `TokenAgg`:

```julia
TokenPredicate{ kind::Symbol,                 # the structured species (FK → a specStructured Species)
                clauses::Vector{Clause} }      # AND-joined; empty ⇒ degenerate "any token of kind"
Clause = (field::Symbol,                       # a token field: protocol (species) or host-struct extra (phase, npv…)
          op::Symbol,                           # ∈ PRED_OP_WHITELIST
          value::ExprNode)                      # RHS: Const | Ref{param|obs} | TimeRef | TokenAgg | arithmetic Call
const PRED_OP_WHITELIST = (:(==), :(!=), :(<), :(<=), :(>), :(>=), :in)
```

Clauses are conjunctive (AND). Disjunction is expressed as two reactant rows or a future `any`/`all` wrapper — kept out of v1 to hold the node closed and JSON-Schema-enumerable (ADR 0005). `field` is read via the same `hasproperty`-guarded accessor as `tokenattr` (ADR 0006 §9.3) so a typo is a `validate()` diagnostic, never AA's silent-`false` swallow (`AlgebraicAgents.jl/src/queries.jl:107-113`). `value` is itself an `ExprNode` over the closed whitelist, so `npv > 1.5 * θ` (θ a param) and `phase == @obs(current_gate)` are expressible; `value` MUST NOT contain a `Sample` node (an RNG draw inside the selection predicate would make the bind set draw-order-dependent — §C measurability).

### (B) The structured-LHS reactant carries an optional predicate

The promoted `ReactantSpec` row (CONTRACT §6.5) gains, for LHS rows on a structured species, an optional `predicate::Union{Nothing, TokenPredicate}`. The runtime `UnfoldedReactant` (`src/state.jl:4-9`) likewise carries it. Binding generalizes the `solvers.jl:202-211` / `:274-283` filter from

```julia
filter(a -> get_species(a) == type && !isblocked(a), structured_token)        # KIND only (today)
```

to

```julia
filter(a -> get_species(a) == type && !isblocked(a) && matches(predicate, a, state, transition),
       structured_token)                                                       # KIND + filtration
```

then the EXISTING sort and take are unchanged: `sort!` by `(priority(a, transName) desc, creation_index)` (ADR 0006 invariant-5 total order — note the creation-index tie-break this ADR depends on), take `allocs[j,i]` instances. **The predicate is just a `filter` in front of the existing `sort`** — it slots into the ADR 0002 allocator (which orders the candidate set the predicate produced) and the ADR 0006 determinism (the sort key is unchanged) with zero new machinery on the allocate/order side.

A `nothing` predicate (or empty `clauses`) is exactly today's behavior, so this is strictly backward-compatible: every existing structured reaction keeps binding by kind.

### (C) The predicate is 𝓕ₜ-measurable — the determinism rule

`matches(pred, token, state, transition)` evaluates each clause's `value` `ExprNode` through `context_eval` against the snapshot AT THE ALLOCATION POINT and compares with the named op. Measurability is enforced by three rules, all checkable in `validate`:

1. **No future, no later-tick draw.** `value` MUST NOT contain `Sample` (RNG) — the bind set cannot depend on a draw made during the same tick's evaluation, else it is order-sensitive and D1 is ill-defined. (`token field op param/obs/time/aggregate` is the legal shape.)
2. **Pinned observation point.** Token attributes and `TokenAgg` aggregates inside `value` are read on the `state.u`-consistent snapshot established by the `update_u_structured!` immediately before `evolve!` (`src/solvers.jl:649` for the genesis pass, `:653` for ongoing). This is the SAME pin ADR 0006 left as an open question for `TokenAgg`; this ADR resolves it for both: **token-reading expressions observe the pre-`evolve!` reflected counts.** A predicate and a rate that read the same sibling pool in the same tick see identical numbers.
3. **Total-order before take.** The filtered candidate set is ordered by the ADR 0006 invariant-5 key `(priority desc, creation_index)` before the integer take, so WHICH matching tokens are bound (when more match than are granted) is reproducible across runs and insertion histories — the predicate narrows the set, the total order makes the choice within it deterministic.

Together: for fixed `(model, seed)`, the set of tokens a transition binds this tick is a deterministic function of the 𝓕ₜ-measurable state — the contract's reproducibility (§4 D1) extends cleanly to predicate selection.

### (D) Phase is an attribute — the canonical model, and the `SetField`/`@advance` write primitive

**Maintainer ruling: a token's lifecycle state is a TOKEN ATTRIBUTE, not a species KIND.** A pipeline program is ONE structured KIND — `Project` — carrying a `phase::Symbol` field (`:Discovery,:Phase1,…,:Market`) plus `npv_estimate`, `cost_to_date`, `pos_remaining`, etc. Selection of "a project in phase 2" is the predicate `@select(Project, phase==:Phase2)` (§A/§B); ADVANCING a project to the next phase is a WRITE to its `phase` field on the transition's RHS. This is now THE structured-token model; the per-phase-KIND style (a distinct species per phase, `@move`-shuttled — the original ADR 0006 north-star sketch) is retired as the canonical pattern.

Why this is a unification and not a new mechanism: `get_species(a)==kind` is just the degenerate predicate clause `(species, ==, kind)` (§A), so KIND selection and attribute selection are one filter. Symmetrically, the existing `@move(:from,:to)` already WRITES a token field — it calls `set_species!(token, to)` then unbinds (`solvers.jl:378,384`). Phase-as-attribute needs exactly the same operation pointed at a DIFFERENT field. So we generalize `@move` into one field-write primitive rather than adding a parallel path:

**`SetField{field, value}` — a new action statement.** It joins the action *type family*, which [ADR 0010](0010-rules-and-conditional-transitions.md) (rules/conditional transitions) independently extends; the reconciled canonical family is `{SetSpecies, SetParams, SetField, AddToken, Activate, Deactivate, Log, Seq}`, defined once in ADR 0010 §C. `SetField` writes a field of the bound token(s) of the firing instance: `field::Symbol` (a token field — `phase`, `npv_estimate`, or the protocol `species` field), `value::ExprNode` (evaluated through `context_eval`; may read the token's own current fields, params, observables, time, and MAY draw via `Sample` since a write may legitimately resample, §F). It is applied in `finish!` to each successful instance's bound tokens, reusing the `@move` write/unbind machinery (`solvers.jl:360-389`). **Context:** because it writes the firing instance's bound token, `SetField` is legal only in a transition **post-action** (and `@move`'s structured-RHS path) — NOT in a standalone Rule action (a Rule has no bound token; ADR 0010 §C restriction).

- `@advance(phase, :Phase3)` ⇒ `SetField{:phase, Const(:Phase3)}` — the canonical phase-advance; writes `phase`, keeps the token's identity, `species` (KIND), `uuid`, `creation_index`, and `past_bonds` history intact (ADR 0006 inv 1, 6).
- `@move(:from, :to)` ⇒ `SetField{:species, Const(:to)}` — the SAME primitive writing the protocol `species` field. `@move` is retained as sugar for the degenerate "write the species field" case (back-compat for any per-KIND model), and its redundant re-`entangle!` (`solvers.jl:376`, flagged for removal in ADR 0006 §D) is dropped — a `SetField` write never changes the `inners` container key (the key is the uuid, ADR 0006 inv 6), so no re-entangle is ever needed.
- `@advance(npv_estimate, npv_estimate * uplift)` ⇒ `SetField{:npv_estimate, Call{:*, [Field(:npv_estimate), Ref{param}(:uplift)]}}` — a general attribute update (e.g. revalue a program on phase transition), the thing a `phase`-only `@move` could never express.

A new `Field{name}` ExprNode leaf (the bound-token field read) joins the closed family so a `SetField.value` can reference the token's own current attributes; it resolves against the firing instance's bound token at apply time and is rejected by `validate` outside a `SetField`/token context. This keeps the whole advance path eval-free and JSON-Schema-enumerable, exactly like `TokenPredicate`.

**Consequences of canonicalizing phase-as-attribute:** (1) no KIND explosion — one `Project` kind instead of `Project_Phase1…Project_Market`; (2) continuous and cross-field predicates (`npv>θ`, `cost_to_date<budget_cap`) are expressible (impossible under per-KIND); (3) `past_bonds` history is uniform on one identity across all phases (ADR 0006); (4) the acquisition lever injects one `Project(phase=:Phase2,…)` rather than choosing among phase-KINDs; (5) the structured `state.u` column per phase is now derived via a `TokenAgg`-style count `nactive(Project, phase==:PhaseN)` rather than one reflected column per KIND — see the §9.5 contract delta and ADR 0006 §A `update_u_structured!`. The per-KIND style is not forbidden (the predicate makes it a legal degenerate), but it is no longer the documented or recommended pattern, and the BD template, examples, and tests use phase-as-attribute exclusively.

### (E) DSL surface

The canonical phase-as-attribute pipeline step — select a Phase-2 Project with NPV over a param θ, consume renewable resources for a cycle, advance it to Phase 3 on success:

```julia
# canonical: one Project kind, phase is a field, @select filters, @advance writes the field
@ct(cycletime), @select(Project, phase == :Phase2 && npv_estimate > θ) + @rate(budget) --> @advance(phase, :Phase3)

# a multi-field advance (revalue on transition) is the same primitive:
@ct(ct), @select(Project, phase == :Phase2) --> @advance(phase, :Phase3); @advance(npv_estimate, npv_estimate * uplift)

# @move is retained as sugar for the degenerate "write the species field" case (per-KIND back-compat):
@ct(ct), Project_Phase2 --> @move(:Project_Phase2, :Project_Phase3)   # ⇒ SetField{:species, …}
```

`@select(kind, clause && clause …)` parses to a `TokenPredicate`: the `&&` tree splits into `clauses`, each `field op value` mapped to a `Clause`, `value` lowered to an `ExprNode` (params/obs/time allowed, `Sample` rejected by `validate`). `@advance(field, value)` parses to a `SetField{field, value}` RHS action (§D); `@move(:from,:to)` is its `field=:species` sugar. Stoichiometry `q * @select(...)` keeps its meaning (bind `q` matching tokens). The parser hooks are `recursive_find_reactants!`'s structured/macrocall branch (`src/interface/reaction_parser.jl:66-76`, `src/interface/create.jl:300-301`), which already special-cases `@structured`/`@move` — `@select` and `@advance` join that closed set — and the structured-RHS path (`structured_rhs`, `solvers.jl:333-397`) that today handles `@move`.

### (F) `validate` additions (ADR 0005 amendment)

`validate(spec)` (CONTRACT §8.3) gains rule 7: for every `TokenPredicate` — `kind` is a declared `specStructured` species; each `Clause.op ∈ PRED_OP_WHITELIST`; each `Clause.value` is a whitelisted `ExprNode` containing NO `Sample` node (§C rule 1); each `Clause.field` is a known token field for that kind where statically knowable (protocol fields always; host-struct extras checked against a registry-declared field list if the kind registers one, else deferred to the `hasproperty` guard at runtime with a warning). For every `SetField{field, value}` (§D): `field` is a known token field of the bound species (same static/`hasproperty` check), and `value` is a whitelisted `ExprNode` — a `Sample` IS permitted in a `SetField.value` (unlike a predicate `value`) since a field WRITE in `finish!` may legitimately draw (e.g. resample an estimate), but it then consumes the seeded RNG stream at a deterministic point (§4 D5); `value` may contain `Field{name}` leaves resolved against the bound token. `@select`/`@advance`/`Field` outside a structured context is a diagnostic (you cannot predicate-filter or field-write a `Float64` count). `perstep` modality on a predicated structured reactant stays illegal (CONTRACT §1.4 — unchanged).

## Consequences

- The maintainer's Q4 is answered with one selection node (`TokenPredicate`) and one write primitive (`SetField`), both eval-free and reproducible: a transition takes "a project in a given state subject to a filtration expression" and advances its state by writing a field.
- Phase-as-attribute is canonical: ONE `Project` kind with a `phase` field replaces a lattice of phase-KINDs, and continuous / cross-field selection criteria (`npv > θ`, `cost_to_date < budget_cap`) become expressible for the first time. The per-KIND style survives only as a legal degenerate (the `species`-field predicate + `@move` sugar), not the recommended pattern.
- `@move` stops being a special case: it is `SetField{:species, …}` sugar, so the engine has ONE field-write path for the protocol `species` field and every host-struct field, and the redundant re-`entangle!` (`solvers.jl:376`) is removed (the uuid container key never changes on a field write).
- It composes with ADR 0006: `TokenPredicate` (bind/consume) and `TokenAgg` (read/aggregate) share the field-access and observation-point discipline; a rate can READ "n active Phase-2 projects" (`TokenAgg`) and a transition can BIND "the highest-priority Phase-2 project with npv>θ" (`TokenPredicate`) consistently in one tick.
- The IR grows by: one selection node (`TokenPredicate`) + `PRED_OP_WHITELIST`; one action statement (`SetField`, extending the ADR 0005 action family to `{SetSpecies, SetParams, SetField, Log, Seq}`); and one ExprNode leaf (`Field{name}`, the bound-token field read). All additive, JSON-Schema-enumerable, eval-free — same discipline as `OP_WHITELIST`/`TOKEN_REDUCER_WHITELIST`.
- Costs: the bind path gains a `matches` filter at four sites (`solvers.jl:202,274` filters; the predicate must thread to `UnfoldedReactant` via `sample_transitions!`, `state.jl:194-205`); the `finish!` RHS path generalizes `@move`'s `set_species!` (`solvers.jl:378`) to a `SetField` field write over bound tokens; `validate` rule 7 (predicate + `SetField`); the `@select`/`@advance` parser branches; and the observation-point pin must be wired (resolving the ADR 0006 open question for `TokenAgg` at the same time). The per-phase reflected `state.u` column (ADR 0006 §A `update_u_structured!`) is replaced for `Project` by phase-filtered counts — a derived-view change to confirm against the engine. Disjunction, cross-token-joins ("a Project AND its paired Patent"), and arg-count-varying clauses are explicitly OUT of v1.

## North-star tie-in (BD)

One `ProjectToken <: AbstractStructuredToken` carrying `phase, npv_estimate, cost_to_date, pos_remaining, acquired, acq_time`. Pipeline advance: `@ct(ct_p2), @select(Project, phase==:Phase2) + 3*@conserved(scientist) + @rate(budget) --> @advance(phase, :Phase3)` on `Binomial(q, PoS)` success — ONE `Project` kind, `phase` is a field, no per-phase kinds. The acquisition lever (ADR 0007 §10.5 / ADR 0006) injects `ProjectToken(phase=:Phase2, acquired=true, npv_estimate=…)`; a `@select(Project, acquired==true && phase==:Phase3)` transition can then apply post-acquisition operational changes — e.g. `@advance(pos_remaining, pos_remaining * 1.1)` — to ONLY the acquired late-stage programs, a filtration + field-write the engine cannot express today. rNPV reductions (`TokenAgg`/`sumattr` over `active(prob, :Project)` filtered by phase) read the same pinned snapshot the bind predicates use, so the pre/post-lever diff (same seed, §4 D5) is well-defined.

## Invariants

1. **Filtration measurability.** A `TokenPredicate.clauses[*].value` contains no `Sample` node and references only token fields, params, observables, `TimeRef`, and `TokenAgg` — all 𝓕ₜ-measurable at the allocation point (validate rule 7).
2. **Pinned observation point.** Token-reading subexpressions in a predicate (and in a `TokenAgg`) read the pre-`evolve!` reflected counts (`update_u_structured!`, `solvers.jl:649/653`); a predicate and a rate reading the same pool in one tick agree.
3. **Deterministic bind set.** For fixed `(model, seed)`, the bound token set per tick is reproducible: predicate narrows, `(priority desc, creation_index)` total order (ADR 0006 inv 5) selects within, integer take bounds.
4. **Degenerate compatibility.** `predicate == nothing` (or empty clauses) is byte-for-byte today's kind-only bind; every existing structured reaction is unaffected.
5. **Bind ≠ read.** `TokenPredicate` consumes/holds (CONTRACT §1 modality applies, `perstep` still illegal); `TokenAgg` reads non-consumingly. The two never alias: a `TokenAgg` in a rate never blocks a token, a `TokenPredicate` on an LHS always can.

## Open questions

- **Field-schema declaration.** Should a structured KIND declare its selectable/writable field list (types) at registration (ADR 0006 C registry) so `validate` can statically type-check both `Clause.value` against `Clause.field` (predicate) and `SetField.value` against `SetField.field` (write), or is the runtime `hasproperty` guard + a soft warning enough? Static is safer for agentic authoring; it costs a richer registry entry. With phase-as-attribute canonical, this matters more — `phase`/`npv_estimate`/etc. are now the primary modeling surface, not protocol internals.
- **Per-phase reflected count.** ADR 0006 §A reflects ONE `state.u` column per structured KIND via `update_u_structured!` (`solvers.jl:630-639`). Under phase-as-attribute there is one `Project` KIND but the model wants per-phase counts (`#Phase2 projects`). Decide whether `update_u_structured!` gains phase-group columns, or whether per-phase counts live purely in `TokenAgg`/`nactive(Project, phase==…)` queries with the single `Project` column counting all active. Affects the §9.2-inv-3 `length(active)==state.u` assertion (it must become per-predicate, not per-KIND).
- **`SetField` ordering within `finish!`.** When one transition writes multiple fields (`@advance(phase,…); @advance(npv,…)`) the writes must apply in a defined order (declaration order) and against a consistent token snapshot; pin it so a `SetField.value` reading `Field{:npv}` sees the pre-write or post-write value deterministically.
- **Disjunction / nesting.** v1 is conjunctive. If real BD models need `(phase==:Phase2 || phase==:Phase3)`, decide between (a) multiple reactant rows, (b) an `in` clause `phase ∈ [:Phase2,:Phase3]` (already covered by `:in`!), or (c) a boolean-tree node — prefer (b)/(a) before (c) to keep the node closed.
- **Cross-token predicates.** "Bind a Project together with its paired Patent token" is a relational join over the token pool, not a per-token predicate; explicitly deferred — likely a separate ADR (a token-join node) if a use case appears.
- **Observation-point ratification.** This ADR pins the pre-`evolve!` snapshot; confirm against a running build that no rate/predicate needs the post-`finish!` count instead (ADR 0006 open question), and add the round-trip + determinism test.

## Contract delta

New CONTRACT §9.5 (Predicate selection & state advance of agentic tokens): the `TokenPredicate{kind, clauses}` node and `PRED_OP_WHITELIST`; the structured-LHS `predicate` field on `ReactantSpec` (§6.5 addition); the `matches`-filter generalization of the bind path; the 𝓕ₜ-measurability rule and pinned observation point (resolving the ADR 0006 §9.4 `TokenAgg` open question); PHASE-AS-ATTRIBUTE as the canonical structured-token model (maintainer ruling), with the per-KIND style retained only as the degenerate `species`-field predicate; the `SetField{field,value}` action statement + `Field{name}` ExprNode leaf, `@move` redefined as `SetField{:species,…}` sugar; the `@select`/`@advance` DSL forms. Amendment to §8.3 (`validate` rule 7 — predicate + `SetField`). Amendment to ADR 0005 (the action family extends to `{SetSpecies,SetParams,SetField,Log,Seq}`; two new closed nodes `TokenPredicate`/`Field` + `PRED_OP_WHITELIST`, all round-trip-tested). Files: `src/solvers.jl` (bind filters `:202-211`, `:274-283`; the `@move`/`structured_rhs` write path `:333-397`, esp. `set_species!` `:378` generalized to `SetField`; `update_u_structured!` `:630-639` per-phase question), `src/state.jl` (`UnfoldedReactant` `:4-9`, `sample_transitions!` `:194-205`), `src/interface/reaction_parser.jl` (`:66-76`), `src/interface/create.jl` (`:300-301`), `src/interface/agents.jl` (`set_species!` `:71`, the protocol the field-write generalizes), a new predicate+field-write evaluator in `src/interface/queries.jl` (beside the ADR 0006 query API).
