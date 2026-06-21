# ADR 0008 — Agentic species under a filtration: predicate-based token selection (`TokenPredicate`, `@select`)

- Status: Proposed 2026-06-21. Generalizes a structured-token LHS reactant from "select by KIND only" to "select the 𝓕ₜ-measurable sub-population a predicate carves out," via one closed eval-free IR node (`TokenPredicate`) and one DSL form (`@select`). Answers the maintainer's Q4 directly: "a transition may want to take a project in a given state / subject to a filtration expression — more complex than a previous species." Unifies the two BD modeling styles (phase-as-species vs phase-as-attribute) by making `get_species(a)==kind` the degenerate predicate, so they stop being two mechanisms.
- Verification note: `file:line` read against `ref-agents` (Julia 1.12.5). The binding filter that this ADR generalizes is literally `filter(a -> get_species(a) == type && !isblocked(a), structured_token)` at the genesis pass (`src/solvers.jl:202-205`) and again at the ongoing pass (`src/solvers.jl:274-277`); both then `sort!` by `priority(a, transName)` descending (`:207-211`, `:279-283`) and take `allocs[j,i]` instances. The selector is KIND + not-blocked; there is no attribute filter anywhere on the bind path. The `priority` hook (`src/interface/agents.jl:67`) lets you ORDER by an attribute but never SELECT by one.
- Relates to: ADR 0006 (structured tokens: the token protocol `get_species`/`isblocked`/`priority`/`past_bonds`, `agents.jl:49-72`; the `(species, creation_index)` total order; the `TokenAgg{reducer,species,field,active_only}` READ node — this ADR adds the BIND analogue), ADR 0005 (closed `ExprNode` IR, `OP_WHITELIST`/`DIST_WHITELIST`, eval-free `validate`; `TokenPredicate` is a new closed node like `TokenAgg`), ADR 0002 (priority-weighted allocation — the predicate filters the candidate set the allocator then orders), ADR 0007 (the initial marking that creates the tokens a predicate selects over). CONTRACT §1 (modalities — `perstep` still illegal for structured species), §4 (D1/D4/D5 — the predicate is 𝓕ₜ-measurable, never reads the future or a later-in-tick RNG draw), §9 (the structured-token subsystem this extends). Becomes a CONTRACT §9.5 + an ADR 0005 amendment.

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

### (D) The phase-as-species vs phase-as-attribute fork — unified, not chosen

A structured token's lifecycle state ("which pipeline phase") can be modeled two ways:

- **Phase-as-species:** each phase is a distinct KIND (`:Discovery,:Phase1,…,:Market`, ADR 0006 north-star); advancing is `@move(:PhaseN,:PhaseN1)` (`solvers.jl:360-389`) re-using the same token object; selection is implicitly by kind.
- **Phase-as-attribute:** ONE `:Project` KIND with a `phase::Symbol` host-struct field; advancing is a `SetAttr(phase=:PhaseN1)` action statement (ADR 0005 action family); selection is `@select(Project, phase==:PhaseN)`.

`TokenPredicate` SUBSUMES both: `get_species(a)==kind` is the degenerate clause `(species, ==, kind)`. So phase-as-species = filtration on the `species` field; phase-as-attribute = filtration on a `phase` field. The engine bind path is one code path either way; the choice becomes a modeling preference, not two mechanisms.

**Recommendation (BD template default): phase-as-attribute.** One `:Project` kind whose `phase` advances by `SetAttr`, selected by `@select`. Rationale: (1) no kind explosion and continuous predicates (`npv > θ`) work; (2) `@move`'s redundant re-`entangle!` (`solvers.jl:376`, flagged for removal in ADR 0006 D) is avoided entirely — the token never changes container key; (3) attribute history is uniform in `past_bonds` (ADR 0006). Phase-as-species stays fully supported for models where phases really are different resource KINDS with different protocol overrides. The predicate node is what makes this a free choice.

### (E) DSL surface

LHS reactant with a filtration, compact:

```julia
# phase-as-attribute (recommended): consume one Phase-2 Project with NPV over a param θ
@ct(cycletime), @select(Project, phase == :Phase2 && npv_estimate > θ) + @rate(budget) --> @move_attr(phase, :Phase3)

# phase-as-species (ADR 0006 style) is the degenerate predicate, written as today:
@ct(cycletime), Project_Phase2 + @rate(budget) --> @move(:Project_Phase2, :Project_Phase3)
```

`@select(kind, clause && clause …)` parses to a `TokenPredicate`: the `&&` tree splits into `clauses`, each `field op value` mapped to a `Clause`, `value` lowered to an `ExprNode` (params/obs/time allowed, `Sample` rejected by `validate`). Stoichiometry `q * @select(...)` keeps its meaning (bind `q` matching tokens). The parser hook is `recursive_find_reactants!`'s structured/macrocall branch (`src/interface/reaction_parser.jl:66-76`, `src/interface/create.jl:300-301`), which already special-cases `@structured`/`@move`; `@select` joins that closed set.

### (F) `validate` additions (ADR 0005 amendment)

`validate(spec)` (CONTRACT §8.3) gains rule 7: for every `TokenPredicate` — `kind` is a declared `specStructured` species; each `Clause.op ∈ PRED_OP_WHITELIST`; each `Clause.value` is a whitelisted `ExprNode` containing NO `Sample` node (§C rule 1); each `Clause.field` is a known token field for that kind where statically knowable (protocol fields always; host-struct extras checked against a registry-declared field list if the kind registers one, else deferred to the `hasproperty` guard at runtime with a warning). `@select` on a non-structured species is a diagnostic (you cannot predicate-filter a `Float64` count). `perstep` modality on a predicated structured reactant stays illegal (CONTRACT §1.4 — unchanged).

## Consequences

- The maintainer's Q4 is answered with one node and one DSL form: a transition can take "a project in a given state subject to a filtration expression," eval-free and reproducible.
- BD models collapse from a lattice of phase-kinds to one `:Project` kind + predicates, and continuous selection criteria (`npv > θ`, `cost_to_date < budget_cap`) become expressible for the first time.
- It composes with ADR 0006: `TokenPredicate` (bind/consume) and `TokenAgg` (read/aggregate) share the field-access and observation-point discipline; a rate can READ "n active Phase-2 projects" (`TokenAgg`) and a transition can BIND "the highest-priority Phase-2 project with npv>θ" (`TokenPredicate`) consistently in one tick.
- The IR grows by exactly one closed node and one closed op-whitelist (`PRED_OP_WHITELIST`), both additive, JSON-Schema-enumerable, eval-free — same discipline as `OP_WHITELIST`/`TOKEN_REDUCER_WHITELIST`.
- Costs: the bind path gains a `matches` filter at four sites (`solvers.jl:202,274` filters; the predicate must thread to `UnfoldedReactant` via `sample_transitions!`, `state.jl:194-205`); `validate` rule 7; the `@select` parser branch; and the observation-point pin must be wired (resolving the ADR 0006 open question for `TokenAgg` at the same time). Disjunction, cross-token-joins ("a Project AND its paired Patent"), and arg-count-varying clauses are explicitly OUT of v1.

## North-star tie-in (BD)

One `ProjectToken <: AbstractStructuredToken` carrying `phase, npv_estimate, cost_to_date, pos_remaining, acquired, acq_time`. Pipeline advance: `@ct(ct_p2), @select(Project, phase==:Phase2) + 3*@conserved(scientist) + @rate(budget) --> SetAttr(phase=:Phase3)` on `Binomial(q, PoS)` success — no per-phase kinds. The acquisition lever (ADR 0007 §C / 0006) injects `ProjectToken(phase=:Phase2, acquired=true, npv_estimate=…)`; a `@select(Project, acquired==true && phase==:Phase3)` transition can then apply post-acquisition operational changes to ONLY the acquired late-stage programs — a filtration the engine cannot express today. rNPV reductions (`TokenAgg`/`sumattr` over `active(prob, :Project)` filtered by phase) read the same pinned snapshot the bind predicates use, so the pre/post-lever diff (same seed, §4 D5) is well-defined.

## Invariants

1. **Filtration measurability.** A `TokenPredicate.clauses[*].value` contains no `Sample` node and references only token fields, params, observables, `TimeRef`, and `TokenAgg` — all 𝓕ₜ-measurable at the allocation point (validate rule 7).
2. **Pinned observation point.** Token-reading subexpressions in a predicate (and in a `TokenAgg`) read the pre-`evolve!` reflected counts (`update_u_structured!`, `solvers.jl:649/653`); a predicate and a rate reading the same pool in one tick agree.
3. **Deterministic bind set.** For fixed `(model, seed)`, the bound token set per tick is reproducible: predicate narrows, `(priority desc, creation_index)` total order (ADR 0006 inv 5) selects within, integer take bounds.
4. **Degenerate compatibility.** `predicate == nothing` (or empty clauses) is byte-for-byte today's kind-only bind; every existing structured reaction is unaffected.
5. **Bind ≠ read.** `TokenPredicate` consumes/holds (CONTRACT §1 modality applies, `perstep` still illegal); `TokenAgg` reads non-consumingly. The two never alias: a `TokenAgg` in a rate never blocks a token, a `TokenPredicate` on an LHS always can.

## Open questions

- **Field-schema declaration.** Should a structured KIND declare its selectable field list (types) at registration (ADR 0006 C registry) so `validate` can statically type-check `Clause.value` against `Clause.field`, or is the runtime `hasproperty` guard + a soft warning enough? Static is safer for agentic authoring; it costs a richer registry entry.
- **Disjunction / nesting.** v1 is conjunctive. If real BD models need `(phase==:Phase2 || phase==:Phase3)`, decide between (a) multiple reactant rows, (b) an `in` clause `phase ∈ [:Phase2,:Phase3]` (already covered by `:in`!), or (c) a boolean-tree node — prefer (b)/(a) before (c) to keep the node closed.
- **Cross-token predicates.** "Bind a Project together with its paired Patent token" is a relational join over the token pool, not a per-token predicate; explicitly deferred — likely a separate ADR (a token-join node) if a use case appears.
- **Observation-point ratification.** This ADR pins the pre-`evolve!` snapshot; confirm against a running build that no rate/predicate needs the post-`finish!` count instead (ADR 0006 open question), and add the round-trip + determinism test.

## Contract delta

New CONTRACT §9.5 (Predicate selection of agentic tokens): the `TokenPredicate{kind, clauses}` node and `PRED_OP_WHITELIST`; the structured-LHS `predicate` field on `ReactantSpec` (§6.5 addition); the `matches`-filter generalization of the bind path; the 𝓕ₜ-measurability rule and pinned observation point (resolving the ADR 0006 §9.4 `TokenAgg` open question); the phase-as-attribute vs phase-as-species unification and BD-template recommendation; the `@select` DSL form. Amendment to §8.3 (`validate` rule 7). Amendment to ADR 0005 (one new closed node + one closed op-whitelist, round-trip-tested). Files: `src/solvers.jl` (bind filters `:202-211`, `:274-283`), `src/state.jl` (`UnfoldedReactant` `:4-9`, `sample_transitions!` `:194-205`), `src/interface/reaction_parser.jl` (`:66-76`), `src/interface/create.jl` (`:300-301`), a new predicate evaluator in `src/interface/queries.jl` (beside the ADR 0006 query API).
