# ADR 0011 — Action callbacks: population-level writes + the general-code escape hatch (`Invoke`), eval-free

- Status: Accepted — maintainer-ratified, 2026-06-21. The maintainer asks that actions/callbacks (transition pre/post-actions and the ADR-0010 Rule actions) be able to (a) WRITE token fields over a selected population, not just the one bound token, and (b) execute GENERAL queries/code, not only the closed declarative verb set. This ADR pins both while preserving the load-bearing eval-free / no-RCE guarantee (CONTRACT §8.4 S4) that ADRs 0005/0006 established. It resolves MVP finding **I** (a Rule cannot write a token field over a `@select`-ed set) and answers the maintainer's "actions should run general queries/code" directive.
- Verification note: all cited engine facts checked against `ref-agents`, Julia 1.12.5. Load-bearing facts this builds on are real and already used by accepted ADRs: the per-network registry indirection lowering a name to `:(registry[fn](...))` (ADR 0006 §C, replacing `@register`'s `@eval`); `context_eval(state, transition, …)` as the single value-eval chokepoint (`src/state.jl:67-71`); the `@move`/`set_species!` token field-write path (`src/solvers.jl:360-389`, `:378`); the deterministic token query API + `(species, creation_index)` total order (ADR 0006 §E, CONTRACT §9.2/§9.3); the append-only live mutation API (ADR 0004).
- Relates to: ADR 0010 (defines the canonical action *type family* in its §C and the Rule channel — this ADR extends that family by two members and is co-located with it), ADR 0008 (defines `SetField` = single-bound-token field write + `TokenPredicate`/`@select` selection + the 𝓕ₜ-measurability rule; `SetTokens` is the population-level generalization of `SetField` over a `TokenPredicate`), ADR 0006 (the per-network function registry + the host-Julia-vs-serializable-data trust boundary — `Invoke` is that exact boundary lifted from value-position to statement-position), ADR 0005 (the eval-free ExprNode/action IR + `validate` this amends), ADR 0004 (append-only mutation — the obligation an `Invoke` body must honor), CONTRACT §4 (determinism — the obligation an `Invoke` body must honor), §9.3/§9.5 (the query/selection surface an action body composes). Maps to a CONTRACT §12.3 amendment + an ADR 0005 amendment.

## Context

After ADR 0010 the action *type family* is `{SetSpecies, SetParams, SetField, AddToken, Activate, Deactivate, Log, Seq}` (canonical definition in ADR 0010 §C). Two expressivity limits remain, both raised by the maintainer:

1. **Single-token vs population write.** `SetField{field, value}` (ADR 0008 §D) writes the field of the FIRING INSTANCE's bound token(s), so it is meaningful only in a transition post-action and only over the tokens that transition just bound. A great many BD decisions are population-level and have no firing transition: "on a competitor's readout, write down EVERY active Phase-2 oncology program's `pos_remaining` by 10%"; "when cash is short, deprioritise all Discovery-phase programs." These are Rule actions over a `@select`-ed SET, which the family cannot express (MVP finding I).

2. **Closed verbs vs general code.** The declarative verbs cover create / write-field / write-pool / write-param / toggle-line / log. Real callbacks sometimes need arbitrary logic — a multi-step query, a bespoke reallocation, a custom valuation update, a conditional branch over the token population — that no fixed verb set will ever fully cover. The maintainer wants actions to run general queries/code.

The hard constraint both must respect: CONTRACT §8.4 S4 — **no field of a `model.rdj.json` is EVER `Meta.parse`d or `eval`d on load** (the RCE that ADR 0005 closed at `loadsave.jl:65,72` and the `Base.convert` hooks at `ReactiveDynamics.jl:99,101,102`). "General code in the file" is exactly what S4 forbids. The resolution is the mechanism ADR 0006 §C ALREADY built for value-helpers and never reopened: code lives in the HOST package (compiled by Julia normally), and the inert JSON references it BY NAME against a host-populated registry. This ADR lifts that boundary from value-position (a registered helper inside a `Call`/`ExprNode`) to statement-position (a registered side-effecting callback as an action). The general-code-RETURNING-A-VALUE case already works today (ADR 0006 §C); what is new here is general-code-AS-A-STATEMENT.

## Decision

Extend the canonical action type family by exactly two members — one declarative, one the escape hatch — keeping the JSON inert and the trust boundary identical to ADR 0006.

### (A) `SetTokens{predicate, assigns}` — declarative population-level field write

```julia
struct SetTokens <: ActionStmt
  predicate :: TokenPredicate                 # ADR 0008 §C; e.g. @select(Project, phase==:Phase2 && area==:onc)
  assigns   :: Vector{Tuple{Symbol,ExprNode}} # field => value, evaluated per matched token
end
```

- **Semantics.** Select every token matching `predicate` (the ADR-0008 `TokenPredicate` selection, reused verbatim), iterate them in the §9.2 `(species, creation_index)` TOTAL ORDER (never raw AA `Dict` order — D4), and for each, write every `(field, value)` pair. `value` is evaluated through `context_eval` in the context of THAT token (a `Field{name}` leaf reads the token's own current field, ADR 0008), so `pos_remaining => Field(:pos_remaining) * 0.9` writes each program down by 10%.
- **Relationship to `SetField`.** `SetField` is the degenerate single-token case (the firing instance's bound token); `SetTokens` is the population case keyed by a predicate. The two are the per-bound-token and per-selected-set forms of the same field-write; `SetField` stays for the post-action advance idiom (`@advance`), `SetTokens` is the new general population write.
- **Context.** Legal in BOTH a Rule action (it carries its own predicate, needs no bound token — this is the finding-I fix) AND a transition pre/post-action. Default DSL: `@set(Project, phase==:Phase2 && area==:onc, pos_remaining *= 0.9)`.
- **Determinism.** A `SetTokens.assigns` `value` MAY contain a `Sample` (a write may legitimately draw, like `SetField`, ADR 0008 §F), consuming the seeded stream in the fixed §9.2 token order, so D1 holds. Selection is 𝓕ₜ-measurable (the predicate is, ADR 0008).
- **Eval-free.** `predicate`/`assigns` are closed `TokenPredicate`/`ExprNode` trees; `validate` checks kind is `specStructured`, fields exist, ops/dists whitelisted (extends ADR 0008 `validate` rule 7). Never an eval.

### (B) `Invoke{fn, args}` — the general-code escape hatch via the registry

```julia
struct Invoke <: ActionStmt
  fn   :: Symbol            # a key into the per-network registry (ADR 0006 §C)
  args :: Vector{ExprNode}  # evaluated and passed positionally
end
```

- **Lowering.** `to_expr(Invoke{fn,args})` → `:(registry[$fn](state, transition, <lowered args>))` — EXACTLY the ADR-0006 §C registry-by-name lowering, except the result is a STATEMENT (return value discarded) rather than a value sub-expression. `transition` is the firing instance, or `nothing` for a Rule action. The function is captured at build so the ADR-0006 §C compile-once property holds; no `invokelatest`, no world-age bump.
- **What the body may do.** `registry[fn]` is an ORDINARY Julia function the host wrote and compiled in its own package. It may run ANY logic: compose the §9.3 query API (`tokens`/`active`/`sumattr`/`@select`), read/write `state.u`/`state.p`, call the append-only live mutation API (`add_structured_token!`/`add_transition!`/`activate!`/…), draw from `state.rng`, branch, loop. This is fully general code.
- **Eval-free / no-RCE — UNCHANGED (S4 holds verbatim).** The JSON carries only the NAME `fn` (a string) and `ExprNode` args. On load, `validate` checks `fn ∈ keys(registry)` and the registered calling-convention/arity (ADR 0006 §C tag); the lowering produces `registry[fn](...)`. NO bytes from the file are `Meta.parse`d or `eval`d — identical to how a registered value-helper `α(...)` already lowers (ADR 0006 §C). The trust boundary is the host program that populated the registry, NOT the file; the file alone stays inert. A malicious model can name `fn` but cannot supply its body, and an unregistered name is a `validate()` diagnostic, never code execution.
- **The honest trade (what `validate` CANNOT guarantee).** A declarative verb is statically checkable: `validate` proves it respects append-only, determinism, and ranges. An `Invoke` body is TRUSTED-BUT-UNVERIFIED: `validate` proves only that the name resolves and the arity matches — it cannot prove the body honors the contract invariants. Therefore the contract places explicit OBLIGATIONS on the host author of an `Invoke` callback (enforced by review/convention, not by `validate`):
  - **O1 (determinism, §4 D5):** any randomness MUST come from `state.rng`; never bare `rand()`/global RNG, never wall-clock/`Date`/`Math.random`, no external I/O except `log`. Else D1 (reproducibility) breaks.
  - **O2 (append-only, ADR 0004):** structural mutation MUST go through the append-only live API; never `rem_parts!`/reorder/reindex on a live model (it would invalidate the frozen varmap).
  - **O3 (tick-boundary, ADR 0010 §12.5):** the body runs at the action's fixed `_step!` slot (post-action in `finish!`/`evolve!`, or step 10 for a Rule); it MUST NOT itself advance the clock or re-enter `_step!`.
  - **O4 (purity):** the body MUST be a function of `(state, transition, state.rng)` only — no captured mutable global state — so a re-run from the same seed reproduces it.
- **Guidance.** PREFER the declarative verbs (validatable, reproducible-by-construction, LLM-authorable, JSON-Schema-enumerable). Reach for `Invoke` ONLY when the logic genuinely exceeds the closed verbs. An `Invoke` model is no longer fully self-validating — it is as trustworthy as the host package, exactly like an ADR-0006 registered behavior function.

### (C) The reconciled canonical action family

ADR 0010 §C remains the single canonical definition; it grows by these two members:

```
{ SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq }
ACTION_VERBS = (:set_species, :set_params, :set_field, :set_tokens,
                :add_token, :activate, :deactivate, :invoke, :log, :seq)
```

Context legality (extends the ADR-0010 §C restriction):

| Verb | Rule action | Transition pre/post-action | Note |
|---|---|---|---|
| `SetSpecies`/`SetParams` | ✓ | ✓ | pool/param writes |
| `AddToken` | ✓ | ✓ (structured RHS) | create a token |
| `Activate`/`Deactivate` | ✓ | ✓ | toggle a line (ADR 0004) |
| `SetField` | ✗ (no bound token) | ✓ | write the firing instance's bound token (ADR 0008) |
| `SetTokens` | ✓ | ✓ | write a `@select`-ed population (this ADR §A) |
| `Invoke` | ✓ | ✓ | general host callback (this ADR §B) |
| `Log`/`Seq` | ✓ | ✓ | log / compose |

`validate` rejects `SetField` in a Rule (no bound token, ADR 0010 §C), an unregistered `Invoke.fn`/`AddToken.kind`, an arity mismatch, an unknown `SetTokens` field/kind, an out-of-set verb, or a `Ref` to an undeclared name — all diagnostics, never eval (§8.4 S4).

## Consequences

- **Finding I resolved.** A Rule can now write a token field over a selected population via `SetTokens` — "write down all Phase-2 oncology `pos_remaining` by 10% on a competitor readout" is one declarative, validatable, reproducible action.
- **General code supported, eval-free.** `Invoke` gives actions arbitrary host logic through the ADR-0006 registry-by-name boundary; the JSON stays inert and S4 holds verbatim. The value-position counterpart (a registered helper in an `ExprNode`) already existed (ADR 0006 §C), so general code is now expressible in BOTH value and statement positions.
- **A two-tier expressivity spectrum, explicit.** Tier 1 = declarative verbs (statically safe, LLM-authorable, the default). Tier 2 = `Invoke` (trusted host code, the escape hatch, obligations O1–O4 on the author). The contract documents which guarantees survive at each tier rather than pretending one size fits all.
- **No new whitelist axis.** `SetTokens` reuses `TokenPredicate`/`Field` (ADR 0008); `Invoke` reuses the ADR-0006 registry key set. The IR grows by two action members only — additive, JSON-Schema-enumerable like the rest.
- **Costs / engine work.** `SetTokens` reuses the ADR-0008 predicate evaluator + the `set_species!`-style field write generalized over a selected set. `Invoke` reuses the ADR-0006 §C registry lowering, retargeted to statement position with the return discarded; the calling-convention tag must mark `Invoke` callbacks as `(state, transition)`-convention. `validate` gains the `SetTokens` and `Invoke` checks. All land with the ADR-0010 action infrastructure.

## North-star tie-in (BD)

```jsonc
// SetTokens: a competitor's Phase-3 readout writes down our same-indication programs (a Rule)
{ "id":"competitor_readout", "fire_mode":"once",
  "guard": { "node":"call","op":">",
             "args":[ {"node":"timeref"}, {"node":"ref","kind":"param","name":"readout_t"} ] },
  "action": { "node":"set_tokens",
              "predicate": { "kind":"Project",
                             "clauses":[ ["phase","==","Phase2"], ["area","==","onc"] ] },
              "assigns":[ ["pos_remaining",
                           {"node":"call","op":"*","args":[ {"node":"field","name":"pos_remaining"},
                                                            {"node":"const","value":0.9} ]}] ] } }
```

```jsonc
// Invoke: a bespoke portfolio-rebalancing callback the host wrote and registered (a Rule)
{ "id":"rebalance", "fire_mode":"every_tick",
  "guard": { "node":"call","op":">=","args":[ {"node":"ref","kind":"species","name":"cash"},
                                              {"node":"const","value":0} ] },
  "action": { "node":"invoke", "fn":"rebalance_portfolio",
              "args":[ {"node":"ref","kind":"param","name":"risk_budget"} ] } }
```
with the host supplying `registry = Dict(:rebalance_portfolio => (state, transition, risk_budget) -> begin … end, …)` — ordinary Julia that may query `active(state, :Project)`, draw via `state.rng`, and mutate append-only. The file references the name; the code is host-compiled.

## Open questions

- **`SetTokens` write order vs allocation.** A `SetTokens` in a transition post-action writes during `finish!`; a Rule `SetTokens` writes at step 10. Both are after the tick's allocation, so they affect the NEXT tick's selection/genesis (consistent with ADR 0010 §12.5). Confirm no use case needs an intra-tick re-selection.
- **`Invoke` reentrancy / nesting.** Should an `Invoke` body be permitted to call `apply_patch` (ADR 0004) or trigger other rules? Current ruling: it may use the append-only API directly but MUST NOT re-enter `_step!` (O3). Confirm whether a guarded sub-action API is wanted later.
- **Static affordances for `Invoke`.** Optional future nicety: let the host DECLARE an `Invoke` callback's effect footprint (which species/params/kinds it writes) so `validate` can at least check the footprint against append-only, recovering partial static safety without constraining the body.

## Contract delta

- **CONTRACT §12.3** (canonical action family, ADR 0010 §C): extend to `{SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}`; add the `SetTokens` population-write semantics, the `Invoke` general-code escape hatch + obligations O1–O4, and the updated context-legality table.
- **CONTRACT §8.3** `validate`: add the `SetTokens` (predicate/field/kind) and `Invoke` (registered name + arity) checks; restate that neither is ever eval'd (§8.4 S4 unchanged).
- **CONTRACT §8.4 S4**: note explicitly that `Invoke`/registered statement callbacks are the ADR-0006 §C boundary in statement position — name-only in the file, body compiled in the host — so S4 holds verbatim; add the trusted-but-unverified caveat and the prefer-declarative guidance.
- **CONTRACT §9.5 / ADR 0008**: `SetTokens` is the population generalization of `SetField` over a `TokenPredicate`; cross-reference.
- **MVP finding I**: RESOLVED by `SetTokens`; the general-code ask resolved by `Invoke`.

Files relevant to implementation: `/Users/bima/ReactiveDynamics-review/src/solvers.jl` (`finish!` post-action + `@move` write path `:360-389`, Rule slot `event_action!` `:316-326`), `/Users/bima/ReactiveDynamics-review/src/state.jl` (`context_eval` `:67-71`, `set_params` `:246`), `/Users/bima/ReactiveDynamics-review/src/interface/queries.jl` (the new predicate/field evaluator the `SetTokens` write reuses, ADR 0008), the ADR-0006 registry plumbing on `ReactionNetworkProblem` (`src/state.jl:40-63`), and the ADR-0005 IR module (action lowering + `validate`).
