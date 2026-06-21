# ADR 0010 — Rules, triggers & conditional transitions: the endogenous decision channel

- Status: Accepted — maintainer-ratified, 2026-06-21. The Business-Development demo ([../MVP_BD_DEMO.md](../MVP_BD_DEMO.md)) requires that DECISIONS — acquire a target, inject capital, toggle a pipeline line, raise a success probability — can be expressed *inside the model* as rules that fire on conditions, including state-contingent conditions ("acquire *if* portfolio value drops below X"), not just as exogenous host interventions at prescribed times. This ADR pins that channel. It is the resolution of MVP finding **B** (no endogenous, state-contingent decision channel) and it repairs the no-op event channel (CONTRACT §3.4 Invariant 7). It extends the Phase-0.5 modeling-language increment (ADRs 0007 interface/initial-state, 0008 token filtration, 0009 refinement/composition) and maps to CONTRACT §12.
- Verification note: all cited engine facts were checked against source on `ref-agents`, Julia 1.12.5. The three load-bearing defects this ADR builds on are real: `event_action!` fetches `:eventAction` but never evaluates it (`src/solvers.jl:323`); the `transActivated` skip-gate exists and is honored (`src/state.jl:179`); `periodic(state, 0.0)` returns `true` (one-shot at t0) and `periodic(state, p)` fires on period boundaries (`src/state.jl:239-244`).
- Relates to: ADR 0001 (native engine, single clock, `_step!` at `src/solvers.jl:642`), ADR 0002 (allocation — guards gate a transition *before* it competes for resources), ADR 0003 (typed IR — Rule and `guard` are new typed records), ADR 0004 (runtime mutation — `AddToken`/`Activate`/`Deactivate` are append-only/index-safe; `fire_mode: once` latch is run-state reset by `_reinit!`), ADR 0005 (eval-free ExprNode IR + the action-statement set this ADR extends), ADR 0006 (`AddToken` lowers to `add_structured_token!`; the acquisition lever), ADR 0007 (interface/initial-state — a `once`-rule latch joins the run state that `reinit!` resets §10.4 and `dump_state` checkpoints §10.5), ADR 0008 (token filtration — a guard/`TokenAgg` reads tokens under the §9.5 measurable-observation-point ruling), ADR 0009 (refinement — `@pipeline`/`@process` fragments carry rules and guarded transitions through compose). Answers maintainer rulings of 2026-06-21 (three forks below).

## Context

ReactiveDynamics models business/R&D processes. A real BD analysis is not a fixed schedule of events — it is a set of *decisions* taken in response to *conditions*: "if cash reserve covers the Phase-3 cost, start the Phase-3 program"; "when an attractive target appears (a proxy threshold is crossed), acquire it and fold in its scientists"; "retire the discovery line once the portfolio is late-stage-heavy." The Phase-0 contract could express *scheduled* interventions (a host-applied `apply_patch` at a tick boundary, §8.6; or the `@deterministic` + `periodic` calendar idiom, §2.8) but had **no in-model, state-contingent decision channel**. Three facts in the current engine made this gap concrete:

1. **The event channel is a no-op.** `event_action!` (`src/solvers.jl:316-326`) iterates `:E`, computes a fire-count `q` from the trigger (`q = v ? 1 : 0` for Bool, `rand(Poisson(v))` for numeric, `:321`), then in the inner loop merely *fetches* `state[i, :eventAction]` (`:323`) — it never evaluates the action. Events do nothing today (CONTRACT §3.4 Inv 7). So the one object that looks like a callback is dead.

2. **There is a latching activation gate but no stateless predicate.** `sample_transitions!` skips a transition whose `transActivated[i]` is `false` (`src/state.jl:179`). `transActivated` is *latching* state (flipped by ADR-0004 `activate!`/`deactivate!`); there is no way to say "fire this transition this tick *iff* a condition holds," re-evaluated every tick.

3. **The boolean/comparison vocabulary already exists, eval-free.** `OP_WHITELIST` already contains `> < >= <= == && || !` (ADR 0005:55). A guard like `cash < threshold && @t() > 90` is therefore *already* a legal, JSON-Schema-enumerable, eval-free `ExprNode` — nothing new is needed to *express* a condition. What is missing is (a) a place to attach a guard, (b) a working action channel, and (c) the action verbs a BD decision needs.

The unifying insight: **a "rule/callback" and a "conditional transition" are the same primitive** — a guard `ExprNode` evaluated per tick, driving a closed action set. They differ only in what the guard gates: a Rule gates an *action block*; a conditional transition gates its own *firing*. Both are deterministic, eval-free, and serializable, so both fit the existing contract with a small, closed extension.

## Decision

Add an **endogenous decision channel** with three additive, eval-free, JSON-Schema-enumerable parts. Nothing in the runtime hot path (compiled closures, allocator, clock) changes shape; the extension is at the IR/authoring boundary plus three localized engine edits.

### (A) A Rule is the repaired Event

The `:E` object becomes a **Rule**: a `(guard, action, fire_mode)` triple evaluated once per tick at a fixed point in `_step!`.

```
Rule {
  id        : Symbol            # closed identity, unique (CONTRACT §6.3)
  guard     : ExprNode          # CONTRACT §1/ADR-0005 IR; resolves to Bool or numeric
  action    : ActionStmt        # the extended closed set, (C) below
  fire_mode : {every_tick, once}   # default every_tick
  enabled   : Bool              # run-state latch; (D)
}
```

- **Guard evaluation.** `guard` is a standard `ExprNode` lowered by `to_expr` and evaluated via the RNG-threaded `context_eval(state, nothing, …)` (CONTRACT §4 D5). A `Bool` guard fires the action once; a numeric guard `v` fires it `rand(state.rng, Poisson(v))` times (preserving the existing `:321` multiplicity semantics, now correctly seeded). A guard MAY read species (`Ref{species}`), params (`Ref{param}`), observables (`Ref{obs}`), and time (`TimeRef`); it consumes RNG iff it contains a `Sample` node.
- **The repair.** `event_action!`'s dead `:323` line — a bare fetch of `:eventAction` — becomes `context_eval(state, nothing, state.wrap_fun(action_i))` inside the `1:q` loop. This is the minimal correctness fix that makes the whole channel live.
- **Determinism / ordering.** Rules iterate `parts(state, :E)` in index order (CONTRACT §4 D4 — already index-ordered, no `Dict` iteration). The fire point is fixed in the `_step!` order (see "Step placement" below), so `(model, seed)` fully determines rule effects (§4 D1).

### (B) A conditional transition carries a stateless `guard`

Add one field to the Transition record (CONTRACT §6.6):

```
guard : ExprNode   # default Const(true); resolves to Bool
```

Evaluated each tick in `sample_transitions!` and AND-ed with the existing latching gate:

```julia
# src/state.jl:179, generalized
fires_this_tick = state.transition_recipes[:transActivated][i] &&
                  context_eval(state, nothing, guard_i)::Bool
fires_this_tick || continue
```

- `guard` is **stateless** (re-evaluated every tick from current state); `transActivated` is **latching** (survives across ticks, flipped only by `activate!`/`deactivate!` or a Rule action). They compose: a line is live iff it is activated AND its guard holds this tick. This is the maintainer-ratified "guard field + rule toggling (both)" fork.
- A `guard` resolving to non-`Bool` is a `validate()` diagnostic (CONTRACT §8.3). A guard consuming RNG (a `Sample` node) is legal but discouraged for transition guards (it makes firing stochastic on top of the genesis Poisson draw); `validate` SHOULD warn.
- Cost: one extra `ExprNode` eval per active transition per tick. Negligible (`Const(true)` short-circuits at lowering — `to_expr(Const(true))` is the literal `true`, compiled away).

### (C) The closed action set, extended to the full endogenous set

The ADR-0005 action-statement set `{SetSpecies, SetParams, Log, Seq}` is extended to the maintainer-ratified **full endogenous set**:

```julia
abstract type ActionStmt end
struct SetSpecies <: ActionStmt; name::Symbol; value::ExprNode; mode::Symbol end  # mode ∈ {set, inc}
struct SetParams  <: ActionStmt; assigns::Vector{Tuple{Symbol,ExprNode}} end
struct SetField   <: ActionStmt; field::Symbol; value::ExprNode end               # write a BOUND token's field (ADR 0008 §D)
struct AddToken   <: ActionStmt; kind::Symbol; fields::Vector{Tuple{Symbol,ExprNode}} end
struct Activate   <: ActionStmt; transition::Symbol end
struct Deactivate <: ActionStmt; transition::Symbol end
struct Log        <: ActionStmt; msg::String end
struct Seq        <: ActionStmt; stmts::Vector{ActionStmt} end

const ACTION_VERBS = (:set_species, :set_params, :set_field, :add_token, :activate, :deactivate, :log, :seq)
```

`SetField` is contributed by [ADR 0008 §D](0008-token-filtration.md) under the maintainer's phase-is-an-attribute ruling: advancing a token's lifecycle state (`phase`) is a write to a token FIELD, and `@move(:from,:to)` is redefined as its degenerate `SetField{:species, Const(:to)}` sugar. It is the per-instance-token counterpart of `SetSpecies` (which writes a plain `state.u` column) and the mutate-existing counterpart of `AddToken` (which creates one). The above struct list is the single canonical definition of the action *type family* `{SetSpecies, SetParams, SetField, AddToken, Activate, Deactivate, Log, Seq}`; ADR 0008 references it rather than declaring a competing set.

**Context restriction (which verbs are legal where).** `SetField` writes the FIRING TRANSITION INSTANCE's bound token(s) (ADR 0008 §D, applied in `finish!` over `Transition.bound_structured_agents`), so it is meaningful ONLY in a transition **post-action** — a standalone Rule (this ADR, fired at step 10 with no associated transition instance and no bound token) has no referent for it. Therefore: a **Rule action** uses `{SetSpecies, SetParams, AddToken, Activate, Deactivate, Log, Seq}` (the system-/pool-level verbs); a **transition pre/post-action** additionally has `SetField` (and `AddToken` for structured RHS). `validate` rejects a `SetField` in a Rule action as a diagnostic. The acquisition lever (a Rule) uses `AddToken` to CREATE a project token; a phase-advance transition uses `SetField`/`@advance` to MUTATE the one it bound — complementary, never the same site.

Lowering (all eval-free, all into the closed allow-list — no `Expr`-head smuggling):

| Stmt | Lowers to | Safety basis |
|---|---|---|
| `SetSpecies{n,v,set}` | `state.u[idx(n)] = eval(v)` | plain-species write; idx frozen at construction (ADR 0004) |
| `SetSpecies{n,v,inc}` | `state.u[idx(n)] += eval(v)` | cash injection / capital lever |
| `SetParams{[(p,v)…]}` | `set_params(state, (p,eval(v))…)` (`src/state.jl:246`) | reserved verb already compiled (`compilers.jl:67`); synergy toggles |
| `SetField{field,v}` | `for a in trans.bound_structured_agents: setproperty!(a, field, eval(v))` (reusing `@move`'s `set_species!` site `solvers.jl:378`) | transition POST-ACTION only; writes the bound token, no `state.u` reindex (ADR 0008 §D); phase-advance = `SetField{:phase,…}` |
| `AddToken{kind,fields}` | `add_structured_token!(state, registry[kind](eval(fields)…))` | `entangle!`, append-only/index-safe (ADR 0006 §9.1); the acquisition |
| `Activate{t}` | `activate!(state, t)` → `transActivated[idx(t)] = true` | ADR-0004 soft gate, no reindex |
| `Deactivate{t}` | `deactivate!(state, t)` → `transActivated[idx(t)] = false` | ADR-0004 soft-deactivate; in-flight instances still finish (`finish!` ignores `transActivated`, `solvers.jl:406`) |
| `Log{msg}` | `log(state, msg)` (`src/state.jl:236`) | reserved verb |
| `Seq{[…]}` | run in order | composition; the acquisition is `Seq[AddToken, SetParams, SetSpecies]` |

`AddToken.kind` is a registry key (ADR 0006 §C); an unregistered kind, an unknown `Activate/Deactivate` target, an unknown `SetField.field`, a `SetField` in a Rule action (vs a transition post-action — the context restriction above), an out-of-set verb, or a `Ref` to an undeclared name is a `validate()` diagnostic — **never an eval**. The JSON stays inert data; the only Julia source produced is `to_expr`/the lowering table over a closed allow-list extended by the registry key set, exactly the ADR-0005/0006 discipline.

### (D) `fire_mode` and the `once` latch

A Rule with `fire_mode = once` fires at most once per run: when its guard first holds and the action runs, the engine sets `enabled = false`, and the rule is skipped thereafter. `every_tick` (the default) re-evaluates the guard every tick (subsuming the legacy periodic/continuous behavior; a periodic financing rule is `every_tick` with a `periodic(state, p)` guard — `src/state.jl:239`).

- `enabled` is **run-state**, part of the seeded state, and `_reinit!` MUST reset every `once` rule's `enabled` back to `true` (CONTRACT §4 D7 — re-run from the same seed reproduces the run, including the lever). This mirrors the §4 D7 obligation already placed on the RNG and the ADR-0006 token creation-counter.
- This is the maintainer-ratified "add `fire_mode {every_tick, once}`" fork. It makes the acquisition lever ("acquire the first time the proxy threshold is crossed") a one-liner instead of requiring a hand-rolled self-latch param.

### Step placement (determinism-critical)

The Rule fire point and the transition-guard evaluation point are fixed in the normative `_step!` order (CONTRACT §3.3):

- **Transition guards** are evaluated inside `sample_transitions!` (§3.3 step 5), the same pass that already reads `transActivated` and realizes per-tick attributes — so a guard sees this tick's freshly-sampled values. A guard gating firing means the transition contributes **no genesis proposal** this tick (it is skipped before `evolve!`), so it never competes for resources (clean interaction with ADR 0002).
- **Rules** fire at §3.3 step 10 (the former event slot, `solvers.jl:657`), i.e. AFTER `finish!` and the structured-count sync, BEFORE the ledger row (step 11) and the clock advance (step 12). Consequence: a Rule's `SetSpecies`/`SetParams`/`AddToken` takes effect on the SAME tick's valuation ledger row and is visible to the NEXT tick's guards and genesis. A token injected by `AddToken` at step 10 is reflected by the step-10-adjacent `update_u_structured!` and becomes bindable from the next tick. This placement is normative; moving it changes results.

## Consequences

- **Finding B resolved.** The acquisition lever, capital injections, synergy toggles, and pipeline on/off control are now expressible *inside* `model.rdj.json` as rules — scheduled (`@t()`-guarded) or state-contingent (state/param/observable-guarded), serializable, reproducible, and LLM-authorable. The BD demo's endogenous-decision requirement is met.
- **Inv 7 (event no-op) is fixed** as a side effect — the channel becomes live with the `:323` repair.
- **The IR grows by closed, enumerable increments only:** one Transition field (`guard`), one object reshape (`:E` → Rule with `fire_mode`/`enabled`), four action verbs (`SetField` from ADR 0008 §D, `AddToken`, `Activate`, `Deactivate`), and a `SetSpecies.mode ∈ {set, inc}` tag. All eval-free, all JSON-Schema-derivable like `OP_WHITELIST`. No new whitelist *axis* beyond the ADR-0006 registry (which `AddToken.kind` reuses).
- **No RCE reopening.** Guards and actions are typed nodes over closed allow-lists; `AddToken.kind` is a registry key (host-populated, ADR 0006). The file stays inert.
- **Determinism preserved.** Guards/rules consume RNG only via explicit `Sample` nodes through `state.rng` (§4 D5); fire points and iteration order are fixed (§4 D4); `once` latches reset on `_reinit!` (§4 D7). `(model, seed)` still fully determines the trajectory (§4 D1).
- **Costs / real engine work.** Three localized edits: (1) repair `event_action!` (`solvers.jl:316-326`) to evaluate the action, RNG-threaded; (2) generalize the `sample_transitions!` gate (`state.jl:179`) to AND in the guard; (3) implement the lowering for `AddToken`/`Activate`/`Deactivate` and the `once`-latch reset in `_reinit!`. The `AddToken` lowering shares code with the ADR-0006 `Construct{kind,args}` node and MUST land with the ADR-0006 token-path bug fixes (`solvers.jl:288,452,455,476,512`). `validate` gains the guard-type, action-verb, and `Activate/Deactivate`-target checks.
- **Authoring guidance.** A condition that must *persist once tripped* (hysteresis / a one-way switch) is a Rule that `Deactivate`s itself or sets a latch param, OR a `once` rule — NOT a bare stateless transition `guard` (which un-fires the moment the condition lapses). The contract (§12) states this so authors pick the right primitive.

## North-star tie-in (BD)

The expressive acquisition lever, entirely in-model:

```jsonc
// rule: acquire the first time the pipeline looks under-valued, fold in the target's capabilities
{ "id": "acquisition_lever", "fire_mode": "once",
  "guard": { "node":"call", "op":"&&", "args": [
              { "node":"call","op":">",  "args":[ {"node":"timeref"}, {"node":"const","value":90} ] },
              { "node":"call","op":"<",  "args":[ {"node":"ref","kind":"obs","name":"portfolio_proxy"},
                                                   {"node":"ref","kind":"param","name":"acq_trigger"} ] } ] },
  "action": { "node":"seq", "stmts": [
              { "node":"add_token", "kind":"ProjectToken",
                "fields":[ ["phase",{"node":"const","value":"Phase2"}], ["acquired",{"node":"const","value":true}],
                           ["acq_time",{"node":"timeref"}], ["npv_peak",{"node":"const","value":1200.0}] ] },
              { "node":"set_species","name":"scientists","mode":"inc","value":{"node":"const","value":30} },
              { "node":"set_params","assigns":[ ["synergy_pos",{"node":"const","value":1}],
                                                 ["synergy_eff",{"node":"const","value":1}] ] } ] } }
```

```jsonc
// conditional transition: only start Phase-3 programs when cash covers the Phase-3 cost
{ "id":"t_p2_to_p3", "genesis":"flow",
  "guard": { "node":"call","op":">=","args":[ {"node":"ref","kind":"species","name":"cash"},
                                               {"node":"ref","kind":"param","name":"phase3_cost"} ] },
  /* …rate, cycletime, prob_of_success, reactants… */ }
```

A scheduled-only deal is the same rule with a pure `@t()` guard (or `fire_mode: once` + a `periodic(state,0.0)` t0 guard). Either way the lever is in the reproducible model document, and pre/post-acquisition rNPV is the ensemble-averaged diff of two same-spec runs differing only by the scenario params that arm the rule (MVP §4.1).

## Open questions

- **Token-field-write action verb (MVP finding I — likely a pre-Phase-1 amendment).** The §9.5/ADR-0008 maintainer ruling makes a token's lifecycle state a FIELD (`phase`, `pos_remaining`, …) advanced by `@advance`/`SetField` as a transition-RHS idiom. But the §12.3 action set `{SetSpecies, SetParams, AddToken, Activate, Deactivate, Log, Seq}` writes only plain species and params — there is NO Rule/action verb that writes a token field, nor one that writes it over a `@select`-ed set (e.g. "on a competitor readout, write down all Phase-2 oncology programs' `pos_remaining` by 10%"). This is a genuine expressivity gap for endogenous decisions that perturb the POPULATION, not just pools/params. Proposed fix: add a closed `SetField{predicate::TokenPredicate, field, value::ExprNode}` (or `SetTokens`) action verb to §12.3, eval-free, lowering through the §9.5 selection + the §9.5 field-write — an ADR-0010/ADR-0008 amendment to settle before Phase-1.
- **Per-entity RNG substreams (MVP finding A).** Deferred: the demo uses ensemble-averaged Δ on a single RNG stream. If a tight per-seed paired Δ is later wanted, substreams keyed by `(root_seed, entity_id)` are the fix and would be a CONTRACT §4 extension (now recorded as CONTRACT §4.6 D-future). The §10.5 `dump_state`/`restore` fork-at-tick checkpoint is the partial mitigation available without the extension.
- **Guard observation point for token-aggregate guards. — now PINNED by ADR 0008 / CONTRACT §9.5.** A guard (or a transition guard) that reads sibling token counts via `TokenAgg` (ADR 0006 §9.4) observes the pre-`evolve!` reflected counts (`update_u_structured!`, `solvers.jl:649/653`), the same measurable observation point §9.5 fixes for `TokenPredicate` and `TokenAgg`. A guard `TokenAgg`/`TokenPredicate.value` MUST be 𝓕ₜ-measurable — no `Sample` node — exactly as §9.5 requires, which is consistent with this ADR's "guards SHOULD NOT consume RNG" guidance (B). So token-aware guards are deterministic and the original concern is closed; MVP finding C is correspondingly downgraded.
- **Action ordering within a tick across multiple rules.** Multiple `every_tick` rules firing on the same tick run in `:E` index order (deterministic). If two rules write the same species/param, last-in-index wins; `validate` MAY warn on write-write conflicts. Confirm this is the desired precedence or whether explicit priority is needed.
- **Per-tick `on_step` action hook (MVP finding D).** Out of scope here, but the action infrastructure built for rules is the natural home if per-tick accrual onto in-flight tokens is later required.

## Contract delta

- New **CONTRACT §12 (Rules, Triggers & Conditional Transitions)** — the contract-level specification of (A)-(D), the step placement, and the determinism obligations. (§10 = ADR 0007 interface/initial-state, §11 = ADR 0009 refinement, §9.5 = ADR 0008 filtration are the sibling Phase-0.5 sections; this is §12.)
- **§3.3** step 10: the event slot is renamed the *rule* slot; its placement (after `finish!`/sync, before ledger/clock) is pinned as normative.
- **§3.4 Invariant 7** updated: the event channel is repaired; the invariant now reads "a Rule whose guard holds fires its action `q` times per tick at step 10," with the `:323` no-op recorded as the corrected defect.
- **§5.3** gains a `guard` row for transitions (TVE? = yes, range = Bool, default `Const(true)`); the §5 TVE column split (finding D-bis) is flagged.
- **§6.6** Transition record gains `guard::ExprNode`; **§6.8** Event record is reshaped to the Rule record (`guard`/`action`/`fire_mode`/`enabled`).
- **§8 (ADR 0005)** amendment: `ACTION_VERBS` extended to the full set; `validate` rules gain guard-type, action-verb, and `Activate/Deactivate`-target checks; the JSON Schema re-derived.
- **§4 D7** obligation extended: `_reinit!` resets every `once` rule's `enabled` latch (joining the RNG state, per-species creation counters, and the §10.4/§10.5 ADR-0007 run-state items that `reinit!`/`dump_state` already govern).
- **§9.5 / §11 interaction:** a transition `guard` composes with the §9.5 `TokenPredicate` LHS selection (a guard gates whether the transition fires at all; a predicate selects which tokens it binds once it does) and travels through §11 `refine`/`@compose` as an ordinary Transition field.

Files relevant to implementation: `/Users/bima/ReactiveDynamics-review/src/solvers.jl` (`event_action!` `:316-326`, `_step!` `:642-673`), `/Users/bima/ReactiveDynamics-review/src/state.jl` (`sample_transitions!` `:174-217`, `set_params` `:246`, `periodic` `:239`, `log` `:236`), `/Users/bima/ReactiveDynamics-review/src/interface/update.jl` (`activate!`/`deactivate!`), `/Users/bima/ReactiveDynamics-review/src/compilers.jl` (action lowering), and the ADR-0005 IR module.
