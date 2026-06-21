# ADR 0007 — Interface & initial-state contract: lifecycle phases, declarative initial marking, and state dump/restore

- Status: Proposed 2026-06-21. Pins the callable API surface as a three-phase lifecycle (authoring → construction → live), adds a DECLARATIVE, SERIALIZABLE initial marking for structured-token populations (today they can only be built imperatively, which breaks the §8.2 S2 "run is determined by `(model.json, seed)`" guarantee), and specifies a state dump/restore (checkpoint) artifact whose schema is a superset of the initial marking. Answers the maintainer's Q3 ("do we have a contract for the interface methods — network definition, initial values/initial state, simulation, joins") and the follow-up requirement that the interface support instantiating structured agents as a list of structs and dumping the system state.
- Verification note: every `file:line` below was read against `ref-agents` (Julia 1.12.5). The two load-bearing gaps — the `"structured"` container is entangled EMPTY at construction (`src/solvers.jl:612`) with no declarative population path, and `_reinit!` (`src/solvers.jl:619-628`) restores neither the token population nor any creation counter — are present-source facts, not hypotheticals.
- Relates to: ADR 0001 (native `ReactionNetworkProblem` engine; `_step!`/`_reinit!`/`_projected_to` at `src/solvers.jl:642/619/675`), ADR 0003 (typed struct-of-columns IR; `specInitVal` is the plain-species initial value), ADR 0004 (append-only runtime mutation; the construction-frozen varmap is the authoring→live boundary), ADR 0005 (eval-free single-JSON + `ExprNode`; `(model.json, seed)` reproducibility, §8.2 S2, §8.5 outputs-are-separate), ADR 0006 (structured tokens: KINDS a priori, INSTANCES live; TYPES/CONSTRUCTORS are host Julia referenced by name through the per-network registry; the `(species, creation_index)` deterministic sort-key resolved 2026-06-20). CONTRACT §4 (D1/D5/D6/D7 determinism) and §5.4 (`specInitVal`) are the obligations this ADR threads the initial state through. Becomes CONTRACT §10.

## Context

CONTRACT §1–§9 pin the SEMANTICS of a network (modalities, time, tick order, determinism, object model, joins, structured tokens) but never state the CALLABLE INTERFACE as one contract: which methods exist, in which order they may be called, and what each requires/guarantees. The pieces are scattered — the append-only index discipline is §6.2/ADR 0004, the live-guard on `rem_parts!` is §7.5/J8, the mutation API is ADR 0004, the query API is §9.3 — but nobody has written the lifecycle state machine they collectively imply. Two concrete holes follow from that omission.

Ground-truth facts, verified against source:

- **The constructor is** `ReactionNetworkProblem(acs::ReactionNetworkSchema, u0=Dict(), p=Dict(); name="reaction_network", kwargs...)` (`src/solvers.jl:536-542`). It resolves time control (`get_tcontrol`, `:552`), captures `structured_token_names` (`:556-557`), compiles attribute closures against a varmap frozen here (ADR 0004; `src/compilers.jl:149`), fills `u0_init` from the `u0` override or `specInitVal` (`:561-569`), merges params (`:571-576`), and at the very end `entangle!(network, FreeAgent("structured"))` (`src/solvers.jl:612`) — an EMPTY structured container. There is no `seed`/`rng` keyword (CONTRACT §4 D6 requires one — unbuilt) and no initial-population keyword.
- **Plain-species initial value is declarative and serializable:** `specInitVal` (a column; `init_u!`, `src/state.jl:85-91`; constructor `:561-569`; `u0` keyword overrides per name). This IS part of `model.rdj.json` (CONTRACT §5.4, §6.4), so a plain run is reproducible from `(model.json, seed)` (§8.2 S2).
- **Structured-species initial population is NOT declarative.** A structured species' `state.u` column is a DERIVED reflected count (`update_u_structured!`, `src/solvers.jl:630-639`), not storage (ADR 0006 A). The only way to create instances is imperative host code: `add_structured_token!(problem, agent) = entangle!(getagent(problem,"structured"), agent)` (`src/interface/agents.jl:42-44`), called in a host `for`-loop (`tutorial/agents-integration/agents.jl:108-109`). So the initial token population lives in a host script, NOT in the model document — and §8.2 S2 ("nothing outside the model document and the seed may influence the trajectory") is VIOLATED for any structured model.
- **Simulation entry.** `simulate(a, max_t=Inf)` is the AlgebraicAgents driver (`AlgebraicAgents.jl/src/interface.jl:160`); it takes a max-TIME (not a step count) and RETURNS the agent — a verified gotcha worth pinning in the interface contract. `_step!` advances one tick (`src/solvers.jl:642-673`); `_projected_to` reports completion once `state.t > state.tspan[2]` (`:675-677`).
- **Re-init is incomplete.** `_reinit!` (`src/solvers.jl:619-628`) restores `u` (from `sol[1]`), `t`, `ongoing_transitions`, `log`, `observables`, `sol` — but NOT the RNG (CONTRACT §4 D7 gap), NOT the token population (the live `inners` of the `"structured"` container are untouched, so a second run sees the END-state token pool, not the initial one), and NOT any per-species creation counter (ADR 0006 resolved the sort-key to `(species, creation_index)`, but no counter exists yet).
- **There is no state-dump path at all.** `as_state(u, t, state)` (`src/state.jl:219-221`) `deepcopy`s the whole live object in-memory for the ensemble summary; the persisted artifacts are the `sol` DataFrame and the `log` (CONTRACT §8.5), neither of which captures the structured-token population, the RNG state, the creation counters, or the in-flight `ongoing_transitions`. You cannot today serialize, halt, and resume a run.

The maintainer's follow-up ("instantiate the structured agents as a list of structures and 'add' them, so the interface should support serialization to dump the state of the system") names both halves of the fix: a declarative LIST form for the initial population, and a state-dump artifact. The key realization is that **these are the same schema**: the initial marking is a zero-tick, no-ongoing-transitions checkpoint, and restoring a checkpoint is constructing with its population as the initial marking plus a clock/RNG/`u` overlay.

## Decision

### (A) The three-phase lifecycle is the interface contract

A model passes through exactly three phases; each method is legal in a defined subset, and the phase boundary is an OBSERVABLE state of the object, not a convention.

| Phase | Object | Mutability | Legal operations | Forbidden |
|---|---|---|---|---|
| **Authoring** | `ModelSpec` / `ReactionNetworkSchema` | freely mutable; reindex allowed | DSL `@ReactionNetworkSchema`; `from_json`; `@push`; `@join`/`union_acs!` (§7); `@equalize`/`equalize!` (§7.4, the `rem_parts!` reindexer `equalize.jl:52`); `refine`/`abstract` (ADR 0009); `validate(spec)` (§8.3) | stepping; live-token ops |
| **Construction** | `ReactionNetworkProblem(spec; seed, u0, registry, population)` | one-shot transition | freezes the varmap (`compilers.jl:149`), seeds the RNG (§4 D6), compiles closures, **instantiates the initial marking** (§B), entangles the `"structured"` container (`solvers.jl:612`) WITH that population | — (atomic) |
| **Live / stepping** | constructed `ReactionNetworkProblem` | APPEND-ONLY (ADR 0004) | `simulate`/`step!`; `add_species!`/`add_transition!`/`add_param!`/`activate!`/`deactivate!`; `add_structured_token!`/`disentangle!`; `apply_patch` (§8.6); the query API (§9.3); `dump_state`/`reinit!` | `@equalize`/any `rem_parts!` reindex (§7.5 J8 live-guard); reordering or deleting object rows |

**The §7.5 live-guard generalizes to a phase guard.** Any operation that reindexes object tables (the sole reindexer is `rem_parts!`, `equalize.jl:52`) MUST refuse once the object is in the Live phase. Construction is the arming point: a `live::Bool` (or "has a constructed problem been built from this spec") flag is checked by `equalize!`, `refine`, and any composition path. Authoring-time composition and refinement are unrestricted; live mutation is append-only and never reindexes (ADR 0004 INV-1/INV-2).

This table IS the answer to Q3: network definition (authoring), initial values/initial state (construction, §B/§C), simulation (live), and joins (authoring) are all placed, with their legal-call ordering made explicit.

### (B) A declarative, serializable INITIAL MARKING for structured tokens

Add a top-level `population[]` array to `model.rdj.json` (ADR 0005 document shape) — the structured analogue of `specInitVal` for plain species. It is INPUT, not output (§8.5), and it is what makes a structured run reproducible from `(model.json, seed)` (§8.2 S2).

```jsonc
"population": [
  { "species": "Project",            // FK → a structured Species.name (validate: must be specStructured)
    "kind":    "ProjectToken",        // registry key → host constructor (ADR 0006 B/C); default = species name
    "count":   12,                    // how many instances (or omit and give an explicit list, below)
    "attributes": {                   // ExprNode per field, sampled once at t=0 through state.rng (§4 D5)
        "phase":        { "Const": "Discovery" },
        "npv_estimate": { "Sample": { "dist": "LogNormal", "args": [ ... ] } },
        "budget":       { "Const": 0.0 } } }
]
```

Two authoring forms, one semantics:
- **Count + attribute distributions** (above): `count` instances of `kind`, each with `attributes` drawn independently through `state.rng`. Compact for "100 inbound deals with sampled NPV."
- **Explicit list of structs** (the maintainer's "list of structures"): `"instances": [ {attributes…}, {attributes…}, … ]` — one row per token, attributes as literal `Const`s or expressions. This is the host-side ergonomic form too: `population = [ProjectToken(:Phase2; npv=3.1), ProjectToken(:Discovery; npv=0.4), …]` passed to the constructor builds the same marking from already-constructed host structs (the registry is not even needed in this path — the host handed over live values).

**Instantiation contract (at construction, before `t=0`, before the first `_step!`):**
1. Iterate `population[]` in DECLARED ORDER; within a `count`-form entry, instantiate `k = 1..count` in order.
2. For each instance: resolve `kind` against the registry (ADR 0006 C) to a host constructor; evaluate each `attributes` `ExprNode` through `context_eval`/`state.rng` (§4 D5) — sampled attributes consume the seeded stream deterministically; `entangle!` the token into the `"structured"` container (`agents.jl:42-44`).
3. Assign the per-species CREATION INDEX: the k-th token instantiated of a given species gets `creation_index = k` (the 2026-06-20 ruling; this is the deterministic tie-break of the §9.2 invariant-5 sort-key `(species, creation_index)`). The counter is part of the seeded run state.
4. After the full marking is built, `update_u_structured!` (`solvers.jl:630-639`) reflects the live unblocked counts into the structured `state.u` columns — so the initial `u` of a structured species equals its initial active population, an assertable invariant.

`validate` (§8.3) gains rule 6: every `population[].species` is a declared species with `specStructured = true`; every `kind` is a registered name (ADR 0006 C dangling-reference diagnostic); every `attributes` field is an `ExprNode` over the closed whitelist; `count`/explicit-list integrality is checked. The marking is eval-free: it references kinds and fields by name, never carrying Julia source (ADR 0006 B trust boundary).

### (C) State dump / restore (checkpoint) — a superset of the initial marking

`dump_state(problem)::StateDump` serializes a live run at a TICK BOUNDARY into an eval-free, JSON-representable artifact; `restore(spec, dump; registry)::ReactionNetworkProblem` reconstructs it. The dump schema is the initial marking (§B) PLUS the dynamic run state:

```jsonc
{ "schema":  "rdj-checkpoint/1",
  "model_hash": "…",                 // matches the dump back to its ModelSpec (S2 / §8.5 run.json)
  "t":   42.0, "tick": 42,           // clock (§2.1)
  "rng": { "type": "Xoshiro", "state": [ … ] },   // §4 D6/D7 — the seeded stream's exact state
  "creation_counters": { "Project": 137, … },     // per-species, for the §9.2 sort-key tie-break
  "u":   [ … ],                      // plain-species columns only; structured columns are re-derived (ADR 0006 A)
  "population": [                     // EVERY live inner token, the §B schema extended with dynamic fields:
     { "species": "Project", "kind": "ProjectToken", "uuid": "…",   // identity (ADR 0006 inv 6)
       "creation_index": 7,
       "attributes": { "phase": "Phase3", "npv_estimate": 3.1, … }, // CURRENT field values (literals, not exprs)
       "bound_transition": "advance_p2" | null,                     // by transition id, not object ref
       "past_bonds": [ ["Phase2", 30.0, "advance_p1"], … ] } ],     // the rNPV/audit history (ADR 0006)
  "ongoing": [                        // in-flight Transition instances (state.jl:14-26)
     { "i": 3, "t": 38.0, "q": 4, "state": 2.5,
       "bound": ["uuid…", …], "nonblock": [ … ] } ],                // bound tokens by uuid
  "log_ref": "runs/<hash>/<seed>/ledger.arrow" }                    // §8.5 — ledger lives in Arrow, by reference
```

**Restore is construction with overlays:** `restore` runs the §B initial-marking machinery using `dump.population` (reconstructing each token's TYPE via the registry, ADR 0006 B/C — the dump carries field VALUES and a kind NAME, never Julia source, so it inherits the §8.4 no-eval guarantee), then overlays `t`, `u` (plain columns), the RNG state, the creation counters, and rebuilds `ongoing` (re-binding bound tokens by uuid). The structured `state.u` columns are NOT restored from the dump — they are re-derived by `update_u_structured!` from the restored population, keeping ADR 0006 invariant 3 (`length(active(prob,s)) == state.u[idx(s)]`) true by construction. `bound_transition`/`bound`/`nonblock` are serialized as ids/uuids and re-linked after both tables exist (the live `Transition` field holds an object ref, `state.jl:13`/`:19-20`, which is not directly serializable — the id/uuid indirection is the eval-free bridge).

**Two artifacts, one schema, distinguished by a flag:** a checkpoint with `t==tspan[1] && ongoing==[]` IS an initial marking (the `population[]` of a `model.rdj.json`); a checkpoint at an arbitrary tick is the full resume artifact. `restore(spec, dump)` with a zero-tick dump is identical to `construct(spec; population=dump.population)`. This is why §B and §C are one decision, not two.

### (D) `reinit!` is completed to honor (B) and (C)

`_reinit!` (`solvers.jl:619-628`) MUST additionally: (1) restore the RNG to the seed-implied initial state (§4 D7 — currently missing); (2) tear down the live token population and rebuild the §B initial marking (currently the END-state tokens survive into the next run); (3) reset every per-species creation counter to the count established by the initial marking. After this, `init → step* → reinit! → step*` reproduces the first trajectory exactly for a structured model, closing D7 for structured runs.

### (E) The canonical interface signatures (the Q3 contract surface)

- `validate(spec; registry_names)::Vector{Diagnostic}` — authoring, eval-free (§8.3, +rule 6 for `population`).
- `from_json(io; registry)::ModelSpec` / `to_json(spec)::String` — authoring (ADR 0005, §8.1–8.2).
- `@join` / `union_acs!` / `@equalize` / `equalize!` — authoring (§7); reindexers refuse when live (§A).
- `refine(spec, transition, sub; ports)` / `abstract(...)` — authoring (ADR 0009).
- `ReactionNetworkProblem(spec; seed::Union{Integer,Nothing}=nothing, u0=Dict(), registry=Dict(), population=spec.population)::ReactionNetworkProblem` — construction; seeds the RNG (§4 D6), instantiates the marking (§B), arms the live-guard (§A). NEW kwargs vs `solvers.jl:536`: `seed`, `registry`, `population`.
- `simulate(problem, max_t=Inf)::problem` — live; max-TIME not step-count, returns the agent (the pinned gotcha).
- `add_structured_token!` / `add_species!` / `add_transition!` / `add_param!` / `activate!` / `deactivate!` / `apply_patch` — live, append-only (ADR 0004, §8.6).
- `tokens`/`active`/`blocked`/`tokenattr`/`ntokens`/`nactive`/`sumattr`/`maxby`/`nbound` — live, deterministic query API (§9.3).
- `dump_state(problem)::StateDump` / `restore(spec, dump; registry)::ReactionNetworkProblem` — live (§C).
- `reinit!(problem)` — live; restores the §B marking + RNG + counters (§D).

## Consequences

- A structured run becomes reproducible from `(model.rdj.json, seed)` exactly as a plain run is — the §8.2 S2 guarantee, currently false for structured models, is restored, because the initial token population is now IN the document (§B) rather than in a host script.
- Halt/resume and "fork a run at a tick to A/B a lever" become first-class: the acquisition-lever counterfactual (ADR 0006 north-star) can `dump_state` immediately before the lever, then `restore` twice and apply the lever to one copy — the cleanest way to isolate the lever's effect under one seed (§4 D5).
- The dump is eval-free and inherits the ADR 0005 RCE-closure: it carries token field VALUES + a kind NAME resolved through the host registry, never Julia source (ADR 0006 B). A malicious checkpoint can no more execute code than a malicious model.
- Costs: real engine work — `seed`/`registry`/`population` constructor kwargs, the §B instantiation pass, per-species creation counters threaded into run state and `_reinit!`, the `dump_state`/`restore` pair with id/uuid relinking for `bound_transition`/`ongoing`, and the §A live-guard generalization. The `ongoing` portion of a mid-run dump is the heaviest piece (it must re-link bound tokens and re-create frozen `trans` snapshots); a zero-tick initial marking (the common case) avoids all of it. None of this disturbs CONTRACT §1–§9 semantics; it is additive.

## North-star tie-in (BD)

The pipeline is seeded declaratively: `population[] = [{species: Project, kind: ProjectToken, count: 30, attributes:{phase: Const(:Discovery), npv_estimate: Sample{LogNormal,…}}}]` — 30 discovery-stage programs with sampled NPV, reproducible from the document + seed. The acquisition counterfactual is `dump_state` at the lever tick → `restore` twice → `add_structured_token!(probA, ProjectToken(:Phase2; acquired=true, …))` on one copy only → step both to horizon → diff the rNPV query-reductions (ADR 0006). The `past_bonds` realized-phase history rides along in both the marking and the dump, so audit survives a halt/resume.

## Invariants

1. **Marking-determinism.** Given `(spec, seed)`, the initial marking — instance count, per-instance sampled attributes, and creation indices — is reproducible (it consumes the seeded stream in declared order, §4 D1/D5).
2. **Reflected-count consistency at t=0.** After the marking is built, `length(active(prob, s)) == state.u[idx(s)]` for every structured species `s` (re-derived by `update_u_structured!`, never set directly).
3. **Dump round-trip.** `restore(spec, dump_state(p))` yields a problem that steps identically to `p` from that point (same RNG state, same population, same ongoing set, same creation counters) — the resume analogue of §8.2 S2.
4. **Eval-free state I/O.** Neither `population[]` nor a `StateDump` is ever `Meta.parse`d or `eval`d on load; both reference kinds/functions by name through the host registry (ADR 0006 B, §8.4 S4).
5. **Phase legality.** A reindexing op (`rem_parts!`) called on a Live object errors; an append-only op called in any phase is safe (ADR 0004 INV-1/2; §7.5 J8 generalized).
6. **`reinit!` completeness.** `init → step* → reinit! → step*` reproduces the first trajectory for BOTH plain and structured models (the §D fix closes §4 D7 for structured runs).

## Open questions

- **`ongoing` snapshot fidelity.** An in-flight `Transition` carries a frozen sampled-attribute snapshot `trans::Dict{Symbol,Any}` (`state.jl:17`) realized from `context_eval` at spawn. Serializing it eval-free requires that snapshot to be a dict of LITERALS (it is, post-evaluation) — confirm no closure/`Function` leaks into it; if one can, mid-run dump must forbid it or the offending attribute must be value-frozen.
- **Model-version binding.** A `StateDump` carries `model_hash`; restoring against a DIFFERENT (e.g. patched, §8.6) spec is ill-defined. Decide whether restore requires an exact hash match or permits an append-only-superset spec (the latter would let a dump from before a mutation be resumed under the post-mutation model — useful, but needs a rule).
- **Checkpoint cadence / size.** Per-tick dumps are O(#live-tokens) each; for long BD runs decide a cadence policy and whether `past_bonds` (unbounded, ADR 0006 retired-token growth) is dumped in full or summarized. No pre-optimization until a real run shows it bites (mirrors the ADR 0006 retired-token ruling).

## Contract delta

New CONTRACT §10 (Interface & initial-state contract): §10.1 the three-phase lifecycle + legal-op table (A); §10.2 the constructor contract and canonical signatures (E); §10.3 the declarative initial marking `population[]`, both authoring forms, and the construction-time instantiation contract (B); §10.4 simulation interface (the `simulate` max-time/returns-agent pin) and the completed `reinit!` (D); §10.5 state dump/restore, the checkpoint schema, and the "initial marking = zero-tick checkpoint" identity (C); §10.6 invariants 1–6. Additions to §8.3 (`validate` rule 6, `population` well-formedness) and §8.5 (the `StateDump` is an OUTPUT-class artifact alongside `sol`/`log`, persisted by reference to the Arrow ledger). Files: `src/solvers.jl` (constructor `:536`, `_reinit!` `:619`, `update_u_structured!` `:630`, structured container `:612`), `src/state.jl` (`as_state` `:219`, run-state fields `:40-63`), `src/interface/agents.jl` (`add_structured_token!` `:42`), a new `src/interface/checkpoint.jl` (`dump_state`/`restore`).
