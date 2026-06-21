# ADR 0009 — Hierarchical refinement, open-port composition, and compact process authoring

- Status: Proposed 2026-06-21. Adds the REFINEMENT axis the contract is otherwise silent on: (1) granularity/substitution — replace a coarse transition with a finer sub-model that is plug-compatible at its boundary, and the inverse abstraction; (2) compositionality — compose fragments by connecting declared open PORTS rather than remembering which names to `@equalize`; (3) compactness — `@pipeline` sugar and reusable parameterized `@process` modules so a business process is authored in a line, not a screenful. Answers the maintainer's framing requirement: "a modeling framework more suitable at modeling business processes — compact/expressive definition, compositionality, and various levels of granularity with more refined dynamics possibly substituted."
- Verification note: `file:line` read against `ref-agents` (Julia 1.12.5). The decisive enabling fact: composition today is name-keyed string surgery (`union_acs!`/`equalize!` over `add_part!`/`incident` loops with `recursively_substitute_vars!`, `operators/joins.jl:10-59`, `operators/equalize.jl:24-66`, `compilers.jl:28-44`), NOT a categorical pushout. CONTRACT §7.4/J7 already established that the ADR-0003 promoted `ReactantSpec` incidence table turns species identification into integer-FK repointing instead of Expr rewriting. Refinement reuses that exact FK-splice — so the headline feature costs almost no new mechanism.
- Relates to: ADR 0003 (the promoted `ReactantSpec` incidence table — the FK-repoint substrate refinement splices on; "structure is implicit, not morphic" is the gap this closes), CONTRACT §7 (composition semantics J1–J9, C1–C6; namespacing `normalize_name` `joins.jl:100`; the dead `prepend_obs` `joins.jl:87-97`; the undefined `include_model` bug J9 `joins.jl:226`), ADR 0004 (refinement is AUTHORING-time, pre-construction — it reindexes, so the §7.5/J8 live-guard forbids it on a stepping model), ADR 0005 (a refinement/module is itself a serializable `ModelSpec` fragment), ADR 0007 (the three-phase lifecycle places refine/compose in the authoring phase; ports are declared in the spec). CONTRACT §6 (object model — ports are a thin annotation on Species; no new table). Becomes CONTRACT §11 + a §7 addendum.

## Context

The signed-off contract pins a FLAT network: one bag of species/transitions, composed by flat name-merge (§7). It says nothing about the VERTICAL axis the maintainer needs — building a process coarsely, then substituting a finer sub-process for one step without disturbing the rest. For a BD/R&D pipeline this is the core workflow: model the portfolio coarsely (`Discovery → … → Market`, one transition per phase) for fast what-if, then ZOOM into the Phase-2 bottleneck (sub-divide into screening / lead-opt / tox / regulatory sub-steps with their own resources) when you need fidelity there — and have the rest of the portfolio model not notice.

Three deficiencies block this today:

- **No substitution.** There is no `refine`/`abstract`. A transition is atomic; to "expand" it you hand-edit the model and re-wire every reactant by name.
- **Flat, fragile composition.** `@join` (`joins.jl:201-238`) namespaces every species `m__X` and merges by name; SHARING is opt-in via `@equalize` equation blocks (§7.4) — but you must REMEMBER which names to identify, there is no declared boundary, and the identification is `recursively_substitute_vars!` string surgery (`equalize.jl:55-63`) that can corrupt a species name colliding inside a subexpression (the §7.4/J7 hazard, recorded in ADR 0003). Also `:E`/`:obs` are silently dropped on join (§7/J4) and the file-include branch calls undefined `include_model` (§7/J9).
- **Verbose authoring.** One reaction per line. A 5-phase pipeline with per-phase cycletime/PoS/resources is ~5 hand-written reactions with repeated `@move` boilerplate; there is no pipeline or module abstraction.

The enabling insight is already in the contract: §7.4/J7 proves that once reactants are first-class `ReactantSpec` FK rows (ADR 0003 Phase 2), species identification becomes "repoint the integer `species` FK from the deleted index to the survivor, drop the row" — a structural O(rows) operation, no Expr rewriting. Refinement is the SAME operation applied at a transition boundary. So the big feature is cheap; the work is mostly defining the boundary and the sugar.

## Decision

### (A) Open ports: a transition/sub-model has a declared boundary

A **port** is a boundary species through which a fragment connects to its environment. Add a thin closed-tag annotation to the Species record (CONTRACT §6.4 — NOT a new table):

```
role::PortRole ∈ {private, input, output, shared}     # default: private
```

- `private` — internal to the fragment; auto-namespaced on compose (today's default `m__X`).
- `input` / `output` — an OPEN port: the fragment expects to be connected here (input = consumed-from boundary, output = produced-into boundary). The directionality is advisory/validation-only; it documents intent and lets `validate` warn on a dangling open port.
- `shared` — a global species identified by bare name across all fragments WITHOUT namespacing (the existing `@catchall` semantics, `joins.jl:117`, made a first-class role).

A **transition's boundary** is just its LHS/RHS `ReactantSpec` rows: LHS species are its input ports, RHS species its output ports. No new structure — the incidence table already IS the boundary. This is what makes (B) free.

### (B) Refinement = boundary-matched splice via FK-repoint (the §7.4/J7 operation)

```julia
refine(spec, transition, submodel; ports = Dict(:in_species => :sub_input, :out_species => :sub_output))
```

A **refinement** of a coarse transition `T` is a sub-model `S` plus a boundary map identifying `S`'s open `input` ports with `T`'s LHS species and `S`'s open `output` ports with `T`'s RHS species. `refine` splices `S` into the parent in four structural moves, all authoring-time:

1. **Namespace** `S`'s `private` species (`normalize_name`, `joins.jl:100`); leave its `input`/`output`/`shared` ports un-prefixed for matching.
2. **Identify ports with boundary species** — for each `ports` entry, identify `S`'s port species with the parent's boundary species by the §7.4/J7 FK-repoint: repoint every `ReactantSpec.species` FK that pointed at the port to the boundary species index, drop the port's `:S` row. NO `recursively_substitute_vars!`, so no collision-corruption.
3. **Append** `S`'s transitions and remaining (private + newly-identified) species/params/obs/events as new rows (append-only, §6.2; this also fixes §7/J4 by merging `:E`/`:obs` uniformly).
4. **Remove the coarse transition** `T` — drop its `:T` row and its `ReactantSpec` rows. (This is a reindex → authoring-time only, §7.5/J8 live-guard; forbidden on a stepping model.)

Because steps 2 leaves the BOUNDARY species unchanged (same indices, same name, same cost/valuation), the coarse and refined models are **plug-compatible**: any other transition referencing the boundary species is untouched, and observables/ledger over the boundary stay meaningful. The portfolio doesn't notice that Phase-2 became four sub-steps.

`abstract(spec, transitions, into; ports)` is the inverse: collapse a connected sub-graph back to one coarse transition, its boundary = the sub-graph's open ports, its `cycletime`/`pos`/`cost` summarized (§D advisory check).

### (C) Boundary consistency — an advisory, not a proof

Refinement does NOT claim the fine model is behaviorally equivalent to the coarse one — that would need a bisimulation the framework can't check. Instead, `validate` (run authoring-time) emits ADVISORY diagnostics comparing the coarse transition's attributes to aggregates of the refinement, where computable:

- `coarse.cycletime ≈ Σ path cycletimes` along the refinement's critical path (warn if off by a tolerance);
- `coarse.prob_of_success ≈ Π sub-PoS` along the success path;
- `coarse.cost ≈ Σ sub-costs`;
- port-balance: every `input` port is consumed by ≥1 sub-transition LHS, every `output` port produced by ≥1 sub-transition RHS (an unconsumed input or unproduced output is a likely modeling error).

These are warnings the author can override (the refinement may legitimately change dynamics — that's the point of zooming in). They make the granularity ladder auditable without overclaiming.

### (D) Compactness — `@pipeline` sugar and `@process` modules

**Pipeline sugar** for the dominant business-process shape (a chain of phases):

```julia
@pipeline Project begin
    Discovery => Phase1 : (ct = 1.0, pos = 0.4, res = 2*@conserved(scientist) + @rate(budget))
    Phase1    => Phase2 : (ct = 2.0, pos = 0.6, res = 3*@conserved(scientist) + @rate(budget))
    Phase2    => Phase3 : (ct = 3.0, pos = 0.5, res = 5*@conserved(scientist))
    Phase3    => Market : (ct = 1.0, pos = 0.9)
end
```

expands to N `flow`-genesis transitions (CONTRACT §2.8 — routing, requires the upstream phase as an upfront-consumed LHS) with the `@select`/`@advance(phase,…)` idiom (ADR 0008 §D — phase-as-attribute is canonical: one `Project` kind, `phase` is a field, advance writes it via `SetField`), each carrying its `(ct, pos, res)`. The whole BD north-star pipeline in one block instead of five hand-wired reactions with repeated boilerplate.

**Process modules** — a named, parameterized fragment instantiated multiple times:

```julia
@process phase_gate(in, out; ct, pos, res) = @ReactionNetworkSchema begin
    @ct(ct), @select($in_species) + $res --> @advance($in_species, $out_species), prob => $pos
end
# instantiate + compose by ports:
gate_a = phase_gate(:Phase1, :Phase2; ct=2.0, pos=0.6, res = 3*@conserved(scientist))
gate_b = phase_gate(:Phase2, :Phase3; ct=3.0, pos=0.5, res = 5*@conserved(scientist))
model  = @compose gate_a gate_b   # connect by matching input/output port names (§E)
```

A module is just a `ModelSpec` fragment with declared ports (A); instantiation substitutes its parameters (eval-free — the params are `ExprNode`/value substitutions into the fragment, not code-gen); reuse + compactness ride the existing namespacing substrate.

### (E) `@compose` — port-connected composition (the §7 join, with a boundary)

`@compose f1 f2 …` is `@join` (§7) PLUS automatic port matching: `output` ports of one fragment are identified with same-named `input` ports of another by the §7.4/J7 FK-repoint (not string surgery), `private` species are namespaced, `shared` species are identified by bare name. It is the explicit-boundary form of `@join`; `@join`/`@equalize` remain available as the manual, no-declared-ports path. `@compose` also CLOSES the §7/J4 and J9 bugs en route: it merges `:E`/`:obs` uniformly (J4) and never takes the undefined `include_model` file branch (J9) — composition is over already-parsed `ModelSpec`s.

### (F) Entity-level refinement (a structured token hosting a sub-network) — scoped out of v1

AlgebraicAgents permits a structured token to itself be a container of inner agents, so a `ProjectToken` could host its OWN sub-network (per-project internal dynamics) — a second, orthogonal granularity locus (entity zoom vs process zoom). This is genuinely useful for "model each acquired program's internal R&D in detail" but is a larger design (per-token sub-engines, their clocks, their coupling to the parent ledger). DEFERRED to a future ADR; v1 refinement is process-structural (B)–(E) only. Noting it so the boundary annotation (A) is designed not to preclude it.

## Consequences

- The granularity/substitution requirement is met by reusing the §7.4/J7 FK-splice: `refine`/`abstract` are structural, collision-safe, and plug-compatible at the boundary — the "coarse portfolio, zoom the bottleneck" workflow works without disturbing the rest of the model.
- Compositionality stops being "remember which names to `@equalize`": fragments declare `input`/`output`/`shared` ports and `@compose` connects them, with the FK-repoint making identification exact.
- Compactness: a pipeline is a block, a reusable process is a parameterized module — the BD north-star is authored in a few lines.
- Three pre-existing §7 defects are closed as a side effect: `:E`/`:obs` merge on compose (J4), the dead `prepend_obs` gains a live caller (the observable-namespacing path), and `@compose` never reaches the undefined `include_model` (J9).
- All of this is AUTHORING-time and additive: it produces a plain `ModelSpec` that constructs/serializes/simulates exactly as a hand-written flat model (no runtime, determinism, or serialization change). The refinement record itself round-trips as ordinary spec structure (ADR 0005).
- Costs: `refine`/`abstract`/`@compose` over the promoted `ReactantSpec` table (depends on ADR 0003 Phase 2 having landed — refinement is gated on the incidence promotion); the `role` annotation + port-matching in compose; the `@pipeline`/`@process` macros; the §C advisory aggregate checks in `validate`. The FK-repoint itself is shared with `equalize` (§7.4/J7), so the core operation is written once. Entity-level refinement (F) and refinement of a STRUCTURED transition (whose LHS selects tokens by predicate, ADR 0008) need care — the sub-model's boundary then includes a token predicate — and are flagged below.

## North-star tie-in (BD)

Coarse portfolio: `@pipeline Project (Discovery⇒Phase1⇒Phase2⇒Phase3⇒Market)` with summary `(ct, pos, cost)` per edge — fast rNPV what-if across the whole book. Then `refine(spec, :Phase2_advance, phase2_detail; ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out))` substitutes a detailed Phase-2 sub-model (screening → lead-opt → tox → filing, each with scientist/budget draws and its own PoS) — the boundary species `:Phase2`/`:Phase3` are unchanged, so the Discovery→Phase1 and Phase3→Market transitions, the acquisition lever (ADR 0007/0006), and the rNPV ledger all keep working untouched. The §C advisory check flags if the detailed sub-model's aggregate PoS drifts from the coarse `0.6` you'd been using — surfacing exactly the fidelity gain (or modeling error) the zoom was meant to expose.

## Invariants

1. **Plug-compatibility.** After `refine(spec, T, S; ports)`, every transition NOT in `{T} ∪ S` is structurally unchanged: same species indices for boundary species, same `ReactantSpec` rows. Refinement is local to the boundary.
2. **FK-exactness.** Port identification and species collapse repoint integer `ReactantSpec.species` FKs (§7.4/J7), never `recursively_substitute_vars!`; a species name colliding inside a subexpression cannot corrupt incidence.
3. **Authoring-only.** `refine`/`abstract`/`@compose`/`@equalize` reindex (drop the coarse `:T`/port `:S` rows) and so are FORBIDDEN on a Live object (§7.5/J8, ADR 0007 §A phase guard); the live mutation API (ADR 0004) is the only runtime change path and never reindexes.
4. **Closure.** `refine`/`abstract`/`@compose` map `ModelSpec`(s) to a `ModelSpec`; the result is a valid input to a further refine/compose/construct (§7 C1 closure, extended to the vertical axis).
5. **Round-trip.** A refined/composed model serializes and re-loads as a flat `model.rdj.json` identical to a hand-authored equivalent (ADR 0005 §8.2 S1); refinement leaves no runtime trace.
6. **Advisory, not equivalence.** The §C consistency diagnostics are WARNINGS; refinement does not assert behavioral equivalence between coarse and fine, and the engine never auto-refines or auto-abstracts.

## Open questions

- **Refining a structured transition.** When the coarse transition selects tokens by an ADR-0008 `TokenPredicate`, the sub-model's input port is a predicated structured species; define how the predicate flows to the sub-model's entry transition (likely: the entry sub-transition inherits the boundary predicate). Needs the ADR-0008 node to exist first.
- **Multi-port / fan-out boundaries.** A coarse transition with multiple LHS or RHS species maps to multiple input/output ports; confirm the `ports` map handles many-to-one and one-to-many boundary identifications, and how `abstract` infers the boundary of a sub-graph with internal branching.
- **Aggregate-check computability.** The §C critical-path cycletime / product-PoS aggregates are well-defined for a linear chain; for a sub-model with branches/loops, decide what "the" aggregate is (or restrict the check to acyclic refinements and warn otherwise).
- **`@process` parameter hygiene.** Parameter substitution into a fragment must not capture or collide with the fragment's private species names; confirm the namespacing order (substitute params, THEN namespace privates) and that a param value referencing a port resolves correctly.
- **Entity-level refinement (F).** Whether/when to let a structured token host a sub-network — a separate ADR; ensure the `role` annotation (A) and the boundary concept don't preclude it.

## Contract delta

New CONTRACT §11 (Refinement & open-port composition): §11.1 ports as a Species `role` annotation (A) — addendum to §6.4; §11.2 a transition's boundary = its `ReactantSpec` LHS/RHS rows; §11.3 `refine`/`abstract` as the §7.4/J7 FK-splice (B), authoring-time, with the four-move algorithm; §11.4 the §C advisory consistency diagnostics (validate addition); §11.5 `@compose` port-matching + the J4/J9 closures; §11.6 `@pipeline`/`@process` compact authoring (D); §11.7 invariants 1–6. §7 addendum: `@compose` is the declared-boundary form of `@join`; `refine`/`abstract` extend the composition algebra to the vertical axis (C1 closure). All gated on ADR 0003 Phase-2 reactant promotion. Files: `src/operators/joins.jl` (`union_acs!` `:10`, `normalize_name` `:100`, dead `prepend_obs` `:87`, `@join` `:201`, `include_model` bug `:226`), `src/operators/equalize.jl` (the FK-repoint home, `:24-66`), `src/compilers.jl` (`recursively_substitute_vars!` `:28` — the string-surgery path refinement replaces), a new `src/operators/refine.jl` and the `@pipeline`/`@process` macros in `src/interface/create.jl`.
