# ReactiveDynamics.jl — Remaining-Work Handoff Plan

> **Purpose.** A self-contained, agentic-handoff plan for the work still open on the `ref-agents` revision as of 2026-07-09. It assumes NO prior conversation context. Read this top to bottom, then verify current state (§0) before touching code — every `file:line` and test tally below drifts and must be re-checked. Source of truth for decisions: `docs/adr/`; normative spec: `docs/CONTRACT_DRAFT.md`; the current-source map: `INVENTORY.md` (note its module table is stale — see WS-4). The Phase-0 sign-off and the already-landed Phase-0.6 increment are recorded in `docs/PHASE0_REVIEW.md`.

## 0. Orient & verify before acting

The revision is tracked as a contract (§1–§15), 14 ADRs (0001–0014), and a phased roadmap. As of this writing the modeling + analysis surface is BUILT and green — suite reported at **500 pass / 7 broken / 507**. Re-establish that baseline first:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'                 # expect ~500 pass / 7 broken
```

The 7 "broken" are deliberate `@test_broken` pins on known engine bugs (see WS-3). AlgebraicAgents is pinned to `Merck/AlgebraicAgents.jl@main` via a `[sources]` block in `Project.toml`; the Manifest records the exact tree-sha. Tests live in `test/semantic/*.jl` under a two-tier scheme: **T1-characterization** (runs against the current engine, pins bugs with `@test_broken`) and **T2-acceptance** (`@test_skip`-wrapped target behavior naming not-yet-built APIs). The acceptance ritual for every task below is: flip the relevant `@test_skip`/`@test_broken` → `@test` as the piece lands, and keep the suite green-when-expected.

Dev-loop gotchas carried from prior work (re-verify): a `dev/` Revise-based test server avoids recompilation; ADDING A NEW `src/` FILE requires clearing the compiled cache (Revise won't pick it up otherwise); when creating git worktrees, use `~/worktrees/<repo-short>/<branch>` per the global convention, never repo-siblings or nested worktrees.

## 1. What is already done (do NOT rebuild)

Contract §1–§15 and ADRs 0001, 0002, 0004, 0005, 0006, 0007, 0008, 0010, 0011, 0012, 0013, 0014 are implemented and tested: the native discrete-event engine (`ReactionNetworkProblem` as an AA `@aagent`), the priority-weighted progressive-fill allocator, runtime append-only mutation, RNG/`seed=` threading, the eval-free JSON `ExprNode` IR with `from_json_model`/`to_json_model` round-trip + `validate`, structured/agentic tokens + the host-function registry, the interface lifecycle with declarative `population[]` and `dump_state`/`restore`/`reinit`, token filtration (`@select`/`@advance`), Rules/conditional transitions, action callbacks (`SetTokens`/`Invoke`), AA coupling (`inputs[]`/`ExternalRef`/`_prestep!` latch, `getobservable`/params), and the whole Phase-0.6 analysis & visualization layer (`token_trajectory`/`representative_token`/`trajectory_envelope`, `ensemble`/`summarize`/`treatment_effect`, `export_run`/`export_ensemble`, `network_graph`/`draw_network`/`exec_map`, the Plots recipes). Four literate demos exist: `demo/core_engine_tour`, `demo/agentic_pipeline`, `demo/aa_integration`, `demo/bd_acquisition`.

## 2. Dependency graph & recommended sequence

```
WS-3 (pinned bugs)  ──┐  independent, do first — small, unblocks confidence
WS-5 (introspection ──┘  independent — everything it needs already ships
       tutorial)

WS-1 Phase 1 (store swap) ──► WS-1 Phase 2 (ReactantSpec promotion) ──► WS-2 (refinement/0009)
                                                    └──► WS-1 Phase 3 (optional acset adapter)
WS-4 (housekeeping) — independent, low-risk, fold in opportunistically
```

**Critical path is WS-1 Phase 2 → WS-2.** ADR 0009 refinement is explicitly gated on the promoted `ReactantSpec` incidence table (the FK-repoint substrate it splices on). Everything else parallelizes. Recommended order for a single agent: WS-3 → WS-5 (quick wins, build context) → WS-1 Phase 1 → WS-1 Phase 2 → WS-2, folding WS-4 in as touched, WS-1 Phase 3 optional.

---

## WS-1 — Data-store migration (ADR 0003) — the largest unbuilt pillar

**Status:** NOT landed. The engine still runs on the Catlab/ACSets `BasicSchema` (`TheoryReactionNetwork` in `src/ReactiveDynamics.jl`), and the transition↔reactant relation is still parsed out of `trans`/`transLHS`/`transRHS` at runtime. `ReactantSpec` currently appears only in a comment in `src/visualize.jl`.

**Decision (ADR 0003, verbatim intent):** replace ACSets with a dependency-free typed-struct-of-columns IR, as a phased migration.

- **Phase 1 — behavior-preserving store swap.** A struct of typed columns per object; a single `const SCHEMA` as the source of truth replacing the `propertynames(acs.subparts)` reflection sites; a signature-preserving shim implementing `getindex`/`setindex!`/`nparts`/`parts`/`add_part!`/`add_parts!`/`incident`/`subpart`/`set_subpart!`/`rem_parts!`. Runtime untouched (`state.jl` already bypasses the store via its own indexing). Guarded by the existing Phase-0 characterization tests — nothing should change behaviorally. If the maintainer stops here, this is a defensible terminal state.
- **Phase 2 — promote the reactant relation** to a first-class typed `ReactantSpec` incidence table with integer FKs (`reactant_trans → :T`, `reactant_species → :S`, plus `stoich`/`side`/modality), WITH an `Expr` escape-hatch for the legitimately dynamic reactants (`@choose`/`@move`/expression-valued stoichiometry). This makes specs FK-checkable for agentic authoring, round-trippable in JSON, and makes `equalize!` species-merge structurally exact (FK-repoint) instead of `recursively_substitute_vars!` string surgery. (`@structured` is no longer in this escape-hatch set — its genesis form is now a typed, registry-resolved `structured{kind, fields}` reactant that round-trips through the JSON IR; the raw inline-constructor form was removed. See ADR 0005 §21.)
- **Phase 3 (optional) — a tested `to_acset`/`from_acset` adapter** behind a weakdep package extension, preserving an optional AlgebraicPetri/Catlab static-spec view at near-zero core cost. Off the critical path; ship only if an interop consumer materializes.

**Files (re-verify lines):** `src/ReactiveDynamics.jl` (schema, `merge_acs!`, `assign_defaults!`, `add_obs!`, `Base.convert` coercions), `src/state.jl` (indexing pass-throughs), `src/compilers.jl` (reactant/attr reads), `src/operators/joins.jl` + `src/operators/equalize.jl` (the string-surgery merge that Phase 2 makes structural), `src/interface/create.jl` + `reaction_parser.jl` (where LHS/RHS get folded — the Phase-2 population point for `ReactantSpec`), `src/visualize.jl` (the `network_graph` Layer-A incidence walk simplifies onto `ReactantSpec` when it lands). `Project.toml`/`Manifest.toml` (drop `ACSets`/`Catlab` from `[deps]` after Phase 1; add the weakdep in Phase 3).

**Acceptance:**
- Phase 1: the full existing suite stays green with ACSets removed from `[deps]` (behavior-preserving). Add a shim-parity test if one is not already implied by the characterization layer.
- Phase 2: species-merge exactness test in `test/semantic/determinism_composition_bugs.jl` (FK-repoint, no string corruption); JSON round-trip of the `ReactantSpec` table (`test/semantic/serialization_ir.jl`); T2-acceptance tests that name `ReactantSpec` FK-repoint flip green.

**Risks:** the shim must be signature-exact or the runtime silently diverges; the `Expr` escape-hatch is essential — do NOT try to make every reactant a static FK row (structured/`@choose`/expression-stoich reactants are legitimately dynamic). Land the store swap (Phase 1) SEPARATELY from the structural promotion (Phase 2) — that separation is the whole point of the phasing.

---

## WS-2 — Hierarchical refinement & multi-granularity (ADR 0009) — the headline unbuilt capability

**Status:** NOT implemented at all (no `refine`/`abstract`, no `@compose`/`@pipeline`/`@process`, no port `role`, no tests). Still "Proposed". **GATED on WS-1 Phase 2** — it reuses the `ReactantSpec` FK-repoint as its splice mechanism. Becomes CONTRACT §11 (already stubbed) + a §7 addendum.

This is the maintainer's "compact/expressive, compositional, various levels of granularity with refined dynamics substituted" ask — model a portfolio coarsely, then zoom one transition (e.g. the Phase-2 bottleneck) into a finer sub-model without disturbing the rest.

**Decision (ADR 0009), build in this order:**
- **(A) Open ports** — a thin closed-tag `role::PortRole ∈ {private, input, output, shared}` annotation on the Species record (CONTRACT §6.4, NOT a new table; default `private`). `private` auto-namespaces on compose (today's `m__X`), `input`/`output` are open ports (directionality advisory), `shared` is bare-name identified (first-class `@catchall`). A transition's boundary IS its LHS/RHS `ReactantSpec` rows — no new structure.
- **(B) `refine(spec, transition, submodel; ports)`** — splice a sub-model into a coarse transition in four authoring-time structural moves: namespace `S`'s private species; identify ports with boundary species by FK-repoint (the WS-1 Phase-2 operation — no `recursively_substitute_vars!`); append `S`'s transitions + remaining species/params/obs/events (this also fixes the §7/J4 `:E`/`:obs`-drop bug uniformly); remove the coarse transition `T`. Boundary species keep the same indices/names/cost/valuation, so the model is plug-compatible — nothing else notices. `abstract(...)` is the inverse (collapse a sub-graph to one coarse transition).
- **(C) Advisory boundary-consistency checks** in `validate` — coarse `cycletime ≈ Σ path cycletimes`, `pos ≈ Π sub-PoS`, `cost ≈ Σ sub-costs`, port-balance (every input consumed, every output produced). WARNINGS the author can override, not equivalence proofs.
- **(D) `@pipeline` sugar** (a chain of phases → N `flow`-genesis routing transitions using the phase-as-attribute `@select`/`@advance` idiom) and **`@process` modules** (named parameterized `ModelSpec` fragments, eval-free param substitution).
- **(E) `@compose f1 f2 …`** — `@join` plus automatic port matching by FK-repoint; closes the §7/J4 (`:E`/`:obs` merge) and J9 (`include_model`) bugs en route because it composes already-parsed `ModelSpec`s. `@join`/`@equalize` remain as the manual no-declared-ports path.
- **(F) Entity-level refinement** (a structured token hosting its own sub-network) is DEFERRED to a future ADR — design (A) so it is not precluded, but do not build it in v1.

**Files:** new authoring-time operators (likely `src/operators/refine.jl` + additions to `joins.jl`); `src/ReactiveDynamics.jl` (the `role` Species field / §6.4 annotation); `src/interface/` new macros for `@pipeline`/`@process`/`@compose`; `validate` (wherever it lives post-WS-1) for the §C advisory checks. Refinement is authoring-only and reindexes — it MUST be forbidden on a live/stepping model (the §7.5/J8 phase guard, ADR 0007 §A).

**Acceptance:** new `test/semantic/` file (e.g. `refinement_composition.jl`): plug-compatibility (transitions outside `{T}∪S` structurally unchanged after `refine`), round-trip of a refined spec through JSON, `@pipeline` expands to the expected `flow` transitions, `@compose` port-matching by FK-repoint, and the §C advisory diagnostics fire on a deliberately-drifted refinement. The `@test_broken` `:E`/`:obs`-merge and `include_model` pins in `determinism_composition_bugs.jl` flip green via `@compose` (coordinate with WS-3).

---

## WS-3 — Close the 7 pinned engine bugs (`@test_broken`)

**Status:** deliberately pinned so the regressions stay visible. Independent, small, do first. Locate each with `grep -rn "@test_broken" test/semantic/*.jl` and re-verify the `file:line` in source before fixing.

- **`@join` drops observables (`:obs`) and events (`:E`)** of joined submodels (`test/semantic/determinism_composition_bugs.jl` ~L216; source in `src/operators/joins.jl`, the dead `prepend_obs`). NOTE: WS-2's `@compose` fixes this uniformly — decide whether to fix `@join` directly here or fold into WS-2. If WS-2 is near, prefer the uniform fix and just flip the pin there.
- **`@join` file branch calls undefined `include_model`** (`determinism_composition_bugs.jl` ~L295; `joins.jl:226,228`). Same note — WS-2 `@compose` never takes this branch.
- **`ceil` breaks dt-invariance** on the non-Poisson genesis path (`modality_genesis.jl` ~L276; `solvers.jl` spawn count). Gate the `ceil` to the Poisson path per CONTRACT §2.3.
- **The remaining pins** in `reference_models.jl` and `allocation_conservation_lifecycle.jl` — read the in-file `# note:` lines; several document lifetime-terminated-instance pruning (`solvers.jl` `finish!` `filter!` keeps `state < cycleTime`, so max-lifetime-terminated instances re-emit RHS every tick — a conservation violation).

Also on the engine-bug list from the audit (verify still present, add tests if fixing): `free_blocked_species!` undefined `q` on the `:nonblock` release path; `set_bound_transition!` passed a `.bound_transition` instead of the token agent (`solvers.jl:288,452`); `delete!` on a Vector with an Int index (should be `deleteat!`, `solvers.jl:455`); `event_action!` fetches `:eventAction` but never evaluates it (though ADR 0010 Rules are the sanctioned decision channel now, so the legacy event no-op may be intentionally superseded — confirm before "fixing").

**Acceptance:** each fixed bug's `@test_broken` → `@test`; add a regression assertion where the pin only documented the target.

---

## WS-4 — Housekeeping (independent, low-risk)

- **Delete or implement exported-but-undefined macros:** `@agentize` (`src/interface/solve.jl:1`), `@export`/`@import` (`src/loadsave.jl:5`), `@prob_role`/`@list_by_role`/`@list_roles` (`src/interface/update.jl:6`). The role exports gesture at a roles/actors ontology that ADR 0009's `PortRole` partially subsumes — decide delete-vs-repurpose in coordination with WS-2. Add a test asserting every `export`ed symbol resolves.
- **Refresh `INVENTORY.md`:** its module table predates `src/exprnode.jl`, `src/predicates.jl`, `src/actions.jl`, `src/ledger.jl`, `src/serialize.jl` and undercounts grown files (see its own 2026-06-30 drift note). Re-run the inventory once WS-1 settles the store, since the schema section will change most.
- **`@prob_check_verbose`** is defined but calls an undefined `check_params` — implement or remove.

---

## WS-5 — Introspection & exec-map tutorial (independent — ships today)

**Status:** the capability is fully built and green (`test/semantic/analysis_observability.jl`, `test/semantic/visualization.jl`) but only surfaces INSIDE `demo/bd_acquisition`'s bespoke figures. The maintainer's headline ask — "a system exec map readable for (in)efficiencies, decorated with simulation results" — has no standalone literate demo. Everything it needs already exists; this is a pure-authoring task with no engine risk.

**Deliverable:** a new literate demo `demo/introspection_tour/` (`.jl` + `README.md`, matching the existing demos' style — block-comment narration, `println` results, an explicit `seed=`, every construct lifted from the passing semantic tests). Cover, on one small self-contained model: reading `prob.sol`/`prob.log`/`program_ledger`; the per-token trajectory log with `representative_token` (medoid) + `trajectory_envelope` (median/IQR); `ensemble`/`summarize`/`treatment_effect` (deterministic `hash((root,k))` seeding); `export_run`/`export_ensemble` bundles; and the three-layer network map — `network_graph` → `draw_network` (Graphviz via AA's `run_graphviz`) → **`exec_map`** with bottleneck/starvation coloring and `@select` token-path highlighting. Note the demo env needs its own `Project.toml` declaring `Plots`/`Arrow` (the Phase-0.6 weakdep demotion means `--project=.` no longer resolves them), mirroring `demo/bd_acquisition/Project.toml`.

**Acceptance:** `julia --project=demo/introspection_tour demo/introspection_tour/introspection_tour.jl` runs end to end and emits a rendered exec map. It invents no API (cross-check every call against `src/analysis.jl`/`src/visualize.jl`/`src/export.jl`).

---

## 3. Deferred / out of scope (gates noted, do NOT start without the gate)

- **Ensemble mode (b)** (reinit-reseed member reuse) — gated on ADR 0007 §D completing `_reinit!` (RNG/counter/population restore). Mode (a) rebuild-per-seed ships now. The mode-(b) acceptance test stays `@test_skip` with its gate note until §D lands.
- **Entity-level refinement** (ADR 0009 §F) — a future ADR; design WS-2(A) to not preclude it.
- **AA `@agentize` sugar + `Opera`-level implicit/fixed-point coupling** — current AA coupling is explicit Jacobi (one-tick lag) only; implicit algebraic-loop coupling is a separate AA-level design.

## 4. Definition of done for this handoff

WS-3 pins all flipped or consciously superseded; WS-5 introspection tour runnable and documented; WS-1 Phase 1 landed (ACSets out of `[deps]`, suite green) and Phase 2 landed (`ReactantSpec` FK table, exact species-merge, JSON round-trip); WS-2 refinement/composition (A)–(E) landed with its semantic tests green and CONTRACT §11 filled in; WS-4 export list authoritative and `INVENTORY.md` refreshed. Phase 3 acset adapter and all §3 items remain explicitly deferred with their gates recorded.
