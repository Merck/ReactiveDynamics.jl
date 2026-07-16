# DRAFT PR — ReactiveDynamics.jl: native discrete-event engine rework

> **Status: DRAFT — not ready to merge.** This is a local draft description (nothing pushed/opened). It intentionally stays draft until the docs/tutorial refinement lands (see "Explicitly out of scope → follow-up PR" below).
>
> **Head:** `rework` (renamed from `ref-agents`) · **Base:** `main` · **111 commits**, 154 files, +23,904 / −2,940. **Suite: 801 pass / 0 broken / 0 skipped**, fully green on Julia ≥ 1.11.

## What this is

`rework` is the full pivot of ReactiveDynamics.jl from a SciML/Catlab embedding to a **native, dependency-light discrete-event engine** for timed, stochastic, resource-constrained business/R&D process modeling (budgeting, ledgers, what-if, rNPV) — NOT chemical reaction networks, despite the Catalyst-derived DSL surface. The whole rework is specified contract-first: a normative operational-semantics contract (`docs/CONTRACT_DRAFT.md` §1–§15) and 15 Architecture Decision Records (`docs/adr/`), with the engine built to satisfy them.

A `ReactionNetworkProblem` IS an AlgebraicAgents `@aagent`, so a network is a node in a larger heterogeneous AA hierarchy. The single entry point for state is `docs/STATUS.md`.

## Headline changes

**Engine & semantics**
- Native discrete-event engine (`ReactionNetworkProblem` stepped via AA's `_step!`); SciML demoted to optional (ADR 0001).
- Priority-weighted progressive-fill (water-filling) resource allocator — work-conserving, deterministic, dependency-free (ADR 0002).
- Append-only + soft-deactivate runtime mutation, so transitions/species/params can be added and transitions retired mid-simulation without breaking position-indexed compiled closures (ADR 0004).
- `AbstractRNG`/`seed=` threaded through every draw; a run is fully determined by `(model, seed)` (CONTRACT §4).

**Data store & serialization**
- ACSets/Catlab **dropped** for a dependency-free typed struct-of-columns IR; the transition↔reactant relation promoted to a first-class typed `ReactantSpec` incidence table (ADR 0003). ACSets interop (a `to_acset` view) is **rejected / will-not-do**, not deferred.
- Single eval-free JSON serialization + typed `ExprNode` IR with `from_json_model`/`to_json_model` round-trip + `validate`; closes the import-time RCE; drops the TOML/CSV/JLD2 zoo (ADR 0005).
- Post-ACSets naming rename (ADR 0015): `@reaction_network` (was `@ReactionNetworkSchema`), `net` (was `acs`), store type `ReactionNetwork`, store verbs renamed to store vocabulary (`nrows`/`row_ids`/`column`/`cell`/`find_rows`/…) **and unexported**; old names survive one release as `@deprecate` shims; `GeneratedExpressions` dropped.

**Modeling language (Phase-0.5)**
- Structured/agentic tokens with live instantiation/query; host-function registry replaces `@register` (eval-free) (ADR 0006).
- Interface lifecycle (authoring → construction → live) + declarative serializable `population[]` initial marking + `dump_state`/`restore` (ADR 0007).
- Token filtration: `TokenPredicate`/`@select` selects tokens by 𝓕ₜ-measurable predicate; phase-as-attribute canonical; `@advance`/`SetField` field-writes (ADR 0008).
- Hierarchical refinement & open-port composition: `refine`/`abstract`/`@compose`/`@pipeline`/`@process` — substitutable granularity via FK-splice (ADR 0009).
- Rules/triggers & conditional transitions — the endogenous decision channel (ADR 0010); action callbacks `SetTokens`/`Invoke` (ADR 0011).
- AlgebraicAgents integration: RD as an AA hierarchy node (outbound `getobservable`/params) + inbound wired external coupling (`inputs[]`/`ExternalRef`/`_prestep!` one-tick Jacobi lag) (ADR 0012). AA moved to the published registry 0.4 release.

**Analysis & visualization (Phase-0.6)**
- Per-token trajectory log + `representative_token`/`trajectory_envelope`; ensemble runner (`ensemble`/`summarize`/`treatment_effect`, `EnsembleProblem` as an AA node); results export bundle (CSV/JSON core, Arrow weakdep) (ADR 0013).
- Model-agnostic Plots recipes behind `RDPlotsExt`; the three-layer network "exec map" (`network_graph`/`draw_network`/`exec_map`) with bottleneck/starvation coloring and `@select` token highlighting (ADR 0014).
- `Plots`/`Arrow` demoted to weakdeps with `RDPlotsExt`/`RDArrowExt` package extensions; `Pluto`/`PlutoUI`/`IJulia`/`DifferentialEquations` dropped from deps.

## This session's increment (on top of the rework)

The most recent work, folded into this branch:
- **Ensemble mode (b) — reinit-reseed member reuse** (`1501c78`, ADR 0013 §14.2): `_reinit!(state; seed=…)` reseed path + `ensemble(...; mode=:reinit)`, reusing one compiled member across seeds. Mode (a)≡mode (b) member-for-member equivalence is tested. The old "gated on ADR 0007 §D" caveat is closed — `_reinit!` fully implements §D.
- **`@agentize` thin constructor sugar** (`ebf7773` + hygiene fix `4f61cfa`, ADR 0001/0012): auto-naming macro over the `ReactionNetworkProblem` constructor; no second construction path.
- **Construction-time modality validator** (`53d1fac`, CONTRACT §1.4): `validate_modalities` rejects the three illegal modality configs (`{:nonblock,:conserved}`; `:rate` with concrete `cycletime==0`; `:rate` on a structured species) at construction with a clear `ArgumentError`. Closed the last 3 `@test_skip` acceptance placeholders → **zero skips remain**.
- **Docs consolidation & records truth-up**: authored repo-root `CLAUDE.md`; retired `REVIEW.md` to a historical snapshot; created `docs/STATUS.md` as the single state index; archived the two completed HANDOFF_PLANs under `docs/history/`; truthed-up all ADR statuses to Implemented with commit citations; reclassified ADR 0003 Phase 3 (ACSets interop) to rejected; fixed a fabricated CI claim; anchored a test's demo `include` to `pkgdir` (was a `homedir()` path — a cross-worktree hazard).

## Tests

`test/semantic/*.jl`, two-tier (T1-characterization / T2-acceptance). **801 pass / 0 broken / 0 skipped.** No `@test_broken` pins, no `@test_skip` placeholders. `test/Project.toml` declares the test-only deps (Plots, Arrow, DataFrames, Distributions) so the weakdep extension paths are exercised, not skipped. No CI is configured in the repo — the suite and JuliaFormatter are run locally.

## Explicitly out of scope → follow-up PR

**Docs/tutorial refinement is deliberately NOT in this PR.** This is why the PR stays a draft: `main` should not receive the rework until the demos/tutorials are polished and coherent against the final (post-ADR-0015) surface. The plan:
- Do the tutorial/demo refinement on a branch off `rework` (e.g. `docs-tutorials`) as its own reviewable PR **targeting `rework`, not `main`**.
- Seven demos are in scope (STATUS.md item 8): `core_engine_tour`, `agentic_pipeline`, `introspection_tour`, `refinement_tour`, `aa_integration`, `wires_viz_tour`, `bd_acquisition`. The old `tutorial/` dir was already removed in favor of these literate demos.
- When that follow-up merges back into `rework`, mark THIS draft PR ready and merge `rework → main` once.

## Genuinely deferred (gates recorded, NOT part of this or the follow-up)

- Entity-level refinement — a structured token hosting its own sub-network (ADR 0009 §F, a future ADR).
- Threaded ensemble backend (`ensemble(...; parallel=true)` is accepted but runs sequentially).
- AA `Opera`-level implicit/fixed-point coupling (current coupling is explicit one-tick-lag Jacobi).
- `dump_state` in-flight limitation — implemented but refuses to dump mid-cycle transitions (a scoped Milestone-1 constraint, not a gap).
- `dt`/`tstep` internal naming pass (flagged in CONTRACT, out of scope for ADR 0015).

## Notes for the reviewer / merge mechanics

- `main` and `rework` share a 2023 merge-base (`d0e195b`); `main` has one trivial commit beyond it (`6a09222` "Update readme.md") that is not on `rework` — reconcile at merge time (rebase/merge `main` or cherry-pick that readme edit) so the readme change is not lost.
- Branch was renamed `ref-agents` → `rework` this session; a stale `origin/ref-agents` remote ref exists (nothing was pushed this session). Push `rework` and update/retire the old remote ref when opening for real.
- Three local `worktree-agent-*` branches are leftovers from merged subagent work and can be deleted.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
