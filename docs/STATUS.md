# ReactiveDynamics.jl — Status & Remaining Work

> The single entry point for "what is the state, and what is left." Last updated 2026-07-15 on branch `ref-agents` (HEAD at the ADR-status truth-up + docs reorg). Re-verify any `file:line` before acting — line numbers drift; grep to confirm.

## Where the durable spec lives

- **[CONTRACT_DRAFT.md](CONTRACT_DRAFT.md)** — the normative operational-semantics spec (§1–§15). This is the durable spine: it pins the time model, per-tick firing/lifecycle rules, modality semantics, determinism/seeding obligations, typed attribute domains, the object model, composition semantics, serialization schema, structured tokens & queries, and the Phase-0.5/0.6 extensions (interface & initial state §10, refinement & composition §11, rules/decisions §12, AA integration §13, analysis & observability §14, visualization §15).
- **[adr/](adr/)** — one Architecture Decision Record per real decision (0001–0015), append-only, with a status table in [adr/README.md](adr/README.md). ADRs record context/decision/consequences/rejected-options; the CONTRACT records the resulting semantics.
- **[../INVENTORY.md](../INVENTORY.md)** — the current-source map (module map, public-API audit, static store, stepping trace, AA touchpoints). Verified against the `ref-agents` tree.
- **[../readme.md](../readme.md)** — the user-facing package README (about, four sketches, demos).
- **[../CLAUDE.md](../CLAUDE.md)** — the agent-facing guide (what the project is, current post-ADR-0015 architecture, build/test commands, dev-loop gotchas). `REVIEW.md` at the repo root is a HISTORICAL 2026-06-15 review snapshot, not a current guide.

Design/record documents (not normative, kept for provenance): **[PHASE0_REVIEW.md](PHASE0_REVIEW.md)** (the Phase-0 sign-off + Phase-0.6 work-package record) and **[MVP_BD_DEMO.md](MVP_BD_DEMO.md)** (the Business-Development acquisition-impact demo design doc that drove §12 and catalogued findings A–I). Superseded planning docs are archived under **[history/](history/)**.

## Overall state

The modeling + analysis + visualization surface is BUILT and green. The full contract (§1–§15) and ADRs 0001–0015 are implemented and tested; the suite runs FULLY green under `test/semantic/` — **801 pass / 0 broken / 801 total**, no `@test_skip` placeholders and no `@test_broken` pins remaining (the last 3 — the CONTRACT §1.4 construction-time modality validators — landed in `53d1fac`). AlgebraicAgents is the published registry release 0.4 (the earlier `Merck/AlgebraicAgents.jl@main` `[sources]` pin was dropped, `eb2ee10`).

### Implemented (with implementing commit)

| Area | ADR / §CONTRACT | Where it lives | Commit |
|---|---|---|---|
| Native discrete-event engine (`ReactionNetworkProblem` as an AA `@aagent`) | 0001 / §1–§3 | `src/solvers.jl`, `src/state.jl` | Accepted 2026-06-17 |
| Priority-weighted progressive-fill allocator | 0002 / §3 | `src/solvers.jl` | Accepted 2026-06-17 |
| Data store: drop ACSets → typed struct-of-columns + `ReactantSpec` promotion | 0003 / §6 | `src/ReactiveDynamics.jl` (`SCHEMA`, `ReactantSpec` `:130`) | Phase 1 `7932111`, Phase 2 `54e3c02` |
| Runtime mutation (append-only + soft-deactivate) | 0004 / §6.10 | `src/solvers.jl`, `src/state.jl` | Accepted 2026-06-18 |
| Single-JSON eval-free serialization + typed `ExprNode` IR | 0005 / §8 | `src/serialize.jl`, `src/exprnode.jl` | (accepted 2026-06-18) |
| Structured/agentic tokens + host-function registry (post-`@register`) | 0006 / §9 | `src/interface/agents.jl` | Accepted 2026-06-20 |
| Interface lifecycle + declarative `population[]` + `dump_state`/`restore` | 0007 / §10 | `src/interface/checkpoint.jl`, `src/solvers.jl` | `44c9d55` |
| Token filtration (`TokenPredicate`/`@select`, `SetField`/`@advance`) | 0008 / §9.5 | `src/predicates.jl`, `src/actions.jl` | `896bc45` |
| Hierarchical refinement + open-port composition (`refine`/`abstract`/`@compose`/`@pipeline`/`@process`) | 0009 / §11 | `src/operators/refine.jl` | `886a610` |
| Rules/triggers & conditional transitions (endogenous decision channel) | 0010 / §12 | `src/actions.jl` | Accepted 2026-06-21 |
| Action callbacks: `SetTokens` population writes + `Invoke` escape hatch | 0011 / §12.3 | `src/actions.jl` | Accepted 2026-06-21 |
| AA integration: RD as a hierarchy node + wired external coupling (`ExternalRef`, `_prestep!`) | 0012 / §13 | `src/interface/aa_coupling.jl`, `src/exprnode.jl` | `d5a367c` |
| `@agentize` thin constructor sugar | 0012 / 0001 | `src/interface/solve.jl` | `ebf7773`/`4f61cfa` |
| Analysis & observability: token trajectory log, ensemble runner, export bundle | 0013 / §14 | `src/analysis.jl`, `src/export.jl`, `ext/RDArrowExt.jl` | `6174ebe` |
| Ensemble mode (b) reinit-reseed (`ensemble(...; mode = :reinit)`) | 0013 / §14.2 | `src/analysis.jl`, `src/solvers.jl` (`_reinit!` `:964`) | `1501c78` |
| Visualization: result-plot recipes + three-layer network exec map | 0014 / §15 | `src/visualize.jl`, `ext/RDPlotsExt.jl` | `6174ebe` |
| Post-ACSets naming rename + drop GeneratedExpressions | 0015 | `src/` (rename), deprecation shims `src/ReactiveDynamics.jl:427-435` | `bc6cc0a` |

## Recently closed

- **Construction-time modality validation (CONTRACT §1.4) — DONE (`53d1fac`).** `validate_modalities` in `src/solvers.jl` runs in the `ReactionNetworkProblem` constructor (before any closure compiles or tick runs) and rejects the three illegal §1.4 configs with a clear `ArgumentError`: `{:nonblock, :conserved}`; `:rate` (perstep) with a concrete `cycletime == 0`; `:rate` on a structured species. It reuses the eval-free static reactant decomposition and mirrors the runtime modality union; Expr/param-valued cycletimes and `@choose`/bidirectional lines are left to run (not statically decidable / the documented escape hatch). The three former `@test_skip` acceptance placeholders in `test/semantic/modality_genesis.jl` are now live passing `@test_throws`, and two collateral test models that had built the now-forbidden `@rate`+ct=0 foot-gun were legalized. Suite is fully green with zero skips.

## Remaining work (all deferred; each with its gate — nothing non-deferred is open)

1. **`dump_state` in-flight limitation (implemented WITH a documented constraint, not a gap).** ADR 0007 §C `dump_state`/`restore` shipped, but `dump_state` REFUSES to dump when any transition is mid-cycle — it errors unless `ongoing_transitions` is empty (`src/interface/checkpoint.jl:37-41`), i.e. dump is a clean-tick-boundary operation only (Milestone-1 scope). Serializing in-flight `ongoing` instances eval-free is the open ADR 0007 §C question; not a bug, an accurately-scoped limitation.
3. **Entity-level refinement (ADR 0009 §F) — deferred.** A structured token hosting its own sub-network. `refine`/`abstract`/`@compose`/`@pipeline`/`@process` (transition-level refinement + composition) all shipped; only the entity-hosts-a-subnetwork axis is deferred to a future ADR. §A of ADR 0009 was designed so as not to preclude it.
4. **Threaded ensemble backend — deferred.** `ensemble(...; parallel = true)` is accepted but currently runs sequentially (`src/analysis.jl:294`); a real threaded backend is future work, orthogonal to mode (b) (members are independent, so it is a safe extension).
5. **AA `Opera`-level implicit/fixed-point (algebraic-loop) coupling — deferred.** Current AA coupling is explicit one-tick-lag Jacobi (`src/interface/aa_coupling.jl:123`); a within-tick fixed point is a separate AA-level `Opera` design (noted in ADR 0012).
6. **ADR 0003 Phase 3 (`to_acset` weakdep interop view) — REJECTED / will-not-do (2026-07-15).** No AlgebraicJulia/ACSets/Catlab interop view will be built; interop is off the BD/rNPV roadmap and impossible for the running stateful engine. This is a definitive rejection, not a parked/deferred option — any future interop would be a new ADR against the then-current engine (the dependency-free typed IR is the permanent core). Recorded in ADR 0003 Resolved and PHASE0_REVIEW §5b. NOT counted as remaining work.
7. **Follow-up naming pass (`dt`/`tstep`) — DONE 2026-07-16.** CONTRACT_DRAFT flagged `tstep` as an internal alias / rename candidate; explicitly out of scope for ADR 0015, folded into this separate pass. Landed: the internal keyword-bag key and the accepted keyword now both use `dt` (matching the `state.dt` field, reconciled by name not position); the legacy `tstep` meta keyword is an accepted-but-deprecated alias mapped to `dt` with a depwarn for one release. Also fixed a latent silent-ignore bug (a user-supplied `tstep` was previously overwritten before use).
8. **Tutorial refinement — separate later pass.** A dedicated pass will refine the demos/tutorials (`demo/core_engine_tour`, `demo/agentic_pipeline`, `demo/introspection_tour`, `demo/refinement_tour`, `demo/aa_integration`, `demo/wires_viz_tour`, `demo/bd_acquisition`). Not part of this docs-consolidation pass.

## Build / test / dev

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install deps
julia --project=. -e 'using Pkg; Pkg.test()'          # run the semantic suite
julia --project=. -e 'using ReactiveDynamics'         # load the package
```

Tests live in `test/semantic/*.jl` (entry `test/runtests.jl` → `test/semantic/runtests.jl`), with their own `test/Project.toml` that declares the test-only deps (`Plots`, `Arrow`, `DataFrames`, `Distributions`, …) that the main project keeps as weakdeps. Formatting is JuliaFormatter with `.JuliaFormatter.toml`. The optional `dev/` Revise-based test server (`dev/run.sh`, `dev/test_server.jl`) avoids recompilation across edits; NB adding a NEW `src/` file requires clearing the compiled cache (Revise won't pick it up otherwise). Weakdep-backed features (Plots recipes, Arrow export) load only when their package is present — several `@test_skip`s in `analysis_observability.jl`/`visualization.jl` are environment guards for that, not gaps.

## History

- **[history/HANDOFF_PLAN.md](history/HANDOFF_PLAN.md)** — the WS-1..WS-5 remaining-work plan (2026-07-09). ALL DONE; archived.
- **[history/HANDOFF_PLAN_2.md](history/HANDOFF_PLAN_2.md)** — Handoff II: ensemble mode (b) + `@agentize` (2026-07-14). ALL DONE; archived.
