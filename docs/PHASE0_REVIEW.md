# ReactiveDynamics.jl — Phase-0 Review & Sign-off Summary

This is the single-page map of everything Phase-0 produced, for maintainer sign-off. Phase 0 is "contract before code": no production `src/` changes land until the contract and ADRs below are approved. Everything here is documentation and tests; the engine is unchanged. All artifacts were adversarially verified against the current `ref-agents` source (several verifier passes ran the engine on Julia 1.12.5), and every `file:line` citation was checked — re-verify before acting, as line numbers drift.

## 1. What Phase 0 decided (the short version)

ReactiveDynamics is a **timed, stochastic, resource-constrained Petri net / discrete-event system** for modeling business / R&D processes (budgeting, ledgers, what-if, rNPV) — not a chemical reaction network. The `ref-agents` branch had already pivoted to a native engine; Phase 0 ratified and specified that pivot rather than re-litigating it. Five decisions, all recorded as ADRs:

| ADR | Decision | Status |
|---|---|---|
| [0001](adr/0001-discrete-event-engine.md) | Keep the native discrete-event engine `ReactionNetworkProblem` (an AlgebraicAgents `@aagent` stepped via `_step!`); SciML demoted to an optional extension. | Accepted 2026-06-17 |
| [0002](adr/0002-priority-weighted-allocation.md) | Priority-weighted allocation via **weighted progressive filling** (water-filling) — work-conserving, conjunctive-consistent, deterministic, dependency-free; replaces the 7-function per-resource-split tangle. | Accepted 2026-06-17 |
| [0003](adr/0003-data-store.md) | **Drop ACSets** for a dependency-free typed-struct-of-columns IR; **promote** the transition↔reactant relation to a first-class typed `ReactantSpec` incidence table; no backward on-disk compatibility. | Accepted 2026-06-18 |
| [0004](adr/0004-runtime-mutation.md) | **Runtime mutation = append-only + soft-deactivate**, so transitions/species/params can be added (and transitions retired) *during* simulation without breaking the position-indexed compiled closures. | Accepted 2026-06-18 |
| [0005](adr/0005-serialization-json-ir.md) | **Single JSON serialization** + typed `ExprNode` IR — eval-free, JSON-Schema-describable for LLM emission/validation; closes the import-time RCE; drops the TOML/CSV/JLD2 zoo. | Accepted 2026-06-18 |

## 2. The modeling contract (`docs/CONTRACT_DRAFT.md`)

Complete, §1–§8. This is the normative specification the engine must satisfy.

| § | Section | Pins |
|---|---|---|
| 1 | Modality Truth Table | Orthogonal axes `allocation∈{upfront,perstep}` × `return∈{consumed,conserved}` × `blocking∈{block,nonblock}`; 5 legal rows + illegal combos; maps to the legacy `Set{Symbol}` and to ADR 0002's `build_requirements!`. |
| 2 | Time Model (+ §2.8 Genesis) | Single clock; Poisson spawn intensity; cycle/lifetime; dt-invariance hazards. §2.8: genesis is two-stage (rate proposal + upfront-LHS gate); closed `{poisson,scheduled,flow,capacity}` intent tag; source vs routing. |
| 3 | Operational Semantics | Instance lifecycle; the 13-step `_step!` tick; 7 contract invariants (non-negativity, conservation, capacity, integrality, determinism, termination-completeness, event-firing). |
| 4 | Determinism & Seeding | D1–D9: thread an `AbstractRNG`, per-trajectory seeding, reinit restores the stream. |
| 5 | Attribute Contract | Per-attribute units / range / default / time-varying policy; replaces the `SampleableValues` catch-all. |
| 6 | Object Model | Typed columnar tables under the append-only index invariant; the promoted `ReactantSpec` incidence table; identity-by-name; structured-token refinement. |
| 7 | Composition Semantics | Join (name-merge S/P/M, disjoint T) + the obs/`:E` merge gap; equalize; FK-repoint vs string-surgery; the `rem_parts!` live-guard; the undefined `include_model` bug. |
| 8 | Serialization Schema | Cross-references ADR 0005: round-trip + `(model, seed)`-determinism guarantees; eval-free `validate`; RCE closure; outputs→Arrow; the append-only mutation-patch form. |

## 3. The Phase-0 semantic test suite (`test/semantic/`)

**63 testsets** across 4 files, wired into `test/runtests.jl`, replacing the assertion-free tutorial smoke tests as the real acceptance layer. All files parse-clean on Julia 1.12.5.

Two tiers:
- **T1-characterization (39 tests)** — run against the *current* engine; they lock in correct behavior or **pin a known bug with `@test_broken`** so the suite is green-when-expected.
- **T2-acceptance (24 tests)** — encode the *target* behavior from the contract/ADRs; `@test_skip`-wrapped with a reference block because they name not-yet-built APIs (`progressive_fill!`, the `seed=` kwarg, `ReactantSpec` FK-repoint). Phase 1 flips each `@test_skip`/`@test_broken` → `@test` as the piece lands.

| File | Testsets | Covers |
|---|---|---|
| `allocation_conservation_lifecycle.jl` | 18 (T1=8, T2=10) | ADR 0002 allocation oracles; conservation & non-negativity (§3.4); instance lifecycle & cycle time. |
| `modality_genesis.jl` | 17 (T1=14, T2=3) | The 5 modality rows + illegal combos (§1); genesis modes poisson/scheduled/flow/capacity (§2.8). |
| `determinism_composition_bugs.jl` | 14 (T1=7, T2=7) | Seeding D1–D9 (§4); join/equalize (§7); the remaining bug-pins. |
| `reference_models.jl` | 14 (T1=10, T2=4) | SIR, toy-pharma, rNPV — the brief's headline acceptance criteria. |

## 4. Known engine bugs the suite pins (fix in Phase 1)

These survive on `ref-agents`; the contract documents them and the tests pin them so the regression is visible. (Re-verify line numbers.)

| Bug | Site | Effect |
|---|---|---|
| Undefined `q` in `free_blocked_species!` | `solvers.jl:512` | `:nonblock` resource-release path crashes whenever a `:nonblock` token is in-flight. |
| Lifetime-terminated instances never pruned | `solvers.jl:501` | `filter!` keeps `state < cycleTime`, so timeout-terminated instances re-emit RHS + re-return conserved tokens every subsequent tick — violates conservation. |
| `event_action!` is a no-op | `solvers.jl:323` | Fetches `:eventAction` but never evaluates it; events do not fire (so the acquisition lever must use the scheduled rate idiom, not events). |
| `add_to_spawn!` deferral broken | `state.jl:251-256` | `findfirst` gets a scalar, and it increments `:transHash` not `:transToSpawn` — capacity-overflow deferral is lost. |
| `resample!` writes nonexistent `o.val` | `state.jl:137` | Range-less observable crashes (field is `sampled`). |
| `@mode` bare `specModality` | `update.jl:108` | `UndefVarError` on any `@mode` call. |
| `ceil` breaks dt-invariance | `solvers.jl:144` | Fractional spawn counts round up every tick off the Poisson path. |
| No RNG seeding | `create.jl:151`, `solvers.jl:413,321`, `state.jl:123` | All draws use the global RNG; determinism-under-seed not achievable (the §4 gap). |
| `@join` file branch | `joins.jl:226,228` | Calls undefined `include_model`. |

## 5. Decisions needed from the maintainer

### 5a. Sign-off (the gate)
Approve the contract (§1–§8) and ADRs 0001–0005 as the Phase-0 baseline. This unblocks Phase 1 implementation. Nothing in `src/` changes until this approval.

### 5b. Open questions that shape Phase 1
These are recorded in the ADRs; a decision now avoids rework later.

- **Resource retirement (ADR 0004).** "Append-only" means transitions can be **retired live** via `deactivate!` (soft, reversible; in-flight instances finish) — but **species have no equivalent**: there's a `transActivated` flag, no `specActivated`. So a resource can only be *orphaned* mid-run (deactivate every transition touching it), not truly retired. Decision: is transition-level deactivation enough, or do you want a symmetric `specActivated` (a small addition to ADR 0004 + §6)?
- **ADR 0003 — Phase 3 interop adapter.** Ship the optional `to_acset` weakdep view for AlgebraicPetri interop, or drop it? (Interop is off the BD/rNPV roadmap.)
- **ADR 0005 — the `@register` user-function path.** Remove outright, or keep behind a closed named-function registry? (It's the remaining eval surface after the RCE fix.)
- **ADR 0005 — action-statement coverage.** Is `{SetSpecies, SetParams, Log, Seq}` enough for real pre/post and event actions, or do tutorials need richer statements?
- **Genesis default (CONTRACT §2.8).** Confirm `flow` (not `poisson`) is the right default for pipeline/routing transitions, with `poisson` reserved for true exogenous sources.

### 5c. Already resolved (for the record)
ADR 0002: priority is dynamic/time-varying, per-tick fairness, `priority=0`=leftover-only. ADR 0003: no on-disk back-compat, promote the reactant, single JSON. These are baked into the contract.

## 6. Provenance

Phase-0 artifacts and commits on `ref-agents`: `b69cfb9` (inventory + ADR 0001/0002), `ebccf9f` (ADR 0003 + contract §1–§5), `3bf2e88` (accept 0003, add 0004/0005, genesis), `c3f6161` (contract §6–§8 + test suite). The ground-truth source map is `INVENTORY.md` (supersedes the stale architecture map in `REVIEW.md`). Workflows that produced these were adversarially verified; corrections from engine-run verifiers were applied (notably: `simulate(a, max_t)` takes a max-*time* and returns the agent; the `seed=` kwarg is currently swallowed; events are non-functional).

## 7. What Phase 1 looks like (after sign-off)

Implement the contract against the native engine: the typed IR + `const SCHEMA`; the weighted-progressive-filling allocator (`build_requirements!`/`progressive_fill!`/`spawn_integer!`); orthogonal typed modalities + construction-time validation; the promoted `ReactantSpec` table; an `AbstractRNG` + `seed=` threaded through every draw; the append-only live mutation API; and the bug fixes in §4. Acceptance = the T2 tests flip green and SIR + toy-pharma reproduce known behavior under seed.
