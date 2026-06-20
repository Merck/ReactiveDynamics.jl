# ReactiveDynamics.jl — Phase-0 Review & Sign-off Summary

> **Status: SIGNED OFF 2026-06-20.** The maintainer approved the contract (§1–§9) and ADRs 0001–0006 and resolved all remaining scoping calls (§5). Phase 1 implementation may begin. Phase 0 made no `src/` changes.

This is the single-page map of everything Phase-0 produced. Phase 0 is "contract before code": no production `src/` changes landed; everything here is documentation and tests; the engine is unchanged. All artifacts were adversarially verified against the current `ref-agents` source (several verifier passes ran the engine on Julia 1.12.5), and every `file:line` citation was checked — re-verify before acting, as line numbers drift.

## 1. What Phase 0 decided (the short version)

ReactiveDynamics is a **timed, stochastic, resource-constrained Petri net / discrete-event system** for modeling business / R&D processes (budgeting, ledgers, what-if, rNPV) — not a chemical reaction network. The `ref-agents` branch had already pivoted to a native engine; Phase 0 ratified and specified that pivot rather than re-litigating it. Five decisions, all recorded as ADRs:

| ADR | Decision | Status |
|---|---|---|
| [0001](adr/0001-discrete-event-engine.md) | Keep the native discrete-event engine `ReactionNetworkProblem` (an AlgebraicAgents `@aagent` stepped via `_step!`); SciML demoted to an optional extension. | Accepted 2026-06-17 |
| [0002](adr/0002-priority-weighted-allocation.md) | Priority-weighted allocation via **weighted progressive filling** (water-filling) — work-conserving, conjunctive-consistent, deterministic, dependency-free; replaces the 7-function per-resource-split tangle. | Accepted 2026-06-17 |
| [0003](adr/0003-data-store.md) | **Drop ACSets** for a dependency-free typed-struct-of-columns IR; **promote** the transition↔reactant relation to a first-class typed `ReactantSpec` incidence table; no backward on-disk compatibility. | Accepted 2026-06-18 |
| [0004](adr/0004-runtime-mutation.md) | **Runtime mutation = append-only + soft-deactivate**, so transitions/species/params can be added (and transitions retired) *during* simulation without breaking the position-indexed compiled closures. | Accepted 2026-06-18 |
| [0005](adr/0005-serialization-json-ir.md) | **Single JSON serialization** + typed `ExprNode` IR — eval-free, JSON-Schema-describable for LLM emission/validation; closes the import-time RCE; drops the TOML/CSV/JLD2 zoo. | Accepted 2026-06-18 |
| [0006](adr/0006-structured-tokens.md) | **Structured/agentic tokens**: live instantiation/query (BD projects as entities), append-only-safe via `entangle!`; **custom-function registry replaces `@register`** (host-Julia vs serializable-data boundary, eval-free). | Accepted 2026-06-20 |

## 2. The modeling contract (`docs/CONTRACT_DRAFT.md`)

Complete, §1–§9. This is the normative specification the engine must satisfy.

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
| 9 | Structured Tokens & Queries | Cross-references ADR 0006: token instance lifecycle (instantiate/bind/move/unbind/retire), append-only-safety via `entangle!`, 7 invariants incl. D4 token total-order; the deterministic query API; the `TokenAgg` ExprNode; the host-Julia-vs-data boundary and the custom-function registry replacing `@register`. |

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

### 5a. Sign-off (the gate) — GRANTED 2026-06-20
The maintainer approved the contract (§1–§9) and ADRs 0001–0006 as the Phase-0 baseline, and accepted all four remaining scoping calls (5b). **Phase 0 is signed off; Phase 1 implementation may begin.** No `src/` changes were made during Phase 0.

### 5b. Scoping calls — RESOLVED 2026-06-20
- **ADR 0003 — Phase 3 interop adapter → DROPPED/deferred.** No `to_acset` weakdep view shipped; interop is off-roadmap and impossible for the running engine. Reversible later.
- **Structured-token determinism source (ADR 0006) → per-species creation counter.** The k-th token of a species gets a monotonic creation index; sort-key = `(species, creation_index)`. The uuid stays as identity; no RNG-threaded uuid generation. Counter is reset by `_reinit!`.
- **Retired-token growth (ADR 0006) → accept for Milestone 1 + escape hatch.** Soft-`:removed` is the default (preserves `past_bonds`); `disentangle!` is the opt-out; periodic archival only if a real run shows the O(#tokens) cost bites.
- **Genesis default (CONTRACT §2.8) → engine default `poisson`, template/authoring default `flow`.** Keeps existing models stable while steering BD/pipeline routing transitions to the correct semantics.

### 5c. Earlier-resolved (for the record)
- **ADR 0002:** priority is dynamic/time-varying; per-tick fairness; `priority=0`=leftover-only.
- **ADR 0003:** no on-disk back-compat; promote the reactant; single JSON.
- **`@register` / custom functions (ADR 0006).** Removed as a model-authoring/eval path; replaced by a per-network **registry** of host-supplied Julia functions referenced BY NAME (closed allow-list, eval-free) — this is why `@register` had to `@eval` into RD's module (the compiled closure's bare call-head, `compilers.jl:28,122`) and why the registry dissolves the need.
- **Action statements (ADR 0005/0006).** `{SetSpecies, SetParams, Log, Seq}` confirmed sufficient for now.
- **Resource retirement (ADR 0004 + 0006).** Species KINDS stay defined a priori (no live species-kind retirement needed); transitions retire live via `deactivate!`; structured-token INSTANCES are created and retired live (soft `:removed` / hard `disentangle!`). No `specActivated` flag is added. (A residual: a *plain* species can still only be orphaned, not retired, mid-run — accepted, since the live-retirement need is on structured-token instances, which ADR 0006 fully supports.)

## 6. Provenance

Phase-0 artifacts and commits on `ref-agents`: `b69cfb9` (inventory + ADR 0001/0002), `ebccf9f` (ADR 0003 + contract §1–§5), `3bf2e88` (accept 0003, add 0004/0005, genesis), `c3f6161` (contract §6–§8 + test suite), `4d2c109` (this review summary), `7867669` (ADR 0006 + contract §9), plus the sign-off commit recording the §5 resolutions. The ground-truth source map is `INVENTORY.md` (supersedes the stale architecture map in `REVIEW.md`). Workflows that produced these were adversarially verified; corrections from engine-run verifiers were applied (notably: `simulate(a, max_t)` takes a max-*time* and returns the agent; the `seed=` kwarg is currently swallowed; events are non-functional).

## 7. What Phase 1 looks like (after sign-off)

Implement the contract against the native engine: the typed IR + `const SCHEMA`; the weighted-progressive-filling allocator (`build_requirements!`/`progressive_fill!`/`spawn_integer!`); orthogonal typed modalities + construction-time validation; the promoted `ReactantSpec` table; an `AbstractRNG` + `seed=` threaded through every draw; the append-only live mutation API; the structured-token subsystem (the per-network function/kind registry replacing `@register`, the deterministic token query API + `token_sortkey` total order, the `Construct{kind,args}`/`TokenAgg` ExprNodes, and the token-binding bug fixes at `solvers.jl:288,452,455,476,512`); and the bug fixes in §4. Acceptance = the T2 tests flip green, SIR + toy-pharma reproduce known behavior under seed, and the BD projects-as-tokens scenario (ADR 0006 north-star) runs with the acquisition lever applied live.
