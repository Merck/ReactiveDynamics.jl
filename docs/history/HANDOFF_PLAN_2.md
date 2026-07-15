# ReactiveDynamics.jl — Handoff II: Ensemble mode (b) + `@agentize` sugar

> **⚠️ SUPERSEDED / HISTORICAL (archived 2026-07-15).** Both work items in this plan have SHIPPED: ensemble mode (b) reinit-reseed (`ensemble(...; mode = :reinit)`, commit `1501c78`) and the `@agentize` thin constructor sugar (commits `ebf7773`/`4f61cfa`). This file is retained only as a record of the design. For the CURRENT state and what genuinely remains, see [../STATUS.md](../STATUS.md).

> **Purpose (as written 2026-07-14).** A self-contained, agentic-handoff plan for the two items that [HANDOFF_PLAN.md](HANDOFF_PLAN.md) left explicitly deferred and that the maintainer now wants built: ensemble **mode (b)** (reinit-reseed member reuse) and the **`@agentize`** authoring sugar. Assumes NO prior conversation context. Read top to bottom, then re-verify current state (§0) before touching code — every `file:line` and test tally below drifts and MUST be re-checked. Source of truth for decisions: `../adr/`; normative spec: `../CONTRACT_DRAFT.md`; current-source map: `../../INVENTORY.md`. The already-landed work (WS-1..WS-5 of [HANDOFF_PLAN.md](HANDOFF_PLAN.md), and the whole Phase-0.6 layer) is recorded in `../PHASE0_REVIEW.md`.

## 0. Orient & verify before acting

The modeling + analysis surface is BUILT and green — the suite is currently **751 pass / 3 broken / 754** (the 3 "broken" are `@test_skip` placeholders for the construction-time modality-validation targets in `modality_genesis.jl`, unrelated to this handoff — Julia's summary counts skips as "Broken"). Re-establish the baseline first:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'                 # expect ~751 pass / 3 broken
```

AlgebraicAgents is on the registry release (0.4) — the old `Merck/AlgebraicAgents.jl@main` GitHub `[sources]` pin was dropped (commit `eb2ee10`). Tests live in `test/semantic/*.jl` under the two-tier scheme (T1-characterization / T2-acceptance). The acceptance ritual for both work items: add the target assertions as real `@test`s as the piece lands, and keep the suite green.

Dev-loop gotchas (re-verify): a `dev/` Revise-based test server avoids recompilation; ADDING A NEW `src/` file requires clearing the compiled cache; git worktrees go under `~/worktrees/<repo-short>/<branch>` per the global convention, never repo-siblings or nested.

## 1. The critical correction to the deferral note

`HANDOFF_PLAN.md §3` and the ADRs describe mode (b) as "gated on ADR 0007 §D completing `_reinit!` (RNG/counter/population restore)." **That gate is now CLOSED.** `AlgebraicAgents._reinit!(::ReactionNetworkProblem)` at `src/solvers.jl:964` already: restores the RNG from `initial_rng` (§4 D7); tears down the live token population (`disentangle!` every inner) and rebuilds the declarative `population[]` marking via `instantiate_population!`; resets `creation_counters`/`creation_index`; clears the per-program ledger and the per-token trajectory log; resets the `once`-rule latches and the external-input buffer. The reproducibility contract `init → step* → reinit! → step*` is a passing live `@test` (`test/semantic/analysis_observability.jl:99`, "reinit! clears the trajectory log and a re-run reproduces it").

So mode (b) is no longer blocked design — it is a bounded extension of an already-working `_reinit!`. What is genuinely missing is (a) a **reseed** path (`_reinit!` restores the SAME seed; mode (b) needs to reinit to t=0 AND install a *different* member seed) and (b) the ensemble wiring + member-retention shape. This handoff builds both.

## 2. Dependency graph

```
WS-B1 (ensemble mode b)  — independent — src/solvers.jl (reseed) + src/analysis.jl (ensemble wiring)
WS-B2 (@agentize sugar)  — independent — src/interface/solve.jl + module export
```

The two are disjoint (mode (b) touches `solvers.jl`/`analysis.jl`; `@agentize` touches `solve.jl`/the module). They can be built fully in parallel; the only shared files are `test/semantic/` (append-only, separate testsets) and `src/ReactiveDynamics.jl` (each may add one export line). Build in isolated worktrees, then integrate.

---

## WS-B1 — Ensemble mode (b): reinit-reseed member reuse (ADR 0013 §14.2, CONTRACT §14.2)

**Status:** mode (a) `:rebuild` ships and is green (`src/analysis.jl:280` `ensemble`, `:256` `EnsembleProblem`). Mode (b) `:reinit` is NOT implemented — `ensemble` hardcodes `mode = :rebuild` and the docstring records mode (b) as gated. The gate is closed (§1). No mode-(b) test exists yet.

**Why it matters.** Mode (a) rebuilds and recompiles a fresh `ReactionNetworkProblem` for every seed — for an `nseed`-large Monte-Carlo that is `nseed ×` construction + closure-compilation cost. Mode (b) builds ONE member, then reinit-reseeds and re-simulates it for each subsequent seed, reusing the compiled closures and the allocated store. It is the cheaper path the maintainer asked for.

**Decision — build in this order.**

- **(1) A reseed entry point.** Extend `_reinit!` (or add a thin `reinit!(prob; seed)` / `reseed!(prob, seed)`) so that, given a new seed, it reinstalls the stream before rebuilding the t=0 marking: set `state.rng = Random.Xoshiro(seed)`, `state.initial_rng = copy(state.rng)`, and store the realized seed on the struct (mirror the constructor at `src/solvers.jl:845-848`). CRITICAL ORDERING: the population is (re)sampled through `state.rng` inside `instantiate_population!` (`_reinit!` currently copies `initial_rng → rng` at `solvers.jl:972` THEN instantiates at `:1003`). The reseed must make the new seed's stream the one `instantiate_population!` draws from, so a `count`+attribute-distribution `population[]` entry yields the NEW seed's sampled initial attributes — identical to what a fresh `build(seed)` would sample. The default (no `seed` kwarg) path must remain byte-identical to today's `_reinit!`.

- **(2) The mode-(a) ≡ mode-(b) equivalence contract (the acceptance gate).** Mode (b) is only a legal substitute for mode (a) if the two produce the SAME ensemble member-for-member (same per-member seed `hash((root_seed, k))` → same trajectory). This holds ONLY when members are structurally homogeneous — same net, same `population[]` schema, differing only in the stochastic stream and seed-sampled initial attributes. State this precondition explicitly and guard it: mode (b) reseeds, it does NOT rebuild, so a `build` that branches structurally on its seed argument is out of contract for mode (b) and must use mode (a). Document the guard; do not attempt to detect arbitrary structural divergence.

- **(3) Member-retention shape (the one real design decision).** Mode (b) reuses ONE problem object, so it CANNOT hold `nseed` live members the way `EnsembleProblem.members::Vector{ReactionNetworkProblem}` does today (`analysis.jl:257`). RECOMMENDED shape: after each member's run, capture a lightweight per-member RESULT SNAPSHOT (copies of the artifacts the read surface needs — `sol`, `log`, `program_ledger` DataFrame, and `token_trajectory` if opted-in) rather than retaining the live agent; the reused problem is reinit-reseeded onward. `summarize`/`treatment_effect` (`analysis.jl:314,336`) take a `metric(member) -> Real`, so the snapshot must expose whatever those metrics read — pin the snapshot type to the artifacts the existing metrics touch (`program_ledger(m)`, `m.sol`) and keep `getobservable`/`observables` (`analysis.jl:366,383`) working off it. If a fully-faithful `metric(::ReactionNetworkProblem)` is required unchanged, the fallback is to apply the caller's metric EAGERLY per member inside the mode-(b) loop and store the scalar — but that changes the `EnsembleProblem` contract, so prefer the snapshot. Whichever is chosen, `length(inners(ens)) == nseed` (§14.2 Invariant 4, the drawable-node contract) must still hold — snapshots entangle as lightweight child nodes.

- **(4) Ensemble wiring.** Add `mode::Symbol = :rebuild` (or `:reinit`) as an `ensemble` kwarg. Mode `:reinit`: build member 1 (`seed = hash((root_seed, 1))`), simulate, snapshot; for `k = 2:nseed` reinit-reseed the SAME problem to `hash((root_seed, k))`, simulate, snapshot. Record the realized mode on the `EnsembleProblem` (the field already exists, `analysis.jl:260`) and `log` which mode ran (the docstring at `:271-274` already promises this — update it to describe both modes).

**Files (re-verify lines):** `src/solvers.jl` (`_reinit!` at `:964`; the reseed; the constructor seeding idiom at `:845-848` to mirror); `src/analysis.jl` (`ensemble` `:280`, `EnsembleProblem` `:256`, `summarize`/`treatment_effect` `:314/:336`, the read surface `:366/:383`); `test/semantic/analysis_observability.jl` (new mode-(b) testset alongside the existing reinit test at `:99`); the `ensemble` docstring (`analysis.jl:263-278`) and CONTRACT §14.2 (drop the "mode (b) gated" language once it lands).

**Acceptance:**
- Mode-(a) ≡ mode-(b) equivalence: `ensemble(build; nseed=N, root_seed=R, mode=:rebuild)` and `... mode=:reinit)` produce member-for-member identical `summarize` (and, on the BD scenario, identical `treatment_effect` Δ) — bit-identical where the metric is exact, otherwise within tight `isapprox`. THIS IS THE HEADLINE TEST.
- Reseed determinism: reinit-reseeding to `hash((R, k))` reproduces exactly what a fresh `build(hash((R, k)))` produces (a single-member A/B).
- The drawable-node invariant survives: `length(inners(ens)) == nseed` for both modes.
- Suite stays green; the mode-(b) `@test_skip`/placeholder (add one if you introduce the testset first) flips to real `@test`s.

**Risks:** the reseed ordering vs `instantiate_population!` is the subtle part — get the "new stream installed BEFORE the marking is re-sampled" order right or the initial attributes silently diverge from mode (a). The explicit-host-token `population[]` form (as in `demo/bd_acquisition`) restores the SAME objects to captured t=0 values (`restore_token_snapshot!`), so those initial markings are seed-INDEPENDENT and only the simulation draws differ — verify both population forms in the equivalence test. Keep the no-`seed` `_reinit!` path byte-identical (it backs the existing green reinit test).

---

## WS-B2 — `@agentize` authoring sugar (ADR 0012 future-work, ADR 0001)

**Status:** NOT implemented. The `export @agentize` was a dangling export (macro never defined) that WS-4 DELETED (`src/interface/solve.jl:1-5` carries the note; `INVENTORY.md:90` records it). ADR 0001 explicitly offers "implement `@agentize` as thin sugar over the constructor OR delete the export" — WS-4 took the delete; the maintainer now wants the sugar built. ADR 0012 / `docs/adr/README.md:25` list it as noted future work.

**Why it matters (and its bound).** Agentization is ALREADY implicit — `ReactionNetworkProblem(acs[, u0, p]; name, kwargs...)` builds the `@aagent` object and entangles the `"structured"` container. So `@agentize` is genuinely thin sugar; its only real value-add over the function call is ERGONOMIC, so scope it tightly and do not let it grow a second construction path.

**Decision — the minimal shape.** Define `@agentize` as a macro that expands to a `ReactionNetworkProblem(...)` call, with the one ergonomic win being AUTO-NAMING: when the acs argument is a plain binding, default the agent `name` to that binding's symbol (so `prob = @agentize mynet` gives `name = "mynet"` instead of the default `"reaction_network"`). Support passing `u0`/`p` positionally and forwarding `kwargs` (`seed=`, `tspan=`, `name=` overriding the auto-name). Suggested surface — confirm against the constructor signature at `src/solvers.jl:817`:

```julia
prob = @agentize net                      # name = "net"
prob = @agentize net u0 p seed=1 tspan=10 # forwards positionals + kwargs
prob = @agentize net name="custom"        # explicit name wins over auto-name
```

It MUST NOT reimplement any construction logic — it lowers to the exact constructor call. If the auto-naming proves to conflict with hygiene or the constructor's own `name` default in a way that adds complexity, drop the auto-name and make `@agentize` a pure pass-through macro; note the decision inline.

**Files:** `src/interface/solve.jl` (define the macro where the note now sits, `:1-5`); `src/ReactiveDynamics.jl` (re-add `export @agentize`); a testset (extend `test/semantic/exports_resolve.jl` — which already asserts every exported symbol resolves — and add a behavioral test that `@agentize` builds a working, simulatable problem and honors auto-naming); `INVENTORY.md` (flip the `@agentize` "REMOVED (WS-4)" rows at `:28/:90/:100/:213/:248/:280` to "implemented — thin constructor sugar"); a short note in ADR 0012 (and/or ADR 0001's alternative) recording that the sugar was implemented.

**Acceptance:**
- `@agentize` is exported and resolves (`exports_resolve.jl` guard stays green — this alone was the WS-4 failure mode).
- A behavioral test: `@agentize`-built problem simulates and matches the equivalent `ReactionNetworkProblem(...)` call (same result under a fixed `seed=`); auto-naming works and is overridable.
- Suite stays green; INVENTORY + the ADR note reflect the implementation.

**Risk:** macro hygiene — the auto-name reads the argument's symbol via `QuoteNode`/`string`; ensure a non-symbol acs expression (e.g. `@agentize build_net()`) still works (falls back to the constructor's default name). Keep it thin: the failure mode to avoid is `@agentize` accreting a divergent second constructor.

---

## 3. Definition of done for this handoff

WS-B1: `ensemble(...; mode=:reinit)` lands with the reseed path in `_reinit!`, the member-retention snapshot, and the mode-(a)≡mode-(b) equivalence test green; the "mode (b) gated" language is removed from the `ensemble` docstring and CONTRACT §14.2. WS-B2: `@agentize` is defined as thin constructor sugar, re-exported, resolves, has a behavioral test, and INVENTORY + the ADR note are updated. The full suite is green (no new `@test_broken`; the 3 pre-existing modality-validation `@test_skip`s are untouched and out of scope here).

## 4. Still deferred after this (unchanged)

- **AA `Opera`-level implicit/fixed-point (algebraic-loop) coupling** — current AA coupling is explicit Jacobi (one-tick lag) only; a within-tick fixed point is a separate AA-level `Opera` design (ADR 0012). NOT part of "`@agentize` sugar" — do not build it here.
- **Threaded ensemble backend** — `ensemble(...; parallel=true)` is accepted but runs sequentially (`analysis.jl:276`); a real threaded backend is future work, orthogonal to mode (b).
- **Entity-level refinement** (ADR 0009 §F) — a structured token hosting its own sub-network; a future ADR.
