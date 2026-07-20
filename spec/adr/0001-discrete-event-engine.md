# ADR 0001 — Native discrete-event engine (keep, harden; do not return to SciML)

- Status: Accepted (maintainer-confirmed 2026-06-17)
- Deciders: maintainer + rework
- Relates to: rework brief Workstream D ("Engine / SciML embedding"); supersedes REVIEW.md's `DiscreteProblem`-era framing; see [INVENTORY.md](../INVENTORY.md) for current source map.

## Context

ReactiveDynamics is a timed, stochastic, resource-constrained Petri net / discrete-event system: transitions are stateful rules that spawn instances at a Poisson rate, occupy resources over a cycle time, and complete with a terminal success probability. Mass-action CRN is one degenerate special case, not the conceptual core. A memoryless SciML `DiscreteProblem`/`FunctionMap` update cannot cleanly carry in-flight instance state (genesis time, accumulated cycle progress, bound resources), which was a structural mismatch and a likely source of latent bugs in the v0.2.x line.

The `ref-agents` branch has ALREADY made this pivot in code (this is not a fresh decision — this ADR ratifies and scopes the hardening of what exists). The simulation output is now the native engine type `ReactionNetworkProblem`, an AlgebraicAgents `@aagent` (`src/state.jl:40`) whose constructor lives at `src/solvers.jl:536`. It owns the instance lifecycle directly: an `ongoing_transitions::Vector{Transition}` heap, a per-step `evolve!`/`finish!` loop, and a single authoritative clock advanced at `src/solvers.jl:668` (`state.t += state.dt`). The old `DiffEqBase.DiscreteProblem` transform, `EnsembleProblem` path, and `optim.jl` are gone; `DiffEqBase`/`DifferentialEquations`/`OrdinaryDiffEq`/`NLopt`/`Catlab` were dropped from `Project.toml`.

## Decision

Keep the native discrete-event engine as the primary and only first-class simulation backbone. Do NOT reintroduce a SciML `DiscreteProblem` as the core. The native step (`AlgebraicAgents._step!(state::ReactionNetworkProblem)`, `src/solvers.jl:642`) is what AlgebraicAgents drives via `simulate`/`step!`.

SciML interop, if ever wanted, is demoted to an OPTIONAL, additive adapter (a package extension) for two narrow purposes only: (a) the mass-action CRN special case, and (b) coupling a continuous sub-dynamics as a sibling agent. It is explicitly out of scope for Milestone 1 and must never become a load-bearing dependency of the core.

## Consequences

Positive: the engine matches the domain semantics (stateful, resource-occupying instances); a single clock removes the dual-clock fragility that REVIEW.md flagged (#19), which is now MOOT; the dependency surface shrinks (no DiffEq/Catlab churn), directly serving the production-readiness and dependency-minimalism goals; and because the engine IS an AA agent, Workstream E (AA integration) is partly in-hand rather than bolted on.

Costs / obligations this ADR creates (tracked, not yet done):
- The engine currently lacks semantic tests. The Phase-0 gate (per the brief) is real `@test`s for conservation, non-negativity, capacity, and determinism-under-seed before any backbone refactor.
- Reproducibility is not yet wired: all stochastic draws use the global RNG (`src/solvers.jl:141,321,413`, `src/state.jl:123`). An `AbstractRNG` + `seed` must be threaded through the state and every draw (REVIEW.md #11). This ADR makes determinism a backbone requirement, formalized in [ADR 0002](0002-priority-weighted-allocation.md).
- Several engine bugs survived the pivot and are not chartered away by this decision — they are tracked in INVENTORY.md ("NEW defects introduced by the pivot" + the REVIEW.md reconciliation): undefined `q` in `free_blocked_species!` (`src/solvers.jl:512`), `add_to_spawn!` arity/`transHash` bug (`src/state.jl:251`), `o.val` vs `o.sampled` (`src/state.jl:137`), no-op `event_action!` (`src/solvers.jl:323`), and lifetime-terminated instances never pruned (`src/solvers.jl:501`).
- `@agentize` is exported but undefined (`src/interface/solve.jl:1`); the real contract is `prob = ReactionNetworkProblem(acs[, u0, p]; …)` then `simulate(prob[, n])`. Either implement `@agentize` as thin sugar over the constructor or delete the export.

## Considered and rejected

- SciML `DiscreteProblem` as core (the v0.2.x design): rejected — memoryless update cannot carry in-flight instance state; forces the `prob.p[:__state__]` smuggling pattern and a second clock.
- A hybrid where SciML owns the clock and the native engine is a callback: rejected — inverts ownership of the lifecycle/allocation logic that is the heart of the system, and re-imports the heavy DiffEq dependency into the core.
