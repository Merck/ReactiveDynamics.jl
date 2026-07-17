# ADR 0002 — Priority-weighted resource allocation under contention

- Status: Accepted (maintainer-confirmed 2026-06-17: "priority-weighted, suggest a cleaner implementation")
- Deciders: maintainer + rework
- Relates to: rework brief design-prior A ("the allocation rule under contention … the centerpiece of the contract"); [ADR 0001](0001-discrete-event-engine.md); REVIEW.md #13 (modalities). Algorithm verified numerically on 10 scenarios (see "Verification").

## Context

Each simulation tick, multiple transition instances compete for shared, finite resource supply `u[s]`. A transition's demand is **conjunctive (Leontief)**: an instance needs ALL of its required species simultaneously to proceed — partial allocation of one species without the others is wasted. The maintainer has fixed the policy as **priority-weighted**: when a resource is scarce, higher-priority transitions get proportionally more. The remaining design question — and the one this ADR answers — is *how* to compute that allocation cleanly and correctly.

### What exists today (the thing we are replacing)

`src/solvers.jl:16-128`: seven entangled functions — `get_reqs_init!` and `get_reqs_ongoing!` (two near-duplicate requirement builders), `get_allocs!` dispatching to `alloc_weighted!` / `alloc_greedy!`, and `get_frac_satisfied` / `get_init_satisfied` (two near-duplicate satisfaction reducers). Confirmed flaws:

1. **Not work-conserving / strands resource.** `alloc_weighted!` splits *each resource row independently* (share ∝ `demand × priority`), ignoring the conjunctive coupling during the split. Conjunctivity is then patched post-hoc by taking `fill[t] = minₛ alloc[s,t]/req[s,t]` and rescaling allocations *down* to that fraction (`get_frac_satisfied:99`). The over-allocated amount of the non-binding resource is **stranded, not redistributed** to other contenders who could use it. The allocation is not Pareto-efficient.
2. **Muddy priority semantics.** Weight is `demand × priority`, so a high-demand transition wins more *even at equal priority*; the number's meaning (rank? weight? rate?) is undocumented.
3. **Non-deterministic.** `alloc_greedy!` uses an unstable `sort` (`solvers.jl:76`); ties break unpredictably — compounding the absence of RNG seeding.
4. **Tangle.** Two reqs builders × two reducers × a dispatcher, dense index code, dead `state` args, misleading `!` (`alloc_weighted!` returns a fresh array), one global `state.p[:strategy]`.
5. **No invariant guards.** The `max(0, u[i])` patch (`solvers.jl:68`) hints negativity has occurred.

### Recognizing the problem

Conjunctive demands + weighted fair sharing of a scarce resource is exactly **weighted max-min fairness**, solved by **weighted progressive filling** (a.k.a. water-filling) — the same primitive behind weighted fair queueing and Dominant Resource Fairness. It is work-conserving and conjunctive-consistent *by construction*, deterministic, and needs **no solver dependency** (it is a simple iterative loop with a closed-form step size).

## Decision

Replace the seven-function tangle with a single **weighted progressive-filling allocator** over per-transition *fill fractions*, plus one parametrized requirements builder and a small integer-rounding wrapper for the spawn phase.

### Priority semantics (the contract)

`priority[t] ≥ 0` (schema attr `transPriority`) is the **fill-rate weight** in progressive filling: every active transition's fill fraction `f[t]` grows at a rate proportional to `priority[t]`. Consequences, documented and testable:
- On a single contended resource, supply splits in proportion to `priority[t] × req[s,t]` among the transitions still contending for it — i.e. at **equal demand, the split equals the priority ratio** (priorities 1 vs 3 → 1:3). This is the headline guarantee practitioners can reason about.
- `priority[t] = 0` means "only run from genuinely idle/leftover resource" (never advances during contention) — maintainer-confirmed reading. Equal priorities → weighted-equal (max-min fair) sharing.

**Priority is dynamic / time-varying (maintainer-confirmed 2026-06-17).** `transPriority` is a `SampleableAttributeT` (`src/ReactiveDynamics.jl:54`) — the same class as `transRate` — so it can be a constant scalar OR an expression evaluated per tick (e.g. a priority that responds to `@t()`, a budget observable, or any parameter). It flows through the same compiled-closure path as rates: the recipe is compiled once at build time (`compile_attrs`, `src/compilers.jl:145`) and evaluated per tick via `context_eval` (`src/state.jl:189`). The allocator consumes the resulting scalar vector and is itself agnostic to where the number comes from. For the agentic-authoring path, a constant priority is a plain serializable `Float64`; a time-varying one is an `ExprNode` like any other sampleable attribute — no special-casing.

**One semantics fix this implies:** priority must be evaluated **fresh per tick for BOTH phases**. Today the spawn phase reads a per-tick value (`state[:, :transPriority]`, `src/solvers.jl:158` → `context_eval`), but the ongoing phase reads a value snapshotted at spawn (`map(t -> t[:transPriority], …)`, `src/solvers.jl:240`, frozen by `get_sampled_transition` at `src/solvers.jl:6-11`). Under time-varying priority these diverge. Decision: priority is a **transition-level** weight, not an instance-level frozen value — all instances of transition `t` (spawning and in-flight) use the current `priority[t]` each tick. The ongoing call site must re-read the fresh per-tick priority rather than the spawn snapshot.

### Algorithm (unified core)

Inputs per call: `req[s,t] ≥ 0` (units of species `s` per unit fill of transition `t`), supply `u[s] ≥ 0`, weights `w[t] = priority[t] ≥ 0`, and a per-transition cap `fmax[t]` (the desired number of new instances for spawn, or `Inf`/desired progress for ongoing). Output: fill fractions `f[t]` and `alloc[s,t] = f[t]·req[s,t]`.

```
f[t]      ← 0
r[s]      ← u[s]                                  # remaining supply
active[t] ← (fmax[t] > 0) and (w[t] > 0) and (req[:,t] has a positive entry)
for t with fmax[t] > 0 and no positive demand:    # demand-free transitions
    f[t] ← fmax[t]                                # fill immediately, consume nothing
while any active:
    # 1. largest fill-step τ that neither overdraws a resource nor overshoots a cap
    dτ ← ∞
    for each species s:  D[s] ← Σ_{active t} w[t]·req[s,t]
                          if D[s] > 0:  dτ ← min(dτ, r[s] / D[s])
    for each active t:    dτ ← min(dτ, (fmax[t] − f[t]) / w[t])
    if dτ not finite: break
    # 2. advance all active fills and debit resources
    for each active t:    f[t] += w[t]·dτ
    for each species s:   r[s] = max(0, r[s] − D[s]·dτ)
    # 3. freeze transitions that hit their cap, or that need a now-saturated resource
    for each active t with f[t] ≥ fmax[t]:                 active[t] ← false; f[t] ← fmax[t]
    for each saturated species s (r[s] ≈ 0):
        for each active t with req[s,t] > 0:               active[t] ← false
alloc[s,t] ← f[t]·req[s,t]
```

This terminates in at most `S + T` iterations (each iteration freezes ≥1 transition or saturates ≥1 resource). It is work-conserving (a resource is only left idle once every transition that could use it is frozen), conjunctive-consistent (allocation is `f[t]·req[:,t]` by construction — nothing is stranded), and deterministic (no sort; ties resolved by the simultaneous freeze rule).

**Zero-priority "leftover only" — two-stage.** The single loop above excludes `w[t] = 0` transitions entirely (they never get the resource). To honor the confirmed `priority = 0` semantics ("runs only from genuinely leftover resource"), run the allocator in two stages: (1) progressive-fill the positive-priority transitions as above; (2) take the supply that remains and progressive-fill the zero-priority transitions among themselves at **equal weight** (max-min fair over the leftover). A zero-priority transition therefore advances iff positive-priority demand did not exhaust its required species, and never competes with positive-priority demand. Verified: with `u=10`, T1(prio 1, cap 4) and T2(prio 0) → `f=[4, 6]`; with T1 uncapped → `f=[10, 0]` (T2 correctly starved); two zero-priority transitions split the leftover equally. This generalizes to a clean **priority-tier** structure if ever needed (descending distinct priority values define tiers; each tier fills from the previous tier's residual) — but for Milestone 1 only the two tiers "positive" and "zero" are specified.

### Spawn phase (integer instances)

Run the core with `fmax[t] = q_desired[t]` to get the real-valued fair fill `f[t]`, then `n[t] = floor(f[t])`. Recover leftover supply `r = u − Σ_t req[:,t]·n[t]` and do a **deterministic priority-ordered integer top-up**: iterate transitions in order `(−priority[t], t)` (priority desc, then index asc for stable tie-break), and while a whole additional instance fits in `r`, grant it; repeat until no instance fits. This keeps spawn counts integral while staying work-conserving and deterministic. (Replaces `get_init_satisfied`.)

### Modality handling (one builder, not two)

A single `build_requirements!(req, transitions, qs, state; counted_modalities, dt_scale)` replaces both `get_reqs_init!` and `get_reqs_ongoing!`:
- **Spawn**: `counted_modalities` = "upfront" tokens only (exclude `:rate` and `:nonblock`), `dt_scale = 1`.
- **Ongoing**: include `:rate` tokens scaled by `state.dt` (only when `transCycleTime > 0`) and `:nonblock` tokens unscaled.

This is also the natural seam to land REVIEW.md #13: replace the ad-hoc `Set{Symbol}` modality with typed orthogonal fields, so "which tokens count and at what scale" is a documented function of the modality type rather than scattered `in(:rate, …)` checks. (Tracked separately; this ADR only fixes allocation.)

### Proposed Julia surface (replaces the 7-function set)

```julia
# Struct-of-arrays scratch reused across ticks (no per-tick heap churn):
struct AllocWorkspace
    req::Matrix{Float64}     # S×T, rebuilt each call
    f::Vector{Float64}       # T
    r::Vector{Float64}       # S
    active::BitVector        # T
    D::Vector{Float64}       # S
end

build_requirements!(ws, state, qs; counted, dt_scale)  # fills ws.req
progressive_fill!(ws, u, w; fmax)                        # core; fills ws.f
spawn_integer!(ws, u, w, q_desired)::Vector{Int}         # spawn wrapper
# alloc[s,t] = ws.f[t]*ws.req[s,t] derived on demand; caller does state.u .-= sum(alloc; dims=2)
```

Deleted: `get_reqs_init!`, `get_reqs_ongoing!`, `get_allocs!`, `alloc_weighted!`, `alloc_greedy!`, `get_frac_satisfied`, `get_init_satisfied`. The `state.p[:strategy]` `:weighted`/`:greedy` switch is removed — strict-greedy is just the `priority → ∞` limit of weighted filling and is no longer a separate code path. (If a strict-priority lexicographic mode is ever wanted, it becomes a documented flag on the single allocator, not a duplicate function.)

### Integration into `evolve!` (`src/solvers.jl:133-313`)

- Spawn call site (`solvers.jl:156-160`): `get_reqs_init!` + `get_allocs!` + `get_init_satisfied` → `build_requirements!(…; counted=upfront, dt_scale=1)` then `spawn_integer!`.
- Ongoing call site (`solvers.jl:235-243`): `get_reqs_ongoing!` + `get_allocs!` + `get_frac_satisfied` → `build_requirements!(…; counted=[:rate,:nonblock], dt_scale=state.dt)` then `progressive_fill!`. **Priority must be re-read fresh here** — `state[t.i, :transPriority]` (per-tick `context_eval`) rather than the spawn-time snapshot `t[:transPriority]` (`solvers.jl:240`), so time-varying priority applies to in-flight instances too.
- `state.u .-= sum(alloc; dims=2)` and the `:allocation`/`:valuation_cost` log rows are unchanged in spirit.
- All `rand` draws (the Poisson spawn count upstream, Binomial PoS at `finish!`) route through the state's `AbstractRNG` per [ADR 0001](0001-discrete-event-engine.md)'s reproducibility obligation; the allocator itself is deterministic and RNG-free.

## Verification (worked numeric examples)

Executed the algorithm independently (Python reference implementation) on 10 scenarios; all invariants held.

| Scenario | Setup | Result | Checks |
|---|---|---|---|
| Priority split | 1 resource u=10, two trans need 1 each, w=1 vs 3 | f=[2.5, 7.5] | ratio exactly 3:1 ✓ |
| Conjunctive 2-res | req T1=[1,2], T2=[1,1], u=[10,10], w=1 | f=[3.33, 3.33], both resources fully used | no stranding ✓ |
| Work-conservation | u=10, w=1, T1 capped at fmax=2 | f=[2, 8], total used 10/10 | leftover 3 flows to T2 (naive split would strand it, using only 7) ✓ |
| Cross-resource redistribute | X shared scarce (u=4), Y ample, T1=[1,0] T2=[1,1] | f=[2,2] | T1 freezes on X, capacity respected ✓ |
| Differing binding res | T1=[2,1], T2=[1,2], u=[10,10], w=1 | f=[3.33,3.33], both resources at 10/10 | symmetric, fully used ✓ |
| Weighted + conjunctive | T1=[1,1] T2=[1,1], X scarce u=6, w=1 vs 2 | f=[2,4], X used 6/6 | fill ratio exactly 2:1 ✓ |
| Integer spawn | u=7, each instance needs 2, want 5 each, eq prio | n=[2,1], leftover 1 | integral; real fair fill was 1.75 each; top-up deterministic ✓ |
| 3-way stranding | u=12, w=1, T1 capped at 2 | f=[2,5,5], total 12/12 | freed amount split evenly among absorbers ✓ |
| Determinism | repeated runs, mixed weights | identical every run | reproducible ✓ |
| Zero supply | u=0 | f=[0,0], alloc=0 | non-negativity, no div-by-zero ✓ |

## Invariant test plan (Phase-0 semantic tests)

- `@test all(alloc .>= 0)` and `@test all(sum(alloc; dims=2) .<= u .+ tol)` — non-negativity & capacity.
- `@test alloc ≈ f' .* req` (elementwise) — conjunctive consistency (no stranding).
- `@test` on the priority-split scenario: `f[2]/f[1] ≈ priority[2]/priority[1]` at equal demand.
- `@test` work-conservation: after allocation, no species has leftover `> 0` while an unfrozen transition still has positive unmet demand for it.
- `@test` determinism: two runs with identical inputs (and seed) give identical `f`/`alloc`.
- `@test` spawn integrality: `spawn_integer!` returns `Int`s and `sum(req[:,t]·n[t]) <= u`.

## Dependencies

**None.** Progressive filling is a closed-form iterative loop; no JuMP/optimization-solver dependency. This satisfies the dependency-minimalism requirement.

## Considered and rejected

- **Keep per-resource independent split + min-rescale (current code):** rejected — not work-conserving (strands resource), muddy priority semantics, non-deterministic greedy path. This is the thing being replaced.
- **LP/convex program** maximizing `Σ priority[t]·f[t]` (linear) or `Σ priority[t]·log f[t]` (Nash/proportional-fair): rejected as the *implementation*, though instructive. The linear objective degenerates to strict greedy-by-priority (no weighted sharing); the log objective is genuinely weighted-fair but requires a solver dependency for a problem that water-filling solves in closed form. Progressive filling gives the weighted max-min fair point directly with no dependency.
- **Strict greedy-by-priority** (`alloc_greedy!`): rejected as default — starves lower-priority transitions entirely and is the `priority → ∞` degenerate limit of the chosen rule, so it needs no separate code path. Retained only as a possible documented flag if a use case demands lexicographic priority.

## Resolved (maintainer, 2026-06-17)

1. **Priority scope → dynamic.** `priority[t]` may be time-varying, like any other parameter (a function of time / observables), evaluated per tick. It is a `SampleableAttributeT` and flows through the same compiled-closure path as rates; the allocator is agnostic to its origin. This obligates the ongoing-phase fresh-read fix noted above (in-flight instances use the current transition priority, not a frozen one).
2. **Cross-tick fairness → per-tick.** Per-tick weighted max-min fairness is sufficient; no carry-over/starvation-memory mechanism for Milestone 1. The allocator stays memoryless and stateless across ticks. Revisit only if scenarios surface starvation artifacts.
3. **`priority = 0` → "runs only from leftover."** Confirmed: a zero-priority transition never advances during contention; it proceeds only from genuinely idle/leftover resource. (Not treated as an epsilon weight.)
