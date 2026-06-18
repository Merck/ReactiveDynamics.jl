# ReactiveDynamics Modeling Contract

> **DRAFT — Phase-0 modeling contract.** The fork-independent sections (§1–§5) are complete. The ADR 0003 data-store decision is now ACCEPTED (drop ACSets for a dependency-free typed IR; promote the reactant relation), with serialization in [ADR 0005](adr/0005-serialization-json-ir.md) (single JSON + ExprNode) and runtime mutation in [ADR 0004](adr/0004-runtime-mutation.md) (append-only). The remaining sections — object model, composition semantics, serialization schema — can now be written against that decision and are the next contract increment (see "Pending sections" at the end).

## Orientation

ReactiveDynamics (RD) is a timed, stochastic, resource-constrained Petri net — a discrete-event system in which transitions periodically spawn in-flight instances that draw on shared, finite resource pools (species) over a service duration, then emit products and return or consume those resources on completion. It is *not* a chemical reaction network and *not* a continuous-time Gillespie process: it is a fixed-step discrete-time engine with priority-weighted resource allocation. This contract pins the *operational semantics* of that engine — the time model, the per-tick firing and lifecycle rules, the modality (resource-claim) semantics, the determinism/seeding obligations, and the typed attribute domains — independently of how the static model is stored. Concrete `file:line` references point at the current `ref-agents` engine and are evidence of present behavior (sometimes evidence of a *bug*, flagged as such), not a commitment to keep that behavior. The object model, composition semantics, and serialization schema are deferred to a later revision pending ADR 0003.

The contract is organized as follows: §1 Modality Truth Table (the resource-claim semantics that the allocator and lifecycle depend on); §2 Time Model (single discrete clock, spawn intensity, cycle/lifetime, dt-invariance); §3 Operational Semantics (instance lifecycle, the ordered step, and the invariants an implementation must uphold); §4 Determinism & Seeding Contract; §5 Attribute Contract (typed attribute domains). ADR 0002 (priority-weighted progressive filling) is the normative specification for the allocator referenced throughout.

---

## 1. Modality Truth Table

A *modality* is the contract attached to each LHS (consumed) token of a transition. It answers three independent questions: **when** is the resource claimed against the pool, **whether** it is returned when the cycle completes, and **whether** the claim blocks the resource while the cycle runs. Today these three concerns are conflated into an unvalidated `Set{Symbol} ⊆ {:nonblock, :conserved, :rate}` (the allowed universe `species_modalities`, `ReactiveDynamics.jl:144`) whose cross-product the solver interprets inconsistently and never documents. This section re-models them as three orthogonal, closed, typed fields and gives the full legal truth table.

### 1.1 The three orthogonal axes

Each consumed token declares a `Modality` value with exactly these three fields. All combinations are legal unless explicitly marked illegal below.

| Field | Type (closed set) | Meaning | Default |
|---|---|---|---|
| `allocation` | `{upfront, perstep}` | **When** the resource is claimed. `upfront`: the full `q·stoich` is reserved once, at spawn time, as a precondition for the cycle to start at all. `perstep`: a slice is claimed each tick while the cycle runs (a flow/throughput draw rather than a one-time seizure). | `upfront` |
| `return` | `{consumed, conserved}` | **Whether** the claimed amount returns to the pool when the cycle finishes. `consumed`: permanently destroyed (true consumption). `conserved`: the reserved amount is credited back at finish (a loan / temporary hold). | `consumed` |
| `blocking` | `{block, nonblock}` | **Whether** the claim holds the resource for the duration. `block`: the resource is unavailable to other transitions until finish. `nonblock`: the resource is touched but immediately released every step, so it never actually constrains anyone (a soft / advisory draw). | `block` |

These map onto the legacy tags as: empty set ⇒ `(upfront, consumed, block)`; `:rate` ⇒ `allocation = perstep`; `:conserved` ⇒ `return = conserved`; `:nonblock` ⇒ `blocking = nonblock`. Note the legacy tags are *not* one-axis-each in isolation — `:nonblock` in the current code simultaneously implies per-step reservation *and* non-blocking release *and* (forced) non-conservation, which is exactly why it cannot be combined with `:conserved`. The re-model separates these so the constraint becomes a single clean rule (below) rather than an ad-hoc `error`.

### 1.2 Effect dimensions each combination resolves to

Every legal `Modality` resolves deterministically to three solver effects:

- **Allocation timing** — *spawn-reserve* (counted by `get_reqs_init!`, must be satisfiable for the transition to spawn) vs *per-step-reserve* (counted by `get_reqs_ongoing!` each tick, optionally `dt`-scaled).
- **Consumption / return semantics** — does `finish!` credit `state.u` back, and by how much.
- **Ledger effect** — what hits the period cost line (every reserved unit is valued at `specCost` in the `:valuation_cost` log entry) and whether it is later offset by a reward/return.

### 1.3 Full truth table (every legal combination)

`s = stoich`, `q = transition multiplicity`, `Δt = state.dt`, `C = transCycleTime`. "Reserve" = subtracted from `state.u` at the indicated time and charged to period cost at `specCost`.

| # | `allocation` | `return` | `blocking` | Legacy `Set{Symbol}` | Allocation timing | Return at finish | Ledger effect | Plain meaning |
|---|---|---|---|---|---|---|---|---|
| 1 | `upfront` | `consumed` | `block` | `{}` (empty) | Reserve `q·s` once at spawn (must fit) | nothing | Net spend of `q·s` valued at `specCost`; gone for good | **Raw consumption.** Classic stoichiometric input — burned to start the cycle. |
| 2 | `upfront` | `conserved` | `block` | `{:conserved}` | Reserve `q·s` once at spawn (must fit) | `+ q·s` | Cost charged at spawn, fully credited back at finish; net-zero on success | **Held capital / equipment.** Seized for the duration, returned intact. |
| 3 | `perstep` | `consumed` | `block` | `{:rate}` | Reserve `q·s·Δt` each tick (only if `C > 0`) | nothing | Continuous spend; period cost accrues per tick | **Metered consumption (flow).** A burn *rate* drawn down over the cycle. |
| 4 | `perstep` | `conserved` | `block` | `{:rate, :conserved}` | Reserve `q·s·Δt` each tick (only if `C > 0`) | `+ q·s·C` | Per-tick cost accrued, then credited back `q·s·C` at finish | **Rented throughput.** A rate-based hold, fully returned at completion. |
| 5 | `perstep` | `consumed` | `nonblock` | `{:nonblock}` | Reserve `q·s` each tick, **freed every step** by `free_blocked_species!` | nothing (already freed) | Touched/measured but not held; no net pool change between ticks | **Soft / advisory draw.** Reads the resource each step without contending for it. |

Rows 1–4 are the four corners of the `allocation × return` plane under `block`. Row 5 is the single coherent `nonblock` form. Every legacy `Set{Symbol}` value in the wild maps to exactly one of these five rows.

### 1.4 Illegal combinations

| `allocation` | `return` | `blocking` | Legacy form | Why illegal |
|---|---|---|---|---|
| any | `conserved` | `nonblock` | `{:nonblock, :conserved}` | Conservation means "credit the held amount back **at finish**"; non-blocking means "release it **every step**." A resource cannot be both held-until-finish and continuously-released. The current code rejects this explicitly (`error` at `solvers.jl:461-465`); the re-model makes it a single validation rule: **`blocking = nonblock` requires `return = consumed`.** |
| `perstep` | any | any (structured species) | `{:rate, ...}` on a structured token | `:rate` is unsupported for structured/agent species (`solvers.jl:38-42`): you cannot reserve a fractional, `dt`-scaled slice of an indivisible agent. Validation rule: **`allocation = perstep` requires a non-structured (countable) species.** |
| `perstep` (`return = consumed`) when `C = 0` | — | — | `{:rate}` with `transCycleTime = 0` | Per-step reservation only fires when `C > 0` (`get_reqs_ongoing!:36`). With `C = 0` a `perstep` token silently reserves nothing — a foot-gun, not a meaning. Validation rule: **`allocation = perstep` requires `transCycleTime > 0`.** |

All other combinations are legal. This collapses the previous 2³ = 8 implicit tag subsets (most of them undefined behavior) into **5 legal rows + 1 illegal rule**, all named and validated at construction.

### 1.5 Connection to ADR 0002 (`build_requirements!`)

ADR 0002's weighted progressive-filling allocator consumes a single requirements matrix `reqs[species, transition]`; the modality axes are precisely what decide **which tokens land in which requirements pass and with what coefficient**. The two legacy entry points `get_reqs_init!` / `get_reqs_ongoing!` unify into one parameterized builder:

```
build_requirements!(reqs, qs, model; counted_modalities, dt_scale)
```

- **`counted_modalities`** selects which `allocation` value this pass reserves for:
  - spawn pass ⇒ `counted_modalities = {allocation = upfront}` (rows 1, 2). This is exactly today's "exclude `:rate` and `:nonblock`" filter (`get_reqs_init!:20`, `get_init_satisfied:115`), now stated positively.
  - ongoing pass ⇒ `counted_modalities = {allocation = perstep}` (rows 3, 4, 5) — i.e. today's `:rate`/`:nonblock` branch (`get_reqs_ongoing!:35,43`).
- **`dt_scale`** is the per-tick coefficient applied to the `perstep` reservation:
  - `dt_scale = Δt` for `return = consumed` *or* `conserved` flow tokens (rows 3, 4 — today's `:rate` path, `solvers.jl:37`), gated on `C > 0`.
  - `dt_scale = 1` for `nonblock` tokens (row 5 — today's unscaled `:nonblock` path, `solvers.jl:43`).

The `return` axis is orthogonal to allocation and is consumed only by `finish!`: `return = conserved` credits `state.u` back by `q·s` (row 2) or `q·s·C` (row 4, the `perstep` flow integral) — exactly the existing `(in(:rate) ? transCycleTime : 1)` factor at `solvers.jl:438-442`, now derived from the typed fields rather than from tag intersection. The `blocking` axis is consumed only by `free_blocked_species!`: `nonblock` tokens are freed each step (`q·s`), and once the undefined-`q` bug (`solvers.jl:512`) is fixed this must credit back exactly what the ongoing pass reserved for that token (`q·s` with `dt_scale = 1`), keeping the per-step reserve/free pair conservative.

Net result: one requirements builder, two call sites distinguished only by `(counted_modalities, dt_scale)`, a single `finish!` return rule keyed on `return`, a single `free` rule keyed on `blocking`, and a closed 5-row truth table that the allocator never has to special-case.

---

## 2. Time Model

This section specifies how the model represents and advances time. It is written against the discrete-event engine that executes the spec; symbol names refer to authoring keywords and their compiled runtime fields. Source references are to the current `ref-agents` tree and are load-bearing — behavior described here is what the code does, not what it should do; known hazards are called out explicitly.

### 2.1 Single discrete clock

The model has exactly one clock: a scalar simulation time `t` (runtime field `state.t::Float64`, `state.jl:48`). There is no per-transition or per-agent local clock; ongoing transitions carry an *age* (`Transition.t`, the spawn time) and a *progress accumulator* (`Transition.state`), but both are read against the single global `t`.

The clock advances by a fixed step `dt` once per tick, as the final action of the step loop: `state.t += state.dt` (`solvers.jl:668`). All within-tick work — freeing blocked resources, updating observables, sampling/spawning new transition instances, advancing ongoing instances, finishing completed instances, firing events, logging the valuation — happens *at the pre-increment value of `t`* and is treated as simultaneous (`_step!`, `solvers.jl:642-672`). The simulation is therefore a fixed-step discrete-time process, not a continuous-time / next-event (Gillespie-style) process: there is no event queue that jumps `t` to the next firing time.

- **Units.** `t` and `dt` are in *time units* (the abstract unit of the model). The constructor normalizes everything by a unit `tunit` (default `oneunit(tspan)`, i.e. `1.0`), so internally `t` runs on a dimensionless `[0, tspan]` axis (`get_tcontrol`, `solvers.jl:526-533`).
- **Default `dt`.** `1.0` time unit (`tunit`) when neither `dt` nor `tstops` is supplied (`solvers.jl:531`).

### 2.2 `tspan` / `dt` / `tstops` resolution

Time control is resolved once at construction by `get_tcontrol(tspan, args)` (`solvers.jl:526-534`):

1. **`tspan`** may be given as a duration (a number) or as a `(t0, t1)` tuple; a tuple is reduced to its length `t1 - t0` (`:527`). The simulation axis is then always normalized to start at `0.0`: the returned span is `(0.0, tspan/tunit)` (`:529, :533`). A nonzero `t0` does not shift the clock — only the duration is honored.
2. **`tunit`** defaults to `oneunit(tspan)` and divides both `tspan` and `dt` (`:528-531`).
3. **`dt`** resolution, in precedence order (`:531`):
   - explicit `dt` keyword → used directly (divided by `tunit`);
   - else, if `tstops` is given → `dt = tspan / tstops` (i.e. `tstops` = number of equal steps to cover the span);
   - else → `dt = tunit` (i.e. `1.0`).

   The number of ticks executed is `tspan / dt`, since `_projected_to` reports completion once `state.t > state.tspan[2]` (`solvers.jl:676`).

**Naming hazard (`dt` vs `tstep`).** The runtime struct field is `dt` (`state.jl:53`), and every consumer reads `state.dt` (`solvers.jl:37, :261, :668`). The constructor, however, stores the resolved step under the keyword key `:tstep` (`keywords[:tspan], keywords[:tstep] = get_tcontrol(...)`, `solvers.jl:552`) and then fills the `dt` field *positionally* from `get(keywords, :tstep, 1)` (`solvers.jl:603`). So the spec keyword is effectively `dt` (or `tstops`), the internal keyword bag uses `tstep`, and the live field is `dt`; the three are reconciled by argument position, not by name. **Contract:** authors set `dt` (or `tstops`); `tstep` is an internal alias that should not be relied upon and is a candidate for renaming.

### 2.3 Spawn intensity (`rate` / `transRate`) — expected instances per unit time

Each transition has a spawn rate, authored as `rate` (aliases: see `prettynames[:transRate] = [:rate]`, `ReactiveDynamics.jl:105`) and stored in the `transRate` column (`ReactiveDynamics.jl:55`). There is no schema default for `transRate`; it is supplied positionally per transition at authoring time (`merge_acs!`, `create.jl:202`).

The rate is interpreted as a **Poisson intensity = expected number of new instances per unit time**. At parse time, `expand_rate` (`create.jl:149-163`) rewrites the authored rate expression `r` into

```julia
rand(Poisson(max(state.dt * r, 0)))
```

so each tick draws the number of fresh instances from `Poisson(dt · r)` (`create.jl:151`). Because `E[Poisson(dt·r)] = dt·r`, the expected spawn count per tick scales linearly with `dt`, and the long-run expected spawn rate is `r` per unit time — independent of `dt`. The `max(·, 0)` guards against negative authored rates yielding an invalid Poisson parameter.

**Rate qualifiers.** Two qualifier macros are recognized inside the rate expression by `expand_rate`:

- **`@deterministic <expr>`** — bypasses the Poisson draw entirely; the spawn count per tick is the bare `<expr>` (no `dt` scaling, no randomness) (`create.jl:150, :153`). Use this when the count is prescribed rather than sampled.
- **`@ct <expr>` / `@cycletime <expr>`** (aliases `prettynames[:transCycleTime] = [:ct, :cycletime]`, `ReactiveDynamics.jl:113`) — rewrites to `1 / <expr>` in-place (`create.jl:157-161`), letting a rate be written as a reciprocal cycle time.

> **Spec note — no `@per_step`.** There is no `@per_step` qualifier in the source (verified across `src/`); the only rate qualifiers are `@deterministic` and `@ct`/`@cycletime`. If a "per-step" (i.e. dt-independent, count-per-tick) semantics is wanted, it is currently expressed via `@deterministic`. A dedicated `@per_step` is a possible future addition but must be documented as such, not assumed.

**Discretization hazard — `ceil` on the spawn count.** In `evolve!`, the per-transition target count is computed as `qs[i] = transRate * transMultiplier` and then `qs .= ceil.(Int, qs)` (`solvers.jl:141-144`). When `transRate` is itself a sampled/numeric value (rather than the Poisson-expanded expression), this `ceil` rounds *up* to the next integer every tick. Combined with halving `dt`, this biases spawning upward: a small expected count like `0.3` becomes `1` each tick, and doing twice as many ticks roughly *doubles* the spawned total instead of preserving it. **Therefore: halving `dt` does NOT in general preserve expected dynamics through this code path.** The Poisson path (`expand_rate`) is dt-invariant in expectation; the `ceil` path is not. Authors who need dt-refinement invariance should use rate expressions that flow through `expand_rate`'s Poisson draw and avoid relying on fractional `transRate` values that get `ceil`-ed.

- **`transMultiplier`** (`ReactiveDynamics.jl:62`, default `1`, `:124`) is a plain integer multiplier on the per-tick target count (`solvers.jl:141`). Units: dimensionless.

### 2.4 Cycle time (`cycletime` / `transCycleTime`) — service duration

`transCycleTime` (column `ReactiveDynamics.jl:56`; aliases `@ct`/`@cycletime`/`cycletime`, `:113`; **default `0.0`**, `:122`) is the amount of *accumulated, resource-weighted progress* an instance must reach before it can complete. Units: time units.

An ongoing instance carries a progress accumulator `Transition.state` (initialized to `0.0` at spawn, `solvers.jl:187`). Each tick it advances by the fraction of its resource demand that was satisfied this tick, times `dt`:

```julia
transition.state += qs[i] * state.dt        # solvers.jl:261
```

where `qs[i] ∈ [0,1]` is the saturation (fraction satisfied) of that instance (`get_frac_satisfied`, `solvers.jl:95-103`). An instance is eligible to finish once `Transition.state >= transCycleTime` (`finish!`, `solvers.jl:409, :412`; survivors retained by `filter!(s -> s.state < s[:transCycleTime], …)`, `:501`).

- **`transCycleTime == 0.0` (default):** the instance completes on the very tick it is created — i.e. instantaneous, single-tick transitions. Note also that `:rate`-modality resource draw for ongoing instances is gated on `transCycleTime > 0` (`solvers.jl:36`), so a zero cycle time means no continuous (rate-scaled) resource consumption.
- **`transCycleTime > 0`:** the instance persists across ticks, accumulating progress at a rate throttled by resource availability; under full saturation (`qs ≡ 1`) it completes after `ceil(transCycleTime / dt)` ticks.

Because completion is tested per tick at the post-step value of `state`, the *realized* service duration is quantized to a whole number of `dt` steps; the residual `state - transCycleTime` overshoot is not carried over. Finer `dt` reduces this quantization error but, per the `ceil` hazard above, does not by itself guarantee invariance of spawn counts.

### 2.5 Max lifetime (`maxlifetime` / `transMaxLifeTime`) — timeout

`transMaxLifeTime` (column `ReactiveDynamics.jl:59`; aliases `@lifetime`/`@maxlifetime`/`@maxtime`/`@timetolive`, `:114`; **default `Inf`**, `:123`) is a wall-clock age cap. Units: time units.

In `finish!`, an instance is forced to terminate once its age exceeds the cap — `(state.t - trans_.t) >= transMaxLifeTime` (`solvers.jl:408`). Lifetime and cycle time interact: an instance terminates on whichever fires first. If it ages out *before* reaching `transCycleTime`, it terminates with success count `q = 0` (no RHS products emitted), because the success draw is gated on `trans_.state >= transCycleTime` (`solvers.jl:412-416`). If it reaches cycle time, the success count is `Binomial(q, transProbOfSuccess)` (`solvers.jl:413`). Default `Inf` means "no timeout": completion is governed solely by cycle time.

### 2.6 Continuous rate ↔ discrete `dt` relationship

The model is a fixed-step Euler-style discretization of an underlying continuous-time intensity model:

- **Spawning** is dt-invariant *in expectation* on the Poisson path: `Poisson(dt·r)` summed over `tspan/dt` ticks has expectation `r·tspan` regardless of `dt`. Refining `dt` refines the *timing resolution* of spawns (and reduces the chance of multiple spawns being collapsed into one tick) without changing the expected total — **provided** the rate flows through `expand_rate` and is not subject to the `ceil` path.
- **Service / resource consumption** is dt-invariant in the limit: progress accrues as `∫ saturation dt` (discretized as `Σ qs·dt`), and `:rate`-modality resource draw is scaled by `dt` (`solvers.jl:37`). Halving `dt` doubles the number of consumption events while halving each, preserving the integral up to quantization.
- **Hazards that break dt-invariance:** (1) the `ceil` on fractional spawn targets (`solvers.jl:144`), which biases spawns upward as `dt → 0`; (2) cycle-time / lifetime quantization to whole `dt` steps; (3) any author-supplied per-step action or rate that is written as a count-per-tick rather than a rate-per-unit-time. **Recommendation:** treat `dt` as a numerical resolution parameter, validate that results are stable under `dt` refinement for the model at hand, and prefer rate-per-unit-time authoring (Poisson path) over count-per-tick.

### 2.7 Determinism contract for time

- **The clock is fully deterministic.** Given `tspan`, `dt`/`tstops`, and `tunit`, the sequence of tick times `t = 0, dt, 2dt, …` and the number of ticks (`tspan/dt`) are fixed and reproducible (`get_tcontrol`, `solvers.jl:526-534`; `state.t += state.dt`, `:668`; `_projected_to`, `:676`). Time does not depend on the RNG and is identical across runs and across seeds.
- **What is stochastic is the *content* of each tick, not its timing.** Randomness enters only through draws evaluated at the fixed tick boundaries: spawn counts `rand(Poisson(dt·r))` (`create.jl:151`), success counts `rand(Binomial(q, p))` (`solvers.jl:413`), event multiplicities `rand(Poisson(v))` (`solvers.jl:321`), and any sampled attributes. Under a fixed RNG seed the entire trajectory — including all draws — is reproducible.
- **`_reinit!` restores the clock** to `state.tspan[1]` (i.e. `0.0`) and clears ongoing transitions, log, and solution, so a re-run from the same seed reproduces the run exactly (`solvers.jl:619-628`).
- **Consequence for the contract:** time advancement is a guaranteed invariant — implementations may change the data store, the allocation strategy, or the RNG stream, but MUST preserve (a) a single scalar clock advanced by a constant `dt`, (b) `tspan/dt` total ticks over a `[0, tspan]` axis, and (c) the property that timing is RNG-independent while per-tick draws are seed-reproducible. (The full RNG/seeding obligations are specified in §4.)

### 2.8 Genesis modes — is Poisson the right primitive for business processes?

This section answers the maintainer's question directly: **Poisson-per-tick is the right primitive for one regime (genuine memoryless exogenous inflow) but the wrong DEFAULT for a business/R&D pipeline.** A drug pipeline is a near-deterministic chain `Phase1 → Phase2 → Phase3 → Market` where a program advances to the next phase *because* the prior phase succeeded — a token flow, not an independent random arrival. Forcing that flow through a memoryless Poisson clock injects spurious variance into rNPV and makes the pre/post-acquisition counterfactual noisy, defeating the comparison the demo exists to make. The good news, established by reading the engine: the better mechanism ALREADY EXISTS — genesis is internally two stages, and the second stage already provides token-flow and capacity gating. No engine redesign is needed; the contract's job is to NAME the regimes and validate intent.

**Genesis is already two stages (verified).** Each tick, for transition `t`: (1) a spawn-count PROPOSAL from the rate expression — `Poisson(dt·rate)` by default (`create.jl:151`) or a bare count under `@deterministic` (`create.jl:153`), times `transMultiplier`, then `ceil`'d (`solvers.jl:140-144`); then (2) an upfront-LHS RESOURCE GATE that clamps the realized count to `floor(allocation / stoich)` over the transition's upfront-consumed LHS tokens, via `get_reqs_init!`→`get_allocs!`→`get_init_satisfied` (`solvers.jl:156-160`, `:110-128`). The gate's `reqs == 0 ⇒ Inf` rule (`solvers.jl:121`) means a transition with NO upfront LHS keeps its full proposal (a pure source), whereas a transition WITH upfront LHS can only spawn as many instances as its input tokens and allocated resources permit (token-flow / capacity-limited). **So flow-triggered and capacity-limited genesis work today with zero new mechanism** — the `toy_pharma_model` already demonstrates both (phase hand-off via `candidate_compound` on discovery's RHS and `dx2market`'s LHS; capacity-limited starts via `3*@conserved(scientist) + @rate(budget) --> candidate_compound`).

**Source vs routing (the latent, correct distinction).** The split between exogenous arrivals and internal routing is already implicit in the data: an EMPTY LHS (`extract_reactants` returns `[]`, `reaction_parser.jl:50-51`) = a SOURCE that bypasses the gate and spawns at the proposal intensity; a NON-EMPTY upfront LHS = ROUTING that fires only when its input tokens exist. The contract makes this explicit rather than accidental.

**The closed genesis-mode tag.** Add ONE append-only transition field `genesis ∈ {poisson, scheduled, flow, capacity}` (a 4-value string tag — eval-free, JSON-Schema-enumerable per ADR 0005; append-only per ADR 0004), an intent declaration over the single shared two-stage execution path:

| Mode | Business meaning | Maps to (existing mechanism) | Source/Routing | New engine surface |
|---|---|---|---|---|
| `poisson` | Random independent exogenous inflow (unsolicited inbound deals) | Today's `expand_rate` `Poisson(dt·rate)` (`create.jl:151`) | Source | none — current default path |
| `scheduled` | Prescribed/calendar genesis (quarterly gate, budget cycle, planned start, the BD acquisition lever) | `@deterministic` bare-count (`create.jl:150,153`) + a documented calendar idiom `@scheduled(period,N) ⇒ @deterministic(N*periodic(period))`, since `periodic(state,period)` already exists (`state.jl:239`) and `period==0.0` returns true (one-shot at t0) | Source | none — documented idiom only |
| `flow` | Token-triggered: a program advances because the upstream phase succeeded | The EXISTING upfront-LHS gate (`solvers.jl:110-128`); author writes a high/`Inf` nominal rate with the upstream species as an upfront-consumed LHS reactant | Routing | none |
| `capacity` | Start as many as resources/headcount allow | Same gate path as `flow`, bounded by a renewable `@conserved`/`@rate` pool + `transCapacity` (`solvers.jl:147-153`) + the ADR-0002 allocator | Routing | none |

`batch` is a PARAMETER, not a mode (`batch::Int`, default 1, applied as post-gate rounding `qs[i] = batch * fld(qs[i], batch)`), and finite-population/one-shot is EXPRESSED (a flow source consuming a finite `specInitVal` pool; one-shot = `scheduled` at t0), not a new mode. Three of the four modes map onto existing mechanisms with zero new engine surface; only the `scheduled` calendar form needs a documented idiom (and `periodic` already compiles, `compilers.jl:67`).

**Validation rules (construction-time).** `flow`/`capacity` REQUIRE ≥1 upfront-consumed LHS species (else they would silently behave as an unbounded Poisson/deterministic source); `poisson`/`scheduled` with a consuming LHS should WARN (the author probably meant routing). This makes the source/routing intent checkable for agentic authoring.

**Rejected (out of scope).** The heavier queueing-lens proposal — redefine routing to drop the Poisson proposal and substitute a `min(1, dt·rate)` service-fraction primitive, splitting genesis into `ArrivalProcess` subtypes + a routing primitive — is rejected: it is a semantic break adding a second intensity code path for marginal fidelity gain, and the existing `rate=Inf + upfront LHS` idiom already realizes token-bounded firing `min(proposal, tokens) = tokens`. The one real artifact it identifies (a routing token left unserved because a small Poisson draw thinned it) is eliminated by the canonical `flow` idiom (`rate=Inf` makes the proposal non-binding), so it is a template/documentation fix, not an engine fix.

**Correctness blockers (independent of genesis design).** Two pre-existing bugs gate the business-process use cases: `event_action!` is a no-op (`solvers.jl:323` fetches `:eventAction` but never evaluates it — §3.4 Invariant 7), so the acquisition lever must use the `scheduled` rate-expression idiom, NOT the event channel, until repaired; and `add_to_spawn!` is doubly broken (`state.jl:251-256` — §3.4 Invariant 3), so `capacity`/`batch` genesis under-delivers under sustained over-demand. Also, because the `ceil` at `solvers.jl:144` breaks dt-invariance for non-Poisson counts (§2.3), `scheduled`/`batch` counts MUST be integer-valued (or `ceil` must be conditioned to the `poisson` path only).

---

## 3. Operational Semantics

This section specifies the per-tick firing and lifecycle rules of the model. It is written against the abstract spec — the *model* is whatever static object the authoring layer produces; the *engine state* is the live runtime object derived from it — and holds regardless of how that model is stored. Concrete `file:line` references point at the current `ref-agents` engine (`src/solvers.jl`, `src/state.jl`) where the rule is implemented or, where flagged, *violated*.

### 3.1 Vocabulary

- **Species** `s`: a resource pool holding a non-negative quantity `u[s]`. A species may be *plain* (a `Float64` count) or *structured* (backed by individual token agents whose live count is reflected into `u[s]`).
- **Transition** `t`: a stateful *recipe* — not an event — that periodically spawns instances. It carries a spawn `rate`, a `priority` (fill-rate weight, see ADR 0002), a `cycleTime`, a `probOfSuccess` (PoS), a `capacity`, a `maxLifeTime`, a `multiplier`, optional pre/post actions, and an LHS/RHS reactant specification with per-token `stoich` and `modality`.
- **Instance**: one in-flight firing of a transition (engine type `Transition`, `state.jl:14-26`), holding its parent recipe index `i`, a frozen sampled-attribute snapshot `trans`, its birth time `t`, its multiplicity `q` (how many concurrent firings this instance object represents), and an accumulated progress `state`.
- **Modality**: a per-LHS-token tag governing *when* a token is debited and *whether/how* it is returned (fully specified in §1; the legacy `Set{Symbol}` form, `{:nonblock, :conserved, :rate}`, is currently unvalidated, `ReactiveDynamics.jl:144`).
- **Tick**: one advance of the simulation clock by `dt`.

### 3.2 Instance lifecycle

Every transition instance passes through the following stages exactly once, in order:

1. **Genesis (Poisson spawn).** Each tick, transition `t` proposes `q_desired = ceil(rate · multiplier)` new firings, where `rate` is authored as `rand(Poisson(max(dt · rate, 0)))` (`create.jl:149-151`, `expand_rate`) so the spawn count is a `dt`-scaled Poisson draw. Cycle-time macros `@ct`/`@cycletime` rewrite to `1/arg` (`create.jl:156-159`).
2. **Capacity gate.** The proposal is clamped so that concurrent instances of `t` never exceed `transCapacity` (`solvers.jl:147-154`): `q = min(capacity − liveCount, q_desired)`, and any overflow is *deferred* to a future tick. (The deferral path `add_to_spawn!`, `state.jl:251-256`, is currently broken — see Invariant 3.)
3. **Resource allocation (genesis).** The clamped demand competes for supply via the allocator of **ADR 0002** (weighted progressive filling). Only *upfront* LHS tokens are debited at genesis — modalities `:rate` and `:nonblock` are excluded from the spawn requirement (`solvers.jl:20`, `get_reqs_init!`). Allocations are floored to whole instances (`get_init_satisfied`, `solvers.jl:110-128`), `u` is debited (`solvers.jl:170`), and an instance object is created with `t = clock`, `q = granted count`, `state = 0.0` (`solvers.jl:178-189`). For structured species the granted integer count of token agents is bound to the instance by descending token `priority` (`solvers.jl:194-225`). The recipe's `transPreAction` runs (`solvers.jl:227`).
4. **Cycle-time accumulation.** Each subsequent tick, in-flight instances compete again for *ongoing* resources (`:rate` tokens scaled by `dt` when `cycleTime > 0`, plus `:nonblock` tokens unscaled — `get_reqs_ongoing!`, `solvers.jl:31-48`). The granted fill fraction `q_frac` advances progress: `instance.state += q_frac · dt` (`solvers.jl:261`). Under contention an instance advances *slower than wall-clock*; with full allocation it advances by exactly `dt` per tick.
5. **Terminal test.** An instance terminates when **either** its cycle completes (`state ≥ cycleTime`) **or** its lifetime is exhausted (`clock − birth ≥ maxLifeTime`) (`solvers.jl:408-410`).
6. **Success draw (Binomial PoS).** On termination, the number of *successful* firings is `q_success = rand(Binomial(q, probOfSuccess))` if the cycle completed, else `0` (a lifetime-only timeout yields no successes) (`solvers.jl:412-416`).
7. **RHS emission.** For each successful firing, the RHS products are emitted into `u` at their stoichiometry; structured products spawn new token agents; the running reward ledger accumulates `specReward` (`solvers.jl:418-435`).
8. **Resource return.** LHS tokens are returned according to modality (`solvers.jl:437-483`): `:conserved` tokens return `q · stoich · (cycleTime if :rate else 1)` (the resource was *held*, not consumed); `:nonblock` tokens return `q · stoich` (held but never blocking). Plain consumed tokens (no return modality) are *not* returned — they were consumed at genesis. `:conserved` together with `:nonblock` is rejected (`solvers.jl:461-465`). Structured tokens are unbound and, for fully-consumed bindings, marked `:removed` (`solvers.jl:443-457, 487-490`).
9. **Termination.** `transPostAction` runs (`solvers.jl:485`) and the instance is pruned from the in-flight set (`solvers.jl:501` — but see Invariant 6).

### 3.3 The ordered step (one tick)

The engine implements one tick as `_step!` (`solvers.jl:642-673`), driven by the host stepping interface. The following order is **normative** — resource accounting, conservation, and determinism all depend on it:

1. **Sync structured counts** — reflect live structured-token agent counts into `u` (`solvers.jl:643`, `update_u_structured!`).
2. **Initial save** — if no history exists yet, record the initial `(t, u)` row (`solvers.jl:644-646`).
3. **Free blocked species** — release `:nonblock` resources held by in-flight instances back into `u` before this tick's allocation (`solvers.jl:648`, `free_blocked_species!`). *(Currently broken — Invariant 1.)*
4. **Update observables** — resample any observable whose `(t − last) ≥ every` (`solvers.jl:650` → `state.jl:144`).
5. **Sample transitions** — clear the per-tick transition table, evaluate each activated recipe's sampleable attributes fresh, and unfold the LHS into reactant records (`solvers.jl:651` → `state.jl:174-217`). This is where per-tick values (rate, priority, stoich, cycleTime) are realized.
6. **Evolve** (`solvers.jl:652` → `solvers.jl:133-313`): **spawn** new instances (genesis + capacity gate + genesis allocation, lifecycle stages 1-3), then **advance** all in-flight instances (ongoing allocation + progress accumulation, stage 4). Both phases allocate via ADR 0002. Per ADR 0002, priority must be re-read fresh per tick in *both* phases; the current ongoing phase reads a spawn-time snapshot (`solvers.jl:240`) and must be fixed to re-read `transPriority` per tick.
7. **Sync structured counts** (`solvers.jl:653`).
8. **Finish** (`solvers.jl:654` → `solvers.jl:400-508`): for every instance past its terminal test, run the success draw, RHS emission, resource return, post-action, and pruning (lifecycle stages 5-9).
9. **Sync structured counts** (`solvers.jl:655`).
10. **Events** — fire scheduled events once per tick (`solvers.jl:657`, `event_action!`). *(Currently a no-op — Invariant 7.)*
11. **Ledger** — push the `:valuation` row: `u' · specValuation` (`solvers.jl:659-666`); cost/reward rows were pushed inside `evolve!`/`finish!`.
12. **Advance clock** — `t += dt` (`solvers.jl:668`). This is the *single* authoritative clock advance; all tick work above occurs at the *old* `t`.
13. **Save** — record the post-increment `(t, u)` row (`solvers.jl:670`).

Termination of the run: the engine reports completion once `t > tspan[2]` (`solvers.jl:675`).

### 3.4 Invariants (the contract)

These must hold at every tick boundary (i.e. after step 13). Each is stated as a contract obligation, followed by where the current code honors or violates it.

**1. Resource non-negativity.** `u[s] ≥ 0` for every species `s`, at all times. No allocation may debit a species below zero, and the order of operations (free blocked → allocate → return) must never transiently require negative supply. *The ADR 0002 allocator guarantees this by construction.* **Violation:** `free_blocked_species!` (`solvers.jl:510-522`) references an undefined variable `q` at **`solvers.jl:512`**, throwing `UndefVarError` on any in-flight `:nonblock` LHS token. The `:nonblock`-release path (step 3) is therefore effectively dead: blocked resources that should re-enter the pool each tick may not, and the function errors whenever a `:nonblock` token is in flight. The presence of a `max(0, u[i])` clamp in the legacy allocator (`solvers.jl:68`) is itself evidence that negativity has occurred in practice.

**2. Conservation (conserved tokens returned exactly).** A token tagged `:conserved` is *held* for the instance's lifetime and returned in full on termination — never consumed. The quantity returned must equal the quantity held: `q · stoich · (cycleTime if :rate else 1)`. The closed system's conserved mass is invariant across spawn→return. *Honored at* `solvers.jl:437-458`. Caveats: (a) `:conserved + :nonblock` is correctly rejected as ill-defined (`solvers.jl:461-465`); (b) conservation is only exact if Invariant 6 holds — a re-emitting un-pruned instance would return conserved tokens repeatedly, inflating the pool.

**3. Capacity.** The number of concurrent in-flight instances of transition `t` never exceeds `transCapacity`. Overflow is deferred to later ticks, not dropped. *Gate present at* `solvers.jl:147-154`. **Violation:** the deferral helper `add_to_spawn!` (`state.jl:251-256`) is doubly broken — `findfirst` is handed a scalar `length(...)` instead of a range, and on match it increments `:transHash` (`+= n`) instead of `:transToSpawn`. Capacity-overflow deferral is non-functional: overflow is silently lost rather than carried forward, so under sustained over-demand the realized spawn rate is below contract.

**4. Integrality of instance counts.** Spawned instance multiplicities `q` are non-negative integers; structured-token stoichiometry must be integer-valued. *Honored:* spawn counts are floored to whole instances (`get_init_satisfied`, `solvers.jl:110-128`); structured stoich is checked and errors on non-integers (`solvers.jl:196-200, 268-272`); the Binomial success draw consumes `Int(trans.q)` (`solvers.jl:413`). Caveat: that `Int(...)` cast throws `InexactError` if `q` is ever non-integral, so integrality is *assumed*, not defensively coerced.

**5. Determinism under seed.** Two runs with the same model, inputs, and RNG seed produce identical trajectories. The allocator is deterministic and RNG-free (ADR 0002); all stochasticity is confined to the Poisson spawn draw, the Binomial PoS draw, and observable resampling, all of which must route through a single seeded `AbstractRNG`. **Violation:** every stochastic draw currently uses the global RNG with no seed threaded through the state — Poisson spawn (`solvers.jl:141` via `create.jl:151`), event Poisson (`solvers.jl:321`), Binomial PoS (`solvers.jl:413`), and observable sampling (`state.jl:123`). Determinism-under-seed is *not* currently achievable. ADR 0001's reproducibility obligation requires threading a state-owned RNG through all three sites. (The full seeding obligations are §4.)

**6. Termination completeness (lifetime prune).** Every instance that passes the terminal test is processed exactly once and then removed from the in-flight set — no instance emits its RHS or returns resources more than once. **Violation:** the prune at **`solvers.jl:501`** keeps instances satisfying `state < cycleTime`, which *retains* instances that terminated solely by `maxLifeTime` (their `state` never reached `cycleTime`). Such instances are re-evaluated every subsequent tick — re-running the success draw, re-emitting RHS products, and re-returning conserved/nonblock resources each tick. This violates Invariants 2 (conservation) and 4 indirectly, and is a correctness defect, not merely a leak. The prune predicate must be "remove every instance that passed the terminal test," matching `solvers.jl:408-410`.

**7. Event firing.** A scheduled event whose trigger holds fires its action `q` times per tick (`q` = 1 for a Bool trigger, `rand(Poisson(v))` for a numeric trigger). **Violation:** `event_action!` (`solvers.jl:316-326`) computes `q` correctly but at **`solvers.jl:323`** merely *fetches* `state[i, :eventAction]` inside the loop without evaluating it — the action expression is never run. Events are currently a complete no-op; any contract behavior depending on events (scheduled budget injections, what-if interventions) does not execute.

### 3.5 Notes for the rework

- Invariants 1, 5, 6, 7 are presently **unmet** by the engine and are blocking for any contract-conformance test suite. Invariants 2 and 4 hold only conditionally (2 depends on 6; 4 assumes upstream integrality).
- The allocation rule referenced throughout step 6 (spawn + ongoing) is fully specified in **ADR 0002** (priority-weighted progressive filling): work-conserving, conjunctive-consistent, deterministic, dependency-free, with `priority = 0` meaning "leftover-only." The invariant test plan in ADR 0002 is the natural home for the allocation-side checks (non-negativity, capacity, conjunctive consistency, priority-split ratio, work-conservation, determinism, spawn integrality).

---

## 4. Determinism & Seeding Contract

### 4.1 The determinism contract

The model MUST satisfy: for a fixed model specification and a fixed seed, every run produces a bit-identical trajectory.

> **D1 (Reproducibility).** Given a model spec `M` and a seed `s`, the pair `(M, s)` fully determines the trajectory: the time series of `state.u`, every `state.log` entry, every spawned/terminated transition, and every observable sample. Re-running `(M, s)` on the same platform MUST yield identical output.

> **D2 (RNG isolation).** A run MUST NOT read from or write to the global RNG (`Random.default_rng()` / `Base.GLOBAL_RNG`). All randomness MUST come from an `AbstractRNG` owned by the run state. Consequently, concurrent or interleaved runs cannot perturb one another, and code outside the model (tests, REPL, other packages) cannot perturb a run.

> **D3 (Allocation determinism).** The resource allocator is a pure deterministic function of `(requirements, available, priorities)` (ADR 0002, weighted progressive filling). It introduces no randomness and MUST remain RNG-free. Determinism of a run therefore reduces entirely to determinism of the stochastic draws enumerated below plus deterministic iteration order.

> **D4 (Order determinism).** Wherever a stochastic draw is made per object (per transition, per species, per event, per observable), the iteration order over those objects MUST be deterministic and stable across runs (e.g. iterate `parts(state, :T)` / `parts(state, :S)` in index order, never `Dict`/`Set` insertion order). The number of draws consumed and the order in which they are consumed are part of the contract.

### 4.2 Where the RNG must be threaded

Today every stochastic draw uses the implicit global RNG. The contract requires a single `rng::AbstractRNG` field on the run state, created at construction, and passed explicitly to every draw. The complete current inventory of global-RNG sites that MUST be converted:

| Draw | Current site | Distribution | Purpose |
| --- | --- | --- | --- |
| Transition spawn count | `expand_rate`, `src/interface/create.jl:151` (`rand(Poisson(max(state.dt * rate, 0)))`); realized in `evolve!`, `src/solvers.jl:140-144` | `Poisson(dt·rate)` | how many new instances of a transition to schedule this tick |
| Probability of success | `finish!`, `src/solvers.jl:413` (`rand(Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))`) | `Binomial(q, PoS)` | how many of `q` completing instances succeed |
| Observable sampling | `sample_range`, `src/state.jl:123` (`rand() * sum(...)`) and `src/state.jl:132` (`rand(r)`) | `Uniform` selector + inner `Sampleable` | choosing and sampling an observable's range entry |
| Event firing count | `event_action!`, `src/solvers.jl:321` (`v isa Number ? rand(Poisson(v)) : 0`) | `Poisson(rate)` | how many times an event's action fires this tick |
| Generic sampleable eval | `context_eval`, `src/state.jl:70` (`o isa Sampleable ? rand(o) : o`) | any `Distributions.Sampleable` | evaluating any attribute that resolves to a distribution |

The last row is the most important: `context_eval` (`src/state.jl:67-71`) is the single chokepoint through which nearly every attribute value (rates, stoichiometries, costs, etc.) flows. Threading the RNG through `context_eval(state, transition, o)` — i.e. `rand(state.rng, o)` — covers the majority of draws by construction. The four explicit `rand(...)` call sites above MUST additionally take `state.rng` as their first argument.

> **D5 (Threading rule).** Add `rng::AbstractRNG` to `ReactionNetworkProblem` (`src/state.jl:40-63`). Every `rand` / `rand(dist)` reachable from a step MUST become `rand(state.rng, ...)`. No code path inside `_step!` (`src/solvers.jl:642-673`) may call `rand` without an explicit RNG argument. A grep for `rand(` lacking an RNG first argument inside `src/` is the enforced invariant.

### 4.3 Construction and re-init

> **D6 (Seed at construction).** `ReactionNetworkProblem(...)` MUST accept a `seed::Union{Integer,Nothing}=nothing` (and/or `rng::AbstractRNG`) keyword. When a seed is given, the state's RNG is initialized deterministically from it (e.g. `Xoshiro(seed)`). When neither is given, a fresh RNG seeded from system entropy is created and its seed recorded in the log so the run can be replayed.

> **D7 (Re-init restores the stream).** `_reinit!` (`src/solvers.jl:619-628`) currently resets `u`, `t`, `ongoing_transitions`, `log`, `observables`, and `sol`, but does NOT reset the RNG — so a second run after `reinit!` diverges. The contract requires `_reinit!` to restore the RNG to the exact state implied by the original `(M, seed)`, so that `init → step* → reinit! → step*` reproduces the first trajectory. The seed (or the initial RNG state) MUST be stored on the state for this purpose.

### 4.4 Ensembles

> **D8 (Per-trajectory seeding from index).** An ensemble of `N` trajectories MUST derive each member's seed deterministically from a single root seed and the member index `k ∈ 1:N`, e.g. `member_seed = hash((root_seed, k))` or `Xoshiro(root_seed)` jumped/split per member. This guarantees: (a) the whole ensemble is reproducible from one root seed; (b) members are mutually independent streams; (c) member `k`'s trajectory does not depend on `N` or on the order members are run (so ensembles may be parallelized or resumed). The per-member seed MUST be recorded in that member's log.

> **D9 (No shared mutable RNG across members).** Ensemble members MUST NOT share one `AbstractRNG` instance; each member owns its own RNG. This is what makes parallel execution (D8c) safe and deterministic.

### 4.5 Out of scope / known non-determinism to fix separately

- Float reduction order (e.g. `sum(allocs; dims=2)` at `src/solvers.jl:170,255`) is deterministic for a fixed array shape and is therefore covered by D1 on a fixed platform; cross-platform bit-identity is explicitly NOT promised by this contract.
- `eval(...)` of authored expressions during compilation (`compileval`, `src/state.jl:119`; `sample_range`, `src/state.jl:123`) is a build-time concern, not a per-step draw, but any randomness it triggers at runtime is still bound by D5.

---

## 5. Attribute Contract (typed attributes)

### 5.1 Motivation and the replacement

Today every authored quantity has the catch-all type `SampleableValues = Union{Expr,Symbol,AbstractString,Float64,Int,Function}` (`src/ReactiveDynamics.jl:10`), used for all of `transRate`, `transPriority`, ..., `specValuation` (`src/ReactiveDynamics.jl:46-62`). This type encodes nothing about units, sign, range, or whether time-variation is permitted, so the only validation is whatever the solver happens to do at runtime. The attribute contract replaces this with a per-attribute specification: each attribute has a name, a **kind**, a **scalar domain** (units + valid range + sign), a **default**, and a **time-variation policy**.

### 5.2 Definitions used below

- **Time-varying expression (TVE):** an attribute value that, instead of a literal number, is an `Expr`/`Function` evaluated each time it is read via `context_eval` (`src/state.jl:67-71`) in the context `(state, transition)`, and which may resolve to a `Distributions.Sampleable` (then sampled per the determinism contract, §4). "May be a TVE = yes" means the solver re-reads it dynamically; "no" means it is read once at construction and frozen.
- **Default source:** literal default applied by `assign_defaults!` via `defargs` (`src/ReactiveDynamics.jl:117-139`), unless noted as "required at authoring."
- **Validation:** the contract REQUIRES validation at construction (`ReactionNetworkProblem`, `src/solvers.jl:536`) for static literals, and at read time for any value resolved from a TVE. Violations MUST error, not silently clamp — except where a clamp is already the documented solver behavior, which is called out explicitly.

### 5.3 Transition attributes

| Attribute | Units | Valid range | Default | TVE? | Notes / current behavior |
| --- | --- | --- | --- | --- | --- |
| `transRate` | instances · time⁻¹ (intensity) | ≥ 0 | **required at authoring** — set from the reaction line (`merge_acs!`, `src/ReactiveDynamics.jl:202`: `transRate = t[1][1]`); has no entry in `defargs`. | **yes** | Spawn count per tick is `Poisson(max(dt·rate, 0))` (`src/interface/create.jl:151`). The `max(...,0)` already clamps negatives; the contract makes "< 0 is invalid" explicit for literals. `@deterministic rate` bypasses the Poisson and uses the value directly (`create.jl:150-154`). |
| `transPriority` | dimensionless fill-weight | ≥ 0 | `1` (`:transPriority => 1`, `src/ReactiveDynamics.jl:119`) | **yes (dynamic, per ADR 0002)** | Per-tick weight in weighted progressive filling. `priority = 0` ⇒ leftover-only (served only after all positive-priority demand). Read dynamically in `evolve!` via `state[:, :transPriority]` (`src/solvers.jl:158`) and `t[:transPriority]` for ongoing (`src/solvers.jl:240`). |
| `transCycleTime` | time | ≥ 0 | `0.0` (`src/ReactiveDynamics.jl:122`) | yes | Duration before an instance may complete. A transition completes when `trans_.state ≥ transCycleTime` (`src/solvers.jl:409,412,501`). With `transCycleTime = 0` the instance can complete on the tick it starts. Also scales `:rate`-modality conserved returns (`src/solvers.jl:442`) and ongoing `:rate` requirements by `dt` (`src/solvers.jl:36-37`). |
| `transProbOfSuccess` | probability | [0, 1] | `1` (`src/ReactiveDynamics.jl:120`) | yes | On completion, successes drawn as `Binomial(q, PoS)` (`src/solvers.jl:413`). `Int(trans_.q)` is required there, so PoS combines with an integer instance count. Out-of-[0,1] MUST error (currently would throw deep inside `Binomial`). |
| `transCapacity` | instances | ≥ 0, or `Inf` | `Inf` (`src/ReactiveDynamics.jl:121`) | yes | Max concurrent live instances of this transition; excess is deferred to `transToSpawn` (`src/solvers.jl:148-153`). `Inf` = uncapped. Must be ≥ 0; non-integer finite capacities are compared against integer instance counts so SHOULD be integer-valued. |
| `transMaxLifeTime` | time | ≥ 0, or `Inf` | `Inf` (`src/ReactiveDynamics.jl:123`) | yes | An instance is force-collected once `(state.t - trans_.t) ≥ transMaxLifeTime` even if it has not reached `transCycleTime` (`src/solvers.jl:408`); on such timeout `q = 0` successes (no products), since success requires `state ≥ transCycleTime` (`src/solvers.jl:412`). `Inf` = never times out. |
| `transMultiplier` | dimensionless | ≥ 0 | `1` (`src/ReactiveDynamics.jl:124`) | yes | Scales the spawn intensity: `qs[i] = transRate · transMultiplier` before `ceil` (`src/solvers.jl:140-144`). `0` disables spawning. |
| `transName` | — (label) | any | `missing` (`src/ReactiveDynamics.jl:127`) | no | Descriptive only; not a quantity. Listed for completeness. |
| `transPreAction` / `transPostAction` | — (code) | callable / `:()` | `:()` (`src/ReactiveDynamics.jl:125-126`) | n/a (action) | Side-effecting code run on spawn / completion (`src/solvers.jl:227,485`); read from `acs` directly, never through `context_eval`'s sampling path (`src/state.jl:73-83`). Not numeric attributes. |

### 5.4 Species attributes

| Attribute | Units | Valid range | Default | TVE? | Notes / current behavior |
| --- | --- | --- | --- | --- | --- |
| `specInitVal` | tokens (count of the species) | ≥ 0 | `0.0` (`src/ReactiveDynamics.jl:131`) | no (read once) | Initial `u`. Used at construction (`src/solvers.jl:561-569`) and re-init (`src/solvers.jl:620`); `u0` keyword overrides per species. Stored as `Float64` but represents a token count, so SHOULD be a nonnegative integer for unstructured species. |
| `specInitUncertainty` | same units as `specInitVal` (absolute) or dimensionless (relative) — MUST be fixed by the contract; recommend **relative, dimensionless ≥ 0** | ≥ 0 | `0.0` (`src/ReactiveDynamics.jl:130`) | no (read once) | Spread applied to `specInitVal` when sampling initial conditions for an ensemble. **Currently declared and defaulted but not consumed by the solver** (no read site in `solvers.jl`/`state.jl`); the contract pins its semantics (interpretation + the distribution used) so the determinism contract §4 governs the draw. |
| `specCost` | value · token⁻¹ (currency per unit consumed) | ≥ 0 (typical); finite | `0.0` (`src/ReactiveDynamics.jl:132`) | yes | Cost charged per allocated unit. Logged each tick as `actual_allocs' · specCost` (`src/solvers.jl:308-311`). For BD/rNPV ledgers this is the spend side. |
| `specReward` | value · token⁻¹ (currency per unit produced) | ≥ 0 (typical); finite | `0.0` (`src/ReactiveDynamics.jl:133`) | yes | Reward credited per produced unit on transition completion (`src/solvers.jl:426,433`). The income side of the ledger. |
| `specValuation` | value · token⁻¹ (currency per unit held) | finite (may be signed) | `0.0` (`src/ReactiveDynamics.jl:134`) | yes | Mark-to-market value of the current holding; logged each tick as `u' · specValuation` (`src/solvers.jl:663-664`). May legitimately be negative (e.g. a liability), so the contract does NOT impose ≥ 0 here. |
| `specModality` | — (set of tags) | ⊆ {`:nonblock`, `:conserved`, `:rate`} | `Set{Symbol}()` (`src/ReactiveDynamics.jl:154`) | no | The species-level modality set, unioned with per-reactant modality during sampling (`src/state.jl:204`). The allowed universe is `species_modalities` (`src/ReactiveDynamics.jl:144`); the contract requires validation against this set (today it is an unvalidated `Set{Symbol}`). `:conserved` + `:nonblock` together is illegal and MUST be rejected at authoring (today only caught at completion, `src/solvers.jl:461-465`). See §1 for the orthogonalized re-model of these tags. |
| `specStructured` | — (flag) | `Bool` | `false` (`src/ReactiveDynamics.jl:135`) | no | Marks a species as a structured token (agent-backed) rather than a plain count. Determines `structured_token_names` (`src/solvers.jl:556-557`). Listed for completeness; not numeric. |

### 5.5 Cross-cutting rules

> **A1 (Sign/range enforcement).** Every range in §5.3–5.4 MUST be enforced: at construction for literals; at each read for TVEs. The only sanctioned silent clamps are the existing `max(dt·rate, 0)` for spawn intensity (`src/interface/create.jl:151`) and the `max(0, u[i])` guard in `alloc_weighted!` (`src/solvers.jl:68`); all other violations MUST error.

> **A2 (Units are advisory but fixed.)** The engine is unit-agnostic at runtime (all quantities are `Float64`), but the contract FIXES the unit interpretation per attribute so that authored models, ledgers, and observables compose consistently. `transRate`'s time unit MUST match the model's `tunit` (`get_tcontrol`, `src/solvers.jl:526-534`); `specCost`/`specReward`/`specValuation` MUST share one currency unit for the valuation logs (`src/solvers.jl:308,663`) to be meaningful.

> **A3 (TVE policy).** Attributes marked "TVE? = yes" MAY be authored as an `Expr`/`Function`/`Sampleable` and are re-read each tick through `context_eval`; attributes marked "no" are read once (construction or re-init) and frozen for the run. Authoring a TVE for a "no" attribute MUST be rejected.

> **A4 (Integrality).** `specInitVal` (unstructured), `transCapacity` (finite), and any stoichiometry feeding a structured species are compared against or converted to integers at runtime (`Int(trans_.q)` at `src/solvers.jl:413`; the `isinteger` checks at `src/solvers.jl:196-200,268-272`). The contract requires these to be integer-valued; non-integer literals MUST be rejected at construction.

> **A5 (Replacement of the catch-all type).** `SampleableValues` (`src/ReactiveDynamics.jl:10`) is replaced as the *semantic* contract by the per-attribute kinds above. The underlying ACSet `AttrType` MAY remain a broad union for storage, but the validation layer defined here is authoritative: a model that type-checks against the ACSet but violates §5.3–5.5 is invalid.

---

## Pending sections (next contract increment)

The ADR 0003 data-store decision is now made, so these sections can be written against the dependency-free typed IR (with the promoted `ReactantSpec` incidence table) and the ADR 0005 JSON serialization:

- **Object model** — the canonical typed representation of species, transitions, the promoted `ReactantSpec` incidence relation (ADR 0003 Phase 2), events, observables, params, and meta; the append-only index invariant (ADR 0004) is part of this section's contract.
- **Composition semantics** — how two models compose. Today this is manual name-matching `add_part!`/`incident` loops (`operators/joins.jl` `union_acs!`, `operators/equalize.jl` `equalize!`), NOT categorical pushouts; the `@join` file branch calls an undefined `include_model` (`joins.jl:226,228`). The promoted reactant table makes species-merge structurally exact (vs today's fragile `recursively_substitute_vars!` string surgery). NOTE: `equalize!`'s `rem_parts!` (`operators/equalize.jl:52`) is the one reindexer that ADR 0004 forbids on a live/stepping model.
- **Serialization schema** — specified in [ADR 0005](adr/0005-serialization-json-ir.md): a single JSON format + typed `ExprNode` IR (eval-free, JSON-Schema-describable for LLM emission, Workstream F), with `from_json`/`to_json`/`validate` and round-trip guarantees. This section will cross-reference rather than restate ADR 0005.
