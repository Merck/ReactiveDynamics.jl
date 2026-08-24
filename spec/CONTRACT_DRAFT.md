# ReactiveDynamics Modeling Contract

> **DRAFT — Phase-0 modeling contract.** The fork-independent sections (§1–§5) are complete. The ADR 0003 data-store decision is now ACCEPTED (drop ACSets for a dependency-free typed IR; promote the arc relation), with serialization in [ADR 0005](adr/0005-serialization-json-ir.md) (single JSON + ExprNode) and runtime mutation in [ADR 0004](adr/0004-runtime-mutation.md) (append-only). The remaining sections — object model, composition semantics, serialization schema — can now be written against that decision and are the next contract increment (see "Pending sections" at the end).

## Orientation

ReactiveDynamics (RD) is a timed, stochastic, resource-constrained Petri net — a discrete-event system in which transitions periodically spawn in-flight instances that draw on shared, finite resource pools — **places**, in Petri-net terms — over a service duration, then emit products and return or consume those resources on completion. It is *not* a chemical reaction network and *not* a continuous-time Gillespie process: it is a fixed-step discrete-time engine with priority-weighted resource allocation. This contract pins the *operational semantics* of that engine — the time model, the per-tick firing and lifecycle rules, the modality (resource-claim) semantics, the determinism/seeding obligations, and the typed attribute domains — independently of how the static model is stored. Concrete `file:line` references point at the current `ref-agents` engine and are evidence of present behavior (sometimes evidence of a *bug*, flagged as such), not a commitment to keep that behavior. The object model, composition semantics, and serialization schema are deferred to a later revision pending ADR 0003.

The contract is organized as follows: §1 Modality Truth Table (the resource-claim semantics that the allocator and lifecycle depend on); §2 Time Model (single discrete clock, spawn intensity, cycle/lifetime, dt-invariance); §3 Operational Semantics (instance lifecycle, the ordered step, and the invariants an implementation must uphold); §4 Determinism & Seeding Contract; §5 Attribute Contract (typed attribute domains). ADR 0002 (priority-weighted progressive filling) is the normative specification for the allocator referenced throughout. The **Glossary** immediately below fixes the vocabulary the rest of the document uses (ADR 0017).

---

## Glossary — the Petri-net term dictionary (see [ADR 0017](adr/0017-petri-net-vocabulary.md))

This contract uses standard Petri-net vocabulary throughout: a resource pool is a **place**, the quantity in it is that place's **marking**, and a transition↔place participation is an **arc** whose multiplicity is its **arc weight**. [ADR 0017](adr/0017-petri-net-vocabulary.md) retired the inherited chemical-reaction-network words (*species*, *reactant*) that earlier ADRs and the pre-0.3 API used; ADRs are append-only, so ADRs 0001–0016 keep their original wording and should be read with this table in hand.

**A token is a discrete unit of resource sitting in a place** — the Petri-net sense, unrelated to language-model tokens. *Token* is canonical and is NOT renamed. Neither is *transition*.

The table below is the durable artifact of ADR 0017: it maps each RD concept onto its published name and the Petri-net extension family that name comes from, so a reader knows which literature answers a question about the engine.

| RD concept | standard term | family / note |
|---|---|---|
| fungible resource pool | **place**; its count is that place's **marking** | classical P/T net. `M₀` is the initial marking (§10.3) |
| structured-token pool | place with a **colour set**; its tokens are distinguishable | Coloured Petri nets (Jensen). Token attributes are colours (§9) |
| LHS / RHS entry | input **arc** / output arc; multiplicity is the **arc weight** (inscription) | classical (§6.5) |
| LHS / RHS multiset | **preset** `•t` / **postset** `t•` | classical |
| transition | **transition** | already canonical — not renamed (§6.6) |
| `rate` expression | **firing rate**, marking-dependent | Stochastic PN / GSPN. RD's Poisson intensity with unbounded concurrency is **infinite-server** semantics (§2.3) |
| `capacity` on a transition | **k-server** semantics | GSPN. Distinct from **place capacity**, an unrelated classical notion RD does not implement |
| `cycletime` | **firing duration** | Timed PN (Ramchandani). NOT Time PN, whose parameter is a firing *interval* (§2.4) |
| `probability` on completion | random switch / probabilistic firing outcome | GSPN (§3.3 step 7) |
| `priority` + the allocator | **priority** and random switches, extended by RD's rationing rule | GSPN + [ADR 0002](adr/0002-priority-weighted-allocation.md) |
| `:nonblock` modality | **read arc** / test arc | classical extension (§1.4) |
| `:conserved` modality | **self-loop** / side condition | classical (§1.4) |
| `:rate` modality | **continuous transition** | hybrid / continuous PN (§1.4) |
| conservation invariant (`S+I+R`) | **P-invariant** (S-invariant) | classical structural analysis |
| resources held over a duration, contended | closest published relative: **Queueing Petri nets** (QPN) | Bause |
| structured token selected by predicate | ≈ **binding element** (a transition plus a variable binding) | Coloured PN (§9.5) |
| the tokens a live firing holds | its **binding** — `bound_tokens`, `nonblock_tokens`, `binding` on the firing. NOT an arc: an arc is static topology with no runtime instance | Coloured PN ([ADR 0017](adr/0017-petri-net-vocabulary.md) A1) |
| the arrow-form lines a model is authored in | **reaction notation** for a Petri net. Under the standard correspondence (species↔places, reactions↔transitions) a reaction network *is* a Petri net, so RD keeps the chemistry register for the notation (`@reaction_network`, "reaction line") and Petri vocabulary for the object model | classical ([ADR 0017](adr/0017-petri-net-vocabulary.md) A3) |

Retired names and their replacements (each old spelling resolves for one release, with a deprecation warning, per [ADR 0015](adr/0015-post-acsets-naming.md)'s shim pattern): `@add_species` → `@add_place`; `SetSpecies` → `SetMarking`; `ReactantSpec` → `ArcSpec`; `reactant_specs` → `arcs`; `specname` → `placename`; `register_structured_species!` → `register_token_kind!`; `get_species`/`set_species!` → `get_place`/`set_place!`; and the `:S` column symbols `specName`/`specModality`/`specCost`/`specValuation`/`specInitVal`/`specReward`/`specStructured`/`specRole`/`specInitUncertainty` → `placeName`/`placeDefaultModality`/`placeCost`/`placeValuation`/`placeInitVal`/`placeReward`/`placeStructured`/`placeRole`/`placeInitUncertainty`. The `SCHEMA` object symbol for places stays `:S` (`:P` is the parameter object). An arc's coefficient is `multiplicity`, retiring `stoich`; it is not spelled `weight`, because `weight` already names `@choose` alternative weights and the allocator's fill-rate weights ([ADR 0017](adr/0017-petri-net-vocabulary.md) A5). The runtime fields carrying a firing's tokens are `bound_tokens`/`nonblock_tokens`/`binding`, retiring `bound_structured_agents`/`nonblock_structured_agents`/`structured_to_agents`, and the per-tick evaluated arc record is `ResolvedArc`, retiring `UnfoldedArc` (A1, A5).

Retired **serialized keys** and their replacements, on the same one-release shim (the writer emits only the new spelling; the loader accepts the old one and warns — §8): the top-level arrays `species[]` → `places[]` and `reactants[]` → `arcs[]`; an arc's `species` → `place`; a `population[]` entry's `species` → `place`; the action verb `set_species` → `set_marking`; a `ref` node's `kind` value `species` → `place` (`REF_KINDS`); the export-bundle manifest key `species` → `places` and a per-token trajectory record's `species` → `place` (§14.3); the result-frame column `:species` → `:place` (§14.1, a `DataFrame` column user analysis code reads); and the `@select`/`@advance` field name `:species` → `:place` (§9.5 — the one alias accepted *silently*, because a deprecation warning there would fire inside the step loop); and an arc's `stoich` → `multiplicity`.

**Two registers, deliberately.** This contract and the API reference say *place*, *marking* and *arc*. The tutorials, case studies, and demo prose say **resource pool** for the same object, glossed once at first use as "a place, in Petri-net terms". Both name one thing; nothing in the engine distinguishes them. A model that wants a third register does not need one of ours: `@aka` renames the objects per model (its own example is `@aka net place = resource transition = reaction`), which is the sanctioned channel for domain drift — canonical names in the API, whatever the readers use in the model file.

**When to adopt a standard term and when to drift** ([ADR 0017](adr/0017-petri-net-vocabulary.md) A2). Adopt the published word when RD's concept is a *superset* of it (an arc carrying a modality is still an arc; a firing with a duration is Ramchandani's timed firing). Drift deliberately, and say so, when RD's concept *contradicts* a commitment the word carries. Never reuse a standard word for a concept at a different *layer* than the one it names — that is the arc-versus-binding trap, and the one failure mode a glossary cannot repair.

---

## 1. Modality Truth Table

A *modality* is the contract attached to each LHS (consumed) token of a transition. It answers three independent questions: **when** is the resource claimed against the pool, **whether** it is returned when the cycle completes, and **whether** the claim blocks the resource while the cycle runs. Today these three concerns are conflated into an unvalidated `Set{Symbol} ⊆ {:nonblock, :conserved, :rate}` (the allowed universe `place_modalities`, `ReactiveDynamics.jl:144`) whose cross-product the solver interprets inconsistently and never documents. This section re-models them as three orthogonal, closed, typed fields and gives the full legal truth table.

### 1.1 The three orthogonal axes

Each consumed token declares a `Modality` value with exactly these three fields. All combinations are legal unless explicitly marked illegal below.

| Field | Type (closed set) | Meaning | Default |
|---|---|---|---|
| `allocation` | `{upfront, perstep}` | **When** the resource is claimed. `upfront`: the full `q·multiplicity` is reserved once, at spawn time, as a precondition for the cycle to start at all. `perstep`: a slice is claimed each tick while the cycle runs (a flow/throughput draw rather than a one-time seizure). | `upfront` |
| `return` | `{consumed, conserved}` | **Whether** the claimed amount returns to the pool when the cycle finishes. `consumed`: permanently destroyed (true consumption). `conserved`: the reserved amount is credited back at finish (a loan / temporary hold). | `consumed` |
| `blocking` | `{block, nonblock}` | **Whether** the claim holds the resource for the duration. `block`: the resource is unavailable to other transitions until finish. `nonblock`: the resource is touched but immediately released every step, so it never actually constrains anyone (a soft / advisory draw). | `block` |

These map onto the legacy tags as: empty set ⇒ `(upfront, consumed, block)`; `:rate` ⇒ `allocation = perstep`; `:conserved` ⇒ `return = conserved`; `:nonblock` ⇒ `blocking = nonblock`. Note the legacy tags are *not* one-axis-each in isolation — `:nonblock` in the current code simultaneously implies per-step reservation *and* non-blocking release *and* (forced) non-conservation, which is exactly why it cannot be combined with `:conserved`. The re-model separates these so the constraint becomes a single clean rule (below) rather than an ad-hoc `error`.

### 1.2 Effect dimensions each combination resolves to

Every legal `Modality` resolves deterministically to three solver effects:

- **Allocation timing** — *spawn-reserve* (counted by `get_reqs_init!`, must be satisfiable for the transition to spawn) vs *per-step-reserve* (counted by `get_reqs_ongoing!` each tick, optionally `dt`-scaled).
- **Consumption / return semantics** — does `finish!` credit `state.u` back, and by how much.
- **Ledger effect** — what hits the period cost line (every reserved unit is valued at `placeCost` in the `:valuation_cost` log entry) and whether it is later offset by a reward/return.

### 1.3 Full truth table (every legal combination)

`s = multiplicity`, `q = transition multiplicity`, `Δt = state.dt`, `C = transCycleTime`. "Reserve" = subtracted from `state.u` at the indicated time and charged to period cost at `placeCost`.

| # | `allocation` | `return` | `blocking` | Legacy `Set{Symbol}` | Allocation timing | Return at finish | Ledger effect | Plain meaning |
|---|---|---|---|---|---|---|---|---|
| 1 | `upfront` | `consumed` | `block` | `{}` (empty) | Reserve `q·s` once at spawn (must fit) | nothing | Net spend of `q·s` valued at `placeCost`; gone for good | **Raw consumption.** Classic consumed input — burned to start the cycle. |
| 2 | `upfront` | `conserved` | `block` | `{:conserved}` | Reserve `q·s` once at spawn (must fit) | `+ q·s` | Cost charged at spawn, fully credited back at finish; net-zero on success | **Held capital / equipment.** Seized for the duration, returned intact. |
| 3 | `perstep` | `consumed` | `block` | `{:rate}` | Reserve `q·s·Δt` each tick (only if `C > 0`) | nothing | Continuous spend; period cost accrues per tick | **Metered consumption (flow).** A burn *rate* drawn down over the cycle. |
| 4 | `perstep` | `conserved` | `block` | `{:rate, :conserved}` | Reserve `q·s·Δt` each tick (only if `C > 0`) | `+ q·s·C` | Per-tick cost accrued, then credited back `q·s·C` at finish | **Rented throughput.** A rate-based hold, fully returned at completion. |
| 5 | `perstep` | `consumed` | `nonblock` | `{:nonblock}` | Reserve `q·s` each tick, **freed every step** by `free_blocked_places!` | nothing (already freed) | Touched/measured but not held; no net pool change between ticks | **Soft / advisory draw.** Reads the resource each step without contending for it. |

Rows 1–4 are the four corners of the `allocation × return` plane under `block`. Row 5 is the single coherent `nonblock` form. Every legacy `Set{Symbol}` value in the wild maps to exactly one of these five rows.

### 1.4 Illegal combinations

| `allocation` | `return` | `blocking` | Legacy form | Why illegal |
|---|---|---|---|---|
| any | `conserved` | `nonblock` | `{:nonblock, :conserved}` | Conservation means "credit the held amount back **at finish**"; non-blocking means "release it **every step**." A resource cannot be both held-until-finish and continuously-released. The current code rejects this explicitly (`error` at `solvers.jl:461-465`); the re-model makes it a single validation rule: **`blocking = nonblock` requires `return = consumed`.** |
| `perstep` | any | any (structured place) | `{:rate, ...}` on a structured token | `:rate` is unsupported for structured/agent places (`solvers.jl:38-42`): you cannot reserve a fractional, `dt`-scaled slice of an indivisible agent. Validation rule: **`allocation = perstep` requires a non-structured (countable) place.** |
| `perstep` (`return = consumed`) when `C = 0` | — | — | `{:rate}` with `transCycleTime = 0` | Per-step reservation only fires when `C > 0` (`get_reqs_ongoing!:36`). With `C = 0` a `perstep` token silently reserves nothing — a foot-gun, not a meaning. Validation rule: **`allocation = perstep` requires `transCycleTime > 0`.** |

All other combinations are legal. This collapses the previous 2³ = 8 implicit tag subsets (most of them undefined behavior) into **5 legal rows + 1 illegal rule**, all named and validated at construction.

### 1.5 Connection to ADR 0002 (`build_requirements!`)

ADR 0002's weighted progressive-filling allocator consumes a single requirements matrix `reqs[place, transition]`; the modality axes are precisely what decide **which tokens land in which requirements pass and with what coefficient**. The two legacy entry points `get_reqs_init!` / `get_reqs_ongoing!` unify into one parameterized builder:

```
build_requirements!(reqs, qs, model; counted_modalities, dt_scale)
```

- **`counted_modalities`** selects which `allocation` value this pass reserves for:
  - spawn pass ⇒ `counted_modalities = {allocation = upfront}` (rows 1, 2). This is exactly today's "exclude `:rate` and `:nonblock`" filter (`get_reqs_init!:20`, `get_init_satisfied:115`), now stated positively.
  - ongoing pass ⇒ `counted_modalities = {allocation = perstep}` (rows 3, 4, 5) — i.e. today's `:rate`/`:nonblock` branch (`get_reqs_ongoing!:35,43`).
- **`dt_scale`** is the per-tick coefficient applied to the `perstep` reservation:
  - `dt_scale = Δt` for `return = consumed` *or* `conserved` flow tokens (rows 3, 4 — today's `:rate` path, `solvers.jl:37`), gated on `C > 0`.
  - `dt_scale = 1` for `nonblock` tokens (row 5 — today's unscaled `:nonblock` path, `solvers.jl:43`).

The `return` axis is orthogonal to allocation and is consumed only by `finish!`: `return = conserved` credits `state.u` back by `q·s` (row 2) or `q·s·C` (row 4, the `perstep` flow integral) — exactly the existing `(in(:rate) ? transCycleTime : 1)` factor at `solvers.jl:438-442`, now derived from the typed fields rather than from tag intersection. The `blocking` axis is consumed only by `free_blocked_places!`: `nonblock` tokens are freed each step (`q·s`), and once the undefined-`q` bug (`solvers.jl:512`) is fixed this must credit back exactly what the ongoing pass reserved for that token (`q·s` with `dt_scale = 1`), keeping the per-step reserve/free pair conservative.

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

**Naming (`dt` vs `tstep`) — resolved 2026-07-16 (separate pass, post ADR 0015).** The runtime struct field is `dt` (`state.jl:76`), and every consumer reads `state.dt`. The constructor now stores the resolved step under the keyword key `:dt` (`keywords[:tspan], keywords[:dt] = get_tcontrol(...)`) and fills the `dt` field from `get(keywords, :dt, 1)` — so the accepted keyword, the internal keyword-bag key, and the live field all use `dt`, reconciled by NAME, not by argument position. **Contract:** authors set `dt` (or `tstops`). The legacy `tstep` meta keyword is an accepted-but-deprecated alias: if supplied (and `dt` is not), it is mapped to `dt` with a `Base.depwarn` before `get_tcontrol` resolves the step, honored for one release (ADR 0015 Tier 1 pattern). Note: this also fixed a latent silent-ignore bug — `get_tcontrol` only ever read `:dt`/`:tstops`/`:tunit`, so a pre-pass user-supplied `tstep` was overwritten (ignored) rather than applied.

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

**Genesis is already two stages (verified).** Each tick, for transition `t`: (1) a spawn-count PROPOSAL from the rate expression — `Poisson(dt·rate)` by default (`create.jl:151`) or a bare count under `@deterministic` (`create.jl:153`), times `transMultiplier`, then `ceil`'d (`solvers.jl:140-144`); then (2) an upfront-LHS RESOURCE GATE that clamps the realized count to `floor(allocation / multiplicity)` over the transition's upfront-consumed LHS tokens, via `get_reqs_init!`→`get_allocs!`→`get_init_satisfied` (`solvers.jl:156-160`, `:110-128`). The gate's `reqs == 0 ⇒ Inf` rule (`solvers.jl:121`) means a transition with NO upfront LHS keeps its full proposal (a pure source), whereas a transition WITH upfront LHS can only spawn as many instances as its input tokens and allocated resources permit (token-flow / capacity-limited). **So flow-triggered and capacity-limited genesis work today with zero new mechanism** — the `toy_pharma_model` already demonstrates both (phase hand-off via `candidate_compound` on discovery's RHS and `dx2market`'s LHS; capacity-limited starts via `3*@conserved(scientist) + @rate(budget) --> candidate_compound`).

**Source vs routing (the latent, correct distinction).** The split between exogenous arrivals and internal routing is already implicit in the data: an EMPTY LHS (`extract_arcs` returns `[]`, `reaction_parser.jl:50-51`) = a SOURCE that bypasses the gate and spawns at the proposal intensity; a NON-EMPTY upfront LHS = ROUTING that fires only when its input tokens exist. The contract makes this explicit rather than accidental.

**The closed genesis-mode tag.** Add ONE append-only transition field `genesis ∈ {poisson, scheduled, flow, capacity}` (a 4-value string tag — eval-free, JSON-Schema-enumerable per ADR 0005; append-only per ADR 0004), an intent declaration over the single shared two-stage execution path:

| Mode | Business meaning | Maps to (existing mechanism) | Source/Routing | New engine surface |
|---|---|---|---|---|
| `poisson` | Random independent exogenous inflow (unsolicited inbound deals) | Today's `expand_rate` `Poisson(dt·rate)` (`create.jl:151`) | Source | none — current default path |
| `scheduled` | Prescribed/calendar genesis (quarterly gate, budget cycle, planned start, the BD acquisition lever) | `@deterministic` bare-count (`create.jl:150,153`) + a documented calendar idiom `@scheduled(period,N) ⇒ @deterministic(N*periodic(period))`, since `periodic(state,period)` already exists (`state.jl:239`) and `period==0.0` returns true (one-shot at t0) | Source | none — documented idiom only |
| `flow` | Token-triggered: a program advances because the upstream phase succeeded | The EXISTING upfront-LHS gate (`solvers.jl:110-128`); author writes a high/`Inf` nominal rate with the upstream place as an upfront-consumed LHS input arc | Routing | none |
| `capacity` | Start as many as resources/headcount allow | Same gate path as `flow`, bounded by a renewable `@conserved`/`@rate` pool + `transCapacity` (`solvers.jl:147-153`) + the ADR-0002 allocator | Routing | none |

`batch` is a PARAMETER, not a mode (`batch::Int`, default 1, applied as post-gate rounding `qs[i] = batch * fld(qs[i], batch)`), and finite-population/one-shot is EXPRESSED (a flow source consuming a finite `placeInitVal` pool; one-shot = `scheduled` at t0), not a new mode. Three of the four modes map onto existing mechanisms with zero new engine surface; only the `scheduled` calendar form needs a documented idiom (and `periodic` already compiles, `compilers.jl:67`).

**Default (maintainer-confirmed 2026-06-20).** The ENGINE default for `genesis` is `poisson` — it is what the current engine already does (`expand_rate` → `Poisson(dt·rate)`, `create.jl:151`), so it is the backward-compatible, least-surprising default for an arbitrary transition; switching the engine default to `flow` would silently change the realized dynamics of every existing model. SEPARATELY, the documented default FOR THE BD/PHARMA-PIPELINE TEMPLATE is `flow`: pipeline/routing transitions (PhaseN → PhaseN+1) are authored explicitly as `flow`, since a program advances because the upstream phase succeeded, not on an independent Poisson clock. So: engine default = `poisson`; template/authoring guidance for routing = `flow` (explicit). This split keeps existing models stable while steering new business-process models toward the correct routing semantics.

**Validation rules (construction-time).** `flow`/`capacity` REQUIRE ≥1 upfront-consumed LHS place (else they would silently behave as an unbounded Poisson/deterministic source); `poisson`/`scheduled` with a consuming LHS should WARN (the author probably meant routing). This makes the source/routing intent checkable for agentic authoring.

**Rejected (out of scope).** The heavier queueing-lens proposal — redefine routing to drop the Poisson proposal and substitute a `min(1, dt·rate)` service-fraction primitive, splitting genesis into `ArrivalProcess` subtypes + a routing primitive — is rejected: it is a semantic break adding a second intensity code path for marginal fidelity gain, and the existing `rate=Inf + upfront LHS` idiom already realizes token-bounded firing `min(proposal, tokens) = tokens`. The one real artifact it identifies (a routing token left unserved because a small Poisson draw thinned it) is eliminated by the canonical `flow` idiom (`rate=Inf` makes the proposal non-binding), so it is a template/documentation fix, not an engine fix.

**Correctness blockers (independent of genesis design).** Two pre-existing bugs gate the business-process use cases: `event_action!` is a no-op (`solvers.jl:323` fetches `:eventAction` but never evaluates it — §3.4 Invariant 7), so the acquisition lever must use the `scheduled` rate-expression idiom, NOT the event channel, until repaired; and `add_to_spawn!` is doubly broken (`state.jl:251-256` — §3.4 Invariant 3), so `capacity`/`batch` genesis under-delivers under sustained over-demand. Also, because the `ceil` at `solvers.jl:144` breaks dt-invariance for non-Poisson counts (§2.3), `scheduled`/`batch` counts MUST be integer-valued (or `ceil` must be conditioned to the `poisson` path only).

---

## 3. Operational Semantics

This section specifies the per-tick firing and lifecycle rules of the model. It is written against the abstract spec — the *model* is whatever static object the authoring layer produces; the *engine state* is the live runtime object derived from it — and holds regardless of how that model is stored. Concrete `file:line` references point at the current `ref-agents` engine (`src/solvers.jl`, `src/state.jl`) where the rule is implemented or, where flagged, *violated*.

### 3.1 Vocabulary

- **Place** `s`: a resource pool holding a non-negative quantity `u[s]` — its *marking*. A place may be *plain* (a `Float64` count) or *structured* (backed by individual token agents whose live count is reflected into `u[s]`).
- **Transition** `t`: a stateful *recipe* — not an event — that periodically spawns instances. It carries a spawn `rate`, a `priority` (fill-rate weight, see ADR 0002), a `cycleTime`, a `probOfSuccess` (PoS), a `capacity`, a `maxLifeTime`, a `multiplier`, optional pre/post actions, and an LHS/RHS arc specification with per-token `multiplicity` and `modality`.
- **Instance**: one in-flight firing of a transition (engine type `Transition`, `state.jl:14-26`), holding its parent recipe index `i`, a frozen sampled-attribute snapshot `trans`, its birth time `t`, its multiplicity `q` (how many concurrent firings this instance object represents), and an accumulated progress `state`.
- **Modality**: a per-LHS-token tag governing *when* a token is debited and *whether/how* it is returned (fully specified in §1; the legacy `Set{Symbol}` form, `{:nonblock, :conserved, :rate}`, is currently unvalidated, `ReactiveDynamics.jl:144`).
- **Tick**: one advance of the simulation clock by `dt`.

### 3.2 Instance lifecycle

Every transition instance passes through the following stages exactly once, in order:

1. **Genesis (Poisson spawn).** Each tick, transition `t` proposes `q_desired = ceil(rate · multiplier)` new firings, where `rate` is authored as `rand(Poisson(max(dt · rate, 0)))` (`create.jl:149-151`, `expand_rate`) so the spawn count is a `dt`-scaled Poisson draw. Cycle-time macros `@ct`/`@cycletime` rewrite to `1/arg` (`create.jl:156-159`).
2. **Capacity gate.** The proposal is clamped so that concurrent instances of `t` never exceed `transCapacity` (`solvers.jl:147-154`): `q = min(capacity − liveCount, q_desired)`, and any overflow is *deferred* to a future tick. (The deferral path `add_to_spawn!`, `state.jl:251-256`, is currently broken — see Invariant 3.)
3. **Resource allocation (genesis).** The clamped demand competes for supply via the allocator of **ADR 0002** (weighted progressive filling). Only *upfront* LHS tokens are debited at genesis — modalities `:rate` and `:nonblock` are excluded from the spawn requirement (`solvers.jl:20`, `get_reqs_init!`). Allocations are floored to whole instances (`get_init_satisfied`, `solvers.jl:110-128`), `u` is debited (`solvers.jl:170`), and an instance object is created with `t = clock`, `q = granted count`, `state = 0.0` (`solvers.jl:178-189`). For structured places the granted integer count of token agents is bound to the instance by descending token `priority` (`solvers.jl:194-225`). The recipe's `transPreAction` runs (`solvers.jl:227`).
4. **Cycle-time accumulation.** Each subsequent tick, in-flight instances compete again for *ongoing* resources (`:rate` tokens scaled by `dt` when `cycleTime > 0`, plus `:nonblock` tokens unscaled — `get_reqs_ongoing!`, `solvers.jl:31-48`). The granted fill fraction `q_frac` advances progress: `instance.state += q_frac · dt` (`solvers.jl:261`). Under contention an instance advances *slower than wall-clock*; with full allocation it advances by exactly `dt` per tick.
5. **Terminal test.** An instance terminates when **either** its cycle completes (`state ≥ cycleTime`) **or** its lifetime is exhausted (`clock − birth ≥ maxLifeTime`) (`solvers.jl:408-410`).
6. **Success draw (Binomial PoS).** On termination, the number of *successful* firings is `q_success = rand(Binomial(q, probOfSuccess))` if the cycle completed, else `0` (a lifetime-only timeout yields no successes) (`solvers.jl:412-416`).
7. **RHS emission.** For each successful firing, the RHS products are emitted into `u` at their multiplicity; structured products spawn new token agents; the running reward ledger accumulates `placeReward` (`solvers.jl:418-435`).
8. **Resource return.** LHS tokens are returned according to modality (`solvers.jl:437-483`): `:conserved` tokens return `q · multiplicity · (cycleTime if :rate else 1)` (the resource was *held*, not consumed); `:nonblock` tokens return `q · multiplicity` (held but never blocking). Plain consumed tokens (no return modality) are *not* returned — they were consumed at genesis. `:conserved` together with `:nonblock` is rejected (`solvers.jl:461-465`). Structured tokens are unbound and, for fully-consumed bindings, marked `:removed` (`solvers.jl:443-457, 487-490`).
9. **Termination.** `transPostAction` runs (`solvers.jl:485`) and the instance is pruned from the in-flight set (`solvers.jl:501` — but see Invariant 6).

### 3.3 The ordered step (one tick)

The engine implements one tick as `_step!` (`solvers.jl:642-673`), driven by the host stepping interface. The following order is **normative** — resource accounting, conservation, and determinism all depend on it:

1. **Sync structured counts** — reflect live structured-token agent counts into `u` (`solvers.jl:643`, `update_u_structured!`).
2. **Initial save** — if no history exists yet, record the initial `(t, u)` row (`solvers.jl:644-646`).
3. **Free blocked places** — release `:nonblock` resources held by in-flight instances back into `u` before this tick's allocation (`solvers.jl:648`, `free_blocked_places!`). *(Currently broken — Invariant 1.)*
4. **Update observables** — resample any observable whose `(t − last) ≥ every` (`solvers.jl:650` → `state.jl:144`).
5. **Sample transitions** — clear the per-tick transition table, evaluate each activated recipe's sampleable attributes fresh, and unfold the LHS into arc records (`solvers.jl:651` → `state.jl:174-217`). This is where per-tick values (rate, priority, multiplicity, cycleTime) are realized.
6. **Evolve** (`solvers.jl:652` → `solvers.jl:133-313`): **spawn** new instances (genesis + capacity gate + genesis allocation, lifecycle stages 1-3), then **advance** all in-flight instances (ongoing allocation + progress accumulation, stage 4). Both phases allocate via ADR 0002. Per ADR 0002, priority must be re-read fresh per tick in *both* phases; the current ongoing phase reads a spawn-time snapshot (`solvers.jl:240`) and must be fixed to re-read `transPriority` per tick.
7. **Sync structured counts** (`solvers.jl:653`).
8. **Finish** (`solvers.jl:654` → `solvers.jl:400-508`): for every instance past its terminal test, run the success draw, RHS emission, resource return, post-action, and pruning (lifecycle stages 5-9).
9. **Sync structured counts** (`solvers.jl:655`).
10. **Rules** — evaluate each Rule's guard and fire its action `q` times once per tick (`solvers.jl:657`, `event_action!`; the repaired event channel — see §12). Placement is normative: AFTER `finish!` and the structured-count sync (steps 8-9), BEFORE the ledger row (step 11) and the clock advance (step 12), so a rule's `SetMarking`/`SetParams`/`AddToken` hits the same tick's valuation row and is visible to the next tick's guards and genesis. *(Currently a no-op — Invariant 7; repaired by [ADR 0010](adr/0010-rules-and-conditional-transitions.md).)*
11. **Ledger** — push the `:valuation` row: `u' · placeValuation` (`solvers.jl:659-666`); cost/reward rows were pushed inside `evolve!`/`finish!`.
12. **Advance clock** — `t += dt` (`solvers.jl:668`). This is the *single* authoritative clock advance; all tick work above occurs at the *old* `t`.
13. **Save** — record the post-increment `(t, u)` row (`solvers.jl:670`).

Termination of the run: the engine reports completion once `t > tspan[2]` (`solvers.jl:675`).

### 3.4 Invariants (the contract)

These must hold at every tick boundary (i.e. after step 13). Each is stated as a contract obligation, followed by where the current code honors or violates it.

**1. Resource non-negativity.** `u[s] ≥ 0` for every place `s`, at all times. No allocation may debit a place below zero, and the order of operations (free blocked → allocate → return) must never transiently require negative supply. *The ADR 0002 allocator guarantees this by construction.* **Violation:** `free_blocked_places!` (`solvers.jl:510-522`) references an undefined variable `q` at **`solvers.jl:512`**, throwing `UndefVarError` on any in-flight `:nonblock` LHS token. The `:nonblock`-release path (step 3) is therefore effectively dead: blocked resources that should re-enter the pool each tick may not, and the function errors whenever a `:nonblock` token is in flight. The presence of a `max(0, u[i])` clamp in the legacy allocator (`solvers.jl:68`) is itself evidence that negativity has occurred in practice.

**2. Conservation (conserved tokens returned exactly).** A token tagged `:conserved` is *held* for the instance's lifetime and returned in full on termination — never consumed. The quantity returned must equal the quantity held: `q · multiplicity · (cycleTime if :rate else 1)`. The closed system's conserved mass is invariant across spawn→return. *Honored at* `solvers.jl:437-458`. Caveats: (a) `:conserved + :nonblock` is correctly rejected as ill-defined (`solvers.jl:461-465`); (b) conservation is only exact if Invariant 6 holds — a re-emitting un-pruned instance would return conserved tokens repeatedly, inflating the pool.

**3. Capacity.** The number of concurrent in-flight instances of transition `t` never exceeds `transCapacity`. Overflow is deferred to later ticks, not dropped. *Gate present at* `solvers.jl:147-154`. **Violation:** the deferral helper `add_to_spawn!` (`state.jl:251-256`) is doubly broken — `findfirst` is handed a scalar `length(...)` instead of a range, and on match it increments `:transHash` (`+= n`) instead of `:transToSpawn`. Capacity-overflow deferral is non-functional: overflow is silently lost rather than carried forward, so under sustained over-demand the realized spawn rate is below contract.

**4. Integrality of instance counts.** Spawned instance multiplicities `q` are non-negative integers; structured-token multiplicity must be integer-valued. *Honored:* spawn counts are floored to whole instances (`get_init_satisfied`, `solvers.jl:110-128`); structured multiplicity is checked and errors on non-integers (`solvers.jl:196-200, 268-272`); the Binomial success draw consumes `Int(trans.q)` (`solvers.jl:413`). Caveat: that `Int(...)` cast throws `InexactError` if `q` is ever non-integral, so integrality is *assumed*, not defensively coerced.

**5. Determinism under seed.** Two runs with the same model, inputs, and RNG seed produce identical trajectories. The allocator is deterministic and RNG-free (ADR 0002); all stochasticity is confined to the Poisson spawn draw, the Binomial PoS draw, and observable resampling, all of which must route through a single seeded `AbstractRNG`. **Violation:** every stochastic draw currently uses the global RNG with no seed threaded through the state — Poisson spawn (`solvers.jl:141` via `create.jl:151`), event Poisson (`solvers.jl:321`), Binomial PoS (`solvers.jl:413`), and observable sampling (`state.jl:123`). Determinism-under-seed is *not* currently achievable. ADR 0001's reproducibility obligation requires threading a state-owned RNG through all three sites. (The full seeding obligations are §4.)

**6. Termination completeness (lifetime prune).** Every instance that passes the terminal test is processed exactly once and then removed from the in-flight set — no instance emits its RHS or returns resources more than once. **Violation:** the prune at **`solvers.jl:501`** keeps instances satisfying `state < cycleTime`, which *retains* instances that terminated solely by `maxLifeTime` (their `state` never reached `cycleTime`). Such instances are re-evaluated every subsequent tick — re-running the success draw, re-emitting RHS products, and re-returning conserved/nonblock resources each tick. This violates Invariants 2 (conservation) and 4 indirectly, and is a correctness defect, not merely a leak. The prune predicate must be "remove every instance that passed the terminal test," matching `solvers.jl:408-410`.

**7. Rule firing (the repaired event channel — see §12).** A Rule whose guard holds fires its action `q` times per tick at step 10 (`q` = 1 for a Bool guard, `rand(state.rng, Poisson(v))` for a numeric guard), where the guard is an eval-free `ExprNode` and the action is from the closed set `{SetMarking, SetParams, AddToken, Activate, Deactivate, Log, Seq}` ([ADR 0010](adr/0010-rules-and-conditional-transitions.md)). **Violation (the defect being corrected):** `event_action!` (`solvers.jl:316-326`) computes `q` correctly but at **`solvers.jl:323`** merely *fetches* `state[i, :eventAction]` inside the loop without evaluating it — the action expression is never run, so events are currently a complete no-op and any behavior depending on them (scheduled budget injections, the endogenous acquisition lever, conditional what-if interventions) does not execute. The repair (ADR 0010 §A) replaces the bare fetch with `context_eval(state, nothing, state.wrap_fun(action_i))`, RNG-threaded per §4 D5. Conditional *transitions* are the sibling mechanism: a stateless `guard::ExprNode` AND-ed with the latching `transActivated` gate (`state.jl:179`) — see §12.

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

> **D4 (Order determinism).** Wherever a stochastic draw is made per object (per transition, per place, per event, per observable), the iteration order over those objects MUST be deterministic and stable across runs (e.g. iterate `parts(state, :T)` / `parts(state, :S)` in index order, never `Dict`/`Set` insertion order). The number of draws consumed and the order in which they are consumed are part of the contract.

### 4.2 Where the RNG must be threaded

Today every stochastic draw uses the implicit global RNG. The contract requires a single `rng::AbstractRNG` field on the run state, created at construction, and passed explicitly to every draw. The complete current inventory of global-RNG sites that MUST be converted:

| Draw | Current site | Distribution | Purpose |
| --- | --- | --- | --- |
| Transition spawn count | `expand_rate`, `src/interface/create.jl:151` (`rand(Poisson(max(state.dt * rate, 0)))`); realized in `evolve!`, `src/solvers.jl:140-144` | `Poisson(dt·rate)` | how many new instances of a transition to schedule this tick |
| Probability of success | `finish!`, `src/solvers.jl:413` (`rand(Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))`) | `Binomial(q, PoS)` | how many of `q` completing instances succeed |
| Observable sampling | `sample_range`, `src/state.jl:123` (`rand() * sum(...)`) and `src/state.jl:132` (`rand(r)`) | `Uniform` selector + inner `Sampleable` | choosing and sampling an observable's range entry |
| Event firing count | `event_action!`, `src/solvers.jl:321` (`v isa Number ? rand(Poisson(v)) : 0`) | `Poisson(rate)` | how many times an event's action fires this tick |
| Generic sampleable eval | `context_eval`, `src/state.jl:70` (`o isa Sampleable ? rand(o) : o`) | any `Distributions.Sampleable` | evaluating any attribute that resolves to a distribution |

The last row is the most important: `context_eval` (`src/state.jl:67-71`) is the single chokepoint through which nearly every attribute value (rates, multiplicities, costs, etc.) flows. Threading the RNG through `context_eval(state, transition, o)` — i.e. `rand(state.rng, o)` — covers the majority of draws by construction. The four explicit `rand(...)` call sites above MUST additionally take `state.rng` as their first argument.

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

### 4.6 Future work — per-entity RNG substreams (the paired-counterfactual gap)

> **D-future (per-entity substreams).** The contract today specifies ONE state-owned RNG stream per trajectory (D2/D5) and independent per-MEMBER streams across an ensemble (D8/D9), but NO per-entity substreams WITHIN a trajectory. This is a known limitation surfaced by the BD acquisition counterfactual ([MVP_BD_DEMO.md](../demo/bd_acquisition/MVP_BD_DEMO.md) finding A): two runs that differ only by an intervention (e.g. a §12 acquisition rule injecting programs) diverge in the NUMBER of draws consumed after the intervention, so the single stream desynchronises and every unrelated entity's draws shift. The mean treatment effect is still unbiased under ensemble averaging (D8), but common-random-numbers variance reduction is unavailable, so a per-seed PAIRED Δ is contaminated and tight CIs need many more seeds. **Decision (2026-06-21):** keep the single stream for now; the demo uses ensemble-averaged Δ. **If a paired counterfactual is later required**, the fix is per-entity substreams keyed by `(root_seed, entity_id)` (per transition and/or per structured token), so injecting an entity perturbs only its own and genuinely-downstream draws. This would extend D5/D8 and interacts with the §10.5 `dump_state`/`restore` fork-at-tick checkpoint (which already enables a partial paired comparison by forking a single run at the lever tick). Recorded as future work, not a Phase-1 obligation.

---

## 5. Attribute Contract (typed attributes)

### 5.1 Motivation and the replacement

Today every authored quantity has the catch-all type `SampleableValues = Union{Expr,Symbol,AbstractString,Float64,Int,Function}` (`src/ReactiveDynamics.jl:10`), used for all of `transRate`, `transPriority`, ..., `placeValuation` (`src/ReactiveDynamics.jl:46-62`). This type encodes nothing about units, sign, range, or whether time-variation is permitted, so the only validation is whatever the solver happens to do at runtime. The attribute contract replaces this with a per-attribute specification: each attribute has a name, a **kind**, a **scalar domain** (units + valid range + sign), a **default**, and a **time-variation policy**.

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

### 5.4 Place attributes

| Attribute | Units | Valid range | Default | TVE? | Notes / current behavior |
| --- | --- | --- | --- | --- | --- |
| `placeInitVal` | tokens (the place's marking) | ≥ 0 | `0.0` (`src/ReactiveDynamics.jl:131`) | no (read once) | Initial `u`. Used at construction (`src/solvers.jl:561-569`) and re-init (`src/solvers.jl:620`); `u0` keyword overrides per place. Stored as `Float64` but represents a token count, so SHOULD be a nonnegative integer for unstructured places. |
| `placeInitUncertainty` | same units as `placeInitVal` (absolute) or dimensionless (relative) — MUST be fixed by the contract; recommend **relative, dimensionless ≥ 0** | ≥ 0 | `0.0` (`src/ReactiveDynamics.jl:130`) | no (read once) | Spread applied to `placeInitVal` when sampling initial conditions for an ensemble. **Currently declared and defaulted but not consumed by the solver** (no read site in `solvers.jl`/`state.jl`); the contract pins its semantics (interpretation + the distribution used) so the determinism contract §4 governs the draw. |
| `placeCost` | value · token⁻¹ (currency per unit consumed) | ≥ 0 (typical); finite | `0.0` (`src/ReactiveDynamics.jl:132`) | yes | Cost charged per allocated unit. Logged each tick as `actual_allocs' · placeCost` (`src/solvers.jl:308-311`). For BD/rNPV ledgers this is the spend side. |
| `placeReward` | value · token⁻¹ (currency per unit produced) | ≥ 0 (typical); finite | `0.0` (`src/ReactiveDynamics.jl:133`) | yes | Reward credited per produced unit on transition completion (`src/solvers.jl:426,433`). The income side of the ledger. |
| `placeValuation` | value · token⁻¹ (currency per unit held) | finite (may be signed) | `0.0` (`src/ReactiveDynamics.jl:134`) | yes | Mark-to-market value of the current holding; logged each tick as `u' · placeValuation` (`src/solvers.jl:663-664`). May legitimately be negative (e.g. a liability), so the contract does NOT impose ≥ 0 here. |
| `placeDefaultModality` | — (set of tags) | ⊆ {`:nonblock`, `:conserved`, `:rate`} | `Set{Symbol}()` (`src/ReactiveDynamics.jl:154`) | no | The place-level modality set, unioned with per-arc modality during sampling (`src/state.jl:204`). The allowed universe is `place_modalities` (`src/ReactiveDynamics.jl:144`); the contract requires validation against this set (today it is an unvalidated `Set{Symbol}`). `:conserved` + `:nonblock` together is illegal and MUST be rejected at authoring (today only caught at completion, `src/solvers.jl:461-465`). See §1 for the orthogonalized re-model of these tags. |
| `placeStructured` | — (flag) | `Bool` | `false` (`src/ReactiveDynamics.jl:135`) | no | Marks a place as a structured token (agent-backed) rather than a plain count. Determines `structured_token_names` (`src/solvers.jl:556-557`). Listed for completeness; not numeric. |

### 5.5 Cross-cutting rules

> **A1 (Sign/range enforcement).** Every range in §5.3–5.4 MUST be enforced: at construction for literals; at each read for TVEs. The only sanctioned silent clamps are the existing `max(dt·rate, 0)` for spawn intensity (`src/interface/create.jl:151`) and the `max(0, u[i])` guard in `alloc_weighted!` (`src/solvers.jl:68`); all other violations MUST error.

> **A2 (Units are advisory but fixed.)** The engine is unit-agnostic at runtime (all quantities are `Float64`), but the contract FIXES the unit interpretation per attribute so that authored models, ledgers, and observables compose consistently. `transRate`'s time unit MUST match the model's `tunit` (`get_tcontrol`, `src/solvers.jl:526-534`); `placeCost`/`placeReward`/`placeValuation` MUST share one currency unit for the valuation logs (`src/solvers.jl:308,663`) to be meaningful.

> **A3 (TVE policy).** Attributes marked "TVE? = yes" MAY be authored as an `Expr`/`Function`/`Sampleable` and are re-read each tick through `context_eval`; attributes marked "no" are read once (construction or re-init) and frozen for the run. Authoring a TVE for a "no" attribute MUST be rejected.

> **A4 (Integrality).** `placeInitVal` (unstructured), `transCapacity` (finite), and any multiplicity feeding a structured place are compared against or converted to integers at runtime (`Int(trans_.q)` at `src/solvers.jl:413`; the `isinteger` checks at `src/solvers.jl:196-200,268-272`). The contract requires these to be integer-valued; non-integer literals MUST be rejected at construction.

> **A5 (Replacement of the catch-all type).** `SampleableValues` (`src/ReactiveDynamics.jl:10`) is replaced as the *semantic* contract by the per-attribute kinds above. The underlying ACSet `AttrType` MAY remain a broad union for storage, but the validation layer defined here is authoritative: a model that type-checks against the ACSet but violates §5.3–5.5 is invalid.

---

## 6. Object Model

This section pins the canonical typed representation of a model under the ADR 0003 dependency-free typed IR. It defines each object as a typed record with named fields, a scalar or [ADR 0005](adr/0005-serialization-json-ir.md) `ExprNode` type per field, and the [ADR 0004](adr/0004-runtime-mutation.md) append-only index discipline. It is the *static authoring/IR* layer only: the live runtime `ReactionNetworkProblem` (`state.jl:40-63`) and its in-flight `Transition`/`Observable` instances (`state.jl:14-38`) are derived from this model and are not part of it. The model is a **struct of columnar object-tables**: each object below is one table whose rows share a single integer index space, exactly mirroring the six current ACSet objects (`:S, :T, :E, :obs, :P, :M`, `ReactiveDynamics.jl:30-76`) plus the one new promoted table (`ArcSpec`). The defining change versus the current source is structural: the transition↔place relation (the net's **arcs**), today *not stored* and re-parsed from the `trans` `Expr` every tick by `extract_arcs` (`interface/reaction_parser.jl:32`, called at `state.jl:194` and `solvers.jl:418`), becomes a first-class typed incidence table with schema-enforced foreign keys.

### 6.1 Field-type vocabulary

Three field categories appear in the records below:

- **Scalar** — a literal stored once and frozen for the run (e.g. `init::Float64`, `structured::Bool`). Never re-read through `context_eval`.
- **`ExprNode`** — a time-varying expression encoded as the closed, eval-free, JSON-representable sum type of [ADR 0005](adr/0005-serialization-json-ir.md) (`Const`/`Ref{place|param|obs}`/`Call{op∈OP_WHITELIST}`/`Sample{dist∈DIST_WHITELIST}`/`TimeRef`/`Choose`). A `Const` ExprNode is the canonical form of a literal that *may* be a TVE; a non-trivial tree is a genuine TVE re-evaluated each tick. This replaces today's catch-all `SampleableValues = Union{Expr,Symbol,AbstractString,Float64,Int,Function}` (`ReactiveDynamics.jl:10`) for every attribute that §5 marks "TVE? = yes."
- **Closed tag** — a `Symbol` constrained to a fixed enumerated set, validated at construction (e.g. `genesis::GenesisKind`, the §1 `Modality` axes, `side∈{lhs,rhs}`). These are JSON-Schema-enumerable per ADR 0005.

Per-field units, ranges, defaults, and TVE policy are normative in **§5** (Attribute Contract) and are not restated here; this section fixes the *record shape*, the *typing*, and the *identity/index* rules. Field names below are the contract's canonical names; the parenthetical is the current schema column they correspond to.

### 6.2 The append-only index invariant (ADR 0004)

Every object table is a positional column store: row `i` of an object is addressed by its integer index, and that index is the object's identity within the live runtime (compiled attribute closures hard-code `state.u[i]` against a varmap frozen at construction, `compilers.jl:149`). Therefore the model obeys [ADR 0004](adr/0004-runtime-mutation.md) INV-1/INV-2/INV-3 verbatim:

- **Append-only.** A new place, transition, param, event, observable, or arc row receives the next free index (`nparts+1`); existing indices NEVER move. Growth is monotone.
- **No mid-run reindex.** Deletion and reordering of rows are forbidden on a live/stepping model. `rem_parts!` (the sole reindexer, `operators/equalize.jl:52`) must refuse under a runtime `live` guard.
- **Soft-deactivate, not delete.** Logical removal of a transition is `transActivated[i]=false` (honored as a skip-gate at `state.jl:179`); the row, its index, and all higher indices stay put, and in-flight instances of a deactivated transition still run to completion (`finish!` iterates `ongoing_transitions` independently of `transActivated`, `solvers.jl:406`).

The live mutation API (`add_place!`/`add_transition!`/`add_param!`/`activate!`/`deactivate!`) and the `refresh_wrap_fun!` re-derivation that keeps pre-existing indices bound to the same `state.u[i]` are specified in ADR 0004 and are not restated here.

### 6.3 Identity-by-name and its risk

Objects are matched and merged **by name**, not by index, at the authoring/composition boundary. The current implementation keys every lookup, merge, and substitution on the name column with no enforced uniqueness:

- `merge_acs!` adds a place only if `incident(acs, r, :placeName)` is empty (`ReactiveDynamics.jl:211-214`) — first-write-wins by name.
- `union_acs!` (`operators/joins.jl:14-26`) merges two networks by looping `incident(acs1, name, :placeName)` + `add_part!`.
- `equalize!` (`operators/equalize.jl`) merges places by name, with a fragile fallback (`equalize.jl:55-63`) that rewrites place names via `recursively_substitute_vars!` over EVERY attribute `Expr`.

**Contract rule (identity).** Within one model table, `name` MUST be unique; the constructor/validator MUST reject duplicate `Place.name`, `Transition.id`, `Param.name`, and `Observable.name`. **Risk this pins.** Because composition is name-keyed string surgery rather than a structural reference repoint, a place name that collides inside a multiplicity/rate subexpression can be silently corrupted by `recursively_substitute_vars!` (ADR 0003 Context). The promoted `ArcSpec` table (§6.5) is the structural fix: composition repoints integer FKs instead of substituting strings, making place-merge exact. The full composition semantics are a separate (pending) contract section; this section only fixes the identity rule the merge operators must honor.

### 6.4 Place

A resource pool. One row per place; the row index is the position of `state.u[i]`, i.e. that place's marking.

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `name` | Symbol (closed identity) | `placeName` (`ReactiveDynamics.jl:44`) | Unique key; the FK target for `ArcSpec.place` (§6.5). |
| `init` | `ExprNode` (read-once → `Const` literal expected; §5 marks "no") | `placeInitVal` (`:46`) | Initial `u[i]`; integer-valued for unstructured places (§5 A4). |
| `init_uncertainty` | `ExprNode` (read-once) | `placeInitUncertainty` (`:47`) | Ensemble initial-condition spread; semantics fixed in §5.4 (declared but currently unconsumed by the solver). |
| `cost` | `ExprNode` (TVE) | `placeCost` (`:48`) | Currency per unit consumed. |
| `reward` | `ExprNode` (TVE) | `placeReward` (`:49`) | Currency per unit produced. |
| `valuation` | `ExprNode` (TVE) | `placeValuation` (`:50`) | Mark-to-market value per unit held; may be signed (§5.4). |
| `modality` | `Modality` (the §1 orthogonal form: `allocation∈{upfront,perstep}` × `return∈{consumed,conserved}` × `blocking∈{block,nonblock}`) | `placeDefaultModality` (`:45`) | Replaces the unvalidated `Set{Symbol} ⊆ {:nonblock,:conserved,:rate}` (`ReactiveDynamics.jl:144`). This place-level modality is the default unioned with each per-arc LHS modality during sampling (`state.jl:202`); see §6.5. Legal rows and the illegal `nonblock⇒consumed` / `perstep` rules are §1.3–§1.4. |
| `structured` | Bool (closed tag) | `placeStructured` (`:51`, default `false`) | Marks the place as agent-backed (a structured token) rather than a plain `Float64` count. See §6.10. |

### 6.5 ArcSpec — the promoted incidence table (ADR 0003 Phase 2)

The model's defining relation: which place each transition consumes (LHS) or produces (RHS), with multiplicity and (LHS only) modality. One row is one **arc**; `multiplicity` is its **arc weight** (inscription). **This table is the structural promotion mandated by ADR 0003 Phase 2** — it replaces the runtime re-parse of the `trans` `Expr` by `extract_arcs` (`interface/reaction_parser.jl:32`), which today reconstructs `FoldedArc(place, multiplicity, modality)` rows (`interface/reaction_parser.jl:5-9`) on the fly at `state.jl:194` and `solvers.jl:418`. Promoting it to a stored typed table delivers what REVIEW.md flagged as missing — "structure is implicit, not morphic" — by making the bipartite place↔transition graph first-class, FK-checkable, and round-trippable (the serialized array is spelled `arcs[]`, ADR 0005 as amended by ADR 0017).

One row per (transition, place, side) participation:

| Field | Type | Notes |
|---|---|---|
| `transition` | FK → `Transition.id` (integer index) | Schema-enforced: validation MUST reject a row whose `transition` is not a live transition index (ADR 0005 validate rule 2, dangling FK). |
| `place` | FK → `Place.name`/index | Likewise schema-enforced; a dangling `place` FK is invalid. |
| `side` | closed tag `∈ {lhs, rhs}` | LHS = consumed input; RHS = produced output. Mirrors `prune_r_line`'s `(l_line, r_line)` split (`state.jl:152-155`). |
| `multiplicity` | `ExprNode` (TVE) | Per-participation multiplicity (the arc weight); integer-valued where it feeds a structured place (§5 A4, checked at `solvers.jl:196-200,268-272`). The `@structured`/`@move` dynamic idioms (`interface/reaction_parser.jl:67`) and the separate `@choose` idiom (`recursively_choose`, `interface/reaction_parser.jl:11-30`) become explicit `ExprNode`/typed-arc variants (`Choose`, structured-token nodes) so this table stays the single source of truth. (`@structured` genesis is DONE: it is now the typed, registry-resolved `structured{kind, fields}` arc that round-trips through the JSON IR — the raw inline-constructor form was removed; see ADR 0005 §21. `@move`/`@choose` remain Expr-level pending the Phase-2 `ArcSpec` promotion.) |
| `modality` | `Modality` (§1 orthogonal form) — **LHS rows only** | Present only when `side = lhs` (a modality is "the contract attached to each LHS consumed token," §1; the current `FoldedArc.modality` is likewise LHS-only, `reaction_parser.jl:8`). RHS rows carry no modality. The effective modality of an LHS token is this per-arc value unioned with the place-level `Place.modality` (today `r.modality ∪ state[j,:placeDefaultModality]`, `state.jl:202`). |

**Why this matters (the three wins ADR 0003 names).** (1) *Schema-enforced FKs* — dangling participations are caught by `validate` at construction/load rather than surfacing as a runtime `find_index` returning `nothing`. (2) *Structural validation* — the LHS/RHS bipartite structure (a real Petri-net incidence) is inspectable without parsing an `Expr`. (3) *Exact composition* — `union_acs!`/`equalize!` merge places by repointing integer `place` FKs, not by `recursively_substitute_vars!` string rewriting (§6.3 risk), so a place-name collision inside a subexpression can no longer corrupt incidence. Append-only (§6.2) applies to this table too: a mutation adds arc rows at the next index and never reorders them.

### 6.6 Transition

A stateful spawning recipe (not an event). One row per transition; the row index/`id` is the FK target for `ArcSpec.transition` and the position the per-tick sampler walks (`state.jl:178`).

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `id` | Symbol/integer (closed identity) | row index of `:T` | Stable unique key; FK target for `ArcSpec`. For reproducible ensembles an injected transition SHOULD get an explicit stable `id`/`name` (ADR 0004: `missing` falls back to `gensym()`). |
| `name` | Symbol/String/`missing` (descriptive scalar) | `transName` (`:63`, default `missing`) | Label only; not a quantity. |
| `rate` | `ExprNode` (TVE) | `transRate` (`:55`) | Spawn intensity (expected instances per unit time); required at authoring (§5.3). Lowers to the genesis draw per `genesis`/`rate_mode` (§6.7). |
| `genesis` | closed tag `GenesisKind ∈ {poisson, scheduled, flow, capacity}` | NEW append-only field (§2.8) | Intent declaration over the shared two-stage genesis path; `poisson`/`scheduled` are sources, `flow`/`capacity` are routing and REQUIRE ≥1 upfront-consumed LHS `ArcSpec` row (§2.8 validation). Engine default `poisson`; BD/pipeline template authors routing as `flow` (§2.8). Subsumes the ADR 0005 `rate_mode∈{poisson,deterministic}` discriminator. |
| `priority` | `ExprNode` (TVE, re-read per tick per ADR 0002) | `transPriority` (`:54`, default `1`) | Per-tick fill-rate weight for weighted progressive filling; `0` = leftover-only. |
| `cycletime` | `ExprNode` (TVE) | `transCycleTime` (`:56`, default `0.0`) | Service duration; gates `perstep` resource draw on `> 0` (§1.4, §2.4). |
| `prob_of_success` | `ExprNode` (TVE, range `[0,1]`) | `transProbOfSuccess` (`:57`, default `1`) | PoS for the `Binomial(q, p)` success draw (`solvers.jl:413`). |
| `capacity` | `ExprNode` (TVE, integer or `Inf`) | `transCapacity` (`:58`, default `Inf`) | Max concurrent live instances; overflow deferred (ADR 0004; deferral currently broken — §3.4 Inv 3). |
| `max_lifetime` | `ExprNode` (TVE, `≥0` or `Inf`) | `transMaxLifeTime` (`:59`, default `Inf`) | Wall-clock age cap (§2.5). |
| `multiplier` | `ExprNode` (TVE, `≥0`) | `transMultiplier` (`:62`, default `1`) | Scales the per-tick spawn target. |
| `guard` | `ExprNode` (TVE, resolves to `Bool`; default `Const(true)`) | NEW append-only field (§12.2, [ADR 0010](adr/0010-rules-and-conditional-transitions.md)) | Stateless per-tick firing condition, AND-ed with the latching `transActivated` gate (`state.jl:179`): `fires = transActivated[i] && eval(guard)`. A guarded-off transition makes no genesis proposal that tick. Non-`Bool` is a `validate()` diagnostic; an RNG-consuming guard SHOULD warn. The conditional-transition half of the §12 endogenous decision channel. |
| `pre_action` | action-statement tree (ADR 0005/§12.3 `{SetMarking,SetParams,AddToken,Activate,Deactivate,Log,Seq}`) | `transPreAction` (`:60`, default `:()`) | Side-effecting code run on spawn (`solvers.jl:227`); not a numeric attribute, never sampled through `context_eval`'s draw path. |
| `post_action` | action-statement tree | `transPostAction` (`:61`, default `:()`) | Run on completion (`solvers.jl:485`). |

The LHS/RHS arc specification is NOT a field here — it lives in the `ArcSpec` table (§6.5), which is the whole point of the Phase-2 promotion. The runtime in-flight instance (engine type `Transition`, `state.jl:14-26`, carrying `i`, frozen `trans` snapshot, birth `t`, multiplicity `q`, progress `state`) is a runtime object derived from this recipe row, not a model object.

### 6.7 Genesis as a recipe property

`genesis` (§6.6) is a closed tag over the single two-stage execution path (rate-proposal then upfront-LHS resource gate; §2.8), not a separate object. It is append-only (ADR 0004) and JSON-Schema-enumerable (ADR 0005). `poisson` lowers `rate` to `Poisson(dt·rate)`; `scheduled` to the `@deterministic` bare-count + calendar idiom; `flow`/`capacity` ride the existing upfront-LHS gate and require a consuming LHS (§2.8 validation, enforced against the `ArcSpec` table). `batch` is a transition parameter, not a genesis mode.

### 6.8 Event → Rule (the endogenous decision channel — see §12, [ADR 0010](adr/0010-rules-and-conditional-transitions.md))

The `:E` object is a **Rule**: a guarded action evaluated once per tick at step 10 (the BD acquisition lever's home, with `event_action!` repaired — §3.4 Inv 7). The rule half of the §12 endogenous decision channel (the conditional-transition half is the §6.6 `guard` field). One row per rule.

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `id` | Symbol (closed identity) | row index of `:E` | Unique key. |
| `guard` | `ExprNode` (TVE; resolves to Bool or numeric) | `eventTrigger` (`:65`) | Bool ⇒ fire once; numeric `v` ⇒ fire `Poisson(v)` times per tick (`solvers.jl:320-321`), RNG-threaded (§4 D5). May be scheduled (`@t()`) or state-contingent (place/param/obs/`TokenAgg` reads under the §9.5 measurable point). |
| `action` | action-statement tree (§12.3 `{SetMarking,SetParams,AddToken,Activate,Deactivate,Log,Seq}`) | `eventAction` (`:66`) | Run when the guard holds. The repair (ADR 0010 §A) evaluates this (`solvers.jl:323` was a no-op fetch — §3.4 Inv 7). `AddToken` is the in-model acquisition; `Activate`/`Deactivate` toggle pipeline lines. |
| `fire_mode` | closed tag `∈ {every_tick, once}` (default `every_tick`) | NEW (§12.4) | `once` self-disables after first firing; `enabled` is seeded run-state reset by `reinit!` (§4 D7). |
| `enabled` | Bool (run-state latch) | NEW (§12.4) | Internal; `true` initially, set `false` by a fired `once` rule. |

### 6.9 Observable, Param, Meta

**Observable** — a periodically-resampled derived quantity. One row per observable.

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `name` | Symbol (closed identity) | `obsName` (`:68`) | Unique key; referenceable via `Ref{obs}` ExprNodes. |
| `range` | `Vector` of weighted `ExprNode` alternatives | `obsOpts.range` (`FoldedObservable`, `ReactiveDynamics.jl:24-28`) | Sampled per §4 determinism. |
| `every` | Float64 (scalar) | `obsOpts.every` (default `Inf`) | Resample period; resampled when `(t − last) ≥ every` (`state.jl:144-148`). |
| `on` | `Vector{ExprNode}` | `obsOpts.on` | Gate expressions. |

NOTE the runtime `Observable` agent (`state.jl:31-38`) carries a `sampled` field; the documented `resample!` bug writes a nonexistent `o.val` (`state.jl:137` vs field `.sampled` at `state.jl:37`) — a model-independent engine defect the test suite must pin.

**Param** — a named scalar constant. One row per param.

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `name` | Symbol (closed identity) | `prmName` (`:71`) | Unique key; FK target for `Ref{param}` ExprNodes; resolved positionally to `state.p[:name]`. |
| `value` | Float64/number (scalar) | `prmVal` (`:72`, default `missing`) | A JSON number under ADR 0005 — directly replacing the `eval(Meta.parseall(prmVal))` import-time RCE (`loadsave.jl:65`). Not an `ExprNode` tree. |

**Meta** — the keyword bag. ADR 0005 collapses the `:M` keyword rows and the solver kwargs into one `meta` object (they already merge into a single `keywords` bag, `solvers.jl:544-552`).

| Field | Type | Corresponds to | Notes |
|---|---|---|---|
| `key` | Symbol (closed identity) | `metaKeyword` (`:73`) | e.g. `name`, `tspan`, `dt`/`tstops`, `tunit`, `seed`, `alloc_strategy`. |
| `value` | scalar (number/string) | `metaVal` (`:74`, default `missing`) | A literal; replaces the `eval`d-string `metaVal` path. The `seed` here is the §4 determinism seed; `tspan`/`dt` feed `get_tcontrol` (§2.2). |

### 6.10 Structured-token refinement of Place

A place with `structured = true` (§6.4) is *not* a separate object table — it is a refinement of a Place row whose pool is backed by individual token agents rather than a `Float64` count. In Coloured-Petri-net terms the token attributes are its *colour set*. The agent type is `AbstractStructuredToken <: AbstractAlgebraicAgent` (`agents.jl:6`), with the default concrete implementation `BaseStructuredToken` (`@aagent FreeAgent`, `agents.jl:11-15`) carrying `place::Union{Nothing,Symbol}`, `bound_transition::Union{Nothing,Transition}`, and `past_bonds`. Live agent counts are reflected into `u[i]` by `update_u_structured!`, and genesis binds granted integer token counts to an instance by descending token `priority` (`agents.jl:67`). For the object model this means: `structured` is a closed-tag flag on the Place record (not a new table); the contract requires structured places to have integer-valued `init` and integer multiplicity on every incident `ArcSpec` row (§5 A4); and `perstep` allocation modality is illegal for a structured place (§1.4 — you cannot reserve a fractional `dt`-scaled slice of an indivisible agent). The append-aware handling of a newly added structured place (`push!(state.structured_token, name)`) is an ADR 0004 open question, not a model-shape question.

Two contract points pinned by [ADR 0006](adr/0006-structured-tokens.md) (full treatment in §9): (1) a structured place's `state.u` column is a DERIVED reflected count (`update_u_structured!`, `solvers.jl:630-639`), not authoritative storage — the source of truth is the set of live unblocked inner-agent instances of that kind, so instances are added/removed live via `entangle!`/`disentangle!` (index-safe, never a `state.u` reindex), while the KIND stays defined a priori. (2) The token TYPE, its CONSTRUCTOR, and any custom protocol-method override are HOST JULIA CODE (compiled in the user's package, never serialized); a JSON model REFERENCES a kind (and any custom function) BY NAME against a per-network registry, and an unregistered name is a `validate()` dangling-reference diagnostic (§8.3) — never an `eval`. This is the eval-free replacement for `@register`'s module-scope eval (`update.jl:513`).

### 6.11 Summary of the object tables

Seven typed columnar tables, all under the §6.2 append-only index discipline: **Place** (§6.4), **Transition** (§6.6, carrying the `genesis` tag §6.7), the promoted **ArcSpec** incidence table (§6.5 — the structural upgrade), **Event** (§6.8), **Observable**, **Param**, **Meta** (§6.9). Identity is by name within each table with enforced uniqueness (§6.3); the cross-table relation is the FK-checked `ArcSpec` bipartite graph. Per-field domains are §5; the orthogonal `Modality` is §1; `GenesisKind` is §2.8; the `ExprNode`/action-statement IR and JSON serialization are [ADR 0005](adr/0005-serialization-json-ir.md); the index/mutation discipline is [ADR 0004](adr/0004-runtime-mutation.md); the drop-ACSets + promote-the-relation decision is [ADR 0003](adr/0003-data-store.md).

---

## 7. Composition Semantics

**Status.** This section specifies how two models are combined into one. It contracts BOTH the current behaviour (`operators/joins.jl` `union_acs!` / `@join`; `operators/equalize.jl` `equalize!` / `@equalize`) and the target behaviour under the ADR 0003 promoted `ArcSpec` incidence table. Composition today is manual name-matching over `add_part!`/`incident` loops, NOT a categorical pushout; the contract keeps that operational model but fixes its gaps and makes place identification structurally exact.

### 7.1 Vocabulary

A **join** combines two (or more) models into one new model, taking the disjoint union of their parts and then *namespacing* one model's places so the two name-spaces do not accidentally collide. **Identification** (a.k.a. *equalize*) is the inverse pressure: it declares that two places in the (possibly already-joined) model are the *same* place and collapses them into one. **Namespacing** rewrites a place name `X` of model `m` to the qualified name `m__X` (`normalize_name`, `joins.jl:100-102`). The composed model is itself a model and is a valid input to a further join (the operation is closed over models).

### 7.2 The JOIN operation (`union_acs!`, `@join`)

`union_acs!(acs1, acs2, name, eqs)` (`joins.jl:10-59`) merges `acs2` into `acs1` in place and returns `acs1`. `@join m1 m2 …` (`joins.jl:201-238`) folds left over a fresh empty `ReactionNetworkSchema()`: it calls `union_acs!(acs_new, mᵢ, :mᵢ, eqs)` once per model in source order, threading the parsed equation blocks `eqs` so identification happens *during* the fold rather than as a separate pass.

**J1 (Places are identified by name; everything else is namespaced).** Before merging, `union_acs!` calls `prepend!(acs2, name, eqs)` (`joins.jl:12`, `64-82`), which renames every place of `acs2` to `name__placeName` UNLESS an equation block in `eqs` maps it to a shared alias. A place in `acs2` is then merged into an EXISTING place of `acs1` iff their (post-namespacing) `placeName` matches (`incident(acs1, …, :placeName)`, `joins.jl:15`); otherwise it is added as a new `:S` part (`joins.jl:17-20`). Modality sets are *unioned*, not overwritten (`union!(…, :placeDefaultModality)`, `joins.jl:22`). Because each model gets its own `m__` prefix, two models' independent places named `A` become `m1__A` and `m2__A` and do NOT merge — sharing is opt-in via §7.4, never accidental.

**J2 (Transitions are always disjoint; never deduplicated).** All of `acs2`'s `:T` parts are appended unconditionally (`add_parts!(acs1, :T, nparts(acs2, :T))`, `joins.jl:30-36`), and each transition is given a namespaced `transName` of the form `name__transName` (`joins.jl:38-44`); an UNNAMED transition (its `transName` is `missing`) is namespaced as `name__<integer-index>` via `normalize_name(Symbol(coalesce(acs1[i,:transName], i)), name)` (`joins.jl:41`). There is no name-based identification of transitions — joining a model with itself yields two copies of every transition. The contract pins this: transitions are *structural* and are NEVER merged by name. (Their arc references to places ARE rewritten so they point at the merged place; see §7.4.)

**J3 (What is merged).** `union_acs!` reflects over `propertynames(acs.subparts)` and merges parts whose attribute names contain the substrings `"place"` (place attributes, `joins.jl:24-27`) and `"trans"` (transition attributes, `joins.jl:31-36`), plus params `:P` by `prmName` (`joins.jl:46-50`) and meta `:M` by `metaKeyword` (`joins.jl:52-56`). The merged part set is therefore exactly **{S, T, P, M}**.

**J4 (What is NOT merged — a gap the rework MUST close).** Events `:E` (`eventTrigger`/`eventAction`, `ReactiveDynamics.jl:65-66`) and observables `:obs` (`obsName`/`obsOpts`, `ReactiveDynamics.jl:68-69`) are silently dropped from the joined model: neither the `"place"`/`"trans"` reflection loops nor the `:P`/`:M` loops touch them, and there is no `:E`/`:obs` loop. `prepend_obs` (`joins.jl:87-97`) — the function that would namespace observable references inside expressions during a join — is **defined but has zero callers** (dead code). The contract flags this as a defect: a join MUST also merge `:E` and `:obs` (events appended disjointly like `:T` per J2; observables identified by namespaced `obsName` like places per J1), and observable references inside merged rate/multiplicity/action expressions MUST be namespaced (the live successor of `prepend_obs`). Until fixed, the documented behaviour is "events and observables are lost on join," and this is a bug, not a feature.

### 7.3 Precedence

**J5 (Last model wins for attribute values; structure accumulates).** For a merged place, every place attribute of `acs2` overwrites `acs1`'s value unless the incoming value is `missing` (`!ismissing(acs2[i, attr]) && (acs1[…] = acs2[i, attr])`, `joins.jl:26`); likewise `prmVal` (`joins.jl:49`) and `metaVal` (`joins.jl:55`) overwrite when present. This matches the `union_acs!` docstring "the attributes in `acs2` taking precedence" (`joins.jl:7-8`) and the `@join` docstring "the last model takes precedence" (`joins.jl:193`). The contract pins this **last-writer-wins on scalar attributes / accumulate on structure (T always; S/P/M by name)** rule. One consequence: `placeDefaultModality` is the exception — it is unioned (J1), so a modality tag set on either side survives regardless of order.

### 7.4 Place IDENTIFICATION / equalize

Identification declares that named places are the *same* and collapses them to a single part. It arrives two ways: inline in a `@join` via equation blocks routed through `prepend!`/`normalize_name` (`joins.jl:104-128`), or post hoc via `@equalize` / `equalize!` (`equalize.jl`). The equation-block grammar (`@catchall(A)`, `@alias`, `m.X = m'.Y`, bare `X = Y`) is parsed by `get_eqs`/`expand_name` (`joins.jl:132-186`) and `get_eqs_ff` (`equalize.jl:3-22`).

**J6 (Identification = collapse to one part + repoint every reference).** `equalize!` (`equalize.jl:24-66`) does this in two moves: (a) for each equation block it finds all matching place indices, copies any `missing` attribute on the surviving (lowest-index) part from the others (`equalize.jl:46-51`), then **deletes the redundant place parts** via `rem_parts!(acs, :S, place_ixs[2:end])` (`equalize.jl:52`); and (b) rewrites every *other* attribute expression in the model so references to the deleted names point at the survivor, via `recursively_substitute_vars!` over each attribute `Expr` (`equalize.jl:55-63`). `prepend!`'s inline path does the analogous rewrite during a join (`joins.jl:71-79`).

**J7 (Why the promoted `ArcSpec` table makes this exact).** The transition↔place relation (the net's arcs) is NOT stored as structure today — it lives inside the `:trans` attribute as an `Expr` (`ReactiveDynamics.jl:53`) and is re-parsed per tick (`extract_arcs`, `interface/reaction_parser.jl:32`). So "repoint place `m2__A` to the survivor `A`" can only be done by *string/Expr surgery* over every attribute: `recursively_substitute_vars!` (`compilers.jl:28-44`) blindly walks every `Expr` arg and replaces any `Symbol` equal to a map key (`compilers.jl:36-37`). This is **fragile by construction**: it cannot distinguish a genuine arc occurrence of place `A` from a coincidentally-equal symbol appearing inside a multiplicity coefficient, a rate subexpression, a parameter name, or an action body. A place name that collides with such a symbol is silently rewritten and the incidence is corrupted (this is the exact hazard recorded in ADR 0003, `0003-data-store.md:21`). Under ADR 0003 Phase 2, arcs are first-class `ArcSpec` rows carrying integer FKs `transition→T` and `place→S` (`0005-serialization-json-ir.md:39`). Identification then becomes: **repoint the `place` FK of every `ArcSpec` row from the deleted place index to the survivor's index, and drop the deleted `:S` row** — a structural, type-checked O(rows) operation with NO expression rewriting and NO possibility of collision-corruption. This is the single concrete structural reason the Phase-2 arc promotion is worth doing; ADR 0003 records this insight in its LENS-C rejection (`0003-data-store.md:50`), and the maintainer confirmed the promotion in the ADR 0003 status note. Expression rewriting survives only for the genuinely Expr-valued escape hatches (`@choose`/`@move`/`@structured`/expression-valued multiplicity), where the place reference is legitimately inside an `ExprNode` and is repointed by `Ref{place}` node identity, not by symbol-name matching.

### 7.5 The ADR 0004 live-guard

**J8 (Identification is forbidden on a stepping model).** `equalize!`'s `rem_parts!` (`equalize.jl:52`) is the one reindexing deletion that ADR 0004 forbids on a live model: collapsing places shifts `:S` indices, and compiled closures hard-code `state.u[i]` against a varmap frozen at construction (`compilers.jl:149`). The contract requires `equalize!` (and any composition path that reaches `rem_parts!`) to **refuse with an error when the model is live/stepping** (i.e. after a `ReactionNetworkProblem` has been constructed / `simulate` has begun). Composition and identification are authoring-time / pre-construction operations; the runtime mutation API (ADR 0004 append-only `add_place!`/`add_transition!`/`deactivate!`) is the *only* sanctioned way to change a live model, and it never reindexes. Under the promoted table (J7) the FK-repoint variant of identification is non-reindexing and could in principle be made live-safe, but the contract does NOT require that — identification stays an authoring-time operation.

### 7.6 Known composition bug to pin

**J9 (`@join` file branch calls an undefined `include_model`).** When a `@join` argument is a file-include macrocall (e.g. `@filename("model.jl")` style), `@join` lowers it to `:(include_model($str_inc))` (`joins.jl:226` and `joins.jl:228`), but `include_model` is **not defined anywhere** in the package (no `function include_model`, no method). Any `@join` that takes the file branch therefore throws `UndefVarError: include_model`. The test suite MUST pin this (it is currently uncovered — only the in-memory model-argument branch is exercised), and the rework MUST either implement `include_model` (load + parse a model file into a `ReactionNetworkSchema`) or remove the file branch. The contract treats the in-memory branch (`@join m1 m2 …` with already-constructed models) as the only currently-functional join entry point.

### 7.7 The composition contract (algebraic properties)

**C1 (Closure).** `union_acs!`/`@join` map model(s) to a model. The result is a valid input to a further join.

**C2 (Identity).** The empty model `ReactionNetworkSchema()` is a left identity for the fold (`@join` starts from it, `joins.jl:204`). Joining a model with the empty model returns a namespaced copy of that model.

**C3 (Commutativity holds only up to precedence and namespacing).** Join is NOT commutative in general. (a) Scalar attributes follow last-writer-wins (J5), so for two models that share an identified place with conflicting values, `@join m1 m2` and `@join m2 m1` differ in the surviving value. (b) The model identifier (`name`) namespacing is order-independent per model, but the *order of `:T`/`:S` part indices* in the result depends on argument order. The contract states: **join is commutative on STRUCTURE up to part reindexing, but its scalar-attribute result is order-dependent (last model wins).** Two joins that differ only in argument order are *isomorphic as models* iff no identified place has conflicting non-`missing` attributes.

**C4 (Associativity up to precedence).** Left-folding `union_acs!` is associative on structure: the disjoint union of parts and the by-name identification of places/params/meta do not depend on the grouping of the fold, because namespacing is per-model and identification is resolved against the global `eqs` blocks, not pairwise. Scalar precedence is the caveat — for a chain `m1, m2, m3` the surviving scalar of an identified place is `m3`'s (then `m2`'s, then `m1`'s) regardless of grouping, so associativity holds **on structure and on the last-writer-wins precedence order**, i.e. up to the same precedence qualifier as C3. (This associativity claim is contingent on J4 being fixed so that `:E`/`:obs` participate uniformly; with the current drop of events/observables, associativity is only guaranteed on {S, T, P, M}.)

**C5 (Identification is idempotent and order-free within a model).** `merge_eqs!` (`joins.jl:174-186`) transitively merges overlapping equation blocks before any rename, so declaring `A = B` and `B = C` collapses `{A, B, C}` to one survivor regardless of the order the equations are written, and re-identifying an already-collapsed set is a no-op. The contract pins identification as an **equivalence-class collapse**: the result depends only on the partition of places induced by the equations, not on their order or grouping. The surviving part is deterministic: its NAME is set by declaration (the `@alias` if present, else the first equation token, `equalize.jl:28`), while the surviving PART is the lowest-index member of the class (the others are removed by `rem_parts!(acs, :S, place_ixs[2:end])`, `equalize.jl:52`); its attributes are the first non-`missing` value across the class (`equalize.jl:46-51`) — a rule the contract pins so identification is reproducible.

**C6 (Meaning of "identification").** Identifying places `X` and `Y` asserts a *modeling* claim: they denote the same resource/quantity, so their token counts, costs, rewards, valuations, and modalities are one. Post-identification, every transition that referenced either now references the survivor; the two pools are summed into one `u` entry. This is a semantic commitment by the author, not an inferred equality — the engine never auto-identifies; sharing is always explicit via an equation block (J1) or `@equalize` (J6).

---

Grounding (all under `/Users/bima/ReactiveDynamics-review`): `src/operators/joins.jl:7-8` (precedence docstring), `:10-59` (`union_acs!`), `:64-82` (`prepend!`), `:87-97` (dead `prepend_obs`), `:100-128` (`normalize_name`), `:174-186` (`merge_eqs!`), `:201-238` (`@join`), `:226,228` (undefined `include_model`); `src/operators/equalize.jl:24-66` (`equalize!`), `:52` (forbidden `rem_parts!`), `:55-63` (Expr surgery); `src/compilers.jl:28-44` (`recursively_substitute_vars!`), `:69-79` (`escape_ref`), `:149` (frozen varmap); `src/ReactiveDynamics.jl:30-76` (schema: parts `:S,:T,:E,:obs,:P,:M`, `:trans` Expr holds arcs at `:53`); `src/interface/reaction_parser.jl:32` (`extract_arcs`); ADR refs `spec/adr/0003-data-store.md:21,39,50`, `spec/adr/0004-runtime-mutation.md`, `spec/adr/0005-serialization-json-ir.md:39`.

---

## 8. Serialization Schema (see ADR 0005)

This section is a CROSS-REFERENCE to [ADR 0005](adr/0005-serialization-json-ir.md) (single JSON serialization + typed `ExprNode` IR) and [ADR 0004](adr/0004-runtime-mutation.md) (append-only runtime mutation). It does NOT restate the `ExprNode` sum type, the `OP_WHITELIST`/`DIST_WHITELIST`/`REF_KINDS` whitelists, the document-shape field mapping, or the concrete JSON example — those are normative in ADR 0005 (§"ExprNode IR", §"Document shape", §"Concrete JSON example"). What this section pins is the *contract obligations* that serialization places on the rest of this document: the canonical artifact, the round-trip guarantee, the eval-free `validate` pass, the no-eval (RCE-closing) guarantee, the outputs-are-separate rule, and the mutation-patch form of the acquisition lever.

### 8.1 Canonical artifact — one JSON document

The canonical serialized form of a model is a SINGLE JSON document, `model.rdj.json` (ADR 0005 "Decision"), serving simultaneously as the on-disk model file and as the agentic-authoring artifact (Workstream F): the unit an LLM emits under structured output and self-validates before it is loaded. There is exactly one model format — the TOML/CSV/JLD2 zoo (`loadsave.jl`) is removed (ADR 0005 "Consequences"), and there is NO backward on-disk compatibility. The document is one JSON object with a `meta` object and the top-level arrays `params[]`, `places[]`, `transitions[]`, `arcs[]`, `observables[]`, `events[]`, mirroring the ADR-0003 typed IR one-to-one; in particular the `arcs[]` array IS the promoted first-class `ArcSpec` incidence table (ADR 0003 Phase 2, ADR 0005 "Document shape"), so the transition↔place relation (the net's arcs) that is today re-parsed per tick from the `trans` Expr (`extract_arcs`, `interface/reaction_parser.jl:32`) becomes structure on disk. The JSON Schema published for LLM structured output is AUTO-DERIVED from the `ModelSpec`/`ExprNode` types (the single-source-of-truth `const SCHEMA`) and round-trip-tested at build time, so schema and loader cannot drift (ADR 0005 "Julia mechanism").

### 8.2 Round-trip guarantee

**S1 (Semantic round-trip).** `to_json(from_json(j))` MUST be semantically equal to `j`: parsing a document and re-serializing it preserves every field's meaning (whitelisted `ExprNode` trees, arc FKs, modality axes, `meta` keywords). Equality is *semantic*, not byte-identical — key order, integer-vs-float JSON spelling (`1` vs `1.0` for an integer attribute, ADR 0005 "Open questions: Integer vs float"), and whitespace need not be preserved, but the constructed `ReactionNetworkProblem` MUST be identical. The `ExprNode` tree lowers via `to_expr` (ADR 0005 "from_json / to_json / validate") to exactly the place-name/param-name `Expr` that today's authoring macros already produce (e.g. the rate tree lowers to the same `:(rand(Poisson(max(state.dt * (0.3 * beta), 0))) * Preclinical)` that `expand_rate` emits at `create.jl:150-155`), so the round-trip is closed at the IR boundary, never through a Julia source string.

**S2 (Run determined by `(model.json, seed)`).** The model JSON is the complete reproducible INPUT: the pair `(model.rdj.json, seed)` fully determines a run, per §4 D1 (reproducibility) and D6 (seed at construction). `from_json(io; seed)::ReactionNetworkProblem` does parse → `validate` → `build_store` → construct, threading `seed` into the state RNG exactly as §4 D6 requires (`Xoshiro(seed)` when given, system-entropy + logged seed when not). Two runs from the same `(model.json, seed)` on the same platform MUST yield bit-identical trajectories (§4 D1). Nothing outside the model document and the seed may influence the trajectory; the `alloc_strategy` recorded in `meta` is deterministic and RNG-free (§4 D3, ADR 0002).

### 8.3 The `validate(spec) -> Vector{Diagnostic}` pass

**S3 (Eval-free validation).** `validate(spec)::Vector{Diagnostic}` is a PURE, eval-free static pass run before construction (ADR 0005 "from_json / to_json / validate"). It MUST NOT `eval`, `Meta.parse`, or execute any field. It is the LLM self-check before `from_json` and the construction-time gate for `validate` failures. It enforces, by walking the typed document:

1. **ExprNode walk** — every `Ref` name resolves to a declared place/param/observable; every `Call.op ∈ OP_WHITELIST`; every `Sample.dist ∈ DIST_WHITELIST`; arities are correct. Unknown names/ops/dists are diagnostics, not exceptions (ADR 0005 validate rule 1).
2. **Dangling arc FKs** — for every arc row, `r.transition` is a declared transition `id` and its place field names a declared place; a dangling FK on either is a diagnostic (ADR 0005 validate rule 2). This is the structural check the promoted incidence table makes possible (§8.1).
3. **§5 ranges / sign / integrality** — `transRate ≥ 0`, `transProbOfSuccess ∈ [0,1]`, `transCycleTime ≥ 0`, `transCapacity ≥ 0`, and integer-valued `capacity`/`init`/`multiplicity` for unstructured places, per the typed attribute domains of §5.3–5.4 and the integrality rule §5 A4. Because JSON has one number type, `validate` MUST check/coerce integrality rather than trust the parsed subtype (ADR 0005 "Open questions: Integer vs float"); emitters (LLMs especially) write `1.0` for integers.
4. **§1.4 illegal modality combos** — the orthogonal-axis modality (§1.1, `allocation ∈ {upfront,perstep}` × `return ∈ {consumed,conserved}` × `blocking ∈ {block,nonblock}`) is checked against the §1.4 rules: `blocking = nonblock` REQUIRES `return = consumed`; `allocation = perstep` REQUIRES a non-structured place AND `transCycleTime > 0`. This is the construction-time replacement for the legacy run-time `error` at `solvers.jl:461-465` and the unvalidated `Set{Symbol}` of `ReactiveDynamics.jl:144`.
5. **§5 A3 TVE policy** — an attribute marked "TVE? = no" in §5.3–5.4 (`placeInitVal`, `placeInitUncertainty`, `placeStructured`, `placeDefaultModality`) MUST be a literal `Const`, not a non-trivial `ExprNode` tree; authoring a time-varying expression for a frozen attribute is a diagnostic (ADR 0005 validate rule 5, enforcing §5 A3).

`validate` returns a `Vector{Diagnostic}` (a possibly-empty list of all violations), NOT a fail-fast exception, so an agentic author gets the complete diagnostic set in one pass to repair before reloading.

### 8.4 No-eval (RCE-closing) guarantee

**S4 (Eval-free by construction).** No field of a `model.rdj.json` document is EVER `Meta.parse`d or `eval`d on load. This closes the import-time remote-code-execution surface that exists today, deleting the eval sites at `loadsave.jl:65` (`eval(Meta.parseall(attrval))` for string-valued params) and `loadsave.jl:72` (`eval(Meta.parseall(row["body"]))` for the `registered` source block), and the eval-on-assignment `Base.convert` hooks at `ReactiveDynamics.jl:99` (`SampleableValues` ← `Meta.parse`), `ReactiveDynamics.jl:101` (`Set{Symbol}` ← `eval(Meta.parse)`), and `ReactiveDynamics.jl:102` (`FoldedObservable`). A model file becomes INERT DATA, not a Julia program: a malicious or malformed model can no longer execute code on load (ADR 0005 "Consequences"). The only place Julia source is produced is `to_expr`, which emits from a closed allow-list of interned symbols (`OP_WHITELIST`/`DIST_WHITELIST`/`REF_KINDS`) and hands the result to the unchanged `wrap_fun`/`compile_attrs` (`compilers.jl:148-180`) for a single construction-time compile; there is no `Expr`-head smuggling and no `apply`/`eval` op. This is the eval-free trade noted in §1, §4.5, and §5 A5: it intentionally REJECTS some currently expressible models (arbitrary `@register`d user bodies, raw call exprs), which the cutover must enumerate (ADR 0005 "Consequences", "Open questions: @register").

**S4 holds verbatim for registry-by-name callables (the general-code path).** A registered value-helper (ADR 0006 §C) and an `Invoke{fn, args}` statement callback ([ADR 0011](adr/0011-action-callbacks-and-general-code.md) §B, §12.3) both let a model reference ARBITRARY host logic, but neither weakens S4: the file carries only the NAME `fn` (a string) plus `ExprNode` args, and lowering produces `registry[fn](...)` — NO bytes from the file are ever `Meta.parse`d or `eval`d. The function body is host Julia, compiled by the user's own package; the trust boundary is the host program that populated the registry, not the file. An unregistered name is a `validate()` diagnostic, never code execution. The one honest caveat: an `Invoke` body is TRUSTED-BUT-UNVERIFIED — `validate` proves the name resolves and the arity matches, but cannot prove the body honors determinism (§4 D5) / append-only (ADR 0004) / tick-boundary (§12.5) — those become author obligations O1–O4 (ADR 0011 §B). Prefer the declarative action verbs (statically validatable); `Invoke` is the escape hatch.

### 8.5 Solutions and ledger are OUTPUTS, serialized separately

**S5 (Outputs are not in the model file).** The `sol` DataFrame (`solvers.jl:588-591`) and the `log` ledger (the `:allocation` and `:valuation_cost` rows in `evolve!`, `solvers.jl:304-312`; the `:valuation_reward` row in `finish!`, `solvers.jl:505`; the per-tick `:valuation` row, `solvers.jl:659-666`) are OUTPUTS of a run and MUST NOT be stored in `model.rdj.json`. They are persisted in a columnar store — Apache Arrow (Parquet, or CSV for human inspection) — replacing the Julia-version-fragile JLD2 blob at `loadsave.jl:205,221` (ADR 0005 "Solutions are a separate concern"). The recommended layout is `runs/<model_content_hash>/<seed>/{trajectory.arrow, ledger.arrow, run.json}`, where `run.json` is a tiny self-describing header `{model_hash, seed, rd_version, tspan, dt}` that matches a run back to its model and seed (S2). The ledger Arrow table is exactly what the rNPV/BD demo reads (the discounted `:valuation`/`:valuation_reward`/`:valuation_cost` streams, §"North-star tie-in" of ADR 0004). Keeping inputs and outputs in separate artifacts is what makes S2 meaningful: the model document is reproducible input, the Arrow tables are the reproduced output.

### 8.6 Mutation-patch form (the acquisition lever — see ADR 0004)

**S6 (Append-only JSON patch = serialized mutation).** A runtime mutation (ADR 0004) has a serialized form: an APPEND-ONLY JSON PATCH document — a `spec_delta` that may ADD places/transitions/arcs (and params/observables/rules) but MUST NEVER reorder or delete existing entries. This is the on-disk shape of the north-star acquisition lever: injecting a candidate program at a pipeline phase, perturbing operational params, or retiring a line (ADR 0004 "North-star tie-in"). The patch is applied by `apply_patch(state, delta)` at a TICK BOUNDARY (never mid-tick, ADR 0004 "Decision"), driving the live mutation API (`add_place!`/`add_param!`/`add_transition!`/`activate!`/`deactivate!`), so it inherits the operational-semantics (§3.3) and determinism (§4) guarantees and consumes no RNG (§4, ADR 0004 INV-4).

The append-only constraint is LOAD-BEARING, not stylistic: compiled attribute closures hard-code each place's position as `state.u[i]` via a varmap frozen at construction (`compilers.jl:148-153`), so a new place/transition MUST take the next free index and existing indices MUST never move (ADR 0004 INV-1, INV-2). A patch that reordered or deleted entries would invalidate every position-indexed closure (the sole mid-run reindexer, `rem_parts!` at `operators/equalize.jl:52`, is forbidden while live — ADR 0004 INV-2). Logical removal in a patch is soft-deactivation (`transActivated[i] = false`, gated at `state.jl:179`), under which in-flight instances still run to completion (`solvers.jl:406`, ADR 0004 INV-3 — stop starting new programs, never vaporize running ones). When a patch introduces a transition referencing a NEW place/param, `apply_patch` must add the place first, then `refresh_wrap_fun!` re-derives the varmap append-safely (pre-existing names keep the same `state.u[i]`, so already-compiled closures stay valid and are not recompiled, ADR 0004 "wrap_fun/varmap refresh mechanism"), then the new transition's expressions are compiled on-add (ADR 0004 "Compile-on-add"). For deterministic ensembles a patched-in transition SHOULD carry an explicit stable `name`; `name = missing` falls back to `gensym()` (`compilers.jl:177`), unique per call but not reproducible across runs (ADR 0004 "North-star tie-in"). The patch's `ExprNode`/action fields use the same closed whitelist and the same eval-free `validate` pass as the full document (§8.3, §8.4), so a mutation is as safe to load as an initial model.

---

## 9. Structured Tokens & Queries (see ADR 0006)

This section pins the structured/agentic-token subsystem — crucial for the Business-Development north-star, where a PROJECT is a structured entity (an agent with attributes, custom behavior, and history). Full treatment and source grounding are in [ADR 0006](adr/0006-structured-tokens.md); this section is the contract-level summary. A structured place is a KIND — a *colour set*, in Coloured-Petri-net terms (`placeStructured = true`) — whose INSTANCES are live AlgebraicAgents inner-agents under the `"structured"` container, distinct from a PLAIN place (a `Float64` `state.u` column). Kinds are defined a priori; instances are created, moved, and retired DURING simulation.

### 9.1 Instance lifecycle

Five phases on the existing token protocol (`interface/agents.jl:49-72`): (1) INSTANTIATE — `add_structured_token!(problem, agent)` = `entangle!(getagent(problem,"structured"), agent)` (`agents.jl:42-44`) at a tick boundary, or mid-RHS via `structured_rhs` (`solvers.jl:333-397`); (2) BIND as a transition resource — priority-ordered allocation in `evolve!` (`solvers.jl:202-215`), after which the token `isblocked` and drops out of re-allocation and the reflected count; (3) MOVE between places — `@move(:from,:to)` reuses the SAME token via `set_place!` (`solvers.jl:378`), preserving identity and `past_bonds` history; (4) UNBIND on completion (`solvers.jl:443-481`); (5) RETIRE — soft (`set_place!(:removed)`, keeps the token queryable for its audit trail) or hard (`disentangle!`, reclaims memory). APPEND-ONLY SAFETY: `entangle!`/`disentangle!` touch only the `inners` dict + AA bookkeeping, never `state.u` positions, so live instantiation cannot break the ADR 0004 frozen-varmap invariant (which governs plain-place columns only); the structured `state.u[i]` is a per-tick reflected count (`update_u_structured!`, `solvers.jl:630-639`), not authoritative storage.

### 9.2 Invariants

(1) exactly one place per token at a time; (2) a bound token is `isblocked`-excluded from re-allocation and the reflected count; (3) `length(active(prob, s)) == state.u[idx(s)]` after each `update_u_structured!`; (4) live instantiation/retirement never touches `state.u` positions; (5) **D4 total order** — any token iteration feeding allocation, an RNG draw, or a log entry MUST be totally ordered by `(place, name, uuid)` before use, never raw AA `Dict`-value order (`inners::Dict{String,…}`), so equal-`priority` binding is reproducible; (6) token `inners` keys derive from the AA `uuid`, not `randstring` names (which can collide and silently evict a live token); (7) token TYPES/CONSTRUCTORS are host Julia, referenced from JSON by name only, an unregistered name being a `validate()` diagnostic.

### 9.3 Query API

A deterministic surface (a new `src/interface/queries.jl`) layered on AA primitives, with `tokens(problem[, place])` as the single chokepoint that erases AA `Dict` order via `token_sortkey = (place, name, uuid)`; everything else composes it: `active`/`blocked` (on `isblocked`), `tokenattr` (guarded attribute read), and reductions `ntokens`/`nactive`/`sumattr`/`maxby`/`nbound`. AA's `@filter`/`f"…"`/raw `inners` iteration is exploratory-only and MUST NOT reach an RNG draw, a log entry, or an allocation order.

### 9.4 Querying tokens from a rate/action

A rate or action that reads the sibling token pool (e.g. "number of active phase-2 projects") uses a NEW closed ExprNode `TokenAgg{reducer ∈ {:count,:sum,:max}, place, field, active_only}`, lowered by `to_expr` to the §9.3 reductions — the ONLY sanctioned way a rate queries tokens (never an arbitrary registered function), keeping the authoring boundary eval-free. This is ratified as an ADR 0005 amendment (a closed `TOKEN_REDUCER_WHITELIST`, with `validate` arity/name checks and a round-trip test). `getobservable` on the structured container surfaces the same aggregates as named AA observables for the Workstream-E coupled-agent read path. The custom-function/kind REGISTRY (ADR 0006 §C) is the closed EXTENSION of the §6.5 `OP_WHITELIST`: `Call{op,…}` resolves builtins inline, registry names by lookup, and anything else is a `validate()` diagnostic — the eval-free replacement for `@register`.

### 9.5 Predicate selection & state advance of agentic tokens (see [ADR 0008](adr/0008-token-filtration.md))

This pins the maintainer's Q4 — "a transition may want to take a project in a given state / subject to a filtration expression, more complex than a previous place." Today a structured-LHS input arc selects tokens by KIND only: `filter(get_place(a)==kind && !isblocked(a))` then `sort!` by `priority`, take the granted count (`solvers.jl:202-211`, `:274-283`). You can ORDER by an attribute (custom `priority`, `agents.jl:67`) but cannot SELECT by one. The fix is ONE closed eval-free ExprNode, `TokenPredicate{kind, clauses}` (sibling of §9.4's `TokenAgg`), where `Clause = (field, op ∈ PRED_OP_WHITELIST{==,≠,<,≤,>,≥,∈}, value::ExprNode)` are AND-joined; the structured-LHS `ArcSpec` row (§6.5) gains an optional `predicate`, and binding generalizes to `filter(kind && !isblocked && matches(predicate,a))` BEFORE the unchanged `sort!`/take — a `filter` in front of the existing order, so it slots into ADR 0002 allocation and the §9.2 invariant-5 `(place, creation_index)` total order with no new machinery. DSL: `@select(Project, phase==:Phase2 && npv>θ)` on the LHS.

**Filtration measurability (determinism).** Honoring the word *filtration* literally: a `TokenPredicate.clauses[*].value` MUST be 𝓕ₜ-measurable — it may read token fields, params, observables, `TimeRef`, and `TokenAgg` aggregates, but NEVER a `Sample` (RNG) node and never the future, else the bind set would depend on draw order and §4 D1 is ill-defined. Token-reading subexpressions (in a predicate AND in a `TokenAgg`) observe the pre-`evolve!` reflected counts (`update_u_structured!`, `solvers.jl:649/653`) — this PINS the ADR 0006 §9.4 `TokenAgg` open observation-point question for both. `matches` is checked at the four bind sites; `validate` gains rule 7 (kind is `placeStructured`; ops/fields/value whitelisted; no `Sample`; `@select` on a plain place is a diagnostic). A `nothing`/empty predicate is exactly today's kind-only bind — strictly backward-compatible.

**Phase is an attribute — the canonical model (maintainer ruling).** A token's lifecycle state (pipeline phase, status) is a TOKEN ATTRIBUTE, not a place KIND: ONE `Project` KIND carries a `phase::Symbol` field; "a project in phase 2" is `@select(Project, phase==:Phase2)`. This unifies with KIND selection because `get_place(a)==kind` is just the degenerate clause `(place,==,kind)`, so the per-KIND-per-phase style (the original ADR 0006 north-star sketch) is the legal degenerate, retained only as `@move` back-compat and not the recommended pattern. Rationale: no KIND explosion; continuous/cross-field predicates (`npv>θ`, `cost_to_date<cap`) become expressible; uniform `past_bonds` history on one identity.

**Advancing state = a field write (`SetField`/`@advance`).** Because phase is a field, advancing a project WRITES that field on the RHS. The ADR 0005 action family extends to `{SetMarking, SetParams, SetField, Log, Seq}`: `SetField{field, value}` writes a field of the firing instance's bound token(s) in `finish!`, reusing the existing `@move` write/unbind machinery (`solvers.jl:378`). `@advance(phase, :Phase3) ⇒ SetField{:phase, Const(:Phase3)}`; `@move(:from,:to) ⇒ SetField{:place, …}` (so `@move` is now sugar for the degenerate place-field write, and its redundant re-`entangle!` `solvers.jl:376` is dropped — a field write never changes the uuid container key); a general `@advance(npv_estimate, npv_estimate*uplift)` revalues. A new `Field{name}` ExprNode leaf reads the bound token's own current field. A `SetField.value` MAY contain a `Sample` (a write may legitimately draw, consuming the seeded stream deterministically, §4 D5) — UNLIKE a predicate `value`, which may not. Note the bind/read/write split: `TokenPredicate` consumes/holds (CONTRACT §1 modality applies, `perstep` still illegal for structured places §1.4); `TokenAgg` reads non-consumingly; `SetField` writes a bound token's field on completion — all share the field/ExprNode vocabulary, never alias. Per-phase counts (`#Phase2`) come from `TokenAgg`/`nactive(Project, phase==…)` rather than one reflected `state.u` column per phase, so the §9.2-inv-3 `length(active)==state.u` assertion becomes per-predicate (an ADR 0008 open question against `update_u_structured!`, `solvers.jl:630-639`).

---

## 10. Interface & Initial-State Contract (see [ADR 0007](adr/0007-interface-and-initial-state.md))

This answers Q3 — "do we have a contract for the interface methods (network definition, initial values / initial state, simulation, joins)" — and the follow-up that the interface must support instantiating structured agents as a list of structures and DUMPING the system state. §1–§9 pin the network's SEMANTICS; this section pins the callable INTERFACE: which methods exist, in which order they may be called, and what each requires/guarantees. It closes two holes: structured-token initial state is not declarative today (it is built by imperative host code, breaking §8.2 S2), and there is no state dump/restore.

### 10.1 The three-phase lifecycle

A model passes through exactly three phases; each method is legal in a defined subset, and the phase boundary is an observable state of the object.

| Phase | Object | Mutability | Legal operations | Forbidden |
|---|---|---|---|---|
| **Authoring** | `ModelSpec`/`ReactionNetworkSchema` | freely mutable; reindex OK | `@ReactionNetworkSchema`; `from_json`; `@push`; `@join`/`union_acs!` (§7); `@equalize`/`equalize!` (§7.4, the `rem_parts!` reindexer); `refine`/`abstract`/`@compose` (§11); `validate(spec)` (§8.3) | stepping; live-token ops |
| **Construction** | `ReactionNetworkProblem(spec; seed, u0, registry, population)` | one-shot transition | freezes the varmap (`compilers.jl:149`), seeds the RNG (§4 D6), compiles closures, INSTANTIATES the initial marking (§10.3), entangles the `"structured"` container (`solvers.jl:612`) with that population | — (atomic) |
| **Live / stepping** | constructed `ReactionNetworkProblem` | APPEND-ONLY (ADR 0004) | `simulate`/`step!`; `add_place!`/`add_transition!`/`add_param!`/`activate!`/`deactivate!`; `add_structured_token!`/`disentangle!`; `apply_patch` (§8.6); query API (§9.3); `dump_state`/`reinit!` | any `rem_parts!` reindex (§7.5 J8); reordering/deleting object rows |

The §7.5/J8 live-guard GENERALIZES to a phase guard: any reindexing op (sole reindexer `rem_parts!`, `equalize.jl:52`) MUST refuse once the object is Live; construction is the arming point. This table IS the Q3 answer — definition (authoring), initial state (construction §10.3), simulation (live §10.4), joins (authoring §7) — with their legal ordering made explicit.

### 10.2 The constructor and canonical signatures

`ReactionNetworkProblem(spec; seed::Union{Integer,Nothing}=nothing, u0=Dict(), registry=Dict(), population=spec.population)::ReactionNetworkProblem` — adds `seed` (§4 D6), `registry` (ADR 0006 C), and `population` (§10.3) to the current head (`solvers.jl:536`). The canonical surface: `validate(spec; registry_names)::Vector{Diagnostic}` (authoring); `from_json(io; registry)`/`to_json(spec)` (authoring, ADR 0005); `@join`/`@equalize`/`refine`/`@compose` (authoring, reindexers refuse when live); the constructor (construction); `simulate(problem, max_t=Inf)::problem` (live); the append-only mutation + query APIs (live); `dump_state`/`restore`/`reinit!` (live).

### 10.3 Declarative initial marking — `population[]`

For a PLAIN place, "initial value" is the literal `placeInitVal` (§5.4), serialized in `model.rdj.json`. For a STRUCTURED place it is NOT a value but a POPULATION of token instances, and today that population is built by imperative host code (`add_structured_token!` in a loop, `agents.jl:42-44`) — so it is NOT in the document and §8.2 S2 ("(model.json, seed) determines the run") is VIOLATED for structured models. The fix is a declarative, serializable top-level `population[]` array (INPUT, not output §8.5), the structured analogue of `placeInitVal`:

- **Two authoring forms, one semantics:** (a) `{place, kind, count, attributes}` — `count` instances of `kind`, each with `attributes` (an `ExprNode` per field) drawn through `state.rng`; (b) `{place, kind, instances:[{attributes…}, …]}` — an explicit LIST of structs (the maintainer's "list of structures"; also the host-ergonomic `population = [ProjectToken(:Phase2; npv=3.1), …]` form).
- **Instantiation contract (at construction, before t=0, before the first `_step!`):** iterate `population[]` in declared order; for each instance resolve `kind` against the registry → host constructor (ADR 0006 B/C), evaluate `attributes` `ExprNode`s through `state.rng` (§4 D5), `entangle!` into the `"structured"` container; assign the per-place CREATION INDEX k=1..count (the §9.2 invariant-5 tie-break, the 2026-06-20 ruling), which is part of the seeded run state; then `update_u_structured!` reflects live counts into the structured `state.u` columns (so initial `u` of a structured place = its initial active population). Eval-free: references kinds/fields by name, never Julia source (ADR 0006 B). `validate` rule 6: every population entry's place is `placeStructured`, every `kind` is registered, every `attributes` field is a whitelisted `ExprNode`.

### 10.4 Simulation interface and `reinit!`

`simulate(problem, max_t=Inf)` is the AlgebraicAgents driver (`AA interface.jl:160`): it takes a max-TIME, not a step count, and RETURNS the agent — a pinned gotcha. `_step!` advances one tick (§3.3, `solvers.jl:642`); `_projected_to` reports completion at `state.t > tspan[2]` (`:675`). `reinit!` (`solvers.jl:619-628`) MUST be COMPLETED to: restore the RNG to the seed-implied initial state (§4 D7 — currently missing); tear down the live token pool and rebuild the §10.3 marking (currently end-state tokens survive into the next run); reset every per-place creation counter. After this, `init → step* → reinit! → step*` reproduces the first trajectory for BOTH plain and structured models.

### 10.5 State dump / restore (checkpoint)

`dump_state(problem)::StateDump` serializes a live run at a TICK BOUNDARY; `restore(spec, dump; registry)::ReactionNetworkProblem` reconstructs it. **The dump schema is the §10.3 initial marking PLUS the dynamic run state:** `t`/`tick`, the RNG state (§4 D7), per-place creation counters, plain-place `u`, the full token population (each token's `kind` name + CURRENT field values + `uuid` + `creation_index` + `bound_transition`-by-id + `past_bonds`), and the in-flight `ongoing` transitions (by recipe index `i`, birth `t`, `q`, progress `state`, bound tokens by uuid). The ledger stays in Arrow by reference (§8.5). **Restore = construction with overlays:** run the §10.3 marking machinery on `dump.population` (reconstructing types via the registry — field VALUES + kind NAME, never source, so it inherits the §8.4 no-eval guarantee), then overlay `t`/`u`/RNG/counters and re-link `ongoing` by uuid; structured `state.u` columns are RE-DERIVED (never restored directly), keeping §9.2 invariant 3 true. **Key identity: a checkpoint with `t==tspan[1] && ongoing==[]` IS an initial marking** — `restore(spec, zero_tick_dump)` ≡ `construct(spec; population=dump.population)`. This is why the initial marking (§10.3) and the dump are one schema; it makes halt/resume and "fork a run at a tick to A/B a lever" first-class (the acquisition counterfactual: `dump_state` at the lever tick → `restore` twice → apply the lever to one copy → diff under one seed, §4 D5).

### 10.6 Interface invariants

(1) **Marking-determinism** — given `(spec, seed)`, the initial marking (counts, sampled attributes, creation indices) is reproducible (§4 D1/D5). (2) **Reflected-count consistency at t=0** — `length(active(prob,s)) == state.u[idx(s)]` after the marking is built. (3) **Dump round-trip** — `restore(spec, dump_state(p))` steps identically to `p` from that point (the resume analogue of S2). (4) **Eval-free state I/O** — neither `population[]` nor a `StateDump` is ever `Meta.parse`d/`eval`d; both reference kinds by name through the host registry (§8.4 S4). (5) **Phase legality** — a reindexing op on a Live object errors; append-only ops are safe in any phase. (6) **`reinit!` completeness** — `init → step* → reinit! → step*` reproduces the first trajectory for plain AND structured models (closes §4 D7 for structured runs).

---

## 11. Refinement & Open-Port Composition (see [ADR 0009](adr/0009-refinement-and-composition.md))

This adds the REFINEMENT axis the contract was otherwise silent on, meeting the maintainer's framing requirement — compact/expressive definition, compositionality, and various levels of granularity with more refined dynamics possibly substituted. It is entirely AUTHORING-time and additive: it produces a plain `ModelSpec` that constructs/serializes/simulates exactly as a hand-written flat model. The enabling insight is §7.4/J7: once arcs are first-class `ArcSpec` FK rows (ADR 0003 Phase 2), place identification is integer-FK repointing, not Expr surgery — refinement reuses that exact splice, so the headline feature costs almost no new mechanism. All of §11 is gated on the ADR 0003 Phase-2 arc promotion.

**§11 is IMPLEMENTED** (ADR 0003 Phase 2 having landed): `src/operators/refine.jl` provides `set_port_role!`/`@port` (§11.1 the `placeRole` closed tag on Place), `refine!`/`refine`/`abstract_transitions` (§11.2 the FK-splice — it delegates namespacing + port identification to `union_acs!`'s equation-alias mechanism, then `rem_parts!`-drops the coarse transition and re-derives the FK-exact `ArcSpec` table), `refinement_diagnostics` (§11.3 advisory port-balance + linear-chain cycletime/PoS aggregate checks), and `@pipeline`/`@process`/`compose`/`@compose` (§11.4). Exercised by `test/semantic/refinement_composition.jl`. `@compose` merges `:E`/`:obs` (J4) via the WS-3 `union_acs!` fix and never reaches `include_model` (J9). Entity-level refinement (§F) remains deferred.

### 11.1 Open ports

A **port** is a boundary place through which a fragment connects to its environment — a thin closed-tag annotation on the Place record (§6.4, NOT a new table): `role ∈ {private, input, output, shared}` (default `private`). `private` = internal, auto-namespaced on compose (`m__X`); `input`/`output` = an OPEN port the fragment expects connected (directionality advisory/validation-only); `shared` = a global place identified by bare name without namespacing (the existing `@catchall` semantics `joins.jl:117` made first-class). A **transition's boundary** is exactly its input/output `ArcSpec` rows — no new structure; the incidence table already IS the boundary.

### 11.2 Refinement = boundary-matched FK-splice

`refine(spec, transition, submodel; ports)` replaces a coarse transition `T` with a finer sub-model `S` plug-compatibly at its boundary, in four authoring-time moves: (1) namespace `S`'s `private` places, leave `input`/`output`/`shared` ports un-prefixed; (2) identify `S`'s open ports with `T`'s boundary places by the §7.4/J7 FK-repoint (repoint `ArcSpec.place` FKs, drop the port `:S` row — NO `recursively_substitute_vars!`, so no collision-corruption); (3) append `S`'s transitions + remaining places/params/obs/EVENTS as new rows (append-only §6.2; also fixes §7/J4 by merging `:E`/`:obs`); (4) remove the coarse `T` (a reindex → authoring-time only, §7.5/J8). Because boundary places keep their indices/names/attributes, coarse and refined models are PLUG-COMPATIBLE: every other transition over the boundary is untouched. `abstract(spec, transitions, into; ports)` is the inverse collapse.

### 11.3 Boundary consistency — advisory, not equivalence

Refinement does NOT claim behavioral equivalence (that would need a bisimulation the framework can't check). `validate` emits ADVISORY diagnostics where computable: `coarse.cycletime ≈ Σ critical-path cycletimes`, `coarse.pos ≈ Π sub-PoS`, `coarse.cost ≈ Σ sub-costs`, and port-balance (every `input` consumed, every `output` produced). These are overridable warnings — the refinement may legitimately change dynamics (the point of zooming in); they make the granularity ladder auditable without overclaiming. The engine never auto-refines/auto-abstracts.

### 11.4 Compact authoring and port-connected composition

**`@pipeline`** sugar for the dominant chain shape expands to N `flow`-genesis transitions (§2.8 routing) with the `@select`/`@advance(phase,…)` idiom (§9.5 phase-as-attribute is canonical: one `Project` kind, advance writes the `phase` field via `SetField`), each carrying per-edge `(ct, pos, res)` — the BD pipeline in one block. **`@process`** = a named parameterized `ModelSpec` fragment with declared ports, instantiated multiple times (eval-free param substitution, not code-gen) and composed. **`@compose f1 f2 …`** is `@join` (§7) PLUS automatic port matching: `output` ports identified with same-named `input` ports by the §7.4/J7 FK-repoint, `private` namespaced, `shared` by bare name. It is the explicit-boundary form of `@join` (which remains the manual path) and CLOSES the §7/J4 (`:E`/`:obs` dropped) and J9 (undefined `include_model`) bugs en route. Entity-level refinement (a structured token hosting its OWN sub-network — a second granularity locus via AlgebraicAgents) is DEFERRED to a future ADR; v1 is process-structural only.

### 11.5 Refinement invariants

(1) **Plug-compatibility** — after `refine`, every transition outside `{T} ∪ S` is structurally unchanged (same boundary-place indices/rows). (2) **FK-exactness** — port identification repoints integer FKs (§7.4/J7), never `recursively_substitute_vars!`. (3) **Authoring-only** — `refine`/`abstract`/`@compose`/`@equalize` reindex and so are FORBIDDEN on a Live object (§7.5/J8, §10.1); the append-only mutation API (ADR 0004) is the only runtime change path. (4) **Closure** — they map `ModelSpec`(s) to a `ModelSpec`, valid input to a further refine/compose/construct (§7 C1 extended to the vertical axis). (5) **Round-trip** — a refined/composed model serializes/re-loads as a flat `model.rdj.json` identical to a hand-authored equivalent (§8.2 S1); refinement leaves no runtime trace. (6) **Advisory** — the §11.3 diagnostics are warnings, not equivalence proofs.

---

## 12. Rules, Triggers & Conditional Transitions (see [ADR 0010](adr/0010-rules-and-conditional-transitions.md))

This pins the ENDOGENOUS DECISION CHANNEL — the ability to express, *inside the model*, decisions that fire on conditions (scheduled OR state-contingent), which a business-process model needs and the Phase-0 contract lacked. A BD analysis is not a fixed schedule; it is decisions taken in response to conditions ("start Phase-3 *if* cash covers its cost"; "acquire *when* the pipeline looks under-valued"; "retire the discovery line once the portfolio is late-stage-heavy"). §8.6 gave a *host*-applied patch (scheduled only); the no-op event channel (§3.4 Inv 7, `solvers.jl:323` fetches but never evaluates `:eventAction`) was the latent in-model channel. This section makes it live and adds a stateless transition guard. The unifying insight: **a rule/callback and a conditional transition are the same primitive — a guard `ExprNode` evaluated per tick, driving a closed action set** — differing only in what the guard gates (a Rule gates an action block; a conditional transition gates its own firing). Both are eval-free (the boolean/comparison ops `> < >= <= == && || !` are already in `OP_WHITELIST`, ADR 0005), deterministic, and serializable.

### 12.1 The Rule object (the repaired Event, §6.8)

The `:E` object is a **Rule**: `{ id, guard::ExprNode, action::ActionStmt, fire_mode ∈ {every_tick, once}, enabled::Bool }`. The guard is a standard `ExprNode` (§1/ADR 0005) evaluated through the RNG-threaded `context_eval(state, nothing, …)` (§4 D5): a `Bool` guard fires the action once, a numeric guard `v` fires it `rand(state.rng, Poisson(v))` times (the existing `:321` multiplicity, now seeded). The repair is one line — `event_action!`'s bare fetch at `solvers.jl:323` becomes `context_eval(state, nothing, state.wrap_fun(action_i))`. Rules iterate `parts(state, :E)` in index order (§4 D4); multiple rules firing on one tick run in that order and last-write-wins on a shared place/param (`validate` MAY warn on write-write conflict).

### 12.2 Conditional transitions — the stateless `guard`

A Transition (§6.6) carries `guard::ExprNode` (default `Const(true)`, resolving to `Bool`). It is evaluated each tick in `sample_transitions!` and AND-ed with the existing latching gate (`state.jl:179`): `fires_this_tick = transActivated[i] && eval(guard_i)`. The two are orthogonal and compose: `guard` is **stateless** (re-evaluated every tick from current state), `transActivated` is **latching** (flipped only by `activate!`/`deactivate!` or a Rule action — ADR 0004). A guarded-off transition contributes NO genesis proposal that tick (skipped before `evolve!`), so it never competes for resources (clean ADR-0002 interaction). A non-`Bool` guard is a `validate()` diagnostic; an RNG-consuming guard (`Sample` node) is legal but `validate` SHOULD warn. **Hysteresis guidance:** a condition that must persist once tripped is a `once` Rule or a Rule that sets a latch param / `Deactivate`s a line — NOT a bare stateless `guard`, which un-fires the moment the condition lapses.

### 12.3 The closed action set

The ADR-0005 action set `{SetMarking, SetParams, Log, Seq}` is extended to the full endogenous action *type family* `{SetMarking, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}` (`ACTION_VERBS = (:set_marking, :set_params, :set_field, :set_tokens, :add_token, :activate, :deactivate, :invoke, :log, :seq)`; canonical definition in [ADR 0010](adr/0010-rules-and-conditional-transitions.md) §C, extended by [ADR 0011](adr/0011-action-callbacks-and-general-code.md)). `SetMarking` carries a `mode ∈ {set, inc}` (the `inc` form is the cash/capital lever). `SetField{field, value}` (§9.5, ADR 0008 §D) writes a field of the firing instance's bound token(s) — the per-token counterpart of `SetMarking`; it is how a program advances phase (`@advance(phase, :Phase3) ⇒ SetField{:phase, …}`), and `@move(:from,:to)` is its degenerate `SetField{:place, …}` sugar. `SetTokens{predicate, assigns}` ([ADR 0011](adr/0011-action-callbacks-and-general-code.md) §A) is the POPULATION generalization of `SetField`: it writes `assigns` (`field => value` pairs, `value` read per-token via a `Field` leaf and MAY draw) over every token matching a §9.5 `TokenPredicate`, iterated in the §9.2 `(place, creation_index)` total order — "write down all Phase-2 oncology `pos_remaining` by 10%" in one action. `AddToken{kind, fields}` lowers to `add_structured_token!(state, registry[kind](…))` — `entangle!`, append-only/index-safe (§9.1) — the in-model acquisition (CREATE a token). `Activate`/`Deactivate{transition}` flip `transActivated` (ADR 0004 soft-deactivate; in-flight instances still finish, `solvers.jl:406`). `Invoke{fn, args}` ([ADR 0011](adr/0011-action-callbacks-and-general-code.md) §B) is the GENERAL-CODE escape hatch: it lowers to `registry[fn](state, transition, <args>)` — the ADR-0006 §C registry-by-name boundary in statement position — calling a host-compiled function that may run arbitrary queries/mutations/logic, subject to author obligations O1–O4 (determinism via `state.rng`, append-only mutation, no clock re-entry, purity). **Context restriction:** `SetField` writes the firing instance's bound token, so it is legal ONLY in a transition **post-action**; a standalone Rule (fired at step 10 with no transition instance) has no bound token, so a **Rule action** uses `{SetMarking, SetParams, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}` and `SetField` in a Rule is a `validate()` diagnostic (`SetTokens` carries its own predicate so it IS legal in a Rule — the finding-I fix). All declarative verbs lower eval-free into the closed allow-list extended by the ADR-0006 registry key set; an unregistered `AddToken.kind`/`Invoke.fn`, an arity mismatch, an unknown `Activate`/`Deactivate` target, an unknown `SetField`/`SetTokens` field/kind, an out-of-set verb, a `SetField` outside a transition post-action, or a `Ref` to an undeclared name is a `validate()` diagnostic, never an eval (§8.4 S4 preserved — `Invoke` carries only a NAME in the file, its body is host Julia). `Invoke` is the one TRUSTED-BUT-UNVERIFIED tier: prefer the declarative verbs (statically validatable, LLM-authorable); reach for `Invoke` only when logic exceeds them.

### 12.4 `fire_mode` and the `once` latch

`every_tick` (default) re-evaluates the guard every tick (a periodic financing rule is `every_tick` with a `periodic(state, p)` guard, `state.jl:239`). `once` fires at most once per run: on first firing the engine sets `enabled = false` and skips the rule thereafter. `enabled` is seeded RUN STATE — **`reinit!` MUST reset every `once` rule's `enabled` to `true`** (§4 D7, joining the RNG state, per-place creation counters §9.2, and the §10.4/§10.5 dump/restore run-state items), so `init → step* → reinit! → step*` reproduces the run including the lever. A `once` rule + a `periodic(state,0.0)` guard is the canonical scheduled one-shot at t0.

### 12.5 Step placement (determinism-critical)

Pinned in the normative `_step!` order (§3.3): **transition guards** are evaluated inside `sample_transitions!` (step 5, against this tick's freshly-sampled attributes); **Rules** fire at step 10 (the former event slot, `solvers.jl:657`) — AFTER `finish!` and the structured-count sync (steps 8-9), BEFORE the ledger row (step 11) and the clock advance (step 12). Consequence: a Rule's `SetMarking`/`SetParams`/`AddToken` hits the SAME tick's valuation row and is visible to the NEXT tick's guards and genesis; a token injected by `AddToken` is reflected by the adjacent `update_u_structured!` and is bindable from the next tick. This placement is normative; moving it changes results.

### 12.6 Determinism obligations

Guards and rules consume RNG only via explicit `Sample` nodes through `state.rng` (§4 D5); fire points and `:E`/transition iteration order are fixed (§4 D4); `once` latches reset on `reinit!` (§4 D7); `(model, seed)` fully determines the trajectory including all rule effects (§4 D1). A guard/`TokenAgg` that reads sibling token counts observes the §9.5 pre-`evolve!` measurable point and MUST be 𝓕ₜ-measurable (no `Sample`), consistent with the §12.2 "guards SHOULD NOT consume RNG" guidance — so token-aware guards are deterministic. The acquisition counterfactual is the ensemble-averaged diff of two same-spec runs differing only by the scenario params that arm a rule (single RNG stream retained; per-entity substreams are §4 future work — see [MVP_BD_DEMO.md](../demo/bd_acquisition/MVP_BD_DEMO.md) §4.1/finding A).

---

## 13. AlgebraicAgents Integration & External Coupling (see [ADR 0012](adr/0012-algebraicagents-integration.md))

A `ReactionNetworkProblem` already IS an AlgebraicAgents (AA) `@aagent` implementing `_step!`/`_reinit!`/`_projected_to` (`solvers.jl:642/619/675`), so `entangle!(parent, rd)` makes a reactive network ONE node in a larger AA hierarchy and `simulate(root)` drives it. This section pins how RD interoperates as a hierarchy node in both directions; the foreign-agent topology lives in host-side `add_wire!` calls, never in the inert RD document (§8.4 S4 preserved).

### 13.1 Clock coordination

AA coordinates HETEROGENEOUS clocks via its least-projected-time gate: `step!` projects the whole hierarchy to the global minimum time and calls an agent's `_step!` only when `_projected_to(a) == t` (`AA interface.jl:191-197`). RD's single clock (§2.1, `state.t += dt`) is unchanged — it is now ONE node's local clock, interleaved with siblings at finer/coarser `dt` with no new code. A coupled net MUST be driven via `simulate(root)`, not `simulate(rd)` (the §10.4 standalone path still works for an un-coupled net); RD's `dt` is the coupling granularity (§2.6).

### 13.2 Outbound — RD as a readable node

RD implements AA's read surface so the hierarchy (and wires) can read it: `observables(rd)` lists exported places, named observables (§9.4), and explicitly-declared `TokenAgg` aggregates (§9.4, e.g. `nactive(Project, phase==:Phase2)`); `getobservable(rd, name|i)` returns the current value (place marking = `state.u[idx]`; observable = `.sampled`; aggregate = the §9.3 ordered reduction); `_getparameters`/`_setparameters!` expose/patch `state.p` (param-only, ADR 0004 index-safe). This subsumes the ADR 0006 §E `getobservable`-on-the-`"structured"`-container TODO and RESOLVES MVP finding G (per-tick token KPIs).

### 13.3 Inbound — `inputs[]` ports + the `ExternalRef` leaf

RD dynamics read external state declaratively (maintainer ruling: wires + `ExternalRef`, NOT a registry walk). The model gains a top-level `inputs[]` array of named read ports (with optional `default`) — the OBSERVABLE-level analogue of the §11.1 open PLACE ports. A new closed `ExprNode` leaf `ExternalRef{port}` (sibling of `Ref{obs}`) reads a port; it lowers to a read of the per-tick external-input buffer and is usable in any rate, guard (§12), or `TokenPredicate.value` (§9.5). Host wiring `add_wire!(sys; from=macro, to=rd, from_var_name="rate", to_var_name="ext_rate")` (`AA wires.jl:54-73`) feeds a foreign agent's observable into a port; the RD document declares only the ports it consumes, never the foreign agents. `validate` rule 8: every `ExternalRef.port` is a declared `inputs[]` port.

### 13.4 The latch point — coupling determinism

External inputs are read ONCE per tick and buffered at a PINNED point: RD implements AA's `_prestep!(rd, t)` to call `retrieve_input_vars(rd)` (`AA wires.jl:26-35`, = `getobservable` over incoming wires) into `state.external_inputs`. This is load-bearing for §4 D4: AA steps siblings in raw `Dict` order (`AA interface.jl:185`), so a LIVE mid-`_step!` cross-agent read would be sibling-order-dependent; reading at `_prestep!` — `prewalk`'d over the whole hierarchy BEFORE any `_step!` (`AA interface.jl:184`) — means RD reads each source's PREVIOUS-boundary projection, order-independent. This is explicit (Jacobi) co-simulation: a one-tick coupling lag, no algebraic loop, reproducible under `(hierarchy, seed)` (§4 D1). Every `ExternalRef` read within a tick returns the same buffered value (the §9.5 same-snapshot discipline). `ExternalRef` is 𝓕ₜ-measurable (reads past/current external state, never the future, consumes no RNG). `state.external_inputs` is transient: `reinit!` clears it (§10.4); `dump_state` need not persist it (§10.5) — it is recomputed on the next `_prestep!`.

### 13.5 Rejected/escape: registry-walk inbound (`Invoke`)

A host-registered function (§12.3 `Invoke`, ADR 0011) that itself walks `getagent(state,"../sibling")` + `getobservable` is REJECTED as the default (imperative, not validatable, re-incurs the ADR 0011 O1–O4 obligations) but RETAINED as the escape hatch for reads exceeding a single wired port (a cross-hierarchy reduction, a glob/regex agent query). `validate` steers authors to `ExternalRef`+wires.

### 13.6 Invariants

(1) RD resolves `getobservable`/`observables` for every exported name (no AA `@error` fall-through). (2) Every `ExternalRef` read in a tick returns the value latched at `_prestep!` (the source's previous-boundary projection), identical across read sites, independent of sibling step order (§4 D4). (3) Explicit coupling — one-tick lag, no within-tick fixed point; `(hierarchy, seed)` determines the coupled trajectory. (4) Eval-free — ports referenced by name, foreign topology in host `add_wire!`, no file bytes parsed/eval'd (§8.4 S4). (5) Read-only inbound; `_setparameters!` writes params only (no reindex/structural change). (6) A coupled RD net is clock-subordinate to `simulate(root)`.

---

## 14. Analysis & Observability (see [ADR 0013](adr/0013-analysis-and-observability.md))

Running a simulation already yields an inspectable result — the plain-place trajectory `state.sol::DataFrame` (`state.jl:85`), the heterogeneous per-tick event log `state.log` (`state.jl:80`), and the per-program ledger `program_ledger(state) → DataFrame` (`ledger.jl:201`). This section pins the missing analysis layer: a TIME-INDEXED history of agentic-token state, an ENSEMBLE runner with cross-run statistics, and a RESULTS export bundle. All three are additive read-only mechanisms over a finished (or in-progress) run; none touches the §1–§9 dynamics. This is the analysis counterpart to §8/§10.5 (which serialize the MODEL and the resumable STATE) — §14 serializes and aggregates the OUTPUTS.

### 14.1 Per-token trajectory log

Structured tokens are passive agents (`_step!(::AbstractStructuredToken) = nothing`, `agents.jl:135`) — they have no self-stepping point, so any per-token history is ORCHESTRATOR-driven. A per-kind hook `log_token_fields(tok)::NamedTuple` (default empty) declares which fields a kind records; the orchestrator appends `(t, token_name, fields)` for each opted-in active token to a columnar `state.token_trajectory` store, iterated in `token_sortkey` (place, creation_index, uuid) order (`predicates.jl:105`) for byte-reproducibility (§4 D4). The new tick step `push_token_trajectory_row!(state)` sits immediately after `push_program_ledger_row!` (`solvers.jl:1037`) and before the clock advance (`:1039`) — the same seam, observation point, and determinism guarantees as the per-program ledger row it generalizes. Read via `token_trajectory(state[, name|pred]) → DataFrame`; a `TokenPredicate` (§9.5) scopes the rows. "Typical" trajectory is made precise as two helpers over this store: `representative_token` (the medoid path — a single typical life) and `trajectory_envelope` (aligned median + IQR band — typical ± spread).

### 14.2 Ensemble runner

`ensemble(build; nseed, root_seed = 2026, max_t, parallel = false, mode = :rebuild) -> EnsembleProblem` runs `nseed` independent members, member `k` seeded `hash((root_seed, k))` (§4 D8); each owns its `state.rng` (`state.jl:91`), so results are order- and parallelism-independent (§4 D9). `EnsembleProblem <: AbstractAlgebraicAgent` is a container whose `inners` are the members and which implements the §13.2 read surface (`observables`/`getobservable` return cross-run reductions), so an ensemble is itself a readable/drawable AA hierarchy node — built ON AA's primitives, not by modifying AA. Statistics: `summarize(ens, run -> metric)` (mean/sem/quantiles) and `treatment_effect(ens_a, ens_b; metric)` (the unpaired Δ with `se = sqrt(var_a/n_a + var_b/n_b)`) — the A/B lever comparison. Two member-production modes, recorded on the result: `mode = :rebuild` (mode a, default) REBUILDS a fresh problem per seed (robust via the §10.3 declarative `population[]`); `mode = :reinit` (mode b) builds ONE member then REINIT-RESEEDS and re-simulates that same problem for each subsequent seed (`reinit!(m; seed)`, the §10.4/ADR 0007 §D `_reinit!` reseed path), reusing its compiled closures + allocated store — the cheaper Monte-Carlo path. After each mode-(b) run a faithful `deepcopy` snapshot is retained, so `ens.members` holds `nseed` independent `ReactionNetworkProblem`s and the read surface + any `metric(member)` are contract-identical to mode (a). **Mode (b) is a legal substitute for mode (a) ONLY when members are structurally HOMOGENEOUS** (same net + `population[]` schema, differing only in the stochastic stream and seed-sampled initial attributes); the reseed installs the new seed's stream BEFORE the t=0 marking is re-sampled, so a reinit-reseeded member is member-for-member identical to a fresh `build(seed)` (the acceptance gate). A `build` that branches structurally on its seed is out of contract for mode (b) and must use `:rebuild` — the guard is documented, not auto-detected.

### 14.3 Results export bundle

`export_run(state, dir)` and `export_ensemble(ens, dir)` write a directory bundle implementing the §8.5 outputs clause, format per artifact: rectangular data (`sol`, `program_ledger`, ensemble summaries) → CSV (human) + Arrow (faithful); heterogeneous/semi-structured data (the `log` event stream, the §14.1 token trajectories, the run manifest) → JSON. The bundle layout is the §8.5 `runs/<model_hash>/<seed>/{trajectory.{csv,arrow}, ledger.{csv,arrow}, events.json, tokens.json, run.json}`, with `ensemble.json` (per-member seeds + the `summarize` table) at the ensemble root. Dependency tiering preserves the ADR 0005 minimalism: JSON always-available (reuses `serialize.jl`'s Dict encoding), CSV via the present `CSV.jl` dep, Arrow behind the `RDArrowExt` weakdep. The `events.json`/`tokens.json` streams round-trip-test like the model JSON (the §8 round-trip guarantee, ADR 0005 symmetry).

### 14.4 Invariants

(1) Trajectory rows are appended in `token_sortkey` order each tick; `(model, seed)` determines the log (§4 D4); the logging hook is 𝓕ₜ-measurable. (2) The trajectory log is bounded by per-kind opt-in (default: log nothing). (3) Ensemble members are independent, seeded `hash((root_seed,k))`; the ensemble is determined by `(root_seed, nseed, build)` (§4 D8/D9). (4) `EnsembleProblem` resolves `getobservable`/`observables` for every exported aggregate (no AA `@error` fall-through, §13.6). (5) Rectangular artifacts round-trip Arrow byte-faithfully; JSON streams round-trip structurally (§8); the manifest pins `(model_hash, seed)` (§4 D6). (6) All of §14 is read-only — no dynamics mutation, no re-indexing.

## 15. Visualization (see [ADR 0014](adr/0014-visualization.md))

The engine ships one live generic plot — `AlgebraicAgents._draw(prob, vars)`, place trajectories from `prob.sol` (`plots.jl:13`) — plus dead SciML plotting code (`solve.jl`'s `@plot`/`plot_summary`/`plot_ensemble_sol`, orphaned with the ADR 0001 SciML demotion) and no way to draw the reactive NETWORK. This section pins two additive capabilities, sharing the Plots/Graphviz infrastructure RD and AA already carry, and removes the dead code.

### 15.1 Result-plot recipes

A model-agnostic set of Plots.jl `@recipe`s behind a `RDPlotsExt` package extension, each keyed off a documented raw artifact so it works on any model: (1) marking trajectory (`prob.sol`, generalizing `_draw`); (2) resource utilization/saturation (`:saturation`/`:allocation` log rows); (3) valuation/rNPV curve (`:valuation` rows + `program_ledger`); (4) per-program ledger (the §14 DataFrame); (5) token trajectory + the §14.1 envelope; (6) ensemble distribution / treatment-effect waterfall (the §14.2 data); (7) throughput (`:new_transitions`/`:terminated_*`). The dead SciML plotting (`solve.jl:7-201`, the `EnsembleSummary`/`EnsembleSolution` paths + `@plot`) and the never-called `plot_df` (`plots.jl:3-10`) are removed; `@agentize` and the generalized `_draw` are kept.

### 15.2 The network "exec map"

Three separable layers. Layer A — `network_graph(prob_or_spec) -> NetworkGraph`: dependency-free structure extraction (place/transition nodes + arcs with multiplicity/modality), built by walking the incidence (`transLHS`/`transRHS` today, `ArcSpec` post-ADR-0003), round-trippable from JSON with no run. Layer B — `to_graphviz`/`draw_network`: emit Petri-net DOT (arc color by §1 modality, transition nodes annotated with rate/priority/cycletime) rendered via AA's `run_graphviz` (`graphviz.jl:18`), and COMPOSED with AA's `wiring_diagram` (`wires.jl:144`) by rendering the RD net as a Graphviz cluster so §13 `inputs[]` coupling wires appear on the same map — one picture of the whole hierarchy. Layer C — `exec_map(prob_or_ensemble; highlight)`: decorate Layer A with §14 statistics — place fill = marking level/trough (starvation), transition color = throughput/saturation (bottlenecks hot), arc thickness = flow, and token-path highlighting from `past_bonds` (`agents.jl:141`) scoped by a `@select` `TokenPredicate` (§9.5). The result→map integration is a styling pass over §14 data, not a new computation.

### 15.3 Primary view

The Petri net is the PRIMARY view (Layers A/B); the AA agent tree is offered as the composed OUTER FRAME via Layer B's cluster embedding, so the two reconcile into one image rather than competing. Tokens are NOT drawn as their own nodes (there can be thousands) — they appear only as Layer C overlay decoration on the place/transition they occupy.

### 15.4 Invariants

(1) `network_graph` is a pure function of the spec — no plotting dep, no simulation (§6). (2) Rendering reuses AA's `run_graphviz`/`wiring_diagram`; no new graph library. (3) `exec_map` only reads finished-run/ensemble statistics (§14) — no mutation, no re-run. (4) Token highlighting is a `TokenPredicate` (§9.5) — a highlighted cohort is a `@select` set. (5) Every recipe consumes a documented raw artifact — none hard-codes BD fields. (6) No dead viz code — every shipped entry point runs on `ref-agents`.

---

## Status

The Phase-0 modeling contract is COMPLETE and SIGNED OFF (§1–§9): §1 modality truth table, §2 time model (+§2.8 genesis modes), §3 operational semantics, §4 determinism & seeding, §5 attribute contract, §6 object model, §7 composition semantics, §8 serialization schema, §9 structured tokens & queries (+§9.5 predicate selection). A Phase-0.5 modeling-language extension increment adds §10 interface & initial-state contract, §11 refinement & open-port composition, §12 rules/triggers & conditional transitions (with the §12.3 action family extended by [ADR 0011](adr/0011-action-callbacks-and-general-code.md)), and §13 AlgebraicAgents integration & external coupling, with §9.5 extending §9 — the business-process modeling refinements (compact/compositional authoring, granularity/substitution, declarative initial state + state dump, agentic-place filtration, the endogenous decision channel, full action expressivity — population writes + a general-code escape hatch — and reactive-network-as-AA-hierarchy-node coupling — that make an expressive BD demo possible). It is backed by [ADR 0001](adr/0001-discrete-event-engine.md) (engine), [ADR 0002](adr/0002-priority-weighted-allocation.md) (allocation), [ADR 0003](adr/0003-data-store.md) (data store), [ADR 0004](adr/0004-runtime-mutation.md) (runtime mutation), [ADR 0005](adr/0005-serialization-json-ir.md) (JSON serialization), [ADR 0006](adr/0006-structured-tokens.md) (structured tokens + custom-function registry), [ADR 0007](adr/0007-interface-and-initial-state.md) (interface & initial state + dump/restore), [ADR 0008](adr/0008-token-filtration.md) (token filtration / predicate selection), [ADR 0009](adr/0009-refinement-and-composition.md) (refinement & open-port composition), [ADR 0010](adr/0010-rules-and-conditional-transitions.md) (rules/triggers & conditional transitions), [ADR 0011](adr/0011-action-callbacks-and-general-code.md) (action callbacks: population writes + general-code `Invoke`), and [ADR 0012](adr/0012-algebraicagents-integration.md) (AlgebraicAgents integration: RD as a hierarchy node + wired external coupling). A Phase-0.6 analysis & visualization increment adds §14 analysis & observability (per-token trajectory log, ensemble runner, results export bundle) and §15 visualization (result-plot recipes + the network exec map), backed by [ADR 0013](adr/0013-analysis-and-observability.md) and [ADR 0014](adr/0014-visualization.md) — the result-inspection layer (inspect an ensemble of runs, gather ledger statistics, typical agentic-token trajectories, export, plot, and a result-decorated network map for (in)efficiency analysis) that sits on top of the §1–§13 modeling layer without disturbing it. **§14/§15 are IMPLEMENTED (2026-06-30, commit `6174ebe`)** — the first `src/` change of the Phase-0 line, landed as `src/analysis.jl`/`src/export.jl`/`src/visualize.jl` + `ext/RDPlotsExt.jl`/`ext/RDArrowExt.jl` with their semantic tests green (suite 500 pass / 7 broken / 507); the §14.2 ensemble reinit-reseed mode (b) landed subsequently once the §10.4/ADR 0007 §D `_reinit!` completion closed its gate (the `reinit!(m; seed)` reseed path + the `mode = :reinit` ensemble wiring, member-for-member equivalent to mode (a)). The [MVP_BD_DEMO.md](../demo/bd_acquisition/MVP_BD_DEMO.md) acquisition-impact demo is the exercise that drove §12 and catalogues the contract weak-points it surfaced; the Phase-0.5 increment resolved findings B (§12), C (§9.5), F (§11), G (§13 `getobservable`), H (§10.3), and I (ADR 0011 `SetTokens`/`Invoke`), leaving open A (counterfactual RNG, §4.6 future work) and D / D-bis (per-program ledger + TVE-freezing semantics); the §14 ensemble runner + export and §15 recipes promote the demo's hand-rolled `analysis.jl`/`figures.jl` into the engine. Every `file:line` citation in §1–§9 was adversarially verified against the current `ref-agents` source (several verifiers ran the engine on Julia 1.12.5); the §10/§11/§12/§13/§9.5 citations are read-verified against the same tree (and against the AlgebraicAgents source for §13) and remain PROPOSED pending maintainer sign-off; the §14/§15 citations were read-verified against the same tree (and the AlgebraicAgents source for the §15 Graphviz/wiring reuse) and are now IMPLEMENTED and exercised by the live engine (the implementation pass corrected several — e.g. the trajectory seam, the AA reuse points against the pinned GitHub source). The Phase-0 semantic test suite that encodes the §1–§9 invariants lives under `test/semantic/`; the §14/§15 acceptance tests have LANDED there (`analysis_observability.jl`, `visualization.jl`), while the §10/§11/§12/§13/§9.5 acceptance tests are still to be added with their ADRs. With maintainer sign-off on the still-proposed sections this closes the Phase-0 gate; broader implementation (Phase 1) may then begin.
