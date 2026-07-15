# Core Engine Tour (classical species)

A single, runnable, literate walkthrough of the ReactiveDynamics engine's core modeling vocabulary, staying entirely in the CLASSICAL regime — every species is a plain counted quantity (a `Float64` stock: molecules, dollars, scientists, jobs). The companion [bd_acquisition](../bd_acquisition) demo covers STRUCTURED / agent tokens (entities with attributes and lifecycle identity); this one is the bare reaction engine, seen end to end without any of that machinery. Each construct is lifted from the package's passing semantic test suite (`test/semantic/*.jl`) — the tour invents no API.

## Run it

```bash
julia --project=. demo/core_engine_tour/core_engine_tour.jl
```

The script is literate: every section is a block comment explaining the modeling idea, then the code, then `println`-narrated results, so running it once tells the whole story top to bottom.

## The mental model in one paragraph

A model is a set of TRANSITIONS. Each transition has a RATE (how often it tries to fire), a left-hand side of REACTANTS it consumes, and a right-hand side of PRODUCTS it emits. Firing can be instantaneous (`cycletime` 0) or take time (`cycletime > 0`, an in-flight instance that may also fail a Binomial success draw). Reactants can be consumed outright, held-and-returned, metered per step, or freed each step — the engine's signature feature is this RESOURCE MODALITY algebra. When several transitions want the same scarce pool in one tick, a priority-weighted ALLOCATOR rations it. All randomness flows through a per-run seeded RNG, so a run is reproducible from `(model, seed)`.

## What each section exercises

| § | Section | Engine capability exercised |
|---|---|---|
| 1 | A first model: SIR | The metalanguage (`@reaction_network`, `@prob_init` / `@prob_params` / `@prob_meta`), mass-action stochastic rates, `simulate`, reading `prob.sol` by column name, `prob.u` / `find_index`, and a conserved-population invariant |
| 2 | Stateful transitions & lifecycle | `cycletime` (in-flight completion delay), Binomial `probability` of success, `capacity` bound on concurrent instances, `maxlifetime` timeout |
| 3 | Resource modalities truth-table | Raw-consumed (bare) vs `@conserved` (held, returned) vs `@rate` (metered per tick) vs `@nonblock` (freed each step) vs stacked `@rate(@conserved(...))` (rented hold) — including the `@rate`-with-`cycletime=0` foot-gun that silently reserves nothing |
| 4 | Priority allocator under contention | `progressive_fill!` (priority-weighted progressive filling) called directly in both the contended and the slack (no-scaling) regime, then genuine in-model contention for a scarce shared pool where the higher-`priority` transition wins more (ADR 0002) |
| 5 | Genesis modes | Poisson SOURCE (`∅` LHS, dt-invariant in expectation), token-gated FLOW / routing (high nominal rate clamped to available upstream tokens), and SCHEDULED batch intake via `@deterministic(N * @periodic(p))` |
| 6 | Registered rates & the ledger | `@register`ed custom rate functions on a toy-pharma pipeline, plus the `@cost` / `@reward` / `@valuation` ledger (`prob.log`) read by tag and reduced to a discounted rNPV (discounting is pure post-processing) |
| 7 | Composition | `@join` (union of species / transitions / params with shared-species identification via `@alias`) and `@equalize` (collapse two species into one and rewrite references); part counts via `nparts` |
| 8 | Determinism & ensembles | `seed=` reproducibility (same seed ⇒ identical trajectory; different seed ⇒ diverges; unseeded ⇒ fresh entropy seed), and a deterministically-seeded ensemble (`hash((root, k))`) with mean ± spread, member-`k` reproducible independent of N and order |
| 9 | Recap | A closing summary of everything shown |

## Honest caveats surfaced in the tour

These are deliberate, narrated in-line so a reader does not trip over them: `@rate` (and `@nonblock`) require `cycletime > 0` or the per-step draw is silently a no-op (§3d); a FRACTIONAL `@deterministic` count on a source is NOT dt-invariant, so use integer deterministic counts for sources (§5a note); and `@join` currently merges species / transitions / params but NOT observables (`:obs`) or events (`:E`) of the joined submodels (§7).
