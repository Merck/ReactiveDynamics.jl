# BD Acquisition-Impact Demo

A reproducible counterfactual on a living pharma pipeline portfolio: run the same pipeline once *without* an acquisition and once *with* it, and read off the change in portfolio value, launches, and capital need *attributable to the deal* — including the deal's effect on the programs the company already owns. This is the [MVP_BD_DEMO.md](../../docs/MVP_BD_DEMO.md) Milestone-1 scope made executable against the Phase-1 engine.

## Run it

```bash
julia --project=. demo/bd_acquisition/run_demo.jl
```

(Or, against the dev test server, `dev/run.sh demo/bd_acquisition/run_demo.jl`.)

## What it shows

The headline is **Δ-rNPV** — the deal's attributable change in risk-adjusted portfolio NPV, computed as an **ensemble-averaged difference of means** across N seeds (not a per-seed paired difference; the contract keeps a single state-owned RNG stream, so the two scenarios desync after the deal — [MVP §4.1](../../docs/MVP_BD_DEMO.md) finding A). A representative 24-seed run:

| Scenario | mean rNPV | mean launches | P(≥1 launch) | Δ-rNPV vs S0 |
|---|--:|--:|--:|--:|
| **S0** Baseline (no deal) | 2117 | 1.17 | 0.75 | — |
| **S1** Deal, pipeline-only | 2622 | 1.79 | 0.88 | +504 |
| **S2** + resource synergy | 2482 | 1.67 | 0.88 | +365 |
| **S3** + capability/PoS synergy | 4625 | 3.42 | 0.92 | **+2508** |
| **S4** + op-efficiency synergy | 2903 | 2.08 | 0.96 | +786 |
| **S5** Full (all synergies) | 4356 | 3.25 | 0.96 | **+2238** |

The decomposition is the BD insight: most of the deal's value comes from the **capability/PoS synergy raising the success probability of the programs already in the pipeline** (S3), not from the acquired programs themselves (S1). A partner can read off "the deal is worth ~+2200 net of a 400 price, and the value is in the platform synergy — so don't pay for it as a pipeline-only bolt-on."

## How it maps to the engine (Phase-1 capabilities)

| Demo element | Engine mechanism |
|---|---|
| Program / asset | a `ProjectToken` structured token (host Julia, [`host.jl`](host.jl)), identity preserved across phases |
| Pipeline phase | a `phase` **attribute** on the single `:Project` kind (phase-as-attribute, ADR 0008) — not a species per phase |
| Phase advance | a transition selecting an in-phase token via `@select(Project, phase==:PhaseN)` (Stage C), consuming `@conserved(scientist)` + `@rate(budget)`, advancing via `@advance(phase, :PhaseN1)` on `Binomial(q, PoS)` success |
| Failure / kill | `Binomial` failure ⇒ the bound token soft-retires (its species flips to `:removed`, ADR 0006); its `phase` records how far it got |
| **The acquisition** | an **endogenous Rule** (ADR 0010, Stage B): `fire_mode: once`, guard `@t() > T_acq`, action `Seq[AddToken(ProjectToken…), SetSpecies(scientist,+Δ), SetParams(synergy…)]`. The lever lives *in the model*, not in host patch code |
| Synergies | param-mediated (MVP §2.1): the acquisition rule flips `synergy_pos`/`synergy_eff`, which the late-phase transitions read in their `probability`/`cycletime` ExprNodes |
| Determinism / counterfactual | every draw routes through the state-owned RNG seeded by `seed=` (Stage A, §4 D1–D9); each ensemble member is seeded `hash((root_seed, k))` (D8); a run is fully determined by `(model, scenario, seed)` |

## Files

- [`host.jl`](host.jl) — the `ProjectToken` kind (host Julia, never serialized), the per-network registry the `AddToken` lever references by name, the coarse pipeline model builder, the initial portfolio, and the acquisition Rule.
- [`analysis.jl`](analysis.jl) — the rNPV roll-up and ensemble/Δ post-processing. Discounting is pure post-processing; the engine does no discounting (MVP §5 / finding D).
- [`run_demo.jl`](run_demo.jl) — the scenario grid (S0–S5) + driver + report.

## Milestone-1 caveats (deliberate)

- **Initial portfolio is imperative host code** (`seed_portfolio!`), not the declarative `population[]` array — that is Phase-1 Stage D (ADR 0007). The starting pipeline is therefore not yet part of the serializable model document.
- **The model is built in-Julia via `@ReactionNetworkSchema`**, not authored as eval-free JSON — that is Phase-1 Stage E (ADR 0005). The demo's *trajectory* does not depend on serialization.
- **Synergies are param-mediated** (a Ref to a param the rule flips); token-pool-mediated synergy (a PoS reading sibling token counts via `TokenAgg`) is a clean Milestone-2 extension.
- **Refinement** (the §11 `refine` granularity-substitution demo) is Milestone-2.
