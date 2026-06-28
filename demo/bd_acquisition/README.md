# BD Acquisition-Impact Demo

A reproducible counterfactual on a living pharma pipeline portfolio: run the same pipeline once *without* an acquisition and once *with* it, and read off the change in portfolio value, launches, and capital need *attributable to the deal* — including the deal's effect on the programs the company already owns. This is the [MVP_BD_DEMO.md](../../docs/MVP_BD_DEMO.md) Milestone-1 scope made executable against the Phase-1 engine.

## Run it

```bash
julia --project=. demo/bd_acquisition/run_demo.jl
```

(Or, against the dev test server, `dev/run.sh demo/bd_acquisition/run_demo.jl`.)

## What it shows

The headline is **Δ-rNPV** — the deal's attributable change in risk-adjusted portfolio NPV, computed as an **ensemble-averaged difference of means** across N seeds (not a per-seed paired difference; the contract keeps a single state-owned RNG stream, so the two scenarios desync after the deal — [MVP §4.1](../../docs/MVP_BD_DEMO.md) finding A). The default run uses 160 seeds, because the per-scenario SE (~±450 at 24 seeds) otherwise swamps the synergy decomposition. The canonical 160-seed run (root seed 2026):

| Scenario | mean rNPV | launches | P(≥1) | cash⌄ | sci⌄ | Δ-rNPV vs S0 (±1 SE) |
|---|--:|--:|--:|--:|--:|--:|
| **S0** Baseline (no deal) | 2282 | 1.21 | 0.81 | 16 | 8 | — |
| **S1** Deal, pipeline-only | 2736 | 1.76 | 0.88 | 16 | 7 | +454 ± 152 |
| **S2** + resource synergy | 3092 | 2.16 | 0.95 | 16 | 9 | +810 ± 168 |
| **S3** + capability/PoS synergy | 3287 | 2.15 | 0.98 | 16 | 7 | **+1005 ± 145** |
| **S4** + op-efficiency synergy | 2750 | 1.83 | 0.89 | 16 | 7 | +468 ± 158 |
| **S5** Full (all synergies) | 3945 | 2.84 | 0.97 | 33 | 16 | **+1663 ± 182** |

`cash⌄`/`sci⌄` = the low-water marks of the budget / scientist pools (near zero ⇒ that resource is binding). The company's own financing dip is exactly `150 − cash⌄` (a 150-unit starting reserve), injection-robust — so a deal that injects capital can only shrink it; it does not appear as a separate column. Marginal synergy contributions over the pipeline-only deal (S1): **resource +356, capability/PoS +551, op-efficiency +14**.

Three reads a BD partner takes away:

1. **The value is in the platform, not the pipeline.** The largest single synergy is **capability/PoS — the target's platform raising the success probability of the programs the company *already owns*** (S3 marginal +551), not the acquired programs themselves (S1 +454). Don't price the deal as a pipeline-only bolt-on.
2. **Cash is the binding constraint, and the deal absorbs its own capital.** The budget pool runs to its floor (`cash⌄ ≈ 16`) in every scenario through S4 — the pipeline is cash-starved. Resource synergy (S2) injects capital, but it is *fully absorbed into running more programs in parallel*: the cash trough stays pinned at the floor and the financing dip doesn't shrink. Only the **full deal (S5)** — where capability synergy also pushes more programs through to launch and out of the pool — finally leaves real slack (cash⌄ 16 → 33, sci⌄ 8 → 16, own-reserve dip 134 → 117). That capital-is-consumed-not-banked dynamic is exactly the system-level effect a spreadsheet rNPV misses. *(The naive "max − min budget" reading would have shown the ask doubling here — an artifact of counting the injected capital itself; the injection-robust dip-below-start is flat-to-down.)*
3. **Operational-efficiency synergy is ≈ 0 here** (+14, within noise). An honest, non-obvious finding: shorter cycle times barely help when **capital, not time, is the binding constraint** (the cash pool floors out while scientists rarely do). The same deal in a time-constrained pipeline would value op-efficiency very differently — which is the point of having a dynamic model.

## How it maps to the engine (Phase-1 capabilities)

| Demo element | Engine mechanism |
|---|---|
| Program / asset | a `ProjectToken` structured token (host Julia, [`host.jl`](host.jl)), identity preserved across phases |
| Initial portfolio | the declarative `population[]` initial marking (ADR 0007 §B, Stage D): `initial_population()` is a list of `ProjectToken` structs passed to the constructor, instantiated before t=0 — reproducible input, not imperative post-construction host code (MVP finding H) |
| Pipeline phase | a `phase` **attribute** on the single `:Project` kind (phase-as-attribute, ADR 0008) — not a species per phase |
| Phase advance | a transition selecting an in-phase token via `@select(Project, phase==:PhaseN)` (Stage C), consuming `@conserved(scientist)` + `@rate(budget)`, advancing via `@advance(phase, :PhaseN1)` on `Binomial(q, PoS)` success |
| Failure / kill | `Binomial` failure ⇒ the bound token soft-retires (its species flips to `:removed`, ADR 0006); its `phase` records how far it got |
| **The acquisition** | an **endogenous Rule** (ADR 0010, Stage B): `fire_mode: once`, guard `@t() > T_acq`, action `Seq[AddToken(ProjectToken…), SetSpecies(scientist,+Δ), SetSpecies(budget,+Δ), SetParams(synergy…)]`. The lever lives *in the model*, not in host patch code |
| Resource contention | the model is calibrated so cash (and headcount) **genuinely bind** — the organic pipeline runs the `budget`/`scientist` pools into single digits, so the ADR-0002 priority allocator actually rations scarce resources (an over-provisioned model would leave the engine's contention machinery idle and make resource synergy inert) |
| Synergies | param-mediated (MVP §2.1): the acquisition rule flips `synergy_pos`/`synergy_eff`, which the late-phase transitions read in their `probability`/`cycletime` ExprNodes |
| **Per-program ledger** | the engine attributes the cost ledger **per program during the run** (MVP finding D, [`src/ledger.jl`](../../src/ledger.jl)): `budget` is priced via `specCost` (which touches only the ledger, not the dynamics), so each program's capital burn is accrued onto it at every advance it sits in. `program_ledger(prob)` returns the per-program cost/reward/valuation summary in deterministic token order; the per-program rows + an `unattributed` bucket reconcile exactly to the aggregate `:valuation_cost` (§8.5). Previously this was reconstructed in post |
| Determinism / counterfactual | every draw routes through the state-owned RNG seeded by `seed=` (Stage A, §4 D1–D9); each ensemble member is seeded `hash((root_seed, k))` (D8); a run is fully determined by `(model, scenario, seed)` |

## Files

- [`host.jl`](host.jl) — the `ProjectToken` kind (host Julia, never serialized), the per-network registry the `AddToken` lever references by name, the coarse pipeline model builder, the initial portfolio, and the acquisition Rule.
- [`analysis.jl`](analysis.jl) — the rNPV roll-up and ensemble/Δ post-processing, plus `program_ledger_summary`/`program_economics`/`ledger_reconciliation` which read the engine-level **per-program ledger** (MVP finding D, now resolved for cost) and cross-check it against the aggregate row and the post-hoc rNPV roll-up. Discounting/valuation stay pure post-processing; the engine does no discounting (MVP §5 / finding D's documented boundary).
- [`run_demo.jl`](run_demo.jl) — the scenario grid (S0–S5) + driver + report (default 160 seeds; the table prints Δ-rNPV ± SE, the financing dip, the cash/scientist low-water marks, and the marginal synergy decomposition).
- [`model.rdj.json`](model.rdj.json) — the same pipeline as an **eval-free JSON model** (ADR 0005, Stage E). `from_json_model(read("model.rdj.json", String); registry = PROJECT_REGISTRY)` builds a model byte-for-byte identical to the in-Julia DSL under the same seed (verified, `serialization_ir.jl::E8`). The host `ProjectToken` type + registry stay host Julia (referenced by name); the JSON carries no code.
- [`figures.jl`](figures.jl) → [`figures/`](figures/) — the presentation charts (Δ-rNPV waterfall, synergy decomposition with ±SE bars, rNPV distribution, binding-cash-constraint view), generated from the same engine run.
- [`BRIEF.md`](BRIEF.md) — the one-page technical-executive brief tying the figures to the engine mechanisms.
- [`presentation.html`](presentation.html) — a self-contained static HTML walkthrough (no server, no build step to view) structured as an HBR-style technical paper: it leads with **the framework as reusable architecture** for rapid valuation/impact modeling, then uses the BD acquisition as one worked case (assume → declare → simulate → read out) with an interactive counterfactual console and charts rendered client-side as SVG from the real 160-seed ensemble. [`export_data.jl`](export_data.jl) dumps that ensemble to `presentation_data.json`; [`build_presentation.jl`](build_presentation.jl) inlines it into the page. Regenerate with `julia --project=. demo/bd_acquisition/export_data.jl && julia demo/bd_acquisition/build_presentation.jl`.

## Milestone-1 caveats (deliberate)

- **Synergies are param-mediated** (a Ref to a param the rule flips); token-pool-mediated synergy (a PoS reading sibling token counts via `TokenAgg`) is a clean Milestone-2 extension.
- **Refinement** (the §11 `refine` granularity-substitution demo) is Milestone-2.
