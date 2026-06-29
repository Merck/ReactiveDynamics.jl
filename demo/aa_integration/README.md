# AlgebraicAgents integration & external coupling

A single, runnable, literate walkthrough of ReactiveDynamics as a NODE inside a larger AlgebraicAgents (AA) hierarchy — coupled to sibling agents in BOTH directions and driven by one `simulate(root)`. This is the executable form of ADR 0012 / CONTRACT §13. The companion [agentic_pipeline](../agentic_pipeline) demo shows the in-model decision channel and the eval-free JSON model on a standalone reaction network; this one takes that same network and drops it into a multi-agent hierarchy, wired to a macro/market source and a downstream finance sink.

## Run it

```bash
julia --project=. demo/aa_integration/aa_integration.jl
```

The script is literate: every section opens with a block comment explaining the integration idea, then runs the code, then `println`-narrates the result, so running it once tells the whole story top to bottom.

## The integration story in one paragraph

A `ReactionNetworkProblem` is ALREADY an AA `@aagent` implementing the three stepping hooks (`_step!`/`_reinit!`/`_projected_to`), so `entangle!(parent, rd)` makes it a hierarchy node and AA's least-projected-time gate interleaves its single clock with sibling clocks for free — no new clock code. What ADR 0012 adds is the COUPLING surface in both directions. OUTBOUND: RD implements `observables`/`getobservable` (so an AA wire can ORIGINATE from an RD net — a sibling can read an RD species or a named observable) and `_getparameters`/`_setparameters!` (so the hierarchy can read or patch RD params, param-only and index-safe). INBOUND: RD reads external state through a DECLARATIVE, eval-free channel — a model-local `inputs[]` port list plus a new closed `ExternalRef` leaf in the expression IR, which a rate, a guard, or an action value can read by NAME. The wiring (which foreign agent feeds which port) lives host-side in `add_wire!`, never in the RD document.

## What the demo exercises

| § | Section | Capability exercised |
|---|---|---|
| 0 | Sibling agents | A `MarketAgent` (source, exports `:sentiment`) and a `FinanceAgent` (sink, reads RD's `cash` at its own `_prestep!`) — ordinary AA `@aagent`s in the host, knowing nothing of RD's internals |
| 1 | The RD net + its `inputs[]` port | The pharma net as an eval-free JSON model declaring one read port `sentiment` (with a pre-wire default); `ExternalRef(sentiment)` read in BOTH a transition RATE and an acquisition RULE GUARD; `validate` rule 8 enforces every `ExternalRef.port` is a declared `inputs[]` port (an undeclared port is a diagnostic, never an eval) |
| 2 | Compose & wire (host-side) | `entangle!` the three agents under one root; `add_wire!` lays `market.sentiment ▶ RD.sentiment` (inbound) and `RD.cash ▶ finance.rd_cash` (outbound). RD's `getobservable` is what makes the outbound wire possible |
| 3 | One `simulate(root)` | AA's least-projected-time gate interleaves the clocks; `_prestep!` latches external reads once/tick; the acquisition lever fires on the first tick BOTH the external sentiment clears its threshold AND internal cash clears its trigger (the ADR North-star: a decision driven by external AND internal state) |
| 4 | The OUTBOUND read | The finance agent reconstructs RD's cash trajectory purely from `getobservable(:cash)` carried over its wire — RD is a first-class wire SOURCE |
| 5 | Determinism | Same `(hierarchy, seed)` ⇒ identical coupled trajectory (the whole reason external reads are pinned at `_prestep!`) |
| 6 | reinit | `_reinit!` restores RD's external-input buffer to its declared default, so latched wire values from a finished run do not leak into the next |
| 7 | Recap | A closing summary of the integration and the ADRs it realizes |

## The load-bearing design points (narrated in-line)

These are the guarantees the demo makes executable, surfaced where they happen:

- **RD is a hierarchy node, clocks for free (§C).** AA's `step!` projects the whole hierarchy to the least projected time and steps each agent only when its `_projected_to` equals that minimum. RD's `_projected_to` returns `state.t` until the horizon, so an RD net interleaves correctly with siblings stepping at any `dt` — the ADR 0001 single-clock contract, now read as ONE node's clock under AA's multi-clock coordinator. A coupled net is driven by `simulate(root)`, not `simulate(rd)` directly.
- **The one-tick Jacobi coupling lag (§B3, the determinism pin).** External inputs are read ONCE per tick at `_prestep!`, which `step!` `prewalk`s over the WHOLE hierarchy in its FIRST phase, before ANY agent's `_step!` runs that tick. So RD always reads each source's value as projected to the PREVIOUS tick boundary — a one-tick lag, no algebraic loop. You can SEE the lag in §3: RD's cash stays 0 at t=1 because the sentiment latched at t=0 was the pre-wire 0.0; the t=1 sentiment is only read at the t=2 boundary. This is explicit (Jacobi-style) co-simulation: fully reproducible under `(hierarchy, seed)`, independent of AA's `Dict`-order sibling stepping.
- **Why `_prestep!`, not a live read (the §4-D4 hazard it closes).** AA's `step!` advances sibling agents in raw `Dict`-iteration order. A LIVE cross-agent `getobservable` mid-`_step!` would see a sibling's pre- or post-step value depending on hash order — a determinism violation. Reading at the pinned `_prestep!` removes that dependence entirely: every `ExternalRef` read in a tick returns the same buffered value, so a rate and a guard that read the same port agree.
- **Eval-free coupling (Invariant 4).** The RD JSON declares the ports it CONSUMES; it never names the foreign agents that fill them. The cross-hierarchy topology lives only in host `add_wire!` calls, so the model document stays portable and no foreign Julia is ever parsed or eval'd.

## Honest scope

The coupling here is EXPLICIT (Jacobi) — a one-tick lag, no within-tick fixed point. An implicit/algebraic coupling (a tight feedback loop solved to a fixed point each tick) is deliberately out of scope; it would be an AA-level `Opera` interaction, a separate design. The `grow` source uses a fractional `@deterministic` rate for didactic clarity at the fixed `dt = 1.0` of this demo; a fractional deterministic count on a source is not `dt`-invariant in general (see the [core_engine_tour](../core_engine_tour) §5 note), so use integer deterministic counts or a stochastic rate for a `dt`-varying study.
