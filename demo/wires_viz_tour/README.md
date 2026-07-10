# Wires & ports, made visible — a visual tour of AlgebraicAgents coupling

A presentation-quality, **visual** walkthrough of ReactiveDynamics as a node inside a larger AlgebraicAgents (AA) hierarchy, coupled to sibling agents through **wires** (open ports carrying one agent's observable into another's input). This is ADR 0012 / CONTRACT §13. The companion [`aa_integration`](../aa_integration) demo teaches the same coupling in prose and `println`; this tour draws it — the wired hierarchy as a diagram, the RD net as a diagram, and the coupled-run data as a table and a chart, intertwined with the narration.

The tutorial is authored as a **Literate.jl** script (`wires_literate.jl`) and rendered to a self-contained HTML page (`wires_literate.html`). It sits on `wires_model.jl` (the wired-hierarchy build + the two Graphviz renders + the data-shaping), lifted verbatim from `aa_integration` — the presentation layer is a thin wrapper over it, so the model-building logic is written once.

## The story

A small pharma R&D portfolio under one `FreeAgent("portfolio")` root:

- a **market** agent (`MarketAgent`) exporting a drifting `sentiment` observable,
- the **RD net** (`ReactionNetworkProblem`) with an `inputs[]` `sentiment` port — growing `cash` at a rate driven by that external sentiment, and firing an acquisition lever when both external sentiment and internal cash clear their thresholds,
- a **finance** agent (`FinanceAgent`) reading the RD net's `cash` back off a wire.

Two wires close the loop: `market.sentiment ▶ RD.sentiment` (inbound) and `RD.cash ▶ finance.rd_cash` (outbound). One `simulate(root)` runs the whole thing. The two headline visuals are `AlgebraicAgents.wiring_diagram(root)` (the cross-agent view — agents, parentship edges, the two labelled wires) and `draw_network(rd)` (the intra-agent view — the pharma net's Petri structure inside the RD box). The coupled-run table and chart show the one-tick **Jacobi lag**: RD `cash` responds to `sentiment` one tick late because external reads are latched once per tick at `_prestep!`.

Every call in the tour already exists and is exercised by `aa_integration` / `introspection_tour` / the tests — no API is invented.

## Rebuild

The HTML is regenerated from the source by `build.jl`:

```bash
julia --project=demo/wires_viz_tour demo/wires_viz_tour/build.jl
```

First time (resolve the demo-local env — a few minutes):

```bash
julia --project=demo/wires_viz_tour -e 'using Pkg; Pkg.instantiate()'
```

`build.jl` renders the tour end-to-end (it executes the tutorial, so `dot` produces real SVGs). Long Julia compiles buffer their output — when driving this from a watchdog'd shell, launch it in the background and poll the log.

### Why a demo-local `Project.toml`

The main ReactiveDynamics project is lean (`Plots` is a weakdep). This tour needs `Plots` (to render the coupling chart and, via `RDPlotsExt`, `draw_network`'s Graphviz path) plus `Literate` for the HTML render. `AlgebraicAgents` (0.4, from the General registry) resolves as an ordinary dependency — no `[sources]` pin. Only `ReactiveDynamics = {path = "../.."}` is a path source.

## How the HTML is built (the exact APIs)

**Literate → markdown → HTML.** Literate.jl's native targets are markdown / notebook / script (no HTML target in this version), so `build.jl` renders Literate → markdown with executed outputs (`Literate.markdown(src, out; execute = true, flavor = CommonMarkFlavor())`). That writes the Plots chart and each network diagram as sidecar `.svg` files referenced by `![](name-N.svg)` — the diagrams are returned as `image/svg+xml`-showable `RawSVG` objects so Literate captures them. `build.jl` then converts the markdown to HTML with the stdlib `Markdown` and **inlines** every sidecar SVG in place of its `<img>` tag, yielding one self-contained HTML with the diagrams embedded. (CommonMark.jl is not in the resolved env; the stdlib `Markdown` renderer is sufficient here.)

## Files

| File | Role |
|---|---|
| `wires_model.jl` | Shared substrate: `MarketAgent` / `FinanceAgent`, `PHARMA_JSON`, `build_portfolio` / `run_coupled`, `wiring_diagram_svg` / `rd_network_svg`, `coupling_table` / `all_wires` (lifted from `aa_integration`) |
| `wires_literate.jl` | The Literate tutorial source |
| `build.jl` | Regenerates the HTML from the source |
| `wires_literate.html` | The rendered tour (committed) |
