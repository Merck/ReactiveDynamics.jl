# Wires & ports, made visible — a visual tour of AlgebraicAgents coupling

A presentation-quality, **visual** walkthrough of ReactiveDynamics as a node inside a larger AlgebraicAgents (AA) hierarchy, coupled to sibling agents through **wires** (open ports carrying one agent's observable into another's input). This is ADR 0012 / CONTRACT §13. The companion [`aa_integration`](../aa_integration) demo teaches the same coupling in prose and `println`; this tour draws it — the wired hierarchy as a diagram, the RD net as a diagram, and the coupled-run data as a table and a chart, intertwined with the narration.

The same tutorial is rendered in **two HTML forms** so you can compare them and pick one:

- **Form A — Pluto notebook** → `wires_viz_tour_pluto.html` (source `wires_viz_tour.jl`)
- **Form B — Literate.jl script** → `wires_literate.html` (source `wires_literate.jl`)

Both share one substrate, `wires_model.jl` (the wired-hierarchy build + the two Graphviz renders + the data-shaping), lifted verbatim from `aa_integration`. The two tutorial sources are thin presentation layers over it — the model-building logic is written once.

## The story (both forms)

A small pharma R&D portfolio under one `FreeAgent("portfolio")` root:

- a **market** agent (`MarketAgent`) exporting a drifting `sentiment` observable,
- the **RD net** (`ReactionNetworkProblem`) with an `inputs[]` `sentiment` port — growing `cash` at a rate driven by that external sentiment, and firing an acquisition lever when both external sentiment and internal cash clear their thresholds,
- a **finance** agent (`FinanceAgent`) reading the RD net's `cash` back off a wire.

Two wires close the loop: `market.sentiment ▶ RD.sentiment` (inbound) and `RD.cash ▶ finance.rd_cash` (outbound). One `simulate(root)` runs the whole thing. The two headline visuals are `AlgebraicAgents.wiring_diagram(root)` (the cross-agent view — agents, parentship edges, the two labelled wires) and `draw_network(rd)` (the intra-agent view — the pharma net's Petri structure inside the RD box). The coupled-run table and chart show the one-tick **Jacobi lag**: RD `cash` responds to `sentiment` one tick late because external reads are latched once per tick at `_prestep!`.

Every call in the tour already exists and is exercised by `aa_integration` / `introspection_tour` / the tests — no API is invented.

## Rebuild

Both HTML files are regenerated from their sources by `build.jl`:

```bash
julia --project=demo/wires_viz_tour demo/wires_viz_tour/build.jl
```

First time (resolve the demo-local env — a few minutes):

```bash
julia --project=demo/wires_viz_tour -e 'using Pkg; Pkg.instantiate()'
```

`build.jl` renders **both** forms end-to-end (each executes the tutorial, so `dot` produces real SVGs). Long Julia compiles buffer their output — when driving this from a watchdog'd shell, launch it in the background and poll the log.

### Why a demo-local `Project.toml`

The main ReactiveDynamics project is lean (`Plots` is a weakdep). This tour needs `Plots` (to render the coupling chart and, via `RDPlotsExt`, `draw_network`'s Graphviz path), plus `Pluto` + `PlutoUI` and `Literate` for the two HTML forms. `AlgebraicAgents` (0.4, from the General registry) resolves as an ordinary dependency — no `[sources]` pin. Only `ReactiveDynamics = {path = "../.."}` is a path source.

## How each form is built (the exact APIs)

**Form A — Pluto → HTML (headless, no server).** `build.jl` opens a `Pluto.ServerSession` with `run_notebook_on_load = true`, runs every cell to completion with `Pluto.SessionActions.open(session, path; run_async = false)`, checks no cell errored, then bakes the executed state into a single self-contained page with `Pluto.generate_html(notebook)`. The notebook's first cell calls `Pkg.activate(@__DIR__)`, which is how Pluto detects manual package management and uses *this* demo env instead of its own nbpkg. The SVGs and the Plots chart are embedded via `HTML(...)` cells.

**Form B — Literate → markdown → HTML.** Literate.jl's native targets are markdown / notebook / script (no HTML target in this version), so `build.jl` renders Literate → markdown with executed outputs (`Literate.markdown(src, out; execute = true, flavor = CommonMarkFlavor())`). That writes the Plots chart and each network diagram as sidecar `.svg` files referenced by `![](name-N.svg)` — the diagrams are returned as `image/svg+xml`-showable `RawSVG` objects so Literate captures them. `build.jl` then converts the markdown to HTML with the stdlib `Markdown` and **inlines** every sidecar SVG in place of its `<img>` tag, yielding one self-contained HTML with the diagrams embedded. (CommonMark.jl is not in the resolved env; the stdlib `Markdown` renderer is sufficient here.)

## Files

| File | Role |
|---|---|
| `wires_model.jl` | Shared substrate: `MarketAgent` / `FinanceAgent`, `PHARMA_JSON`, `build_portfolio` / `run_coupled`, `wiring_diagram_svg` / `rd_network_svg`, `coupling_table` / `all_wires` (lifted from `aa_integration`) |
| `wires_viz_tour.jl` | Form A source — the Pluto notebook |
| `wires_literate.jl` | Form B source — the Literate script |
| `build.jl` | Regenerates both HTML files from the sources |
| `wires_viz_tour_pluto.html` | Form A output (committed) |
| `wires_literate.html` | Form B output (committed) |
