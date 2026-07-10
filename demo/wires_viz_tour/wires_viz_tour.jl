### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 441fd186-7c61-11f1-bb7e-a5c123e2e57e
begin
    import Pkg
    # Run this notebook in the demo-local env (path-dev'd ReactiveDynamics + AlgebraicAgents 0.4 +
    # Plots), NOT Pluto's own package manager. Referencing `Pkg.activate` here is exactly how Pluto
    # detects "the author manages packages manually" and stands its own nbpkg down.
    Pkg.activate(@__DIR__)
    # The shared substrate: the wired-hierarchy build + the two Graphviz renders + the data-shaping.
    # Lifted verbatim from demo/aa_integration; the notebook is a thin presentation layer over it.
    include(joinpath(@__DIR__, "wires_model.jl"))
    using .WiresModel
    using ReactiveDynamics, AlgebraicAgents, DataFrames, Plots, PlutoUI
    TableOfContents(; title = "Wires & ports")
end

# ╔═╡ 4420bf4a-7c61-11f1-b365-4be8aa004f63
md"""
# ReactiveDynamics as an AlgebraicAgents node — the wires, made visible

A ReactiveDynamics reaction network is a first-class node in a larger **AlgebraicAgents (AA)** hierarchy. Coupled to sibling agents through **wires** — open ports carrying one agent's observable into another's input — it runs under a single `simulate(root)`. This is ADR 0012 / CONTRACT §13.

The [`aa_integration`](../aa_integration) demo teaches this coupling in prose and `println`. This tour makes it **visual**: it draws the wired hierarchy as a diagram and the RD net as a diagram, side by side with the coupled-run data. The subject is the *ports and wires you could not see before* — and *what lives inside the RD box they connect*.

The system is a small **pharma R&D portfolio** under one root:

- a **market** agent (a `MarketAgent`) that exports a drifting `sentiment` signal,
- the **RD net** (a `ReactionNetworkProblem`) with an `inputs[]` `sentiment` port, growing cash at a rate driven by that external sentiment and firing an acquisition lever when both external sentiment and internal cash clear their thresholds,
- a **finance** agent (a `FinanceAgent`) that reads the RD net's `cash` back off a wire.

Two wires close the loop: `market.sentiment ▶ RD.sentiment` (inbound) and `RD.cash ▶ finance.rd_cash` (outbound).
"""

# ╔═╡ 4420bfb0-7c61-11f1-b8e0-d1971bc7efd7
md"""
## 1 · Build the wired hierarchy

`build_portfolio` entangles the three agents under one `FreeAgent("portfolio")` root and lays the two wires with `add_wire!(root; from, to, from_var_name, to_var_name)` — the *only* place the cross-agent topology lives (the RD JSON names only its `sentiment` port, never the market agent that feeds it: Invariant 4, eval-free coupling). We build and run the whole thing with one `simulate(root)`, wrapped as `run_coupled`.
"""

# ╔═╡ 4420bfbc-7c61-11f1-8d1a-9d469ceb89b3
sys = run_coupled(seed = 11);

# ╔═╡ 4420bfc4-7c61-11f1-b5ab-c780881c3d7f
md"""
The hierarchy is `portfolio (root) ⊇ {market, RD net, finance}`. RD exports observables **$(AlgebraicAgents.observables(sys.rd))**, and the finance agent reads `cash` off its incoming wire. Here are the two wires, read straight off the shared `Opera`:
"""

# ╔═╡ 4420bfce-7c61-11f1-9aab-d3b02316b037
DataFrame(all_wires(sys.root))

# ╔═╡ 4420bfe2-7c61-11f1-92cb-b7a95081a069
md"""
## 2 · Visualize the wiring — the headline diagram

`AlgebraicAgents.wiring_diagram(root)` returns a Graphviz **DOT** string of the whole hierarchy — agent nodes, parentship edges (dashed), and the annotated **wires** — which `run_graphviz` renders to SVG (the same backend RD's `draw_network` uses). This is the cross-agent view: *who is wired to whom*.

Read it: `portfolio` (node 1) parents `market`, the RD net (`reaction_network`, which itself parents its `structured` token store), and `finance`. The two solid arrows are the wires — `market → reaction_network` labelled `sentiment` at both ends (the inbound port), and `reaction_network → finance` labelled `cash` at the tail / `rd_cash` at the head (the outbound wire, tail = the source observable, head = the target input).
"""

# ╔═╡ 4420bfee-7c61-11f1-b6b0-7107b3a6345d
HTML(wiring_diagram_svg(sys.root))

# ╔═╡ 4420bff6-7c61-11f1-85ab-6b56a961d342
md"""
## 3 · Visualize the RD node's internals

The wiring diagram shows the RD net as one box. `draw_network(rd)` opens that box: it renders the pharma net's **internal Petri structure** — the `cash` and `acquired` places (circles), the `grow` transition (box), and the arc `grow → cash`. Together with the wiring diagram above, this is the full coupling picture at two granularities: the cross-agent wiring, and the intra-agent structure the `sentiment` port feeds into.
"""

# ╔═╡ 4420c000-7c61-11f1-8269-f1fee2135a5c
HTML(rd_network_svg(sys.rd))

# ╔═╡ 4420c00a-7c61-11f1-b7d1-81ab333e526a
md"""
## 4 · Run the coupled model + show the data

`simulate(root)` already ran inside `run_coupled`. AA's least-projected-time gate interleaves the clocks; every tick, `_prestep!` latches each agent's incoming wires **once** — so RD reads the market's *previous* tick-boundary sentiment (an explicit, one-tick **Jacobi lag**, no algebraic loop). The table aligns, per RD tick: the market `sentiment`, RD `cash`, the finance agent's reconstructed `cash` (RD's cash read purely over the outbound wire), and the `acquired` lever.
"""

# ╔═╡ 4420c014-7c61-11f1-bde2-4193fe6f111f
tbl = coupling_table(sys)

# ╔═╡ 4420c020-7c61-11f1-b469-bb0a725d78c9
md"""
The **one-tick lag is visible in the data**: sentiment rises to `0.12` at `t=1`, but RD `cash` stays `0` until `t=2` — because the sentiment RD latched at `t=1` was the pre-wire default from the `t=0` boundary. And `finance_cash_seen` tracks `rd_cash` exactly: the finance agent reconstructed RD's trajectory knowing nothing of RD's internals, purely through `getobservable` on wire 2.
"""

# ╔═╡ 4420c028-7c61-11f1-b710-13f2ff6512ca
let
    plot(tbl.t, tbl.rd_cash; label = "RD cash", lw = 3, marker = :circle,
        color = :steelblue, legend = :topleft,
        xlabel = "tick t", ylabel = "cash", title = "The coupled trajectory (one-tick Jacobi lag)")
    plot!(tbl.t, coalesce.(tbl.finance_cash_seen, NaN); label = "finance saw (wire 2)",
        lw = 0, marker = :xcross, markersize = 7, color = :orange)
    sc = maximum(tbl.rd_cash) / maximum(tbl.sentiment)
    plot!(tbl.t, tbl.sentiment .* sc; label = "market sentiment (right-scaled)",
        lw = 2, ls = :dash, color = :seagreen)
    acq = findfirst(==(1), tbl.acquired)
    isnothing(acq) || vline!([tbl.t[acq]]; label = "acquisition fires", color = :firebrick, ls = :dot, lw = 2)
end

# ╔═╡ 4420c032-7c61-11f1-a731-47d2faf493ee
md"""
## 5 · Recap — the coupling surface (ADR 0012 / §13)

- **RD is a first-class hierarchy node.** `entangle!` makes it a child; AA's least-projected-time gate interleaves its single clock with the siblings' for free (no new clock code). A coupled net is driven by `simulate(root)`, never `simulate(rd)`.
- **Inbound coupling** is a declarative, eval-free port: the RD JSON's `inputs[]` list plus a closed `ExternalRef` leaf that a rate or a guard reads by name. Here `sentiment` drives both the `grow` rate and the `acquire` guard.
- **Outbound coupling** is `getobservable`: any sibling can read an RD species or named observable off a wire, as `finance` reads `cash`.
- **The wiring lives host-side** in `add_wire!` — never in the portable RD document (Invariant 4).
- **The coupling is explicit and deterministic**: external reads are latched once per tick at `_prestep!` (the one-tick Jacobi lag), so a coupled run is a function of `(hierarchy, seed)` alone.

The two diagrams above are the two halves of that story: `wiring_diagram` shows the ports and wires *between* agents; `draw_network` shows the Petri structure *inside* the RD node they connect.
"""

# ╔═╡ Cell order:
# ╟─4420bf4a-7c61-11f1-b365-4be8aa004f63
# ╠═441fd186-7c61-11f1-bb7e-a5c123e2e57e
# ╟─4420bfb0-7c61-11f1-b8e0-d1971bc7efd7
# ╠═4420bfbc-7c61-11f1-8d1a-9d469ceb89b3
# ╟─4420bfc4-7c61-11f1-b5ab-c780881c3d7f
# ╠═4420bfce-7c61-11f1-9aab-d3b02316b037
# ╟─4420bfe2-7c61-11f1-92cb-b7a95081a069
# ╠═4420bfee-7c61-11f1-b6b0-7107b3a6345d
# ╟─4420bff6-7c61-11f1-85ab-6b56a961d342
# ╠═4420c000-7c61-11f1-8269-f1fee2135a5c
# ╟─4420c00a-7c61-11f1-b7d1-81ab333e526a
# ╠═4420c014-7c61-11f1-bde2-4193fe6f111f
# ╟─4420c020-7c61-11f1-b469-bb0a725d78c9
# ╠═4420c028-7c61-11f1-b710-13f2ff6512ca
# ╟─4420c032-7c61-11f1-a731-47d2faf493ee
