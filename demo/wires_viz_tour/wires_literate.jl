# # ReactiveDynamics as an AlgebraicAgents node — the wires, made visible
#
# A ReactiveDynamics reaction network is a first-class node in a larger **AlgebraicAgents (AA)**
# hierarchy. Coupled to sibling agents through **wires** — open ports carrying one agent's
# observable into another's input — it runs under a single `simulate(root)`. This is ADR 0012 /
# CONTRACT §13.
#
# The [`aa_integration`](../aa_integration) demo teaches this coupling in prose and `println`. This
# tour makes it **visual**: it draws the wired hierarchy as a diagram and the RD net as a diagram,
# side by side with the coupled-run data. The subject is the *ports and wires you could not see
# before* — and *what lives inside the RD box they connect*.
#
# The system is a small **pharma R&D portfolio** under one root: a **market** agent that exports a
# drifting `sentiment` signal; the **RD net** with an `inputs[]` `sentiment` port, growing cash at a
# rate driven by that external sentiment and firing an acquisition lever when both external sentiment
# and internal cash clear their thresholds; and a **finance** agent that reads the RD net's `cash`
# back off a wire. Two wires close the loop: `market.sentiment ▶ RD.sentiment` (inbound) and
# `RD.cash ▶ finance.rd_cash` (outbound).

# We load the shared substrate — the wired-hierarchy build, the two Graphviz renders, and the
# data-shaping — lifted verbatim from `demo/aa_integration`. This script is a thin presentation
# layer over it (the SAME `wires_model.jl` the Pluto notebook uses). `build.jl` sets
# `WIRES_VIZ_DIR` because Literate executes with the working directory changed to its output folder;
# when the script is run directly, `@__DIR__` is the right fallback.
const WIRES_VIZ_DIR = get(ENV, "WIRES_VIZ_DIR", @__DIR__)
include(joinpath(WIRES_VIZ_DIR, "wires_model.jl"))
using .WiresModel
using ReactiveDynamics, AlgebraicAgents, DataFrames, Plots

# A tiny helper so a raw SVG string renders as an image in the generated page: an object whose
# `show(::MIME"image/svg+xml")` prints the SVG. Literate writes it out as a sidecar `.svg` and the
# HTML build step inlines it — so the network diagrams appear inline, not as code.
struct RawSVG
    s::String
end
Base.show(io::IO, ::MIME"image/svg+xml", x::RawSVG) = print(io, x.s)

# ## 1 · Build the wired hierarchy
#
# `build_portfolio` entangles the three agents under one `FreeAgent("portfolio")` root and lays the
# two wires with `add_wire!(root; from, to, from_var_name, to_var_name)` — the *only* place the
# cross-agent topology lives (the RD JSON names only its `sentiment` port, never the market agent
# that feeds it: Invariant 4, eval-free coupling). We build and run the whole thing with one
# `simulate(root)`, wrapped as `run_coupled`.

sys = run_coupled(seed = 11)

# RD exports these observables (the outbound-wire sources a sibling can read):

AlgebraicAgents.observables(sys.rd)

# And here are the two wires, read straight off the shared `Opera` — `from` agent / `from_var`
# observable ▶ `to` agent / `to_var` input:

DataFrame(all_wires(sys.root))

# ## 2 · Visualize the wiring — the headline diagram
#
# `AlgebraicAgents.wiring_diagram(root)` returns a Graphviz **DOT** string of the whole hierarchy —
# agent nodes, parentship edges (dashed), and the annotated **wires** — which `run_graphviz` renders
# to SVG (the same backend RD's `draw_network` uses). This is the cross-agent view: *who is wired to
# whom*.
#
# Read it: `portfolio` (node 1) parents `market`, the RD net (`reaction_network`, which itself
# parents its `structured` token store), and `finance`. The two solid arrows are the wires —
# `market → reaction_network` labelled `sentiment` at both ends (the inbound port), and
# `reaction_network → finance` labelled `cash` at the tail / `rd_cash` at the head (the outbound
# wire: tail = the source observable, head = the target input).

RawSVG(wiring_diagram_svg(sys.root))

# ## 3 · Visualize the RD node's internals
#
# The wiring diagram shows the RD net as one box. `draw_network(rd)` opens that box: it renders the
# pharma net's **internal Petri structure** — the `cash` and `acquired` places (circles), the `grow`
# transition (box), and the arc `grow → cash`. Together with the wiring diagram above, this is the
# full coupling picture at two granularities: the cross-agent wiring, and the intra-agent structure
# the `sentiment` port feeds into.

RawSVG(rd_network_svg(sys.rd))

# ## 4 · Run the coupled model + show the data
#
# `simulate(root)` already ran inside `run_coupled`. AA's least-projected-time gate interleaves the
# clocks; every tick, `_prestep!` latches each agent's incoming wires **once** — so RD reads the
# market's *previous* tick-boundary sentiment (an explicit, one-tick **Jacobi lag**, no algebraic
# loop). The table aligns, per RD tick: the market `sentiment`, RD `cash`, the finance agent's
# reconstructed `cash` (RD's cash read purely over the outbound wire), and the `acquired` lever.

tbl = coupling_table(sys)

# The **one-tick lag is visible in the data**: sentiment rises to `0.12` at `t=1`, but RD `cash`
# stays `0` until `t=2` — because the sentiment RD latched at `t=1` was the pre-wire default from the
# `t=0` boundary. And `finance_cash_seen` tracks `rd_cash` exactly: the finance agent reconstructed
# RD's trajectory knowing nothing of RD's internals, purely through `getobservable` on wire 2. The
# chart shows all three series and marks the tick the acquisition lever fires:

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
    current()
end

# ## 5 · Recap — the coupling surface (ADR 0012 / §13)
#
# - **RD is a first-class hierarchy node.** `entangle!` makes it a child; AA's least-projected-time
#   gate interleaves its single clock with the siblings' for free (no new clock code). A coupled net
#   is driven by `simulate(root)`, never `simulate(rd)`.
# - **Inbound coupling** is a declarative, eval-free port: the RD JSON's `inputs[]` list plus a
#   closed `ExternalRef` leaf that a rate or a guard reads by name. Here `sentiment` drives both the
#   `grow` rate and the `acquire` guard.
# - **Outbound coupling** is `getobservable`: any sibling can read an RD species or named observable
#   off a wire, as `finance` reads `cash`.
# - **The wiring lives host-side** in `add_wire!` — never in the portable RD document (Invariant 4).
# - **The coupling is explicit and deterministic**: external reads are latched once per tick at
#   `_prestep!` (the one-tick Jacobi lag), so a coupled run is a function of `(hierarchy, seed)` alone.
#
# The two diagrams above are the two halves of that story: `wiring_diagram` shows the ports and wires
# *between* agents; `draw_network` shows the Petri structure *inside* the RD node they connect.
