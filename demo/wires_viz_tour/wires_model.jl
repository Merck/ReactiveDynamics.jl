# ════════════════════════════════════════════════════════════════════════════════════════
# wires_model.jl — the SHARED model-building + visualization substrate for the wires_viz_tour.
#
# Both HTML forms (the Pluto notebook `wires_viz_tour.jl` and the Literate script
# `wires_literate.jl`) `include` this file and drive it — so the wired-hierarchy code, the two
# Graphviz renders, and the coupled-run data-shaping live in EXACTLY ONE place. The tutorials
# themselves are thin presentation layers (prose + a call each) over this module.
#
# Every construct here is lifted verbatim from demo/aa_integration/aa_integration.jl (the
# MarketAgent / FinanceAgent @aagents, the PHARMA_JSON model, build_portfolio) and from the
# exported visualization surface (wiring_diagram / run_graphviz from AlgebraicAgents, draw_network
# from src/visualize.jl). No API is invented.
# ════════════════════════════════════════════════════════════════════════════════════════

module WiresModel

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, from_json_model, validate
using AlgebraicAgents          # reexported by ReactiveDynamics; entangle!/add_wire!/simulate/wiring_diagram/run_graphviz
using DataFrames
import JSON

export MarketAgent, FinanceAgent, PHARMA_JSON, build_portfolio, run_coupled,
    HORIZON, wiring_diagram_svg, rd_network_svg, coupling_table, sentiment_path, all_wires

const RD = ReactiveDynamics
const HORIZON = 8.0

# ── §0. The sibling agents — a MARKET source and a FINANCE sink (verbatim from aa_integration) ──

@aagent struct MarketAgent
    sentiment::Float64   # the exported signal (deal-appetite index, 0..1)
    drift::Float64       # per-tick increment
    dt::Float64
    t::Float64
    horizon::Float64
end
MarketAgent(name::AbstractString, s0::Real, drift::Real, dt::Real, horizon::Real) =
    MarketAgent(name, Float64(s0), Float64(drift), Float64(dt), 0.0, Float64(horizon))

AlgebraicAgents.observables(m::MarketAgent) = [:sentiment]
AlgebraicAgents.getobservable(m::MarketAgent, ::Union{Symbol, AbstractString}) = m.sentiment
AlgebraicAgents.getobservable(m::MarketAgent, ::Int) = m.sentiment
AlgebraicAgents._step!(m::MarketAgent) = (m.sentiment += m.drift; m.t += m.dt; m.t)
AlgebraicAgents._projected_to(m::MarketAgent) = m.t > m.horizon ? true : m.t

@aagent struct FinanceAgent
    cash_seen::Vector{Float64}   # the wired-in RD cash, recorded each tick (the latched lag)
    dt::Float64
    t::Float64
    horizon::Float64
end
FinanceAgent(name::AbstractString, dt::Real, horizon::Real) =
    FinanceAgent(name, Float64[], Float64(dt), 0.0, Float64(horizon))

function AlgebraicAgents._prestep!(f::FinanceAgent, t)
    inputs = AlgebraicAgents.retrieve_input_vars(f)   # Dict(to_var_name => value)
    haskey(inputs, "rd_cash") && push!(f.cash_seen, Float64(inputs["rd_cash"]))
    return f
end
AlgebraicAgents._step!(f::FinanceAgent) = (f.t += f.dt; f.t)
AlgebraicAgents._projected_to(f::FinanceAgent) = f.t > f.horizon ? true : f.t

# ── §1. The RD pharma net as an eval-free JSON model (verbatim from aa_integration) ──

const PHARMA_JSON = """
{ "rd_format":"reactive-dynamics-model", "version":"1.0",
  "meta":{ "tspan":8.0, "dt":1.0 },
  "params":[ { "name":"base_inflow", "value":30.0 },
             { "name":"acq_cash_trigger", "value":50.0 },
             { "name":"sentiment_threshold", "value":0.45 } ],
  "inputs":[ { "port":"sentiment", "default":{ "node":"const", "value":0.0 } } ],
  "species":[ { "name":"cash", "init":0 },
              { "name":"acquired", "init":0 } ],
  "transitions":[
    { "id":"grow", "name":"grow", "rate_mode":"deterministic",
      "rate":{ "node":"call", "op":"*",
               "args":[ { "node":"externalref", "port":"sentiment" },
                        { "node":"ref", "kind":"param", "name":"base_inflow" } ] } } ],
  "reactants":[
    { "transition":"grow", "side":"rhs", "species":"cash", "stoich":1 } ],
  "rules":[
    { "id":"acquire", "fire_mode":"once",
      "guard":{ "node":"call", "op":"&&",
                "args":[ { "node":"call", "op":">",
                           "args":[ { "node":"externalref", "port":"sentiment" },
                                    { "node":"ref", "kind":"param", "name":"sentiment_threshold" } ] },
                         { "node":"call", "op":">=",
                           "args":[ { "node":"ref", "kind":"species", "name":"cash" },
                                    { "node":"ref", "kind":"param", "name":"acq_cash_trigger" } ] } ] },
      "action":{ "verb":"seq", "stmts":[
                   { "verb":"set_species", "name":"acquired", "mode":"inc",
                     "value":{ "node":"const", "value":1 } },
                   { "verb":"log", "msg":"ACQUISITION: external sentiment cleared + cash trigger met" } ] } } ] }
"""

# ── §2. Compose the hierarchy & WIRE it — both directions, host-side (verbatim from aa_integration) ──

function build_portfolio(; seed)
    rd = from_json_model(PHARMA_JSON; seed = seed)
    market = MarketAgent("market", 0.0, 0.12, 1.0, HORIZON)   # sentiment 0.0, +0.12/tick
    finance = FinanceAgent("finance", 1.0, HORIZON)
    root = FreeAgent("portfolio")
    entangle!(root, market)
    entangle!(root, rd)
    entangle!(root, finance)
    # wire 1 — INBOUND: the market's sentiment feeds RD's `sentiment` port.
    add_wire!(root; from = market, to = rd, from_var_name = "sentiment", to_var_name = "sentiment")
    # wire 2 — OUTBOUND: RD's cash (read via getobservable) feeds the finance agent.
    add_wire!(root; from = rd, to = finance, from_var_name = "cash", to_var_name = "rd_cash")
    return (; root, rd, market, finance)
end

# Every wire in the hierarchy, as (from, from_var, to, to_var) name-tuples for a tidy printout. The
# wires live on the shared Opera (getopera(root).wires); each agent's get_wires_from/get_wires_to
# is the per-agent filtered view of that same list.
function all_wires(root)
    ws = AlgebraicAgents.getopera(root).wires
    return [
        (
                from = AlgebraicAgents.getname(w.from), from_var = String(w.from_var_name),
                to = AlgebraicAgents.getname(w.to), to_var = String(w.to_var_name),
            ) for w in ws
    ]
end

# Build + run the whole coupled system under a single simulate(root).
function run_coupled(; seed = 11)
    sys = build_portfolio(; seed)
    simulate(sys.root)
    return sys
end

# ── Visualization helpers — the two Graphviz renders, each returning a self-contained SVG STRING ──

# Strip the XML/DOCTYPE preamble Graphviz emits so the `<svg>` can be inlined directly into an HTML
# document (safe for both Pluto's HTML() and the Literate→HTML inliner), then scale it up for
# on-screen readability. Graphviz encodes size as `width="Xpt" height="Ypt"` alongside a `viewBox`;
# multiplying the pt dimensions (viewBox untouched) makes the browser render it larger with no
# quality loss — pure presentation, the vector content is unchanged.
function _clean_svg(s::AbstractString; scale::Real = 1.0)
    i = findfirst("<svg", s)
    svg = i === nothing ? String(s) : String(s[first(i):end])
    scale == 1.0 && return svg
    return replace(
        svg, r"width=\"([\d.]+)pt\" height=\"([\d.]+)pt\"" =>
            m -> begin
            mm = match(r"width=\"([\d.]+)pt\" height=\"([\d.]+)pt\"", m)
            w = round(parse(Float64, mm.captures[1]) * scale; digits = 1)
            h = round(parse(Float64, mm.captures[2]) * scale; digits = 1)
            "width=\"$(w)pt\" height=\"$(h)pt\""
        end; count = 1
    )
end

# (2) The HEADLINE visual: AA's wiring_diagram(root) → DOT → run_graphviz → SVG. This is the
# cross-agent view — the three agents, the parentship edges, and the two LABELED wires.
function wiring_diagram_svg(root; scale = 2.6)
    dot = AlgebraicAgents.wiring_diagram(root)         # returns a Graphviz DOT string
    io = IOBuffer()
    AlgebraicAgents.run_graphviz(io, dot; prog = :dot, format = "svg")
    return _clean_svg(String(take!(io)); scale = scale)
end

# (3) The intra-agent view: RD's own draw_network of the pharma net INSIDE the RD box. With no
# `path`, draw_network returns the rendered SVG bytes as a String (src/visualize.jl).
rd_network_svg(rd; scale = 2.6) = _clean_svg(draw_network(rd; format = "svg"); scale = scale)

# ── Data-shaping — the coupled trajectory as one aligned DataFrame (the one-tick Jacobi lag) ──

# The market's sentiment at each RD tick boundary. The market starts at s0=0.0 and drifts up by
# `drift` each tick (MarketAgent._step!), so the boundary value at tick k is s0 + k*drift — the
# per-tick source series the inbound wire carries into RD (one tick before RD reads it).
function sentiment_path(sys)
    n = length(sys.rd.sol.t)
    s0, drift = 0.0, sys.market.drift
    return round.([s0 + k * drift for k in 0:(n - 1)]; digits = 3)
end

# One tidy table aligned on the RD clock: sentiment(t), RD cash(t), the finance agent's
# reconstructed cash (RD's cash read over the outbound wire), and the acquisition lever.
function coupling_table(sys)
    t = Int.(sys.rd.sol.t)
    n = length(t)
    cash = round.(sys.rd.sol.cash; digits = 1)
    acquired = Int.(sys.rd.sol.acquired)
    sent = sentiment_path(sys)
    # finance.cash_seen is latched once per tick at _prestep!; align it to the RD clock, padding
    # with `missing` if the lengths differ (the lag can leave it one shorter/longer).
    seen = sys.finance.cash_seen
    finance_seen = Union{Float64, Missing}[
        k <= length(seen) ? round(seen[k]; digits = 1) : missing
            for k in 1:n
    ]
    return DataFrame(
        t = t,
        sentiment = sent,
        rd_cash = cash,
        finance_cash_seen = finance_seen,
        acquired = acquired,
    )
end

end # module WiresModel
