# ════════════════════════════════════════════════════════════════════════════════════════
# ReactiveDynamics.jl — a literate tour of AlgebraicAgents INTEGRATION & external coupling
# ════════════════════════════════════════════════════════════════════════════════════════
#
# Run it:   julia --project=. demo/aa_integration/aa_integration.jl
#
# This script is an essay-with-code (the same literate style as demo/core_engine_tour and
# demo/agentic_pipeline). It demonstrates ADR 0012 / CONTRACT §13: a ReactiveDynamics reaction
# network sitting as ONE node inside a larger AlgebraicAgents (AA) hierarchy, COUPLED to sibling
# agents in BOTH directions, and driven by a single `simulate(root)`.
#
# The story (the ADR North-star): a pharma R&D portfolio (the RD net) lives beside a MACRO/MARKET
# agent and a downstream FINANCE agent under one portfolio-level root:
#
#                        ┌─────────────── root (FreeAgent) ───────────────┐
#                        │                                                │
#              market ───┤  sentiment ─wire─▶ [ RD pharma net ]           │
#            (a sibling  │                          │ getobservable       │
#             @aagent)   │                          ▼                     │
#                        │            cash ─wire─▶ finance (downstream)    │
#                        └────────────────────────────────────────────────┘
#
#   (1) INBOUND  — the RD net READS the market agent's `sentiment` through a declared `inputs[]`
#       port + an `ExternalRef` leaf, used BOTH in a transition RATE and in an acquisition RULE
#       GUARD (the lever fires only once external sentiment clears a threshold). This is the
#       endogenous decision channel (ADR 0010) driven by EXTERNAL state.
#   (2) OUTBOUND — a downstream FINANCE agent READS the RD net's `cash` via AA `getobservable`
#       carried on a wire, so RD is a wire SOURCE as well as a sink.
#
# The load-bearing guarantees this demo makes executable:
#   • RD is a first-class hierarchy node: AA's least-projected-time gate interleaves RD's clock
#     with the sibling clocks for free — no new clock code (§C).
#   • The coupling is EXPLICIT (Jacobi): external reads are latched ONCE per tick at `_prestep!`,
#     before any sibling `_step!`, so RD always reads each source's PREVIOUS-tick-boundary value
#     — a one-tick lag, no algebraic loop, and fully reproducible under (hierarchy, seed). That
#     `_prestep!` latch is the determinism pin (Invariants 2-3; it closes the §4-D4 sibling-order
#     hazard a live mid-step cross-agent read would open).
#   • The RD model document stays self-contained: it declares the PORTS it reads; the foreign-
#     agent topology (which agent feeds which port) lives host-side in `add_wire!`, never in the
#     JSON (Invariant 4, eval-free coupling).
#
# Everything below uses only constructs from the engine's PASSING semantic test suite
# (test/semantic/{aa_integration,serialization_ir,rules_decisions}.jl). No invented API.

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, from_json_model, validate, find_index, ExternalRef
using AlgebraicAgents          # reexported by ReactiveDynamics; entangle!/add_wire!/simulate/step!
using Random, DataFrames
import JSON

const RD = ReactiveDynamics

banner(title) = (println(); println("="^88); println(title); println("="^88))

# ════════════════════════════════════════════════════════════════════════════════════════
# §0. The sibling agents — a MARKET source and a FINANCE sink (minimal AA @aagents)
# ════════════════════════════════════════════════════════════════════════════════════════
#
# To couple the RD net to "the rest of the hierarchy" we need at least one other agent. We define
# two tiny ones in the host (NOT in RD's scope — these are ordinary AA agents):
#
#   MarketAgent — a SOURCE. It carries a `sentiment` signal that drifts upward each tick (a stand-
#     in for a macro/market model: rising deal appetite). It EXPORTS `sentiment` as an observable
#     so a wire can carry it into the RD net.
#   FinanceAgent — a SINK. Each tick it READS a wired-in value (the RD net's `cash`) and records
#     the running maximum it has seen. It demonstrates RD as a wire SOURCE: the finance agent
#     consumes RD's `getobservable(:cash)` with no knowledge of RD's internals.
#
# Both implement only the AA surface they need. NB the custom constructor types its first arg
# `::AbstractString` so dispatch picks it over the @aagent-generated positional `T(name, args...)`.

@aagent struct MarketAgent
    sentiment::Float64   # the exported signal (deal-appetite index, 0..1)
    drift::Float64       # per-tick increment
    dt::Float64
    t::Float64
    horizon::Float64
end
MarketAgent(name::AbstractString, s0::Real, drift::Real, dt::Real, horizon::Real) =
    MarketAgent(name, Float64(s0), Float64(drift), Float64(dt), 0.0, Float64(horizon))

# OUTBOUND surface for the market agent: it exports one observable, `sentiment`.
AlgebraicAgents.observables(m::MarketAgent) = [:sentiment]
AlgebraicAgents.getobservable(m::MarketAgent, ::Union{Symbol, AbstractString}) = m.sentiment
AlgebraicAgents.getobservable(m::MarketAgent, ::Int) = m.sentiment
# Stepping: drift the sentiment upward, advance the clock.
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

# INBOUND for the finance agent: each tick latch the wired-in RD cash. We read the agent's own
# incoming wires the SAME way RD does — `retrieve_input_vars` (AA), at `_prestep!` — so the
# finance read carries the same one-tick lag and the same determinism guarantee as the RD read.
function AlgebraicAgents._prestep!(f::FinanceAgent, t)
    inputs = AlgebraicAgents.retrieve_input_vars(f)   # Dict(to_var_name => value)
    haskey(inputs, "rd_cash") && push!(f.cash_seen, Float64(inputs["rd_cash"]))
    return f
end
AlgebraicAgents._step!(f::FinanceAgent) = (f.t += f.dt; f.t)
AlgebraicAgents._projected_to(f::FinanceAgent) = f.t > f.horizon ? true : f.t

banner("§0. The sibling agents (a MARKET source + a FINANCE sink)")
println("MarketAgent  — exports observable :sentiment (a drifting deal-appetite index).")
println("FinanceAgent — reads RD's :cash through a wire at _prestep! and records what it sees.")
println("Both are ordinary AA @aagents in the host; neither knows anything about RD's internals.")

# ════════════════════════════════════════════════════════════════════════════════════════
# §1. The RD pharma net as an eval-free JSON model — declaring its external INPUT port
# ════════════════════════════════════════════════════════════════════════════════════════
#
# The reactive network is a small portfolio cash model authored as a JSON document (ADR 0005). The
# NEW piece (ADR 0012 §B1) is the top-level `inputs[]` array: it declares the OBSERVABLE-level read
# ports the network consumes — here a single port `sentiment`, with a pre-wire `default` of 0.0
# (the value a read returns before any wire has delivered, or in a standalone run with no wires).
#
# The model reads that port through the `ExternalRef` leaf in TWO places:
#   • a transition RATE — the `grow` line accrues cash at a rate proportional to external
#     sentiment (`ExternalRef(sentiment) * base`), so a hotter market funds the portfolio faster;
#   • a RULE GUARD — the `acquire` lever (a once-rule) fires only when `ExternalRef(sentiment)`
#     clears a threshold AND the portfolio has built enough cash. That is the ADR North-star:
#     an endogenous decision driven by BOTH internal portfolio state and EXTERNAL market state.
#
# The JSON names only the PORT. Which agent feeds it is decided host-side in §2 via `add_wire!` —
# the document is portable and carries no foreign-agent topology (Invariant 4).

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

banner("§1. The RD pharma net as JSON — its declared external INPUT port")
# A pure, eval-free self-check BEFORE building (validate rule 8: every ExternalRef.port is a
# declared inputs[] port). [] means clean.
diags = validate(JSON.parse(PHARMA_JSON))
println("validate(model) -> ", isempty(diags) ? "OK (no diagnostics)" : diags)
println("Declared input port: :sentiment (default 0.0). It is read by an ExternalRef in BOTH")
println("the `grow` transition rate AND the `acquire` rule guard — a rate and a decision driven")
println("by external state. The JSON names only the PORT; the wiring is host-side (§2).")

# A demonstration of validate rule 8: an UNDECLARED port is a diagnostic, never an eval.
bad = JSON.parse(PHARMA_JSON)
bad["transitions"][1]["rate"]["args"][1]["port"] = "undeclared_signal"
bad_diags = validate(bad)
println("validate(model with an UNDECLARED port) -> ", length(bad_diags), " diagnostic(s):")
for d in bad_diags
    println("    ", string(d))
end

# ════════════════════════════════════════════════════════════════════════════════════════
# §2. Compose the hierarchy & WIRE it — both directions, host-side
# ════════════════════════════════════════════════════════════════════════════════════════
#
# We entangle the three agents under one portfolio root and lay TWO wires. A wire
# `add_wire!(root; from, to, from_var_name, to_var_name)` says "carry the source's observable
# `from_var_name` into the target's input `to_var_name`." This is the ONLY place the cross-agent
# topology lives — keeping it out of the inert RD document (Invariant 4).
#
#   wire 1 (INBOUND to RD) : market.sentiment ─▶ RD.sentiment port
#   wire 2 (OUTBOUND from RD): RD.cash ─▶ finance.rd_cash
#
# `getobservable` on the RD net is what makes wire 2 possible: the finance agent reads RD's `cash`
# with no knowledge of how RD computes it.

const HORIZON = 8.0

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

banner("§2. Compose & wire the hierarchy (both directions)")
sys = build_portfolio(seed = 11)
println("Hierarchy: portfolio (root) ⊇ {market, RD pharma net, finance}.")
println("RD exports observables: ", AlgebraicAgents.observables(sys.rd))
println("Wire 1 (inbound) : market.sentiment ──▶ RD.sentiment port")
println("Wire 2 (outbound): RD.cash (getobservable) ──▶ finance.rd_cash")
println(
    "RD reads sentiment=", AlgebraicAgents.getobservable(sys.market, :sentiment),
    " at t=0 (the pre-wire default in RD's buffer is ", sys.rd.external_inputs[:sentiment], ")."
)

# ════════════════════════════════════════════════════════════════════════════════════════
# §3. Run the WHOLE coupled model under one simulate(root)
# ════════════════════════════════════════════════════════════════════════════════════════
#
# A coupled RD net is driven by `simulate` on the ROOT, not `simulate(rd)` directly (§C). AA's
# `step!` (a) `prewalk`s every agent's `_prestep!` FIRST (the latch phase — this is where RD and
# the finance agent read their incoming wires, once, for the whole tick), then (b) steps each agent
# whose projected time equals the hierarchy minimum. So every external read in a tick is the
# source's PREVIOUS-tick-boundary value (the one-tick Jacobi lag).

banner("§3. Run the coupled model — one simulate(root)")
simulate(sys.root)

sentiment_path = round.(0.0:0.12:(0.12 * 8); digits = 3)
println("market sentiment over the run (drifts up):   ", sentiment_path[1:9])
println()
println("RD cash trajectory (grow rate = sentiment × base_inflow, latched one tick late):")
println("   t    : ", Int.(sys.rd.sol.t))
println("   cash : ", round.(sys.rd.sol.cash; digits = 1))
println("   acquired (the lever): ", Int.(sys.rd.sol.acquired))
println()
acq_tick = findfirst(>(0.0), sys.rd.sol.acquired)
println(
    "The acquisition lever fired at t = ", isnothing(acq_tick) ? "never" : Int(sys.rd.sol.t[acq_tick]),
    " — the first tick BOTH the external sentiment cleared its threshold AND cash ≥ trigger."
)
println("rule `acquire` enabled? ", sys.rd.rules[1].enabled, "  (false ⇒ the :once lever has fired)")

# ════════════════════════════════════════════════════════════════════════════════════════
# §4. The OUTBOUND read — the finance agent saw RD's cash through a wire
# ════════════════════════════════════════════════════════════════════════════════════════
#
# The finance agent never touched RD; it only read `rd_cash` off its incoming wire each `_prestep!`.
# Because that read is pinned at the latch (the source's previous boundary), the values it recorded
# are RD's cash trajectory shifted by the one-tick coupling lag — the explicit-coupling signature.

banner("§4. The downstream OUTBOUND read (finance read RD's cash via getobservable)")
println("finance.cash_seen (RD cash as the finance agent latched it each tick):")
println("   ", round.(sys.finance.cash_seen; digits = 1))
println("RD cash (the source):")
println("   ", round.(sys.rd.sol.cash; digits = 1))
println("The finance series is RD's cash carried over wire 2 — RD is a wire SOURCE, read purely")
println("through getobservable, with the same one-tick lag the latch guarantees.")

# ════════════════════════════════════════════════════════════════════════════════════════
# §5. Determinism — the coupled trajectory is reproducible under (hierarchy, seed)
# ════════════════════════════════════════════════════════════════════════════════════════
#
# The whole point of pinning external reads at `_prestep!` is that the coupled run is a function of
# (hierarchy, seed) ALONE — independent of AA's Dict-order sibling stepping (the §4-D4 hazard a
# live mid-step cross-agent read would expose). We build and run the SAME system twice and compare.

banner("§5. Determinism — same (hierarchy, seed) ⇒ identical coupled trajectory")
a = build_portfolio(seed = 11); simulate(a.root)
b = build_portfolio(seed = 11); simulate(b.root)
println("two independent coupled runs, same seed:")
println("   RD trajectories identical (sol == sol)? ", a.rd.sol == b.rd.sol)
println("   finance reads identical?                ", a.finance.cash_seen == b.finance.cash_seen)
# a DIFFERENT seed leaves THIS deterministic model identical (no RNG draws here), but the machinery
# is seed-threaded: any stochastic rate/lever would diverge. We assert reproducibility, the contract.
println("The latch makes the coupling explicit (Jacobi, one-tick lag) and fully reproducible —")
println("no algebraic loop, no sibling-order dependence (Invariants 2-3).")

# ════════════════════════════════════════════════════════════════════════════════════════
# §6. reinit — a hierarchy-wide reset clears the external buffer back to its declared default
# ════════════════════════════════════════════════════════════════════════════════════════
#
# `_reinit!` on the RD node restores its run-state AND re-seeds `external_inputs` from the declared
# input defaults (§B3 / §4 D7), so the latched wire values from the finished run do not leak into
# the next one — the first post-reinit `_prestep!` re-latches from a clean seed.

banner("§6. reinit clears the external-input buffer to its declared default")
println(
    "after the run, RD's latched sentiment = ", round(sys.rd.external_inputs[:sentiment]; digits = 3),
    " (a stale wire value)"
)
AlgebraicAgents._reinit!(sys.rd)
println(
    "after _reinit!, RD's sentiment buffer  = ", sys.rd.external_inputs[:sentiment],
    "  (restored to the declared inputs[] default)"
)
println("equal to the declared defaults snapshot? ", sys.rd.external_inputs == sys.rd.external_input_defaults)

# ════════════════════════════════════════════════════════════════════════════════════════
# §7. Recap
# ════════════════════════════════════════════════════════════════════════════════════════
banner("§7. Recap")
println(
    """
    We ran a pharma portfolio RD net as ONE node inside an AlgebraicAgents hierarchy, coupled in BOTH
    directions and driven by a single simulate(root):

      §0  sibling agents      a MarketAgent (source) + a FinanceAgent (sink), ordinary AA @aagents
      §1  inputs[] port       the RD JSON declares a `sentiment` read port; ExternalRef reads it in a
                              RATE and a RULE GUARD; validate rule 8 enforces the port is declared
      §2  wiring (host-side)  add_wire! lays market.sentiment ▶ RD, and RD.cash ▶ finance — topology
                              lives in the host, never in the RD document (Invariant 4)
      §3  simulate(root)      AA's least-projected-time gate interleaves the clocks; _prestep! latches
                              external reads once/tick; the acquisition lever fires on external state
      §4  OUTBOUND read       finance consumed RD's cash purely via getobservable on a wire
      §5  determinism         same (hierarchy, seed) ⇒ identical coupled trajectory (Invariants 2-3)
      §6  reinit              the external buffer resets to its declared default (§B3 / §4 D7)

    The through-line (ADR 0012 / CONTRACT §13): a reactive network is a first-class AA hierarchy node
    in BOTH directions — readable by the hierarchy (getobservable/observables) and able to read it
    (inputs[] + ExternalRef) — and the coupling is explicit, eval-free, and deterministic: one new
    closed ExternalRef leaf, a model-local inputs[] port list, and a pinned _prestep! latch.
    """
)
println("Done.")
