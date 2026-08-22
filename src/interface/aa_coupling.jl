# AlgebraicAgents coupling surface (ADR 0012 / CONTRACT §13).
#
# A `ReactionNetworkProblem` is already an AA `@aagent` implementing the three stepping hooks
# (`_step!`/`_reinit!`/`_projected_to`, solvers.jl), so `entangle!(parent, rd)` makes it a node in
# a larger AlgebraicAgents hierarchy and AA's least-projected-time gate interleaves its single
# clock with sibling clocks for free (§C). What it did NOT have was AA's READ / COUPLING surface:
# nothing in the hierarchy could read an RD place/observable (so no AA wire could ORIGINATE from
# an RD net), the hierarchy could not read or set RD params, and there was no pinned point at
# which RD reads its incoming wires. This file adds exactly that, in BOTH directions:
#
#   (A) OUTBOUND — `getobservable`/`observables` (RD as a wire SOURCE) + `_getparameters`/
#       `_setparameters!` (param read/patch). All pure / RNG-free / index-safe.
#   (B) INBOUND  — `_prestep!`, the per-tick latch that reads RD's incoming wires ONCE, before any
#       sibling `_step!`, into `state.external_inputs` (which `ExternalRef` lowers to read).
#
# Placed in its own interface/ file (auto-included via the readdir(interface/) line in
# ReactiveDynamics.jl) so the AA method overrides stay OUT of the heavily-edited solvers.jl —
# they need nothing from it (ExternalRef lowering is self-contained, exprnode.jl). The AA method
# signatures here were verified against the installed source (AlgebraicAgents/ovDs5/src):
# `getobservable(a, args...)` (interface.jl:299), `observables(a)` (:312), `_getparameters(a)`
# (:101), `_setparameters!(a, parameters)` (:113), `_prestep!(a, _)` (:242), and the incoming-wire
# reader `retrieve_input_vars(a)` (wires.jl:26-35).

# ════════════════════════════════════════════════════════════════════════════════════════
# (A) OUTBOUND — RD as a first-class readable hierarchy node
# ════════════════════════════════════════════════════════════════════════════════════════

"""
    observables(rd::ReactionNetworkProblem)

The ordered list of names this network exports to the AlgebraicAgents hierarchy (ADR 0012 §A,
Invariant 1): every SPECIES name (`net[:,:placeName]`) followed by every NAMED observable
(`keys(state.observables)`, §9.4). This is the canonical order `getobservable(rd, i::Int)` indexes.

Token aggregates are surfaced the LEAN-EXPLICIT way the ADR open question settles on: an author
declares the aggregates worth exporting as NAMED observables (an `observables[]` entry, e.g. a
count of Phase-2 tokens), which then appear here automatically — rather than auto-enumerating a
combinatorial `nactive × kind × phase` set. Returns `Vector{Symbol}`.
"""
function AlgebraicAgents.observables(rd::ReactionNetworkProblem)
    return Symbol[collect(rd.network[:, :placeName]); collect(keys(rd.observables))]
end

"""
    getobservable(rd::ReactionNetworkProblem, name)
    getobservable(rd::ReactionNetworkProblem, i::Int)

The current value of an exported observable (ADR 0012 §A, Invariant 1). Reads are PURE and
RNG-free — they never advance `state.rng`, so a coupled read does not perturb the trajectory:

  - a SPECIES count is `state.u[idx(name)]` (the live stock, structured or classical);
  - a NAMED observable is its last-sampled `.sampled` value (§9.4);
  - an `Int` indexes `observables(rd)` (the §A canonical order).

`name` may be a `Symbol` or a `String` (AA wires carry the `from_var_name` as a string through
`retrieve_input_vars`, so both must resolve). An unknown name is a hard `error` — a diagnostic,
never AA's silent `@error` fall-through (Invariant 1).
"""
function AlgebraicAgents.getobservable(rd::ReactionNetworkProblem, name::Symbol)
    # place count (classical stock or structured `state.u`-consistent count, §9.5 observation point)
    i = find_index(name, rd)
    isnothing(i) || return rd.u[i]
    # named observable → its last-sampled value (§9.4)
    haskey(rd.observables, name) && return rd.observables[name].sampled
    return error(
        "getobservable: `$name` is not an exported observable of $(getname(rd)); " *
            "exported names are $(AlgebraicAgents.observables(rd)) (ADR 0012 §A)",
    )
end

AlgebraicAgents.getobservable(rd::ReactionNetworkProblem, name::AbstractString) =
    AlgebraicAgents.getobservable(rd, Symbol(name))

function AlgebraicAgents.getobservable(rd::ReactionNetworkProblem, i::Int)
    names = AlgebraicAgents.observables(rd)
    checkbounds(Bool, names, i) ||
        error("getobservable: index $i out of range 1:$(length(names)) (ADR 0012 §A)")
    return AlgebraicAgents.getobservable(rd, names[i])
end

"""
    _getparameters(rd::ReactionNetworkProblem)

Expose the network's parameter space `state.p` to the hierarchy (ADR 0012 §A). `state.p` is a
`Dict{Symbol,Any}`, so AA's `getparameters` walk and a coupled controller can both read it.
"""
AlgebraicAgents._getparameters(rd::ReactionNetworkProblem) = rd.p

"""
    _setparameters!(rd::ReactionNetworkProblem, parameters)

Patch the network's parameters from a `Symbol=>value` dict (ADR 0012 §A, Invariant 5). Writes are
PARAM-ONLY: they `merge!` into `state.p` and NEVER touch structure (place/transitions/the network
index), so they are ADR-0004 index-safe — structural change stays on the append-only mutation API.

Caveat (ADR 0012 open question / §5 A3): a param feeding an attribute FROZEN at spawn (a token's
`cycletime`/`prob_of_success` read from its spawn snapshot) only affects instances spawned AFTER
the set — the same mid-run-TVE caveat as the synergy levers.
"""
function AlgebraicAgents._setparameters!(rd::ReactionNetworkProblem, parameters)
    merge!(rd.p, parameters)
    return rd.p
end

# ════════════════════════════════════════════════════════════════════════════════════════
# (B) INBOUND — the `_prestep!` latch (the determinism pin, ADR 0012 §B3, Invariants 2-3)
# ════════════════════════════════════════════════════════════════════════════════════════

"""
    _prestep!(rd::ReactionNetworkProblem, t)

Latch this network's declared external inputs ONCE per tick, at a pinned point (ADR 0012 §B3).
`step!` `prewalk`s `_prestep!` over the WHOLE hierarchy in its FIRST phase, before ANY agent's
`_step!` runs that tick (AA interface.jl:184), so reading the incoming wires HERE — rather than
live mid-`_step!` — means RD always sees each source's value as projected to the PREVIOUS tick
boundary, independent of AA's `Dict`-order sibling stepping (the §4-D4 hazard this closes).

Mechanics: `retrieve_input_vars(rd)` (AA wires.jl) returns `to_var_name => getobservable(from,
from_var_name)` over RD's incoming wires. We merge that OVER a fresh copy of the declared input
DEFAULTS into `state.external_inputs`, so (i) wire reads override defaults, (ii) a defaulted-but-
unwired port keeps its default, and (iii) no value latched in a prior tick lingers. Every
`ExternalRef(port)` read during this tick's `sample_transitions!`/`evolve!`/rule-firing then
returns the SAME buffered value (Invariant 2), and the coupling is explicit/Jacobi — a one-tick
lag, no algebraic loop, reproducible under `(hierarchy, seed)` (Invariant 3).

The `to_var_name` keys arrive as the strings passed to `add_wire!`; we `Symbol`-ize them so they
match the `inputs[]` port symbols `ExternalRef` reads. RNG-free, structure-free.
"""
function AlgebraicAgents._prestep!(rd::ReactionNetworkProblem, t)
    incoming = AlgebraicAgents.retrieve_input_vars(rd)   # Dict(to_var_name => value) over in-wires
    # Re-seed from the immutable defaults, then overlay this tick's wire reads (string keys → Symbol).
    empty!(rd.external_inputs)
    merge!(rd.external_inputs, rd.external_input_defaults)
    for (k, v) in incoming
        rd.external_inputs[Symbol(k)] = v
    end
    return rd
end
