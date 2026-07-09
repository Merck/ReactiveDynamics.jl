module ReactiveDynamics

using Reexport
using MacroTools
using ComponentArrays

@reexport using GeneratedExpressions

# ADR 0003 Phase 1: the static authoring/IR store is a dependency-free typed-struct-of-columns
# (no ACSets). The store type keeps the name `ReactionNetworkSchema` and a public `.subparts`
# NamedTuple of typed columns, so the ~70 indexing sites and the 8 `propertynames(acs.subparts)`
# reflection loops compile unchanged. These verbs are the signature-preserving shim that replaces
# the ACSets API surface RD used; they are exported because callers (incl. tests) used them bare
# from ACSets before. See the SCHEMA + shim block below.
export nparts, parts, dom_parts, incident, subpart, set_subpart!, add_part!, add_parts!, rem_parts!

const SampleableValues = Union{Expr,Symbol,AbstractString,Float64,Int,Function}
const ActionableValues = Union{Function,Symbol,Float64,Int}

const SampleableRange = Union{
    Float64,
    Int64,
    AbstractString,
    Expr,
    Symbol,
    Tuple{Float64,Union{Float64,Int64,AbstractString,Expr,Symbol}},
}

Base.convert(::Type{SampleableRange}, x::Tuple) = (Float64(x[1]), x[2])

Base.@kwdef mutable struct FoldedObservable
    range::Vector{SampleableRange} = SampleableRange[]
    every::Float64 = Inf
    on::Vector{SampleableValues} = SampleableValues[]
end

# ── ADR 0003 Phase 1: the typed-struct-of-columns static store ────────────────────────────────
#
# `const SCHEMA` is the single source of truth for the object model — six objects (:S species,
# :T transitions, :E events, :obs observables, :P params, :M meta), ZERO homs — replacing the old
# ACSets `BasicSchema`. Each object maps to a NamedTuple of `column ⇒ element-type`. The declaration
# ORDER (S-cols, then T, E, obs, P, M) is load-bearing: the store's `.subparts` NamedTuple is built
# in this order, so `propertynames(acs.subparts)` reproduces the exact ACSets column order the eight
# reflection loops (compilers.jl, solvers.jl, joins.jl, equalize.jl) filter over by substring.
const SCHEMA = (
    S = (
        specName = Symbol,
        specModality = Set{Symbol},
        specInitVal = SampleableValues,
        specInitUncertainty = SampleableValues,
        specCost = SampleableValues,
        specReward = SampleableValues,
        specValuation = SampleableValues,
        specStructured = Bool,
    ),
    T = (
        trans = SampleableValues,
        transPriority = SampleableValues,
        transRate = SampleableValues,
        transCycleTime = SampleableValues,
        transProbOfSuccess = SampleableValues,
        transCapacity = SampleableValues,
        transMaxLifeTime = SampleableValues,
        transPreAction = SampleableValues,
        transPostAction = SampleableValues,
        transMultiplier = SampleableValues,
        transName = Union{String,Symbol,Missing},
    ),
    E = (eventTrigger = SampleableValues, eventAction = SampleableValues),
    obs = (obsName = Symbol, obsOpts = FoldedObservable),
    P = (prmName = Symbol, prmVal = Any),
    M = (metaKeyword = Symbol, metaVal = SampleableValues),
)

# attr → owning object, and the flat ordered attr list (matches ACSets `propertynames(subparts)`).
const ATTR2OBJ = Dict{Symbol,Symbol}(
    a => obj for obj in keys(SCHEMA) for a in keys(SCHEMA[obj])
)
const ALLATTRS = Tuple(a for obj in keys(SCHEMA) for a in keys(SCHEMA[obj]))

# `columns(SCHEMA)` = all attrs (ordered); `columns(SCHEMA, obj)` = that object's attrs. These
# replace the `propertynames(acs.subparts)` reflection where a schema-driven list is wanted; the
# in-place `.subparts` loops keep using `propertynames` since `.subparts` IS a NamedTuple.
columns(::typeof(SCHEMA)) = ALLATTRS
columns(::typeof(SCHEMA), obj::Symbol) = keys(SCHEMA[obj])

# A single typed column: values plus a `defined` bitmap. An unassigned cell reads back as `nothing`
# (cloning ACSets' `get(col, i, default=nothing)`), so the `isnothing(acs[i,k])` guards in
# assign_defaults!/solvers keep working, and `grow!` need not fabricate a typed default value.
mutable struct AttrColumn{T}
    v::Vector{T}
    def::Vector{Bool}
end
AttrColumn{T}() where {T} = AttrColumn{T}(T[], Bool[])

@inline getcell(c::AttrColumn, i::Int) = @inbounds(c.def[i]) ? (@inbounds c.v[i]) : nothing
# Assigning into a `Vector{T}` invokes `convert(T, x)` — this is what preserves the String→Symbol
# (and String→SampleableValues, parse-not-eval) coercion the DSL/loader relied on (hooks below).
@inline setcell!(c::AttrColumn{T}, i::Int, x) where {T} =
    (@inbounds c.v[i] = x; @inbounds c.def[i] = true; x)
@inline grow!(c::AttrColumn{T}) where {T} = (resize!(c.v, length(c.v) + 1); push!(c.def, false))

# The static network container. Keeps the name `ReactionNetworkSchema` so every existing signature
# and `state.acs::ReactionNetworkSchema` annotation compiles unchanged. `parts` counts rows per
# object; `subparts` is the NamedTuple of typed columns in ALLATTRS order.
struct ReactionNetworkSchema
    parts::Dict{Symbol,Int}
    subparts::NamedTuple
end

function ReactionNetworkSchema()
    parts = Dict{Symbol,Int}(obj => 0 for obj in keys(SCHEMA))
    cols = NamedTuple{ALLATTRS}(AttrColumn{SCHEMA[ATTR2OBJ[a]][a]}() for a in ALLATTRS)
    return ReactionNetworkSchema(parts, cols)
end

# ACSets hashed a static model by CONTENT (so export.jl `_model_hash = hash(prob.acs)` names a
# stable bundle dir); a struct's default hash is object-identity. Preserve content-hashing.
function Base.hash(acs::ReactionNetworkSchema, h::UInt)
    h = hash(:ReactionNetworkSchema, h)
    for a in ALLATTRS
        col = acs.subparts[a]
        for i = 1:acs.parts[ATTR2OBJ[a]]
            h = hash(getcell(col, i), h)
        end
    end
    return h
end

# ── the signature-preserving shim (replaces the ACSets API surface RD used) ───────────────────
@inline _col(acs::ReactionNetworkSchema, attr::Symbol) = getfield(acs, :subparts)[attr]

nparts(acs::ReactionNetworkSchema, obj::Symbol) = acs.parts[obj]
parts(acs::ReactionNetworkSchema, obj::Symbol) = Base.OneTo(acs.parts[obj])
dom_parts(acs::ReactionNetworkSchema, attr::Symbol) = Base.OneTo(acs.parts[ATTR2OBJ[attr]])

# scalar get/set
Base.getindex(acs::ReactionNetworkSchema, i::Int, attr::Symbol) = getcell(_col(acs, attr), i)
Base.setindex!(acs::ReactionNetworkSchema, v, i::Int, attr::Symbol) = setcell!(_col(acs, attr), i, v)
# whole-column and row-subset reads return COPIES (matching ACSets `collect_column`/`map(identity)`);
# `map(identity, …)` narrows e.g. `acs[:, :specName]` back to `Vector{Symbol}`.
Base.getindex(acs::ReactionNetworkSchema, ::Colon, attr::Symbol) =
    map(identity, [getcell(_col(acs, attr), i) for i = 1:acs.parts[ATTR2OBJ[attr]]])
Base.getindex(acs::ReactionNetworkSchema, rows::AbstractVector, attr::Symbol) =
    [getcell(_col(acs, attr), i) for i in rows]

subpart(acs::ReactionNetworkSchema, attr::Symbol) = acs[:, attr]           # COPY (see note)
subpart(acs::ReactionNetworkSchema, i::Int, attr::Symbol) = acs[i, attr]
set_subpart!(acs::ReactionNetworkSchema, i::Int, attr::Symbol, v) = (acs[i, attr] = v)

# `incident(acs, val, attr)` = findall over the attr column (never an FK follow — no homs). `isequal`
# (not `==`) so `nothing`/`missing` cells compare `false`, never poison the result with `missing`.
incident(acs::ReactionNetworkSchema, val, attr::Symbol) =
    findall(i -> isequal(getcell(_col(acs, attr), i), val), 1:acs.parts[ATTR2OBJ[attr]])

function add_part!(acs::ReactionNetworkSchema, obj::Symbol; kwargs...)
    n = (acs.parts[obj] += 1)
    for a in keys(SCHEMA[obj])
        grow!(_col(acs, a))
    end
    for (k, v) in kwargs
        setcell!(_col(acs, k), n, v)
    end
    return n
end

function add_parts!(acs::ReactionNetworkSchema, obj::Symbol, m::Int)
    n0 = acs.parts[obj]
    for _ = 1:m
        add_part!(acs, obj)
    end
    return (n0+1):(n0+m)
end

# rem_parts! is SWAP-AND-POP (verified against ACSets 0.2.29: it moves the LAST row into each freed
# slot and shrinks, iterating the sorted victims in REVERSE), NOT shift-down. equalize.jl:52 is the
# sole reindexer and the surviving-row ORDER it produces feeds species→state.u indexing / the sol
# DataFrame columns / valuation dot-products, so this must clone the swap-and-pop order exactly.
function rem_parts!(acs::ReactionNetworkSchema, obj::Symbol, idxs)
    idxs = issorted(idxs) ? idxs : sort(idxs)
    for p in Iterators.reverse(idxs)
        last = acs.parts[obj]
        for a in keys(SCHEMA[obj])
            c = _col(acs, a)
            if p != last
                # Move the last row into the freed slot. Guard on `def[last]`: for a non-bits column
                # (e.g. Vector{Set{Symbol}}, Vector{FoldedObservable}) an UNDEFINED last cell is a
                # `#undef` slot, and reading `c.v[last]` would throw UndefRefError — so only copy the
                # value when it is defined; otherwise just carry the (un)defined flag.
                if c.def[last]
                    c.v[p] = c.v[last]
                end
                c.def[p] = c.def[last]
            end
            resize!(c.v, last - 1)
            pop!(c.def)
        end
        acs.parts[obj] -= 1
    end
    return acs
end

Base.convert(::Type{Symbol}, ex::String) = Symbol(ex)

Base.convert(::Type{Union{String,Symbol,Missing}}, ex::String) =
    try
        Symbol(ex)
    catch
        string(ex)
    end

Base.convert(::Type{SampleableValues}, ex::String) = MacroTools.striplines(Meta.parse(ex))

# The Set{Symbol}/FoldedObservable string→eval convert hooks were removed (ADR 0005): they
# `eval`'d attribute strings on assignment (an import-time RCE vector). The JSON loader builds
# these as typed values directly (modality 3-axis → Set via to_set, observables structurally),
# so no string-eval path remains. (The SampleableValues parse above is parse-only — no eval —
# and is retained for the legacy DSL string-attr assignment.)

prettynames = Dict(
    :transRate => [:rate],
    :specInitUncertainty => [:uncertainty, :stoch, :stochasticity],
    :transPreAction => [:preAction, :action, :pre],
    :transPostAction => [:postAction, :post],
    :transName => [:name, :interpretation],
    :transPriority => [:priority],
    :transProbOfSuccess => [:probability, :prob, :pos],
    :transCapacity => [:cap, :capacity],
    :transCycleTime => [:ct, :cycletime],
    :transMaxLifeTime => [:lifetime, :maxlifetime, :maxtime, :timetolive],
)

defargs = Dict(
    :T => Dict{Symbol,Any}(
        :transPriority => 1,
        :transProbOfSuccess => 1,
        :transCapacity => Inf,
        :transCycleTime => 0.0,
        :transMaxLifeTime => Inf,
        :transMultiplier => 1,
        :transPreAction => :(),
        :transPostAction => :(),
        :transName => missing,
    ),
    :S => Dict{Symbol,Any}(
        :specInitUncertainty => 0.0,
        :specInitVal => 0.0,
        :specCost => 0.0,
        :specReward => 0.0,
        :specValuation => 0.0,
        :specStructured => false,
    ),
    :P => Dict{Symbol,Any}(:prmVal => missing),
    :M => Dict{Symbol,Any}(:metaVal => missing),
)

# (`compilable_attrs` removed with the ACSets swap — it was dead: `eltype(::Symbol)==SampleableValues`
# is never true, so the filter was always empty, and it had zero references anywhere.)

species_modalities = [:nonblock, :conserved, :rate]

function assign_defaults!(acs::ReactionNetworkSchema)
    for (_, v_) in defargs, (k, v) in v_
        for i in dom_parts(acs, k)
            isnothing(acs[i, k]) && (acs[i, k] = v)
        end
    end

    foreach(
        i -> !isnothing(acs[i, :specModality]) || (acs[i, :specModality] = Set{Symbol}()),
        parts(acs, :S),
    )
    k = [:specCost, :specReward, :specValuation]
    foreach(
        k -> foreach(i -> !isnothing(acs[i, k]) || (acs[i, k] = 0.0), parts(acs, :S)),
        k,
    )

    return acs
end

function ReactionNetworkSchema(transitions, reactants, obs, events)
    return merge_acs!(ReactionNetworkSchema(), transitions, reactants, obs, events)
end

function ReactionNetworkSchema(transitions, reactants, obs)
    return merge_acs!(ReactionNetworkSchema(), transitions, reactants, obs, [])
end

function add_obs!(acs, obs)
    for p in obs
        sym = p.args[3].value
        i = incident(acs, sym, :obsName)
        i = if isempty(incident(acs, sym, :obsName))
            add_part!(acs, :obs; obsName = sym, obsOpts = FoldedObservable())
        else
            i[1]
        end
        for opt in p.args[4:end]
            if isexpr(opt, :(=)) && (opt.args[1] ∈ fieldnames(FoldedObservable))
                opt.args[1] == :every &&
                    (acs[i, :obsOpts].every = min(acs[i, :obsOpts].every, opt.args[2]))
                opt.args[1] == :on && union!(acs[i, :obsOpts].on, [opt.args[2]])
            elseif isexpr(opt, :tuple) || opt isa SampleableValues
                push!(
                    acs[i, :obsOpts].range,
                    isexpr(opt, :tuple) ? tuple(opt.args...) : opt,
                )
            end
        end
    end

    return acs
end

function merge_acs!(acs::ReactionNetworkSchema, transitions, reactants, obs, events)
    foreach(
        t -> add_part!(acs, :T; trans = t[1][2], transRate = t[1][1], t[2]...),
        transitions,
    )
    add_obs!(acs, obs)
    unique!(reactants)
    foreach(
        ev -> add_part!(acs, :E; eventTrigger = ev.trigger, eventAction = ev.action),
        events,
    )
    foreach(
        r -> isempty(incident(acs, r, :specName)) && add_part!(acs, :S; specName = r),
        reactants,
    )

    return assign_defaults!(acs)
end

include("state.jl")
include("exprnode.jl")
include("compilers.jl")
include.(readdir(joinpath(@__DIR__, "interface"); join = true))
include.(readdir(joinpath(@__DIR__, "utils"); join = true))
include.(readdir(joinpath(@__DIR__, "operators"); join = true))
include("solvers.jl")
include("predicates.jl")
include("actions.jl")
# Per-program (per-structured-token) ledger (MVP finding D) — placed after actions.jl so it sees
# the ReactionNetworkProblem/Transition types and the token accessors; the hooks in solvers.jl
# call into it. Kept in its own file to minimize merge surface with the allocator rewrite.
include("ledger.jl")
include("serialize.jl")
#include("optim.jl")
include("loadsave.jl")

# Phase-0.6 analysis & visualization layer (ADR 0013/0014, CONTRACT §14/§15). All read-only over a
# finished run; placed last so they see the run-state, ledger, predicate/sortkey, and serializer.
#   analysis.jl  — §14.1 per-token trajectory log helpers + §14.2 ensemble runner / EnsembleProblem.
#   export.jl    — §14.3 results export bundle (JSON+CSV core; Arrow via RDArrowExt weakdep).
#   visualize.jl — §15.2 network exec map (Layer A network_graph, Layer B to_graphviz/draw_network,
#                  Layer C exec_map). Result-plot recipes + `_draw` live in ext/RDPlotsExt.jl.
include("analysis.jl")
include("export.jl")
include("visualize.jl")

end
