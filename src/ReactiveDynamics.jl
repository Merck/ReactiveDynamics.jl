module ReactiveDynamics

using Reexport
using MacroTools
using ComponentArrays

# ADR 0003 Phase 1: the static authoring/IR store is a dependency-free typed-struct-of-columns
# (no ACSets). The store type is `ReactionNetwork` (ADR 0015) with a public `.columns` NamedTuple of
# typed columns, so the ~70 indexing sites and the 8 `propertynames(net.columns)` reflection loops
# work uniformly. The store verbs (`nrows`/`row_ids`/`column`/`cell`/`find_rows`/`add_row!`/… — see
# the SCHEMA + shim block below) are an INTERNAL store shim, NOT exported (ADR 0015 Tier 2 retired
# the old exported ACSets vocabulary `nparts`/`subpart`/`incident`/…; deprecated aliases live at the
# bottom of this file for one release).
# ADR 0003 Phase 2: the promoted transition↔reactant incidence table + its accessors.
export ReactantSpec, reactant_specs, specname

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
# ORDER (S-cols, then T, E, obs, P, M) is load-bearing: the store's `.columns` NamedTuple is built
# in this order, so `propertynames(net.columns)` reproduces the exact ACSets column order the eight
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
        # ADR 0009 §A / CONTRACT §11.1 — open-port role: a thin closed tag on the Species record
        # (NOT a new table). role ∈ {:private (default, auto-namespaced m__X on compose), :input,
        # :output (open ports; directionality advisory), :shared (bare-name identified, the
        # first-class @catchall)}. Drives @compose port-matching and refine boundary identification.
        specRole = Symbol,
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
# replace the `propertynames(net.columns)` reflection where a schema-driven list is wanted; the
# in-place `.columns` loops keep using `propertynames` since `.columns` IS a NamedTuple.
columns(::typeof(SCHEMA)) = ALLATTRS
columns(::typeof(SCHEMA), obj::Symbol) = keys(SCHEMA[obj])

# A single typed column: values plus a `defined` bitmap. An unassigned cell reads back as `nothing`
# (cloning ACSets' `get(col, i, default=nothing)`), so the `isnothing(net[i,k])` guards in
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

# ── ADR 0003 Phase 2: the promoted transition↔reactant incidence table ────────────────────────
#
# A `ReactantSpec` is one row of the bipartite transition↔species relation, promoted from the
# re-parsed `:trans` Expr into a first-class typed record with INTEGER foreign keys: `trans` → a :T
# index, `species` → a :S index. This makes the model's defining relation FK-checkable for agentic
# authoring and — the headline win — lets `equalize!` merge species by structurally REPOINTING the
# `species` FK (§7.4/J7) instead of `recursively_substitute_vars!` string surgery that can corrupt a
# name colliding inside a subexpression.
#
# The `expr` field is the ADR-mandated ESCAPE-HATCH for the legitimately dynamic reactants that are
# NOT a static (species, stoich) pair — `@choose` (random species per call), `@move`/`@structured`/
# `@advance` (RHS macrocalls), `@select(Kind, …)` (a token PREDICATE whose "species" is a structured
# kind, not a :S row), and expression-valued stoichiometry. Such a row carries `species = 0` (the
# "no static FK" sentinel — the store reads an undefined Int cell as `nothing`, but here the table is
# a plain Vector so we use 0 explicitly) and stashes the original term Expr in `expr`; a static row
# has `species ≥ 1` and `expr === nothing`.
#
# The table is DERIVED from the authoritative `:trans` column (see `populate_reactant_specs!`); the
# runtime engine still parses `:trans` per tick (state.jl), so promoting the table is additive and
# behavior-preserving. It lives as a struct field (NOT a 7th SCHEMA object) precisely so it never
# enters `propertynames(net.columns)` — the eight reflection loops in compilers/solvers/joins/
# equalize keep iterating exactly the six original objects' columns, untouched.
struct ReactantSpec
    trans::Int              # FK → :T
    species::Int            # FK → :S, or 0 for a dynamic (expr-carried) reactant
    stoich::SampleableValues
    side::Symbol            # :lhs or :rhs
    modality::Set{Symbol}
    expr::Union{Nothing,Expr,Symbol}   # escape-hatch term for a dynamic reactant, else nothing
end

# The static network container `ReactionNetwork` (ADR 0015: renamed from the ACSets-lineage
# `ReactionNetworkSchema` — it is a populated network INSTANCE, not the schema; the type-level
# object model is `const SCHEMA`). `counts` counts rows per object; `columns` is the NamedTuple of
# typed columns in ALLATTRS order; `reactants` is the promoted ReactantSpec incidence table (ADR
# 0003 Phase 2), populated lazily/on-merge (empty for a freshly-constructed or not-yet-promoted
# model — the runtime never reads it).
struct ReactionNetwork
    counts::Dict{Symbol,Int}
    columns::NamedTuple
    reactants::Vector{ReactantSpec}
    # Explicit TYPED inner constructor. Without it Julia auto-generates an untyped
    # `ReactionNetwork(::Any,::Any,::Any)` field constructor, which collides with the legacy
    # semantic 3-arg outer constructor `ReactionNetwork(transitions, reactants, obs)` below
    # (method overwrite → precompile error). The typed inner ctor is the only field-init path.
    ReactionNetwork(counts::Dict{Symbol,Int}, columns::NamedTuple,
                          reactants::Vector{ReactantSpec}) = new(counts, columns, reactants)
end

function ReactionNetwork()
    counts = Dict{Symbol,Int}(obj => 0 for obj in keys(SCHEMA))
    cols = NamedTuple{ALLATTRS}(AttrColumn{SCHEMA[ATTR2OBJ[a]][a]}() for a in ALLATTRS)
    return ReactionNetwork(counts, cols, ReactantSpec[])
end

# Deprecated compatibility alias (ADR 0015 Tier 1): the store type was renamed
# `ReactionNetworkSchema` → `ReactionNetwork`. Old caller code / annotations keep working (with a
# depwarn) for one release; removed in a follow-up.
Base.@deprecate_binding ReactionNetworkSchema ReactionNetwork

# ── ReactantSpec accessors (ADR 0003 Phase 2 public surface) ──────────────────────────────────
# The promoted incidence table. `populate_reactant_specs!` (serialize.jl) fills it from `:trans`;
# `equalize!` keeps it FK-exact across a species merge. A caller that wants the table on a model
# authored before promotion can call `populate_reactant_specs!(net)` first (equalize! does).
reactant_specs(net::ReactionNetwork) = net.reactants

# The :S species name at index `i` (the inverse of `find_index`), used to check FK targets.
specname(net::ReactionNetwork, i::Integer) = net[i, :specName]

# The :S index of a species name on the STATIC schema (the ReactionNetworkProblem overload lives in
# state.jl:263). Returns nothing if absent. Used by equalize!'s FK-repoint and the acceptance tests.
function find_index(species::Symbol, net::ReactionNetwork)
    inc = find_rows(net, species, :specName)
    return isempty(inc) ? nothing : first(inc)
end

# ── ADR 0009 §A / CONTRACT §11.1 — open-port roles (a closed tag on the Species record) ────────
const PORT_ROLES = (:private, :input, :output, :shared)

# The role of species `i`, defaulting to :private (a species authored before roles existed, or one
# whose specRole cell is unset, is internal/namespaced). Assign-defaults seeds :private, but read
# defensively so `port_role` is correct on a not-yet-defaulted schema too.
function port_role(net::ReactionNetwork, i::Integer)
    r = net[i, :specRole]
    return (r === nothing || r === missing) ? :private : r
end
port_role(net::ReactionNetwork, name::Symbol) =
    (i = find_index(name, net); i === nothing ? nothing : port_role(net, i))

is_open_port(role::Symbol) = role === :input || role === :output

# ACSets hashed a static model by CONTENT (so export.jl `_model_hash = hash(prob.network)` names a
# stable bundle dir); a struct's default hash is object-identity. Preserve content-hashing. The
# promoted reactant table is DERIVED from `:trans`, so it is intentionally NOT hashed — a model and
# its post-`populate_reactant_specs!` self must hash identically (the table adds no new information),
# keeping `_model_hash` stable across the Phase-2 promotion.
function Base.hash(net::ReactionNetwork, h::UInt)
    h = hash(:ReactionNetwork, h)
    for a in ALLATTRS
        col = net.columns[a]
        for i = 1:net.counts[ATTR2OBJ[a]]
            h = hash(getcell(col, i), h)
        end
    end
    return h
end

# ── the signature-preserving shim (replaces the ACSets API surface RD used) ───────────────────
@inline _col(net::ReactionNetwork, attr::Symbol) = getfield(net, :columns)[attr]

nrows(net::ReactionNetwork, obj::Symbol) = net.counts[obj]
row_ids(net::ReactionNetwork, obj::Symbol) = Base.OneTo(net.counts[obj])
col_row_ids(net::ReactionNetwork, attr::Symbol) = Base.OneTo(net.counts[ATTR2OBJ[attr]])

# scalar get/set
Base.getindex(net::ReactionNetwork, i::Int, attr::Symbol) = getcell(_col(net, attr), i)
Base.setindex!(net::ReactionNetwork, v, i::Int, attr::Symbol) = setcell!(_col(net, attr), i, v)
# whole-column and row-subset reads return COPIES (matching ACSets `collect_column`/`map(identity)`);
# `map(identity, …)` narrows e.g. `net[:, :specName]` back to `Vector{Symbol}`.
Base.getindex(net::ReactionNetwork, ::Colon, attr::Symbol) =
    map(identity, [getcell(_col(net, attr), i) for i = 1:net.counts[ATTR2OBJ[attr]]])
Base.getindex(net::ReactionNetwork, rows::AbstractVector, attr::Symbol) =
    [getcell(_col(net, attr), i) for i in rows]

column(net::ReactionNetwork, attr::Symbol) = net[:, attr]           # COPY (see note)
cell(net::ReactionNetwork, i::Int, attr::Symbol) = net[i, attr]
set_cell!(net::ReactionNetwork, i::Int, attr::Symbol, v) = (net[i, attr] = v)

# `find_rows(net, val, attr)` = findall over the attr column (never an FK follow — no homs). `isequal`
# (not `==`) so `nothing`/`missing` cells compare `false`, never poison the result with `missing`.
find_rows(net::ReactionNetwork, val, attr::Symbol) =
    findall(i -> isequal(getcell(_col(net, attr), i), val), 1:net.counts[ATTR2OBJ[attr]])

function add_row!(net::ReactionNetwork, obj::Symbol; kwargs...)
    n = (net.counts[obj] += 1)
    for a in keys(SCHEMA[obj])
        grow!(_col(net, a))
    end
    for (k, v) in kwargs
        setcell!(_col(net, k), n, v)
    end
    return n
end

function add_rows!(net::ReactionNetwork, obj::Symbol, m::Int)
    n0 = net.counts[obj]
    for _ = 1:m
        add_row!(net, obj)
    end
    return (n0+1):(n0+m)
end

# rem_rows! is SWAP-AND-POP (verified against ACSets 0.2.29: it moves the LAST row into each freed
# slot and shrinks, iterating the sorted victims in REVERSE), NOT shift-down. equalize.jl:52 is the
# sole reindexer and the surviving-row ORDER it produces feeds species→state.u indexing / the sol
# DataFrame columns / valuation dot-products, so this must clone the swap-and-pop order exactly.
function rem_rows!(net::ReactionNetwork, obj::Symbol, idxs)
    idxs = issorted(idxs) ? idxs : sort(idxs)
    for p in Iterators.reverse(idxs)
        last = net.counts[obj]
        for a in keys(SCHEMA[obj])
            c = _col(net, a)
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
        net.counts[obj] -= 1
    end
    return net
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
        :specRole => :private,          # ADR 0009 §A — default port role (internal, namespaced)
    ),
    :P => Dict{Symbol,Any}(:prmVal => missing),
    :M => Dict{Symbol,Any}(:metaVal => missing),
)

# (`compilable_attrs` removed with the ACSets swap — it was dead: `eltype(::Symbol)==SampleableValues`
# is never true, so the filter was always empty, and it had zero references anywhere.)

species_modalities = [:nonblock, :conserved, :rate]

function assign_defaults!(net::ReactionNetwork)
    for (_, v_) in defargs, (k, v) in v_
        for i in col_row_ids(net, k)
            isnothing(net[i, k]) && (net[i, k] = v)
        end
    end

    foreach(
        i -> !isnothing(net[i, :specModality]) || (net[i, :specModality] = Set{Symbol}()),
        row_ids(net, :S),
    )
    k = [:specCost, :specReward, :specValuation]
    foreach(
        k -> foreach(i -> !isnothing(net[i, k]) || (net[i, k] = 0.0), row_ids(net, :S)),
        k,
    )

    return net
end

function ReactionNetwork(transitions, reactants, obs, events)
    return merge_network!(ReactionNetwork(), transitions, reactants, obs, events)
end

function ReactionNetwork(transitions, reactants, obs)
    return merge_network!(ReactionNetwork(), transitions, reactants, obs, [])
end

function add_obs!(net, obs)
    for p in obs
        sym = p.args[3].value
        i = find_rows(net, sym, :obsName)
        i = if isempty(find_rows(net, sym, :obsName))
            add_row!(net, :obs; obsName = sym, obsOpts = FoldedObservable())
        else
            i[1]
        end
        for opt in p.args[4:end]
            if isexpr(opt, :(=)) && (opt.args[1] ∈ fieldnames(FoldedObservable))
                opt.args[1] == :every &&
                    (net[i, :obsOpts].every = min(net[i, :obsOpts].every, opt.args[2]))
                opt.args[1] == :on && union!(net[i, :obsOpts].on, [opt.args[2]])
            elseif isexpr(opt, :tuple) || opt isa SampleableValues
                push!(
                    net[i, :obsOpts].range,
                    isexpr(opt, :tuple) ? tuple(opt.args...) : opt,
                )
            end
        end
    end

    return net
end

function merge_network!(net::ReactionNetwork, transitions, reactants, obs, events)
    foreach(
        t -> add_row!(net, :T; trans = t[1][2], transRate = t[1][1], t[2]...),
        transitions,
    )
    add_obs!(net, obs)
    unique!(reactants)
    foreach(
        ev -> add_row!(net, :E; eventTrigger = ev.trigger, eventAction = ev.action),
        events,
    )
    foreach(
        r -> isempty(find_rows(net, r, :specName)) && add_row!(net, :S; specName = r),
        reactants,
    )

    return assign_defaults!(net)
end

# ── Deprecated ACSets-vocabulary aliases (ADR 0015 Tier 2) ────────────────────────────────────
# The store verbs were renamed to store vocabulary AND unexported (they are an internal store shim).
# These thin, UN-exported aliases keep prior `RD.nparts(...)`-style calls working for one release
# with a depwarn; removed in a follow-up. They forward to the new names, which dispatch on
# ReactionNetwork / ReactionNetworkProblem (the RNP overloads live in state.jl). `add_part!` is
# written by hand because `@deprecate` cannot express its `kwargs`.
@deprecate nparts(x, obj) nrows(x, obj) false
@deprecate parts(x, obj) row_ids(x, obj) false
@deprecate dom_parts(x, attr) col_row_ids(x, attr) false
@deprecate incident(x, val, attr) find_rows(x, val, attr) false
@deprecate subpart(x, attr) column(x, attr) false
@deprecate subpart(x, i::Int, attr) cell(x, i, attr) false
@deprecate set_subpart!(x, i, attr, v) set_cell!(x, i, attr, v) false
@deprecate add_parts!(x, obj, m) add_rows!(x, obj, m) false
@deprecate rem_parts!(x, obj, idxs) rem_rows!(x, obj, idxs) false
function add_part!(x, obj::Symbol; kwargs...)
    Base.depwarn("`add_part!` is deprecated (ADR 0015); use `ReactiveDynamics.add_row!`.", :add_part!)
    return add_row!(x, obj; kwargs...)
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
