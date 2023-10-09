module ReactiveDynamics

using ACSets
using Reexport
using MacroTools
using ComponentArrays

@reexport using GeneratedExpressions

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

TheoryReactionNetwork = BasicSchema(
    [:S, :T, :E, :obs, :P, :M], # species, transitions, events, processes (observables), model params, solver args
    [], # no homs
    [
        :SymbolicAttributeT,
        :DescriptiveAttributeT,
        :SampleableAttributeT,
        :ModalityAttributeT,
        :PcsOptT,
        :PrmAttributeT,
    ], # AttrTypes
    [
        # species
        (:specName, :S, :SymbolicAttributeT),
        (:specModality, :S, :ModalityAttributeT),
        (:specInitVal, :S, :SampleableAttributeT),
        (:specInitUncertainty, :S, :SampleableAttributeT),
        (:specCost, :S, :SampleableAttributeT),
        (:specReward, :S, :SampleableAttributeT),
        (:specValuation, :S, :SampleableAttributeT),
        # transitions
        (:trans, :T, :SampleableAttributeT),
        (:transPriority, :T, :SampleableAttributeT),
        (:transRate, :T, :SampleableAttributeT),
        (:transCycleTime, :T, :SampleableAttributeT),
        (:transProbOfSuccess, :T, :SampleableAttributeT),
        (:transCapacity, :T, :SampleableAttributeT),
        (:transMaxLifeTime, :T, :SampleableAttributeT),
        (:transPostAction, :T, :SampleableAttributeT),
        (:transMultiplier, :T, :SampleableAttributeT),
        (:transName, :T, :DescriptiveAttributeT),
        # events
        (:eventTrigger, :E, :SampleableAttributeT),
        (:eventAction, :E, :SampleableAttributeT),
        # observables
        (:obsName, :obs, :SymbolicAttributeT),
        (:obsOpts, :obs, :PcsOptT),
        # params, args
        (:prmName, :P, :SymbolicAttributeT),
        (:prmVal, :P, :PrmAttributeT),
        (:metaKeyword, :M, :SymbolicAttributeT),
        (:metaVal, :M, :SampleableAttributeT),
    ],
)

@acset_type FoldedReactionNetworkType(TheoryReactionNetwork)

const ReactionNetworkSchema = FoldedReactionNetworkType{
    Symbol,
    Union{String,Symbol,Missing},
    SampleableValues,
    Set{Symbol},
    FoldedObservable,
    Any,
}

Base.convert(::Type{Symbol}, ex::String) = Symbol(ex)

Base.convert(::Type{Union{String,Symbol,Missing}}, ex::String) =
    try
        Symbol(ex)
    catch
        string(ex)
    end

Base.convert(::Type{SampleableValues}, ex::String) = MacroTools.striplines(Meta.parse(ex))
Base.convert(::Type{Set{Symbol}}, ex::String) = eval(Meta.parse(ex))
Base.convert(::Type{FoldedObservable}, ex::String) = eval(Meta.parse(ex))

prettynames = Dict(
    :transRate => [:rate],
    :specInitUncertainty => [:uncertainty, :stoch, :stochasticity],
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
        :transPostAction => :(),
        :transName => missing,
    ),
    :S => Dict{Symbol,Any}(
        :specInitUncertainty => 0.0,
        :specInitVal => 0.0,
        :specCost => 0.0,
        :specReward => 0.0,
        :specValuation => 0.0,
    ),
    :P => Dict{Symbol,Any}(:prmVal => missing),
    :M => Dict{Symbol,Any}(:metaVal => missing),
)

compilable_attrs =
    filter(attr -> eltype(attr) == SampleableValues, propertynames(ReactionNetworkSchema()))

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
include("compilers.jl")
include.(readdir(joinpath(@__DIR__, "interface"); join = true))
include.(readdir(joinpath(@__DIR__, "utils"); join = true))
include.(readdir(joinpath(@__DIR__, "operators"); join = true))
include("solvers.jl")
#include("optim.jl")
include("loadsave.jl")

end
