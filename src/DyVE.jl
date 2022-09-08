module DyVE

using Catlab, Catlab.CategoricalAlgebra, Catlab.Present
using Reexport
using MacroTools
using NLopt
using ComponentArrays

@reexport using GeneratedExpressions

const SampleableValues = Union{Expr, Symbol, AbstractString, Float64, Int, Function}
const ActionableValues = Union{Function, Symbol, Float64, Int}

const SampleableRange = Union{Float64, Int64, AbstractString, Expr, Symbol, Tuple{Float64, Union{Float64, Int64, AbstractString, Expr, Symbol}}}

Base.convert(::Type{SampleableRange}, x::Tuple) = (Float64(x[1]), x[2])

Base.@kwdef mutable struct FoldedObservable
    range::Vector{SampleableRange} = SampleableRange[]
    every::Float64 = Inf
    on::Vector{SampleableValues} = SampleableValues[]
end

@present TheoryReactionNetwork(FreeSchema) begin
    (S, T)::Ob # species, transitions

    (SymbolicAttributeT, DescriptiveAttributeT, SampleableAttributeT, 
        ModalityAttributeT, PcsOptT, PrmAttributeT)::AttrType

    specName::Attr(S, SymbolicAttributeT)
    specModality::Attr(S, ModalityAttributeT)
    specInitVal::Attr(S, SampleableAttributeT)
    specInitUncertainty::Attr(S, SampleableAttributeT)
    (specCost, specReward, specValuation)::Attr(S, SampleableAttributeT)

    trans::Attr(T, SampleableAttributeT)
    transPriority::Attr(T, SampleableAttributeT)
    transRate::Attr(T, SampleableAttributeT)
    transCycleTime::Attr(T, SampleableAttributeT)
    transProbOfSuccess::Attr(T, SampleableAttributeT)
    transCapacity::Attr(T, SampleableAttributeT)
    transMaxLifeTime::Attr(T, SampleableAttributeT)
    transPostAction::Attr(T, SampleableAttributeT)
    transMultiplier::Attr(T, SampleableAttributeT)
    transName::Attr(T, DescriptiveAttributeT)

    E::Ob # events
    (eventTrigger, eventAction)::Attr(E, SampleableAttributeT)

    obs::Ob # processes (observables)
    obsName::Attr(obs, SymbolicAttributeT)
    obsOpts::Attr(obs, PcsOptT)

    (P, M)::Ob # model params, solver args

    prmName::Attr(P, SymbolicAttributeT)
    prmVal::Attr(P, PrmAttributeT)

    metaKeyword::Attr(M, SymbolicAttributeT)
    metaVal::Attr(M, SampleableAttributeT)
end

@acset_type FoldedReactionNetworkType(TheoryReactionNetwork)

const ReactionNetwork = FoldedReactionNetworkType{Symbol, Union{String, Symbol, Missing}, SampleableValues, Set{Symbol}, FoldedObservable, Any}

Base.convert(::Type{Symbol}, ex::String) = Symbol(ex)
Base.convert(::Type{Union{String, Symbol, Missing}}, ex::String) = try Symbol(ex) catch; string(ex) end
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
    :transMaxLifeTime => [:lifetime, :maxlifetime, :maxtime, :timetolive]
)

defargs = Dict(
    :T => Dict{Symbol, Any}(:transPriority => 1, :transProbOfSuccess => 1, :transCapacity => Inf, :transCycleTime => 1,
        :transMaxLifeTime => Inf, :transMultiplier=>1, :transPostAction => :(), :transName => missing),
    :S => Dict{Symbol, Any}(:specInitUncertainty => .0, :specInitVal => .0, :specCost => .0, :specReward => .0, :specValuation => .0),
    :P => Dict{Symbol, Any}(:prmVal => missing),
    :M => Dict{Symbol, Any}(:metaVal => missing)
)

compilable_attrs = filter(attr -> eltype(attr) == SampleableValues, propertynames(ReactionNetwork()))

species_modalities = [:nonblock, :conserved, :rate]

function assign_defaults!(acs::ReactionNetwork)
    for (_, v_) in defargs, (k, v) in v_
        for i in eachindex(acs.attrs[k])
            !isassigned(acs.attrs[k], i) && (acs.attrs[k][i] = v)
        end
    end
    
    foreach(i -> isassigned(acs.attrs.specModality, i) || (acs.attrs.specModality[i] = Set{Symbol}()), 1:nparts(acs, :S))
    k = [:specCost, :specReward, :specValuation]
    foreach(k -> foreach(i -> isassigned(getproperty(acs.attrs, k), i) || (getproperty(acs.attrs, k)[i] = .0), 1:nparts(acs, :S)), k)

    acs
end

ReactionNetwork(transitions, reactants, obs, events) = merge_acs!(ReactionNetwork(), transitions, reactants, obs, events)
ReactionNetwork(transitions, reactants, obs) = merge_acs!(ReactionNetwork(), transitions, reactants, obs, [])

function add_obs!(acs, obs)
    for p in obs 
        sym = p.args[3].value; i = incident(acs, sym, :obsName)
        i = isempty(incident(acs, sym, :obsName)) ? add_part!(acs, :obs; obsName=sym, obsOpts=FoldedObservable()) : i[1]
        for opt in p.args[4:end]
            if isexpr(opt, :(=)) && (opt.args[1] ∈ fieldnames(FoldedObservable))
                opt.args[1] == :every && (acs[i, :obsOpts].every = min(acs[i, :obsOpts].every, opt.args[2]))
                opt.args[1] == :on && union!(acs[i, :obsOpts].on, [opt.args[2]])
            elseif isexpr(opt, :tuple) || opt isa SampleableValues
                push!(acs[i, :obsOpts].range, isexpr(opt, :tuple) ? tuple(opt.args...) : opt)
            end
        end
    end

    acs
end

function merge_acs!(acs::ReactionNetwork, transitions, reactants, obs, events)
    foreach(t -> add_part!(acs, :T; trans=t[1][2], transRate=t[1][1], t[2]...), transitions)
    add_obs!(acs, obs); unique!(reactants)
    foreach(ev -> add_part!(acs, :E; eventTrigger=ev.trigger, eventAction=ev.action), events)
    foreach(r -> isempty(incident(acs, r, :specName)) && add_part!(acs, :S; specName=r), reactants)

    assign_defaults!(acs)
end



include("state.jl"); include("compilers.jl")
include.(readdir(joinpath(@__DIR__, "interface"), join=true))
include.(readdir(joinpath(@__DIR__, "utils"), join=true))
include.(readdir(joinpath(@__DIR__, "operators"), join=true))
include("solvers.jl"); include("optim.jl")
include("loadsave.jl")

end