# default parameters
const species_defaults = Dict{Symbol, Any}(
    :speciesInitUncertainty => .0, :speciesInitVal => .0,
    :speciesCost => .0, :speciesReward => .0,
    :speciesValuation => .0, :speciesModality => Set{Symbol}()
    )

@rnagent ReactionNetworkSpecies ReactionNetworkAgent begin
    u; u0 # current quantity of species, initial quantity
    params_expression; params_interpreted; params_sampled # parameters
end

@doc "(Parametric) species in a reaction network which is associated with quantity."

"""
    ReactionNetworkSpecies(name, u0, p)
Initialize a parametric species in a reaction network.
"""
function ReactionNetworkSpecies(name::AbstractString, u0=missing,
    params=deepcopy(species_defaults), all_species=String[])
    
    o = ReactionNetworkSpecies(); o.name = name
    
    o.u0 = u0; o.u = deepcopy(o.u0)
    o.params_expression = merge(params, species_defaults)
    o.params_interpreted = Dict{Symbol, Any}()
    o.params_sampled = Dict{Symbol, Any}()

    interpret_params!(o, all_species)

    o
end

# implement AlgebraicAgents.jl interface
## pretty printing
function AlgebraicAgents.print_custom(io::IO, mime::MIME"text/plain", o::ReactionNetworkSpecies)
    indent = get(io, :indent, 0)
    print(io, "\n", " "^(indent+3), "custom properties:\n")
    print(io, " "^(indent+3), crayon"italics", "u", ": ", crayon"reset", "\n")
    show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.u)

    print(io, "\n" * " "^(indent+3), crayon"italics", "parameters", ": ", crayon"reset", "\n")
    show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.params_expression)
end

## (de)reference
function Base.getindex(obj::ReactionNetworkSpecies, keys...)
    if length(obj.u) > 1
        isempty(keys) ? obj.u : obj.u[keys...]
    else first(obj.u) end
end

function Base.setindex!(obj::ReactionNetworkSpecies, val, keys...)
    if isempty(keys) obj.u .= val
    else obj.u[keys...] .= val end
end

AlgebraicAgents._step!(::ReactionNetworkSpecies, _) = nothing
AlgebraicAgents._projected_to(::ReactionNetworkSpecies) = nothing
AlgebraicAgents._reinit!(o::ReactionNetworkSpecies) = o.u = deepcopy(o.u0)
AlgebraicAgents.getobservable(o::ReactionNetworkSpecies, obs...) = o[obs...]

@params_interface ReactionNetworkSpecies # default parameter interface