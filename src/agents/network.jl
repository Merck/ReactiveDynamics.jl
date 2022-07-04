# supertype of transition network agents, this subtypes abstract algebraic agent
abstract type ReactionNetworkAgent <: AbstractAlgebraicAgent end
# supertype of transitions: stateful, instantious actions
abstract type ReactionNetworkAction <: ReactionNetworkAgent end

# a transition complex consisting of reactants, transitions, and events
@rnagent ReactionNetwork begin
    log::Vector{Tuple}
end

@doc "A transition complex consisting of reactants, transitions, and events."

"""
    ReactionNetwork(name, expression=:())
Initialize an (empty, by default) transition network.
If `expression` is provided, parse the transition dynamics and populate the transition network.
"""
function ReactionNetwork(name::AbstractString, expression::Expr; _...)
    o = ReactionNetwork(); o.name = name
    o.log = Tuple[]

    transitions, species = parse_transitions(expression)
    network_append!(o; species, transitions)

    o
end

"Append objects to a reation network."
function network_append!(rn::ReactionNetwork; species=[], transitions=[], events=[], sampleables=[])
    # add species
    for r in unique!(species)
        parent = mkpath!(rn, dirname(r)) # create intermediate hierachy

        if !haskey(inners(parent), basename(r))
            entangle!(parent, ReactionNetworkSpecies(basename(r)))
        end
    end

    # add transitions
    all_species = get_all_species(rn)

    for t in transitions
        name = t[2][:name] = string(get(t[2], :name, gensym(:transition)))

        params = push!(t[2], :trans => t[1][2], :rate => t[1][1])
        if !haskey(inners(rn), name)
            entangle!(rn, ReactionNetworkTransition(name, params, all_species))
        else
            t = inners(rn)[name]
            merge!(t.expression, params); interpret_params!(t, all_species)
        end
    end

    # add events
    for e in events
        name = e.params[:name] = string(get(e.params, :name, gensym(:event)))

        entangle!(rn, ReactionNetworkEvent(name, e.event_body, e.params, all_species))
    end

    # add sampleables
    for o in sampleables
        name = o.params[:name] = string(get(o.params, :name, gensym(:event)))
        parent = mkpath!(rn, dirname(name)) # create intermediate hierachy

        entangle!(parent, ReactionNetworkSampleable(basename(name), o.range, o.params, all_species))
    end

    rn
end

"Return parent reaction network."
getnetwork(a::ReactionNetworkAgent) = a isa ReactionNetwork ? a : getnetwork(parent(a))

"Translate `a` into `(a, getnetwork(a))`."
macro args(a) :(($(esc(a)), getnetwork($(esc(a))))) end

"Retrieve species names as symbols (ignores nested hierarchies)."
function get_all_species(rn::ReactionNetwork)
    [Symbol(getname(v)) for (_, v) in inners(rn) if v isa ReactionNetworkSpecies]
end