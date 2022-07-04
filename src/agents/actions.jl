# default parameters
const transition_defaults = Dict{Symbol, Any}(
    :priority => 1, :probability_of_success => 1,
    :capacity => Inf, :cycle_time => 0,
    :max_life_time => Inf
    )

# contained within `ReactionNetworkTransition`
"Ongoing transition; an instance of a `ReactionNetworkTransition`."
mutable struct Transition
    uuid::UUID # transition uuid
    params_sampled # sampled attrs

    t::Float64; q::Float64 # time spawn, quantity
    state::Float64 # [0, 1] ∋ transition progress
end

# parametric transition
@rnagent ReactionNetworkTransition ReactionNetworkAction begin
    params_expression; params_interpreted; params_sampled

    ongoing_transitions::Vector{Transition}
end

@doc "Stateful transition in a reaction network." ReactionNetworkTransition 

"""
    ReactionNetworkTransition(name, parameters)
Initialize a parametric transition in a reaction network.

By default, parameters will assume default values from `transition_defaults`.
"""
function ReactionNetworkTransition(name, parameters=deepcopy(transition_defaults), all_species=[])
    o = ReactionNetworkTransition(); o.name = name

    # set up properties
    o.params_expression = merge!(parameters, transition_defaults)
    o.params_interpreted = Dict{Symbol, Any}()
    o.params_sampled = Dict{Symbol, Any}(:toSpawn => 0.)
    interpret_params!(o, all_species) # turn expression into functions, escape refs

    o.ongoing_transitions = Transition[]

    o
end

@params_interface ReactionNetworkTransition # default parameter interface

# instantious action, triggered at each integration step
@rnagent ReactionNetworkEvent ReactionNetworkAction begin
    event_expression; event_interpreted

    params_expression; params_interpreted
end

@doc "Instantious action in a reaction network, triggered at each integration step." ReactionNetworkEvent

"""
    ReactionNetworkEvent(name, event_body, params, species=[])
Initialize a parametric, instantious event in a reaction network.

Expression will translate into a function `(event, network) -> \$event_body`.
Optionally, provide a list of names of variables in the system;
this will resolve corresponding references in `expression`.
"""
function ReactionNetworkEvent(name, event_body, params, all_species=[])
    o = ReactionNetworkEvent(); o.name = name
    
    o.event_expression = event_body
    o.event_interpreted = interpret_eval(event_body, all_species)

    o.params_expression = params; o.params_interpreted = Dict{Symbol, Any}()
    interpret_params!(o, all_species)

    o
end

@params_interface ReactionNetworkEvent # default parameter interface

# implement AlgebraicAgents.jl interface - delegated to the containing reaction network
AlgebraicAgents._step!(o::ReactionNetworkEvent, _) = nothing
AlgebraicAgents._projected_to(::ReactionNetworkEvent) = nothing
AlgebraicAgents.getobservable(o::ReactionNetworkEvent, _...) = nothing

"Execute event in a reaction network."
AlgebraicAgents._interact!(o::ReactionNetworkEvent) = o.event_interpreted(getnetwork(o), o)

## pretty printing
function AlgebraicAgents.print_custom(io::IO, mime::MIME"text/plain", o::ReactionNetworkEvent)
    indent = get(io, :indent, 0)
    print(io, "\n", " "^(indent+3), "custom properties:\n")
    print(io, " "^(indent+3), crayon"italics", "event", ": ", crayon"reset", "\n")
    show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.event_expression) 
end