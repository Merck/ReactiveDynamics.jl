# default parameters
const sampleable_defaults = Dict{Symbol, Any}(:every => 1., :on => false)

@rnagent ReactionNetworkSampleable ReactionNetworkAgent begin
    range_expression; range_interpreted

    params_expression; params_interpreted

    last_sampled::Float64 # last sampling time
    sampled # last sampled value
end

@doc "Sampleable (dynamically resampled process)." ReactionNetworkSampleable

"""
    ReactionNetworkSampleable(name, range, params, species=[])
Initialize a parametric sampleable in a reaction network.

Range will translate into `rng = (event, network) -> \$range`;
the sampling range will equal `rng(event, network)`.

Optionally, resample the sampleable's value 
whenever parameter `on` evaluates to true, or after `every` time steps.

Optionally, provide a list of names of variables in the system;
this will resolve corresponding references in `expression`.
"""
function ReactionNetworkSampleable(name, range, params, all_species=[])
    o = ReactionNetworkSampleable(); o.name = name
   
    o.range_expression = range
    o.range_interpreted = interpret_eval(range, all_species)

    o.params_expression = merge!(params, sampleable_defaults)
    o.params_interpreted = Dict{Symbol, Any}()
    interpret_params!(o, all_species)

    o.last_sampled = -Inf; o.sampled = missing

    o
end

@params_interface ReactionNetworkSampleable # default parameter interface

"Resample a sampleable agent."
function resample!(o::ReactionNetworkSampleable, t)
    o.last_sampled = t + sample_params(o)[:every]

    o.sampled = take_sample(o)
end

"Take sample of sampleable agent's stochastic process."
function take_sample(o::ReactionNetworkSampleable, rn=getnetwork(o))
    r = o.range_interpreted(rn, o) 
    r isa Sampleable ? rand(r) : r
end

# implement AlgebraicAgents.jl interface
function AlgebraicAgents._step!(o::ReactionNetworkSampleable, t)
    (projected_to(o) <= t) && resample!(o, t)

    o.last_sampled
end

AlgebraicAgents._projected_to(::ReactionNetworkSampleable) = o.last_sampled

## (de)reference
AlgebraicAgents.getobservable(o::ReactionNetworkSampleable) = o.sampled

## pretty printing
function AlgebraicAgents.print_custom(io::IO, mime::MIME"text/plain", o::ReactionNetworkSampleable)
    indent = get(io, :indent, 0)
    print(io, "\n", " "^(indent+3), "custom properties:\n")
    print(io, " "^(indent+3), crayon"italics", "sampleable", ": ", crayon"reset", "\n")
    show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.range_expression)

    print(io, "\n" * " "^(indent+3), crayon"italics", "resampled every", ": ", crayon"reset", "\n")
    show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.params_expression[:every])

    if (o.params_expression[:on] !== false)
        print(io, "\n" * " "^(indent+3), crayon"italics", "resampled on", ": ", crayon"reset", "\n")
        show(IOContext(io, :indent=>get(io, :indent, 0)+4), mime, o.params_expression[:on])
    end    
end