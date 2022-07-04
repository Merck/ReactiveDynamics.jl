"""
Run `GeneratedExpressions.generate` on expressions at input, and return a flat representation of the result.
Optionally specify top-level options for `generate` and a module to evaluate in.
"""
function expand(expr; tl_opts=[], eval_module=@__MODULE__)
    flatten!(generate(Expr(:braces, expr, tl_opts...); eval_module))
end

"""
Macro that takes an expression corresponding to a reaction network and outputs an instance of reaction network agent.

Most arrows are accepted (both right, left, and bi-drectional arrows). Use 0 or ∅ for annihilation/creation to/from nothing.

Custom functions and sampleable objects can be used as numeric parameters. Note these have to be accessible from DyVE's module.

# Examples
```julia
acs = @ReactionNetwork begin
    1.0, X ⟶ Y
    1.0, X ⟶ Y, priority=>6., prob=>.7, capacity=>3.
    1.0, ∅ --> (Poisson(.3γ)X, Poisson(.5)Y)
end
```
"""
macro ReactionNetwork(expr, tl_opts...)
    name, generated = if typeof(expr) ∈ [Symbol, String]
        string(expr), :() 
    else
        "reaction_network", expand(expr; tl_opts=collect(tl_opts), eval_module=__module__)
    end

    :(ReactionNetwork($name, $(QuoteNode(generated)); eval_module=$__module__))
end

"Append transitions to a reaction network."
macro transitions(network, expr, tl_opts...)
    generated = expand(expr; tl_opts=collect(tl_opts), eval_module=__module__)

    quote
        transitions, species = parse_transitions($(QuoteNode(generated)))
        network_append!($(esc(network)); transitions, species)
    end
end

"""
Add species to a reaction network.

# Examples
```julia
@species network a"A/B" "C" D
```
"""
macro species(network, exs...)
    species = []
    for expr in exs
        generated = expand(expr; eval_module=__module__)
        for s in (generated isa Expr && isexpr(generated, :block) ? generated.args : [generated])
            push!(species, to_string(s))
        end
    end
    
    quote
        network_append!($(esc(network)); species=$species)
    end
end

macro events(network, expr, tl_opts...)
    generated = expand(expr; tl_opts=collect(tl_opts), eval_module=__module__)

    quote
        events = parse_events($(QuoteNode(generated)))
        network_append!($(esc(network)); events)
    end
end

macro sampleables(network, expr, tl_opts...)
    generated = expand(expr; tl_opts=collect(tl_opts), eval_module=__module__)

    quote
        sampleables = parse_observables($(QuoteNode(generated)))
        network_append!($(esc(network)); sampleables)
    end
end