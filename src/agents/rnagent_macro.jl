# same as AlgebraicAgents.jl's `@aagent` macro, but doesn't set `:name` field

"""
    @rnagent agent_name begin
        extra_fields...
    end

Create a custom algebraic agent type, and include fields expected by default interface methods (see [`FreeAgent`](@ref)).

# Example 
```julia
@rnagent ReactionNetworkSpecies begin
    u; u0 # current quantity of species, initial quantity
    params_expression; params_interpreted; params_sampled # parameters
end
```
"""
macro rnagent(new_name, extra_fields)
    AlgebraicAgents.define_agent(new_name, FreeAgent, extra_fields, __module__, quote
            function $(new_name)()
                m = new()
                m.uuid = AlgebraicAgents.uuid4()
                m.parent = nothing; m.inners = Dict{String, AbstractAlgebraicAgent}()
                m.relpathrefs = Dict{AbstractString, AlgebraicAgents.UUID}()
                m.opera = AlgebraicAgents.Opera(m.uuid => m)

                m
            end
        end
    )
end

"""
    @rnagent agent_name supertype begin
        extra_fields...
    end

Create a custom algebraic agent type, and include fields expected by default interface methods (see [`FreeAgent`](@ref)).

# Example 
```julia
@rnagent ReactionNetworkSpecies ReactionNetworkAgent begin
    u; u0 # current quantity of species, initial quantity
    params_expression; params_interpreted; params_sampled # parameters
end
```
"""
macro rnagent(new_name, super_type, extra_fields)
    AlgebraicAgents.define_agent_with_supertype(new_name, FreeAgent, super_type, extra_fields, __module__, quote
            function $(new_name)()
                m = new()
                m.uuid = AlgebraicAgents.uuid4()
                m.parent = nothing; m.inners = Dict{String, AbstractAlgebraicAgent}()
                m.relpathrefs = Dict{AbstractString, AlgebraicAgents.UUID}()
                m.opera = AlgebraicAgents.Opera(m.uuid => m)

                m
            end
        end
    )
end