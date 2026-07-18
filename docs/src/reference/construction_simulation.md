```@meta
CurrentModule = ReactiveDynamics
```

# Construction & simulation

An authored network becomes a runnable problem via [`ReactionNetworkProblem`](@ref), which owns the live simulation state and the per-run RNG (a run is determined by `(model, seed)`). `@agentize` builds the problem inside an AlgebraicAgents hierarchy.

```@docs
ReactionNetworkProblem
@agentize
```

## Stepping the problem

A `ReactionNetworkProblem` *is* an [AlgebraicAgents.jl](https://github.com/Merck/AlgebraicAgents.jl) agent, so it is stepped and reset with AA's verbs — `simulate` (advance the run to its horizon) and `reinit!` (restore the exact RNG stream for a fresh run, taking a `seed` kwarg). These are AlgebraicAgents functions, reexported by ReactiveDynamics; see the AlgebraicAgents documentation for their full signatures.

```julia
prob = ReactionNetworkProblem(net; seed = 1)
simulate(prob)     # AlgebraicAgents verb, reexported by RD
reinit!(prob)      # AlgebraicAgents verb, reexported by RD
```
