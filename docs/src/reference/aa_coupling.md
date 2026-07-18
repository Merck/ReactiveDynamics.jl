```@meta
CurrentModule = ReactiveDynamics
```

# AlgebraicAgents coupling

A `ReactionNetworkProblem` *is* an [AlgebraicAgents.jl](https://github.com/Merck/AlgebraicAgents.jl) (AA) agent, so a network is a node in a larger heterogeneous AA hierarchy: it can be entangled with sibling agents, expose observables, and read a sibling's output off a wire (with a one-tick Jacobi lag). This page also covers checkpointing a live run.

## AA verbs (reexported)

`entangle!`, `add_wire!`, and `getobservable` are AlgebraicAgents functions — `getobservable` is overloaded by ReactiveDynamics for `ReactionNetworkProblem`, and all three are reexported (via `@reexport using AlgebraicAgents`); they are documented in the AlgebraicAgents documentation. `entangle!(parent, rd)` makes the network a child node of `parent`; `add_wire!(root; from, to, from_var_name, to_var_name)` connects one agent's observable to another agent's input port; `getobservable(rd, name)` reads a network's exported observable by name (or `getobservable(rd, i::Int)` by canonical index). A cross-agent read is expressed on the RD side with an [`ExternalRef`](@ref) (see [Serialization](@ref) for its IR docstring): the referenced value resolves one tick late, giving a deterministic Jacobi coupling.

```julia
using ReactiveDynamics
using AlgebraicAgents            # reexported by ReactiveDynamics

# rd is a ReactionNetworkProblem; finance is a sibling AA agent under root
entangle!(root, rd)                                                         # AA verb — make rd a node in the hierarchy
add_wire!(root; from = rd, to = finance,                                    # AA verb — wire rd's :cash observable
    from_var_name = "cash", to_var_name = "rd_cash")                        #   into finance's :rd_cash input port
getobservable(rd, :cash)                                                    # AA verb — read the observable by name
```

## Checkpointing

`dump_state`/`restore` snapshot and reload a live run through a `StateDump`; these are ReactiveDynamics exports. (`dump_state` requires a clean tick boundary.)

```@docs
StateDump
dump_state
restore
```
