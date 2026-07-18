```@meta
CurrentModule = ReactiveDynamics
```

# Authoring

The authoring surface is the `@reaction_network` DSL and the macros that populate a network with transitions, species, attributes, and the initial-state / parameter / solver metadata a run needs. A network authored here is a pure typed data artifact — no host code is captured — that is later handed to [`ReactionNetworkProblem`](@ref) for simulation.

```@docs
@reaction_network
@push
@add_species
@mode
@aka
@name_transition
@append_transitions
@cost
@reward
@valuation
@prob_init
@prob_uncertainty
@prob_params
@prob_meta
@periodic
@jump
@register
```

## Attributes & shorthands

Common transition attributes. When set inside `@reaction_network` (or the update macros) they may be referred to by any of their shorthand names.

| attribute | shorthand | interpretation |
| :----- | :----- | :----- |
| `transPriority` | `priority` | priority of a transition (influences resource allocation) |
| `transProbOfSuccess` | `probability` `prob` `pos` | probability that a transition terminates successfully |
| `transCapacity` | `cap` `capacity` | maximum number of concurrent instances of the transition |
| `transCycleTime` | `ct` `cycletime` | duration of a transition's instance (adjusted by resource allocation) |
| `transMaxLifeTime` | `lifetime` `maxlifetime` `maxtime` `timetolive` | maximal duration of a transition's instance |
| `transPostAction` | `postAction` `post` | action to be executed once a transition's instance terminates |
| `transName` | `name` `interpretation` | name of a transition, either a string or unquoted text |

Common species attributes.

| attribute | shorthand | interpretation |
| :----- | :----- | :----- |
| `specInitUncertainty` | `uncertainty` `stoch` `stochasticity` | uncertainty about a variable's initial state (modelled as a Gaussian standard deviation) |
| `specInitVal` | | initial value of a variable |

## Rate semantics

The "rate" term of a transition governs how many instances spawn per step. By default a bare numeric rate is a stochastic (Poisson) intensity: at each step `n ~ Poisson(rate * dt)` instances are spawned. To specify the rate as a cycle time instead, use `@ct(cycle_time)` — e.g. `@ct(ex), A --> B, ...`, a shorthand for `1/ex, A --> B, ...`. For a deterministic "rate", use `@deterministic(ex)`, where `ex` evaluates to a deterministic count (floored to whole instances) spawned per integrator step. Note that a deterministic count does **not** scale with the step length.
