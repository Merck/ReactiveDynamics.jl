# --------------------------------------------------------------------------------
# Structured (agent-based) species

#=
import Pkg
Pkg.activate(".")
Pkg.dev("../..")
Pkg.add(["AlgebraicAgents"])
=#

using ReactiveDynamics
using AlgebraicAgents

# Define the "symbolic" reaction network.
network = @ReactionNetworkSchema

# Below, we combine
# - "classical species" (continuous or discrete; considered as pure quantities);
# - "structured agents" (possibly with custom evolutionary function; these can appear both on LHS and RHS).

@push network begin
    # With specified intensities, generate experimental resources.
    ρ1, ∅ --> R1
    ρ2, ∅ --> R2

    # Generate "Molecule 1" (where the integer corresponds to a "state" of, e.g., experimental triage).
    ρ3, ∅ --> M1(@t(), rand(4))

    # Based on properties of particular "structured agent" assigned to the transition,
    # we can update the attributes of the instance of a transition (such as probability of success).

    # Transition "Molecule 1" into "Molecule 2."
    # Update transition probability based on properties of "M1," 
    # which was assigned as a "resource" to the transition.
    ρ4,
    R1 + M1 --> M2(@t(), rand(4)),
    preAction => update_prob_transition(state, transition)
end

@prob_init network R1 = 10 R2 = 15

# As for structured agents, we will need to instantiate the instances
# and add them to the instance of a network. But first, we still need to define these types.
@prob_init network M1 = 2 M2 = 0

@prob_params network ρ1 = 2 ρ2 = 1 ρ3 = 3 ρ4 = 4

@prob_meta network tspan = 100 dt = 1.0

# We use `@structured` macro, which is a convenience wrapper around `@aagent`),
# defined in ReactiveDynamics.jl
@structured network struct M1
    descriptor::Any
    time_created::Any
end

using Random

# Type `M1` lives in the scope of ReactiveDynamics.
# Accordingly, we have to explicitly declare the scope.
using ReactiveDynamics: M1

ReactiveDynamics.M1(time, descriptor) = M1("M1" * randstring(4), nothing, descriptor, time)

# We define the function which updates the transition probability.
# This has to be accessible from within the name scope of ReactiveDynamics.
@register begin
    update_prob_transition = function (state, transition)
        if !isnothing(transition) && !isempty(transition.bound_structured_agents)
            bound_agent = first(transition.bound_structured_agents)

            transition[:transProbOfSuccess] = min(1.0, sum(bound_agent.descriptor))
        end
    end
end

# Alternatively, we can define a structured agent type using
# the usual `@aagent` macro. This must be evaluated inside the scope
# of ReactiveDynamics.
@register begin
    @aagent BaseStructuredSpecies AbstractStructuredSpecies struct M2
        descriptor::Any
        time_created::Any
    end

    using Random
    M2(time, descriptor) = M2("M2" * randstring(4), nothing, descriptor, time)
end

# Let the network know that the species is structured.
ReactiveDynamics.register_structured_species!(network, :M2)

# --------------------------------------------------------------------------------
# Instantiate the network.
network_instance = ReactionNetworkProblem(network)

for i = 1:2
    add_structured_species!(network_instance, ReactiveDynamics.M1(0.0, rand(4)))
end

# --------------------------------------------------------------------------------
# Simulate the network.
simulate(network_instance, 10)
