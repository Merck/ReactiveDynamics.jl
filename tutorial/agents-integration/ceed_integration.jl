using CEEDesigns, CEEDesigns.GenerativeDesigns
using DataFrames
using ScientificTypes
using Statistics, Copulas
import POMDPs, POMDPTools, MCTS

import Distributions

# ----- Experimental Setup -----

# We generate a synthetic dataset.
# This is taken from https://github.com/Merck/CEEDesigns.jl/blob/34588ae0e5563cb93f6818e3a9c8b3a77c5e3c47/tutorials/SimpleGenerative.jl

include("experimental_setup.jl")

# ----- Get a sampling function -----

(; sampler, uncertainty, weights) = DistanceBased(
    data;
    target = "y",
    uncertainty = Variance(),
    similarity = GenerativeDesigns.Exponential(; λ = 5),
);

# ----- Set up a reaction network -----

#=
Pkg.activate(".")
Pkg.develop(path = "../..")
=#

using ReactiveDynamics

# Set up parameters that will be used to define a network.

# Experiments and costs
features_experiments = Dict(["x$i" => "e$i" for i = 1:4])

experiments_costs = Dict([
    features_experiments[e] => (i, i) => [e] for (i, e) in enumerate(names(data)[1:4])
])

experiments_costs["ey"] = (100, 100) => ["y"]

# Experimental resources
experimental_resources = [:ER1, :ER2, :ER3]
resources_quantities = [rand(1:3, length(experimental_resources)) for _ = 1:5]

# "Compound," which is a structured token.
# We use `@structured` macro, which is a convenience wrapper around `@aagent`,
# defined in ReactiveDynamics.jl
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct Compound
        state::Any
        history::Vector{Symbol}
    end

    get_cmpds = function (transition::Transition) end
    run_experiment = function (agent::Compound, experiment::Symbol, rng) end
    assign_to_places = function (state::ReactionNetworkProblem, threshold) end
end

# Provide a constructor for `Compound` and define functions that will
# "execute" the experiments.
import ReactiveDynamics: Compound, get_cmpds, assign_to_places, sample
using ReactiveDynamics: ReactionNetworkProblem, Transition

function Compound(id::AbstractString, predictions::Dict)
    state = State((Evidence(predictions...), Tuple(zeros(2))))

    return Compound("Compound $id", :pool, nothing, [], state, String[])
end

using Random: default_rng

ReactiveDynamics.get_cmpds = function (transition::Transition)
    if !isnothing(transition) && !isempty(transition.bound_structured_agents)
        agent_ix = findall(x -> x isa Compound, transition.bound_structured_agents)

        return transition.bound_structured_agents[agent_ix]
    end
end

ReactiveDynamics.run_experiment =
    function (agents::Vector, experiment::Symbol, rng = default_rng())
        println("running experiment $experiment")
        for agent in agents
            push!(agent.history, experiment)
            experiment = String(experiment)

            observation = sampler(
                agent.state["preclinical_state"].evidence,
                getindex(experiments_costs[experiment], 2),
                rng,
            )

            agent.state["preclinical_state"] = merge(
                agent.state["preclinical_state"],
                observation,
                first(experiments_costs[experiment]),
            )
        end
        return agents
    end

transitions_experiments = String[]
for i = 1:5
    experiment = i < 5 ? "e$i" : "ey"

    resource_part = join(
        [
            "$(resources_quantities[i][res_i]) * $res" for
            (res_i, res) in enumerate(experimental_resources)
        ],
        " + ",
    )

    push!(
        transitions_experiments,
        """@deterministic($experiment), $experiment + $resource_part --> @move(:$experiment, :pool),
            action => run_experiment(get_cmpds(@transition), :$experiment)
        """,
    )
end

# Set up a reaction network

network = @ReactionNetworkSchema

# Resource generation part
@push network begin
    p1, ∅ --> ER1
    p2, ∅ --> ER2
    p3, ∅ --> ER3
end

@prob_init network ER1 = 1200 ER2 = 1500 ER3 = 1300

@prob_params network p1 = 1 p2 = 1 p3 = 1

# Experiments part
for species in union(Symbol.(keys(experiments_costs)), [:pool])
    ReactiveDynamics.register_structured_species!(network, species)
end

str_network_def = """
    begin
        @push network begin
            $(join(transitions_experiments, '\n'))
        end
    end 
"""

using MacroTools: striplines
expr_network_def = striplines(Meta.parseall(str_network_def))

eval(expr_network_def)

@prob_meta network tspan = 100

problem = ReactionNetworkProblem(network)

using Random: randstring

# Simplified setup: we assume that compounds are already assigned
# to the "experimental" places.
for _ = 1:10
    cmpd = Compound(randstring(4), Dict())
    cmpd.species = :pool

    add_structured_token!(problem, cmpd)
end

simulate(problem, 10)

# To allow for dynamic assignment of compounds to places, we need to create an agent
# that will move the agent to the corresponding places.
# This can be expressed in two ways:
# - Create an "algebraic agent" which will modify the reaction network's state,
# - Create a "placeholder" transition which will run the mutating function.

# In either case, we need to define the function that will facilitate the assignments.
evidence = Evidence()

solver = GenerativeDesigns.DPWSolver(; n_iterations = 500, tree_in_info = true)
repetitions = 5
mdp_options = (; max_parallel = 1, discount = 1.0, costs_tradeoff = (0.5, 0.5))

ReactiveDynamics.assign_to_places =
    function (state::ReactionNetworkProblem, threshold = 0.1)
        compounds = filter(
            x -> ReactiveDynamics.get_species(x) == :pool,
            collect(values(inners(inners(state)["structured"]))),
        )

        for cmpd in compounds
            e = get_next_experiment(cmpd.state["preclinical_state"].evidence, threshold)
            if !isnothing(e)
                cmpd.species = first(e)
            end
        end
    end

function get_next_experiment(evidence::Evidence, threshold = 0.1)
    design = efficient_design(
        experiments_costs;
        sampler = sampler,
        uncertainty = uncertainty,
        threshold,
        evidence = evidence,
        solver = solver,
        repetitions = repetitions,
        mdp_options = mdp_options,
    )

    if !haskey(design[2], :arrangement) || isempty(design[2].arrangement)
        return nothing
    else
        return first(design[2].arrangement)
    end
end

@push network begin
    @deterministic(2), ∅ --> @structured(Compound(randstring(4), Dict()))
    @deterministic(1), ∅ --> ∅, preAction => assign_to_places(@state)
end

problem = ReactionNetworkProblem(network)

# Simplified setup: we assume that compounds are already assigned
# to the "experimental" places.
for _ = 1:10
    cmpd = Compound(randstring(4), Dict())
    cmpd.species = Symbol(:pool)

    add_structured_token!(problem, cmpd)
end

simulate(problem, 10)
