import Pkg

# Manually add dependencies
#=
Pkg.activate(".")
Pkg.develop(path = "../..")
Pkg.add(["MacroTools", "AlgebraicAgents", "JSON3", "Distributions"])
=#

# --------------------------------------------------------------------------------
# Parametrize the reaction network

# We will have two classes of entities, each of which can be in a number of different states.
# States "S" and "F" will internally denote some terminal, "success" and "failure" states, respectively.

terminal_states = ["S", "F"]
M1_states = "M1_" .* ["A", "B", "C", "D", terminal_states...]
M2_states = "M2_" .* ["A", "B", "C", "D", "E", terminal_states...]

experimental_resources = ["ER1", "ER2", "ER3"]

M1_transition_probs = round.(rand(length(M1_states), length(M1_states)); digits=2)
M2_transition_probs = round.(rand(length(M2_states), length(M2_states)); digits=2)

# Ensure forward flow between states and set zero transition prob to terminal states
for probs in [M1_transition_probs, M2_transition_probs]
    for i in axes(probs, 1)
        probs[i, i:end] .= 0
    end

    probs[:, end-1] .= 0
end

M1_priorities = round.(rand(length(M1_states)); digits=2)
M2_priorities = 2 * round.(rand(length(M2_states)); digits=2)

M1_resources = [rand(1:5, length(experimental_resources)) for _ in M1_states]
M2_resources = [2 * rand(1:5, length(experimental_resources)) for _ in M2_states]

M1_durations = rand(2:5, length(M1_states))
M2_durations = rand(2:5, length(M2_states))

# --------------------------------------------------------------------------------
# Add initial quantities

M1_initial = [rand(1:10, length(M1_states)-2)..., 0, 0]
M2_initial = [rand(1:10, length(M2_states)-2)..., 0, 0]

experimental_resources_initial = rand(1e2:3e2, size(experimental_resources))

# --------------------------------------------------------------------------------
# Export to JSON

# First, we build a dictionary.
data = Dict("species" => [], "transitions" => [])

species, transitions = data["species"], data["transitions"]

# --------------------------------------------------------------------------------
# Add species with initial quantities.
# Later, we may add a couple more attributes (conserved resources, etc.).

for (state, q) in zip(M1_states, M1_initial)
    push!(species, Dict("name" => "$state", "initial" => q))
end

for (state, q) in zip(M2_states, M2_initial)
    push!(species, Dict("name" => "$state", "initial" => q))
end

for (res, q) in zip(experimental_resources, experimental_resources_initial)
    push!(species, Dict("name" => "$res", "initial" => q))
end

# --------------------------------------------------------------------------------
# Add transitions for entity class M1.

for (i, state) in enumerate(M1_states[begin:end-2])
    t = Dict{String, Any}("priority" => M1_priorities[i], "duration" => M1_durations[i], "rate" => state, "name" => "M1")

    # from (left-hand side)
    push!(t, "from" => [Dict("name" => res, "q" => q) for (res, q) in zip(experimental_resources, M1_resources[i])])
    push!(t["from"], Dict("name" => state, "q" => 1))

    # to (right-hand side)
    push!(t, "to" => [])
    to = t["to"]

    for (j, state2) in enumerate(M1_states[i:end])
        if M1_transition_probs[j+i-1, i] > 0
            push!(to, Dict("probability" => M1_transition_probs[j+i-1, i], "reactants" => [Dict("name" => state2, "q" => 1)]))
        end
    end

    push!(transitions, t)
end

# Add transitions for entity class M2.

for (i, state) in enumerate(M2_states[begin:end-2])
    t = Dict{String, Any}("priority" => M2_priorities[i], "duration" => M2_durations[i], "rate" => state, "name" => "M2")

    # from (left-hand side)
    push!(t, "from" => [Dict("name" => res, "q" => q) for (res, q) in zip(experimental_resources, M2_resources[i])])
    push!(t["from"], Dict("name" => state, "q" => 1))

    # to (right-hand side)
    push!(t, "to" => [])
    to = t["to"]

    for (j, state2) in enumerate(M2_states[i:end])
        if M2_transition_probs[j+i-1, i] > 0
            push!(to, Dict("probability" => M2_transition_probs[j+i-1, i], "reactants" => [Dict("name" => state2, "q" => 1)]))
        end
    end

    push!(transitions, t)
end

using JSON3
open("reaction_network.json", "w") do io
    #JSON3.pretty(io, data)
end

# --------------------------------------------------------------------------------
# Import from JSON
# Eventually move into the module (refactor current CSV interface).

# Species and initial values.
name = "reaction_network"
str_init = "@prob_init $name"
for s in data["species"]
    str_init *= " $(s["name"])=$(s["initial"])"
end

function get_subline(d::Vector)
    if isempty(d)
        return "∅"
    elseif haskey(first(d), "probability")
        sublines = ["($(sd["probability"]), " * join(["$(s["q"]) * $(s["name"])" for s in sd["reactants"]], " + ") * ")" for sd in d]
        return "@choose(" * join(sublines, ", ") * ")"
    else
        return join(["$(s["q"]) * $(s["name"])" for s in d], " + ")
    end
end

# Transitions.
str_transitions = []
for t in data["transitions"]
    line = t["rate"] * ", " * get_subline(t["from"]) * " --> " * get_subline(t["to"])
    line *= ", name => $(t["name"]), priority => $(t["priority"]), cycletime => $(t["duration"])"
    push!(str_transitions, line)
end

push!(str_transitions, "i1, ∅ --> M1_A, name => M1_creation")
push!(str_transitions, "i2, ∅ --> M2_A, name => M2_creation")

str_params = "@prob_params $name i1 = .3 i2 = .2"

str_network_def = """
    begin
        $name = @ReactionNetworkSchema
        @push $name begin
            $(join(str_transitions, '\n'))
        end
        $str_init
        $str_params
    end
"""

using MacroTools: striplines
expr_network_def = striplines(Meta.parseall(str_network_def))

using ReactiveDynamics

eval(expr_network_def)

@isdefined reaction_network

# --------------------------------------------------------------------------------
# Solve problem

@prob_meta reaction_network tspan = 100 dt = 1.0

# Convert network into an AlgAgents hierarchy.
problem = ReactionNetworkProblem(reaction_network; name = "network")

# AlgAgents: "periodic" callback.
using AlgebraicAgents

@aagent struct Controller
    M1_λ::Float64
    M2_λ::Float64

    log::Vector{String}

    current_time::Float64
    time_step::Float64

    initial_time::Float64
end

function Controller(name::String, M1_λ::Float64, M2_λ::Float64, time::T, time_step::T) where {T<:Real}
    return Controller(name, M1_λ, M2_λ, String[], time, time_step, time)
end

using Distributions: Poisson

function AlgebraicAgents._step!(c::Controller)
    n_removed_M1 = rand(Poisson(c.time_step * c.M1_λ))
    n_removed_M2 = rand(Poisson(c.time_step * (c.M2_λ + c.M1_λ)))

    transitions = getagent(c, "../network").ongoing_transitions
    M1_transitions = filter(x -> x[:transName] === :M1, transitions)
    M2_transitions = filter(x -> x[:transName] === :M2, transitions)

    M1_transitions_delete = isempty(M1_transitions) ? [] : unique(rand(M1_transitions, n_removed_M1))
    for trans in M1_transitions_delete
        trans.state = trans[:transCycleTime]
        trans.trans[:transProbOfSuccess] = 0
    end

    M2_transitions_delete = isempty(M2_transitions) ? [] : unique(rand(M2_transitions, n_removed_M2))
    for trans in M2_transitions_delete
        trans.state = trans[:transCycleTime]
        trans.trans[:transProbOfSuccess] = 0
    end

    push!(c.log, "t = $(c.current_time) removed compounds: " * join(getname.(union(M1_transitions_delete, M2_transitions_delete)), ", "))

    return c.current_time += c.time_step
end

AlgebraicAgents._projected_to(c::Controller) = c.current_time

c = Controller("controller", 1e-1, 2e-1, 0., 1.);

compound_problem = ⊕(problem, c; name = "compound problem");

# Simulate
simulate(compound_problem, 100);

# Access solution
compound_problem.inners["network"].sol
compound_problem.inners["controller"].log

# Plot solution
draw(compound_problem.inners["network"])

# --------------------------------------------------------------------------------
# Add resource making part to the reaction network

eval(expr_network_def)

n_primary_resources = 5
primary_resources = ["R$i" for i in 1:n_primary_resources]

str_resource_making_transitions = [
    "p_primary_$i, ∅ --> $res, name => primary_resource_maker_$i" for (i, res) in enumerate(primary_resources)
] 

for (i, res) in enumerate(experimental_resources)
    stoich = rand(1:5, n_primary_resources)
    rate = "p_$i * (" * join(primary_resources, " + ") * ")"
    transition_lhs = join(["$(stoich[i])" * "$primary_res" for (i, primary_res) in enumerate(primary_resources)], " + ")
    
    push!(str_resource_making_transitions, "$rate, $transition_lhs --> $res, name => resource_maker_$i")
end

str_resource_making_params = 
    "@prob_params resource_making_network " *
    join(["p_$i = $(rand(1:3))" for i in 1:length(experimental_resources)], " ") * " " *
    join(["p_primary_$i = $(rand(1:5))" for i in 1:length(primary_resources)], " ")

str_network_def = """
    begin
        resource_making_network = @ReactionNetworkSchema
        @push resource_making_network begin
            $(join(str_resource_making_transitions, '\n'))
        end
        $str_resource_making_params
    end
"""

expr_network_resource_making_def = striplines(Meta.parseall(str_network_def))

eval(expr_network_resource_making_def)

extended_network = union_acs!(reaction_network, resource_making_network, "resource_making_network")

#equalize!(extended_network, [Meta.parseall("$res=resource_making_network.$res") for res in experimental_resources])
str_equalize = "@equalize extended_network " * join(["$res=resource_making_network.$res" for res in experimental_resources], " ")

equalize_expr = striplines(Meta.parseall(str_equalize))
eval(equalize_expr)

@prob_meta extended_network tspan = 100 dt = 1.0

# Convert network into an AlgAgents hierarchy.
extended_problem = ReactionNetworkProblem(extended_network; name = "extended_network")

simulate(extended_problem, 100);