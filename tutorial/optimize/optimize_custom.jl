using ReactiveDynamics

# fit unknown dynamics to empirical data
## some embedded neural network, etc.
@register begin
    function function_to_learn(A, B, C, params)
        return [A, B, C]' * params # params: 3-element vector
    end
end

acs = @ReactionNetwork begin
    function_to_learn(A, B, C, params), A --> B + C
    1.0, B --> C
    2.0, C --> B
end

# initial values, check params
@prob_init acs A = 60.0 B = 10.0 C = 150.0
@prob_params acs params = [0.01, 0.01, 0.01]
@prob_meta acs tspan = 100.0

prob = @problematize acs
sol = @solve prob

time_points = [1, 50, 100]
data = [60 30 5]

@fit_and_plot acs data time_points vars = [A] params α maxeval = 200 lower_bounds = 0 upper_bounds =
    0.01

# a slightly more extended example: export solver as a function of given params
## some embedded neural network, etc.
@register begin
    function learnt_function(A, B, C, params, α)
        return [A, B, C]' * params + α # params: 3-element vector
    end
end
acs = @ReactionNetwork begin
    learnt_function(A, B, C, params, α), A --> B + C, priority => 0.6
    1.0, B --> C
    2.0, C --> B
end

# initial values, check params
@prob_init acs A = 60.0 B = 10.0 C = 150.0
@prob_params acs params = [1, 2, 3] α = 1
@prob_meta acs tspan = 100.0

prob = @problematize acs
sol = @solve prob
@optimize acs abs(C) params α maxeval = 200 lower_bounds = 0 upper_bounds = 200 final_only =
    true

## export solution as a function of params
parametrized_solver = @build_solver acs A params α trajectories = 10
parametrized_solver = @build_solver acs A params α
parametrized_solver([1.0, 0.0, 0.0, 0.0, 1.0])

using Statistics # for mean
my_little_objective(params) =
    mean(1:10) do _ # take avg over 10 samples
        sol = parametrized_solver(params) # compute solution
        return sol[3, end] # some objective
    end

## now some optimizer...
using NLopt # https://github.com/JuliaOpt/NLopt.jl
obj_for_nlopt = (vec, _) -> my_little_objective(vec) # the second term corresponds to a gradient - won't be used, but is a part of the NLopt interface

opt = Opt(:GN_DIRECT, 5)
opt.lower_bounds = 0
opt.upper_bounds = 100

opt.max_objective = obj_for_nlopt
opt.maxeval = 100

(minf, minx, ret) = optimize(opt, rand(5))
