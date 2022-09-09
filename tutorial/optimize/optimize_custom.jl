using DyVE

## some embedded neural network, etc.
@register begin
    function function_to_learn(A, B, C, params)
        [A, B, C]' * params # params: 3-element vector
    end
end

acs = @ReactionNetwork begin
    function_to_learn(A, B, C, params), A --> B+C
    1., B --> C
    2., C --> B
end

# initial values, check params
@prob_init acs A=60. B=10. C=150.
@prob_params acs params=[.01, .01, .01]
@prob_meta acs tspan=100.

prob = @problematize acs
sol = @solve prob

time_points = [1, 50, 100]
data = [60 30 5]

@fit_and_plot acs data time_points vars=[A] params α maxeval=200 lower_bounds=0 upper_bounds=.01

## export solution as a function of params
parametrized_solver = @build_solver acs A params α trajectories=10
parametrized_solver = @build_solver acs A params α
parametrized_solver([1., 0., 0., 0., 1.])

using Statistics # for mean
function my_little_objective(params)
    mean(1:10) do _ # take avg over 10 samples
        sol = parametrized_solver(params) # compute solution
        sol[3, end] # some objective
    end
end

## now some optimizer...
using NLopt # https://github.com/JuliaOpt/NLopt.jl
obj_for_nlopt = (vec, _) -> my_little_objective(vec) # the second term corresponds to a gradient - won't be used, but is a part of the NLopt interface

opt = Opt(:GN_DIRECT, 5)
opt.lower_bounds = 0
opt.upper_bounds = 100

opt.max_objective = obj_for_nlopt
opt.maxeval = 100

(minf,minx,ret) = optimize(opt, rand(5))