acs = @ReactionNetwork begin
    3.0, A --> A, priority=>.6, name=>aa
    1.0, B + .2*A --> 2*α*B, prob=>.7, priority=>.6, name=>bb
    3.0, A + 2*B--> 2*C, prob=>.7, priority=>.7, name=>cc
end

# initial values, check params
@prob_init acs A=60. B=10. C=150.
@prob_params acs α=10.
@prob_meta acs tspan=100.

prob = @problematize acs
sol = @solve prob

# check https://github.com/JuliaOpt/NLopt.jl for solver opts
@optimize acs abs(A-B) A B=20. α=2. lower_bounds=0 upper_bounds=200 min_t=50 max_t=100 maxeval=1000 final_only=true

t = [1, 50, 100]
data = [50 50 50]

@fit acs data t vars=[A] B=20 A α lower_bounds=0 upper_bounds=100 maxeval=1000

@fit_and_plot acs data t vars=[A] B=20 A α lower_bounds=0 upper_bounds=100 maxeval=1000