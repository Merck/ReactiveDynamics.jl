using ReactionDynamics

# solve for steady state
acss = @ReactionNetwork begin
    3.0, A --> A, priority=>.6, name=>aa
    1.0, B + .2*A --> 2*α*B, prob=>.7, priority=>.6, name=>bb
    3.0, A + 2*B--> 2*C, prob=>.7, priority=>.7, name=>cc
end

# initial values, check params
@prob_init acss A=100. B=200 C=150.
@prob_params acss α=10
@prob_meta acss tspan=100.

prob = @problematize acss
sol = @solve prob trajectories=10

# check https://github.com/JuliaOpt/NLopt.jl for solver opts
@optimize acss abs(A-B) A B=20. α=2. lower_bounds=0 upper_bounds=200 min_t=50 max_t=100 maxeval=20 final_only=true

t_ = [1, 50, 100]
data = [70 50 50]

@fit acss data t_ vars=[A] B=20 A α lower_bounds=0 upper_bounds=200 maxeval=20

@fit_and_plot acss data t_ vars=[A] B=20 A α lower_bounds=0 upper_bounds=200 maxeval=2 trajectories=10