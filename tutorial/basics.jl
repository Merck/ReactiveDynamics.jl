using ReactiveDynamics

# define the network
acs = @ReactionNetworkSchema begin
    1.0, X --> Y, name => "transition1"
end

@prob_init acs X = 10 Y = 20
@prob_params acs
@prob_meta acs tspan = 25 dt = 0.10

# convert network into an AlgAgents hierarchy
prob = ReactionNetworkProblem(acs)

# simulate
simulate(prob)

# access solution
prob.sol

draw(prob)
