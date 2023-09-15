using ReactiveDynamics, AlgebraicAgents

# acs = @ReactionNetwork begin
#     1.0, X ⟺ Y
# end

acs = @ReactionNetwork begin
    1.0, X --> Y, name => "transition1"
end

@prob_init acs X = 10 Y = 20
@prob_params acs
@prob_meta acs tspan = 250 dt = 0.1


# sol = ReactiveDynamics.solve(prob)

#sol = @solve prob

prob = @agentize acs

for i in 1:30
    AlgebraicAgents.step!(prob)
end
prob.history_u
using Plots

@plot sol plot_type = summary
∪

prob = @problematize acs
@solve prob
