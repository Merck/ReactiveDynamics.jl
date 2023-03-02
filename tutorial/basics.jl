using ReactiveDynamics

acs = @ReactionNetwork begin
    1.0, X ⟺ Y
end

acs = @ReactionNetwork begin
    1.0, X ⟺ Y, name=>"transition1"
end

@prob_init acs X=10 Y=20
@prob_params acs
@prob_meta acs tspan=250 dt=0.1

prob = @problematize acs

# sol = @solve prob trajectories=50

ReactiveDynamics.solve(prob)