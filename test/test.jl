using ReactiveDynamics

acs = @ReactionNetwork begin
    α * S * I, S + I --> 2I, cycle_time => 0, name => I2R
    β * I, I --> R, cycle_time => 0, name => R2S
end

acs = @ReactionNetwork begin
    α * S * I, S + a"C/I" --> 2 * MyResource("I", randn()), cycle_time => 0, name => I2R
    β * I, I --> R, cycle_time => 0, name => R2S
end

reaction_network = @ReactionNetwork "network"

@transitions reaction_network begin
    α * S * I, S + a"C/I" --> 2 * MyResource("I", randn()), cycle_time => 0, name => I2R
    β * I, I --> R, cycle_time => 0, name => R2S
end

@species reaction_network "X" "Y" "Z"

@events reaction_network begin (X += 1), name => "xinc" end

@sampleables reaction_network begin [X, a"C/I"], name => "obs", every => 1 end
