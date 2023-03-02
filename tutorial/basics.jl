using ReactiveDynamics

acs = @ReactionNetwork begin
    1.0, X ⟺ Y
end

acs = @ReactionNetwork begin
    1.0, X ⟺ Y, name=>"transition1"
end