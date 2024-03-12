using ReactiveDynamics

rd_model = @ReactionNetwork

# need to bring data into ReactiveDynamics's (solver) scope
@register begin
    # number of phases
    n_phase = 2
    n_workforce = 5
    n_equipment = 4
    # transition matrix: ones at the upper diagonal + rand noise
    prob_success = [i / (i + 1) for i = 1:n_phase]
    transitions = zeros(n_phase, n_phase)
    foreach(i -> transitions[i, i+1] = 1, 1:(n_phase-1))
    transitions .+= 0.1 * rand(n_phase, n_phase)
end

rd_model = @ReactionNetworkSchema begin
    $i, phase[$i] --> phase[$j], cycle_time => $i * $j
end i = 1:3 j = 1:($i)

using ReactiveDynamics: nparts
u0 = rand(1:1000, nparts(rd_model, :S))
@prob_init rd_model u0

@prob_meta rd_model tspan = 100

prob = @problematize rd_model
sol = @solve prob trajectories = 20
@plot sol plot_type

# need to bring data into ReactiveDynamics's (solver) scope
@register begin
    k = 5
    r = 3
    M = zeros(k, k)
    for i = 1:k
        for conn_in in rand(1:k, rand(1:4))
            M[conn_in, i] += 0.1 * rand(1:10)
        end
        for conn_out in rand(1:k, rand(1:4))
            M[i, conn_out] += 0.1 * rand(1:10)
        end
    end
    cycle_times = rand(1:5, k, k)
    resource = rand(1:10, k, k, r)
end

rd_model = @ReactionNetworkSchema begin
    M[$i, $j],
    mod[$i] +
    {resource[$i, $j, $k] * resource[$k], k = rand(1:(ReactiveDynamics.r)), dlm = +} -->
    mod[$j],
    cycle_time => cycle_times[$i, $j]
end i = 1:(ReactiveDynamics.k) j = 1:(ReactiveDynamics.k)

using ReactiveDynamics: nparts
u0 = rand(1:1000, nparts(rd_model, :S))
@prob_init rd_model u0

@prob_meta rd_model tspan = 100

prob = @problematize rd_model
sol = @solve prob
@plot sol plot_type = allocation
