# substitute $r as the global number of resources, $i as the submodel identifier
@register begin
    push!(ns, rand(1:5))
    ϵ = 10e-2
    push!(M, rand(ns[$i], ns[$i]))
    foreach(i -> M[$i][i, i] += ϵ, 1:ns[$i])
    foreach(i -> M[$i][i, :] /= sum(M[$i][i, :]), 1:ns[$i])
    push!(cycle_times, rand(1:5, ns[$i], ns[$i]))
    push!(demand, rand(1:10, ns[$i], ns[$i], $r))
    push!(production, rand(1:10, ns[$i], ns[$i], $r))
end

# generate submodel dynamics
push!(
    rd_models,
    @ReactionNetwork begin
        M[$i][$m, $n],
        state[$m] + {demand[$i][$m, $n, $l] * resource[$l], l = 1:($r), dlm = +} -->
        state[$n] + {production[$i][$m, $n, $l] * resource[$l], l = 1:($r), dlm = +},
        cycle_time => cycle_times[$i][$m, $n],
        probability_of_success => $m * $n / (n[$i])^2
    end m = 1:ReactiveDynamics.ns[$i] n = 1:ReactiveDynamics.ns[$i]
)
