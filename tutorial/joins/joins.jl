using ReactiveDynamics
## setup the environment
n_models = 5;
r = 2; # number of submodels, resources
rd_models = ReactiveDynamics.ReactionNetworkSchema[] # submodels

@register begin
    ns = Int[] # size of submodels
    M = Array[] # transition intensities of submodels
    cycle_times = Array[] # cycle times of transitions in submodels
    demand = Array[]
    production = Array[] # resource production / generation for transitions in submodels
end

# submodels: dense interactions
@generate {@fileval(submodel.jl, i = $i, r = r), i = 1:n_models}

# batch join over the submodels
rd_model = @generate "@join {rd_models[\$i], i=1:n_models, dlm=' '}"

# identify resources
@generate {
    @equalize(
        rd_model,
        @alias(resource[$j]) = {rd_models[$i].resource[$j], i = 1:n_models, dlm = :(=)}
    ),
    j = 1:r,
}

# sparse off-diagonal interactions, sparse declaration
sparse_off_diagonal = zeros(sum(ReactiveDynamics.ns), sum(ReactiveDynamics.ns))
for i = 1:n_models
    j = rand(setdiff(1:n_models, (i,)))
    i_ix = rand(1:ReactiveDynamics.ns[i])
    j_ix = rand(1:ReactiveDynamics.ns[j])
    sparse_off_diagonal[
        i_ix+sum(ReactiveDynamics.ns[1:(i-1)]),
        j_ix+sum(ReactiveDynamics.ns[1:(j-1)]),
    ] += 1
    interaction_ex = """@push rd_model begin 1., var"rd_models[$i].state[$i_ix]" --> var"rd_models[$j]__state[$j_ix]" end"""
    eval(Meta.parseall(interaction_ex))
end

sparse_off_diagonal += cat(ReactiveDynamics.M...; dims = (1, 2))
using Plots;
heatmap(1 .- sparse_off_diagonal; color = :greys, legend = false);

using ReactiveDynamics: nparts
u0 = rand(1:1000, nparts(rd_model, :S))
@prob_init rd_model u0

@prob_meta rd_model tspan = 10

prob = ReactionNetworkProblem(sir_acs)
sol = simulate(prob)
draw(sol)
