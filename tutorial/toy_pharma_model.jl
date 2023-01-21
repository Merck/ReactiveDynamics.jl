using ReactiveDynamics

# model dynamics
toy_pharma_model = @ReactionNetwork begin
    α(candidate_compound, marketed_drug, κ),
    3 * @conserved(scientist) + @rate(budget) --> candidate_compound, name => discovery,
    probability => 0.3, cycletime => 10.0, priority => 0.5
    β(candidate_compound, marketed_drug),
    candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) -->
    marketed_drug + 5 * budget, name => dx2market, probability => 0.5 + 0.001 * @t(),
    cycletime => 4
    γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
end

@periodic toy_pharma_model 1.0 budget+=11 * marketed_drug

@register function α(number_candidate_compounds, number_marketed_drugs, κ)
    κ + exp(-number_candidate_compounds) + exp(-number_marketed_drugs)
end
@register function β(number_candidate_compounds, number_marketed_drugs)
    number_candidate_compounds + exp(-number_marketed_drugs)
end

# simulation parameters
## initial values
@prob_init toy_pharma_model candidate_compound=5 marketed_drug=6 scientist=20 budget=100
## parameters
@prob_params toy_pharma_model κ=4 γ=0.1
## other arguments passed to the solver
@prob_meta toy_pharma_model tspan=250 dt=0.1

prob = @problematize toy_pharma_model

sol = @solve prob trajectories=20

using Plots

@plot sol plot_type=summary

@plot sol plot_type=summary show=:marketed_drug

## for deterministic rates 

# model dynamics
toy_pharma_model = @ReactionNetwork begin
    @per_step(α(candidate_compound, marketed_drug, κ)),
    3 * @conserved(scientist) + @rate(budget) --> candidate_compound, name => discovery,
    probability => 0.3, cycletime => 10.0, priority => 0.5
    @per_step(β(candidate_compound, marketed_drug)),
    candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) -->
    marketed_drug + 5 * budget, name => dx2market, probability => 0.5 + 0.001 * @t(),
    cycletime => 4
    @per_step(γ*marketed_drug), marketed_drug --> ∅, name => drug_killed
end

@periodic toy_pharma_model 0.0 budget+=11 * marketed_drug

@register function α(number_candidate_compounds, number_marketed_drugs, κ)
    κ + exp(-number_candidate_compounds) + exp(-number_marketed_drugs)
end
@register function β(number_candidate_compounds, number_marketed_drugs)
    number_candidate_compounds + exp(-number_marketed_drugs)
end

# simulation parameters
## initial values
@prob_init toy_pharma_model candidate_compound=5 marketed_drug=6 scientist=20 budget=100
## parameters
@prob_params toy_pharma_model κ=4 γ=0.1
## other arguments passed to the solver
@prob_meta toy_pharma_model tspan=250

prob = @problematize toy_pharma_model

sol = @solve prob trajectories=20

using Plots

@plot sol plot_type=summary

@plot sol plot_type=summary show=:marketed_drug
