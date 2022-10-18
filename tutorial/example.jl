using ReactiveDynamics

# acs as a model : incomplete dynamics
sir_acs = @ReactionNetwork

# set up ontology: try ?@aka
@aka sir_acs transition=reaction species=population_group
# add single reaction
@push sir_acs β*S*I*tdecay(@t()) S+I --> @choose((1.5, 2*@choose(1, 2)*I), (1.6, @register(:RATE, 1, every=2)*I)) name=>SI2I

# additional dynamics
@push sir_acs begin 
        ν*I, I --> R, name=>I2R
        γ, R --> S, name=>R2S
end

# register custom function
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t/10^3)

@valuation sir_acs I=.1 # for testing purpose only

# atomic data: initial values, parameter values
@prob_init sir_acs S=999 I=100 R=100
u0 = [999, 10, 0] # alternative specification
@prob_init sir_acs u0
@prob_params sir_acs β=0.0001 ν=0.01 γ=5
@prob_meta sir_acs tspan=100

#prob = @problematize sir_acs
prob = @problematize sir_acs tspan=200

sol = @solve prob trajectories=20
@plot sol plot_type=summary
sol = @solve prob
@plot sol plot_type=allocation
@plot sol plot_type=valuations
@plot sol plot_type=new_transitions

# acs as a model : incomplete dynamics
sir_acs = @ReactionNetwork

# set up ontology: try ?@aka
@aka sir_acs transition=reaction species=population_group
# add single reaction

# additional dynamics
@push sir_acs begin 
        ν*I*tdecay(@t()), I --> @choose(E, R, S), name=>I2R
        γ, R --> S, name=>R2S
end
@jump sir_acs 3 (S > 0 && bool_cond(@t()) && (I += 1; S -= 1)) # drift term, support for event conditioning
@periodic sir_acs 5. (β += 1; println(β))
# register custom function
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t/10^3)

# atomic data: initial values, parameter values
@prob_init sir_acs S=999 I=10 R=0
u0 = [999, 10, 0] # alternative specification
@prob_init sir_acs u0
@prob_params sir_acs β=0.0001 ν=0.01 γ=5

# model dynamics
sir_acs = @ReactionNetwork begin
        α*S*I, S+I --> 2I, cycle_time=>0, name=>I2R
        β*I, I --> R, cycle_time=>0, name=>R2S 
end

# simulation parameters
@prob_init sir_acs S=999 I=10 R=0
@prob_uncertainty sir_acs S=10. I=5.
@prob_params sir_acs α=0.0001 β=0.01
@prob_meta sir_acs tspan=250 dt=.1

# batch simulation
prob = @problematize sir_acs
sol = @solve prob trajectories=20
@plot sol plot_type=summary
@plot sol plot_type=summary show=:S
@plot sol plot_type=summary c=:green xlimits=(.0, 100.)

acs_1 = @ReactionNetwork begin
        1., A --> B + C
end
acs_2 = @ReactionNetwork begin
        1., A --> B + C
end

acs_union = @join acs_1 acs_2 acs_1.A=acs_2.A=@alias(A) @catchall(B) # @catchall: equalize species "B" from the joined submodels
@equalize acs_union acs_1.C=acs_2.C=@alias(C) # the species can be set equal either within the call to @join the acsets, but you may always equalize a set of species of any given acset