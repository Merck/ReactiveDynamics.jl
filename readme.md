# ReactiveDynamics.jl <br> 

<p align="center">
  <img src="docs/assets/diagram1.png" alt="wiring diagram"> <br>
  <a href="#about">About</a> |
  <a href="#context-dynamics-of-value-evolution-dyve">Context</a> |
  <a href="#three-sketches">Three Sketches</a> |
  <a href="https://merck.github.io/ReactiveDynamics.jl">Documentation</a>
</p>

## About

The package provides a category of reaction (transportation) network-type problems formalized on top of the **[generalized algebraic theory](https://ncatlab.org/nlab/show/generalized+algebraic+theory)**, and is compatible with the **[SciML](https://sciml.ai/)** ecosystem.

Our motivation stems from the area of **[system dynamics](https://www.youtube.com/watch?v=o-Yp8A7BPE8)**, which is a mathematical modeling approach to frame, analyze, and optimize complex (nonlinear) dynamical systems, to augment the strategy and policy design.

<img src="docs/assets/diagram2.png" align="right" alt="wiring diagram"></a>
<p>The central concept is of a <b>transition</b> (transport, flow, rule, reaction - a subcategory of general algebraic action). Generally, a transition prescribes a stochastic rule which repeatedly transforms the modeled system's <b>resources</b>. An elementary instance of such an ontology is provided by chemical reaction networks.

A <b>reaction network</b> (system modeled) is then a tuple $(T, R)$, where $T$ is a set of the transitions and $R$ is a set of the network's resource classes (aka species). The simultaneous action of transitions on the resources evolves the dynamical system.

The transitions are generally **stateful** (i.e., act over a period of time). Moreover, at each time step a quantum of the a transition's instances is brought into the scope, where the size of the batch is drived by a Poisson counting process. A transition takes the from `rate, a*A + b*B + ... --> c*C + ...`, where `rate` gives the expected batch size per time unit. `A`, `B`, etc., are the resources, and `a`, `b`, etc., are the generalized stoichiometry coefficients. Note that both `rate` and the "coefficients" can in fact be given by a function which depends on the system's instantaneous state (stochastic, in general). In particular, even the structural form of a transition can be stochastic, as will be demonstrated shortly.
</p>

An instance of a stateful transition evolves gradually from its genesis up to the terminal point, where the products on the instance's right hand-side are put into the system. An instance advances proportionally to the quantity of resources allocated. To understand this behavior, we note that the system's resource classes (or occurences of a resource class in a transition) may be assigned a **modality**; a modality governs the interaction between the resource and the transition, as well as the interpretation of the generalized stoichiometry coefficient.

In particular, a resource can be allocated to an instance either for the instance's lifetime or a single time step of the model's evolution (after each time step, the resource may be reallocated based on the global demand). Moreover, the coefficient of a resource on the left hand-side can either be interpreted as a total amount of the resource required or as an amount required per time unit. Similarly, it is possible to declare a resource **conserved**, in which case it is returned into the scope once the instance terminates.

<img src="docs/assets/diagram3.png" align="left" alt="attributes diagram"></a>

The transitions are <b>parametric</b>. That is, it is possible to set the period over which an instance of a transition acts in the system (as well as the maximal period of this action), the total number of transition's instances allowed to exist in the system, etc. An annotated transition takes the form `rate, a*A + b*B + ... --> c*C + ..., prm => val, ...`, where the numerical values can be given by a function which depends on the system's state. Internally, the reaction network is represented as an <a href=https://algebraicjulia.github.io/Catlab.jl/dev/generated/wiring_diagrams/wd_cset/><b>attributed C-set</b></a>.

A network's dynamics is specified using a compact **modeling metalanguage**. Moreover, we have integrated another expression comprehension metalanguage which makes it easy to generate arbitrarily complex dynamics from a single template transition!

Taking **unions** of reaction networks is fully supported, and it is possible to identify the resource classes as apropriate.

Moreover, it is possible to **export and import** reaction network dynamics using the [TOML](https://toml.io/) format.

Once a network's dynamics is specified, it can be converted to a problem and simulated. The exported problem is a **`DiscreteProblem`** compatible with **[DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)** ecosystem, and hence the latter package's all powerful capabilities are available. For better user experience, we have tailored and exported many of the functionalities within the modeling metalanguage, including ensemble analysis, parameter optimization, parameter inference, etc. Would you guess that **[universal differential equations](https://arxiv.org/abs/2001.04385)** are supported? If even the dynamics is unknown, you may just infer it!

## Context: Dynamics of Value Evolution (DyVE)
 
The package is an integral part of the **Dynamics of Value Evolution (DyVE)** computational framework for learning, designing, integrating, simulating, and optimizing R&D process models, to better inform strategic decisions in science and business.
 
As the framework evolves, multiple functionalities have matured enough to become standalone packages.
 
This includes **[GeneratedExpressions.jl](https://github.com/Merck/GeneratedExpressions.jl)**, a metalanguage to support code-less expression comprehensions. In the present context, expression comprehensions are used to generate complex dynamics from user-specified template transitions.
 
Another package is **[AlgebraicAgents.jl](https://github.com/Merck/AlgebraicAgents.jl)**, a lightweight package to enable hierarchical, heterogeneous dynamical systems co-integration. It implements a highly scalable, fully customizable interface featuring sums and compositions of dynamical systems. In present context, we note it can be used to co-integrate a reaction network problem with, e.g., a stochastic ordinary differential problem!

## Three Sketches

For other examples, see the **[tutorials](tutorial)**.

### SIR Model

The acronym SIR stands for susceptible, infected, and recovered, and as such the SIR model attempts to capture the dynamics of disease spread. We express the SIR dynamics as a reaction network using the compact modeling metalanguage. 

Follow the SIR model's reactions:

<p align="center">
  <img src="docs/assets/sir_reactions.png" alt="SIR reactions"> <br>
</p>

```julia
using ReactiveDynamics

# model dynamics
sir_acs = @ReactionNetwork begin
        α*S*I, S+I --> 2I, name=>I2R
        β*I, I --> R, name=>R2S 
end

# simulation parameters
## initial values
@prob_init sir_acs S=999 I=10 R=0
## uncertainty in initial values (Gaussian)
@prob_uncertainty sir_acs S=10. I=5.
## parameters
@prob_params sir_acs α=0.0001 β=0.01
## other arguments passed to the solver
@prob_meta sir_acs tspan=250 dt=.1
```

The resulting reaction network is represented as an attributed C-set:

![sir acs](docs/assets/sir_acs.png)

Next we solve the problem.

```
# turn model into a problem
prob = @problematize sir_acs

# solve the problem over multiple trajectories
sol = @solve prob trajectories=20

# plot the solution
@plot sol plot_type=summary
## show only species S
@plot sol plot_type=summary show=:S
## plot evolution over (0., 100.) in green (propagates to Plots.jl)
@plot sol plot_type=summary c=:green xlimits=(.0, 100.)
```

![sir plots](docs/assets/sir_plot.png)

### Sparse Interactions

We introduce a complex reaction network as a union of $n_{\text{models}}$ reaction networks, where the off-diagonal interactions are sparse.

To harness the capabilities of **GeneratedExpressions.jl**, let us first declare a template atomic model.

```julia
# submodel.jl
# substitute $r as the global number of resources, $i as the submodel identifier
@register begin 
        push!(ns, rand(1:5)); ϵ = 10e-2
        push!(M, rand(ns[$i], ns[$i])); foreach(i -> M[$i][i, i] += ϵ, 1:ns[$i])
        foreach(i -> M[$i][i, :] /= sum(M[$i][i, :]), 1:ns[$i])
        push!(cycle_times, rand(1:5, ns[$i], ns[$i]))
        push!(demand, rand(1:10, ns[$i], ns[$i], $r)); push!(production, rand(1:10, ns[$i], ns[$i], $r))
end
    
# generate submodel dynamics
push!(rd_models, @ReactionNetwork begin
                M[$i][$m, $n], state[$m] + {demand[$i][$m, $n, $l]*resource[$l], l=1:$r, dlm=+} --> state[$n] + 
                        {production[$i][$m, $n, $l]*resource[$l], l=1:$r, dlm=+}, cycle_time=>cycle_times[$i][$m, $n], probability_of_success=>$m*$n/(n[$i])^2
        end m=1:ReactiveDynamics.ns[$i] n=1:ReactiveDynamics.ns[$i]
)
```

The next step is to instantiate the atomic models (submodels).

```julia
using ReactiveDynamics
## setup the environment
rd_models = ReactiveDynamics.ReactionNetwork[] # submodels

# needs to live within ReactiveDynamics's scope
# the arrays will contain the submodel
@register begin
    ns = Int[] # size of submodels
    M = Array[] # transition intensities of submodels
    cycle_times = Array[] # cycle times of transitions in submodels
    demand = Array[]; production = Array[] # resource production / generation for transitions in submodels
end
```

Load $n_{\text{models}}$ the atomic models (substituting into `submodel.jl`):

```julia
n_models = 5; r = 2 # number of submodels, resources

# submodels: dense interactions
@generate {@fileval(submodel.jl, i=$i, r=r), i=1:n_models}
```

We take union of the atomic models, and we identify common resources.

```julia
# batch join over the submodels
rd_model = @generate "@join {rd_models[\$i], i=1:n_models, dlm=' '}"

# identify resources
@generate {@equalize(rd_model, @alias(resource[$j])={rd_models[$i].resource[$j], i=1:n_models, dlm=:(=)}), j=1:r}
```

Next step is to add some off-diagonal interactions.

```julia
# sparse off-diagonal interactions, sparse declaration
# again, we use GeneratedExpressions.jl
sparse_off_diagonal = zeros(sum(ReactiveDynamics.ns), sum(ReactiveDynamics.ns))
for i in 1:n_models
    j = rand(setdiff(1:n_models, (i, )))
    i_ix = rand(1:ReactiveDynamics.ns[i]); j_ix = rand(1:ReactiveDynamics.ns[j])
    sparse_off_diagonal[i_ix+sum(ReactiveDynamics.ns[1:i-1]), j_ix+sum(ReactiveDynamics.ns[1:j-1])] += 1
    interaction_ex = """@push rd_model begin 1., var"rd_models[$i].state[$i_ix]" --> var"rd_models[$j]__state[$j_ix]" end"""
    eval(Meta.parseall(interaction_ex))
end
```

Let's plot the interactions:

![interactions](docs/assets/interactions.png)

The resulting model can then be conveniently simulated using the modeling metalanguage.

```julia
using ReactiveDynamics: nparts
u0 = rand(1:1000, nparts(rd_model, :S))
@prob_init rd_model u0

@prob_meta rd_model tspan=10

prob = @problematize rd_model
sol = @solve prob trajectories=2

# plot "state" species only
@plot sol plot_type=summary show=r"state"
```

### Universal Differential Equations: Fitting Unknown Dynamics

We demonstrate how to fit unknown part of dynamics to empirical data.

We use `@register` to define a simple linear function within the scope of module `ReactiveDynamics`; parameters of the function will be then optimized for. Note that `function_to_learn` can be generally replaced with a neural network (Flux chain), etc.

```julia
## some embedded function (neural network, etc.)
@register begin
    function function_to_learn(A, B, C, params)
        [A, B, C]' * params # params: 3-element vector
    end
end
```

Next we set up a simple dynamics and supply initial parameters.

```julia
acs = @ReactionNetwork begin
    function_to_learn(A, B, C, params), A --> B+C
    1., B --> C
    2., C --> B
end

# initial values, check params
@prob_init acs A=60. B=10. C=150.
@prob_params acs params=[.01, .01, .01]
@prob_meta acs tspan=100.
```

Let's next see the numerical results for the initial guess.

```julia
sol = @solve acs
@plot sol
```

![plot](docs/assets/optim1.png)

Next we supply empirical data and fit `params`.

```julia
time_points = [1, 50, 100]
data = [60 30 5]

@fit_and_plot acs data time_points vars=[A] params α maxeval=200 lower_bounds=0 upper_bounds=.01
```

![plot](docs/assets/optim2.png)