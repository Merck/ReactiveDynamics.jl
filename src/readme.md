# Interface documentation

## Flowchart

![macros_overview](macros_overview/macros_small.png)

<a id='DyVE-Interface-Documentation'></a>

<a id='DyVE-Interface-Documentation-1'></a>

## DyVE Interface Documentation

<a id='DyVE.@ReactionNetwork' href='#DyVE.@ReactionNetwork'>#</a>
**`DyVE.@ReactionNetwork`** &mdash; *Macro*.



Macro that takes an expression corresponding to a reaction network and outputs an instance of `TheoryReactionNetwork` that can be converted to a `DEProblem` or solved directly.

Most arrows accepted (both right, left, and bi-drectional arrows). Use 0 or ∅ for annihilation/creation to/from nothing.

Custom functions and sampleable objects can be used as numeric parameters. Note that these have to be accessible from DyVE's source code.

**Examples**

```julia
acs = @ReactionNetwork begin
    1.0, X ⟶ Y
    1.0, X ⟶ Y, priority=>6., prob=>.7, capacity=>3.
    1.0, ∅ --> (Poisson(.3γ)X, Poisson(.5)Y)
    (XY > 100) && (XY -= 1)
end
@push acs 1.0 X ⟶ Y 
@prob_init acs X=1 Y=2 XY=α
@prob_params acs γ=1 α=4
@solve_and_plot acs
```

<a id='DyVE.@add_species-Tuple{Any, Vararg{Any}}' href='#DyVE.@add_species-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@add_species`** &mdash; *Macro*.



Add new species to a model.

**Examples**

```julia
@add_species acs S I R
```

<a id='DyVE.@aka-Tuple{Any, Vararg{Any}}' href='#DyVE.@aka-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@aka`** &mdash; *Macro*.



Alias object name in an acs.

Followed with an ACS*name, transition=transition*alias species=species_alias.

**Default names**

| name       | short name |
|:---------- |:---------- |
| species    | S          |
| transition | T          |
| action     | A          |
| event      | E          |
| param      | P          |
| meta       | M          |

**Examples**

```julia
@aka acs species=resource transition=reaction
```

<a id='DyVE.@cost-Tuple{Any, Vararg{Any}}' href='#DyVE.@cost-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@cost`** &mdash; *Macro*.



Set cost.

**Examples**

```julia
@cost model experimental1=2 experimental2=3
```

<a id='DyVE.@equalize-Tuple{Any, Vararg{Any}}' href='#DyVE.@equalize-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@equalize`** &mdash; *Macro*.



Identify (collapse) a set of species in a model.

**Examples**

```julia
@join acs acs1.A=acs2.A B=C
```

<a id='DyVE.@export_as_table' href='#DyVE.@export_as_table'>#</a>
**`DyVE.@export_as_table`** &mdash; *Macro*.



```julia
@export_as_table sol
```

Export a solution as a `DataFrame`.

**Examples**

```julia
@export_as_table sol
```

<a id='DyVE.@export_csv' href='#DyVE.@export_csv'>#</a>
**`DyVE.@export_csv`** &mdash; *Macro*.



```julia
@export_csv sol
@export_csv sol "sol.csv"
```

Export a solution to a file.

**Examples**

```julia
@export_csv sol "sol.csv"
```

<a id='DyVE.@export_model-Tuple{Any, Any}' href='#DyVE.@export_model-Tuple{Any, Any}'>#</a>
**`DyVE.@export_model`** &mdash; *Macro*.



Export model to a file.

**Examples**

```julia
@export_model acs "acs_data.toml"
```

<a id='DyVE.@export_solution' href='#DyVE.@export_solution'>#</a>
**`DyVE.@export_solution`** &mdash; *Macro*.



```julia
@export sol
@export sol "sol.jld2"
```

Export a solution to a file.

**Examples**

```julia
@export_solution sol "sol.jdl2"
```

<a id='DyVE.@fit-Tuple{Any, Any, Any, Vararg{Any}}' href='#DyVE.@fit-Tuple{Any, Any, Any, Vararg{Any}}'>#</a>
**`DyVE.@fit`** &mdash; *Macro*.



```julia
@fit acset data_points time_steps empiric_variables <free_var=[init_val]>... <free_prm=[init_val]>... opts...
```

Take an acset and fit initial values and parameters to empirical data.

The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model. The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

**Examples**

```julia
t = [1, 50, 100]
data = [80 30 20]
@fit acs data t vars=A B=20 A α # fit B, A, α; empirical data is for variable A
```

<a id='DyVE.@fit_and_plot-Tuple{Any, Any, Any, Vararg{Any}}' href='#DyVE.@fit_and_plot-Tuple{Any, Any, Any, Vararg{Any}}'>#</a>
**`DyVE.@fit_and_plot`** &mdash; *Macro*.



```julia
@fit acset data_points time_steps empiric_variables <free_var=[init_val]>... <free_prm=[init_val]>... opts...
```

Take an acset, fit initial values and parameters to empirical data, and plot the result.

The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model. The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

**Examples**

```julia
t = [1, 50, 100]
data = [80 30 20]
@fit acs data t vars=A B=20 A α # fit B, A, α; empirical data is for variable A
```

<a id='DyVE.@import_model' href='#DyVE.@import_model'>#</a>
**`DyVE.@import_model`** &mdash; *Macro*.



Import a model from a file.

**Examples**

```julia
@import_model "model.toml"
```

<a id='DyVE.@import_solution' href='#DyVE.@import_solution'>#</a>
**`DyVE.@import_solution`** &mdash; *Macro*.



```julia
@import "sol.jld2"
@import "sol.jld2" sol
```

Import a solution from a file.

**Examples**

```julia
@import_solution "sir_acs_sol/serialized/sol.jld2"
```

<a id='DyVE.@join-Tuple' href='#DyVE.@join-Tuple'>#</a>
**`DyVE.@join`** &mdash; *Macro*.



Performs join of models and identifies model variables, as specified.

Model variables / parameter values and metadata are propagated; the last model takes precedence.

**Examples**

```julia
@join acs1 acs2 @catchall(A)=acs2.Z @catchall(XY) @catchall(B)
```

<a id='DyVE.@jump-Tuple{Any, Any, Any}' href='#DyVE.@jump-Tuple{Any, Any, Any}'>#</a>
**`DyVE.@jump`** &mdash; *Macro*.



Add a jump process (with specified Poisson intensity per unit time step) to a model.

Followed with ACS*name Poisson*rate (conditiuon1 && condition2 && … && (action1; action2; …)).

**Examples**

```julia
@jump acs λ Z += rand(Poisson(1.))
```

<a id='DyVE.@mode-Tuple{Any, Any, Any}' href='#DyVE.@mode-Tuple{Any, Any, Any}'>#</a>
**`DyVE.@mode`** &mdash; *Macro*.



Set species modality.

**Supported modalities**

  * nonblock
  * conserved
  * rate

**Examples**

@mode acs (r"proj\w+", r"experimental\w+") conserved @mode acs (S, I) conserved @mode acs S conserved

<a id='DyVE.@name_transition-Tuple{Any, Vararg{Any}}' href='#DyVE.@name_transition-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@name_transition`** &mdash; *Macro*.



Set name of a transition in the model.

**Examples**

@name*transition acs 1="name" @name*transition acs name="transition*name" @name*transition acs "name"="transition_name"

<a id='DyVE.@optimize-Tuple{Any, Any, Vararg{Any}}' href='#DyVE.@optimize-Tuple{Any, Any, Vararg{Any}}'>#</a>
**`DyVE.@optimize`** &mdash; *Macro*.



```julia
@optimize acset objective <free_var=[init_val]>... <free_prm=[init_val]>... opts...
```

Take an acset and optimize given functional.

Objective is an expression which may reference the model's variables and parameters, i.e., `A+β`. The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model. The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

By default, the functional is minimized. Specify `objective=max` to perform maximization. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

**Examples**

```julia
@optimize acs abs(A-B) A B=20. α=2. lower_bounds=0 upper_bounds=100
@optimize acss abs(A-B) A B=20. α=2. upper_bounds=[200,300,400] maxeval=200 objective=min
```

<a id='DyVE.@periodic-Tuple{Any, Any, Any}' href='#DyVE.@periodic-Tuple{Any, Any, Any}'>#</a>
**`DyVE.@periodic`** &mdash; *Macro*.



Add a periodic callback to a model.

**Examples**

```julia
@periodic acs 1. X += 1
```

<a id='DyVE.@plot-Tuple{Any, Vararg{Any}}' href='#DyVE.@plot-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@plot`** &mdash; *Macro*.



Plot the solution (summary).

**Examples**

@plot sol plot*type=summary @plot sol plot*type=allocation # not supported for ensemble solutions! @plot sol plot*type=valuations # not supported for ensemble solutions! @plot sol plot*type=new_transitions # not supported for ensemble solutions!

<a id='DyVE.@prob_check_verbose-Tuple{Any}' href='#DyVE.@prob_check_verbose-Tuple{Any}'>#</a>
**`DyVE.@prob_check_verbose`** &mdash; *Macro*.



Check model parameters have been set.

<a id='DyVE.@prob_init-Tuple{Any, Vararg{Any}}' href='#DyVE.@prob_init-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@prob_init`** &mdash; *Macro*.



Set initial values of species in an acset.

Two different formats can be applied. Followed with ACS*name initial*values.

**Examples**

```julia
@prob_init acs X=1 Y=2 Z=h(α)
@prob_init acs [1., 2., 3.]
```

<a id='DyVE.@prob_meta-Tuple{Any, Vararg{Any}}' href='#DyVE.@prob_meta-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@prob_meta`** &mdash; *Macro*.



Set model metadata (i.e., solver arguments)

Followed with arguments, e.g. ACS_name tspan tstep.

**Examples**

```julia
@prob_meta acs tspan=(0, 100.) schedule=schedule_weighted!
@prob_meta sir_acs tspan=250 tstep=1
```

<a id='DyVE.@prob_params-Tuple{Any, Vararg{Any}}' href='#DyVE.@prob_params-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@prob_params`** &mdash; *Macro*.



Set parameter values in an acset.

Followed with ACS*name parameter*values.

**Examples**

```julia
@prob_params acs α=1. β=2.
```

<a id='DyVE.@prob_uncertainty-Tuple{Any, Vararg{Any}}' href='#DyVE.@prob_uncertainty-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@prob_uncertainty`** &mdash; *Macro*.



Set uncertainty in initial values of species in an acset (stderr).

Followed with ACS*name standard*error.

**Examples**

```julia
@prob_uncertainty acs X=.1 Y=.2
@prob_uncertainty acs [.1, .2,]
```

<a id='DyVE.@problematize-Tuple{Any, Vararg{Any}}' href='#DyVE.@problematize-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@problematize`** &mdash; *Macro*.



Convert a model to a `DiscreteProblem`. If passed a problem instance, return the instance.

**Examples**

@problematize acs tspan=1:100

<a id='DyVE.@push-Tuple{Any, Vararg{Any}}' href='#DyVE.@push-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@push`** &mdash; *Macro*.



Add reactions to an acset.

**Examples**

```julia
@push sir_acs β*S*I*tdecay(@time()) S+I --> 2I name=>SI2I
@push sir_acs begin 
    ν*I, I --> R, name=>I2R
    γ, R --> S, name=>R2S
end
```

<a id='DyVE.@register-Tuple{Any}' href='#DyVE.@register-Tuple{Any}'>#</a>
**`DyVE.@register`** &mdash; *Macro*.



Evaluate expression in DyVE scope.

**Examples**

```julia
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t/10^3)
```

<a id='DyVE.@reward-Tuple{Any, Vararg{Any}}' href='#DyVE.@reward-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@reward`** &mdash; *Macro*.



Set reward.

**Examples**

```julia
@reward model experimental1=2 experimental2=3
```

<a id='DyVE.@solve-Tuple{Any, Vararg{Any}}' href='#DyVE.@solve-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@solve`** &mdash; *Macro*.



Solve the problem. Solverargs passed at the calltime take precedence.

**Examples**

@solve prob @solve prob tspan=1:100 @solve prob tspan=100 trajectories=20

<a id='DyVE.@valuation-Tuple{Any, Vararg{Any}}' href='#DyVE.@valuation-Tuple{Any, Vararg{Any}}'>#</a>
**`DyVE.@valuation`** &mdash; *Macro*.



Set valuation.

**Examples**

```julia
@valuation model experimental1=2 experimental2=3
```