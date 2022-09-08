export @problematize, @solve, @plot
export @optimize, @fit, @fit_and_plot, @build_solver

using DifferentialEquations: DiscreteProblem, EnsembleProblem, FunctionMap, EnsembleSolution
import MacroTools
import Plots

"""
Convert a model to a `DiscreteProblem`. If passed a problem instance, return the instance.

# Examples
@problematize acs tspan=1:100
"""
macro problematize(acsex, args...)
    args, kwargs = args_kwargs(args)
    quote
        $(esc(acsex)) isa DiscreteProblem ? $(esc(acsex)) : 
            DiscreteProblem($(esc(acsex)), $(args...); $(kwargs...))
    end
end

"""
Solve the problem. Solverargs passed at the calltime take precedence.

# Examples
@solve prob
@solve prob tspan=1:100
@solve prob tspan=100 trajectories=20
"""
macro solve(probex, args...)
    args, kwargs = args_kwargs(args)
    mode = find_kwargex_delete!(kwargs, :mode, nothing)
    !isnothing(findfirst(el -> el.args[1] == :trajectories, kwargs)) && (mode = :ensemble)

    quote
        prob = $(esc(probex)) isa DiscreteProblem ? $(esc(probex)) : DiscreteProblem($(esc(probex)), $(args...); $(kwargs...))
        if $(preserve_sym(mode)) == :ensemble
            solve(EnsembleProblem(prob; prob_func=get_prob_func(prob)), FunctionMap(), $(kwargs...))
        else solve(prob) end
    end
end

# auxiliary plotting functions
function plot_summary(s, labels, ixs; kwargs...) 
    isempty(ixs) && return @warn "Set of species to plot must be non-empty!"
    s = EnsembleSummary(s)
    f_ix = first(ixs); p = Plots.plot(s.t, s.u[f_ix,:], ribbon=(-s.qlow[f_ix,:]+s.u[f_ix,:], s.qhigh[f_ix,:]-s.u[f_ix,:]), label=labels[f_ix], fillalpha=.2, w=2.; kwargs...)
    foreach(i -> Plots.plot!(p, s.t, s.u[i,:], ribbon=(-s.qlow[i,:]+s.u[i,:], s.qhigh[i,:]-s.u[i,:]), label=labels[i], fillalpha=.2, w=2.; kwargs...),
        ixs[2:end]
    )

    p
end

function plot_ensemble_sol(sol, label, ixs; kwargs...)
    !(sol isa EnsembleSolution) ? Plots.plot(sol, idxs=ixs, label=reshape(label[ixs], 1, :); kwargs...) : Plots.plot(sol, idxs=ixs, alpha=.7; kwargs...)
end

first_sol(sol) = sol isa EnsembleSolution ? first(sol) : sol

function match_names(selector, names)
    if isnothing(selector) 1:length(names)
    else 
        selector = selector isa Union{AbstractString, Regex, Symbol} ? [selector] : selector
        ixs = Int[]; for s in selector
            s isa Symbol && (s = string(s))
            append!(ixs, findall(name -> s isa Regex ? occursin(s, name) : (string(s) == string(name)), names))
        end
        unique!(ixs)
    end
end

"""
Plot the solution (summary).

# Examples
@plot sol plot_type=summary
@plot sol plot_type=allocation # not supported for ensemble solutions!
@plot sol plot_type=valuations # not supported for ensemble solutions!
@plot sol plot_type=new_transitions # not supported for ensemble solutions!
"""
macro plot(solex, args...)
    _, kwargs = args_kwargs(args)
    plot = find_kwargex_delete!(kwargs, :plot_type, nothing)
    selector = find_kwargex_delete!(kwargs, :show, nothing)

    quote
        plot_type = $(preserve_sym(plot)); sol = $(esc(solex))
        if plot_type ∈ [nothing, :ensemble]
            names = string.(first_sol(sol).prob.p[:__state__][:, :specName])
            plot_ensemble_sol(sol, names, match_names($selector, names); $(kwargs...))
        elseif plot_type == :summary
            names = string.(first_sol(sol).prob.p[:__state__][:, :specName])[:]
            plot_summary(sol, names, match_names($selector, names); $(kwargs...))
        else
            names = string.(sol.prob.p[:__state__][:, :specName])
            plot_from_log(sol.prob.p[:__state__], plot_type, match_names($selector, names); $(kwargs...)) end
    end
end

function get_transitions(log)
    transitions = Symbol[]; for r in log
        r[1] ∈ [:new_transitions, :saturation] && union!(transitions, map(r -> r[1], r[3:end]))
    end

    transitions
end

function complete_log_transitions(log, hashes, record_type)
    pos = Dict(hashes[ix] => ix for ix in eachindex(hashes))
    steps = unique(map(r -> r[2], log)); vals = zeros(length(steps), length(hashes))
    ix = 1; for r in log
        if r[1] == record_type
            for (hash, val) in r[3:end]; vals[ix, pos[hash]] = val end
            ix += 1
        end
    end

    vals, steps
end

function complete_log_species(log, n_species, record_type)
    steps = unique(map(r -> r[2], log)); vals = zeros(length(steps), n_species)
    ix = 1; for r in log
        if r[1] == record_type
            vals[ix, :] .= r[3]; ix += 1
        end
    end

    vals, steps
end

function complete_log_valuations(log)
    steps = unique(map(r -> r[2], log)); vals = zeros(length(steps), 4)
    ix = 0; time_last = -Inf; for r in log
        r[2] > time_last && (time_last = r[2]; ix += 1)
        if r[1] == :valuation; vals[ix, 1] = r[3]
        elseif r[1] == :valuation_cost; vals[ix, 2] = r[3]
        elseif r[1] == :valuation_reward; vals[ix, 3] = r[3]  end
    end

    vals, steps
end

function get_names(state, hashes)
    names = String[]
    for hash in hashes
        ix = findfirst(==(hash), state[:, :transHash])
        push!(names, assigned(state.transition_recipes[:transName], ix) ? string(state[ix, :transName]) : "transition $ix")
    end
 
    names
end

function plot_from_log(state, record_type, ixs; kwargs...)
    if record_type ∈ [:new_transitions, :terminated_transitions, :saturation]
        hashes = get_transitions(state.log)
        label = get_names(state, hashes)
        vals, steps = complete_log_transitions(state.log, hashes, record_type)
    elseif record_type ∈ [:allocation]
        label = string.(state[:, :specName])
        vals, steps = complete_log_species(state.log, length(label), record_type)
        vals = vals[:, ixs]; label = label[ixs]
    elseif record_type ∈ [:valuations]
        label = ["portfolio valuation", "costs", "rewards", "total balance"]
        vals, steps = complete_log_valuations(state.log)
    end

    Plots.plot(steps, vals, title=string(record_type), label=reshape(label, 1, :); kwargs...)
end

"""
    @optimize acset objective <free_var=[init_val]>... <free_prm=[init_val]>... opts...

Take an acset and optimize given functional.

Objective is an expression which may reference the model's variables and parameters, i.e., `A+β`.
The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model.
The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

By default, the functional is minimized. Specify `objective=max` to perform maximization. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

# Examples
```julia
@optimize acs abs(A-B) A B=20. α=2. lower_bounds=0 upper_bounds=100
@optimize acss abs(A-B) A B=20. α=2. upper_bounds=[200,300,400] maxeval=200 objective=min
```
"""
macro optimize(acsex, obex, args...)
    args_all = args; args, kwargs = args_kwargs(args)
    min_t = find_kwargex_delete!(kwargs, :min_t, -Inf)
    max_t = find_kwargex_delete!(kwargs, :max_t, Inf)
    final_only = find_kwargex_delete!(kwargs, :final_only, false)
    okwargs = filter(ex -> ex.args[1] in [:loss, :trajectories], kwargs)
    
    quote
        u0, p = get_free_vars($(esc(acsex)), $(QuoteNode(args_all)))
        prob_ = DiscreteProblem($(esc(acsex)))
        prep_u0!(u0, prob_); prep_params!(p, prob_)

        init_p = [k => v for (k, v) in p]
        init_vec = if length(u0) > 0
                ComponentVector{Float64}(; species=collect(wvalues(u0)), init_p...)
            else
                ComponentVector{Float64}(; init_p...)
            end

        o = build_loss_objective($(esc(acsex)), init_vec, u0, p, $(QuoteNode(obex)); min_t=$min_t, max_t=$max_t, final_only=$final_only, $(okwargs...))

        optim!(o, init_vec; $(kwargs...))
    end
end

"""
    @fit acset data_points time_steps empiric_variables <free_var=[init_val]>... <free_prm=[init_val]>... opts...

Take an acset and fit initial values and parameters to empirical data.

The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model.
The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

# Examples
```julia
t = [1, 50, 100]
data = [80 30 20]
@fit acs data t vars=A B=20 A α # fit B, A, α; empirical data is for variable A
```
"""
macro fit(acsex, data, t, args...)
    args_all = args; args, kwargs = args_kwargs(args)
    okwargs = filter(ex -> ex.args[1] in [:loss, :trajectories], kwargs)
    vars = (ix = findfirst(ex -> ex.args[1] == :vars, kwargs); !isnothing(ix) ? (v=kwargs[ix].args[2]; deleteat!(kwargs, ix); v) : :())

    quote
        u0, p = get_free_vars($(esc(acsex)), $(QuoteNode(args_all)))
        vars = get_vars($(esc(acsex)), $(QuoteNode(vars)))
        prob_ = DiscreteProblem($(esc(acsex)))
        prep_u0!(u0, prob_); prep_params!(p, prob_)
        
        init_p = [k => v for (k, v) in p]
        init_vec = if length(u0) > 0
            ComponentVector{Float64}(; species=collect(wvalues(u0)), init_p...)
        else
            ComponentVector{Float64}(; init_p...)
        end

        o = build_loss_objective_datapoints($(esc(acsex)), init_vec, u0, p, $(esc(t)), $(esc(data)), vars; $(okwargs...))

        optim!(o, init_vec; $(kwargs...))
    end
end

"""
    @fit acset data_points time_steps empiric_variables <free_var=[init_val]>... <free_prm=[init_val]>... opts...

Take an acset, fit initial values and parameters to empirical data, and plot the result.

The values to optimized are listed using their symbolic names; unless specified, the initial value is inferred from the model.
The vector of free variables passed to the `NLopt` solver has the form `[free_vars; free_params]`; order of vars and params, respectively, is preserved. 

Propagates `NLopt` solver arguments; see [NLopt documentation](https://github.com/JuliaOpt/NLopt.jl).

# Examples
```julia
t = [1, 50, 100]
data = [80 30 20]
@fit acs data t vars=A B=20 A α # fit B, A, α; empirical data is for variable A
```
"""
macro fit_and_plot(acsex, data, t, args...)
    args_all = args
    trajectories = get_kwarg(args, :trajectories, 1)
    args, kwargs = args_kwargs(args)
    okwargs = filter(ex -> ex.args[1] in [:loss, :trajectories], kwargs)
    vars = (ix = findfirst(ex -> ex.args[1] == :vars, kwargs); !isnothing(ix) ? (v=kwargs[ix].args[2]; deleteat!(kwargs, ix); v) : :())

    quote
        u0, p = get_free_vars($(esc(acsex)), $(QuoteNode(args_all)))
        vars = get_vars($(esc(acsex)), $(QuoteNode(vars)))
        prob_ = DiscreteProblem($(esc(acsex)); suppress_warning=true)
        prep_u0!(u0, prob_); prep_params!(p, prob_)
        
        init_p = [k => v for (k, v) in p]
        init_vec = if length(u0) > 0
            ComponentVector{Float64}(; species=collect(wvalues(u0)), init_p...)
        else
            ComponentVector{Float64}(; init_p...)
        end

        o = build_loss_objective_datapoints($(esc(acsex)), init_vec, u0, p, $(esc(t)), $(esc(data)), vars; $(okwargs...))

        r = optim!(o, init_vec; $(kwargs...))
        if r[3] != :FORCED_STOP
            s_ = build_parametrized_solver($(esc(acsex)), init_vec, u0, p; trajectories=$trajectories)
            sol = first(s_(init_vec)); sol_ = first(s_(r[2]))

            p = Plots.plot(sol; idxs=vars, label="(initial) " .* reshape(String.($(esc(acsex))[:, :specName])[vars], 1, :))
            Plots.plot!(p, $(esc(t)), transpose($(esc(data))), label="(empirical) " .* reshape(String.($(esc(acsex))[:, :specName])[vars], 1, :))
            Plots.plot!(p, sol_; idxs=vars, label="(fitted) " .* reshape(String.($(esc(acsex))[:, :specName])[vars], 1, :))       
            p
        else :FORCED_STOP end
    end
end

"""
    @build_solver acset <free_var=[init_val]>... <free_prm=[init_val]>... opts...

Take an acset and export a solution as a function of free vars and free parameters.

# Examples
```julia
solver = @build_solver acs S α β # function of variable S and parameters α, β
solver([S, α, β])
```
"""
macro build_solver(acsex, args...)
    args_all = args#; args, kwargs = args_kwargs(args)
    trajectories = get_kwarg(args, :trajectories, 1)

    quote
        u0, p = get_free_vars($(esc(acsex)), $(QuoteNode(args_all)))
        prob_ = DiscreteProblem($(esc(acsex)))
        prep_u0!(u0, prob_); prep_params!(p, prob_)

        init_p = [k => v for (k, v) in p]
        init_vec = if length(u0) > 0
            ComponentVector{Float64}(; species=collect(wvalues(u0)), init_p...)
        else
            ComponentVector{Float64}(; init_p...)
        end

        build_parametrized_solver_($(esc(acsex)), init_vec, u0, p; trajectories=$trajectories)
    end
end