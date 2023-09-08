export @agentize, @solve, @plot
export @optimize, @fit, @fit_and_plot, @build_solver

import MacroTools
import Plots

"""
Convert a model to a `ReactiveNetwork`. If passed a problem instance, return the instance.

# Examples

```julia
@agentize acs tspan = 1:100
```
"""
macro agentize(acsex, args...)
    args, kwargs = args_kwargs(args)
    quote
        if $(esc(acsex)) isa ReactiveNetwork
            $(esc(acsex))
        else
            ReactiveNetwork($(esc(acsex)), $(args...); $(kwargs...))
        end
    end
end

"""
Solve the problem. Solverargs passed at the calltime take precedence.

# Examples

```julia
@solve prob
@solve prob tspan = 1:100
@solve prob tspan = 100
```
"""
macro solve(probex, args...)
    args, kwargs = args_kwargs(args)
    mode = find_kwargex_delete!(kwargs, :mode, nothing)
    !isnothing(findfirst(el -> el.args[1] == :trajectories, kwargs)) && (mode = :ensemble)

    quote
        prob = if $(esc(probex)) isa ReactiveNetwork
            $(esc(probex))
        else
            ReactiveNetwork($(esc(probex)), $(args...); $(kwargs...))
        end
        
        simulate(prob)
    end
end

# auxiliary plotting functions
function plot_summary(s, labels, ixs; kwargs...)
    isempty(ixs) && return @warn "Set of species to plot must be non-empty!"
    s = EnsembleSummary(s)
    f_ix = first(ixs)
    p = Plots.plot(
        s.t,
        s.u[f_ix, :];
        ribbon = (-s.qlow[f_ix, :] + s.u[f_ix, :], s.qhigh[f_ix, :] - s.u[f_ix, :]),
        label = labels[f_ix],
        fillalpha = 0.2,
        w = 2.0,
        kwargs...,
    )
    foreach(
        i -> Plots.plot!(
            p,
            s.t,
            s.u[i, :];
            ribbon = (-s.qlow[i, :] + s.u[i, :], s.qhigh[i, :] - s.u[i, :]),
            label = labels[i],
            fillalpha = 0.2,
            w = 2.0,
            kwargs...,
        ),
        ixs[2:end],
    )

    return p
end

function plot_ensemble_sol(sol, label, ixs; kwargs...)
    return if !(sol isa EnsembleSolution)
        Plots.plot(sol; idxs = ixs, label = reshape(label[ixs], 1, :), kwargs...)
    else
        Plots.plot(sol; idxs = ixs, alpha = 0.7, kwargs...)
    end
end

first_sol(sol) = sol isa EnsembleSolution ? first(sol) : sol

function match_names(selector, names)
    if isnothing(selector)
        1:length(names)
    else
        selector = selector isa Union{AbstractString,Regex,Symbol} ? [selector] : selector
        ixs = Int[]
        for s in selector
            s isa Symbol && (s = string(s))
            append!(
                ixs,
                findall(
                    name -> s isa Regex ? occursin(s, name) : (string(s) == string(name)),
                    names,
                ),
            )
        end
        unique!(ixs)
    end
end

"""
Plot the solution (summary).

# Examples

```julia
@plot sol plot_type = summary
@plot sol plot_type = allocation # not supported for ensemble solutions!
@plot sol plot_type = valuations # not supported for ensemble solutions!
@plot sol plot_type = new_transitions # not supported for ensemble solutions!
```
"""
macro plot(solex, args...)
    _, kwargs = args_kwargs(args)
    plot = find_kwargex_delete!(kwargs, :plot_type, nothing)
    selector = find_kwargex_delete!(kwargs, :show, nothing)

    quote
        plot_type = $(preserve_sym(plot))
        sol = $(esc(solex))
        if plot_type ∈ [nothing, :ensemble]
            names = string.(first_sol(sol).prob.p[:__state__][:, :specName])
            plot_ensemble_sol(sol, names, match_names($selector, names); $(kwargs...))
        elseif plot_type == :summary
            names = string.(first_sol(sol).prob.p[:__state__][:, :specName])[:]
            plot_summary(sol, names, match_names($selector, names); $(kwargs...))
        else
            names = string.(sol.prob.p[:__state__][:, :specName])
            plot_from_log(
                sol.prob.p[:__state__],
                plot_type,
                match_names($selector, names);
                $(kwargs...),
            )
        end
    end
end

function get_transitions(log)
    transitions = Symbol[]
    for r in log
        r[1] ∈ [:new_transitions, :saturation] &&
            union!(transitions, map(r -> r[1], r[3:end]))
    end

    return transitions
end

function complete_log_transitions(log, hashes, record_type)
    pos = Dict(hashes[ix] => ix for ix in eachindex(hashes))
    steps = unique(map(r -> r[2], log))
    vals = zeros(length(steps), length(hashes))
    ix = 1
    for r in log
        if r[1] == record_type
            for (hash, val) in r[3:end]
                vals[ix, pos[hash]] = val
            end
            ix += 1
        end
    end

    return vals, steps
end

function complete_log_species(log, n_species, record_type)
    steps = unique(map(r -> r[2], log))
    vals = zeros(length(steps), n_species)
    ix = 1
    for r in log
        if r[1] == record_type
            vals[ix, :] .= r[3]
            ix += 1
        end
    end

    return vals, steps
end

function complete_log_valuations(log)
    steps = unique(map(r -> r[2], log))
    vals = zeros(length(steps), 4)
    ix = 0
    time_last = -Inf
    for r in log
        r[2] > time_last && (time_last = r[2]; ix += 1)
        if r[1] == :valuation
            vals[ix, 1] = r[3]
        elseif r[1] == :valuation_cost
            vals[ix, 2] = r[3]
        elseif r[1] == :valuation_reward
            vals[ix, 3] = r[3]
        end
    end

    return vals, steps
end

function get_names(state, hashes)
    names = String[]
    for hash in hashes
        ix = findfirst(==(hash), state[:, :transHash])
        push!(names, if assigned(state.transition_recipes[:transName], ix)
            string(state[ix, :transName])
        else
            "transition $ix"
        end)
    end

    return names
end

function plot_from_log(state, record_type, ixs; kwargs...)
    if record_type ∈ [:new_transitions, :terminated_transitions, :saturation]
        hashes = get_transitions(state.log)
        label = get_names(state, hashes)
        vals, steps = complete_log_transitions(state.log, hashes, record_type)
    elseif record_type ∈ [:allocation]
        label = string.(state[:, :specName])
        vals, steps = complete_log_species(state.log, length(label), record_type)
        vals = vals[:, ixs]
        label = label[ixs]
    elseif record_type ∈ [:valuations]
        label = ["portfolio valuation", "costs", "rewards", "total balance"]
        vals, steps = complete_log_valuations(state.log)
    end

    return Plots.plot(
        steps,
        vals;
        title = string(record_type),
        label = reshape(label, 1, :),
        kwargs...,
    )
end