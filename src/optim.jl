function build_parametrized_solver(acs, init_vec, u0, params; trajectories=1)
    prob = DiscreteProblem(acs)
    vars = prob.p[:__state__][:, :specInitUncertainty]
    init_vec = deepcopy(init_vec)

    function (vec)
        vec = vec isa ComponentVector ? vec : (init_vec .= vec)
        data = []
        for _ in 1:trajectories
            prob.p[:__state__] = deepcopy(prob.p[:__state0__])
            for i in eachindex(prob.u0)
                rv = randn()*vars[i]
                prob.u0[i] = (sign(rv + prob.u0[i]) == sign(prob.u0[i])) ? rv + prob.u0[i] : prob.u0[i]
            end

            for (i, k) in enumerate(wkeys(u0)) prob.u0[k] = vec.species[i] end
            for k in wkeys(params) prob.p[k] = vec[k] end

            sync!(prob.p[:__state__], prob.u0, prob.p)
            push!(data, solve(prob))
        end

        data
    end
end

function build_parametrized_solver_(acs, init_vec, u0, params; trajectories=1)
    prob = DiscreteProblem(acs)
    vars = prob.p[:__state__][:, :specInitUncertainty]
    init_vec = deepcopy(init_vec)

    function (vec)
        vec = vec isa ComponentVector ? vec : (init_vec .= vec; init_vec)
        data = map(1:trajectories) do _
            prob.p[:__state__] = deepcopy(prob.p[:__state0__])
            for i in eachindex(prob.u0)
                rv = randn()*vars[i]
                prob.u0[i] = (sign(rv + prob.u0[i]) == sign(prob.u0[i])) ? rv + prob.u0[i] : prob.u0[i]
            end

            for (i, k) in enumerate(wkeys(u0)) prob.u0[k] = vec.species[i] end
            for k in wkeys(params) prob.p[k] = vec[k] end

            sync!(prob.p[:__state__], prob.u0, prob.p)
            
            solve(prob)
        end

        trajectories == 1 ? data[1] : EnsembleSolution(data, .0, true)
    end
end

## optimization part

BOUND_DEFAULT = 5000

function optim!(obj, init; nlopt_kwargs...)
    nlopt_kwargs = Dict(nlopt_kwargs)
    alg = pop!(nlopt_kwargs, :algorithm, :GN_DIRECT)

    opt = Opt(alg, length(init))

    # match to a ComponentVector
    foreach(o -> setproperty!(opt, o...), filter(x -> x[1] in propertynames(opt), nlopt_kwargs))
    get(nlopt_kwargs, :objective, min) == min ? (opt.min_objective = obj) : (opt.max_objective = obj)

    optimize(opt, deepcopy(init))
end

const n_steps = 100

# loss objective given an objective expression
function build_loss_objective(acs, init_vec, u0, params, obex; loss=identity, trajectories=1, min_t=-Inf, max_t=Inf, final_only=false)
    ob = eval(get_wrap_fun(acs)(obex))
    obj_ = build_parametrized_solver(acs, init_vec, u0, params; trajectories)

    function (vec, _)
        ls = []
        for sol in obj_(vec) 
            t_points = 
                if final_only; [last(sol.t)]
                else 
                    min_t = max(min_t, sol.prob.tspan[1])
                    max_t = min(max_t, sol.prob.tspan[2])

                    range(min_t, max_t; length=n_steps)
                end

            push!(ls, mean(t -> loss(ob(as_state(sol(t), t, sol.prob.p[:__state__]))), t_points)) 
        end
    
        mean(ls)
    end
end

# loss objective given empirical data
function build_loss_objective_datapoints(acs, init_vec, u0, params, t, data, vars; loss=abs2, trajectories=1)
    obj_ = build_parametrized_solver(acs, init_vec, u0, params; trajectories)

    function (vec, _)
        ls = []
        for sol in obj_(vec)
            push!(ls, mean(t -> sum(i -> loss(sol(t[2])[i[2]] - data[i[1], t[1]]), enumerate(vars)), enumerate(t))) 
        end
        
        mean(ls)
    end
end

# set initial model parameter values in an optimization problem
function prep_params!(params, prob)
    for (k, v) in params; (v === NaN) && wset!(params, k, get(prob.p, k, NaN)) end
    any(p -> (p[2] === NaN) && @warn("Uninitialized prm: $p"), params)

    params
end

# set initial model variable values in an optimization problem
function prep_u0!(u0, prob)
    for (k, v) in u0; (v === NaN) && wset!(u0, k, get(prob.u0, k, NaN)) end
    any(u -> (u[2] === NaN) && @warn("Uninitialized prm: $(u[1])"), u0)

    u0
end

"Extract symbolic variables referenced in `acs`, `args`."
function get_free_vars(acs, args)
    u0_syms = collect(acs[:, :specName]); p_syms = collect(acs[:, :prmName])
    u0 = []; p = []

    for arg in args
        if arg isa Symbol; (k, v) = (arg, NaN)
        elseif isexpr(arg, :(=)); (k, v) = (arg.args[1], arg.args[2])
        else continue end
        
        if ((k in u0_syms || k isa Number) && !in(k, wkeys(u0))) push!(u0, k => v) 
        elseif (k in p_syms && !in(k, wkeys(p))) push!(p, k => v) end
    end

    u0_ = []; for (k, v) in u0;
        if k isa Number push!(u0_, Int(k) => v)
        else 
            for i in 1:length(subpart(acs, :specName))
                (acs[i, :specName] == k) && (push!(u0_, i => v); break)
            end
        end
    end

    u0_, p
end

"Resolve symbolic / positional model variable names to positional."
function get_vars(acs, args)
    (args == :()) && return args
    args_ = []

    for arg in (MacroTools.isexpr(args, :vect, :tuple) ? args.args : [args])
        arg = recursively_expand_dots(arg)
        if arg isa Number push!(args_, Int(arg))
        else 
            for i in 1:length(subpart(acs, :specName))
                isassigned(subpart(acs, :specName), i) && (acs[i, :specName] == arg) &&
                    (push!(args_, i); break)
            end
        end
    end

    args_
end