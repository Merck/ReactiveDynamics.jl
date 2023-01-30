# reaction network DSL: UPDATE part; add species, name, add modalities, set model variables, set solver arguments 

export @push, @name_transition, @mode, @add_species
export @periodic, @jump
export @prob_init, @prob_uncertainty, @prob_params, @prob_meta
export @prob_role, @list_by_role, @list_roles
export @prob_check_verbose
export @aka
export @register

using DataFrames
using MacroTools: striplines

function push_to_acs!(acsex, exs...)
    if isexpr(exs[1], :block)
        ex = striplines(exs[1])
    else
        args = Any[]
        map(el -> isexpr(el, :(=)) ? push!(args, :($(el.args[1]) => $(el.args[2]))) :
                  push!(args, el), exs)
        ex = Expr(:tuple, args...)
    end

    quote
        ex = blockize($(QuoteNode(ex)))
        merge_acs!($(esc(acsex)), get_data(ex)...)
    end
end

"""
Add reactions to an acset.

# Examples
```julia
@push sir_acs β*S*I*tdecay(@time()) S+I --> 2I name=>SI2I
@push sir_acs begin 
    ν*I, I --> R, name=>I2R
    γ, R --> S, name=>R2S
end
```
"""
macro push(acsex, exs...)
    push_to_acs!(acsex, exs...)
end

"""
Set name of a transition in the model.

# Examples
```julia
@name_transition acs 1="name"
@name_transition acs name="transition_name"
@name_transition acs "name"="transition_name"
```
"""
macro name_transition(acsex, exs...)
    call = :(begin end)
    for ex in exs
        call_ = if ex.args[1] isa Number
            :($(esc(acsex))[$(ex.args[1]), :transName] = $(QuoteNode(ex.args[2])))
        else
            quote
                acs = $(esc(acsex))
                ixs = findall(i -> string(acs[i, :transName]) == $(string(ex.args[1])),
                              1:nparts(acs, :T))
                foreach(i -> acs[i, :transName] = $(string(ex.args[2])), ixs)
            end
        end

        push!(call.args, call_)
    end

    call
end

function incident_pattern(pattern, attr)
    ix = []
    for i in 1:length(attr)
        !isnothing(attr[i]) &&
            (m = match(pattern, string(attr[i]));
             !isnothing(m) &&
                 (string(attr[i]) == m.match)) && push!(ix, i)
    end

    ix
end

function mode!(acs, dict)
    for (spex, mods) in dict
        i = if spex isa Regex
            incident_pattern(spex, acs[:, :specName])
        else
            incident(acs, Symbol(spex), :specName)
        end

        for ix in i
            isnothing(acs[ix, specModality]) &&
                (acs[ix, specModality] = Set{Symbol}())
            union!(acs[ix, :specModality], mods)
        end
    end
end

"""
Set species modality.

# Supported modalities
 - nonblock
 - conserved
 - rate

# Examples
```julia
@mode acs (r"proj\\w+", r"experimental\\w+") conserved
@mode acs (S, I) conserved
@mode acs S conserved
```
"""
macro mode(acsex, spexs, mexs)
    mods = !isa(mexs, Expr) ? [mexs] : collect(mexs.args)
    spexs = isexpr(spexs, :tuple) ? spexs.args : [spexs]
    exs = map(ex -> striplines(ex), spexs)

    quote
        dictcall = Dict()
        exs_ = []
        foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
        exs__ = []
        foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
        foreach(ex -> push!(dictcall, get_pattern(recursively_expand_dots(ex)) => $mods),
                exs__)

        mode!($(esc(acsex)), dictcall)
    end
end

function set_valuation!(acs, dict, valuation_type)
    for (spex, val) in dict
        i = if spex isa Regex
            incident_pattern(spex, subpart(acs, :specName))
        else
            incident(acs, Symbol(spex), :specName)
        end

        foreach(ix -> acs[ix, Symbol(:spec, Symbol(uppercasefirst(string(valuation_type))))] = eval(val),
                i)
    end
end

export @cost, @reward, @valuation
for valuation_type in (:cost, :reward, :valuation)
    eval(quote
             export $(Symbol(Symbol("@"), valuation_type))
             @doc """
             Set $($(string(valuation_type))).

             # Examples
             ```julia
             @$($(string(valuation_type))) model experimental1=2 experimental2=3
             ```
             """
             macro $valuation_type(acsex, exs...)
                 dictcall = Dict()
                 exs_ = []
                 foreach(s -> push!(exs_, striplines(blockize(s))), exs)
                 exs__ = []
                 foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
                 foreach(ex -> push!(dictcall,
                                     get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2]),
                         exs__)
                 :(set_valuation!($(esc(acsex)), $dictcall,
                                  $(QuoteNode($(QuoteNode(valuation_type))))))
             end
         end)
end

"""
Add new species to a model.

# Examples
```julia
@add_species acs S I R
```
"""
macro add_species(acsex, exs...)
    call = :(begin end)
    spexs_ = []
    foreach(s -> push!(spexs_, s), exs)

    for ex in recursively_expand_dots.(spexs_)
        push!(call.args, :(add_part!($(esc(acsex)), :S; specName = $(QuoteNode(ex)))))
    end

    push!(call.args, :(assign_defaults!($(esc(acsex)))))
    call
end

get_pattern(ex) = ex isa Expr && (macroname(ex) == :r_str) ? eval(ex) : ex

"""
Set initial values of species in an acset.

# Examples
```julia
@prob_init acs X=1 Y=2 Z=h(α)
@prob_init acs [1., 2., 3.]
```
"""
macro prob_init(acsex, exs...)
    exs = map(ex -> striplines(ex), exs)

    return if length(exs) == 1 && (isexpr(exs[1], :vect) || (exs[1] isa Symbol))
        :(init!($(esc(acsex)), $(esc(exs[1]))))
    else
        quote
            dictcall = Dict()
            exs_ = []
            foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
            exs__ = []
            foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
            foreach(ex -> push!(dictcall,
                                get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2]),
                    exs__)

            init!($(esc(acsex)), dictcall)
        end
    end
end

# deprecate
macro prob_init_from_vec(acsex, vecex)
    :(init!($(esc(acsex)), $(esc(vecex))))
end

function init!(acs, inits)
    inits isa AbstractVector && length(inits) == nparts(acs, :S) &&
        (subpart(acs, :specInitVal) .= inits; return)
    inits isa AbstractDict && for (k, init_val) in inits
        k isa Number ? acs[k, :specInitVal] = init_val :
        begin
            i = k isa Regex ? incident_pattern(k, subpart(acs, :specName)) :
                incident(acs, k, :specName)
            foreach(ix -> (acs[ix, :specInitVal] = init_val), i)
        end
    end

    acs
end

"""
Set uncertainty in initial values of species in an acset (stderr).

# Examples
```julia
@prob_uncertainty acs X=.1 Y=.2
@prob_uncertainty acs [.1, .2,]
```
"""
macro prob_uncertainty(acsex, exs...)
    exs = map(ex -> striplines(ex), exs)

    return if length(exs) == 1 && (isexpr(exs[1], :vect) || (exs[1] isa Symbol))
        :(uncinit!($(esc(acsex)), $(esc(exs[1]))))
    else
        quote
            dictcall = Dict()
            exs_ = []
            foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
            exs__ = []
            foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
            foreach(ex -> push!(dictcall,
                                get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2]),
                    exs__)

            uncinit!($(esc(acsex)), dictcall)
        end
    end
end

function uncinit!(acs, inits)
    inits isa AbstractVector && length(inits) == nparts(acs, :S) &&
        (subpart(acs, :specInitUncertainty) .= inits; return)
    inits isa AbstractDict && for (k, init_val) in inits
        k isa Number ? acs[k, :specInitUncertainty] = init_val :
        begin
            i = k isa Regex ? incident_pattern(k, subpart(acs, :specName)) :
                incident(acs, k, :specName)
            foreach(ix -> (acs[ix, :specInitUncertainty] = init_val), i)
        end
    end

    acs
end

function set_params!(acs, params)
    params isa AbstractDict && for (k, init_val) in params
        k = get_pattern(k)
        if k isa Regex
            i = incident_pattern(k, subpart(acs, :prmName))
        else
            i = incident(acs, k, :prmName)
            isempty(i) && (i = add_part!(acs, :P, prmName = k))
        end

        foreach(ix -> acs[ix, :prmVal] = eval(init_val), i)
    end
end

"""
Set parameter values in an acset.

# Examples
```julia
@prob_params acs α=1. β=2.
```
"""
macro prob_params(acsex, exs...)
    exs = map(ex -> striplines(ex), exs)

    quote
        dictcall = Dict()
        exs_ = []
        foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
        exs__ = []
        foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
        foreach(ex -> push!(dictcall,
                            get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2]),
                exs__)

        set_params!($(esc(acsex)), dictcall)
    end
end

function meta!(acs, metas)
    for (k, metaval) in metas
        i = incident(acs, k, :metaKeyword)
        isempty(i) && (i = add_part!(acs, :M, metaKeyword = k))
        set_subpart!(acs, first(i), :metaVal, metaval)
    end
end

"""
Set model metadata (e.g. solver arguments)

# Examples
```julia
@prob_meta acs tspan=(0, 100.) schedule=schedule_weighted!
@prob_meta sir_acs tspan=250 tstep=1
```
"""
macro prob_meta(acsex, exs...)
    dictcall = :(Dict([]))
    foreach(ex -> push!(dictcall.args[2].args, (ex.args[1] => eval(ex.args[2]))), exs)

    :(meta!($(esc(acsex)), $dictcall))
end

"""
Alias object name in an acs.

# Default names
| name | short name |
| :--- | :--- |
| species | S |
| transition | T |
| action | A |
| event | E |
| param | P |
| meta | M |

# Examples
```julia
@aka acs species=resource transition=reaction
```
"""
macro aka(acsex, exs...)
    dictcall = :(Dict([]))
    foreach(ex -> push!(dictcall.args[2].args,
                        Symbol("alias_", findfirst(==(ex.args[1]), alias_default)) => ex.args[2]),
            exs)
    :(meta!($(esc(acsex)), $dictcall))
end

alias_default = Dict(:S => :species, :T => :transition, :A => :action, :E => :event,
                     :P => :param, :M => :meta)

function get_alias(acs, ob)
    (i = incident(acs, Symbol(:alias_, ob), :metaKeyword);
     !isempty(i) ?
     acs[first(i), :metaVal] :
     alias_default[ob])
end

"Check model parameters have been set."
macro prob_check_verbose(acsex) # msg as return value
    :(missing_params = check_params($(esc(acsex)));
      isempty(missing_params) ? "Params OK." :
      "Missing params: $missing_params")
end

"""
Add a periodic callback to a model.

# Examples
```julia
@periodic acs 1. X += 1
```
"""
macro periodic(acsex, pex, acex)
    push_to_acs!(acsex, Expr(:&&, :(@periodic($pex)), acex))
end

"""
Add a jump process (with specified Poisson intensity per unit time step) to a model.

# Examples
```julia
@jump acs λ Z += rand(Poisson(1.))
```
"""
macro jump(acsex, inex, acex)
    push_to_acs!(acsex,
                 Expr(:&&,
                      Expr(:call, :rand, :(Poisson(max(@solverarg(:tstep) * $inex, 0)))),
                      acex))
end

"""
Evaluate expression in ReactiveDynamics scope.

# Examples
```julia
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t/10^3)
```
"""
macro register(ex)
    :(@eval ReactiveDynamics $ex)
end
