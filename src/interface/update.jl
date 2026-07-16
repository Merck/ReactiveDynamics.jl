# reaction network DSL: UPDATE part; add species, name, add modalities, set model variables, set solver arguments 

export @push, @name_transition, @mode, @add_species
export @periodic, @jump
export @prob_init, @prob_uncertainty, @prob_params, @prob_meta
export @aka
export @register

# NOTE (WS-4 housekeeping): three dangling exports were DELETED from here — `@prob_role`,
# `@list_by_role`, `@list_roles` (a legacy roles/actors ontology that was never implemented: no
# macro definitions, no `specRole` schema attribute), and `@prob_check_verbose` (see below). The
# role concept was dropped; ADR 0009's `PortRole` is an unrelated per-Species `role` field authored
# inside `@reaction_network`, not a `@prob_role`-style config macro, so nothing is repurposed.

using DataFrames
using MacroTools: striplines

function push_to_acs!(acsex, exs...)
    if isexpr(exs[1], :block)
        ex = striplines(exs[1])
    else
        args = Any[]
        map(el -> if isexpr(el, :(=))
            push!(args, :($(el.args[1]) => $(el.args[2])))
        else
            push!(args, el)
        end, exs)
        ex = Expr(:tuple, args...)
    end

    return quote
        ex = blockize($(QuoteNode(ex)))
        merge_network!($(esc(acsex)), get_data(ex)...)
    end
end

"""
Add reactions to an network.

# Examples

```julia
@push sir_acs β * S * I * tdecay(@time()) S + I --> 2I name => SI2I
@push sir_acs begin
    ν * I, I --> R, name => I2R
    γ, R --> S, name => R2S
end
```
"""
macro push(acsex, exs...)
    return push_to_acs!(acsex, exs...)
end

"""
Set name of a transition in the model.

# Examples

```julia
@name_transition net 1 = "name"
@name_transition net name = "transition_name"
@name_transition net "name" = "transition_name"
```
"""
macro name_transition(acsex, exs...)
    call = :(
        begin end
    )
    for ex in exs
        call_ = if ex.args[1] isa Number
            :($(esc(acsex))[$(ex.args[1]), :transName] = $(QuoteNode(ex.args[2])))
        else
            quote
                net = $(esc(acsex))
                ixs = findall(
                    i -> string(net[i, :transName]) == $(string(ex.args[1])),
                    row_ids(net, :T),
                )
                foreach(i -> net[i, :transName] = $(string(ex.args[2])), ixs)
            end
        end

        push!(call.args, call_)
    end

    return call
end

function incident_pattern(pattern, attr)
    ix = []
    for i = 1:length(attr)
        !isnothing(attr[i]) &&
            (
                m = match(pattern, string(attr[i]));
                !isnothing(m) && (string(attr[i]) == m.match)
            ) &&
            push!(ix, i)
    end

    return ix
end

function mode!(net, dict)
    for (spex, mods) in dict
        i = if spex isa Regex
            incident_pattern(spex, net[:, :specName])
        else
            find_rows(net, Symbol(spex), :specName)
        end

        for ix in i
            isnothing(net[ix, :specModality]) && (net[ix, :specModality] = Set{Symbol}())
            union!(net[ix, :specModality], mods)
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
@mode net (r"proj\\w+", r"experimental\\w+") conserved
@mode net (S, I) conserved
@mode net S conserved
```
"""
macro mode(acsex, spexs, mexs)
    mods = !isa(mexs, Expr) ? [mexs] : collect(mexs.args)
    spexs = isexpr(spexs, :tuple) ? spexs.args : [spexs]
    exs = map(ex -> striplines(ex), spexs)

    return quote
        dictcall = Dict()
        exs_ = []
        foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
        exs__ = []
        foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
        foreach(
            ex -> push!(dictcall, get_pattern(recursively_expand_dots(ex)) => $mods),
            exs__,
        )

        mode!($(esc(acsex)), dictcall)
    end
end

function set_valuation!(net, dict, valuation_type)
    for (spex, val) in dict
        i = if spex isa Regex
            incident_pattern(spex, column(net, :specName))
        else
            find_rows(net, Symbol(spex), :specName)
        end

        foreach(
            ix ->
                net[ix, Symbol(:spec, Symbol(uppercasefirst(string(valuation_type))))] =
                    eval(val),
            i,
        )
    end
end

export @cost, @reward, @valuation
for valuation_type in (:cost, :reward, :valuation)
    eval(
        quote
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
                foreach(
                    ex -> push!(
                        dictcall,
                        get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2],
                    ),
                    exs__,
                )
                return :(set_valuation!(
                    $(esc(acsex)),
                    $dictcall,
                    $(QuoteNode($(QuoteNode(valuation_type)))),
                ))
            end
        end,
    )
end

"""
Add new species to a model.

# Examples

```julia
@add_species net S I R
```
"""
macro add_species(acsex, exs...)
    call = :(
        begin end
    )
    spexs_ = []
    foreach(s -> push!(spexs_, s), exs)

    for ex in recursively_expand_dots.(spexs_)
        push!(call.args, :(add_row!($(esc(acsex)), :S; specName = $(QuoteNode(ex)))))
    end

    push!(call.args, :(assign_defaults!($(esc(acsex)))))
    return call
end

get_pattern(ex) = ex isa Expr && (macroname(ex) == :r_str) ? eval(ex) : ex

"""
Set initial values of species in an network.

# Examples

```julia
@prob_init net X = 1 Y = 2 Z = h(α)
@prob_init net [1.0, 2.0, 3.0]
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
            foreach(
                ex -> push!(
                    dictcall,
                    get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2],
                ),
                exs__,
            )

            init!($(esc(acsex)), dictcall)
        end
    end
end

# deprecate
macro prob_init_from_vec(acsex, vecex)
    return :(init!($(esc(acsex)), $(esc(vecex))))
end

function init!(net, inits)
    if inits isa AbstractVector && length(inits) == nrows(net, :S)
        column(net, :specInitVal) .= inits
    elseif inits isa AbstractDict
        for (k, init_val) in inits
            if k isa Number
                net[k, :specInitVal] = init_val
            else
                begin
                    i = if k isa Regex
                        incident_pattern(k, column(net, :specName))
                    else
                        find_rows(net, k, :specName)
                    end
                    foreach(ix -> (net[ix, :specInitVal] = init_val), i)
                end
            end
        end
    end

    return net
end

"""
Set uncertainty in initial values of species in an network (stderr).

# Examples

```julia
@prob_uncertainty net X = 0.1 Y = 0.2
@prob_uncertainty net [0.1, 0.2]
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
            foreach(
                ex -> push!(
                    dictcall,
                    get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2],
                ),
                exs__,
            )

            uncinit!($(esc(acsex)), dictcall)
        end
    end
end

function uncinit!(net, inits)
    inits isa AbstractVector &&
        length(inits) == nrows(net, :S) &&
        (column(net, :specInitUncertainty) .= inits; return)
    inits isa AbstractDict && for (k, init_val) in inits
        if k isa Number
            net[k, :specInitUncertainty] = init_val
        else
            begin
                i = if k isa Regex
                    incident_pattern(k, column(net, :specName))
                else
                    find_rows(net, k, :specName)
                end
                foreach(ix -> (net[ix, :specInitUncertainty] = init_val), i)
            end
        end
    end

    return net
end

function set_params!(net, params)
    return params isa AbstractDict && for (k, init_val) in params
        k = get_pattern(k)
        if k isa Regex
            i = incident_pattern(k, column(net, :prmName))
        else
            i = find_rows(net, k, :prmName)
            isempty(i) && (i = add_row!(net, :P; prmName = k))
        end

        foreach(ix -> net[ix, :prmVal] = eval(init_val), i)
    end
end

"""
Set parameter values in an network.

# Examples

```julia
@prob_params net α = 1.0 β = 2.0
```
"""
macro prob_params(acsex, exs...)
    exs = map(ex -> striplines(ex), exs)

    return quote
        dictcall = Dict()
        exs_ = []
        foreach(s -> push!(exs_, striplines(blockize(s))), $(QuoteNode(exs)))
        exs__ = []
        foreach(s -> foreach(s -> push!(exs__, s), s.args), exs_)
        foreach(
            ex -> push!(
                dictcall,
                get_pattern(recursively_expand_dots(ex.args[1])) => ex.args[2],
            ),
            exs__,
        )

        set_params!($(esc(acsex)), dictcall)
    end
end

meta!(net, metas) =
    for (k, metaval) in metas
        i = find_rows(net, k, :metaKeyword)
        isempty(i) && (i = add_row!(net, :M; metaKeyword = k))
        set_cell!(net, first(i), :metaVal, metaval)
    end

"""
Set model metadata (e.g. solver arguments)

# Examples

```julia
@prob_meta net tspan = (0, 100.0) schedule = schedule_weighted!
@prob_meta sir_acs tspan = 250 dt = 1   # `tstep` is a deprecated alias for `dt`
```
"""
macro prob_meta(acsex, exs...)
    dictcall = :(Dict([]))
    foreach(ex -> push!(dictcall.args[2].args, (ex.args[1] => eval(ex.args[2]))), exs)

    return :(meta!($(esc(acsex)), $dictcall))
end

"""
Alias object name in an net.

# Default names

| name       | short name |
|:---------- |:---------- |
| species    | S          |
| transition | T          |
| action     | A          |
| event      | E          |
| param      | P          |
| meta       | M          |

# Examples

```julia
@aka net species = resource transition = reaction
```
"""
macro aka(acsex, exs...)
    dictcall = :(Dict([]))
    foreach(
        ex -> push!(
            dictcall.args[2].args,
            Symbol("alias_", findfirst(==(ex.args[1]), alias_default)) => ex.args[2],
        ),
        exs,
    )
    return :(meta!($(esc(acsex)), $dictcall))
end

alias_default = Dict(
    :S => :species,
    :T => :transition,
    :A => :action,
    :E => :event,
    :P => :param,
    :M => :meta,
)

function get_alias(net, ob)
    return (
        i = find_rows(net, Symbol(:alias_, ob), :metaKeyword);
        !isempty(i) ? net[first(i), :metaVal] : alias_default[ob]
    )
end

# NOTE (WS-4 housekeeping): `@prob_check_verbose` was REMOVED here (definition + export above). It
# called an undefined `check_params`, so it threw at call time; it was vestigial — referenced by no
# code, test, demo, or doc, and the CONTRACT specifies no param-completeness check macro. Removed
# rather than back-filled with a `check_params` implementation nothing consumes.

"""
Add a periodic callback to a model.

# Examples

```julia
@periodic net 1.0 X += 1
```
"""
macro periodic(acsex, pex, acex)
    return push_to_acs!(acsex, Expr(:&&, :(@periodic($pex)), acex))
end

"""
Add a jump process (with specified Poisson intensity per unit time step) to a model.

# Examples

```julia
@jump net λ Z += rand(state.rng, Poisson(1.0))
```
"""
macro jump(acsex, inex, acex)
    # The Poisson intensity draws from the state-owned RNG (§4 D2/D5); `state` is in scope
    # because the generated trigger is compiled into a (state, transition) closure.
    return push_to_acs!(
        acsex,
        Expr(:&&, Expr(:call, :rand, :(state.rng), :(Poisson(max(state.dt * $inex, 0)))), acex),
    )
end

"""
Evaluate expression in ReactiveDynamics scope.

# Examples

```julia
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t / 10^3)
```
"""
macro register(ex)
    return :(@eval ReactiveDynamics $ex)
end
