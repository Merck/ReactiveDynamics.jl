# reaction network DSL: UPDATE part; add place, name, add modalities, set model variables, set solver arguments

export @push, @name_transition, @mode, @add_place
export @periodic, @jump
export @prob_init, @prob_uncertainty, @prob_params, @prob_meta
export @aka
export @register

# NOTE (WS-4 housekeeping): three dangling exports were DELETED from here — `@prob_role`,
# `@list_by_role`, `@list_roles` (a legacy roles/actors ontology that was never implemented: no
# macro definitions, no `placeRole` schema attribute), and `@prob_check_verbose` (see below). The
# role concept was dropped; ADR 0009's `PortRole` is an unrelated per-Place `role` field authored
# inside `@reaction_network`, not a `@prob_role`-style config macro, so nothing is repurposed.

using DataFrames
using MacroTools: striplines

"""
Build the expression that appends reactions/attributes to the network bound to `netex`. Accepts either a
`begin … end` block of reaction lines or a single line with trailing `key = value` attributes (folded
into `key => value` pairs), then splices the parsed data into a `merge_network!` call. The shared engine
behind [`@push`](@ref) and [`@periodic`](@ref).
"""
function push_to_network!(netex, exs...)
    if isexpr(exs[1], :block)
        ex = striplines(exs[1])
    else
        args = Any[]
        map(
            el -> if isexpr(el, :(=))
                push!(args, :($(el.args[1]) => $(el.args[2])))
            else
                push!(args, el)
            end, exs
        )
        ex = Expr(:tuple, args...)
    end

    return quote
        ex = blockize($(QuoteNode(ex)))
        merge_network!($(esc(netex)), get_data(ex)...)
    end
end

"""
Add reactions to a network.

# Examples

```julia
@push sir β * S * I * tdecay(@time()) S + I --> 2I name => SI2I
@push sir begin
    ν * I, I --> R, name => I2R
    γ, R --> S, name => R2S
end
```
"""
macro push(netex, exs...)
    return push_to_network!(netex, exs...)
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
macro name_transition(netex, exs...)
    call = :(
        begin end
    )
    for ex in exs
        call_ = if ex.args[1] isa Number
            :($(esc(netex))[$(ex.args[1]), :transName] = $(QuoteNode(ex.args[2])))
        else
            quote
                net = $(esc(netex))
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

"""
Row indices of `attr` whose stringified value matches `pattern` in FULL (the regex must span the whole
name, not just a substring). The regex-valued analogue of [`find_rows`](@ref), used by [`mode!`](@ref) to
apply a modality to every place whose name matches a pattern.
"""
function incident_pattern(pattern, attr)
    ix = []
    for i in 1:length(attr)
        !isnothing(attr[i]) &&
            (
            m = match(pattern, string(attr[i]));
            !isnothing(m) && (string(attr[i]) == m.match)
        ) &&
            push!(ix, i)
    end

    return ix
end

"""
Union the modality tags in `dict` (`place-or-regex => modalities`) into each matching place'
`placeDefaultModality` set, in place. A plain key matches one place by name ([`find_rows`](@ref)); a `Regex`
key matches every place whose name matches ([`incident_pattern`](@ref)). The runtime behind
[`@mode`](@ref).
"""
function mode!(net, dict)
    for (spex, mods) in dict
        i = if spex isa Regex
            incident_pattern(spex, net[:, :placeName])
        else
            find_rows(net, Symbol(spex), :placeName)
        end

        for ix in i
            isnothing(net[ix, :placeDefaultModality]) && (net[ix, :placeDefaultModality] = Set{Symbol}())
            union!(net[ix, :placeDefaultModality], mods)
        end
    end
    return
end

"""
Set place modality.

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
macro mode(netex, spexs, mexs)
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

        mode!($(esc(netex)), dictcall)
    end
end

"""
Set the `valuation_type` economic attribute (`:cost`/`:reward`/`:valuation` → the `placeCost`/`placeReward`/
`placeValuation` column) of each place named in `dict` (`place-or-regex => value`), in place. Matches
by name ([`find_rows`](@ref)) or, for a `Regex` key, by pattern ([`incident_pattern`](@ref)). The runtime
behind the `@cost`/`@reward`/`@valuation` macros.
"""
function set_valuation!(net, dict, valuation_type)
    for (spex, val) in dict
        i = if spex isa Regex
            incident_pattern(spex, column(net, :placeName))
        else
            find_rows(net, Symbol(spex), :placeName)
        end

        foreach(
            ix ->
            net[ix, Symbol(:place, Symbol(uppercasefirst(string(valuation_type))))] =
                eval(val),
            i,
        )
    end
    return
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
            macro $valuation_type(netex, exs...)
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
                return :(
                    set_valuation!(
                        $(esc(netex)),
                        $dictcall,
                        $(QuoteNode($(QuoteNode(valuation_type)))),
                    )
                )
            end
        end,
    )
end

"""
Add new place to a model.

# Examples

```julia
@add_place net S I R
```
"""
macro add_place(netex, exs...)
    call = :(
        begin end
    )
    spexs_ = []
    foreach(s -> push!(spexs_, s), exs)

    for ex in recursively_expand_dots.(spexs_)
        push!(call.args, :(add_row!($(esc(netex)), :S; placeName = $(QuoteNode(ex)))))
    end

    push!(call.args, :(assign_defaults!($(esc(netex)))))
    return call
end

# Resolve a place/param selector to a matcher: an `r"…"` regex-string macrocall becomes the compiled
# `Regex` (so callers can pattern-match names), anything else passes through as the literal name.
get_pattern(ex) = ex isa Expr && (macroname(ex) == :r_str) ? eval(ex) : ex

"""
Set initial values of place in a network.

# Examples

```julia
@prob_init net X = 1 Y = 2 Z = h(α)
@prob_init net [1.0, 2.0, 3.0]
```
"""
macro prob_init(netex, exs...)
    exs = map(ex -> striplines(ex), exs)

    return if length(exs) == 1 && (isexpr(exs[1], :vect) || (exs[1] isa Symbol))
        :(init!($(esc(netex)), $(esc(exs[1]))))
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

            init!($(esc(netex)), dictcall)
        end
    end
end

# deprecate
macro prob_init_from_vec(netex, vecex)
    return :(init!($(esc(netex)), $(esc(vecex))))
end

"""
Set place initial values (`placeInitVal`) from `inits`, in place. A vector matching the place count is
assigned positionally; a dict maps `index-or-name-or-regex => value` (a `Regex` key sets every matching
place). The runtime behind [`@prob_init`](@ref).
"""
function init!(net, inits)
    if inits isa AbstractVector && length(inits) == nrows(net, :S)
        column(net, :placeInitVal) .= inits
    elseif inits isa AbstractDict
        for (k, init_val) in inits
            if k isa Number
                net[k, :placeInitVal] = init_val
            else
                begin
                    i = if k isa Regex
                        incident_pattern(k, column(net, :placeName))
                    else
                        find_rows(net, k, :placeName)
                    end
                    foreach(ix -> (net[ix, :placeInitVal] = init_val), i)
                end
            end
        end
    end

    return net
end

"""
Set uncertainty in initial values of place in a network (stderr).

# Examples

```julia
@prob_uncertainty net X = 0.1 Y = 0.2
@prob_uncertainty net [0.1, 0.2]
```
"""
macro prob_uncertainty(netex, exs...)
    exs = map(ex -> striplines(ex), exs)

    return if length(exs) == 1 && (isexpr(exs[1], :vect) || (exs[1] isa Symbol))
        :(uncinit!($(esc(netex)), $(esc(exs[1]))))
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

            uncinit!($(esc(netex)), dictcall)
        end
    end
end

"""
Set place initial-value uncertainty (`placeInitUncertainty`, a stderr) from `inits`, in place — the
uncertainty counterpart of [`init!`](@ref), with the same vector/dict/regex handling. The runtime behind
[`@prob_uncertainty`](@ref).
"""
function uncinit!(net, inits)
    inits isa AbstractVector &&
        length(inits) == nrows(net, :S) &&
        (column(net, :placeInitUncertainty) .= inits; return)
    inits isa AbstractDict && for (k, init_val) in inits
        if k isa Number
            net[k, :placeInitUncertainty] = init_val
        else
            begin
                i = if k isa Regex
                    incident_pattern(k, column(net, :placeName))
                else
                    find_rows(net, k, :placeName)
                end
                foreach(ix -> (net[ix, :placeInitUncertainty] = init_val), i)
            end
        end
    end

    return net
end

"""
Set parameter values (`prmVal`) from the dict `params` (`name-or-regex => value`), in place, ADDING a
`:P` row for a plain name that does not yet exist. A `Regex` key sets every matching existing param. The
runtime behind [`@prob_params`](@ref).
"""
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
Set parameter values in a network.

# Examples

```julia
@prob_params net α = 1.0 β = 2.0
```
"""
macro prob_params(netex, exs...)
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

        set_params!($(esc(netex)), dictcall)
    end
end

"""
Set network metadata (`:M` rows) from the dict `metas` (`keyword => value`), in place, adding a row for
an unseen keyword and overwriting an existing one. The runtime behind [`@prob_meta`](@ref) and
[`@aka`](@ref).
"""
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
@prob_meta sir tspan = 250 dt = 1   # `tstep` is a deprecated alias for `dt`
```
"""
macro prob_meta(netex, exs...)
    dictcall = :(Dict([]))
    foreach(ex -> push!(dictcall.args[2].args, (ex.args[1] => eval(ex.args[2]))), exs)

    return :(meta!($(esc(netex)), $dictcall))
end

"""
Alias an object name in a network.

# Default names

| name       | short name |
|:---------- |:---------- |
| place      | S          |
| transition | T          |
| action     | A          |
| event      | E          |
| param      | P          |
| meta       | M          |

The pre-ADR-0017 spelling `species` is still accepted on the left of `=` for one release (it selects
the `:S` object exactly as `place` does) and then removed.

# Examples

```julia
@aka net place = resource transition = reaction
```
"""
macro aka(netex, exs...)
    dictcall = :(Dict([]))
    foreach(
        ex -> push!(
            dictcall.args[2].args,
            Symbol("alias_", _aka_object(ex.args[1])) => ex.args[2],
        ),
        exs,
    )
    return :(meta!($(esc(netex)), $dictcall))
end

alias_default = Dict(
    :S => :place,
    :T => :transition,
    :A => :action,
    :E => :event,
    :P => :param,
    :M => :meta,
)

# The object an `@aka` left-hand name selects. Normally the inverse of `alias_default`; the extra
# entry keeps the retired ADR-0017 spelling working for ONE release, so `@aka net species = resource`
# still targets `:S` instead of silently becoming `alias_nothing`.
const _AKA_LEGACY_NAMES = Dict(:species => :S)
_aka_object(name) = get(_AKA_LEGACY_NAMES, name, findfirst(==(name), alias_default))

"""
The display alias for object `ob` (`:S`/`:T`/`:A`/`:E`/`:P`/`:M`) — a user-set `alias_<ob>` metadata
value if one was authored via [`@aka`](@ref), else the built-in default from `alias_default`.
"""
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
macro periodic(netex, pex, acex)
    return push_to_network!(netex, Expr(:&&, :(@periodic($pex)), acex))
end

"""
Add a jump process (with specified Poisson intensity per unit time step) to a model.

# Examples

```julia
@jump net λ Z += rand(state.rng, Poisson(1.0))
```
"""
macro jump(netex, inex, acex)
    # The Poisson intensity draws from the state-owned RNG (§4 D2/D5); `state` is in scope
    # because the generated trigger is compiled into a (state, transition) closure.
    return push_to_network!(
        netex,
        Expr(:&&, Expr(:call, :rand, :(state.rng), :(Poisson(max(state.dt * $inex, 0)))), acex),
    )
end

"""
Register a host function/definition into `ReactiveDynamics` scope so model expressions may call it by
name (e.g. a custom rate/guard helper). Evaluates `ex` at macro-expansion time — an AUTHORING-time escape
hatch, distinct from the eval-free runtime; do not use it to inject per-run data.

# Examples

```julia
@register bool_cond(t) = (100 < t < 200) || (400 < t < 500)
@register tdecay(t) = exp(-t / 10^3)
```
"""
macro register(ex)
    return :(@eval ReactiveDynamics $ex)
end
