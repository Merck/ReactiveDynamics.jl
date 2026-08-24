# Assortment of expression handling utilities

using MacroTools: striplines

"""
Collect the `key = value` expressions in `args` into a vector of `:kw` Exprs (`Expr(:kw, key, value)`),
dropping any non-assignment entries. Used to lift a macro's trailing `k = v` arguments into real keyword
arguments for a generated call.
"""
function kwargize(args)
    kwargs = Any[]
    map(el -> isexpr(el, :(=)) && push!(kwargs, Expr(:kw, el.args[1], el.args[2])), args)

    return kwargs
end

"""
Partition `args` into `(positional, keywords)`: `key = value` entries become `:kw` Exprs (the second
return), everything else is `esc`aped and kept as a positional argument (the first). The escaping means
this is intended for building a call inside a macro. Returns the two vectors as a tuple.
"""
function args_kwargs(args)
    args_ = Any[]
    kwargs = Any[]
    map(
        el -> if isexpr(el, :(=))
            push!(kwargs, Expr(:kw, el.args[1], el.args[2]))
        else
            push!(args_, esc(el))
        end, args
    )

    return args_, kwargs
end

"""
Return the value expression stored for `key` in a collection of `key = value` keyword expressions AND
delete that entry, or return `default` if the key is absent. The mutating counterpart of
[`get_kwarg`](@ref) — used to consume a keyword out of a macro's argument list.
"""
function find_kwargex_delete!(collection, key, default = :())
    ix = findfirst(ex -> ex.args[1] == key, collection)
    return if !isnothing(ix)
        (v = collection[ix].args[2]; deleteat!(collection, ix); v)
    else
        default
    end
end

"""
The bare macro name of a macrocall Expr, `@` stripped (`:(@foo x)` → `:foo`); errors on a non-macrocall.
"""
macroname(ex) =
if isexpr(ex, :macrocall)
    (str = string(ex.args[1]); Symbol(strip(str, '@')))
else
    error("expr $ex is not a macrocall")
end
"Strip a leading `@` from a Symbol (`Symbol(\"@foo\")` → `:foo`); the Symbol-input twin of [`macroname`](@ref)."
strip_sym(sym) = (str = string(sym); Symbol(strip(str, '@')))

"Wrap `ex` in a `:block` (unless it already is one) and strip line-number nodes — normalizes a macro body to a clean block."
blockize(ex) = striplines(isexpr(ex, :block) ? ex : Expr(:block, ex))

"Quote a Symbol as a `QuoteNode` so it survives interpolation as data (a literal `:sym`), passing non-Symbols through unchanged."
preserve_sym(el) = el isa Symbol ? QuoteNode(el) : el

"`true` if `collection[ix]` is both assigned (not `#undef`) and not `missing` — a guard for sparse/optional columns."
assigned(collection, ix) = isassigned(collection, ix) && !ismissing(collection[ix])

"""
Flatten one level of nested `:block` expressions into `ex`'s own args (a shallow, non-recursive splice),
so a block of blocks reads as a single flat statement list. Mutates and returns `ex`.
"""
function unblock_shallow!(ex)
    isexpr(ex, :block) || return ex
    args = []
    for ex in ex.args
        isexpr(ex, :block) ? append!(args, ex.args) : push!(args, ex)
    end
    ex.args = args
    return ex
end

"""
Flatten a compound name to a single symbol: rewrite `.` → `__` and drop parentheses (`a.b` → `:a__b`),
passing numbers through unchanged. The primitive behind [`recursively_expand_dots`](@ref) — turns dotted
place notation into one atomic place name.
"""
function underscorize(ex)
    return if ex isa Number
        ex
    else
        (str = string(ex); replace(str, '.' => "__", '(' => "", ')' => "") |> Symbol)
    end
end

"Keys of a mixed iterable of `Pair`s and bare values — a `Pair`'s first element, else the element itself."
wkeys(itr) = map(x -> x isa Pair ? x[1] : x, itr)
"Values of a mixed iterable of `Pair`s and bare values — a `Pair`'s second element, else the element itself."
wvalues(itr) = map(x -> x isa Pair ? x[2] : x, itr)

"""
Rewrite the `key => _` entry of a pair-collection `col` in place to `key => val`, matching keys via
[`wkeys`](@ref); a no-op if `key` is absent. Returns `col`.
"""
function wset!(col, key, val)
    ix = findfirst(==(key), wkeys(col))
    !isnothing(ix) && (col[ix] = (key => val))

    return col
end

"""
Return the (expression) value stored for `key` in a collection of keyword expressions, or `default` if
the key is absent. The non-mutating counterpart of [`find_kwargex_delete!`](@ref): it reads the mapping
in place, leaving `collection` untouched.
"""
function get_kwarg(collection, key, default = :())
    ix = findfirst(
        ex -> ex isa Expr && hasproperty(ex, :args) && (ex.args[1] == key),
        collection,
    )

    return !isnothing(ix) ? collection[ix].args[2] : default
end

"""
The structured tokens bound to `place` by a firing `transition` — filters the transition's
`bound_tokens` to those whose `place` is `place`. Returns `nothing` when there is no
transition or nothing is bound. Used by action/predicate evaluation to resolve `@field`-style token reads.
"""
function get_bound_agent(transition, place)
    return if !isnothing(transition) && !isempty(transition.bound_tokens)
        bound_agents = filter(x -> x.place == place, transition.bound_tokens)

        bound_agents
    end
end
