# Assortment of expression handling utilities

using MacroTools: striplines

"Return array of keyword expressions in `args`."
function kwargize(args)
    kwargs = Any[]
    map(el -> isexpr(el, :(=)) && push!(kwargs, Expr(:kw, el.args[1], el.args[2])), args)

    kwargs
end

"""Split an array of expressions into arrays of "argument" expressions and keyword expressions."""
function args_kwargs(args)
    args_ = Any[]; kwargs = Any[]
    map(el -> isexpr(el, :(=)) ? push!(kwargs, Expr(:kw, el.args[1], el.args[2])) : push!(args_, esc(el)), args)

    args_, kwargs
end

"Return the (expression) value stored for the given key in a collection of keyword expression, or the given default value if no mapping for the key is present."
function find_kwargex_delete!(collection, key, default=:())
    ix = findfirst(ex -> ex.args[1] == key, collection)
    !isnothing(ix) ? (v = collection[ix].args[2]; deleteat!(collection, ix); v) : default
end

macroname(ex) = isexpr(ex, :macrocall) ? (str = string(ex.args[1]); Symbol(strip(str, '@'))) : error("expr $ex is not a macrocall")
strip_sym(sym) = (str = string(sym); Symbol(strip(str, '@')))

blockize(ex) = striplines(isexpr(ex, :block) ? ex : Expr(:block, ex))

preserve_sym(el) = el isa Symbol ? QuoteNode(el) : el

assigned(collection, ix) = isassigned(collection, ix) && !ismissing(collection[ix])

function unblock_shallow!(ex)
    isexpr(ex, :block) || return ex
    args = []; for ex in ex.args
        isexpr(ex, :block) ? append!(args, ex.args) : push!(args, ex)
    end
    ex.args = args; ex
end

underscorize(ex) = ex isa Number ? ex : (str = string(ex); replace(str, '.' => "__", '(' => "", ')' => "") |> Symbol)

wkeys(itr) = map(x -> x isa Pair ? x[1] : x, itr)
wvalues(itr) = map(x -> x isa Pair ? x[2] : x, itr)

function wset!(col, key, val)
    ix = findfirst(==(key), wkeys(col))
    !isnothing(ix) && (col[ix] = (key => val))

    col
end

"Return the (expression) value stored for the given key in a collection of keyword expression, or the given default value if no mapping for the key is present."
function get_kwarg(collection, key, default=:())
    ix = findfirst(ex -> ex isa Expr && hasproperty(ex, :args) && (ex.args[1] == key), collection)
    
    !isnothing(ix) ? collection[ix].args[2] : default
end