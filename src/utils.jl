# various utility functions

"If `ex` is a macrocall, return the macro's name, else return `nothing`."
macroname(e) = isexpr(e, :macrocall) ? Symbol(strip(string(e.args[1]), '@')) : nothing

"Check if a object is a valid agent reference (string, symbol, `a`-string)."
function ispath(o)
    isa(o, Union{Symbol, AbstractString}) || 
        (isexpr(o, :macrocall) && (macroname(o) == :a_str))
end

"Join nested block expressions in `ex`."
function flatten!(ex)
    isexpr(ex, :block) || return ex

    args = []; for ex in ex.args
        isexpr(ex, :block) ? append!(args, ex.args) : push!(args, ex)
    end

    ex.args = args; ex
end

"Return a copy of `expr` where, for each pair `old=>new` in `old_new`, all occurrences of `old` are replaced by `new`."
function replace_in_expr(expr, old_new...)
    old_new = Dict(old_new...)

    prewalk(ex -> haskey(old_new, ex) ? old_new[ex] : ex, expr)
end

"Convert path refs in an expression to strings."
function make_paths(ex, names)
    replacements = map(filter(ispath, names)) do sym
        sym => to_string(sym)
    end

    replace_in_expr(ex, replacements...)
end

"Wrap an expression into a block expression."
blockize(ex) = striplines(isexpr(ex, :block) ? ex : Expr(:block, ex))

"""
Take a random sample from a weighted range. If an item is a `Pair`,
the first argument provides the sampling weight (defaults to 1).
"""
function rand_polyrange(rng)
    isempty(rng) && return missing
    r = rand() * sum(r -> r isa Tuple ? r[1] : 1, rng)
    
    ix = 0; s = 0; while s <= r && (ix < length(rng))
        ix += 1
        s += rng[ix] isa Tuple ? rng[ix][1] : 1
    end

    r = rng[ix] isa Tuple ? rng[ix][2] : rng[ix]
    r isa Sampleable ? rand(r) : r
end

"""
Take a random sample from a weighted range.

If an item is a `Pair`, the first argument provides
the sampling weight (defaults to 1) evaluated as a function of `(network, this)`.
"""
function rand_polyrange(rng, network, this)
    isempty(rng) && return missing
    ws = map(r -> r isa Tuple ? interpret_eval(r[1])(network, this) : 1, rng)
    r = rand() * sum(ws)
    
    ix = 0; s = 0; while s <= r && ix < length(rng)
        ix += 1; s += ws[ix]
    end

    r = rng[ix] isa Tuple ? rng[ix][2] : rng[ix]
    r isa Sampleable ? rand(r) : r
end

"Return a string representation of a path."
function to_string end

to_string(o::Symbol)::AbstractString = string(o)
to_string(o::AbstractString)::AbstractString = o

function to_string(o::Expr)
    if isexpr(o, :macrocall) && (macroname(o) == :a_str)
        o.args[end]
    else @error("could not parse agent's path $o") end
end

"Decorated operadic agent's path."
macro a_str(p::AbstractString) normpath(p) end

"""
Set up nested hieararchy of agents according to `path`.
If an intermediate agent doesn't exist, create an instance of `FreeAgent`.
"""
function mkpath!(a::AbstractAlgebraicAgent, path::String)
    path = normpath(path)

    for p in splitpath(path)
        (p ∈ [".", "/", ""]) && continue
        (p == "..") && @error "Can't backtrack outside of the reaction network!"

        a = haskey(inners(a), p) ? inners(a)[p] : entangle!(a, FreeAgent(p))
    end

    a
end