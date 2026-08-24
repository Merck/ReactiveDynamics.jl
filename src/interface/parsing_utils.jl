# Reaction-line parsing helpers. Parts of this file were originally adapted from Catalyst.jl.

"""
Flatten dotted notation `A.B` into the single place symbol `Symbol("A__B")` (an alias for
[`underscorize`](@ref)). Applied to arc terms so a dotted place name becomes one atomic symbol.
"""
recursively_expand_dots(ex) = underscorize(ex)

"""
The number of elements in a tuple expression, or `1` for a non-tuple `ex` (a bare Symbol/number is
treated as a length-1 "tuple"). Paired with [`get_tup_arg`](@ref) to read multiplicity/multiplicity
terms uniformly whether or not they are written as a tuple.
"""
function tup_leng(ex::SampleableValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

"""
The `i`-th element of a tuple expression, or `ex` itself if it is not a tuple (a bare Symbol/number).
The accessor counterpart of [`tup_leng`](@ref).
"""
function get_tup_arg(ex::SampleableValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

"""
Combine a base multiplicity `mult` with additional factors `mults...` into a single multiplicity term,
folding all numeric factors into one constant and preserving symbolic ones as a `*` product. A purely
numeric set multiplies to a number; otherwise returns the simplified product Expr (dropping a redundant
`1` coefficient). Flattens nested `*` products via [`recursively_find_mults!`](@ref).
"""
function multiplex(mult, mults...)
    all(m -> isa(m, Number), [mult] ∪ mults) && return mult * prod(mults; init = 1.0)
    multarray = SampleableValues[]
    recursively_find_mults!(multarray, mults...)
    mults_numeric =
        prod(filter(m -> isa(m, Number), multarray); init = 1.0) *
        (mult isa Number ? mult : 1.0)
    mults_expr = filter(m -> !isa(m, Number), multarray)
    mult = mult isa Expr ? deepcopy(mult) : :(*())
    mults_numeric == 1 && length(mults_expr) == 1 && return mults_expr[1]
    mults_numeric != 1.0 && push!(mult.args, mults_numeric)
    append!(mult.args, mults_expr)

    return mult
end

"""
Collect the leaf factors of a (possibly nested) `*` product into `multarray`, descending into any
`*`-headed subexpression so `a * (b * c)` flattens to `[a, b, c]`. The recursion behind [`multiplex`](@ref).
"""
function recursively_find_mults!(multarray, mults...)
    for m in mults
        if isa(m, Expr) && m.args[1] == :*
            recursively_find_mults!(multarray, m.args[2:end]...)
        else
            push!(multarray, m)
        end
    end
    return
end
