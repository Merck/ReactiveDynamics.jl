# parts of the code were taken from Catalyst.jl and adapted

recursively_expand_dots(ex) = underscorize(ex)

# Returns the length of a expression tuple, or 1 if it is not an expression tuple (probably a  Symbol/Numerical).
function tup_leng(ex::SampleableValues)
    (typeof(ex) == Expr && ex.head == :tuple) && (return length(ex.args))
    return 1
end

#Gets the ith element in a expression tuple, or returns the input itself if it is not an expression tuple (probably a  Symbol/Numerical).
function get_tup_arg(ex::SampleableValues, i::Int)
    (tup_leng(ex) == 1) && (return ex)
    return ex.args[i]
end

function multiplex(mult, mults...)
    all(m -> isa(m, Number), [mult] ∪ mults) && return mult * prod(mults; init = 1.0)
    multarray = SampleableValues[]
    recursively_find_mults!(multarray, mults...)
    mults_numeric = prod(filter(m -> isa(m, Number), multarray); init = 1.0) *
                    (mult isa Number ? mult : 1.0)
    mults_expr = filter(m -> !isa(m, Number), multarray)
    mult = mult isa Expr ? deepcopy(mult) : :(*())
    mults_numeric == 1 && length(mults_expr) == 1 && return mults_expr[1]
    mults_numeric != 1.0 && push!(mult.args, mults_numeric)
    append!(mult.args, mults_expr)

    mult
end

function recursively_find_mults!(multarray, mults...)
    for m in mults
        isa(m, Expr) && m.args[1] == :* ?
        recursively_find_mults!(multarray, m.args[2:end]...) : push!(multarray, m)
    end
end
