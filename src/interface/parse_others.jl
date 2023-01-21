"""
Parse events; expects `event_body` or `(event_body, kwargs)` tuples.

Return a named tuple `(event_body, params)`.
"""
function parse_events(expr)
    objects = []
    if isexpr(expr, :block)
        expr = striplines(expr)
        esc_dollars!(expr)
        append!(objects, expr.args)
    elseif expr != :()
        push!(objects, expr)
    end

    map(objects) do event
        exs = isexpr(event, :tuple) ? event.args : [event]

        params = Dict()
        for prm in exs[2:end]
            (prm isa Expr && isexpr(prm, :call)) || continue
            key = findfirst(k -> prm.args[2] ∈ k, prettynames[:event]) # localize parameters

            !isnothing(key) && push!(params, key => prm.args[3])
        end

        (; event_body = exs[1], params)
    end
end

"""
Parse sampleables; expects `range, params...` tuples.

Return a named tuple `(range, params)`.
"""
function parse_observables(expr)
    objects = []
    if isexpr(expr, :block)
        expr = striplines(expr)
        esc_dollars!(expr)
        append!(objects, expr.args)
    elseif expr != :()
        push!(objects, expr)
    end

    map(objects) do sampleable
        exs = isexpr(sampleable, :tuple) ? sampleable.args : [sampleable]
        params = Dict()
        for prm in exs[2:end]
            (prm isa Expr && isexpr(prm, :call)) || continue
            key = findfirst(k -> prm.args[2] ∈ k, prettynames[:sampleable]) # localize parameters

            push!(params, (isnothing(key) ? prm.args[2] : key) => prm.args[3])
        end

        (; range = exs[1], params)
    end
end
