# in-macro relative reference
const refs = [:network, :this]

"Interprets references within an expression."
function interpret(expr, all_species = [])
    expr = prewalk(expr) do e
        if e isa Symbol && (string(e) ∈ all_species)
            Expr(:macrocall, :a_str, string(e))
        else
            e
        end
    end

    expr = prewalk(expr) do e
        # resolve symbolic references (macro form)
        return if macroname(e) ∈ refs
            macroname(e)
        elseif macroname(e) == :a_str
            :(getagent(network, $(e.args[3]))[]) # ref
        elseif macroname(e) == :getagent
            :(getagent(network, $(e.args[3])))
        else
            e
        end
    end

    ref_aagents(expr, :network)
end

"Interpret a-string as agent references."
function ref_aagents(expr, agent = :network)
    prewalk(expr) do e
        if macroname(e) == :a_str
            :(getpath($agent, $(e.args[3])))
        else
            e
        end
    end
end

"Interprets references within an expression and returns an evaluated function of `(network, this)`."
function interpret_eval(expr, all_species = [])
    expr = interpret(expr, all_species)

    eval(:((network, this) -> $expr))
end
