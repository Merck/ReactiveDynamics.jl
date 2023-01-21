# `_getparameters`, `_setparameters!`, `interpret_params!` default constructors
"""
Provide `_getparameters`, `_setparameters!`, `interpret_params!`,
assuming interpretable parameters under `params_expression`.
"""
macro params_interface(T)
    quote
        # parameter interpretation
        """
        Interpret agent's `a` parameters (in expression form) stored under key `params_expression`,
        and store the resulting objects under `params_interpreted`.
        """
        function interpret_params!(a::$T, all_species = String[])
            for (k, ex) in a.params_expression
                push!(a.params_interpreted, k => interpret_eval(ex, all_species))
            end

            a
        end

        # sample params
        "Evaluate (interpreted) parameters under `params_interpreted`."
        function sample_params(a::$T)
            rn = getnetwork(a)

            Dict(k => v(rn, a) for (k, v) in a.params_interpreted)
        end

        # AlgebraicAgents.jl parameter get/set interface
        function AlgebraicAgents._getparameters(o::$T)
            Dict(k => v(@args(o)...) for (k, v) in o.params_interpreted)
        end

        function AlgebraicAgents._setparameters!(o::$T, parameters)
            merge!(o.params_expression, parameters)

            interpret_params!(o, all_species(getnetwork(o)))
        end
    end |> esc
end
