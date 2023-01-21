# reaction network DSL: UPDATE part; add species, name, add modalities, set model variables, set solver arguments 

"""
Add a periodic callback to a model.

# Examples
```julia
@periodic acs 1. X += 1
```
"""
macro periodic(acsex, pex, acex)
    push_to_acs!(acsex, Expr(:&&, :(@periodic($pex)), acex))
end

"""
Add a jump process (with specified Poisson intensity per unit time step) to a model.

# Examples
```julia
@jump acs λ Z += rand(Poisson(1.))
```
"""
macro jump(acsex, inex, acex)
    push_to_acs!(acsex,
                 Expr(:&&,
                      Expr(:call, :rand, :(Poisson(max(@solverarg(:tstep) * $inex, 0)))),
                      acex))
end

"Evaluate expression in DyVE scope."
macro register(ex)
    :(@eval DyVE $ex)
end
