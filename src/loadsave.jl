# Model + solution persistence. The model format is the single eval-free JSON (`model.rdj.json`,
# ADR 0005): `@import_model`/`@export_model` delegate to the typed-IR loader (src/serialize.jl),
# which NEVER `Meta.parse`s or `eval`s a model field — closing the import-time RCE that the old
# TOML/CSV loader had (it `eval`'d attribute strings and whole `registered` function bodies). The
# TOML/CSV/JLD2 model zoo is removed per ADR 0005. Solution OUTPUTS (the `sol` DataFrame) are a
# separate concern and may still be written to CSV for human inspection.

export @import_model, @export_model
export @export_solution_as_table, @export_solution_as_csv

using DataFrames
using CSV

"""
    @import_model "model.rdj.json"
    @import_model "model.rdj.json" prob   registry=Dict(:Kind=>ctor)

Load a model from an eval-free `model.rdj.json` document (ADR 0005) and construct a
`ReactionNetworkProblem`. The document is inert data: it is validated and lowered through the
typed ExprNode IR, never `eval`'d. Host token kinds / callbacks are supplied BY NAME through the
`registry` (ADR 0006 §C); the file alone cannot execute code.

# Examples

```julia
@import_model "pipeline.rdj.json" prob
prob = from_json_model(read("pipeline.rdj.json", String); seed = 1, registry = REG)
```
"""
macro import_model(pathex, name = gensym(), kwargs...)
    return :(
        $(esc(name)) =
            from_json_model(read($(esc(pathex)), String); $(map(esc, kwargs)...))
    )
end

"""
    @export_model prob "model.rdj.json"

Serialize a model to an eval-free `model.rdj.json` document (ADR 0005).

# Examples

```julia
@export_model prob "pipeline.rdj.json"
```
"""
macro export_model(probex, pathex)
    return :(write($(esc(pathex)), to_json_model($(esc(probex)))))
end

# ── Solution outputs (CSV; the model is the reproducible input, the solution is output) ──

"""
    @export_solution_as_table sol

Export a solution's trajectory as a `DataFrame`.
"""
macro export_solution_as_table(solex, pathex = "sol")
    return :(DataFrame($(esc(solex)).sol))
end

get_DataFrame(sol) = sol.sol

"""
    @export_solution_as_csv sol "sol.csv"

Export a solution's trajectory to a CSV file.
"""
macro export_solution_as_csv(solex, pathex = "sol.csv")
    return :(CSV.write($(string(pathex)), get_DataFrame($(esc(solex)))))
end
