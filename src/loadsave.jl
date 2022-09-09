export @import_model, @export_model
export @import_solution, @export_solution
export @export_as_table, @export_csv
export @export, @import

using TOML, JLD2, CSV
using DataFrames

const objects_aliases = Dict(:S => "spec", :T => "trans", :P => "prm", :M => "meta", :E => "event", :obs => "obs")
const RN_attrs = string.(propertynames(ReactionNetwork().attrs))

get_attrs(object) = (object = object isa Symbol ? objects_aliases[object] : object; filter(x -> occursin(object, x), RN_attrs))

function serialize_to_toml(acs::ReactionNetwork, io::IO=stdout)
    data = Dict()
    for (key, val) in objects_aliases
        push!(data, val => []);
        for i in 1:nparts(acs, key)
            data_ = Dict()
            for attr in get_attrs(val)
                attr_val = acs[i, Symbol(attr)]
                ismissing(attr_val) && continue
                attr_val = attr_val isa Number ? attr_val : string(attr_val)
                push!(data_, string(attr) => attr_val)
            end
            push!(data[val], data_)
        end
    end

    TOML.print(io, data)
end

function deserialize_from_dict(data::Dict)
    acs = ReactionNetwork()
    for (key, val) in objects_aliases
        for row in data[val]
            i = add_part!(acs, key)
            for (attr, attrval) in row
                set_subpart!(acs, i, Symbol(attr), attrval)
            end
        end
    end

    for row in get(data, "registered", []) eval(Meta.parseall(row["body"])) end

    assign_defaults!(acs)
end

export_model(acs::ReactionNetwork, path::AbstractString) = open(path, "w") do io; serialize_to_toml(acs, io) end
import_model(path::AbstractString) = (dict = TOML.parsefile(path); deserialize_from_dict(dict))

function load_models(io::IO)
    for line in eachline(io)
        name, path = split(line, ';')
        name = isempty(name) ? gensym() : name
        
    end
end

"""
Export model to a file.

# Examples
```julia
@export_model acs "acs_data.toml"
```
"""
macro export_model(acsex, pathex) :(export_model($(esc(acsex)), $(string(pathex)))) end

"""
Import a model from a file.

# Examples
```julia
@import_model "model.toml"
```
"""
macro import_model(pathex, name=gensym()) :($(esc(name)) = import_model($(string(pathex)))) end

macro load_models(pathex)
    callex = :(begin end)
    for line in readlines(string(pathex))
        name, pathex = split(line, ';')
        name = isempty(name) ? gensym() : Symbol(name)
        push!(callex.args, :($(esc(name)) = import_model($(string(pathex)))))
    end

    callex
end

"""
    @import_solution "sol.jld2"
    @import_solution "sol.jld2" sol
Import a solution from a file.

# Examples
```julia
@import_solution "sir_acs_sol/serialized/sol.jld2"
```
"""
macro import_solution(pathex, varname="sol") :(JLD2.load($(string(pathex)), $(string(varname)))) end

"""
    @export_solution sol
    @export_solution sol "sol.jld2"
Export a solution to a file.

# Examples
```julia
@export_solution sol "sol.jdl2"
```
"""
macro export_solution(solex, pathex="sol.jld2") :(JLD2.save($(string(pathex)), $(string(solex)), $(esc(solex)))) end

"""
    @export_as_table sol
Export a solution as a `DataFrame`.

# Examples
```julia
@export_as_table sol
```
"""
macro export_as_table(solex, pathex="sol.jld2") :(DataFrame($(esc(solex)))) end

get_DataFrame(sol) = sol isa EnsembleSolution ? DataFrame(sol)[!, [:u, :t]] : DataFrame(sol)

"""
    @export_csv sol
    @export_csv sol "sol.csv"
Export a solution to a file.

# Examples
```julia
@export_csv sol "sol.csv"
```
"""
macro export_csv(solex, pathex="sol.csv") :(CSV.write($(string(pathex)), get_DataFrame($(esc(solex))))) end