export @import_network, @export_network
export @load_models
export @import_solution, @export_solution
export @export_solution_as_table, @export_solution_as_csv
export @export, @import

using TOML, JLD2, CSV
using DataFrames

const objects_aliases = Dict(:S => "spec", :T => "trans", :P => "prm", :M => "meta",
                             :E => "event", :obs => "obs")

const RN_attrs = string.(propertynames(ReactionNetwork().subparts))

function get_attrs(object)
    object = object isa Symbol ? objects_aliases[object] : object

    filter(x -> occursin(object, x), RN_attrs)
end

function export_network(acs::ReactionNetwork)
    dict = Dict()
    for (key, val) in objects_aliases
        push!(dict, val => [])
        for i in 1:nparts(acs, key)
            dict_ = Dict()
            for attr in get_attrs(val)
                attr_val = acs[i, Symbol(attr)]
                ismissing(attr_val) && continue
                attr_val = attr_val isa Number ? attr_val : string(attr_val)
                push!(dict_, string(attr) => attr_val)
            end
            push!(dict[val], dict_)
        end
    end

    dict
end

function load_network(dict::Dict)
    acs = ReactionNetwork()
    for (key, val) in objects_aliases
        val == "prm" && continue
        for row in get(dict, val, [])
            i = add_part!(acs, key)
            for (attr, attrval) in row
                set_subpart!(acs, i, Symbol(attr), attrval)
            end
        end
    end

    for row in get(dict, "prm", [])
        i = add_part!(acs, :P)
        for (attr, attrval) in row
            if attr == "prmVal"
                attrval = attrval isa String ? eval(Meta.parseall(attrval)) : attrval
            end
            set_subpart!(acs, i, Symbol(attr), attrval)
        end
    end

    for row in get(dict, "registered", [])
        eval(Meta.parseall(row["body"]))
    end

    assign_defaults!(acs)
end

function import_network_csv(pathmap)
    dict = Dict()
    for (key, paths) in pathmap
        push!(dict, key => [])
        for path in paths
            data = DataFrame(CSV.File(path; delim = ";;", types = String,
                                      stripwhitespace = true, comment = "#"))
            for row in eachrow(data)
                object = Dict()
                for (attr, val) in Iterators.zip(keys(row), values(row))
                    !ismissing(val) && push!(object, string(attr) => val)
                end
                push!(dict[key], object)
            end
        end
    end

    load_network(dict)
end

function import_network(path::AbstractString)
    if splitext(path)[2] == ".csv"
        pathmap = Dict(val => []
                       for val in [collect(values(objects_aliases)); "registered"])
        for row in CSV.File(path; delim = ";;", stripwhitespace = true, comment = "#")
            push!(pathmap[row.type], joinpath(dirname(path), row.path))
        end

        import_network_csv(pathmap)
    else
        load_network(TOML.parsefile(path))
    end
end

function export_network(acs::ReactionNetwork, path::AbstractString)
    if splitext(path)[2] == ".csv"
        exported_network = export_network(acs)
        paths = DataFrame(type = [], path = [])
        for (key, objs) in exported_network
            push!(paths, (key, "export-$key.csv"))
            objs_exported = DataFrame(Dict(attr => [] for attr in get_attrs(key)))
            for obj in objs
                push!(objs_exported,
                      [get(obj, key, missing) for key in names(objs_exported)])
            end

            CSV.write(joinpath(dirname(path), "export-$key.csv"), objs_exported,
                      delim = ";;")
        end
        CSV.write(path, paths, delim = ";;")
    else
        open(io -> TOML.print(io, export_network(acs)), path, "w+")
    end
end

"""
Export model to a file: this can be either a single TOML file encoding the entire model,
or a batch of CSV files (a root file and a number of files, each per a class of objects).

See `tutorials/loadsave` for an example.

# Examples
```julia
@export_network acs "acs_data.toml" # as a TOML
@export_network acs "csv/model.csv" # as a CSV
```
"""
macro export_network(acsex, pathex)
    :(export_network($(esc(acsex)), $(string(pathex))))
end

"""
Import a model from a file: this can be either a single TOML file encoding the entire model,
or a batch of CSV files (a root file and a number of files, each per a class of objects).

See `tutorials/loadsave` for an example.

# Examples
```julia
@import_network "model.toml"
@import_network "csv/model.toml"
```
"""
macro import_network(pathex, name = gensym())
    :($(esc(name)) = import_network($(string(pathex))))
end

macro load_models(pathex)
    callex = :(begin end)
    for line in readlines(string(pathex))
        name, pathex = split(line, ';')
        name = isempty(name) ? gensym() : Symbol(name)
        push!(callex.args, :($(esc(name)) = import_network($(string(pathex)))))
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
macro import_solution(pathex, varname = "sol")
    :(JLD2.load($(string(pathex)), $(string(varname))))
end

"""
    @export_solution sol
    @export_solution sol "sol.jld2"
Export a solution to a file.

# Examples
```julia
@export_solution sol "sol.jdl2"
```
"""
macro export_solution(solex, pathex = "sol.jld2")
    :(JLD2.save($(string(pathex)), $(string(solex)), $(esc(solex))))
end

"""
    @export_solution_as_table sol
Export a solution as a `DataFrame`.

# Examples
```julia
@export_solution_as_table sol
```
"""
macro export_solution_as_table(solex, pathex = "sol.jld2")
    :(DataFrame($(esc(solex))))
end

get_DataFrame(sol) = sol isa EnsembleSolution ? DataFrame(sol)[!, [:u, :t]] : DataFrame(sol)

"""
    @export_solution_as_csv sol
    @export_solution_as_csv sol "sol.csv"
Export a solution to a file.

# Examples
```julia
@export_solution_as_csv sol "sol.csv"
```
"""
macro export_solution_as_csv(solex, pathex = "sol.csv")
    :(CSV.write($(string(pathex)), get_DataFrame($(esc(solex)))))
end
