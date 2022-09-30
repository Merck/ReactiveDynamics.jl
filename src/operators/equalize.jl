export @equalize

expand_name_ff(ex) = ex isa Expr && isexpr(ex, :macrocall) ? (macroname(ex), underscorize(ex.args[end])) : (nothing, underscorize(ex))

"Parse species equation blocks."
function get_eqs_ff(eq)
    if eq isa Expr && isexpr(eq, :(=))
        [expand_name_ff(eq.args[1]); isexpr(eq.args[2], :(=)) ? get_eqs_ff(eq.args[2]) : expand_name_ff(eq.args[2])]
    else [expand_name_ff(eq)] end
end

function equalize!(acs::ReactionNetwork,  eqs=[])
    specmap = Dict()
    for block in eqs
        block_alias = findfirst(e -> e[1] == :alias, block)
        block_alias = !isnothing(block_alias) ? block[block_alias][2] : first(block)[2]
        species_ixs = Int64[]
        for e in block, i in 1:nparts(acs, :S)
            ((i == e[2]) || (e[1] == :catchall && occursin(Regex("(__$(e[2])|$(e[2]))\$"), string(acs[i, :specName]))) ||
                (e[2] == acs[i, :specName])) &&
                    (push!(species_ixs, i); push!(specmap, acs[i, :specName] => (acs[i, :specName] = block_alias)))
        end
        isempty(species_ixs) && continue
        species_ixs = sort(unique!(species_ixs))
        lix = first(species_ixs); for attr in propertynames(acs.subparts)
            !occursin("spec", string(attr)) && continue
            for i in species_ixs
                ismissing(acs[lix, attr]) && (acs[lix, attr] = acs[i, attr])
            end
        end
        rem_parts!(acs, :S, species_ixs[2:end])
    end

    for attr in propertynames(acs.subparts)
        attr == :specName && continue
        attr_ = getproperty(acs.subparts, attr)
        for i in 1:length(attr_)
            attr_[i] = escape_ref(attr_[i], collect(keys(specmap)))
            attr_[i] = recursively_substitute_vars!(specmap, attr_[i])
        end
    end

    acs
end

"""
Identify (collapse) a set of species in a model.

# Examples
```julia
@join acs acs1.A=acs2.A B=C
```
"""
macro equalize(acsex, exs...)
    exs = collect(exs); foreach(i -> (exs[i] = MacroTools.striplines(exs[i])), 1:length(exs))
    eqs = []; foreach(ex -> ex isa Expr && merge_eqs!(eqs, get_eqs_ff(ex)), exs)

    :(equalize!($(esc(acsex)), $eqs))
end