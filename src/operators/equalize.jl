export equalize!, @equalize

expand_name_ff(ex) =
    if ex isa Expr && isexpr(ex, :macrocall)
        (macroname(ex), underscorize(ex.args[end]))
    else
        (nothing, underscorize(ex))
    end

"""
Parse species equation blocks.
"""
function get_eqs_ff(eq)
    if eq isa Expr && isexpr(eq, :(=))
        [
            expand_name_ff(eq.args[1])
            isexpr(eq.args[2], :(=)) ? get_eqs_ff(eq.args[2]) : expand_name_ff(eq.args[2])
        ]
    else
        [expand_name_ff(eq)]
    end
end

function equalize!(acs::ReactionNetworkSchema, eqs = [])
    specmap = Dict()
    for block in eqs
        block_alias = findfirst(e -> e[1] == :alias, block)
        block_alias = !isnothing(block_alias) ? block[block_alias][2] : first(block)[2]
        species_ixs = Int64[]
        for e in block, i in parts(acs, :S)
            (
                (i == e[2]) ||
                (
                    e[1] == :catchall &&
                    occursin(Regex("(__$(e[2])|$(e[2]))\$"), string(acs[i, :specName]))
                ) ||
                (e[2] == acs[i, :specName])
            ) && (
                push!(species_ixs, i);
                push!(specmap, acs[i, :specName] => (acs[i, :specName] = block_alias))
            )
        end
        isempty(species_ixs) && continue
        species_ixs = sort(unique!(species_ixs))
        lix = first(species_ixs)
        for attr in propertynames(acs.subparts)
            !occursin("spec", string(attr)) && continue
            for i in species_ixs
                ismissing(acs[lix, attr]) && (acs[lix, attr] = acs[i, attr])
            end
        end
        rem_parts!(acs, :S, species_ixs[2:end])
    end

    for attr in propertynames(acs.subparts)
        attr == :specName && continue
        attr_ = acs[:, attr]
        for i in eachindex(attr_)
            attr_[i] = escape_ref(attr_[i], collect(keys(specmap)))
            attr_[i] = recursively_substitute_vars!(specmap, attr_[i])
            acs[i, attr] = attr_[i]
        end
    end

    # ADR 0003 Phase 2 (§7.4/J7): promote the transition↔reactant relation to the FK-exact
    # ReactantSpec table AFTER the merge. Because the merge above already collapsed the identified
    # species to a single surviving `:S` row and rewrote every reference to the survivor's name,
    # rebuilding the typed table from the post-merge `:trans` lines repoints every reactant's integer
    # `species` FK onto the survivor STRUCTURALLY — no dangling FK to a removed row, and no reactant
    # still names an eliminated alias. This is the collision-safe replacement for the string surgery
    # above at the STRUCTURAL grain (the string rewrite of `:trans` is retained only because the
    # runtime engine still parses `:trans` per tick, ADR 0003 Phase 1's behavior-preserving contract).
    populate_reactant_specs!(acs)

    return acs
end

"""
Identify (collapse) a set of species in a model.

# Examples

```julia
@join acs acs1.A = acs2.A B = C
```
"""
macro equalize(acsex, exs...)
    exs = collect(exs)
    foreach(i -> (exs[i] = MacroTools.striplines(exs[i])), 1:length(exs))
    eqs = []
    foreach(ex -> ex isa Expr && merge_eqs!(eqs, get_eqs_ff(ex)), exs)

    return :(equalize!($(esc(acsex)), $eqs))
end
