export equalize!, @equalize

# Parse one entry of an `@equalize` equation into a `(qualifier, name)` pair: a `@catchall`/`@alias`
# macrocall yields `(:catchall, name)`/`(:alias, name)`, a bare species yields `(nothing, name)`.
# Dotted names are flattened by `underscorize`. The "_ff" suffix distinguishes it from the join
# operators' `get_eqs`.
expand_name_ff(ex) =
if ex isa Expr && isexpr(ex, :macrocall)
    (macroname(ex), underscorize(ex.args[end]))
else
    (nothing, underscorize(ex))
end

"""
Flatten one `@equalize` equation block into the list of `(qualifier, name)` pairs to be merged.
Recurses through a chained `A = B = C` (right-nested `:(=)` Exprs) so every member of the chain is
collected, tagging each via [`expand_name_ff`](@ref). Consumed by the [`@equalize`](@ref) macro.
"""
function get_eqs_ff(eq)
    return if eq isa Expr && isexpr(eq, :(=))
        [
            expand_name_ff(eq.args[1])
            isexpr(eq.args[2], :(=)) ? get_eqs_ff(eq.args[2]) : expand_name_ff(eq.args[2])
        ]
    else
        [expand_name_ff(eq)]
    end
end

"""
    equalize!(net, eqs = []) -> ReactionNetwork

Identify (collapse) sets of species in the static `net`, in place. Each block in `eqs` names species to merge — by exact `:S` index, by bare name, or by a `:catchall` name match (matching a `__`-namespaced suffix) — into a single surviving row aliased to the block's `:alias` (or its first entry). Missing attribute cells on the survivor are filled from the merged rows, every reference to a removed name is rewritten, and the removed rows are dropped by swap-and-pop; the promoted [`ReactantSpec`](@ref) table is then rebuilt so each reactant's `species` FK is repointed structurally onto the survivor (§7.4/J7, ADR 0003 Phase 2). Authoring-time only — FORBIDDEN on a live/stepping model, since it reindexes. [`@equalize`](@ref) is the declarative macro form.
"""
function equalize!(net::ReactionNetwork, eqs = [])
    specmap = Dict()
    for block in eqs
        block_alias = findfirst(e -> e[1] == :alias, block)
        block_alias = !isnothing(block_alias) ? block[block_alias][2] : first(block)[2]
        species_ixs = Int64[]
        for e in block, i in row_ids(net, :S)
            (
                (i == e[2]) ||
                    (
                    e[1] == :catchall &&
                        occursin(Regex("(__$(e[2])|$(e[2]))\$"), string(net[i, :specName]))
                ) ||
                    (e[2] == net[i, :specName])
            ) && (
                push!(species_ixs, i);
                push!(specmap, net[i, :specName] => (net[i, :specName] = block_alias))
            )
        end
        isempty(species_ixs) && continue
        species_ixs = sort(unique!(species_ixs))
        lix = first(species_ixs)
        for attr in propertynames(net.columns)
            !occursin("spec", string(attr)) && continue
            for i in species_ixs
                ismissing(net[lix, attr]) && (net[lix, attr] = net[i, attr])
            end
        end
        rem_rows!(net, :S, species_ixs[2:end])
    end

    for attr in propertynames(net.columns)
        attr == :specName && continue
        attr_ = net[:, attr]
        for i in eachindex(attr_)
            attr_[i] = escape_ref(attr_[i], collect(keys(specmap)))
            attr_[i] = recursively_substitute_vars!(specmap, attr_[i])
            net[i, attr] = attr_[i]
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
    populate_reactant_specs!(net)

    return net
end

"""
Identify (collapse) sets of species in a model — each equation names species to merge into one, so a
downstream species can be fused with an upstream one (the FK-repoint of §7.4/J7, not string surgery).

# Examples

```julia
@equalize net A = B C = D
```
"""
macro equalize(netex, exs...)
    exs = collect(exs)
    foreach(i -> (exs[i] = MacroTools.striplines(exs[i])), 1:length(exs))
    eqs = []
    foreach(ex -> ex isa Expr && merge_eqs!(eqs, get_eqs_ff(ex)), exs)

    return :(equalize!($(esc(netex)), $eqs))
end
