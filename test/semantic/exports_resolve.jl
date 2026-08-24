# WS-4 housekeeping test — every RD-owned exported symbol must resolve.
#
# The module root `src/ReactiveDynamics.jl` declares no `export` of its own; every public symbol
# comes from an included file's `export` line (see INVENTORY.md → Public-API Audit). Historically a
# handful of those exports were DANGLING — exported names with no definition anywhere (`@agentize`,
# `@prob_role`, `@list_by_role`, `@list_roles`, plus a `@prob_check_verbose` whose body called an
# undefined `check_params`). Those threw only at call time (or, for the macros, resolved to nothing),
# so `using ReactiveDynamics` stayed green and the rot went unnoticed. This testset is the guard: it
# iterates `names(ReactiveDynamics)` and asserts each RD-OWNED export is `isdefined` and resolvable.
#
# Scope note: `names(ReactiveDynamics)` also surfaces symbols pulled in by `@reexport using
# AlgebraicAgents` (state.jl:1). Some of those are themselves dangling *upstream* (e.g. AA's
# `@derived`, `@integration`, `AgentCall`) — not RD's bug to fix and not RD's export to police. So we
# test only the symbols RD itself defines-or-exports, identified as: the symbol is exported by RD AND
# (it is defined in a source file owned by the ReactiveDynamics module, i.e. `parentmodule` of the
# binding is `ReactiveDynamics`, OR it is NOT exported by the reexported upstream module). This
# isolates the invariant WS-4 established.

using ReactiveDynamics, Test

RD = ReactiveDynamics

# Upstream reexport surface we do NOT police (its danglers are upstream's concern). We reach it
# THROUGH RD (`@reexport using AlgebraicAgents` at state.jl:1) rather than `import`ing it directly, so
# the test needs no extra test/Project.toml deps — it is an RD-visible submodule. (ADR 0015 dropped
# the GeneratedExpressions reexport, so AlgebraicAgents is now the sole upstream surface here.)
const _REEXPORT_MODULES =
    filter(m -> m isa Module, (getproperty(RD, :AlgebraicAgents),))

"Is `sym` merely reexported into RD from an upstream package (not authored by RD)?"
function _is_upstream_reexport(sym)
    return any(m -> sym in names(m), _REEXPORT_MODULES)
end

"Resolve the binding behind an exported name (handles the `@macro` → `Symbol(\"@macro\")` case)."
_resolves(sym) = isdefined(RD, sym)

@testset "Every RD-owned export resolves" begin
    exported = names(RD)                       # includes the module's own name + reexports
    # Drop the module's self-reference and the upstream reexports; keep RD-authored public symbols.
    rd_owned = filter(exported) do sym
        sym === :ReactiveDynamics && return false
        _is_upstream_reexport(sym) && return false
        return true
    end

    @test !isempty(rd_owned)                   # sanity: we actually found RD's public API

    # Report-friendly: collect any danglers before asserting, so a failure names them all.
    dangling = filter(sym -> !_resolves(sym), rd_owned)
    @test isempty(dangling)
    isempty(dangling) || @error "Dangling RD exports (exported but undefined)" dangling

    # And assert positively on each so the count contributes and failures pinpoint the symbol.
    for sym in rd_owned
        @test _resolves(sym)
    end

    # Macros specifically: an exported `@foo` shows up in `names` as `Symbol("@foo")`; assert the
    # macro binding is a genuine macro (a `getproperty` on the module must not throw).
    for sym in rd_owned
        startswith(string(sym), "@") || continue
        @test getproperty(RD, sym) isa Function   # macros are `Function`s (methods on `var"@foo"`)
    end

    # WS-B2 regression pin: `@agentize` was the original WS-4 dangling export (exported, never
    # defined → this very guard would have failed). Assert it is BOTH exported and resolvable now
    # that it ships as thin constructor sugar (src/interface/solve.jl). This is a named guard so a
    # future re-deletion of the macro without the export fails here, not silently at call time.
    @test Symbol("@agentize") in rd_owned
    @test _resolves(Symbol("@agentize"))
    @test getproperty(RD, Symbol("@agentize")) isa Function
end

# ── ADR 0017 vocabulary gate: the retired chemistry spellings live only where documented ────────
#
# The rename retired `species`, `reactant`, the `spec<Attr>` store columns and `stoich`. Each still
# resolves for ONE release, but only from a small, documented set of files: the Tier-1/Tier-5
# deprecation shims, the `@aka` legacy object name, the `:species` selector field, and the
# serializer's wire-key read aliases. This testset is the ADR's step-6 grep made executable — a
# retired spelling reappearing anywhere else in `src/` fails here, rather than quietly rebuilding
# the two-vocabulary problem the rename exists to remove.
@testset "ADR 0017: retired vocabulary confined to the documented legacy sites" begin
    retired = r"\bspecies\b|reactant|spec[A-Z]|\bstoich"
    # Each entry carries a one-release shim that a forwarding binding cannot express.
    legacy_files = Set(
        [
            "src/ReactiveDynamics.jl",   # Tier 1/5 @deprecate shims + the @add_species macro shim
            "src/interface/update.jl",   # @aka net species = resource
            "src/predicates.jl",         # @select field :species → :place
            "src/solvers.jl",            # @advance field :species → :place
            "src/serialize.jl",          # wire-key aliases: places/arcs/place/multiplicity
        ]
    )
    root = pkgdir(RD)
    srcfiles = String[]
    for (dir, _, files) in walkdir(joinpath(root, "src")), fn in files
        endswith(fn, ".jl") && push!(srcfiles, relpath(joinpath(dir, fn), root))
    end
    @test !isempty(srcfiles)                      # sanity: we actually walked the source tree

    offenders =
        filter(f -> f ∉ legacy_files && occursin(retired, read(joinpath(root, f), String)), srcfiles)
    @test isempty(offenders)
    isempty(offenders) ||
        @error "retired ADR-0017 vocabulary outside the documented legacy sites" offenders

    # And the allowlist is not stale: every file named above really does still carry legacy handling,
    # so removing a shim forces this list to shrink with it.
    for f in legacy_files
        @test occursin(retired, read(joinpath(root, f), String))
    end
end

# ── ADR 0018 gate: the swapped firing/transition names do not come back ──────────────────────────
#
# ADR 0018 unswapped the inversion ADR 0017 left standing — the in-flight instance type is `Firing`,
# the static rows are `state.transitions`, the per-tick snapshot is `state.sampled_transitions`, and a
# token's back-pointer is `bound_firing`. The retired spellings have no shim except the type alias, so
# the only thing that keeps them from creeping back into new code is this grep.
@testset "ADR 0018: the pre-rename firing vocabulary is gone from src/" begin
    retired = r"\btransition_recipes\b|\bongoing_transitions\b|\bbound_transition\b"
    root = pkgdir(RD)
    # `src/ReactiveDynamics.jl` names all three in the comment above the type alias, documenting what
    # was renamed and why the fields get no shim — the one place the old spellings may appear.
    legacy_files = Set(["src/ReactiveDynamics.jl"])
    srcfiles = String[]
    for (dir, _, files) in walkdir(joinpath(root, "src")), fn in files
        endswith(fn, ".jl") && push!(srcfiles, relpath(joinpath(dir, fn), root))
    end
    offenders =
        filter(f -> f ∉ legacy_files && occursin(retired, read(joinpath(root, f), String)), srcfiles)
    @test isempty(offenders)
    isempty(offenders) || @error "pre-ADR-0018 firing vocabulary is back in src/" offenders

    # The two live tables are distinct fields, not one aliased dict (the collision that made the
    # first sweep attempt fail to load), and both carry the `trans*` column family.
    @test :transitions ∈ fieldnames(RD.ReactionNetworkProblem)
    @test :sampled_transitions ∈ fieldnames(RD.ReactionNetworkProblem)
    @test :ongoing_firings ∈ fieldnames(RD.ReactionNetworkProblem)
    @test :bound_firing ∈ fieldnames(RD.BaseStructuredToken)

    # The type rename keeps a one-release forwarding alias, since every structured-token example
    # writes `ReactiveDynamics.Transition` by hand in a `past_bonds` element type.
    @test RD.Transition === RD.Firing
end
