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
# AlgebraicAgents` / `GeneratedExpressions` (state.jl:1, ReactiveDynamics.jl:9). Some of those are
# themselves dangling *upstream* (e.g. AA's `@derived`, `@integration`, `AgentCall`) — not RD's bug
# to fix and not RD's export to police. So we test only the symbols RD itself defines-or-exports,
# identified as: the symbol is exported by RD AND (it is defined in a source file owned by the
# ReactiveDynamics module, i.e. `parentmodule` of the binding is `ReactiveDynamics`, OR it is NOT
# exported by either reexported upstream module). This isolates the invariant WS-4 established.

using ReactiveDynamics, Test

RD = ReactiveDynamics

# Upstream reexport surfaces we do NOT police (their danglers are upstream's concern). We reach
# them THROUGH RD (`@reexport using AlgebraicAgents` at state.jl:1; `@reexport using
# GeneratedExpressions` at ReactiveDynamics.jl:9) rather than `import`ing them directly, so the
# test needs no extra test/Project.toml deps — both are RD-visible submodules.
const _REEXPORT_MODULES =
    filter(m -> m isa Module, (getproperty(RD, :AlgebraicAgents), getproperty(RD, :GeneratedExpressions)))

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
end
