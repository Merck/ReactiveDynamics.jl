# WS-B2: `@agentize` is now a thin authoring sugar over the `ReactionNetworkProblem` constructor
# (ADR 0001 offered "implement as thin sugar OR delete the export"; WS-4 deleted it, WS-B2 re-adds
# the sugar). It reimplements NO construction logic — it expands to exactly one
# `ReactionNetworkProblem(net[, u0, p]; …)` call, the same public contract as
# `ReactionNetworkProblem(net) |> simulate` (ADR 0012). See the macro at the foot of this file.
#
# EXPORT SITE (deviation flagged): the export is co-located here at solve.jl:1 — the exact line the
# WS-4 dangling `export @agentize` was deleted from, and matching the file-local export convention
# (`create.jl` exports `@reaction_network`; `update.jl` exports its role macros next to their
# definitions). The WS-B2 brief said "re-add in `src/ReactiveDynamics.jl`", but that module root
# centralizes only the ADR-0003 store-shim exports; every interface macro is exported from its own
# file. Co-location keeps the export beside the definition.
export @agentize
#
# NOTE (ADR 0014 §A1 / CONTRACT §15.1): this file previously carried ~195 lines of dead SciML
# plotting — `plot_summary`/`plot_ensemble_sol`/`first_sol`/`plot_from_log` and the `@plot` macro —
# all written against SciML `EnsembleSummary`/`EnsembleSolution` objects that ADR 0001 removed when
# it demoted SciML. Those types were never imported, so the macro was dead-on-arrival on
# `ref-agents`. The result-plotting story is now the model-agnostic `@recipe` set in the `RDPlotsExt`
# package extension (ext/RDPlotsExt.jl); the live generic `_draw` reduction moved there too (it
# needs `Plots`, which is now a weakdep). Nothing here imports `Plots` anymore.

"""
    @agentize net [u0] [p] [seed=…] [tspan=…] [name=…] …

Thin authoring sugar over the [`ReactionNetworkProblem`](@ref) constructor: expand to exactly one
`ReactionNetworkProblem(net[, u0, p]; kwargs...)` call (ADR 0001 / ADR 0012). Agentization is already
implicit in that constructor — it builds the `@aagent` and entangles the `"structured"` container —
so this macro adds NO second construction path; it only lowers to the public constructor call.

The one ergonomic win is AUTO-NAMING: when `net` is a plain binding (a `Symbol`), the agent `name`
defaults to that binding's name (`@agentize net` ⇒ `name = "net"`). An explicit `name=` kwarg always
wins, and a non-symbol `net` expression (`@agentize build_net()`) falls back to the constructor's own
`name` default — the auto-name is a compile-time `String` literal, so it never fights macro hygiene.

# Examples

```julia
prob = @agentize net                          # name = "net"
prob = @agentize net u0 p seed=1 tspan=10     # positional u0/p + forwarded kwargs
prob = @agentize net name="custom"            # explicit name overrides the auto-name
prob = @agentize build_net()                  # non-symbol acs → constructor's default name
```
"""
macro agentize(netex, args...)
    # `args_kwargs` esc's the POSITIONAL args (correct — u0/p are caller expressions) but leaves the
    # kwarg VALUES unescaped in `Expr(:kw, key, value)`. Forwarded verbatim into the emitted call,
    # such a value would resolve in RD's module scope, so a caller-local (`@agentize net seed = s`)
    # would UndefVarError. Esc each kwarg value here — NOT in the shared helper (@problematize et al.
    # rely on its current behavior) — remapping `Expr(:kw, k, v)` → `Expr(:kw, k, esc(v))`.
    args, kwargs = args_kwargs(args)
    kwargs = [Expr(:kw, kw.args[1], esc(kw.args[2])) for kw in kwargs]
    # Auto-name only for a bare binding, and only when the author did not pass an explicit `name=`.
    # `netex isa Symbol` ⇒ inject a `String` literal, so the name is fixed at expansion time and needs
    # no `esc` (leave it bare — it is a constant, not a caller binding); anything else defers to the
    # constructor's `name` default. Appended AFTER the esc-remap so it stays an unescaped literal.
    if netex isa Symbol && isnothing(findfirst(ex -> ex.args[1] == :name, kwargs))
        push!(kwargs, Expr(:kw, :name, string(netex)))
    end

    quote
        ReactionNetworkProblem($(esc(netex)), $(args...); $(kwargs...))
    end
end
