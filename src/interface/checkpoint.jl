# State dump / restore (ADR 0007 §C, CONTRACT §10.5) — serialize a live run at a TICK BOUNDARY
# into an eval-free, JSON-representable artifact, and reconstruct it. The dump schema is the
# declarative initial marking (§B) PLUS the dynamic run state (clock, RNG, creation counters,
# plain u, the full token population with CURRENT field values, rule latches). A zero-tick dump
# with no in-flight transitions IS an initial marking — the "initial marking = checkpoint"
# identity (§C). Eval-free: tokens carry field VALUES + a kind NAME resolved through the host
# registry on restore, never Julia source (§8.4 S4 / ADR 0006 §B).
#
# SCOPE (Milestone-1): the common tick-boundary case — t, rng, counters, u, the token population,
# and `once`-rule latches. Mid-run in-flight `ongoing` Transition snapshots (the §C open question —
# frozen sampled-attr dicts + bound-token relink by uuid) are NOT serialized; `dump_state` requires
# an empty `ongoing` set (a clean tick boundary) and says so if not. This covers halt/resume and the
# zero-tick "dump == initial marking" identity; the heavier mid-cycle resume is deferred.

export StateDump, dump_state, restore

"""
    StateDump

An eval-free, JSON-representable snapshot of a live run at a TICK BOUNDARY (ADR 0007 §C / CONTRACT §10.5), produced by [`dump_state`](@ref) and consumed by [`restore`](@ref). It holds the declarative initial marking (§B) PLUS the dynamic run state: the model hash, the clock (`t`, `tspan`, `dt`), the RNG state (Xoshiro `s0..s4`), the creation counters and `name → creation_index` map, the plain-place column vector `u`, the full token population as `(place, name, creation_index, fields, bound)` tuples with CURRENT field values, and the `once`-rule enabled latches. Eval-free by construction: tokens carry field VALUES + a kind NAME resolved through the host registry on restore, never Julia source (§8.4 S4 / ADR 0006 §B). A zero-tick dump with nothing in flight IS an initial marking — the "initial marking = checkpoint" identity (§C).
"""
struct StateDump
    model_hash::UInt64
    t::Float64
    tspan::Tuple{Float64, Float64}    # time control, so restore is self-contained
    dt::Float64
    rng_state::NTuple{5, UInt64}      # Xoshiro s0..s4 (eval-free, JSON-representable)
    creation_counters::Dict{Symbol, Int}
    creation_index::Dict{String, Int}
    u::Vector{Float64}               # plain-place columns; structured columns are re-derived
    tokens::Vector{NamedTuple}       # (place, name, creation_index, fields::Dict, bound::Bool)
    rule_latches::Dict{Symbol, Bool}  # `once`-rule enabled state
end

# The host-struct fields of a token beyond the @aagent/protocol injected ones — the modeling
# attributes (phase, npv, …) we snapshot as current literal values.
const _PROTOCOL_FIELDS = (
    :uuid, :name, :parent, :inners, :relpathrefs, :opera,
    :place, :bound_transition, :past_bonds,
)
_token_attr_fields(tok) = filter(f -> !(f in _PROTOCOL_FIELDS), fieldnames(typeof(tok)))

"""
    dump_state(problem) -> StateDump

Serialize a live `problem` into an eval-free [`StateDump`](@ref) for halt/resume or the zero-tick "dump == initial marking" identity (ADR 0007 §C / CONTRACT §10.5). Captures the clock, RNG state, creation counters, plain-place `u`, the token population with each token's CURRENT field values, and the `once`-rule latches; pair with [`restore`](@ref) to reconstruct the run.

DELIBERATE DEFERRAL (Milestone-1): `dump_state` requires a CLEAN TICK BOUNDARY — an empty `ongoing` transition set — and `error`s otherwise. Mid-cycle in-flight `Transition` instances (their frozen sampled-attr dicts and bound-token relink-by-uuid — the §C open question) are NOT serialized; the heavier mid-cycle resume is deferred. Step to a boundary where no instance is mid-cycle (or `reinit!`) before dumping.
"""
function dump_state(problem::ReactionNetworkProblem)
    isempty(problem.ongoing_transitions) || error(
        "dump_state: $(length(problem.ongoing_transitions)) in-flight transition(s) — dump is " *
            "only supported at a clean tick boundary (empty `ongoing`) in Milestone-1 (ADR 0007 §C " *
            "open question). Step to a boundary where no instance is mid-cycle, or use reinit! to reset.",
    )
    toks = NamedTuple[]
    for tok in values(inners(getagent(problem, "structured")))
        fields = Dict{Symbol, Any}(f => getproperty(tok, f) for f in _token_attr_fields(tok))
        push!(
            toks,
            (
                place = get_place(tok),
                name = AlgebraicAgents.getname(tok),
                creation_index = get(problem.creation_index, AlgebraicAgents.getname(tok), 0),
                fields = fields,
                bound = isblocked(tok),
            ),
        )
    end
    rng = problem.rng
    return StateDump(
        hash(problem.network),
        problem.t,
        problem.tspan,
        problem.dt,
        (rng.s0, rng.s1, rng.s2, rng.s3, rng.s4),
        copy(problem.creation_counters),
        copy(problem.creation_index),
        copy(problem.u),
        toks,
        Dict{Symbol, Bool}(r.id => r.enabled for r in problem.rules if r.fire_mode === :once),
    )
end

"""
    restore(spec, dump::StateDump; registry = Dict{Symbol, Any}(), kwargs...) -> ReactionNetworkProblem

Reconstruct a live run from a [`StateDump`](@ref) — construction with overlays (ADR 0007 §C). Builds a fresh `ReactionNetworkProblem` from `spec`, then overlays the dumped clock/RNG/creation-counters/`u` and rebuilds the token population in creation-index order — each token's kind resolved via the `registry` to a host constructor, its fields set to the dumped literal values (eval-free; no Julia source is carried). Structured `u` columns are re-derived from the restored population, not copied. Warns on a model-hash mismatch (restoring against a different `spec` is ill-defined). `kwargs` are forwarded to the `ReactionNetworkProblem` constructor.
"""
function restore(spec, dump::StateDump; registry = Dict{Symbol, Any}(), kwargs...)
    hash(spec) == dump.model_hash || @warn "restore: model hash mismatch — restoring a dump " *
        "against a different spec is ill-defined (ADR 0007 §C open question)."
    problem = ReactionNetworkProblem(
        spec;
        registry = registry,
        population = [],
        tspan = dump.tspan,
        dt = dump.dt,
        kwargs...,
    )
    # rebuild the token population from the dump, in creation-index order, via the registry
    container = getagent(problem, "structured")
    for tok in collect(values(inners(container)))
        disentangle!(tok)
    end
    empty!(problem.creation_counters)
    empty!(problem.creation_index)
    for td in sort(dump.tokens; by = t -> t.creation_index)
        haskey(registry, td.place) ||
            error("restore: no registry constructor for kind $(td.place)")
        tok = registry[td.place](problem, td.fields)
        for (f, v) in td.fields
            hasproperty(tok, f) && setproperty!(tok, f, v)
        end
        add_structured_token!(problem, tok)
    end
    # overlay the dynamic run state
    problem.t = dump.t
    problem.u .= dump.u
    problem.rng = Random.Xoshiro(dump.rng_state...)
    # Restore the creation counters AND the realized (name → creation_index) map from the dump,
    # rather than relying on the rebuild order to reproduce them — the dump is the source of truth
    # for the (place, creation_index) selection order (defensive against future rebuild changes).
    merge!(empty!(problem.creation_counters), dump.creation_counters)
    merge!(empty!(problem.creation_index), dump.creation_index)
    for r in problem.rules
        haskey(dump.rule_latches, r.id) && (r.enabled = dump.rule_latches[r.id])
    end
    update_u_structured!(problem)   # re-derive structured columns from the restored population
    return problem
end
