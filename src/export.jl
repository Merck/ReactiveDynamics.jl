# Results export bundle (ADR 0013 §C / CONTRACT §14.3).
#
# Writes a finished run (or ensemble) to a directory bundle implementing the §8.5 outputs clause,
# format PER ARTIFACT (one format does not fit all — they have different shapes):
#
#   | artifact                              | shape                         | format             |
#   | prob.sol trajectory                   | rectangular, all-Float64      | CSV (+ Arrow)      |
#   | program_ledger(prob)                  | rectangular DataFrame         | CSV (+ Arrow)      |
#   | token_trajectory(prob)                | semi-structured (kind fields) | CSV (+ Arrow)      |
#   | ensemble summary                      | rectangular                   | CSV (+ Arrow)      |
#   | prob.log event stream                 | heterogeneous tuple-per-tick  | JSON               |
#   | tokens.json (per-token histories)     | semi-structured               | JSON               |
#   | run.json manifest (hash/seed/tspan)   | nested                        | JSON               |
#
# Dependency tiering preserves the ADR 0005 minimalism: JSON is always-available (reuses the
# JSON.jl dep `serialize.jl` already uses), CSV via the present CSV.jl dep, and Arrow lives behind
# the `RDArrowExt` weakdep — `_arrow_write` is a no-op here and is overridden by the extension when
# `Arrow` is loaded, so an Arrow-less run silently emits only the CSV/JSON core. Layout is the §8.5
# `runs/<model_hash>/<seed>/{trajectory.{csv,arrow}, ledger.{csv,arrow}, tokens.{csv,arrow},
# events.json, tokens.json, run.json}`, with `ensemble.json` at the ensemble root.

using CSV

export export_run, export_ensemble

# Arrow seam (ADR 0013 §C). Declared as an EMPTY generic (no method) so RDArrowExt can supply the
# SOLE method without the "method overwriting during precompilation" error a core default body would
# cause. `_arrow_available()` reports whether the extension has loaded a method, so `export_run`
# emits the `.arrow` sibling exactly when the user has `Arrow` and otherwise silently writes only the
# CSV/JSON core.
function _arrow_write end
_arrow_available() = !isempty(methods(_arrow_write))

# Write a rectangular DataFrame as the CSV core + the optional faithful Arrow sibling. `stem` is the
# path WITHOUT extension; emits `stem.csv` always and `stem.arrow` when RDArrowExt is loaded.
function _write_rectangular(stem::AbstractString, df::DataFrame)
    CSV.write(stem * ".csv", df)
    _arrow_available() && _arrow_write(stem * ".arrow", df)
    return stem
end

# The model hash that names the bundle directory (the same `hash(net)` the checkpoint manifest uses,
# interface/checkpoint.jl) — pins the bundle to a replayable (model, seed) pair (§4 D6, Invariant 5).
_model_hash(prob::ReactionNetworkProblem) = hash(prob.network)

# Encode the heterogeneous `prob.log` event stream as JSON-friendly records. Each row is a tuple
# `(tag::Symbol, t, payload…)`; we emit `{"event": tag, "t": t, "data": [payload…]}` with the
# payload rendered structurally (Dicts/vectors pass through JSON.jl; everything else stringifies so
# the stream always round-trips structurally — the §8 guarantee for a non-rectangular artifact).
function _log_to_records(log)
    recs = Vector{Dict{String, Any}}(undef, length(log))
    for (i, row) in enumerate(log)
        tag = row[1]
        t = length(row) >= 2 ? row[2] : nothing
        payload = length(row) >= 3 ? collect(row[3:end]) : Any[]
        recs[i] = Dict{String, Any}(
            "event" => string(tag),
            "t" => t,
            "data" => map(_jsonable, payload),
        )
    end
    return recs
end

# Render an arbitrary log payload element into something JSON.jl emits faithfully. Dicts and vectors
# recurse; numbers/strings/bools pass through; anything else (a Symbol, a NamedTuple) stringifies.
_jsonable(x::Union{Real, AbstractString, Bool, Nothing}) = x
_jsonable(x::Symbol) = string(x)
_jsonable(x::AbstractDict) = Dict{String, Any}(string(k) => _jsonable(v) for (k, v) in x)
_jsonable(x::AbstractVector) = map(_jsonable, x)
_jsonable(x::NamedTuple) = Dict{String, Any}(string(k) => _jsonable(getfield(x, k)) for k in keys(x))
_jsonable(x::Tuple) = map(_jsonable, collect(x))
_jsonable(x) = string(x)

# Per-token histories as nested JSON: token name → ordered list of {t, species, fields…} records,
# from the §14.1 trajectory store. Round-trips structurally like the model JSON (Invariant 5).
function _tokens_to_records(prob::ReactionNetworkProblem)
    out = Dict{String, Vector{Dict{String, Any}}}()
    for (t, name, species, fields) in prob.token_trajectory
        rec = Dict{String, Any}("t" => t, "species" => string(species))
        for k in keys(fields)
            rec[string(k)] = _jsonable(getfield(fields, k))
        end
        push!(get!(out, name, Dict{String, Any}[]), rec)
    end
    return out
end

"""
    export_run(prob, dir; with_model = true) -> dir

Write the finished run `prob` to `dir` as the §14.3 bundle: `trajectory.{csv,arrow}` (`prob.sol`),
`ledger.{csv,arrow}` (`program_ledger`), `tokens.{csv,arrow}` (the §14.1 `token_trajectory` long
form, when any token opted in), `events.json` (the `prob.log` stream), `tokens.json` (per-token
histories), and `run.json` (the manifest: model hash, seed, tspan, dt, schema version, and — when
`with_model` — the embedded JSON model so the bundle is self-describing/replayable). Arrow siblings
appear only when `RDArrowExt` is loaded (the user has `Arrow`); otherwise the CSV/JSON core is
written alone. Returns `dir`.
"""
function export_run(prob::ReactionNetworkProblem, dir::AbstractString; with_model::Bool = true)
    mkpath(dir)

    # Track exactly what is written so the manifest index reflects the REAL bundle contents (the
    # `.arrow` siblings appear only when RDArrowExt is loaded; `tokens.*` only when a kind opted in).
    artifacts = String[]
    arrow = _arrow_available()
    _emit_rect(stem) = (push!(artifacts, stem * ".csv"); arrow && push!(artifacts, stem * ".arrow"))

    # Rectangular artifacts → CSV core (+ Arrow sibling via the weakdep).
    _write_rectangular(joinpath(dir, "trajectory"), prob.sol); _emit_rect("trajectory")
    _write_rectangular(joinpath(dir, "ledger"), program_ledger(prob)); _emit_rect("ledger")
    traj = token_trajectory(prob)
    if nrow(traj) > 0
        _write_rectangular(joinpath(dir, "tokens"), traj); _emit_rect("tokens")
    end

    # Heterogeneous artifacts → JSON.
    open(joinpath(dir, "events.json"), "w") do io
        JSON.print(io, _log_to_records(prob.log))
    end
    open(joinpath(dir, "tokens.json"), "w") do io
        JSON.print(io, _tokens_to_records(prob))
    end
    append!(artifacts, ["events.json", "tokens.json"])

    # The manifest (nested) — pins (model_hash, seed) so the bundle is traceable to a replayable run.
    manifest = Dict{String, Any}(
        "schema" => "rd-run/1",
        "model_hash" => string(_model_hash(prob)),
        "seed" => prob.seed === nothing ? nothing : string(prob.seed),
        "tspan" => collect(prob.tspan),
        "dt" => prob.dt,
        "species" => string.(prob.network[:, :specName]),
        "artifacts" => artifacts,
    )
    with_model && (manifest["model"] = JSON.parse(to_json_model(prob)))
    open(joinpath(dir, "run.json"), "w") do io
        JSON.print(io, manifest)
    end
    return dir
end

"""
    export_ensemble(ens, dir; metric = nothing, with_model = true) -> dir

Write the ensemble `ens` to `dir` as the §14.3 layout: one subdirectory `member_<k>/` per member
(each an `export_run` bundle), plus a top-level `ensemble.json` recording the per-member seeds, the
root seed, the run mode, and — when a `metric::(member -> Real)` is supplied — the `summarize` table
over it. Returns `dir`.
"""
function export_ensemble(
        ens::EnsembleProblem, dir::AbstractString;
        metric = nothing, with_model::Bool = true
    )
    mkpath(dir)
    for (k, m) in enumerate(ens.members)
        export_run(m, joinpath(dir, "member_$k"); with_model = with_model)
    end
    top = Dict{String, Any}(
        "schema" => "rd-ensemble/1",
        "root_seed" => ens.root_seed,
        "mode" => string(ens.mode),
        "nseed" => length(ens.members),
        "seeds" => string.(ens.seeds),
    )
    if metric !== nothing
        s = summarize(ens, metric)
        top["summary"] = Dict{String, Any}(string(k) => getfield(s, k) for k in keys(s))
    end
    open(joinpath(dir, "ensemble.json"), "w") do io
        JSON.print(io, top)
    end
    return dir
end
