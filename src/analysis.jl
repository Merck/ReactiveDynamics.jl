# Analysis & observability (ADR 0013 / CONTRACT §14).
#
# Three additive, READ-ONLY pieces over a finished (or in-progress) run — none touches the §1–§9
# dynamics:
#   §14.1  the per-token trajectory log    — `push_token_trajectory_row!` (the tick step, called
#          from solvers.jl `_step!`), the `token_trajectory` read API, and the "typical" helpers
#          `representative_token` (medoid) / `trajectory_envelope` (median + IQR band).
#   §14.2  the ensemble runner             — `ensemble`/`summarize`/`treatment_effect`, with
#          `EnsembleProblem <: AbstractAlgebraicAgent` implementing the §13.2 read surface so an
#          ensemble is itself a readable/drawable AA hierarchy node (promotes demo/.../analysis.jl).
#
# The `log_token_fields` per-kind hook (default empty) lives host-side in interface/agents.jl, beside
# the token-kind structs. The store field `state.token_trajectory` is defined in state.jl.

using Statistics

export token_trajectory, representative_token, trajectory_envelope
export ensemble, summarize, treatment_effect, EnsembleProblem

# ════════════════════════════════════════════════════════════════════════════════════════
# §14.1 — the per-token trajectory log
# ════════════════════════════════════════════════════════════════════════════════════════

"""
    push_token_trajectory_row!(state)

Append this tick's per-token snapshot to `state.token_trajectory` (ADR 0013 §A / CONTRACT §14.1).
For every structured token whose kind opts in (`log_token_fields(tok)` returns a non-empty
NamedTuple), record `(t, token_name, place, fields)`. Tokens are iterated in `token_sortkey`
order (place, creation_index, uuid) — the SAME deterministic order, and the SAME tick seam (right
after `push_program_ledger_row!`), as the per-program ledger row this generalizes, so the log is
byte-reproducible under `(model, seed)` (§4 D4). The hook is a pure 𝓕ₜ-measurable field read — no
RNG, no future (Invariant 1). A kind that doesn't opt in contributes no rows (Invariant 2).
"""
function push_token_trajectory_row!(state::ReactionNetworkProblem)
    container = getagent(state, "structured")
    toks = collect(values(inners(container)))
    sort!(toks; by = a -> token_sortkey(state, a))
    for tok in toks
        fields = log_token_fields(tok)
        isempty(fields) && continue
        push!(
            state.token_trajectory,
            (state.t, String(AlgebraicAgents.getname(tok)), get_place(tok), fields),
        )
    end
    return state
end

# The union of field names appearing in the stored rows, in first-seen order — the columns the long
# DataFrame carries beyond (t, program, place). A field absent from a given row is `missing`.
function _trajectory_field_names(rows)
    names = Symbol[]
    for (_, _, _, fields) in rows
        for k in keys(fields)
            k in names || push!(names, k)
        end
    end
    return names
end

"""
    token_trajectory(state) -> DataFrame
    token_trajectory(state, name::AbstractString) -> DataFrame
    token_trajectory(state, pred::TokenPredicate) -> DataFrame

The per-token trajectory log in long form: columns `t, program, place, <field>…`, one row per
(tick, opted-in token), in append (= `token_sortkey`-per-tick) order (§14.1). With a `name` it is one
token's life; with a `TokenPredicate` (ADR 0008 / §9.5) it is the rows of the tokens that CURRENTLY
match the predicate (the same selection machinery the model dynamics use — `select_tokens`), so e.g.
`token_trajectory(state, @select(Project, phase==:Market))` is the launched cohort's history.
Heterogeneous field sets across kinds are unioned; an absent field is `missing`.
"""
function token_trajectory(state::ReactionNetworkProblem)
    return _trajectory_dataframe(state.token_trajectory)
end

function token_trajectory(state::ReactionNetworkProblem, name::AbstractString)
    rows = filter(r -> r[2] == String(name), state.token_trajectory)
    return _trajectory_dataframe(rows)
end

function token_trajectory(state::ReactionNetworkProblem, pred::TokenPredicate)
    names = Set(String(AlgebraicAgents.getname(t)) for t in select_tokens(state, pred))
    rows = filter(r -> r[2] in names, state.token_trajectory)
    return _trajectory_dataframe(rows)
end

function _trajectory_dataframe(rows)
    fieldnames = _trajectory_field_names(rows)
    df = DataFrame(t = Float64[], program = String[], species = Symbol[])
    for f in fieldnames
        df[!, f] = Vector{Any}()
    end
    for (t, name, place, fields) in rows
        row = Dict{Symbol, Any}(:t => t, :program => name, :species => place)
        for f in fieldnames
            row[f] = haskey(fields, f) ? fields[f] : missing
        end
        push!(df, NamedTuple(row); cols = :subset)
    end
    return df
end

# ── "Typical" trajectory helpers (§14.1, two precise senses over the store) ──────────────

# Numeric coordinates of one token's logged path: each NUMERIC field's value at each logged tick,
# flattened in (tick, field) order. Non-numeric fields (a Symbol `phase`) are skipped for the metric
# but still available in the DataFrame view. Returns (name => Vector{Float64}); paths of differing
# length are compared on their common prefix by the metric below (tokens aligned by tick index).
function _numeric_paths(state::ReactionNetworkProblem, names)
    paths = Dict{String, Vector{Float64}}()
    by_name = Dict{String, Vector{Tuple{Float64, String, Symbol, NamedTuple}}}()
    for r in state.token_trajectory
        r[2] in names || continue
        push!(get!(by_name, r[2], eltype(state.token_trajectory)[]), r)
    end
    for (name, rows) in by_name
        coords = Float64[]
        for (_, _, _, fields) in rows
            for v in values(fields)
                v isa Real && push!(coords, Float64(v))
            end
        end
        paths[name] = coords
    end
    return paths
end

# Per-field normalized L2 distance between two equal-or-truncated coordinate vectors (the medoid
# heuristic of ADR 0013's open question — documented as a heuristic, not a canonical "typical").
function _path_distance(a::Vector{Float64}, b::Vector{Float64})
    n = min(length(a), length(b))
    n == 0 && return Inf
    s = 0.0
    @inbounds for i in 1:n
        d = a[i] - b[i]
        s += d * d
    end
    return sqrt(s / n)
end

"""
    representative_token(state[, pred]) -> String

The MEDOID program (ADR 0013 §A5, sense (i) — "show me a typical program's life"): the token whose
numeric logged path is closest to the per-cohort mean path, by per-field normalized distance. The
cohort is all opted-in tokens, or those matching `pred` (a `TokenPredicate`). Returns the token name
(index into `token_trajectory(state, name)`), or `nothing` if the cohort logged nothing. The metric
is a documented heuristic (open question), uniform-weighted over numeric fields.
"""
function representative_token(state::ReactionNetworkProblem)
    names = unique(r[2] for r in state.token_trajectory)
    return _medoid(state, names)
end

function representative_token(state::ReactionNetworkProblem, pred::TokenPredicate)
    names = [String(AlgebraicAgents.getname(t)) for t in select_tokens(state, pred)]
    return _medoid(state, names)
end

function _medoid(state::ReactionNetworkProblem, names)
    paths = _numeric_paths(state, Set(names))
    isempty(paths) && return nothing
    ks = sort!(collect(keys(paths)))                      # deterministic order (§4 D4)
    maxlen = maximum(length(paths[k]) for k in ks)
    # cohort mean path over the common (max) length, treating short paths as covering their prefix
    meanpath = zeros(Float64, maxlen)
    counts = zeros(Int, maxlen)
    for k in ks, (i, v) in enumerate(paths[k])
        meanpath[i] += v
        counts[i] += 1
    end
    for i in 1:maxlen
        counts[i] > 0 && (meanpath[i] /= counts[i])
    end
    best, bestd = ks[1], Inf
    for k in ks
        d = _path_distance(paths[k], meanpath)
        d < bestd && (best, bestd = k, d)
    end
    return best
end

"""
    trajectory_envelope(state[, pred]; align_on = :t) -> DataFrame

The ENVELOPE (ADR 0013 §A5, sense (ii) — "typical ± spread"): for each numeric logged field, the
per-alignment-index median and inter-quartile band (q25/q75) across the cohort. The cohort is all
opted-in tokens or those matching `pred`. With `align_on = :t` (default) tokens are aligned by
absolute tick; pass another field symbol to align by an event index (e.g. ticks-since-first-row).
Columns: `align, field, median, q25, q75, n`. This is also the single-run input to the ensemble
band (ADR 0014 recipe 6). Pure post-processing over the store; adds no engine state.
"""
function trajectory_envelope(state::ReactionNetworkProblem; align_on::Symbol = :t)
    names = unique(r[2] for r in state.token_trajectory)
    return _envelope(state, names; align_on = align_on)
end

function trajectory_envelope(state::ReactionNetworkProblem, pred::TokenPredicate; align_on::Symbol = :t)
    names = [String(AlgebraicAgents.getname(t)) for t in select_tokens(state, pred)]
    return _envelope(state, names; align_on = align_on)
end

function _envelope(state::ReactionNetworkProblem, names; align_on::Symbol = :t)
    nameset = Set(names)
    # Group rows by token, in append order, so per-token alignment index is well-defined.
    by_name = Dict{String, Vector{Tuple{Float64, String, Symbol, NamedTuple}}}()
    for r in state.token_trajectory
        r[2] in nameset || continue
        push!(get!(by_name, r[2], eltype(state.token_trajectory)[]), r)
    end
    fieldnames = filter(
        f -> any(haskey(fields, f) && fields[f] isa Real for (_, _, _, fields) in state.token_trajectory),
        _trajectory_field_names(state.token_trajectory),
    )
    # alignment index → field → collected values
    bucket = Dict{Any, Dict{Symbol, Vector{Float64}}}()
    for (_, rows) in by_name
        for (k, (t, _, _, fields)) in enumerate(rows)
            align = align_on === :t ? t : (k - 1)            # absolute tick or per-token event index
            fb = get!(bucket, align, Dict{Symbol, Vector{Float64}}())
            for f in fieldnames
                if haskey(fields, f) && fields[f] isa Real
                    push!(get!(fb, f, Float64[]), Float64(fields[f]))
                end
            end
        end
    end
    df = DataFrame(
        align = Float64[], field = Symbol[], median = Float64[],
        q25 = Float64[], q75 = Float64[], n = Int[]
    )
    for align in sort!(collect(keys(bucket)))
        for f in fieldnames
            vals = get(bucket[align], f, Float64[])
            isempty(vals) && continue
            push!(
                df, (
                    Float64(align), f, median(vals),
                    quantile(vals, 0.25), quantile(vals, 0.75), length(vals),
                )
            )
        end
    end
    return df
end

# ════════════════════════════════════════════════════════════════════════════════════════
# §14.2 — the ensemble runner + EnsembleProblem (an AA-readable result node)
# ════════════════════════════════════════════════════════════════════════════════════════

"""
    EnsembleProblem <: AbstractAlgebraicAgent

The result of `ensemble` — a `FreeAgent`-style container whose `inners` are the member
`ReactionNetworkProblem`s (entangled), plus the realized per-member seeds and the run mode. Because
it implements the ADR 0012 read surface (`observables`/`getobservable`), an ensemble is itself a
readable/drawable AA hierarchy node (CONTRACT §14.2 Invariant 4) — built ON AA's primitives, not by
modifying AA. `getobservable(ens, name)` returns a cross-run reduction (see `getobservable` below).
"""
@aagent struct EnsembleProblem
    members::Vector{ReactionNetworkProblem}
    seeds::Vector{UInt64}
    root_seed::Int
    mode::Symbol            # :rebuild (mode a) or :reinit (mode b, reinit-reseed member reuse)
end

"""
    ensemble(build; nseed, root_seed = 2026, max_t = nothing, parallel = false, mode = :rebuild) -> EnsembleProblem

Run `nseed` INDEPENDENT members of a scenario. `build(seed) -> ReactionNetworkProblem` constructs ONE
member; member `k ∈ 1:nseed` is seeded `hash((root_seed, k))` (§4 D8) and, if `max_t` is given, run
to `max_t` (otherwise `build` is assumed to already simulate, as the demo's `run_scenario` does).
Each member owns its own `state.rng` (no shared global RNG), so the result is order- and
parallelism-independent (§4 D9) and `(root_seed, nseed, build)` fully determines the ensemble
(Invariant 3).

Two member-production modes, recorded on the result as `ens.mode`:

  - `mode = :rebuild` (mode a, the default): `build(seed)` constructs a FRESH problem per member —
    robust because the structured token population is re-instantiated cleanly from the declarative
    `population[]` (ADR 0007 §B). Pays `nseed ×` construction + closure-compilation cost.
  - `mode = :reinit` (mode b, reinit-reseed member reuse): build ONE member, then reinit-reseed and
    re-simulate that SAME problem for each subsequent seed (`AlgebraicAgents.reinit!(m; seed)`, ADR
    0007 §D + the reseed path), reusing its compiled closures and allocated store — the cheaper Monte
    Carlo path. After each member's run a faithful `deepcopy` SNAPSHOT of the finished problem is
    retained (so `ens.members` still holds `nseed` independent `ReactionNetworkProblem`s the read
    surface + any `metric(member)` see unchanged); the reused problem is reinit-reseeded onward.

**Mode-(a) ≡ mode-(b) equivalence + its precondition.** Mode (b) is a legal substitute for mode (a)
ONLY when members are structurally HOMOGENEOUS — same net, same `population[]` schema, differing only
in the stochastic stream and the seed-sampled initial attributes. The reseed installs the new seed's
stream BEFORE the t=0 marking is re-sampled (`_reinit!`), so a reseeded member is identical to a fresh
`build(seed)` and the two modes produce the ensemble member-for-member. A `build` that BRANCHES
STRUCTURALLY on its seed argument (a different net/population per seed) is out of contract for mode
(b) and must use `mode = :rebuild`: mode (b) reseeds, it does NOT rebuild. This guard is documented,
not auto-detected — `build(root_seed)` is called once for member 1 and reused thereafter.

`parallel = true` is accepted (members are independent) but currently runs sequentially; a threaded
backend is future work (the API is fixed so callers need not change) — and it is only meaningful for
`:rebuild` (mode (b) serially reuses one object). Subsumes `demo/bd_acquisition/analysis.jl`'s
hand-rolled `ensemble`.
"""
function ensemble(
        build; nseed::Integer, root_seed::Integer = 2026,
        max_t = nothing, parallel::Bool = false, name = "ensemble",
        mode::Symbol = :rebuild
    )
    mode in (:rebuild, :reinit) ||
        error("ensemble: mode must be :rebuild (mode a) or :reinit (mode b), got $(repr(mode)).")
    seeds = UInt64[hash((root_seed, k)) for k in 1:nseed]
    members = _ensemble_members(build, seeds, max_t, mode)
    ens = EnsembleProblem(name, members, seeds, Int(root_seed), mode)
    # Make each member a child so the ensemble is a genuine AA hierarchy node (drawable/walkable,
    # §14.2 Invariant 4 / ADR 0014 composition). `entangle!` keys `inners` by `getname`, and members
    # built by the same `build` closure all carry the default name "reaction_network" — so they MUST
    # be renamed to a unique `member_<k>` first, else each entangle! overwrites the previous key and
    # the hierarchy collapses to one child (the drawable-node invariant `length(inners) == nseed`).
    # The per-run statistics read `ens.members` (the vector), so renaming is purely for the AA view.
    for (k, m) in enumerate(members)
        m.name = "member_$k"
        entangle!(ens, m)
    end
    return ens
end

# Mode (a): REBUILD a fresh problem per seed (the closure re-runs construction + compilation each
# time). `build` may already simulate; `max_t` runs it if given.
function _ensemble_members(build, seeds, max_t, mode::Val{:rebuild})
    members = ReactionNetworkProblem[]
    for s in seeds
        push!(members, _build_member(build, s, max_t))
    end
    return members
end

# Mode (b): build ONE member, then reinit-reseed the SAME problem for each subsequent seed and keep a
# faithful deepcopy snapshot of each finished run. The retained snapshots are detached, fully-formed
# `ReactionNetworkProblem`s — `sol`/`log`/`program_ledger`/`token_trajectory`/`observables` all read
# off them exactly as off a rebuilt member — so `ens.members`, `summarize`, `treatment_effect`, and
# any `metric(::ReactionNetworkProblem)` are contract-identical to mode (a). `build(seeds[1])` is
# called ONCE (the reused object); a `build` that branches structurally on its seed is out of
# contract here (see the `ensemble` docstring's equivalence precondition).
function _ensemble_members(build, seeds, max_t, mode::Val{:reinit})
    isempty(seeds) && return ReactionNetworkProblem[]
    prob = _build_member(build, first(seeds), max_t)
    members = ReactionNetworkProblem[deepcopy(prob)]
    for s in seeds[2:end]
        AlgebraicAgents.reinit!(prob; seed = s)
        max_t === nothing ? simulate(prob) : simulate(prob, max_t)
        push!(members, deepcopy(prob))
    end
    return members
end

_ensemble_members(build, seeds, max_t, mode::Symbol) =
    _ensemble_members(build, seeds, max_t, Val(mode))

function _build_member(build, s, max_t)
    m = build(s)
    m isa ReactionNetworkProblem ||
        error(
        "ensemble: build(seed) must return a ReactionNetworkProblem (got $(typeof(m))); " *
            "for a coupled member, extract the RD child before returning (ADR 0013 open Q)."
    )
    max_t === nothing || simulate(m, max_t)
    return m
end

"""
    summarize(ens, metric) -> NamedTuple

Reduce a per-run scalar `metric(member::ReactionNetworkProblem) -> Real` across the ensemble:
`(; mean, sem, q25, median, q75, n)` (CONTRACT §14.2). The standard error of the mean is the honest
spread the ensemble exists to report (`sem = std/√n`, 0 for n ≤ 1). Ledger statistics are e.g.
`summarize(ens, p -> last(program_ledger(p).valuation))`.
"""
function summarize(ens::EnsembleProblem, metric)
    xs = Float64[Float64(metric(m)) for m in ens.members]
    n = length(xs)
    return (
        mean = n == 0 ? NaN : mean(xs),
        sem = n <= 1 ? 0.0 : std(xs) / sqrt(n),
        q25 = n == 0 ? NaN : quantile(xs, 0.25),
        median = n == 0 ? NaN : median(xs),
        q75 = n == 0 ? NaN : quantile(xs, 0.75),
        n = n,
    )
end

"""
    treatment_effect(ens_baseline, ens_deal, metric) -> NamedTuple

The unpaired treatment effect on `metric` between two INDEPENDENT ensembles (CONTRACT §14.2): the
difference of means with the unpaired SE `se = sqrt(var_b/n_b + var_d/n_d)` (the two arms desync the
shared-RNG-free streams, so this is the unpaired estimator — MVP §4.1 finding A). Returns
`(; delta, se, baseline, deal, n_baseline, n_deal)`. This is the A/B lever comparison that drove the
BD Δ-rNPV (`treatment_effect(ensemble(baseline), ensemble(deal), rnpv)`).
"""
function treatment_effect(ens_baseline::EnsembleProblem, ens_deal::EnsembleProblem, metric)
    b = Float64[Float64(metric(m)) for m in ens_baseline.members]
    d = Float64[Float64(metric(m)) for m in ens_deal.members]
    nb, nd = length(b), length(d)
    se = (nb <= 1 || nd <= 1) ? 0.0 : sqrt(var(b) / nb + var(d) / nd)
    return (
        delta = mean(d) - mean(b),
        se = se,
        baseline = mean(b),
        deal = mean(d),
        n_baseline = nb,
        n_deal = nd,
    )
end

# ── EnsembleProblem as an AA-readable node (ADR 0012 read surface, §14.2 Invariant 4) ────

# AA stepping hooks: an ensemble is a finished RESULT, not a steppable dynamical system — it is a
# read-only container. `_step!` is a no-op and it never projects forward (mirrors the passive-token
# pattern, agents.jl).
AlgebraicAgents._step!(::EnsembleProblem) = nothing
AlgebraicAgents._projected_to(::EnsembleProblem) = nothing

"""
    observables(ens::EnsembleProblem)

The cross-run aggregate names an ensemble exports (ADR 0012 §A / §14.2 Invariant 4): the union of
its members' `observables`, each surfaced as the across-member MEAN by `getobservable`. Stable,
sorted order.
"""
function AlgebraicAgents.observables(ens::EnsembleProblem)
    names = Symbol[]
    for m in ens.members
        for o in AlgebraicAgents.observables(m)
            o in names || push!(names, o)
        end
    end
    return sort!(names)
end

"""
    getobservable(ens::EnsembleProblem, name)

The across-member MEAN of member observable `name` (ADR 0012 §A / §14.2 Invariant 4): reads
`getobservable(member, name)` on each member and averages. An unknown name is a hard `error` (a
diagnostic, never AA's silent `@error` fall-through). `Int` indexes `observables(ens)`.
"""
function AlgebraicAgents.getobservable(ens::EnsembleProblem, name::Symbol)
    name in AlgebraicAgents.observables(ens) || error(
        "getobservable: `$name` is not an exported aggregate of $(getname(ens)); " *
            "exported names are $(AlgebraicAgents.observables(ens)) (§14.2 Invariant 4)",
    )
    vals = Float64[]
    for m in ens.members
        name in AlgebraicAgents.observables(m) || continue
        push!(vals, Float64(AlgebraicAgents.getobservable(m, name)))
    end
    return isempty(vals) ? NaN : mean(vals)
end

AlgebraicAgents.getobservable(ens::EnsembleProblem, name::AbstractString) =
    AlgebraicAgents.getobservable(ens, Symbol(name))

function AlgebraicAgents.getobservable(ens::EnsembleProblem, i::Int)
    names = AlgebraicAgents.observables(ens)
    checkbounds(Bool, names, i) ||
        error("getobservable: index $i out of range 1:$(length(names)) (§14.2)")
    return AlgebraicAgents.getobservable(ens, names[i])
end

# ════════════════════════════════════════════════════════════════════════════════════════
# Result-plot specs (ADR 0014 §15.1) — the wrapper types the RDPlotsExt @recipes dispatch on
# ════════════════════════════════════════════════════════════════════════════════════════
#
# Plots.jl recipes need a TYPE to dispatch on; dispatching on raw `prob.sol` (a DataFrame) would
# hijack every DataFrame plot. So each recipe is keyed off a thin wrapper constructed here in the
# CORE (no Plots dependency — a user can name `MarkingPlot(prob)` with Plots absent) and given its
# `@recipe` in ext/RDPlotsExt.jl. Each wraps a raw result artifact (Invariant 5: model-agnostic, no
# hard-coded BD field), so `plot(MarkingPlot(prob))` works on ANY model. The four bespoke BD figures
# collapse onto `EnsembleBar`/`TreatmentEffectPlot` (recipe 6) and `SaturationPlot` (recipe 2).

export MarkingPlot, SaturationPlot, ValuationPlot, LedgerPlot,
    TokenTrajectoryPlot, EnsembleBar, TreatmentEffectPlot, ThroughputPlot

"""
    MarkingPlot(prob; vars = all places)

Plot spec (ADR 0014 recipe 1) for per-place token COUNTS over time — the marking trajectory of the named `vars` across a finished run's `prob.sol`. Realized by a `@recipe` in `RDPlotsExt`; `plot(MarkingPlot(prob))` needs `Plots` loaded.
"""
struct MarkingPlot              # recipe 1 — per-place token counts over time (generalizes _draw)
    prob::ReactionNetworkProblem
    vars::Vector{String}
end
MarkingPlot(prob::ReactionNetworkProblem; vars = string.(prob.network[:, :placeName])) =
    MarkingPlot(prob, collect(String.(vars)))

"""
    SaturationPlot(prob; vars = all place)

Plot spec (ADR 0014 recipe 2) for RESOURCE UTILIZATION over time — the troughs of the named resource pools `vars` across a finished run, showing when a `@conserved`/`@rate` resource is drawn down (saturated). Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct SaturationPlot           # recipe 2 — resource utilization / pool troughs
    prob::ReactionNetworkProblem
    vars::Vector{String}
end
SaturationPlot(prob::ReactionNetworkProblem; vars = string.(prob.network[:, :placeName])) =
    SaturationPlot(prob, collect(String.(vars)))

"""
    ValuationPlot(prob)

Plot spec (ADR 0014 recipe 3) for the portfolio VALUATION curve — the cumulative cost/reward/valuation series the ledger logged over a finished run. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct ValuationPlot            # recipe 3 — portfolio valuation / cost / reward curve
    prob::ReactionNetworkProblem
end

"""
    LedgerPlot(prob)

Plot spec (ADR 0014 recipe 4) for PER-PROGRAM cost/reward bars — the final per-program totals from `program_ledger(prob)`, one bar group per program. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct LedgerPlot               # recipe 4 — per-program cost/reward bars from program_ledger
    prob::ReactionNetworkProblem
end

"""
    TokenTrajectoryPlot(prob, field; pred = nothing)

Plot spec (ADR 0014 recipe 5) for one logged `field`'s PER-TOKEN paths over time, overlaid with the cohort's typical envelope band (median + IQR, from [`trajectory_envelope`](@ref)). The cohort is all opted-in tokens, or those matching the `TokenPredicate` `pred`. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct TokenTrajectoryPlot      # recipe 5 — a field's per-token paths + the typical envelope band
    prob::ReactionNetworkProblem
    field::Symbol
    pred::Union{Nothing, TokenPredicate}
end
TokenTrajectoryPlot(prob::ReactionNetworkProblem, field::Symbol; pred = nothing) =
    TokenTrajectoryPlot(prob, field, pred)

"""
    EnsembleBar(ens, metric)

Plot spec (ADR 0014 recipe 6a) for the DISTRIBUTION of a per-run scalar `metric(member) -> Real` across an ensemble's members — a histogram of the metric over the [`ensemble`](@ref) runs. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct EnsembleBar              # recipe 6a — per-run metric distribution (histogram) across members
    ens::EnsembleProblem
    metric::Any
end

"""
    TreatmentEffectPlot(baseline, deal, metric)

Plot spec (ADR 0014 recipe 6b) for an A/B comparison: the `metric` distributions of the `baseline` and `deal` ensembles side by side, with the treatment effect Δ (see [`treatment_effect`](@ref)) annotated. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct TreatmentEffectPlot      # recipe 6b — baseline vs deal metric distributions, Δ annotated
    baseline::EnsembleProblem
    deal::EnsembleProblem
    metric::Any
end

"""
    ThroughputPlot(prob)

Plot spec (ADR 0014 recipe 7) for THROUGHPUT over time — the number of transition firings and terminations per tick across a finished run. Realized by a `@recipe` in `RDPlotsExt`; needs `Plots` loaded.
"""
struct ThroughputPlot           # recipe 7 — firings / terminations per tick
    prob::ReactionNetworkProblem
end

# Helper the recipes (and callers) use: pull a scalar-per-tick series tagged `tag` out of the log
# (`:valuation`/`:valuation_cost`/`:valuation_reward` are scalar rows). Returns (t::Vector, v::Vector).
function log_scalar_series(prob::ReactionNetworkProblem, tag::Symbol)
    ts = Float64[]
    vs = Float64[]
    for r in prob.log
        r[1] === tag || continue
        push!(ts, Float64(r[2]))
        push!(vs, Float64(r[3]))
    end
    return ts, vs
end

# Helper: per-tick count of a (key, q)-style row family — the number of firing/terminating
# instances summed per tick. The element shapes DIFFER by tag: `:new_transitions` splats `(hash, q)`
# Tuples (solvers.jl), while `:terminated_all`/`:terminated_success` splat `Symbol=>Float64` Pairs
# (a Dict, solvers.jl) — a `Pair` is NOT `<: Tuple`, so we must accept both or terminations read as 0.
function log_count_series(prob::ReactionNetworkProblem, tag::Symbol)
    ts = Float64[]
    vs = Float64[]
    _qty(x) = (x isa Pair || x isa Tuple) ? Float64(last(x)) : 0.0
    for r in prob.log
        r[1] === tag || continue
        push!(ts, Float64(r[2]))
        push!(vs, sum(_qty, r[3:end]; init = 0.0))
    end
    return ts, vs
end
