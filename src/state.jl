@reexport using AlgebraicAgents
using DataFrames
using Random

# Per-program (per-structured-token) ledger accumulator (MVP finding D — the logic lives in
# src/ledger.jl; only this small data struct is here so the ReactionNetworkProblem field below can
# name its element type, since state.jl is `include`d before ledger.jl). `cost_incurred` is capital
# burned on behalf of this program; `reward_realized` is reward credited when a transition it was
# bound to finished successfully; `valuation` is its current mark-to-market; `entries` is the
# append-only `(t, kind, amount, transition_name)` audit trail. See ledger.jl for the attribution
# rule and its documented boundary.
"""
    ProgramLedger(species, creation_index)

Per-program (per-structured-token) ledger accumulator (CONTRACT §12, MVP finding D — the attribution logic lives in `src/ledger.jl`). Tracks one structured token's running economics: `cost_incurred` (capital burned on its behalf), `reward_realized` (reward credited when a transition it was bound to finished successfully), `valuation` (current mark-to-market), and `entries` — the append-only `(t, kind, amount, transition_name)` audit trail. `species`/`creation_index` identify the program. Per-program rows plus the state's `unattributed_cost`/`unattributed_reward` buckets sum exactly to the aggregate ledger rows.
"""
mutable struct ProgramLedger
    species::Symbol
    creation_index::Int
    cost_incurred::Float64
    reward_realized::Float64
    valuation::Float64
    entries::Vector{Tuple{Float64, Symbol, Float64, String}}
end
ProgramLedger(species::Symbol, creation_index::Int) =
    ProgramLedger(species, creation_index, 0.0, 0.0, 0.0, Tuple{Float64, Symbol, Float64, String}[])

struct UnfoldedReactant
    index::Int
    species::Symbol
    stoich::ActionableValues
    modality::Set{Symbol}
    predicate::Any   # nothing (kind-only bind, the default) or a TokenPredicate (ADR 0008 §B)
end
# Backward-compatible constructor: no predicate ⇒ today's kind-only binding.
UnfoldedReactant(index, species, stoich, modality) =
    UnfoldedReactant(index, species, stoich, modality, nothing)

"""
Ongoing transition auxiliary structure.
"""
@aagent struct Transition
    i::Int

    trans::Dict{Symbol, Any}

    bound_structured_agents::Vector{AbstractAlgebraicAgent}
    nonblock_structured_agents::Vector{AbstractAlgebraicAgent}
    structured_to_agents::Vector

    t::Float64
    q::Float64
    state::Float64
end

Base.getindex(state::Transition, key) = state.trans[key]
Base.setindex!(state::Transition, val, key) = state.trans[key] = val

@aagent struct Observable
    last::Float64 # last sampling time
    range::Vector{Union{Tuple{Float64, SampleableValues}, SampleableValues}}
    every::Float64
    on::Vector{ActionableValues}

    sampled::Any
end

"""
The live simulation state — an AlgebraicAgents `@aagent`, so a running network is itself a node in a larger heterogeneous AA hierarchy (an `AbstractAlgebraicAgent`). It is constructed from a static authoring store by the `ReactionNetworkProblem(net; …)` outer constructor and advanced by the `_step!` loop. Key fields: `.network` (the static `ReactionNetwork` IR store the run was compiled from), `.u` (the current plain-species marking vector), `.p` (parameters), `.t`/`.tspan`/`.dt` (time control), `.sol` (the per-step marking log as a `DataFrame`), `.log` (the event/message log), `.observables`, `.ongoing_transitions` (in-flight transition instances), `.program_ledgers` (per-program economics, §12), and `.token_trajectory` (per-token trajectory log, §14.1). Determinism is contractual (§4): `.rng` is the state-owned RNG that is the SOLE source of randomness in the step loop, `.seed` records the realized construction seed, and `.initial_rng` snapshots the stream at t=0 so `_reinit!` restores it exactly. The endogenous decision channel lives in `.rules`/`.registry`; the runtime store is append-only (ADR 0004), so compiled attribute closures may position-index it safely.
"""
@aagent struct ReactionNetworkProblem
    network::ReactionNetwork

    attrs::Dict{Symbol, Vector}
    transition_recipes::Dict{Symbol, Vector}

    u::Vector{Float64}
    p::Any
    t::Float64

    structured_token::Vector{Symbol}

    tspan::Tuple{Float64, Float64}
    dt::Float64

    transitions::Dict{Symbol, Vector}
    ongoing_transitions::Vector{Transition}
    log::Vector{Tuple}

    observables::Dict{Symbol, Observable}

    wrap_fun::Any
    sol::DataFrame

    # Determinism (§4 D1–D9): a state-owned RNG is the SOLE source of randomness in the
    # step loop. `seed` records the REALIZED construction seed (the explicit `seed=` kwarg, or
    # the entropy-drawn one when none was given — always concrete, so any run is replayable, D6).
    # `initial_rng` is a snapshot of `rng` at t=0 so `_reinit!` restores the exact stream (D7).
    rng::Random.AbstractRNG
    seed::Union{Integer, Nothing}
    initial_rng::Random.AbstractRNG

    # Endogenous decision channel (ADR 0010/0011, §12). `rules` are guard/action/fire_mode
    # triples fired at _step! step 10. `registry` is the per-network host-function/kind
    # allow-list (ADR 0006 §C) keyed by name — used by AddToken/Invoke; never eval'd.
    rules::Vector
    registry::Dict{Symbol, Any}

    # Structured-token determinism (ADR 0006 §E / ADR 0008): per-species monotonic creation
    # counter, and the realized (token-name → creation_index) map that fixes the (species,
    # creation_index) total order tokens are selected in. Reset by _reinit! (§4 D7).
    creation_counters::Dict{Symbol, Int}
    creation_index::Dict{String, Int}

    # Declarative initial marking (ADR 0007 §B). `population` is the structured-token initial
    # state (the analogue of specInitVal for plain species): either a vector of declarative
    # PopulationEntry specs (count + seeded attribute exprs) OR already-constructed host token
    # agents. Stored so _reinit! can rebuild the exact t=0 marking (§D, closing §4 D7 for
    # structured runs). For the explicit-host-token form `init_snapshot` records each token's
    # initial field values (by token name) so _reinit! can restore the SAME objects to their t=0
    # attributes (species/phase/…), not just reset their bonds. `live` arms the §A phase guard:
    # once constructed, reindexers (rem_parts!) refuse.
    population::Vector
    init_snapshot::Dict{String, Dict{Symbol, Any}}
    live::Bool

    # Per-program (per-structured-token) ledger (MVP finding D — src/ledger.jl). `program_ledgers`
    # maps a token's stable network identity (`AlgebraicAgents.getname`, the same key as
    # `creation_index`) to its running cost/reward/valuation accumulator + append-only audit trail.
    # `unattributed_cost`/`unattributed_reward` collect spend/reward from transitions with NO bound
    # structured token (plain reactions), so the per-program rows + these buckets SUM EXACTLY to the
    # aggregate `:valuation_cost`/`:valuation_reward` rows (the invariant in program_ledger.jl).
    # All three are reset by _reinit! (§4 D7), mirroring `creation_counters`.
    program_ledgers::Dict{String, ProgramLedger}
    unattributed_cost::Float64
    unattributed_reward::Float64

    # External coupling buffer (ADR 0012 §B3). The per-tick latch of declared `inputs[]` ports:
    # `_prestep!` reads RD's incoming AA wires ONCE per tick (before any sibling `_step!`) and
    # merges them over the declared input DEFAULTS into this Dict, so every `ExternalRef(port)`
    # read within a tick returns the SAME value — the source's previous-tick-boundary projection
    # (Invariant 2, the §4-D4 sibling-order hazard the latch closes). Transient run-state: seeded
    # from the declared defaults at construction, re-seeded by `_reinit!` (§10.4), recomputed every
    # `_prestep!`, and NOT persisted by `dump_state` (recovered on the next prestep).
    external_inputs::Dict{Symbol, Any}
    # The immutable declared-default snapshot (the §B1 `inputs[]` defaults), paired with
    # `external_inputs` exactly as `initial_rng` is paired with `rng` (§4 D7): `_reinit!` restores
    # `external_inputs` to a copy of this so a re-run drops stale latched wire values but KEEPS the
    # pre-wire defaults, and every `_prestep!` merges wire reads OVER a copy of it.
    external_input_defaults::Dict{Symbol, Any}

    # Per-token trajectory log (ADR 0013 §A / CONTRACT §14.1). The time-indexed companion to the
    # per-program ledger: each tick `push_token_trajectory_row!` appends `(t, token_name, species,
    # fields)` for every token whose KIND opts in via `log_token_fields(tok)::NamedTuple` (default
    # empty), iterated in `token_sortkey` order — the same seam (`solvers.jl`, right after
    # `push_program_ledger_row!`), observation point, and determinism guarantee as the ledger row it
    # generalizes (§4 D4). `species` is captured at log time (it can change under soft-retire). Sibling
    # of `log`; bounded by per-kind opt-in (Invariant 2); reset by `_reinit!` like the ledger (§4 D7).
    token_trajectory::Vector{Tuple{Float64, String, Symbol, NamedTuple}}
end

# get value of a numeric expression
# evaluate compiled numeric expression in context of (u, p, t)
function context_eval(state::ReactionNetworkProblem, transition, o)
    o = o isa Function ? Base.invokelatest(o, state, transition) : o

    return o isa Sampleable ? rand(state.rng, o) : o
end

function Base.getindex(state::ReactionNetworkProblem, keys...)
    if any(occursin.(["transPreAction", "transPostAction"], Ref(string(keys[2]))))
        return state.network[keys[1], keys[2]]
    else
        return context_eval(
            state,
            nothing,
            (contains(string(keys[2]), "trans") ? state.transitions : state.attrs)[keys[2]][keys[1]],
        )
    end
end

function init_u!(state::ReactionNetworkProblem)
    return (
        u = fill(0.0, nrows(state, :S));
        foreach(i -> u[i] = state[i, :specInitVal], row_ids(state, :S));
        state.u = u
    )
end
save!(state::ReactionNetworkProblem) = push!(state.sol, (state.t, state.u[:]...))

function compile_observables(net::ReactionNetwork)
    observables = Dict{Symbol, Observable}()
    species_names = collect(net[:, :specName])
    prm_names = collect(net[:, :prmName])
    varmap = Dict([name => :(state.u[$i]) for (i, name) in enumerate(species_names)])

    for (name, opts) in Iterators.zip(net[:, :obsName], net[:, :obsOpts])
        on = map(on -> wrap_expr(on, species_names, prm_names, varmap), opts.on)
        range = map(
            r -> begin
                r = r isa Tuple ? r : (1.0, r)
                (r[1], wrap_expr(r[2], species_names, prm_names, varmap))
            end,
            opts.range,
        )

        push!(
            observables,
            name => Observable(string(name), -Inf, range, opts.every, on, missing),
        )
    end

    return observables
end

compileval(ex, state) = !isa(ex, Expr) ? ex : eval(state.wrap_fun(ex))

# `rng` here is the range vector (a historical misnomer); randomness is drawn from the
# state-owned `state.rng` (§4 D2/D5), never the global RNG.
function sample_range(rng, state)
    isempty(rng) && return missing
    r = rand(state.rng) * sum(r -> r isa Tuple ? eval(r[1]) : 1, rng)
    ix = 0
    s = 0
    while s <= r && (ix < length(rng))
        ix += 1
        s += rng[ix] isa Tuple ? compileval(rng[ix][1], state) : 1
    end

    r = rng[ix] isa Tuple ? rng[ix][2] : rng[ix]
    return r isa Sampleable ? rand(state.rng, r) : r
end

function resample!(state::ReactionNetworkProblem, o::Observable)
    o.last = state.t
    # Range-less observable resamples to `missing` on the `.sampled` field (the struct has
    # no `.val` field — that was the bug pinned by determinism_composition_bugs.jl).
    isempty(o.range) && (return o.sampled = missing)

    return o.sampled = context_eval(state, nothing, sample_range(o.range, state))
end

resample(state::ReactionNetworkProblem, o::Symbol) = resample!(state, state.observables[o])

function update_observables(state::ReactionNetworkProblem)
    return foreach(
        o -> (state.t - o.last) >= o.every && resample!(state, o),
        values(state.observables),
    )
end

function prune_r_line(r_line)
    return if r_line isa Expr && r_line.args[1] ∈ fwd_arrows
        r_line.args[[2, 3]]
    elseif r_line isa Expr && r_line.args[1] ∈ bwd_arrows
        r_line.args[[3, 2]]
    elseif isexpr(r_line, :macrocall) && (macroname(r_line) == :choose)
        sample_range(
            [
                (
                        if isexpr(r, :tuple)
                            (r.args[1], prune_r_line(r.args[2]))
                    else
                            prune_r_line(r)
                    end
                    ) for r in r_line.args[3:end]
            ],
            state,
        )
    end
end

function find_index(species::Symbol, state::ReactionNetworkProblem)
    return findfirst(i -> state[i, :specName] == species, row_ids(state, :S))
end

function sample_transitions!(state::ReactionNetworkProblem)
    for (_, v) in state.transitions
        empty!(v)
    end
    for i in 1:length(state.transition_recipes[:trans])
        # A transition fires this tick iff it is activated (latching gate, ADR 0004) AND its
        # stateless guard holds this tick (ADR 0010 §B). We ALWAYS realize and push every
        # transition's attributes so `state.transitions[attr]` stays parallel to the `:T`
        # part-index that `evolve!`/`get_allocs!` index by — a non-firing transition is
        # recorded with `transFiring = false` and `evolve!` zeroes its genesis quantity, so
        # it makes no proposal and never competes for resources (clean ADR-0002 interaction).
        fires =
            state.transition_recipes[:transActivated][i] &&
            (context_eval(state, nothing, state.transition_recipes[:transGuard][i]) != false)
        l_line, r_line = prune_r_line(state.transition_recipes[:trans][i])

        for attr in keys(state.transition_recipes)
            (
                attr ∈ [
                    :trans,
                    :transPreAction,
                    :transPostAction,
                    :transActivated,
                    :transHash,
                    :transGuard,
                ]
            ) && continue
            push!(
                state.transitions[attr],
                context_eval(state, nothing, state.transition_recipes[attr][i]),
            )
        end

        reactants = []
        for r in extract_reactants(l_line, state)
            j = find_index(r.species, state)
            push!(
                reactants,
                UnfoldedReactant(
                    j,
                    r.species,
                    context_eval(state, nothing, state.wrap_fun(r.stoich)),
                    r.modality ∪ state[j, :specModality],
                    r.predicate,
                ),
            )
        end

        push!(state.transitions[:transLHS], reactants)
        push!(state.transitions[:transRHS], r_line)
        push!(state.transitions[:transFiring], fires)

        foreach(
            k -> push!(state.transitions[k], state.transition_recipes[k][i]),
            [:transPreAction, :transPostAction, :transToSpawn, :transHash],
        )

        state.transition_recipes[:transToSpawn] .= 0
    end
    return
end

function as_state(u, t, state::ReactionNetworkProblem)
    return (state = deepcopy(state); state.u .= u; state.t = t; state)
end

# Extend RD's own nrows/row_ids store generics (ADR 0003 Phase 1 — no longer ACSets'; renamed in
# ADR 0015) for the live state, delegating to the static store. Defined in ReactiveDynamics.jl.
function nrows(state::ReactionNetworkProblem, obj::Symbol)
    return nrows(state.network, obj)
end

function row_ids(state::ReactionNetworkProblem, obj::Symbol)
    return row_ids(state.network, obj)
end

## query the state

t(state::ReactionNetworkProblem) = state.t
solverarg(state::ReactionNetworkProblem, arg) = state.p[arg]
take(state::ReactionNetworkProblem, pcs::Symbol) = state.observables[pcs].sampled
log(state::ReactionNetworkProblem, msg) = (println(msg); push!(state.log, (:log, msg)))
state(state::ReactionNetworkProblem) = state

function periodic(state::ReactionNetworkProblem, period)
    return period == 0.0 || (
        length(state.sol.t) > 1 &&
            (fld(state.t, period) - fld(state.sol.t[end - 1], period) > 0)
    )
end

set_params(state::ReactionNetworkProblem, vals...) =
    for (p, v) in vals
    state.p[p] = v
end

function add_to_spawn!(state::ReactionNetworkProblem, hash, n)
    ix = findfirst(
        ix -> state.transition_recipes[:transHash][ix] == hash,
        1:length(state.transition_recipes[:transHash]),
    )
    return !isnothing(ix) && (state.transition_recipes[:transToSpawn][ix] += n)
end
