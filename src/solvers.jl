using Distributions
using Random

export ReactionNetworkProblem

function get_sampled_transition(state, i)
    transition = Dict{Symbol, Any}()
    foreach(k -> push!(transition, k => state[i, k]), keys(state.transitions))

    return transition
end

# The token-selection predicate (ADR 0008) for structured place `type` in an LHS arc
# list, or `nothing` (kind-only bind) when that arc carries none.
function lhs_predicate(lhs, type::Symbol)
    for r in lhs
        r isa ResolvedArc && r.place == type && return r.predicate
    end
    return nothing
end

isinteger(x::Number) = x == trunc(x)

# ── Priority-weighted progressive-filling allocator (ADR 0002) ───────────────────────────
#
# Each tick, transition instances compete for the finite shared supply `u[s]`. Demand is
# conjunctive (Leontief): an instance needs ALL of its required place at once, so a partial
# share of one place without the others is wasted. ADR 0002 fixes the policy as priority-
# weighted max-min fairness, computed by weighted progressive filling (water-filling). The
# single allocator below replaces the old seven-function tangle (`get_reqs_init!`,
# `get_reqs_ongoing!`, `get_allocs!`/`alloc_weighted!`/`alloc_greedy!`, `get_frac_satisfied`,
# `get_init_satisfied`) and the `state.p[:strategy]` `:weighted`/`:greedy` switch — strict
# greedy is just the `priority → ∞` limit of weighted filling and needs no separate path.
#
# The allocator is RNG-free and deterministic: no `sort` with unstable ties; the only tie-break
# (integer top-up) is a stable `(-priority, index)` order. It is work-conserving (a resource is
# left idle only once every transition that could use it is frozen) and conjunctive-consistent
# (`alloc[s,t] = f[t]·req[s,t]` by construction — nothing is stranded).

"""
Struct-of-arrays scratch for the progressive-filling allocator (ADR 0002 "Proposed Julia
surface"). One workspace is built per `evolve!` call site each tick; `req` is filled by
`build_requirements!`. The other fields are the iteration state:

  - `req::Matrix` — S×T, units of place `s` per unit fill of transition `t`.
  - `f::Vector`   — T, the per-transition fill fraction (the allocator's output).
  - `r::Vector`   — S, remaining supply during filling.
  - `active`      — T, which transitions are still being filled.
  - `D::Vector`   — S, per-place weighted demand of the currently active transitions.

`AllocWorkspace(req::Matrix)` builds a workspace sized to a given requirement matrix (the
unit-test entry point); `AllocWorkspace(nS, nT)` allocates a zeroed workspace of a given shape.
"""
struct AllocWorkspace
    req::Matrix{Float64}     # S×T, rebuilt each call
    f::Vector{Float64}       # T
    r::Vector{Float64}       # S
    active::BitVector        # T
    D::Vector{Float64}       # S
end

function AllocWorkspace(nS::Integer, nT::Integer)
    return AllocWorkspace(
        zeros(Float64, nS, nT),
        zeros(Float64, nT),
        zeros(Float64, nS),
        falses(nT),
        zeros(Float64, nS),
    )
end

function AllocWorkspace(req::AbstractMatrix)
    return AllocWorkspace(copy(Matrix{Float64}(req)), size(req)...)
end
function AllocWorkspace(req::Matrix{Float64}, nS::Integer, nT::Integer)
    return AllocWorkspace(
        req,
        zeros(Float64, nT),
        zeros(Float64, nS),
        falses(nT),
        zeros(Float64, nS),
    )
end

"""
    build_requirements!(ws, state, qs; counted, dt_scale, ongoing)

Fill `ws.req[s,t]` — the per-unit-fill demand of each transition — from the transitions' LHS
tokens. One parametrized builder for both phases (replaces `get_reqs_init!`/`get_reqs_ongoing!`):

  - **Spawn** (`ongoing=false`): `counted` excludes `:rate` and `:nonblock` (upfront tokens only),
    `dt_scale = 1`. `qs[t]` is the desired spawn count of recipe row `t`.
  - **Ongoing** (`ongoing=true`): include `:rate` tokens scaled by `dt_scale = state.dt` (only when
    the in-flight transition's `transCycleTime > 0`) and `:nonblock` tokens unscaled. `qs[i]` is the
    in-flight instance count of `state.ongoing_transitions[i]`.

The modality rules match the deleted builders exactly; the conjunctive `req[s,t]` here is the
demand per unit fill, so `alloc[s,t] = f[t]·req[s,t]`.
"""
function build_requirements!(
        ws::AllocWorkspace,
        state,
        qs;
        counted = nothing,
        dt_scale = 1,
        ongoing = false,
    )
    reqs = ws.req
    reqs .= 0.0
    if ongoing
        for i in eachindex(state.ongoing_transitions)
            for tok in state.ongoing_transitions[i][:transLHS]
                if in(:rate, tok.modality)
                    in(tok.place, state.structured_token) && error(
                        "Modality `:rate` is not supported for structured place in transition $(state.ongoing_transitions[i][:transName]).",
                    )
                    (state.ongoing_transitions[i][:transCycleTime] > 0) &&
                        (reqs[tok.index, i] += qs[i] * tok.multiplicity * dt_scale)
                end
                in(:nonblock, tok.modality) && (reqs[tok.index, i] += qs[i] * tok.multiplicity)
            end
        end
    else
        for i in axes(reqs, 2)
            for tok in state[i, :transLHS]
                # Spawn counts only "upfront" tokens: everything that is neither `:rate`
                # (consumed continuously while in-flight) nor `:nonblock` (claimed but not held).
                !any(m -> m in tok.modality, [:rate, :nonblock]) &&
                    (reqs[tok.index, i] += qs[i] * tok.multiplicity)
            end
        end
    end

    return reqs
end

"""
    progressive_fill!(ws, u, w; fmax = fill(Inf, length(w))) -> f

Weighted progressive filling (water-filling) — the unified allocator core (ADR 0002 "Algorithm").
Returns the per-transition fill-fraction vector `f` (also stored in `ws.f`); the caller derives
`alloc[s,t] = ws.req[s,t] * f[t]` and debits `state.u .-= sum(alloc; dims=2)`.

Inputs: `ws.req[s,t] ≥ 0` (rebuilt by `build_requirements!`), supply `u[s] ≥ 0`, fill-rate weights
`w[t] = priority[t] ≥ 0`, and per-transition caps `fmax[t]` (desired instances/progress; `Inf`
= fill until a required resource exhausts).

Every active transition's fill grows at rate `w[t]`; at each step the largest step `dτ` is taken
that neither overdraws a resource nor overshoots a cap, then any transition that hit its cap or
that needs a now-saturated resource is frozen. Properties (ADR 0002):

  - **work-conserving** — a resource is idle only once every transition wanting it is frozen;
  - **conjunctive-consistent** — `alloc = f' .* req`, so nothing is stranded;
  - **deterministic** — no RNG, no unstable sort; simultaneous freezes are order-independent.

**Zero-priority "leftover only" — two-stage** (ADR 0002 Resolved #3). `w[t] = 0` means "run only
from genuinely leftover resource": stage 1 fills the positive-priority transitions; stage 2 fills
the zero-priority transitions, at equal weight, from the supply that stage 1 left behind. A
zero-priority transition therefore never competes with positive-priority demand and advances only
if positive-priority demand did not exhaust its required place.
"""
function progressive_fill!(ws::AllocWorkspace, u, w; fmax = fill(Inf, length(w)))
    nT = length(w)
    f = ws.f
    fill!(f, 0.0)

    # Stage 1: positive-priority tier fills from the full supply.
    pos = [t for t in 1:nT if w[t] > 0]
    _fill_tier!(ws, u, w, fmax, pos)

    # Stage 2: zero-priority tier fills the leftover, at equal weight (max-min fair). Recover the
    # residual supply after stage 1, then fill the w[t]==0 transitions with a flat unit weight.
    zero_tier = [t for t in 1:nT if w[t] == 0]
    if !isempty(zero_tier)
        r = ws.r
        copyto!(r, u)
        for t in 1:nT, s in axes(ws.req, 1)
            r[s] -= ws.req[s, t] * f[t]
        end
        @. r = max(0.0, r)
        _fill_tier!(ws, r, fill(1.0, nT), fmax, zero_tier)
    end

    return f
end

# Core water-filling loop restricted to a tier of transition indices `tier`, drawing from supply
# `supply` with weights `weights` and caps `fmax`. Advances `ws.f[t]` for t in the tier in place;
# transitions outside the tier are untouched. Terminates in at most S+|tier| iterations because
# each iteration freezes ≥1 transition or saturates ≥1 resource.
function _fill_tier!(ws::AllocWorkspace, supply, weights, fmax, tier)
    isempty(tier) && return ws.f
    f = ws.f
    r = ws.r
    active = ws.active
    D = ws.D
    nS = size(ws.req, 1)

    copyto!(r, supply)
    fill!(active, false)

    # Demand-free transitions (positive cap, no positive resource demand) fill immediately and
    # consume nothing; others with a positive cap and positive weight start active.
    for t in tier
        has_demand = any(s -> ws.req[s, t] > 0, 1:nS)
        if fmax[t] > 0 && !has_demand
            f[t] = fmax[t]
        elseif fmax[t] > f[t] && weights[t] > 0 && has_demand
            active[t] = true
        end
    end

    tol = 1.0e-12
    while any(active)
        # 1. largest fill-step τ that neither overdraws a resource nor overshoots a cap.
        dτ = Inf
        fill!(D, 0.0)
        for t in tier
            active[t] || continue
            for s in 1:nS
                D[s] += weights[t] * ws.req[s, t]
            end
        end
        for s in 1:nS
            D[s] > 0 && (dτ = min(dτ, r[s] / D[s]))
        end
        for t in tier
            active[t] && (dτ = min(dτ, (fmax[t] - f[t]) / weights[t]))
        end
        isfinite(dτ) || break

        # 2. advance all active fills and debit resources.
        for t in tier
            active[t] && (f[t] += weights[t] * dτ)
        end
        for s in 1:nS
            r[s] = max(0.0, r[s] - D[s] * dτ)
        end

        # 3. freeze transitions that hit their cap, or that need a now-saturated resource.
        for t in tier
            active[t] && f[t] >= fmax[t] - tol && (active[t] = false; f[t] = fmax[t])
        end
        for s in 1:nS
            if r[s] <= tol
                for t in tier
                    active[t] && ws.req[s, t] > 0 && (active[t] = false)
                end
            end
        end
    end

    return f
end

"""
    spawn_integer!(ws, u, w, q_desired) -> Vector{Int}

Integer spawn wrapper (ADR 0002 "Spawn phase"). Runs `progressive_fill!` with `fmax = q_desired`
to get the real-valued fair fill `f`, floors it to whole instances `n = floor.(f)`, then does a
**deterministic priority-ordered integer top-up**: with the supply left after the floored grants,
iterate transitions in `(-priority, index)` order (priority descending, index ascending as the
stable tie-break) and, while a whole additional instance fits, grant it; repeat until none fits.

Keeps spawn counts integral while staying work-conserving and deterministic. `ws.req` must already
hold the spawn requirements (call `build_requirements!` first). Replaces `get_init_satisfied`.
"""
function spawn_integer!(ws::AllocWorkspace, u, w, q_desired)
    nT = length(w)
    fmax = Float64.(q_desired)
    progressive_fill!(ws, u, w; fmax = fmax)

    n = floor.(Int, ws.f)

    # Residual supply after the floored grants.
    r = ws.r
    copyto!(r, u)
    for t in 1:nT, s in axes(ws.req, 1)
        r[s] -= ws.req[s, t] * n[t]
    end

    # Deterministic top-up: (-priority, index) order, grant whole instances while they fit.
    order = sort(1:nT; by = t -> (-w[t], t))
    progress = true
    while progress
        progress = false
        for t in order
            n[t] >= q_desired[t] && continue
            fits = all(s -> r[s] >= ws.req[s, t] - 1.0e-12, axes(ws.req, 1))
            fits || continue
            for s in axes(ws.req, 1)
                r[s] -= ws.req[s, t]
            end
            n[t] += 1
            progress = true
        end
    end

    return n
end

"""
Evolve transitions, spawn new transitions.
"""
function evolve!(state)
    actual_allocs = zero(state.u)

    ## schedule new transitions
    qs = zeros(nrows(state, :T))

    foreach(
        i -> qs[i] = state[i, :transRate] * state[i, :transMultiplier],
        row_ids(state, :T),
    )
    # Integerize the spawn-count proposal by FLOORING, not `ceil` (CONTRACT §2.3, §2.8). On the
    # Poisson genesis path `transRate` is already a whole `rand(Poisson(dt·rate))` draw (create.jl:153),
    # so `floor` is a no-op there and spawning stays dt-invariant in expectation. The ONLY source of a
    # fractional `qs` is the `@deterministic` bare-count path: a fractional count (e.g. 0.3) previously
    # got `ceil`'d up to 1 EVERY tick, so halving `dt` (doubling the tick count) roughly doubled the
    # spawned total — an upward bias as dt→0. Flooring truncates a fractional deterministic count to
    # whole instances (0.3 → 0), removing that discretization hazard; integer `@deterministic` counts
    # (the only in-contract form, §2.8) are unaffected. This is the "condition ceil to the poisson path"
    # fix given there is no genesis-mode tag to branch on.
    qs .= floor.(Ref(Int), qs)
    # A transition gated off this tick (deactivated or guard false, ADR 0010 §B) proposes no
    # new instances — zero its genesis quantity so it never competes for resources.
    foreach(i -> state.transitions[:transFiring][i] || (qs[i] = 0), row_ids(state, :T))

    for i in row_ids(state, :T)
        new_instances = qs[i] + state[i, :transToSpawn]
        capacity =
            state[i, :transCapacity] -
            count(t -> t[:transHash] == state[i, :transHash], state.ongoing_transitions)
        (capacity < new_instances) &&
            add_to_spawn!(state, state[i, :transHash], new_instances - capacity)
        qs[i] = min(capacity, new_instances)
    end

    # Spawn allocation (ADR 0002). Build PER-INSTANCE upfront requirements (qs=1) so a unit of
    # fill is one instance, then `spawn_integer!` returns whole instance counts under priority-
    # weighted progressive filling + a deterministic integer top-up. `allocs[s,t]` is the realized
    # demand (per-instance req × granted count) consumed below and read by the structured bind loop.
    ws = AllocWorkspace(nrows(state, :S), nrows(state, :T))
    build_requirements!(ws, state, ones(nrows(state, :T)); ongoing = false)
    n = spawn_integer!(ws, state.u, state[:, :transPriority], qs)
    allocs = ws.req .* reshape(Float64.(n), 1, :)
    qs .= n

    push!(
        state.log,
        (
            :new_transitions,
            state.t,
            [(hash, q) for (hash, q) in zip(state[:, :transHash], qs)]...,
        ),
    )
    state.u .-= sum(allocs; dims = 2)
    actual_allocs .+= sum(allocs; dims = 2)

    # MVP finding D — per-program ledger: snapshot the spawn-phase per-transition consumption
    # BEFORE the bind loop below mutates `allocs[j,i]` (it decrements it per bound token). Each
    # spawned transition i consumes `spawn_allocs[:, i]`; we attribute that cost to the tokens it
    # binds (src/ledger.jl::attribute_cost!), once the bind list is populated, just below.
    spawn_allocs = copy(allocs)

    structured_token = collect(values(inners(getagent(state, "structured"))))

    # add spawned transitions to the heap
    for i in row_ids(state, :T)
        if qs[i] != 0
            transition = Transition(
                string(state[i, :transName]) * "_@$(state.t)",
                i,
                get_sampled_transition(state, i),
                AbstractAlgebraicAgent[],
                AbstractAlgebraicAgent[],
                [],
                state.t,
                qs[i],
                0.0,
            )
            push!(state.ongoing_transitions, transition)

            bound = transition.bound_tokens
            binding = transition.binding

            for (j, type) in enumerate(state.network[:, :placeName])
                if type ∈ state.structured_token
                    if !isinteger(allocs[j, i])
                        error(
                            "For structured place, multiplicity coefficient must be integer in transition $i.",
                        )
                    end

                    # ADR 0008 §B: narrow the candidate set by the LHS arc's predicate
                    # (kind-only when none), then the unchanged priority sort + integer take.
                    pred = lhs_predicate(state.transitions[:transLHS][i], type)
                    available_places = filter(
                        a ->
                        get_place(a) == type &&
                            !isblocked(a) &&
                            matches(pred, a, state, transition),
                        structured_token,
                    )

                    # Total order (ADR 0008 inv 3): highest priority first, ties broken by the
                    # deterministic (place, creation_index) key — NOT the AA Dict / random-name
                    # order, which would make WHICH equal-priority token binds non-reproducible.
                    sort!(
                        available_places;
                        by = a -> (
                            -priority(a, state.network[i, :transName]),
                            token_sortkey(state, a),
                        ),
                    )

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(available_places)
                        set_bound_transition!(available_places[ix], transition)

                        push!(bound, available_places[ix])
                        push!(binding, type => available_places[ix])
                        add_to_log!(available_places[ix], type, state.t, transition)

                        allocs[j, i] -= 1
                        ix += 1
                    end
                end
            end

            # MVP finding D — attribute this spawned transition's upfront resource cost
            # (spawn_allocs[:, i]) to the program(s) it just bound (src/ledger.jl). Done here,
            # AFTER the bind loop, so `transition.bound_tokens` is populated.
            attribute_cost!(state, transition, @view spawn_allocs[:, i])

            context_eval(
                state,
                transition,
                state.wrap_fun(state.network[i, :transPreAction]),
            )
        end
    end

    ## evolve ongoing transitions
    # Ongoing allocation (ADR 0002). `req` is the full per-tick demand of each in-flight instance
    # group (instance count × multiplicity; `:rate` scaled by dt when transCycleTime>0, `:nonblock`
    # unscaled). Each group fills at most fraction 1.0 of its requested progress this tick
    # (`fmax = 1`), so the fill fraction `f[i]` ∈ [0,1] IS the saturation `qs[i]` and progress
    # advances by `qs[i]*dt`. Priority is re-read FRESH per tick from the recipe row `t.i`
    # (`state[t.i, :transPriority]` → per-tick context_eval), NOT the spawn-time snapshot — so a
    # time-varying transPriority applies to in-flight instances too (ADR 0002 "fresh per tick").
    nong = length(state.ongoing_transitions)
    ws = AllocWorkspace(nrows(state, :S), nong)
    qs = map(t -> t.q, state.ongoing_transitions)
    build_requirements!(ws, state, qs; ongoing = true, dt_scale = state.dt)
    w = [state[t.i, :transPriority] for t in state.ongoing_transitions]
    progressive_fill!(ws, state.u, w; fmax = fill(1.0, nong))
    qs = copy(ws.f)
    allocs = ws.req .* reshape(qs, 1, :)
    push!(
        state.log,
        (
            :saturation,
            state.t,
            [
                (state.ongoing_transitions[i][:transHash], qs[i]) for
                    i in eachindex(state.ongoing_transitions)
            ]...,
        ),
    )
    state.u .-= sum(allocs; dims = 2)
    actual_allocs .+= sum(allocs; dims = 2)

    # MVP finding D — snapshot the ongoing-phase per-transition consumption (mostly the @rate burn,
    # e.g. the BD demo's per-tick `budget` spend) BEFORE the nonblock-bind loop mutates allocs.
    ongoing_allocs = copy(allocs)

    for i in eachindex(state.ongoing_transitions)
        transition = state.ongoing_transitions[i]
        if qs[i] != 0
            transition.state += qs[i] * state.dt

            bound = transition.nonblock_tokens
            binding = transition.binding

            for (j, type) in enumerate(state.network[:, :placeName])
                if type ∈ state.structured_token
                    if !isinteger(allocs[j, i])
                        error(
                            "For structured place, multiplicity coefficient must be integer in transition $i.",
                        )
                    end

                    # ADR 0008 §B: narrow by the in-flight transition's LHS predicate.
                    pred = lhs_predicate(transition[:transLHS], type)
                    available_places = filter(
                        a ->
                        get_place(a) == type &&
                            !isblocked(a) &&
                            matches(pred, a, state, transition),
                        structured_token,
                    )

                    # Total order (ADR 0008 inv 3): highest priority first, ties broken by the
                    # deterministic (place, creation_index) key. NB use `transition.i` (the recipe
                    # index stored at spawn), NOT the loop var `i` — here `i` indexes the
                    # ongoing_transitions array, not the :T schema row, so `state.network[i, …]` would
                    # read the wrong transition's priority (latent: harmless only while priority is
                    # the default 0.0 for all tokens; a per-transition priority override would hit it).
                    sort!(
                        available_places;
                        by = a -> (
                            -priority(a, state.network[transition.i, :transName]),
                            token_sortkey(state, a),
                        ),
                    )

                    ix = 1
                    while allocs[j, i] > 0 && ix <= length(available_places)
                        set_bound_transition!(available_places[ix], transition)

                        push!(bound, available_places[ix])
                        push!(binding, type => available_places[ix])
                        add_to_log!(available_places[ix], type, state.t, transition)

                        allocs[j, i] -= 1
                        ix += 1
                    end
                end
            end

            # MVP finding D — attribute this in-flight transition's per-tick (rate/nonblock)
            # resource cost to its bound program(s). `bound_tokens` (the @select'd token,
            # bound at spawn) plus any nonblock tokens just bound are charged; an instance with no
            # bound program books its burn against the unattributed bucket (src/ledger.jl).
            attribute_cost!(state, transition, @view ongoing_allocs[:, i])
        end
    end

    push!(state.log, (:allocation, state.t, actual_allocs))
    return push!(
        state.log,
        (
            :valuation_cost,
            state.t,
            actual_allocs' * [state[i, :placeCost] for i in row_ids(state, :S)],
        ),
    )
end

# The legacy `event_action!` (a no-op fetch of :eventAction, the repaired CONTRACT §3.4 Inv 7
# defect) is superseded by the endogenous decision channel `fire_rules!` (ADR 0010, src/actions.jl),
# which evaluates each Rule's guard and runs its action at _step! step 10. :E rows are lifted to
# Rules at construction.

function allocate_for_move(t::Transition, s::Symbol)
    return t.bound_tokens ∩
        map(x -> x[2], filter(x -> x[1] == s, t.binding))
end

function structured_rhs(expr::Expr, state, transition)
    if isexpr(expr, :macrocall) && macroname(expr) == :structured
        if length(expr.args) >= 3 && expr.args[3] isa QuoteNode
            # NAMED form: `@structured(:Kind, field = node, …)` — the eval-free-serializable genesis
            # product (ADR 0005 §39: the structured escape-hatch promoted to a typed node). It is the
            # RHS-product twin of the `AddToken` rule action (actions.jl:203) and shares its host
            # contract EXACTLY: the document/reaction-line carries the registry KEY plus field-value
            # nodes, never the constructor — the host `state.registry[:Kind]` (a `(state, fields::Dict)
            # -> token` ctor, ADR 0006 §C) is resolved BY NAME at firing time. Field values evaluate
            # through the same seeded closure as everywhere else (`_eval_value`), so `@t()` / a
            # `rand(state.rng, …)` draw are reproducible under the run's seed. Discriminated from the
            # raw forms below by `args[3] isa QuoteNode` — a bare quoted kind symbol is never a valid
            # raw token body (a Symbol is not a token), so the two never collide.
            kind = expr.args[3].value
            haskey(state.registry, kind) ||
                error("@structured: kind $kind not in the network registry (ADR 0006 §C)")
            ctor = state.registry[kind]
            fieldvals = Dict{Symbol, Any}()
            for kw in @view expr.args[4:end]
                (isexpr(kw, :(=)) && kw.args[1] isa Symbol) || error(
                    "@structured($kind, …): each field must be `name = value`, got `$(kw)`",
                )
                fieldvals[kw.args[1]] = _eval_value(state, transition, kw.args[2])
            end
            token = ctor(state, fieldvals)
            entangle!(getagent(state, "structured"), token)
            return token, get_place(token)
        else
            # The raw `@structured(Ctor(…))` / `@structured(token, place)` forms were removed —
            # the named, registry-resolved form above is the only supported genesis product (it is
            # the sole one that serializes eval-free; ADR 0005 §39 / 0006 §C). Construction rejects
            # a raw line (recursively_find_arcs!, create.jl), so reaching here means a
            # hand-built :trans Expr bypassed that check — surface it rather than eval host code.
            error(
                "@structured: only the named form `@structured(:Kind, field = value, …)` is " *
                    "supported; the raw constructor form was removed (see create.jl). Got: $expr",
            )
        end
    elseif isexpr(expr, :macrocall) && macroname(expr) == :move
        expr = quote
            place_from = $(expr.args[end - 1])
            place_to = $(expr.args[end])

            return place_from, place_to
        end

        place_from, place_to =
            Symbol.(context_eval(state, transition, state.wrap_fun(expr)))

        tokens =
            filter(x -> get_place(x) == place_from, transition.bound_tokens)

        if !isempty(tokens)
            token = first(tokens)
            entangle!(getagent(state, "structured"), token)

            set_place!(token, place_to)
            ix = findfirst(
                i -> transition.bound_tokens[i] == token,
                eachindex(transition.bound_tokens),
            )
            deleteat!(transition.bound_tokens, ix)
            set_bound_transition!(token, nothing)

            return token, place_to
        else
            # No bound token of place_from to move — a graceful no-op (finish! skips a nothing
            # place), consistent with @advance; do NOT fall through to an implicit nothing that
            # would crash the (token, place) unpack at the call site.
            @error "Not enough tokens to allocate for a move."
            return nothing, nothing
        end

    elseif isexpr(expr, :macrocall) && macroname(expr) == :advance
        # @advance(field, value): advance a bound token's lifecycle by writing one field, keeping
        # its identity/kind/uuid/creation_index/past_bonds (ADR 0008 §D). The phase-as-attribute
        # generalization of @move (which writes the `place` field). `value` may read the token's
        # own current fields via @field(name), and MAY draw (§F). The advanced token is then
        # released. It is consumed from the first bound token of this transition.
        field = expr.args[3]
        field isa Symbol ||
            error("@advance: first argument must be a field name, got $field")
        valex = expr.args[4]
        # No bound token to advance (the @select predicate matched nothing this firing) — a
        # silent no-op: the instance produced no advance. finish! skips the nothing return.
        isempty(transition.bound_tokens) && return nothing, nothing
        token = first(transition.bound_tokens)
        # Evaluate the value with the bound token in scope so @field(name) reads its attributes.
        val = eval_with_token(state, transition, token, valex)
        # `:species` is the retired ADR-0017 spelling of `:place`, accepted (silently — this is
        # inside the step loop) for one release.
        if field === :place || field === :species
            set_place!(token, Symbol(val))
        else
            setproperty!(token, field, val)
        end
        deleteat!(transition.bound_tokens, 1)
        set_bound_transition!(token, nothing)
        return token, get_place(token)

    else
        token = context_eval(state, transition, state.wrap_fun(expr))
        entangle!(getagent(state, "structured"), token)

        return token, get_place(token)
    end
end

# collect terminated transitions
function finish!(state)
    val_reward = 0
    terminated_all = Dict{Symbol, Float64}()
    terminated_success = Dict{Symbol, Float64}()

    ix = 1
    while ix <= length(state.ongoing_transitions)
        trans_ = state.ongoing_transitions[ix]
        ((state.t - trans_.t) < trans_.trans[:transMaxLifeTime]) &&
            (trans_.state < trans_[:transCycleTime]) &&
            (ix += 1; continue)

        q = if trans_.state >= trans_[:transCycleTime]
            rand(state.rng, Distributions.Binomial(Int(trans_.q), trans_[:transProbOfSuccess]))
        else
            0
        end

        # MVP finding D — per-program ledger: snapshot the program(s) bound to this finishing
        # transition (an @advance/@move RHS op moves the bound token OUT of bound_tokens
        # during the loop below, so we must capture them first) and the running reward BEFORE its
        # RHS emission, to attribute this transition's realized reward to its program(s) afterward.
        reward_before = val_reward
        finishing_tokens = _all_bound_tokens(trans_)

        for r in extract_arcs(trans_[:transRHS], state)
            if r.place isa Expr
                multiplicity = context_eval(state, trans_, state.wrap_fun(r.multiplicity))

                for _ in 1:(q * multiplicity)
                    token, place = structured_rhs(r.place, state, trans_)
                    # A structured-RHS op may legitimately produce nothing (e.g. @advance with no
                    # bound token to advance) — skip the count/reward in that case.
                    place === nothing && continue
                    i = find_index(place, state)
                    state.u[i] += 1
                    val_reward += state[i, :placeReward]
                end
            else
                i = find_index(r.place, state)
                multiplicity = context_eval(state, trans_, state.wrap_fun(r.multiplicity))

                state.u[i] += q * multiplicity
                val_reward += state[i, :placeReward] * q * multiplicity
            end
        end

        # MVP finding D — attribute this transition's realized reward (the delta it just emitted) to
        # the program(s) it was bound to, split evenly; an unbound (plain) transition's reward goes to
        # the unattributed bucket (src/ledger.jl). Use the pre-RHS snapshot so an advanced/moved
        # program (already removed from bound_tokens) still receives its reward.
        attribute_reward!(state, trans_, finishing_tokens, val_reward - reward_before)

        for tok in trans_[:transLHS]
            if in(:conserved, tok.modality)
                state.u[tok.index] +=
                    trans_.q *
                    tok.multiplicity *
                    (in(:rate, tok.modality) ? trans_[:transCycleTime] : 1)
                if tok.place ∈ state.structured_token
                    for _ in 1:(trans_.q * tok.multiplicity)
                        agent_ix = findfirst(
                            a -> get_place(a) == tok.place,
                            trans_.bound_tokens,
                        )
                        # No more bound tokens of this place to release (a multi-place
                        # transition may exhaust one place before the q*multiplicity count) — stop.
                        isnothing(agent_ix) && break

                        set_bound_transition!(
                            trans_.bound_tokens[agent_ix],
                            nothing,
                        )
                        deleteat!(trans_.bound_tokens, agent_ix)
                    end
                end
            end

            if in(:nonblock, tok.modality)
                if in(:conserved, tok.modality)
                    error(
                        "Modalities `:conserved` and `:nonblock` cannot be specified at the same time.",
                    )
                end

                state.u[tok.index] += trans_.q * tok.multiplicity
                if tok.place ∈ state.structured_token
                    for _ in 1:(trans_.q * tok.multiplicity)
                        agent_ix = findfirst(
                            a -> get_place(a) == tok.place,
                            trans_.nonblock_tokens,
                        )
                        isnothing(agent_ix) && break   # no more nonblock tokens of this place

                        set_bound_transition!(
                            trans_.nonblock_tokens[agent_ix],
                            nothing,
                        )
                        deleteat!(trans_.nonblock_tokens, agent_ix)
                    end
                end
            end
        end

        context_eval(
            state,
            trans_,
            state.wrap_fun(state.network[trans_.i, :transPostAction]),
        )

        for agent in trans_.bound_tokens
            set_place!(agent, :removed)
            set_bound_transition!(agent, nothing)
        end

        terminated_all[Symbol(trans_[:transHash])] =
            get(terminated_all, Symbol(trans_[:transHash]), 0) + trans_.q

        terminated_success[Symbol(trans_[:transHash])] =
            get(terminated_success, Symbol(trans_[:transHash]), 0) + q

        ix += 1
    end

    # Prune every instance that passed the terminal test above — i.e. keep only those that
    # have neither completed their cycle NOR aged out. This must mirror the skip condition at
    # the top of the loop; the old predicate ignored max-lifetime, so timed-out instances
    # (state < cycleTime but age ≥ maxLifeTime) were retained and re-emitted/​re-credited every
    # subsequent tick (violating conservation + termination-completeness, §3.4 INV2/INV6).
    filter!(
        s -> ((state.t - s.t) < s[:transMaxLifeTime]) && (s.state < s[:transCycleTime]),
        state.ongoing_transitions,
    )

    push!(state.log, (:terminated_all, state.t, terminated_all...))
    push!(state.log, (:terminated_success, state.t, terminated_success...))
    push!(state.log, (:valuation_reward, state.t, val_reward))

    return state.u
end

function free_blocked_places!(state)
    for trans in state.ongoing_transitions, tok in trans[:transLHS]
        in(:nonblock, tok.modality) && (state.u[tok.index] += trans.q * tok.multiplicity)
    end

    for trans in state.ongoing_transitions
        for a in trans.nonblock_tokens
            a.bound_transition = nothing
        end

        empty!(trans.nonblock_tokens)
    end
    return
end

## resolve tspan, dt

function get_tcontrol(tspan, args)
    tspan isa Tuple && (tspan = tspan[2] - tspan[1])
    tunit = get(args, :tunit, oneunit(tspan))
    tspan = tspan / tunit

    dt = get(args, :dt, haskey(args, :tstops) ? tspan / args[:tstops] : tunit) / tunit

    return ((0.0, tspan), dt)
end

# ── CONTRACT §1.4 — construction-time modality validation ─────────────────────────────────────
# The three orthogonal modality axes (§1.1) admit exactly five legal rows (§1.3); three cross-
# products are illegal (§1.4) and, left unchecked, fail LATE or SILENTLY deep in the stepper. We
# reject all three HERE — before any closure compiles or any tick runs — with a clear ArgumentError
# naming the transition, the offending LHS token, and the §1.4 rule. The corresponding deep-path
# errors are thereby UNREACHABLE for these cases (the construction check fires first) but are left in
# place as defensive belt-and-suspenders:
#   1. {:nonblock, :conserved}    — else errors in finish! (solvers.jl ~735) / crashes on `q` in
#                                    free_blocked_places! on the 2nd tick.
#   2. :rate (perstep) with C==0  — else constructs and runs SILENTLY (build_requirements! gates the
#                                    per-step draw on C>0, ~117, so the token never meters).
#   3. :rate (perstep) on a       — else errors deep in build_requirements! (~114).
#      structured/agent place
#
# Arcs are read via the eval-free static decomposition (`_split_reaction_line` +
# `_static_arcs`, serialize.jl) — the SAME parse the runtime/exporter use — so the checked
# modality Set matches what the engine forms per tick. The effective per-token modality unions the
# arc's wrapper tags with the place's `:placeDefaultModality` (the `@mode` channel), exactly as the
# runtime does at state.jl:309. Lines the static splitter cannot handle (`@choose`/bidirectional)
# are the escape hatch and are left un-validated (they are un-validatable statically).
function validate_modalities(net::ReactionNetwork)
    for t in row_ids(net, :T)
        lhs = try
            first(_split_reaction_line(net[t, :trans]))
        catch
            continue    # @choose / bidirectional / non-standard line: escape hatch — skip
        end
        tname = net[t, :transName]
        tlabel = (tname === missing || tname === nothing) ? "t$t" : string(tname)
        ct = net[t, :transCycleTime]
        for r in _static_arcs(lhs)
            sname = string(r.place)
            i = r.place isa Symbol ? find_index(r.place, net) : nothing
            mod = i === nothing ? r.modality : (r.modality ∪ net[i, :placeDefaultModality])

            # Rule 1 — blocking = nonblock requires return = consumed.
            if in(:nonblock, mod) && in(:conserved, mod)
                throw(
                    ArgumentError(
                        "Transition `$tlabel`, LHS token `$sname`: modality {:nonblock, :conserved} is " *
                            "illegal (CONTRACT §1.4) — a resource cannot be both held-until-finish " *
                            "(:conserved) and released-every-step (:nonblock). `blocking = nonblock` " *
                            "requires `return = consumed`.",
                    ),
                )
            end

            # Rule 2 — allocation = perstep requires transCycleTime > 0. Only a concrete numeric
            # cycletime is statically checkable; an Expr/param-valued cycletime is left to run.
            if in(:rate, mod) && ct isa Real && iszero(ct)
                throw(
                    ArgumentError(
                        "Transition `$tlabel`, LHS token `$sname`: modality :rate (perstep) with " *
                            "cycletime == 0 is illegal (CONTRACT §1.4) — a per-step reservation only fires " *
                            "when cycletime > 0, so with cycletime == 0 it silently reserves nothing. " *
                            "`allocation = perstep` requires `transCycleTime > 0`.",
                    ),
                )
            end

            # Rule 3 — allocation = perstep requires a non-structured (countable) place.
            if in(:rate, mod) && i !== nothing && net[i, :placeStructured] === true
                throw(
                    ArgumentError(
                        "Transition `$tlabel`, LHS token `$sname`: modality :rate (perstep) on a " *
                            "structured/agent place is illegal (CONTRACT §1.4) — you cannot reserve a " *
                            "fractional, dt-scaled slice of an indivisible agent. `allocation = perstep` " *
                            "requires a non-structured (countable) place.",
                    ),
                )
            end
        end
    end
    return net
end

"""
    ReactionNetworkProblem(net::ReactionNetwork, u0 = Dict(), p = Dict(); name = "reaction_network", seed = nothing, tspan, dt = 1, kwargs...)

Construct a live simulation state (`ReactionNetworkProblem`) from a static authoring/IR store `net` — the central entry point that turns an authored `@reaction_network` into a runnable, steppable AA node. `u0` overrides plain-place initial markings by name (defaulting to each place's `placeInitVal`); `p` supplies/overrides parameters (merged over the store's declared params); `name` is the agent name. Meta keywords declared in the store (e.g. `tspan`, `dt`, `tunit`) are read as defaults and may be overridden by the matching kwargs. The constructor validates modalities up front (CONTRACT §1.4), compiles the attribute/transition closures against the frozen store positions (ADR 0004), builds the `rules`/`registry` endogenous-decision channel, and instantiates the declarative initial token population before arming the live phase guard.

The `seed` kwarg owns the per-run RNG (CONTRACT §4): it fixes the state-owned stream so a run is fully determined by `(model, seed)`; absent, a fresh seed is drawn from system entropy and the REALIZED value stored on `.seed`, so any run stays replayable. `initial_rng` snapshots the stream at t=0 for `_reinit!`.
"""
function ReactionNetworkProblem(
        net::ReactionNetwork,
        u0 = Dict(),
        p = Dict();
        name = "reaction_network",
        kwargs...,
    )
    assign_defaults!(net)
    # CONTRACT §1.4: reject the three illegal modality configurations up front, before any closure
    # compiles or any tick runs (the T2 acceptance tests require the throw from the constructor
    # itself, not deep in the stepper). Runs after assign_defaults! so :placeDefaultModality is materialized.
    validate_modalities(net)
    keywords = Dict{Symbol, Any}(
        [
            net[i, :metaKeyword] => net[i, :metaVal] for i in row_ids(net, :M) if
                !isnothing(net[i, :metaKeyword]) && !isnothing(net[i, :metaVal])
        ]
    )

    merge!(keywords, Dict(collect(kwargs)))
    # `alloc_strategy` (legacy :weighted/:greedy switch) is accepted and IGNORED: ADR 0002 makes
    # the priority-weighted progressive-filling allocator the single policy (strict greedy is its
    # priority→∞ limit, not a separate path), so the allocator no longer reads `state.p[:strategy]`.
    # A JSON/DSL model may still carry the key; it is a harmless no-op rather than a hard error.

    # Deprecation shim (ADR 0015 follow-up dt/tstep naming pass): the integrator step is authored
    # as `dt` (the settled struct-field name); the legacy `tstep` meta keyword is honored for one
    # release with a depwarn, following the ADR 0015 Tier 1 alias pattern. Map it BEFORE
    # get_tcontrol resolves the step — get_tcontrol reads only :dt/:tstops/:tunit, so a user-supplied
    # `tstep` was previously SILENTLY IGNORED (immediately overwritten by the resolved value below);
    # mapping it here both renames it and fixes that latent silent-ignore bug.
    if haskey(keywords, :tstep) && !haskey(keywords, :dt)
        Base.depwarn(
            "The `tstep` meta keyword is deprecated (dt/tstep naming pass); use `dt` instead.",
            :ReactionNetworkProblem,
        )
        keywords[:dt] = pop!(keywords, :tstep)
    end

    keywords[:tspan], keywords[:dt] = get_tcontrol(keywords[:tspan], keywords)

    # Determinism (§4 D5/D6): build the state-owned RNG. A `seed` kwarg fixes the stream;
    # absent it, draw a fresh seed from system entropy so a default run is still self-contained.
    # The REALIZED seed is stored on the struct (so an entropy-seeded run is replayable, D6)
    # and `initial_rng` snapshots the stream at t=0 so `_reinit!` restores it exactly (D7).
    # Any `Integer` seed is accepted verbatim (e.g. a `hash((root, k))` ensemble member key, §4 D8).
    seed = get(keywords, :seed, nothing)
    seed = isnothing(seed) ? rand(Random.RandomDevice(), UInt64) : seed
    rng = Random.Xoshiro(seed)
    initial_rng = copy(rng)

    net = register_observables(net)

    structured_token_names =
        net[filter(i -> net[i, :placeStructured], 1:nrows(net, :S)), :placeName]

    attrs, transitions, wrap_fun = compile_attrs(net, structured_token_names)
    transition_recipes = transitions
    u0_init = zeros(nrows(net, :S))

    for i in row_ids(net, :S)
        if !isnothing(net[i, :placeName]) && haskey(u0, net[i, :placeName])
            u0_init[i] = u0[net[i, :placeName]]
        else
            u0_init[i] = net[i, :placeInitVal]
        end
    end

    prms = Dict{Symbol, Any}(
        (
            net[i, :prmName] => net[i, :prmVal] for
                i in Iterators.filter(i -> !isnothing(net[i, :prmVal]), 1:nrows(net, :P))
        )
    )

    merge!(p, prms)

    ongoing_transitions = Transition[]
    log = NamedTuple[]
    observables = compile_observables(net)
    transitions_attrs =
        setdiff(
        filter(a -> contains(string(a), "trans"), propertynames(net.columns)),
        (:trans,),
    ) ∪ [:transLHS, :transRHS, :transToSpawn, :transHash, :transFiring]
    transitions = Dict{Symbol, Vector}(a => [] for a in transitions_attrs)

    sol = DataFrame(
        "t" => Float64[],
        (string(name) => Float64[] for name in net[:, :placeName])...,
    )

    # Endogenous decision channel (ADR 0010 §12). Per-network host registry for AddToken/Invoke
    # (ADR 0006 §C — by-name, never eval'd). Rules are built from the :E rows: a legacy event
    # `trigger && action` becomes a Rule{guard=trigger, action=RawExpr(action), every_tick}.
    # Typed Rules can also be supplied directly via the `rules=` kwarg / @rule authoring.
    registry = Dict{Symbol, Any}(get(keywords, :registry, Dict{Symbol, Any}()))
    rules = Any[
        Rule(Symbol("rule_", i), net[i, :eventTrigger], RawExpr(net[i, :eventAction]))
            for i in row_ids(net, :E) if
            !isnothing(net[i, :eventTrigger]) && !isnothing(net[i, :eventAction])
    ]
    append!(rules, get(keywords, :rules, Any[]))

    # External coupling (ADR 0012 §B). Seed the per-tick input buffer from the declared `inputs[]`
    # defaults (passed by from_json_model as `external_inputs=`), so a port read before any AA wire
    # delivers — or in a standalone wire-less run — has a well-defined fallback (§B3). The buffer is
    # re-derived every `_prestep!` (merged over a copy of the defaults); the defaults snapshot is
    # kept immutable so `_reinit!` can restore the pre-wire seed (§4 D7).
    external_input_defaults =
        Dict{Symbol, Any}(get(keywords, :external_inputs, Dict{Symbol, Any}()))
    external_inputs = copy(external_input_defaults)

    network = ReactionNetworkProblem(
        name,
        net,
        attrs,
        transition_recipes,
        u0_init,
        p,
        keywords[:tspan][1],
        structured_token_names,
        keywords[:tspan],
        get(keywords, :dt, 1),
        transitions,
        ongoing_transitions,
        log,
        observables,
        wrap_fun,
        sol,
        rng,
        seed,
        initial_rng,
        rules,
        registry,
        Dict{Symbol, Int}(),
        Dict{String, Int}(),
        collect(get(keywords, :population, [])),
        Dict{String, Dict{Symbol, Any}}(),
        false,
        # Per-program ledger (MVP finding D, src/ledger.jl): empty at construction, populated at the
        # bind/finish sites in evolve!/finish!, reset by _reinit!.
        Dict{String, ProgramLedger}(),
        0.0,
        0.0,
        external_inputs,
        external_input_defaults,
        # Per-token trajectory log (ADR 0013 §14.1): empty at construction, appended each tick by
        # push_token_trajectory_row! for opted-in kinds, reset by _reinit!.
        Tuple{Float64, String, Symbol, NamedTuple}[],
    )

    entangle!(network, FreeAgent("structured"))

    # Instantiate the declarative initial marking (ADR 0007 §B) into the structured container,
    # in declared order, with seeded attribute draws + creation indices — BEFORE t=0/the first
    # step, so a structured run is reproducible from (model, seed, population). Then arm the §A
    # phase guard: the model is now Live and reindexers (rem_parts!) must refuse.
    instantiate_population!(network)
    snapshot_population!(network)
    update_u_structured!(network)
    network.live = true

    # save!(network)

    return network
end

function AlgebraicAgents._reinit!(state::ReactionNetworkProblem; seed = nothing)
    state.u .= isempty(state.sol) ? state.u : Vector(state.sol[1, 2:end])
    state.t = state.tspan[1]
    empty!(state.ongoing_transitions)
    empty!(state.log)
    state.observables = compile_observables(state.network)
    empty!(state.sol)
    # RNG restore (§4 D7) — OR reseed (ensemble mode b, ADR 0013 §14.2). With NO `seed` (the default
    # AA `reinit!(a)` path, byte-identical to before) restore the construction stream so the second
    # run reproduces the first. With a `seed`, INSTALL that seed's stream as the NEW construction
    # stream — mirroring the constructor (`:845-848`): `Xoshiro(seed)`, snapshot `initial_rng`, store
    # the realized seed. This makes the reseeded state equivalent to a fresh `build(seed)`. CRITICAL
    # ORDERING: it MUST precede `instantiate_population!` below — the t=0 marking's `count`+attribute
    # draws are sampled through `state.rng`, so a reseeded member's initial attributes match what a
    # fresh `build(seed)` samples (the mode-(a) ≡ mode-(b) equivalence contract, §14.2).
    if seed === nothing
        state.rng = copy(state.initial_rng)
    else
        state.rng = Random.Xoshiro(seed)
        state.initial_rng = copy(state.rng)
        state.seed = seed
    end
    # Reset every `once` rule's latch so a re-run from the same seed reproduces the lever (§4 D7).
    for r in state.rules
        r.fire_mode === :once && (r.enabled = true)
    end
    # Tear down the live token population and rebuild the declarative initial marking (ADR 0007
    # §D): without this the END-state tokens survive into the next run. Reset the creation
    # counters, then re-instantiate population[] in declared order through the restored RNG — so
    # `init → step* → reinit! → step*` reproduces the first trajectory for structured models too
    # (closing §4 D7 for structured runs). NB the explicit-host-token population form holds the
    # SAME agent objects; re-entangling them after a run re-uses them with reset bonds.
    container = getagent(state, "structured")
    for tok in collect(values(inners(container)))
        disentangle!(tok)
    end
    empty!(state.creation_counters)
    empty!(state.creation_index)
    # Clear the per-program ledger so a re-run from the same seed rebuilds it identically (§4 D7,
    # MVP finding D — mirrors the creation_counters reset above).
    reset_program_ledger!(state)
    # Clear the per-token trajectory log too (ADR 0013 §14.1 — same §4 D7 rebuild contract as the
    # ledger: `init → step* → reinit! → step*` must reproduce the first run's trajectory rows).
    empty!(state.token_trajectory)
    for entry in state.population
        if !(entry isa PopulationEntry)
            set_bound_transition!(entry, nothing)
            empty!(entry.past_bonds)
            # restore the SAME object's attributes (place/phase/…) to their captured t=0 values
            restore_token_snapshot!(state, entry)
        end
    end
    instantiate_population!(state)
    update_u_structured!(state)

    # Drop stale latched external-input wire values and restore the pre-wire defaults (ADR 0012
    # §B3 / §10.4 / §4 D7), so the first post-reinit `_prestep!` re-latches from a clean seed and
    # `init → step* → reinit! → step*` reproduces the first coupled trajectory.
    empty!(state.external_inputs)
    merge!(state.external_inputs, state.external_input_defaults)

    return state
end

# The `reinit!` hierarchy walker (AA `interface.jl:256`) takes NO kwargs — it calls `_reinit!(a)`
# positionally then recurses into `inners`. Ensemble mode (b) needs to forward a `seed` to the
# top-level `_reinit!` (the reseed path above), so we override `reinit!` for our type with a `seed`
# kwarg while otherwise reproducing AA's walk EXACTLY. The `seed = nothing` default makes the plain
# `reinit!(prob)` call byte-identical to AA's (same ops, same order) — it backs the existing green
# reproducibility test — so this override only ADDS the reseed capability, it changes nothing else.
function AlgebraicAgents.reinit!(state::ReactionNetworkProblem; seed = nothing)
    AlgebraicAgents._reinit!(state; seed = seed)
    for a in values(AlgebraicAgents.inners(state))
        AlgebraicAgents.reinit!(a)
    end
    return state
end

function update_u_structured!(state)
    structured_tokens = collect(values(inners(getagent(state, "structured"))))
    for (i, place) in enumerate(state.network[:, :placeName])
        if state.network[i, :placeStructured]
            state.u[i] =
                count(a -> get_place(a) == place && !isblocked(a), structured_tokens)
        end
    end

    return state.u
end

function AlgebraicAgents._step!(state::ReactionNetworkProblem)
    update_u_structured!(state)
    if isempty(state.sol)
        save!(state)
    end

    free_blocked_places!(state)
    update_u_structured!(state)
    update_observables(state)
    sample_transitions!(state)
    evolve!(state)
    update_u_structured!(state)
    finish!(state)
    update_u_structured!(state)

    # Step 10 (§3.3): fire the endogenous decision channel (ADR 0010). Rules see this tick's
    # post-finish state; their writes land on this tick's ledger row and the next tick's
    # genesis/guards. Replaces the old no-op event_action! slot.
    fire_rules!(state)
    update_u_structured!(state)

    push!(
        state.log,
        (
            :valuation,
            state.t,
            state.u' * [state[i, :placeValuation] for i in row_ids(state, :S)],
        ),
    )

    # MVP finding D — mark each live program to market (its place's placeValuation) and push the
    # per-tick per-program ledger row, in deterministic token order, right after the aggregate
    # :valuation row so the per-program and aggregate views are consistent (src/ledger.jl).
    attribute_valuation!(state)
    push_program_ledger_row!(state)

    # ADR 0013 §14.1: snapshot each opted-in token's declared fields into the per-token trajectory
    # log, at the SAME seam as the per-program ledger row (all physics done, tokens in final
    # positions, timestamp still this tick's `t`) and in the SAME deterministic token order — so the
    # trajectory log and the ledger share one observation point and one determinism guarantee.
    push_token_trajectory_row!(state)

    state.t += state.dt

    save!(state)

    return state.t
end

function AlgebraicAgents._projected_to(state::ReactionNetworkProblem)
    return state.t > state.tspan[2] ? true : state.t
end

function fetch_params(net::ReactionNetwork)
    return Dict{Symbol, Any}(
        (
            net[i, :prmName] => net[i, :prmVal] for
                i in Iterators.filter(i -> !isnothing(net[i, :prmVal]), row_ids(net, :P))
        )
    )
end
