# Per-program (per-structured-token) ledger — MVP finding D / D-bis.
#
# ── Why this exists ─────────────────────────────────────────────────────────────────────
# RD's built-in ledger (CONTRACT §8.5) is AGGREGATE only: each tick `evolve!`/`finish!`/`_step!`
# push pool-level rows into `state.log` — `(:valuation_cost, t, scalar)`, `(:valuation_reward, t,
# scalar)`, `(:valuation, t, scalar)`, `(:allocation, t, vector)`. There is NO attribution of
# cost/reward to an individual STRUCTURED TOKEN (a "program"). MVP finding D records exactly this
# gap: "two non-composing accounting systems; no per-tick hook to accrue economics onto an in-flight
# token", forcing the BD acquisition demo to RECONSTRUCT per-program economics in post-processing
# (demo/bd_acquisition/analysis.jl::portfolio_rnpv walks the final token population and reads host
# fields like `npv_peak`/`pos_remaining`, NOT an engine ledger). This file gives the engine a
# first-class per-program ledger, attributed DURING the run, that the demo can read or cross-check.
#
# NB this is DISTINCT from CONTRACT §4.6 (per-entity RNG substreams). That is about RNG determinism
# (which draws shift when an entity is injected); this is about cost/reward ATTRIBUTION. The two do
# not interact: attribution here is eval-free, append-only, and consumes NO randomness.
#
# ── The attribution rule (chosen, with its boundary documented) ─────────────────────────
# A transition's per-tick resource COST is `Σ_s allocs[s] · placeCost[s]` over the species it
# consumed this tick — the SAME quantity the aggregate `:valuation_cost` row sums (solvers.jl
# `evolve!`). We attribute that transition's cost to the structured token(s) BOUND to it this tick,
# split EVENLY across them:
#
#     each bound token's share  =  (transition cost this tick) / (number of bound tokens)
#
# Boundary (documented, deliberate):
#   • A transition with EXACTLY ONE bound structured token (the BD case — genesis is token-gated by
#     @select, so each `@advance` instance binds one Project) receives that transition's FULL
#     consumed cost. This is the simplest correct rule and makes the BD per-program rNPV directly
#     reconstructable: a program's `cost_incurred` is the sum of the burn of every advance it sat in.
#   • A transition with SEVERAL bound tokens splits its cost evenly among them. (Alternative rules —
#     split by per-token stoichiometric weight, or by a token "size" field — are equally defensible;
#     even split is chosen for being the simplest order-independent rule and because the demo's
#     transitions each bind one token, where every split rule coincides.)
#   • A PLAIN-species-only transition (no bound structured token — e.g. the `financing` inflow, or a
#     classic non-agentic reaction) has NO program to attribute to. Its cost is recorded against the
#     network-level UNATTRIBUTED bucket (`state.unattributed_cost`), NOT silently dropped — so the
#     per-program rows + the unattributed bucket SUM EXACTLY to the aggregate `:valuation_cost` row
#     (the invariant asserted in test/semantic/program_ledger.jl). This is the honest scope of
#     finding D's "what is cleanly attributable": cost at the bind site, reward at the finishing
#     transition; pool-level spend with no bound program stays pool-level.
#
# REWARD is attributed at `finish!`: when a transition completes, its realized reward
# (`Σ placeReward · q · stoich` over RHS products, the per-transition contribution to the aggregate
# `:valuation_reward` row) is split EVENLY across the tokens that were bound to that finishing
# transition (same rule as cost). For an @advance/@move pipeline the produced token IS the bound
# program (identity preserved, ADR 0008 §D), so reward lands on the program that advanced. Reward
# from a finishing transition with no bound program (a plain reaction) goes to the unattributed
# bucket, preserving the same sum invariant against `:valuation_reward`.
#
# ── Determinism (§4 D4) ─────────────────────────────────────────────────────────────────
# Attribution draws NO randomness and is append-only. The per-tick `:program_ledger` log row
# iterates tokens in the deterministic (species, creation_index) total order (`token_sortkey`,
# src/predicates.jl), so the ledger is byte-for-byte reproducible under a fixed seed and is rebuilt
# (cleared) by `_reinit!` exactly like `creation_counters` (closing §4 D7 for the ledger).

export ProgramLedger, program_ledger, program_ledger_entries

# (The `ProgramLedger` accumulator struct is defined in src/state.jl — it must precede the
# ReactionNetworkProblem field that names its element type. All of its LOGIC — attribution, the
# per-tick row, the query API, and the reinit reset — lives here.)

# Fetch (or lazily create) the ledger for a token, keyed by its stable network identity
# (`AlgebraicAgents.getname`, the same key as `creation_index`/`init_snapshot`). Lazy creation
# means a token that is never bound to a costed transition still gets a (zeroed) ledger the first
# time it is touched, and tokens injected mid-run (ADR 0010 AddToken) are tracked from first bind.
function _program_ledger!(state::ReactionNetworkProblem, token)
    name = AlgebraicAgents.getname(token)
    led = get(state.program_ledgers, name, nothing)
    if led === nothing
        ci = get(state.creation_index, name, 0)
        led = ProgramLedger(get_species(token), ci)
        state.program_ledgers[name] = led
    end
    return led
end

# Cost of one transition this tick: `Σ_s consumed[s] · placeCost[s]` over the species column
# `consumed` (one transition's allocation). This is the per-transition decomposition of the
# aggregate `actual_allocs' · placeCost` the `:valuation_cost` row sums.
function _transition_cost(state::ReactionNetworkProblem, consumed::AbstractVector)
    c = 0.0
    for s in row_ids(state, :S)
        @inbounds c += consumed[s] * state[s, :placeCost]
    end
    return c
end

# The structured tokens bound to a transition this tick — both the upfront/blocking binds
# (`bound_structured_agents`) and the nonblock binds (`nonblock_structured_agents`). The cost a
# transition consumed is split EVENLY across these (see the attribution rule at the top).
function _bound_tokens(transition::Transition)
    return vcat(transition.bound_structured_agents, transition.nonblock_structured_agents)
end

"""
    attribute_cost!(state, transition, consumed)

Attribute the COST a `transition` consumed this tick (`Σ consumed[s]·placeCost[s]`) to the
structured token(s) bound to it, split evenly (MVP finding D attribution rule, see this file's
header). A transition with no bound program books its cost against the network UNATTRIBUTED
bucket. Append-only; draws no RNG. Returns the cost it accounted for (so the caller can assert the
per-transition sum equals the tick aggregate).
"""
function attribute_cost!(state::ReactionNetworkProblem, transition::Transition, consumed::AbstractVector)
    cost = _transition_cost(state, consumed)
    cost == 0.0 && return 0.0
    toks = _bound_tokens(transition)
    if isempty(toks)
        state.unattributed_cost += cost
        return cost
    end
    share = cost / length(toks)
    tname = AlgebraicAgents.getname(transition)
    for tok in toks
        led = _program_ledger!(state, tok)
        led.cost_incurred += share
        push!(led.entries, (state.t, :cost, share, tname))
    end
    return cost
end

"""
    attribute_reward!(state, transition, tokens, reward)

Attribute the REWARD a finishing `transition` realized this tick (`Σ placeReward·q·stoich` over its
RHS products) to `tokens` — the programs that were bound to it — split evenly (same rule as cost).
For an @advance/@move pipeline the produced token IS the bound program (identity preserved, ADR
0008 §D), so the reward lands on the advancing program. A finishing transition with no bound program
(`tokens` empty) books its reward against the UNATTRIBUTED bucket, preserving the sum invariant
against `:valuation_reward`.

`tokens` MUST be the bound list snapshotted at `finish!` BEFORE the RHS emission: an @advance/@move
RHS op moves its token out of `transition.bound_structured_agents` mid-loop, so reading the bind
list afterward would lose exactly the program that earned the reward.
"""
function attribute_reward!(
        state::ReactionNetworkProblem,
        transition::Transition,
        tokens::AbstractVector,
        reward::Real,
    )
    reward == 0.0 && return 0.0
    if isempty(tokens)
        state.unattributed_reward += reward
        return reward
    end
    share = reward / length(tokens)
    tname = AlgebraicAgents.getname(transition)
    for tok in tokens
        led = _program_ledger!(state, tok)
        led.reward_realized += share
        push!(led.entries, (state.t, :reward, share, tname))
    end
    return reward
end

# Recompute each program's mark-to-market valuation as its species' `placeValuation` unit value
# (so the per-program valuations of live tokens sum to the structured part of the aggregate
# `:valuation` row). Overwrites `valuation` (it is a STOCK, not a flow — unlike cost/reward which
# accumulate), so it is NOT appended to `entries`. Iterated in deterministic token order. Tokens
# whose species carries no `placeValuation` (the BD case, where valuation is a post-hoc rNPV roll-up,
# MVP finding D) keep valuation 0.0 here — the demo reads cost/reward from this ledger and computes
# rNPV itself.
function attribute_valuation!(state::ReactionNetworkProblem)
    container = getagent(state, "structured")
    for tok in collect(values(inners(container)))
        sp = get_species(tok)
        sp === nothing && continue
        i = find_index(sp, state)
        i === nothing && continue
        led = _program_ledger!(state, tok)
        led.valuation = isblocked(tok) ? led.valuation : state[i, :placeValuation]
    end
    return state
end

# Push the per-tick per-program ledger row into state.log, consistent with the aggregate rows.
# A `(:program_ledger, t, Dict(token_name => (cost, reward, valuation)))` row whose per-program
# `cost`/`reward` entries (PLUS the tick's unattributed deltas, carried on the row) reconcile to the
# aggregate `:valuation_cost`/`:valuation_reward` rows of the same tick. Tokens are iterated in the
# deterministic (species, creation_index) order so the row is reproducible (§4 D4). `cost`/`reward`
# here are the RUNNING totals (matching the running `cost_incurred`/`reward_realized` fields); a
# consumer wanting per-tick flow diffs successive rows.
function push_program_ledger_row!(state::ReactionNetworkProblem)
    snapshot = Dict{String, NamedTuple{(:cost, :reward, :valuation), Tuple{Float64, Float64, Float64}}}()
    container = getagent(state, "structured")
    toks = collect(values(inners(container)))
    sort!(toks; by = a -> token_sortkey(state, a))
    for tok in toks
        led = _program_ledger!(state, tok)
        snapshot[AlgebraicAgents.getname(tok)] =
            (cost = led.cost_incurred, reward = led.reward_realized, valuation = led.valuation)
    end
    push!(state.log, (:program_ledger, state.t, snapshot))
    return state
end

# ── User-facing query API (the seam the BD demo reads) ──────────────────────────────────

"""
    program_ledger(state) -> DataFrame

Per-program (per-structured-token) cost/reward/valuation summary for a finished (or in-progress)
run, in the deterministic (species, creation_index) token order (§4 D4) — the engine-level
replacement for the BD demo's post-hoc reconstruction (MVP finding D). Columns:

  `program`         the token's stable network identity (`AlgebraicAgents.getname`)
  `species`         the token's CURRENT species/kind (`:removed` if soft-retired)
  `creation_index`  the per-species monotonic creation index (ADR 0006 §E) — the order key
  `cost_incurred`   total capital burned on behalf of this program (sum of its bind-cost shares)
  `reward_realized` total reward credited when a transition it was bound to finished successfully
  `valuation`       current mark-to-market = the species' `placeValuation` (0 when none — see header)
  `net`             reward_realized − cost_incurred (the realized economics to date)

The per-program `cost_incurred` summed over ALL programs PLUS `state.unattributed_cost` equals the
sum of the aggregate `:valuation_cost` rows (likewise reward). Pass the live `state`; only the
structured tokens that have ever existed appear.
"""
function program_ledger(state::ReactionNetworkProblem)
    container = getagent(state, "structured")
    live = collect(values(inners(container)))
    # Order: live tokens by deterministic sortkey; then any retired/disentangled programs that have
    # a ledger but are no longer in the container (so their economics are not lost), by stored key.
    sort!(live; by = a -> token_sortkey(state, a))
    rows_name = String[]
    seen = Set{String}()
    for tok in live
        push!(rows_name, AlgebraicAgents.getname(tok))
        push!(seen, AlgebraicAgents.getname(tok))
    end
    extra = sort(
        [k for k in keys(state.program_ledgers) if !(k in seen)];
        by = k -> (string(state.program_ledgers[k].species), state.program_ledgers[k].creation_index, k),
    )
    append!(rows_name, extra)

    live_by_name = Dict(AlgebraicAgents.getname(t) => t for t in live)
    df = DataFrame(
        program = String[],
        species = Symbol[],
        creation_index = Int[],
        cost_incurred = Float64[],
        reward_realized = Float64[],
        valuation = Float64[],
        net = Float64[],
    )
    for name in rows_name
        led = get(state.program_ledgers, name, nothing)
        led === nothing && continue
        # current species: prefer the live token (it may have advanced/retired since first bind)
        sp = haskey(live_by_name, name) ? get_species(live_by_name[name]) : led.species
        push!(
            df,
            (
                name,
                sp,
                led.creation_index,
                led.cost_incurred,
                led.reward_realized,
                led.valuation,
                led.reward_realized - led.cost_incurred,
            ),
        )
    end
    return df
end

"""
    program_ledger_entries(state, token_name) -> Vector{Tuple{Float64,Symbol,Float64,String}}

The append-only per-event audit trail for one program: `(t, kind, amount, transition_name)` rows
(`kind ∈ (:cost, :reward)`), in attribution order. Empty for an unknown/never-bound program.
"""
function program_ledger_entries(state::ReactionNetworkProblem, token_name::AbstractString)
    led = get(state.program_ledgers, String(token_name), nothing)
    return led === nothing ? Tuple{Float64, Symbol, Float64, String}[] : led.entries
end

# Reset the per-program ledger to its t=0 (empty) state — mirrors how `_reinit!` clears
# `creation_counters`. Called from `_reinit!` so `init → step* → reinit! → step*` reproduces the
# first run's program ledger (§4 D7).
function reset_program_ledger!(state::ReactionNetworkProblem)
    empty!(state.program_ledgers)
    state.unattributed_cost = 0.0
    state.unattributed_reward = 0.0
    return state
end
