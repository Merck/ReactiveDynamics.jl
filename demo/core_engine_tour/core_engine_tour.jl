# =============================================================================
# ReactiveDynamics.jl — CORE ENGINE TOUR (classical / plain-Float64 place)
# =============================================================================
#
# This is a single, runnable, literate walkthrough of the ReactiveDynamics
# engine's core modeling vocabulary. It deliberately stays in the CLASSICAL
# regime — every place is a plain counted quantity (a Float64 stock, like a
# population of molecules, dollars, scientists, or jobs). The engine also
# supports STRUCTURED / agent tokens (programs with attributes that move through
# a lifecycle and carry identity); those are a separate demo (see
# demo/bd_acquisition). Here we want to see, end to end, what the bare reaction
# engine can express WITHOUT any of that machinery.
#
# Run it:   julia --project=. demo/core_engine_tour/core_engine_tour.jl
#
# Every construct below is copied from the package's passing semantic test
# suite (test/semantic/*.jl) — this file invents no API. As you read, the code
# is preceded by an explanation of the MODELING IDEA and followed by narrated
# `println` output, so running it tells a story.
#
# The mental model in one paragraph: a model is a set of TRANSITIONS. Each
# transition has a RATE (how often it tries to fire), a left-hand side of
# INPUTS it consumes from its resource pools ("places", in Petri-net terms),
# and a right-hand side of PRODUCTS it emits. Firing can
# be instantaneous (cycletime 0) or take time (cycletime > 0, an "in-flight
# instance"). Arcs can be consumed outright, held-and-returned, metered
# per-step, and so on — the engine's signature feature is this RESOURCE MODALITY
# system. When several transitions want the same scarce pool in the same tick,
# a priority-weighted ALLOCATOR rations it. Randomness (Poisson genesis,
# Binomial success) flows through a per-run seeded RNG, so a run is reproducible
# from (model, seed). That is the whole tour.

using ReactiveDynamics
using ReactiveDynamics: nrows           # row-counting store accessor for composed networks
using Statistics                        # mean/std for the ensemble section

banner(title) = (println(); println("="^74); println(title); println("="^74))


# =============================================================================
banner("§1. A first model: SIR — the metalanguage, simulate, and an invariant")
# =============================================================================
#
# The classic susceptible–infected–recovered epidemic. It introduces every part
# of the authoring surface you will reuse for the rest of the tour:
#
#   * @reaction_network begin ... end  — the model DSL. Each line is
#       `rate, LHS --> RHS, name => ...`. Here both rates are mass-action
#       expressions in the places and parameters (α·S·I, β·I): a bare numeric
#       expression is a STOCHASTIC (Poisson) rate.
#   * @prob_init  — initial counts (the marking at t=0).
#   * @prob_params — the named parameters the rate expressions reference.
#   * @prob_meta  — the simulation horizon `tspan` and the time step `dt`.
#   * ReactionNetworkProblem(...; seed=) — compile the schema into a runnable
#       problem. The `seed=` kwarg owns a per-run RNG, which is the ONLY route
#       to reproducibility (more on that in §8).
#   * simulate(prob) — advance to tspan. The solution lands in `prob.sol`, a
#       DataFrame with a "t" column plus one column per place.
#
# S+I→2I converts one S into one I (net −1 S, +1 I); I→R converts one I into
# one R. No place is created or destroyed outright, so S+I+R is a structural
# INVARIANT — a sanity check the engine should preserve exactly.

sir = @reaction_network begin
    α * S * I, S + I --> 2I, name => I2R
    β * I, I --> R, name => R2S
end
@prob_init sir S = 999 I = 10 R = 0
@prob_params sir α = 0.0001 β = 0.01
@prob_meta sir tspan = 250 dt = 0.1

sir_prob = ReactionNetworkProblem(sir; seed = 1)
simulate(sir_prob)

# IMPORTANT: read solution columns BY NAME. Column order is CONSTRUCTION order,
# not the order you wrote the place, so positional indexing is a foot-gun.
S = sir_prob.sol[!, "S"]
I = sir_prob.sol[!, "I"]
R = sir_prob.sol[!, "R"]
total = S .+ I .+ R

peak_ix = argmax(I)
println("Initial population S+I+R           : ", total[1])
println(
    "Population invariant holds (max-min): ",
    round(maximum(total) - minimum(total); digits = 9), "  (≈ 0 ⇒ conserved)"
)
println(
    "Epidemic peak |I|                   : ", round(maximum(I); digits = 1),
    "  vs initial I = ", I[1]
)
println(
    "Peak occurs at t                    : ",
    round(sir_prob.sol[!, "t"][peak_ix]; digits = 1), " (interior ⇒ a genuine outbreak)"
)
println(
    "Infected at horizon                 : ", round(I[end]; digits = 1),
    "  (declines after the peak)"
)
println("`prob.u` final state vector         : ", round.(sir_prob.u; digits = 1))
println("Index of :S in the state vector     : ", ReactiveDynamics.find_index(:S, sir_prob))


# =============================================================================
banner("§2. Stateful transitions & lifecycle: cycletime, probability, capacity, maxlifetime")
# =============================================================================
#
# In §1 every reaction completed in the same tick it fired (cycletime 0). Real
# processes take TIME and can FAIL. The engine models this with per-transition
# attributes attached after the reaction (key => value, comma-separated):
#
#   * cycletime => c   — a fired instance stays "in-flight" and completes after
#       ceil(c/dt) ticks. Until then it occupies the transition and holds any
#       conserved resources (§3).
#   * probability => p — on completion, success is a Binomial(q, p) draw; only
#       successful instances emit the RHS. (alias: `prob`)
#   * capacity => k    — at most k instances may be concurrently in-flight;
#       proposals beyond k are DEFERRED to later ticks, never dropped.
#   * maxlifetime => L — an instance that has not completed by age L TIMES OUT;
#       a timed-out instance yields zero successes and is pruned cleanly.
#
# We model a small pipeline: jobs start from an ample `budget`, each taking 3
# ticks, succeeding 60% of the time, with at most 4 running concurrently.
# Because cycletime = 3 and dt = 1, the FIRST product cannot appear before t = 3.

pipeline = @reaction_network begin
    @deterministic(1.0), budget --> product,
        name => job, cycletime => 3.0, probability => 0.6, capacity => 4
end
@prob_init pipeline budget = 1000 product = 0
@prob_params pipeline
@prob_meta pipeline tspan = 20 dt = 1.0

pipe_prob = ReactionNetworkProblem(pipeline; seed = 7)
simulate(pipe_prob)

t = pipe_prob.sol[!, "t"]
prod = pipe_prob.sol[!, "product"]
first_completion = t[findfirst(>(0.0), prod)]
# Live concurrency = instances currently in-flight; the capacity gate bounds it.
h = pipe_prob[1, :transHash]
inflight = count(tr -> tr[:transHash] == h, pipe_prob.ongoing_transitions)

println("cycletime = 3.0, dt = 1.0  ⇒ no product can appear before t = 3")
println("First product appears at t : ", first_completion, "  (cycletime delay, as expected)")
println(
    "Products completed by end   : ", Int(prod[end]),
    "  (~60% of started jobs succeed; the rest fail the Binomial draw)"
)
println("In-flight instances at end  : ", inflight, "  (bounded by capacity => 4)")


# =============================================================================
banner("§3. Resource modalities — the engine's signature feature (a truth-table tour)")
# =============================================================================
#
# A arc is not just "consumed". The engine has a small ALGEBRA of resource
# behaviors, set by wrapping the place in a modality macro on the LHS. The
# behavior depends on WHEN the resource is drawn and WHETHER it comes back:
#
#   bare  X        — RAW CONSUMED: debited at spawn, never returned (mass burned).
#   @conserved(X)  — HELD then RETURNED in full at finish (a reusable resource
#                     that is merely occupied during the instance's cycle).
#   @rate(X)       — METERED per ongoing tick: draws q·stoich·Δt each step the
#                     instance is alive (a flow that is consumed continuously).
#                     REQUIRES cycletime > 0, or it silently reserves NOTHING.
#   @nonblock(X)   — held but FREED every step (a soft hold). Also needs ct > 0.
#   @rate(@conserved(X)) — stacked: BOTH tags. A "rented hold" — drawn per tick
#                     but fully credited back at finish (net occupancy is small).
#
# We exercise each in its own tiny model and read the pool's trajectory. The
# tell is the SHAPE of the pool over time: raw drains monotonically; conserved
# plateaus above zero; rate keeps draining (gated on ct); rented plateaus high.

# --- 3a. RAW CONSUMED: 2 material per firing, never returned -----------------
raw = @reaction_network begin
    @deterministic(1.0), 2 * material --> widget, name => build
end
@prob_init raw material = 100 widget = 0
@prob_params raw
raw_prob = ReactionNetworkProblem(raw, Dict(); tspan = 3, dt = 1.0)
simulate(raw_prob)
println(
    "3a. raw  `2*material --> widget`    : material ",
    raw_prob.sol[!, "material"][1], " → ", raw_prob.sol[!, "material"][end],
    "  (monotone drain; consumed mass never comes back)"
)

# --- 3b. @conserved: 3 cash held for 3 ticks, then returned ------------------
# One holder spawns per tick; each ties up 3 cash for its 3-tick cycle, then
# returns it. The in-flight backlog is bounded, so the pool settles to a steady
# FLOOR strictly above zero — proof the resource was held, not consumed.
cons = @reaction_network begin
    @deterministic(1.0), 3 * @conserved(cash) --> product, name => hold, cycletime => 3.0
end
@prob_init cons cash = 100 product = 0
@prob_params cons
cons_prob = ReactionNetworkProblem(cons, Dict(); tspan = 12, dt = 1.0)
simulate(cons_prob)
cash_tail = cons_prob.sol[!, "cash"][(end - 3):end]
println(
    "3b. @conserved(cash)                : cash steady floor ", cash_tail[end],
    "  (held during cycle, returned in full ⇒ plateaus above 0)"
)

# --- 3c. @rate: fuel metered per ongoing tick (gated on cycletime) -----------
# One instance starts per tick and lives 3 ticks, each drawing 1 fuel/tick. As
# the in-flight population builds to 3, the per-tick draw ramps 1, 2, 3, then
# saturates — a continuously consumed FLOW, not a one-shot debit.
rate = @reaction_network begin
    @deterministic(1.0), @rate(fuel) --> trip, name => drive, cycletime => 3.0
end
@prob_init rate fuel = 1000 trip = 0
@prob_params rate
rate_prob = ReactionNetworkProblem(rate, Dict(); tspan = 6, dt = 1.0)
simulate(rate_prob)
draws = -diff(rate_prob.sol[!, "fuel"])
println(
    "3c. @rate(fuel) (ct=3)               : per-tick draws ", Int.(draws[1:4]),
    "...  (metered q·s·Δt each ongoing tick, ramps then saturates at 3 concurrent)"
)

# --- 3d. The @rate cycletime=0 FOOT-GUN --------------------------------------
# @rate's per-step draw is GATED on cycletime > 0. With the default cycletime 0,
# an instance never persists across a tick boundary, so the @rate resource is
# NEVER touched — the token becomes a silent free input. A real trap worth
# seeing explicitly: `out` still grows while `fuel` never moves.
footgun = @reaction_network begin
    @deterministic(1.0), @rate(fuel) --> out, name => r0
end
@prob_init footgun fuel = 100 out = 0
@prob_params footgun
fg_prob = ReactionNetworkProblem(footgun, Dict(); tspan = 4, dt = 1.0)
simulate(fg_prob)
println(
    "3d. @rate FOOT-GUN (ct defaults 0)  : fuel ",
    fg_prob.sol[!, "fuel"][1], " → ", fg_prob.sol[!, "fuel"][end],
    " (UNTOUCHED!) while out → ", Int(fg_prob.sol[!, "out"][end]),
    "  ⇒ @rate needs cycletime > 0"
)

# --- 3e. @nonblock: held but freed every step --------------------------------
# A soft hold: the sensor is reserved while the instance runs but credited back
# each step, so the pool stays non-negative and finite (it does not drain away).
nb = @reaction_network begin
    @deterministic(1.0), @nonblock(sensor) --> reading, name => measure, cycletime => 3.0
end
@prob_init nb sensor = 10 reading = 0
@prob_params nb
nb_prob = ReactionNetworkProblem(nb, Dict(); tspan = 5, dt = 1.0)
simulate(nb_prob)
println(
    "3e. @nonblock(sensor) (ct=3)        : sensor ",
    nb_prob.sol[!, "sensor"][1], " → ", nb_prob.sol[!, "sensor"][end],
    "  (freed every step ⇒ stays non-negative, does not drain away)"
)

# --- 3f. Stacked @rate(@conserved(...)): a rented hold -----------------------
# Both tags attach. The resource is drawn per tick like @rate, but the full
# per-tick integral is credited back at finish like @conserved. Net effect: a
# rented throughput whose pool plateaus HIGH — only one cohort's reservation is
# ever outstanding.
rented = @reaction_network begin
    @deterministic(1.0), 2 * @rate(@conserved(fuel)) --> made, name => rc, cycletime => 2.0
end
@prob_init rented fuel = 1000 made = 0
@prob_params rented
rented_prob = ReactionNetworkProblem(rented, Dict(); tspan = 8, dt = 1.0)
simulate(rented_prob)
println(
    "3f. 2*@rate(@conserved(fuel))       : fuel steady floor ",
    rented_prob.sol[!, "fuel"][end],
    "  (drawn per tick BUT fully returned at finish ⇒ small net hold)"
)


# =============================================================================
banner("§4. The priority-weighted allocator under genuine contention (ADR 0002)")
# =============================================================================
#
# When several transitions want the same scarce pool in one tick, the engine
# rations it with a PRIORITY-WEIGHTED progressive-filling allocator (ADR 0002).
# The core routine is `progressive_fill!(ws, supply, weights; fmax)`, returning a
# per-transition fill fraction `f`; the realized allocation is `ws.req .* f'`.
# Each transition's fill grows at a rate proportional to its priority weight; a
# transition freezes when it hits its cap (`fmax`) or a resource it needs runs
# out. It is WORK-CONSERVING (nothing usable is left idle) and conjunctive-
# consistent (nothing is stranded). Let's see both regimes directly.
#
# Two requesters each demand 5 from a supply of 8 (uncapped, fmax=Inf ⇒ fill
# until the resource exhausts), with priority weights 1 and 3. At equal demand
# the split equals the priority ratio, and the whole supply is used.

reqs = reshape([5.0, 5.0], 1, 2)        # req[resource, transition]
supply = [8.0]
weights = [1.0, 3.0]
ws = ReactiveDynamics.AllocWorkspace(reqs)
f = ReactiveDynamics.progressive_fill!(ws, supply, weights; fmax = [Inf, Inf])
allocs = vec(ws.req .* f')
println("Direct call — contended (supply 8 < uncapped demand), weights 1:3")
println(
    "  allocation          : ", allocs, "  (ratio ", round(allocs[2] / allocs[1]; digits = 2),
    " ≈ 3.0, the priority ratio)"
)
println("  sum allocated       : ", sum(allocs), "  (= supply 8 ⇒ work-conserving)")

# The no-contention regime: cap each requester at its full demand (fmax = 1 unit
# of fill each) over an ample supply — both granted in full, priority irrelevant.
ws_slack = ReactiveDynamics.AllocWorkspace(reqs)
f_slack = ReactiveDynamics.progressive_fill!(ws_slack, [20.0], weights; fmax = [1.0, 1.0])
allocs_slack = vec(ws_slack.req .* f_slack')
println("Direct call — SLACK (supply 20 ≥ demand 10)")
println("  allocation          : ", allocs_slack, "  (granted in full; priority is irrelevant when nothing is scarce)")

# --- The same rationing INSIDE a running model -------------------------------
# Now make two transitions compete for a scarce shared `cash` pool, each held
# via @conserved with cycletime so the reservation persists and the pools draw
# down. They draw the SAME amount but carry different priorities (1 vs 3). With
# financing calibrated to keep cash genuinely scarce, the higher-priority
# transition should win more instances — visible in the product counts.

contend = @reaction_network begin
    @deterministic(3.0), 4 * @conserved(cash) --> lowprod,
        name => low, cycletime => 2.0, priority => 1.0
    @deterministic(3.0), 4 * @conserved(cash) --> highprod,
        name => high, cycletime => 2.0, priority => 3.0
    @deterministic(6.0), ∅ --> cash, name => financing   # steady but insufficient inflow
end
@prob_init contend cash = 12 lowprod = 0 highprod = 0
@prob_params contend
@prob_meta contend tspan = 30 dt = 1.0
contend_prob = ReactionNetworkProblem(contend; seed = 4)
simulate(contend_prob)
lo = contend_prob.sol[!, "lowprod"][end]
hi = contend_prob.sol[!, "highprod"][end]
println("In-model contention for a scarce `cash` pool (both demand 4, priority 1 vs 3):")
println("  low-priority output : ", Int(lo))
println(
    "  high-priority output: ", Int(hi),
    hi > lo ? "  ⇒ the higher-priority transition won more of the scarce resource" :
        "  (priority allocator active)"
)


# =============================================================================
banner("§5. Genesis modes: source (∅), flow/routing, and scheduled (@periodic)")
# =============================================================================
#
# "Genesis" is how new tokens enter a model. Three distinct idioms:
#
#   * SOURCE (∅ LHS)  — an empty left-hand side bypasses the input gate, so the
#       rate IS the realized count. A Poisson source `2.0, ∅ --> arrival` is
#       dt-invariant in expectation (E = rate·tspan). CAVEAT: a FRACTIONAL
#       @deterministic count on a source is NOT dt-invariant (a known ceil
#       rounding behavior) — use INTEGER deterministic counts for sources.
#   * FLOW / ROUTING  — a non-empty LHS with a high nominal rate is TOKEN-GATED:
#       it fires bounded by the available upstream tokens, not the nominal rate.
#       The idiom for "route whatever is available": high rate + an upstream
#       place on the LHS.
#   * SCHEDULED (@periodic) — `@deterministic(N * @periodic(p))` fires N spawns
#       at each multiple of period p, nothing in between — a calendar/batch
#       intake. NOTE: the MACRO `@periodic(p)` must appear INSIDE the rate; a
#       bare `periodic(p)` call errors.

# --- 5a. Poisson source: empty LHS, dt-invariant in expectation --------------
function source_total(dt; seed)
    src = @reaction_network begin
        2.0, ∅ --> arrival, name => inflow
    end
    @prob_init src arrival = 0
    @prob_params src
    p = ReactionNetworkProblem(src, Dict(); tspan = 50.0, dt = dt, seed = seed)
    simulate(p)
    return p.u[1]
end
m_dt1 = mean(source_total(1.0; seed = s) for s in 1:60)
m_dt2 = mean(source_total(0.5; seed = s) for s in 1:60)
println("5a. Poisson source  `2.0, ∅ --> arrival`  (E = rate·tspan = 2·50 = 100)")
println("    ensemble mean total @ dt=1.0 : ", round(m_dt1; digits = 1))
println(
    "    ensemble mean total @ dt=0.5 : ", round(m_dt2; digits = 1),
    "  ⇒ halving dt preserves the expected total (dt-invariant)"
)

# --- 5b. Flow / routing: non-empty LHS is token-gated ------------------------
# `upstream` deposits 2 feed/tick; `router` has nominal rate 100 but can only
# route what feed actually holds, so realized routing tracks the 2/tick deposit
# — NOT the nominal 100.
flow = @reaction_network begin
    @deterministic(2.0), ∅ --> feed, name => upstream
    @deterministic(100.0), feed --> product, name => router
end
@prob_init flow feed = 0 product = 0
@prob_params flow
flow_prob = ReactionNetworkProblem(flow, Dict(); tspan = 5, dt = 1.0)
simulate(flow_prob)
println("5b. Flow/routing (nominal rate 100, fed 2/tick)")
println(
    "    product at t=0,1 : ", flow_prob.sol[!, "product"][1], ", ",
    flow_prob.sol[!, "product"][2], "  (0 while feed empty)"
)
println(
    "    max routed/tick  : ", maximum(diff(flow_prob.sol[!, "product"])),
    "  ⇒ token-gated to the 2/tick supply, NOT the nominal 100"
)

# --- 5c. Scheduled: @periodic fires N at each calendar boundary --------------
sched = @reaction_network begin
    @deterministic(3 * @periodic(2.0)), ∅ --> cohort, name => intake
end
@prob_init sched cohort = 0
@prob_params sched
sched_prob = ReactionNetworkProblem(sched, Dict(); tspan = 7, dt = 1.0)
simulate(sched_prob)
deltas = diff(sched_prob.sol[!, "cohort"])
println("5c. Scheduled  `@deterministic(3 * @periodic(2.0))`  (period 2.0, 3 per boundary)")
println(
    "    cohort end total : ", Int(sched_prob.sol[!, "cohort"][end]),
    "  (3 boundaries at t=2,4,6 × 3 each = 9)"
)
println("    per-tick deltas  : ", Int.(deltas), "  (spawns only at period boundaries, flat between)")


# =============================================================================
banner("§6. Custom registered rate functions + the cost / reward / valuation ledger")
# =============================================================================
#
# Rates need not be closed-form mass-action expressions: you can @register an
# ordinary Julia function and call it as a rate. Registered functions live in
# the engine's scope and must be defined BEFORE the problem is built. Here a
# toy-pharma pipeline: discovery turns scientists+budget into candidate
# compounds at a registered rate α(...), and dx2market converts a candidate into
# a marketed drug at a registered rate β(...).
#
# This section also tours the VALUATION LEDGER. Attach @cost / @reward /
# @valuation to place and the engine records financial events to `prob.log`,
# a vector of NamedTuple-like rows. Read it by tag:
#   costs  = [r[3] for r in prob.log if r[1] == :valuation_cost]
#   rewards = [r[3] for r in prob.log if r[1] == :valuation_reward]
# Discounting/NPV is pure post-processing — the engine does no discounting.

@register function α(n1, n2, κ)
    return κ + exp(-n1) + exp(-n2)
end
@register function β(n1, n2)
    return n1 + exp(-n2)
end

toy = @reaction_network begin
    α(candidate_compound, marketed_drug, κ),
        3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
        name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
    β(candidate_compound, marketed_drug),
        candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
        name => dx2market, probability => 0.5, cycletime => 4
    γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
end
@periodic toy 1.0 budget += 11 * marketed_drug          # financing tied to revenue
@prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
@prob_params toy κ = 4 γ = 0.1
@cost toy budget = 1.0 scientist = 2.0                  # cost per unit of resource used
@reward toy marketed_drug = 50.0                        # reward credited on a successful launch
@valuation toy marketed_drug = 100.0
@prob_meta toy tspan = 50 dt = 0.1

toy_prob = ReactionNetworkProblem(toy; seed = 1)
simulate(toy_prob)

costs = [r[3] for r in toy_prob.log if r[1] == :valuation_cost]
rewards = [r[3] for r in toy_prob.log if r[1] == :valuation_reward]
# A simple discounted rNPV reduction over the ledger (post-processing only).
rnpv = sum(
    (
            row[1] == :valuation_reward ? row[3] :
            row[1] == :valuation_cost ? -row[3] : 0.0
        ) / (1 + 0.1)^row[2]
        for row in toy_prob.log
)
println("Toy-pharma pipeline with @register'd rates α, β")
println("  solution columns   : ", names(toy_prob.sol), "  (construction order ⇒ index by name)")
println(
    "  scientists held     : min ", round(minimum(toy_prob.sol[!, "scientist"]); digits = 1),
    " / max ", round(maximum(toy_prob.sol[!, "scientist"]); digits = 1),
    "  (≤ 20 initial ⇒ @conserved never overruns its holding)"
)
println("  ledger cost rows    : ", length(costs), "  totalling ", round(sum(costs); digits = 1))
println("  ledger reward rows  : ", length(rewards), "  totalling ", round(sum(rewards); digits = 1))
println(
    "  discounted rNPV     : ", round(rnpv; digits = 1),
    "  (a post-processing reduction over the ledger; the engine does no discounting)"
)


# =============================================================================
banner("§7. Composition: @join two submodels and @equalize place")
# =============================================================================
#
# Models compose. `@join` takes the UNION of two schemas' places, transitions,
# and parameters, optionally IDENTIFYING shared places across the two via
# equations. `@equalize` collapses two places WITHIN one schema into a single
# pool and rewrites every reference. Both operate at AUTHORING time (on a
# schema), before construction.
#
# `@join` merges place / transitions / params AND (since WS-3) events (:E) and
# observables (:obs) too — `merge_networks!` walks all six objects and appends :E/:obs
# structurally, so nothing is silently dropped on a join. `@join` / `@equalize` are
# the MANUAL, no-declared-ports path (you name the place to identify); the
# declared-port counterpart is `@compose`, which matches open input/output ports
# automatically (see demo/refinement_tour).
#
# Here two reaction sub-systems each consume a shared resource A; we join them,
# identifying the two A's as one pool, and count parts of the merged schema.

acs1 = @reaction_network begin
    1.0, A --> B, name => t1
end
acs2 = @reaction_network begin
    1.0, A --> C, name => t2
end
joined = @join acs1 acs2 acs1.A = acs2.A = @alias(A)
println("@join acs1 acs2 (identifying the shared place A)")
println(
    "  place in join : ", nrows(joined, :S),
    "  (union {A,B,C} ⇒ 3; the two A's merged into one)"
)
println("  transitions     : ", nrows(joined, :T), "  (1 + 1, none lost)")

# @equalize: two conceptually-identical place A and A2 collapse to one.
eqacs = @reaction_network begin
    1.0, A  --> B, name => t1
    1.0, A2 --> B, name => t2
end
before_S = nrows(eqacs, :S)
equalized = @equalize eqacs A = A2
println("@equalize eqacs A = A2 (collapse A and A2 into one pool)")
println(
    "  place before  : ", before_S, "  → after : ", nrows(equalized, :S),
    "  (dropped by exactly 1; references rewritten)"
)
println("  transitions     : ", nrows(equalized, :T), "  (preserved; only :S was touched)")


# =============================================================================
banner("§8. Determinism & seeded ensembles")
# =============================================================================
#
# Reproducibility comes from the `seed=` construction kwarg ALONE. The state
# owns its own RNG, isolated from Julia's global RNG — so `Random.seed!(n)`
# followed by constructing does NOT pin a run, and a seeded run neither reads
# nor perturbs the global stream. The contract:
#
#   * same (model, seed)      ⇒ byte-identical trajectory and ledger,
#   * different seed          ⇒ (almost surely) a different trajectory,
#   * unseeded                ⇒ a fresh entropy seed each construction (diverges).
#
# For ensembles, derive each member's seed deterministically from a single root
# seed + member index — `hash((root, k))`. Then member k is reproducible from
# (root, k) regardless of how many members you run or in what order — the basis
# for stable Monte-Carlo statistics.

function build_birth()
    net = @reaction_network begin
        3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
    end
    @prob_init net A = 100 B = 0
    @prob_params net
    @prob_meta net tspan = 20 dt = 1.0
    return net
end

a = (p = ReactionNetworkProblem(build_birth(); seed = 42); simulate(p); p.sol[!, "B"])
b = (p = ReactionNetworkProblem(build_birth(); seed = 42); simulate(p); p.sol[!, "B"])
c = (p = ReactionNetworkProblem(build_birth(); seed = 7); simulate(p); p.sol[!, "B"])
println("Same seed (42 vs 42) identical : ", a == b)
println("Different seed (42 vs 7) differ: ", a != c)

# A seeded ensemble: member k is seeded from (root, k); collect the final B.
member_seed(root, k) = hash((root, k))
function run_member(root, k)
    p = ReactionNetworkProblem(build_birth(); seed = member_seed(root, k))
    simulate(p)
    return p.sol[!, "B"][end]
end
ens = [run_member(2026, k) for k in 1:200]
# Reproducible from (root, k), independent of N and order:
solo = run_member(2026, 3)
ens_member_3 = ens[3]
println("Ensemble of 200 members, final B (cycletime=2, probability=0.5):")
println("  mean ± std         : ", round(mean(ens); digits = 2), " ± ", round(std(ens); digits = 2))
println("  range [min, max]   : [", Int(minimum(ens)), ", ", Int(maximum(ens)), "]")
println(
    "  member 3 in-ensemble vs standalone equal : ", ens_member_3 == solo,
    "  (reproducible from (root,k), independent of N/order)"
)


# =============================================================================
banner("§9. Recap — what this tour exercised")
# =============================================================================
println(
    """
      §1  The metalanguage: @reaction_network / @prob_init / @prob_params /
          @prob_meta, mass-action rates, simulate, reading prob.sol by name, and a
          conserved-population invariant on an SIR epidemic.
      §2  Stateful lifecycle: cycletime (in-flight delay), Binomial `probability`,
          `capacity` bound on concurrency, and `maxlifetime` timeout.
      §3  Resource modalities: raw-consumed vs @conserved vs @rate vs @nonblock vs
          the stacked rented hold — plus the @rate-with-cycletime=0 foot-gun.
      §4  The priority-weighted allocator: progressive_fill! directly (contended and
          slack), then genuine in-model contention where higher priority wins.
      §5  Genesis: Poisson source (∅, dt-invariant), token-gated flow/routing, and
          scheduled @periodic batch intake.
      §6  @register'd custom rate functions on a toy-pharma pipeline, and the
          cost / reward / valuation ledger with a discounted-rNPV reduction.
      §7  Composition: @join (union + shared-place identification) and @equalize
          (collapse + rewrite).
      §8  Determinism: seed= reproducibility, divergence on different seeds, and a
          deterministically-seeded ensemble with mean ± spread.

      Everything above used CLASSICAL (plain Float64) place only. Structured /
      agent tokens with attributes and lifecycle identity are a separate demo:
      see demo/bd_acquisition.
    """
)
