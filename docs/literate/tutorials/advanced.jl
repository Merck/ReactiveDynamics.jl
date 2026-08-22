# # Advanced tutorial — a structured-token R&D portfolio under contention
#
# **What you will build.** A portfolio of R&D programs modeled as *structured tokens* — first-class
# entities that each carry a lifecycle `phase` and a net-present value `npv`, keep a stable identity
# as they advance, and compete for a scarce shared resource. We author the lifecycle as a
# phase-as-attribute pipeline, select programs by a value threshold, ration a contended capital pool
# through the priority allocator, add an in-model management lever, and close on the marginal value
# of a financing decision computed across a seeded ensemble.
#
# **Who this is for.** Readers who have done the [introductory tutorial](introductory.md) and are
# comfortable with `@reaction_network`, `simulate`, and reading `prob.sol` by name. The introductory
# tier stayed entirely in the *classical* regime — every quantity was a plain counted stock. Here we
# take the step that distinguishes ReactiveDynamics from a plain reaction engine: tokens with
# attributes and identity, moved through a lifecycle by predicate-selected transitions, under a
# resource algebra with a priority allocator.
#
# **The through-line.** One small R&D portfolio, growing in sophistication section by section:
# structured tokens (§1) → a phase lifecycle (§2) → value-qualified selection (§3) → the resource
# modalities and the allocator that make capital genuinely scarce (§4) → an in-model financing lever
# (§5) → the marginal value of that lever, with a standard error (§6).

using ReactiveDynamics
using ReactiveDynamics: register_token_kind!, get_place, inners, getagent,
    Rule, Seq, SetMarking, AddToken, Log, PopulationEntry
using Statistics                # mean / std for the ensemble reductions
using Distributions             # Normal, for a sampled starting portfolio
using Plots                     # inline figures

# The structured-token TYPE and several helpers live in the ReactiveDynamics module (that is where
# the selection / advancement machinery can see them), so we keep a short qualified alias.
const RD = ReactiveDynamics

# ## 1. Structured tokens: a program is an entity, not a count
#
# A *classical* place is a scalar — a single `Float64` saying "how many `A` there are." That is
# exactly right for indistinguishable molecules, but a program in a portfolio is not a molecule. We
# want each program to carry **attributes** (its current phase, its value) and a stable **identity**
# — the *same* object as it advances Phase1 → Phase2 → …, so a downstream report can follow it.
#
# A *structured token* gives us that: a first-class agent carrying host-Julia fields, whose identity
# (uuid / kind / creation index) is preserved as the engine mutates its fields. We define the kind
# in ReactiveDynamics' own scope with the `@register` / `@aagent` idiom, because the engine's
# selection and advancement machinery lives there and must see the type. The four leading
# constructor arguments are the `@aagent` protocol fields, in order — a unique `name`, the `place`
# kind tag (here `:Project`; every program shares one kind), a `bound_transition` (`nothing` until a
# transition binds the token), and an empty `past_bonds` history — followed by our two modeling
# attributes, `phase` and `npv`.

@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
        phase::Symbol     # lifecycle stage — the canonical "phase-as-attribute"
        npv::Float64      # the program's net-present value (a plain descriptor field)
    end
    function ProjectToken(phase, npv)
        return ProjectToken(
            "Proj" * string(rand(1:(10^9))),                       # name
            :Project,                                              # kind (one kind for all phases)
            nothing,                                               # bound_transition
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[], # past_bonds
            phase,
            npv,
        )
    end
end

# The **registry** maps a kind symbol to a constructor `(state, fields::Dict) -> token`. A serialized
# model, a declarative population, or an injected token all reference host token kinds *by name*, and
# the registry is how those names resolve to real Julia constructors without any data file carrying
# code. The same registry serves the declarative population here and the `AddToken` lever in §5.

const REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RD.ProjectToken(get(f, :phase, :Phase1), get(f, :npv, 100.0)),
)

# A few reductions over the live token pool (the structured container is a `Dict` keyed by token
# name; we usually want the values):

livetokens(p) = collect(values(inners(getagent(p, "structured"))))
nphase(p, ph) = count(t -> get_place(t) == :Project && t.phase == ph, livetokens(p))
nlaunched(p) = nphase(p, :Launched)
nretired(p) = count(t -> get_place(t) == :removed, livetokens(p))

# We seed the starting portfolio with the **declarative** `population` initial marking, passed to the
# constructor and instantiated before `t = 0`. This is preferred over an imperative post-construction
# loop because the run is then reproducible from `(model, population, seed)` — the portfolio is
# reproducible *input*, not host code that runs after the model exists. Two authoring forms exist.
#
# Form A is an explicit host-token list — a fixed, hand-authored opening portfolio:

explicit_portfolio() = [
    RD.ProjectToken(:Phase1, 120.0),
    RD.ProjectToken(:Phase1, 90.0),
    RD.ProjectToken(:Phase2, 200.0),
    RD.ProjectToken(:Phase2, 150.0),
    RD.ProjectToken(:Phase3, 300.0),
]

# Form B is the `PopulationEntry` "count + attribute expressions" form — "N programs with these
# attributes." A symbol literal is wrapped with `QuoteNode`; a value may be a *sampled* expression
# drawn from the run's seeded RNG (`state.rng`), so a same-seed construction is reproducible:

sampled_portfolio() = [
    PopulationEntry(
        :Project, :Project; count = 4,
        attributes = Dict(
            :phase => QuoteNode(:Phase1),
            :npv => :(rand(state.rng, Normal(100.0, 15.0))),
        ),
    ),
    PopulationEntry(
        :Project, :Project; count = 3,
        attributes = Dict(:phase => QuoteNode(:Phase2), :npv => 200.0),
    ),
]

# To *see* a token before we build a model around it, construct one directly and read its fields.
# A demo token, and the two population forms materialized under a seed:

demo_tok = RD.ProjectToken(:Phase1, 120.0)
println("A single ProjectToken : phase=", demo_tok.phase, "  npv=", demo_tok.npv, "  kind=", get_place(demo_tok))
println(
    "Form A explicit portfolio: ", length(explicit_portfolio()), " programs, phases = ",
    sort(string.([t.phase for t in explicit_portfolio()]))
)

# ## 2. A phase-as-attribute lifecycle
#
# A naive design would make a *place per phase* (`Phase1`, `Phase2`, …) and "advance" by destroying
# a `Phase1` token and creating a `Phase2` token. That breaks identity — the new token is a different
# object — and multiplies the place count. The canonical ReactiveDynamics design is
# **phase-as-attribute**: there is one `:Project` kind, and `phase` is a field. A pipeline step reads
#
#     @select(Project, <clause>) --> @advance(phase, :NextPhase)
#
# `@select(Project, clauses)` binds only tokens of kind `Project` whose attributes satisfy the
# conjunctive `&&` clause (operators `== != < <= > >= in`). `@advance(phase, :Phase2)` writes the
# bound token's `phase` field **in place** — the same object, identity preserved.
#
# A stage gate can also *fail*: `probability => q` makes each advance a `Binomial(·, q)` trial. On
# failure the bound token **soft-retires** — its place flips to `:removed` and its `phase` records
# how far it got (a killed Phase2 program stays at `phase == :Phase2` but `place == :removed`).
#
# We build the portfolio's lifecycle as three timed advances. Each advance also holds a shared
# `capital` pool via `@conserved` — capital is *occupied* for the duration of an in-flight advance
# and returned when it completes — so the number of programs that can advance at once is bounded by
# capital. That is the contention we ration in §4 and relieve in §5; here we just watch the pipeline
# run. `priority` orders who wins capital when it is scarce (late-stage programs first). We read one
# tick as one quarter, so `tspan = 12` is a three-year horizon.

function portfolio_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) + 2 * @conserved(capital) --> @advance(phase, :Phase2),
            name => adv12, cycletime => 1.0, probability => 0.9, priority => 1.0
        @deterministic(1.0),
            @select(Project, phase == :Phase2) + 3 * @conserved(capital) --> @advance(phase, :Phase3),
            name => adv23, cycletime => 2.0, probability => 0.6, priority => 2.0
        @deterministic(1.0),
            @select(Project, phase == :Phase3) + 4 * @conserved(capital) --> @advance(phase, :Launched),
            name => adv3L, cycletime => 2.0, probability => 0.9, priority => 3.0
    end
    @prob_init net capital = 18
    register_token_kind!(net, :Project)
    @prob_meta net tspan = 12 dt = 1.0
    return net
end

# Build with the explicit Form-A portfolio, thread a seed, and simulate. The starting portfolio and
# the launch count at the horizon:

p_life = ReactionNetworkProblem(
    portfolio_model(); seed = 1, registry = REGISTRY, population = explicit_portfolio(),
)
println(
    "t=0 by phase : Phase1=", nphase(p_life, :Phase1), "  Phase2=", nphase(p_life, :Phase2),
    "  Phase3=", nphase(p_life, :Phase3), "  Launched=", nphase(p_life, :Launched)
)
simulate(p_life)
println("After 12 quarters:")
println("  Launched              : ", nlaunched(p_life))
println("  still in Phase2/Phase3: ", nphase(p_life, :Phase2), " / ", nphase(p_life, :Phase3))
println("  soft-retired (failed a gate, :removed): ", nretired(p_life))

# The Phase2 → Phase3 gate succeeds only 60% of the time, so some programs soft-retire while others
# reach `:Launched`. A launched program is the *same object* that started in Phase1 — `@advance`
# rewrote its field in place; it never became a different token.

# Form B, materialized: the `PopulationEntry` form builds the same kind from a count plus seeded
# attribute expressions. The sampled Phase1 NPVs are reproducible under the seed.

p_formB = ReactionNetworkProblem(
    portfolio_model(); seed = 42, registry = REGISTRY, population = sampled_portfolio(),
)
println("Form B population: built ", length(livetokens(p_formB)), " programs (4 Phase1 sampled + 3 Phase2).")
println(
    "  sampled Phase1 NPVs (seeded, reproducible): ",
    round.(sort([t.npv for t in livetokens(p_formB) if t.phase == :Phase1]); digits = 1)
)

# ## 3. Value-qualified selection
#
# The power of `@select` is that its predicate is a filter over token attributes: only the matching
# subset is bindable. A continuous clause like `npv > θ` lets a transition act on a value threshold —
# "fast-track only the high-value Phase2 programs." When several tokens match but the transition can
# fire on only a few per tick, *which* bind first is deterministic: equal-priority ties break by
# creation order (the earlier-added token wins), not by dictionary hash order or the tokens' random
# names, so a predicate-selected pipeline reproduces exactly under `(model, seed)`.

function fasttrack_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(phase, :Phase3),
            name => fasttrack, cycletime => 1.0, probability => 1.0
    end
    register_token_kind!(net, :Project)
    @prob_meta net tspan = 3 dt = 1.0
    return net
end

p_sel = ReactionNetworkProblem(
    fasttrack_model(); seed = 1, registry = REGISTRY,
    population = [
        RD.ProjectToken(:Phase2, 100.0),   # below θ = 150 — stays in Phase2
        RD.ProjectToken(:Phase2, 220.0),   # above θ — fast-tracked
        RD.ProjectToken(:Phase2, 180.0),   # above θ — fast-tracked
        RD.ProjectToken(:Phase1, 999.0),   # wrong phase — the && clause gates BOTH phase and npv
    ],
)
println("Predicate: @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(:Phase3)")
println("Before: Phase1=", nphase(p_sel, :Phase1), " Phase2=", nphase(p_sel, :Phase2), " Phase3=", nphase(p_sel, :Phase3))
simulate(p_sel)
println("After : Phase1=", nphase(p_sel, :Phase1), " Phase2=", nphase(p_sel, :Phase2), " Phase3=", nphase(p_sel, :Phase3))

# Only the two high-NPV Phase2 programs advanced; the `npv = 100` program stayed (below θ), and the
# Phase1 program was never eligible (the conjunctive clause gates phase *and* npv). The same
# selection logic drives the dynamics, the analysis queries, and the population-write actions — one
# predicate machinery throughout.

# ## 4. Resource modalities and the allocator under contention
#
# A arc is not simply "consumed." The engine has a small **algebra** of resource behaviors, set
# by wrapping a place in a modality macro on the left-hand side. The behavior depends on *when* the
# resource is drawn and *whether* it comes back — this is the engine's signature feature, and the
# truth table lands here:
#
# | LHS form | meaning | pool shape over time |
# |---|---|---|
# | `X` (bare) | **raw consumed** — debited at spawn, never returned | drains monotonically |
# | `@conserved(X)` | **held then returned** in full at finish | plateaus above zero |
# | `@rate(X)` | **metered** per ongoing tick (needs `cycletime > 0`) | keeps draining while instances live |
# | `@nonblock(X)` | held but **freed every step** (a soft hold) | stays non-negative, does not drain |
# | `@rate(@conserved(X))` | drawn per tick **and** credited back at finish (a rented hold) | plateaus high |
#
# The foot-gun to remember: `@rate`'s per-step draw is *gated on* `cycletime > 0`. With the default
# `cycletime = 0`, an instance never persists across a tick boundary, so a `@rate` reservation would
# never fire. Rather than let that reserve nothing silently, the engine's construction validator
# (CONTRACT §1.4) now *rejects* the combination outright — a modeling mistake caught at build time
# rather than a silently-wrong run. We exercise three legal rows in their own tiny models and read
# each pool's trajectory (the tell is the *shape*), then show the validator refusing the illegal one.

raw = @reaction_network begin
    @deterministic(1.0), 2 * material --> widget, name => build
end
@prob_init raw material = 100 widget = 0
raw_prob = ReactionNetworkProblem(raw, Dict(); tspan = 3, dt = 1.0)
simulate(raw_prob)
println(
    "raw  2*material --> widget : material ", raw_prob.sol[!, "material"][1], " → ",
    raw_prob.sol[!, "material"][end], "  (monotone drain; consumed mass never returns)"
)

cons = @reaction_network begin
    @deterministic(1.0), 3 * @conserved(cash) --> product, name => hold, cycletime => 3.0
end
@prob_init cons cash = 100 product = 0
cons_prob = ReactionNetworkProblem(cons, Dict(); tspan = 12, dt = 1.0)
simulate(cons_prob)
println(
    "@conserved(cash)           : cash steady floor ", cons_prob.sol[!, "cash"][end],
    "  (held during the cycle, returned in full ⇒ plateaus above 0)"
)

rate = @reaction_network begin
    @deterministic(1.0), @rate(fuel) --> trip, name => drive, cycletime => 3.0
end
@prob_init rate fuel = 1000 trip = 0
rate_prob = ReactionNetworkProblem(rate, Dict(); tspan = 6, dt = 1.0)
simulate(rate_prob)
println(
    "@rate(fuel) (ct=3)         : per-tick draws ", Int.((-diff(rate_prob.sol[!, "fuel"]))[1:4]),
    "...  (metered each ongoing tick; ramps then saturates at 3 concurrent)"
)

footgun = @reaction_network begin
    @deterministic(1.0), @rate(fuel) --> out, name => r0
end
@prob_init footgun fuel = 100 out = 0
try
    ReactionNetworkProblem(footgun, Dict(); tspan = 4, dt = 1.0)
    println("@rate FOOT-GUN (ct=0)      : constructed (unexpected)")
catch err
    msg = sprint(showerror, err)
    println(
        "@rate FOOT-GUN (ct=0)      : REJECTED at construction ⇒ ",
        occursin("cycletime == 0 is illegal", msg) ? "CONTRACT §1.4 validator fired (@rate needs cycletime > 0)" : msg
    )
end

# ### The priority-weighted allocator
#
# When several transitions want the same scarce pool in one tick, the engine rations it with a
# priority-weighted progressive-filling allocator. Each transition's fill grows at a rate
# proportional to its priority weight; it freezes when it hits its cap or a resource it needs runs
# out. The routine is work-conserving (nothing usable is left idle). We can call it directly: two
# requesters each demand 5 from a supply of 8, with priority weights 1 and 3.

reqs = reshape([5.0, 5.0], 1, 2)        # req[resource, transition]
ws = RD.AllocWorkspace(reqs)
f = RD.progressive_fill!(ws, [8.0], [1.0, 3.0]; fmax = [Inf, Inf])
allocs = vec(ws.req .* f')
println("progressive_fill! — contended (supply 8 < demand 10), weights 1:3")
println(
    "  allocation : ", round.(allocs; digits = 2), "  (ratio ≈ ",
    round(allocs[2] / allocs[1]; digits = 2), ", the priority ratio; Σ = ", sum(allocs), " = supply)"
)

# The same rationing happens *inside* a running model. Two transitions compete for a scarce shared
# `cash` pool, each holding it via `@conserved` over a cycle so the reservation persists. They demand
# the same amount but carry priorities 1 and 3; with financing calibrated to keep cash genuinely
# scarce, the higher-priority transition should win more instances — visible in the output counts.

contend = @reaction_network begin
    @deterministic(3.0), 4 * @conserved(cash) --> lowprod,
        name => low, cycletime => 2.0, priority => 1.0
    @deterministic(3.0), 4 * @conserved(cash) --> highprod,
        name => high, cycletime => 2.0, priority => 3.0
    @deterministic(6.0), ∅ --> cash, name => financing   # steady but insufficient inflow
end
@prob_init contend cash = 12 lowprod = 0 highprod = 0
@prob_meta contend tspan = 30 dt = 1.0
p_contend = ReactionNetworkProblem(contend; seed = 4)
simulate(p_contend)
lo = Int(p_contend.sol[!, "lowprod"][end]); hi = Int(p_contend.sol[!, "highprod"][end])
println("In-model contention for a scarce `cash` pool (both demand 4; priority 1 vs 3):")
println(
    "  low-priority output : ", lo, "    high-priority output: ", hi,
    hi > lo ? "  ⇒ the higher-priority transition won more" : "  (allocator active)"
)

# This is exactly the mechanism at work in our portfolio model: the three advances hold `@conserved`
# capital with ascending priority, so when capital is tight, late-stage programs (nearest to launch)
# win it first. That is what makes the financing lever in the next section a real decision — capital
# genuinely binds.

# ## 5. An in-model decision rule as a management lever
#
# A management decision — "if conditions hold, take an action" — is itself part of the model, a typed
# `Rule`, not host patch code that reaches in and mutates the run from outside. A rule has
#
#     Rule(id, guard::Expr, action; fire_mode = :once | :every_tick)
#
# The guard is evaluated against the live state (`@t()` is the clock; places and params are in
# scope). `fire_mode = :once` fires the action the first tick its guard holds, then latches off
# (`_reinit!` re-arms it). Actions compose via `Seq`: `SetMarking` injects into a resource pool,
# `SetParams` flips a model parameter, `AddToken` injects a fresh token *by kind* through the
# registry, and `Log` annotates.
#
# Our lever models a **Series-B raise**: once `t > 3`, inject capital into the contended pool and add
# one fresh Phase2 program to the portfolio — a recapitalization that both relieves the binding
# constraint and expands the pipeline, all *in* the model.

raise_lever() = Rule(
    :series_b, :(@t() > 3.0),
    Seq(
        [
            SetMarking(:capital, 40, :inc),                                   # +40 capital into the pool
            AddToken(:Project, [:phase => QuoteNode(:Phase2), :npv => 200.0]), # seed one more Phase2 program
            Log("Series-B raised: +40 capital, +1 Phase2 program"),
        ],
    );
    fire_mode = :once,
)

# Run the portfolio with the lever armed, over the same starting portfolio as §2:

p_lever = ReactionNetworkProblem(
    portfolio_model(); seed = 1, registry = REGISTRY,
    population = explicit_portfolio(), rules = [raise_lever()],
)
simulate(p_lever)
capital_series = p_lever.sol[!, "capital"]
println("Rule: once @t() > 3, Seq[ +40 capital, AddToken(Phase2 npv=200) ]  (fire_mode = :once)")
println("  positive capital jumps : ", count(>(0.0), diff(capital_series)), "  (a :once rule fires exactly once)")
println("  peak capital reached   : ", round(maximum(capital_series); digits = 1), "  (was 18 at t=0)")
println("  rule latched off       : enabled = ", p_lever.rules[1].enabled)
println("  Launched with the raise: ", nlaunched(p_lever), "  (vs ", nlaunched(p_life), " without, same seed)")

# The decision lives in the *model* — the driver only set the seed and armed the rule; no host code
# reached in mid-run. That is what makes the scenario a reproducible `(model, rules, seed)` triple
# rather than an imperative script, and what lets us treat "raise or not" as a clean A/B in the next
# section.

# ## 6. The marginal value of the raise
#
# One run is a single sample of a stochastic process. To value the raise we need a *distribution*, so
# we build two **ensembles** over the same seeds — a baseline (no lever) and a deal arm (lever armed)
# — and compare launches. `ensemble(build; nseed, root_seed)` runs `nseed` independent members,
# member `k` seeded deterministically from `hash((root_seed, k))`, so member `k` is reproducible
# regardless of how many members you run or in what order. The result is an `EnsembleProblem`.

launches(p) = nlaunched(p)

baseline_arm(s) = (
    p = ReactionNetworkProblem(
        portfolio_model(); seed = s, registry = REGISTRY, population = explicit_portfolio(),
    );
    simulate(p); p
)
deal_arm(s) = (
    p = ReactionNetworkProblem(
        portfolio_model(); seed = s, registry = REGISTRY,
        population = explicit_portfolio(), rules = [raise_lever()],
    );
    simulate(p); p
)

base_ens = ensemble(baseline_arm; nseed = 40, root_seed = 2026)
deal_ens = ensemble(deal_arm; nseed = 40, root_seed = 2026)

# `treatment_effect(base, deal, metric)` reduces a per-member scalar metric across both arms to the
# unpaired difference of means with a standard error `√(var_b/n_b + var_d/n_d)`:

te = treatment_effect(base_ens, deal_ens, launches)
println("Launches over the horizon (40-seed ensemble):")
println("  no raise    : ", round(te.baseline; digits = 2))
println("  Series-B    : ", round(te.deal; digits = 2))
println(
    "  the raise buys: +", round(te.delta; digits = 2), " expected launches  (± ",
    round(te.se; digits = 2), " SE)"
)

# The two launch distributions, with their means, make the shift visible — the deal-arm distribution
# sits to the right of the baseline:

base_launches = [launches(m) for m in base_ens.members]
deal_launches = [launches(m) for m in deal_ens.members]
histogram(
    base_launches; bins = -0.5:1:10.5, alpha = 0.5, label = "no raise",
    xlabel = "launches over horizon", ylabel = "ensemble members",
    title = "Marginal value of a Series-B raise",
)
histogram!(deal_launches; bins = -0.5:1:10.5, alpha = 0.5, label = "Series-B raise")
vline!([te.baseline, te.deal]; label = "means", lw = 2, color = :black, ls = :dash)

# ### Reading the result
#
# **The Series-B raise buys the reported +Δ expected launches over the three-year horizon**, with a
# standard error that is small relative to the effect — so it is a real gain, not sampling noise.
# That converts directly into a financing rule: the raise clears its hurdle when the value of those
# additional launches exceeds its dilution and cost.
#
# What matters is not the specific number but its *kind*. The gain is not "one added program advances
# to launch" arithmetic — it is a *system* effect. The extra capital relaxes the `@conserved`
# constraint that the priority allocator was rationing (§4), so programs already in the portfolio
# advance sooner and clear more gates within the horizon, on top of the one program the raise added.
# A static spreadsheet that adds a standalone NPV cannot see that interaction; only a timed,
# stochastic, resource-aware model surfaces it. The applied [case studies](../case_studies/marginal_scientist.md)
# scale exactly this reasoning up to the shadow price of a scarce resource and the value of an
# in-licensing deal, where "value is not additive under contention" is the headline.

# ## Recap and what comes next
#
# You have, on one growing portfolio model:
#
# 1. defined a structured-token kind (`@register` / `@aagent`) and seeded a portfolio with the
#    declarative `population` marking (both the explicit list and the `PopulationEntry` count forms);
# 2. authored a phase-as-attribute lifecycle with `@select` / `@advance`, including a stochastic gate
#    that soft-retires failures while preserving token identity;
# 3. selected a value-qualified subset with a continuous `npv > θ` predicate;
# 4. toured the resource-modality truth table and watched the priority allocator ration a contended
#    pool, in isolation and inside the portfolio;
# 5. added an in-model `Rule` — a Series-B raise — as a management lever, not host patch code;
# 6. computed the marginal value of that lever across a seeded ensemble with `treatment_effect`,
#    ending on +Δ expected launches ± SE.
#
# Two deep-dives go further on the machinery touched here: the [serialization deep-dive](../deep_dives/serialization.md)
# shows how this whole model — place, pipeline, lever, and portfolio — becomes an eval-free JSON
# document that round-trips loss-free, and the [composition deep-dive](../deep_dives/composition.md)
# covers `@join` / `@compose` / `refine` for building a portfolio out of fragments and moving between
# granularities. The [expert tutorial](expert.md) then places the portfolio as a node in a larger
# heterogeneous system — coupled to sibling agents over wires, checkpointed with `dump_state` /
# `restore`, and read through the exec-map analysis layer.
