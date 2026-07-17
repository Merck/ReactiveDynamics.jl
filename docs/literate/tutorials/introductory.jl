# # Introductory tutorial — your first ReactiveDynamics model
#
# **What you will build.** A timed, stochastic R&D pipeline, end to end: author it in
# the modeling metalanguage, simulate it reproducibly, read the results by name, and
# — the point of the whole exercise — turn the run into *a number a manager acts on*.
#
# **Who this is for.** Anyone new to ReactiveDynamics. You need no prior exposure to
# discrete-event simulation or to the package. We stay entirely in the *classical*
# regime — every quantity is a plain counted stock (molecules, dollars, scientists,
# jobs). Structured/agentic tokens (projects with attributes and lifecycle identity)
# are the subject of the [advanced tutorial](advanced.md); here we learn the bare engine.
#
# **The mental model in one paragraph.** A model is a set of **transitions**. Each
# transition has a **rate** (how often it tries to fire), a left-hand side of
# **reactants** it consumes, and a right-hand side of **products** it emits. Firing can
# be instantaneous or take time (a `cycletime`, during which an in-flight instance may
# also fail a success draw). All randomness flows through a per-run **seeded RNG**, so a
# run is fully reproducible from `(model, seed)`. That is the whole engine.
#
# Every construct below is drawn from the package's passing semantic test suite — this
# tutorial invents no API.

using ReactiveDynamics
using Statistics                # mean / std / quantile for the closing decision number

# ## 1. A first model: the SIR epidemic
#
# We start with the textbook susceptible–infected–recovered epidemic, because it
# introduces every part of the authoring surface you will reuse for the rest of the
# tutorial. A model is written with the `@reaction_network` macro; each line reads
# `rate, LHS --> RHS, name => …`.
#
# - `S + I --> 2I` — a susceptible meets an infected and *becomes* infected (net −1 `S`,
#   +1 `I`). Its rate `α*S*I` is a **mass-action** expression: a bare numeric expression
#   is a *stochastic (Poisson)* intensity.
# - `I --> R` — an infected recovers, at rate `β*I`.
#
# No individual is created or destroyed outright, so `S + I + R` is a structural
# **invariant** — a built-in sanity check the engine must preserve exactly.

sir = @reaction_network begin
    α * S * I, S + I --> 2I, name => infection
    β * I, I --> R, name => recovery
end

# The model is authored; now we attach its numbers. Each of these companion macros
# takes **literal** right-hand sides (they evaluate in module scope, so a loop variable
# would not resolve — use literals):
#
# - `@prob_init` — the initial counts (the *marking* at t = 0),
# - `@prob_params` — the named parameters the rate expressions reference,
# - `@prob_meta` — the simulation horizon `tspan` and the time step `dt`.

@prob_init sir S = 999 I = 10 R = 0
@prob_params sir α = 0.0001 β = 0.01
@prob_meta sir tspan = 250 dt = 0.1

# `ReactionNetworkProblem(model; seed = …)` compiles the authored network into a runnable
# problem. **The `seed=` kwarg owns a per-run RNG — it is the only route to
# reproducibility** (§8 returns to this). `simulate(prob)` then advances to `tspan`; the
# solution lands in `prob.sol`, a `DataFrame` with a `"t"` column plus one column per
# species.

sir_prob = ReactionNetworkProblem(sir; seed = 1)
simulate(sir_prob)

# **Read solution columns by name.** Column order is *construction* order, not the order
# you wrote the species, so positional indexing is a foot-gun — always index by name.

S = sir_prob.sol[!, "S"]
I = sir_prob.sol[!, "I"]
R = sir_prob.sol[!, "R"]
total = S .+ I .+ R

println("Initial population S+I+R      : ", total[1])
println(
    "Invariant drift (max − min)   : ",
    round(maximum(total) - minimum(total); digits = 9), "  (≈ 0 ⇒ conserved)"
)
println("Epidemic peak |I|             : ", round(maximum(I); digits = 1), "  (started at ", I[1], ")")
println("Infected at the horizon       : ", round(I[end]; digits = 1), "  (declines after the peak)")

# The population is conserved to floating-point exactness, and `I` rises to a genuine
# interior peak before burning out — an outbreak, reproduced from `(sir, seed=1)`.

# ## 2. From epidemic to pipeline: timed, probabilistic transitions
#
# The SIR reactions completed in the same tick they fired. Real processes take **time**
# and can **fail**. We now model the object we actually care about: a three-phase R&D
# pipeline. Candidate programs enter *discovery*, are *screened* into preclinical,
# *advanced* into the clinic, and run a clinical *trial* to approval. Each phase is a
# plain counted pool; each transition carries lifecycle attributes:
#
# - `cycletime => c` — a fired instance stays *in-flight* and completes after `ceil(c/dt)`
#   ticks (nothing appears downstream before then).
# - `probability => p` — on completion, success is a `Binomial(q, p)` draw; only successes
#   emit the RHS. (alias: `prob`)
# - `capacity => k` — at most `k` instances may be in-flight at once; proposals beyond `k`
#   are *deferred* to later ticks, never dropped. This is the scarce clinical trial slot.
#
# The routing rates are high (`@deterministic(50.0)`) so that each phase moves whatever
# its upstream pool holds — the flow is **token-gated**, clamped to available programs,
# not to the nominal rate. That makes the long, capacity-limited clinical `trial` the
# real bottleneck, which is exactly the lever we will price in §4.
#
# We read one tick as one month, so `tspan = 60` is a five-year horizon (`dt = 1.0` ⇒ 60
# monthly ticks). Time units are whatever you choose; the engine only sees ticks.

pipeline = @reaction_network begin
    2.0, ∅ --> discovery, name => intake
    @deterministic(50.0), discovery --> preclinical,
        name => screen, cycletime => 1.0, probability => 0.6
    @deterministic(50.0), preclinical --> clinical,
        name => advance, cycletime => 1.0, probability => 0.7
    @deterministic(50.0), clinical --> approved,
        name => trial, cycletime => 5.0, probability => 0.5, capacity => 3
end
@prob_init pipeline discovery = 0 preclinical = 0 clinical = 0 approved = 0
@prob_params pipeline
@prob_meta pipeline tspan = 60 dt = 1.0

pipe_prob = ReactionNetworkProblem(pipeline; seed = 1)
simulate(pipe_prob)

approved = pipe_prob.sol[!, "approved"]
clinical = pipe_prob.sol[!, "clinical"]
println("Solution columns              : ", names(pipe_prob.sol), "  (index by name!)")
println("Programs approved by horizon  : ", Int(approved[end]))
println("Clinical queue depth (max)    : ", Int(maximum(clinical)), "  (programs waiting on a trial slot)")

# A queue builds in front of the clinical trial: more programs are ready than the three
# slots can run. That backlog is the signature of a **binding constraint** — and the
# reason the next slot might be worth paying for.

# ## 3. Reproducibility: the seed *is* the run
#
# Reproducibility comes from the `seed=` construction kwarg **alone**. The state owns its
# own RNG, isolated from Julia's global stream, so `Random.seed!` does *not* pin a run and
# a seeded run neither reads nor perturbs the global RNG. The contract:
#
# - same `(model, seed)` ⇒ byte-identical trajectory,
# - different seed ⇒ (almost surely) a different trajectory,
# - unseeded ⇒ a fresh entropy seed each construction.

run_once(seed) = (p = ReactionNetworkProblem(pipeline; seed); simulate(p); p.sol[!, "approved"][end])
println("seed 1 vs seed 1 identical    : ", run_once(1) == run_once(1))
println("seed 1 vs seed 2 differ       : ", run_once(1) != run_once(2))

# A single seeded run is one sample of a stochastic process — never report it as *the*
# answer. To get a distribution, run an **ensemble**: derive each member's seed
# deterministically from a single root seed plus the member index, `hash((root, k))`.
# Member `k` is then reproducible from `(root, k)` regardless of how many members you run
# or in what order — the basis for stable Monte-Carlo statistics.

member_seed(root, k) = hash((root, k))
function approvals(model, root, k)
    p = ReactionNetworkProblem(model; seed = member_seed(root, k))
    simulate(p)
    return p.sol[!, "approved"][end]
end

ens = [approvals(pipeline, 2026, k) for k in 1:200]
println("Ensemble of 200 runs, approvals by horizon:")
println("  mean ± std       : ", round(mean(ens); digits = 2), " ± ", round(std(ens); digits = 2))
println("  p10 / p90        : ", quantile(ens, 0.1), " / ", quantile(ens, 0.9))

# ## 4. The decision: is a fourth clinical trial slot worth it?
#
# Here is the discipline this framework is built for: a run is not a chart, it is an
# **input to a decision**. Our pipeline is cash-constrained at three concurrent clinical
# trials, and a backlog is queueing behind them (§2). *Does adding a fourth slot move the
# number a portfolio manager cares about — approvals over the horizon — enough to justify
# its cost?*
#
# We answer it as a **counterfactual**: the identical pipeline with `capacity => 4` on the
# clinical trial, run over the same 200 ensemble seeds, and compare the mean approvals.
# (Model attributes are literals, so the two capacities are two literal model builders —
# the honest, idiomatic way to hold everything else fixed.)

pipeline_4slots = @reaction_network begin
    2.0, ∅ --> discovery, name => intake
    @deterministic(50.0), discovery --> preclinical,
        name => screen, cycletime => 1.0, probability => 0.6
    @deterministic(50.0), preclinical --> clinical,
        name => advance, cycletime => 1.0, probability => 0.7
    @deterministic(50.0), clinical --> approved,
        name => trial, cycletime => 5.0, probability => 0.5, capacity => 4
end
@prob_init pipeline_4slots discovery = 0 preclinical = 0 clinical = 0 approved = 0
@prob_params pipeline_4slots
@prob_meta pipeline_4slots tspan = 60 dt = 1.0

ens_3 = ens                                                    # the capacity-3 baseline from §3
ens_4 = [approvals(pipeline_4slots, 2026, k) for k in 1:200]   # same seeds, one more slot

marginal = mean(ens_4) - mean(ens_3)
se = sqrt(var(ens_3) / length(ens_3) + var(ens_4) / length(ens_4))

println("Approvals over the horizon (200-seed ensemble):")
println("  3 clinical slots : ", round(mean(ens_3); digits = 2))
println("  4 clinical slots : ", round(mean(ens_4); digits = 2))
println("  marginal 4th slot: +", round(marginal; digits = 2), " approvals  (± ", round(se; digits = 2), " SE)")

# ### The number, and what it means
#
# **A fourth clinical trial slot buys roughly +4 additional approvals over the five-year
# horizon** — well outside the ensemble noise. The verdict a manager acts on: *fund the
# fourth slot if the incremental value of ~4 approvals exceeds the cost of standing it up.*
#
# The point is not the specific number; it is that the framework produced a *marginal,
# system-level* answer that a static spreadsheet cannot. The fourth slot's value comes
# from relieving a contended bottleneck, and only a timed, stochastic, resource-aware
# model surfaces that — the same machinery that, scaled up, prices the marginal scientist
# and values an in-licensing deal in the [applied case studies](../case_studies/marginal_scientist.md).

# ## Recap
#
# You have, end to end:
#
# 1. authored a model in the `@reaction_network` metalanguage and attached its numbers
#    with `@prob_init` / `@prob_params` / `@prob_meta`;
# 2. simulated it with `ReactionNetworkProblem(…; seed=)` + `simulate`, and read
#    `prob.sol` **by column name**;
# 3. seen mass-action (Poisson) vs `@deterministic` rates and the timed lifecycle
#    (`cycletime`, `probability`, `capacity`);
# 4. made a run *reproducible* via `seed=` and turned an ensemble into a **decision
#    number** — the marginal value of one more clinical trial slot.
#
# Next: the [advanced tutorial](advanced.md) replaces plain counted pools with *structured
# tokens* — programs that carry attributes and identity through their lifecycle — and adds
# resource modalities, the priority allocator, and in-model decision rules.
