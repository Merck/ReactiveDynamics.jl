# # What is this in-licensing asset worth to *this* pipeline?
#
# **The answer, up front.** A late-stage asset we can in-license is not worth a number we look up in
# isolation and add to our balance sheet. Its worth is a *system* effect on the programs we already
# own: it competes for the same scientists and the same cash, and — if the target brings a platform —
# it raises the odds on our organic programs too. When we run our real pipeline once *without* the deal
# and once *with* it, over a seeded ensemble, the fully-integrated deal is clearly value-accretive: it
# adds on the order of **+1,900 value units of risk-adjusted NPV, net of the acquisition price**, at a
# standard error several times smaller than the effect. The non-obvious part is the second-order read:
# **rNPV is not additive under contention.** The full deal is worth materially *more* than the sum of
# its synergies priced one at a time, and the single largest lever is the platform lifting the success
# probability of programs *already in our pipeline* — not the acquired programs themselves. The value is
# in the platform, not the pipeline.
#
# We build that answer here from a living portfolio: structured-token programs advancing through a
# timed, resource-contended pipeline, with the acquisition itself expressed as an in-model rule that
# fires once. Everything below runs at build time from an explicit seed. (We run a modest ensemble here
# so the page builds quickly; the full study — [`demo/bd_acquisition`](https://github.com/Merck/ReactiveDynamics.jl/tree/main/demo/bd_acquisition)
# — runs 160 seeds, which is what it takes to resolve the *fine* ranking of the individual synergies.)
#
# The underlying semantics — structured tokens, the endogenous decision channel, resource modalities,
# the per-program ledger — are specified in the normative [operational-semantics contract](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/CONTRACT_DRAFT.md); here we use them.

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, register_structured_species!,
    Rule, Seq, SetSpecies, SetParams, AddToken, get_species, inners, getagent, program_ledger
using Random, Distributions, DataFrames, Statistics, Printf
using Plots

const RD = ReactiveDynamics

# ## 1. The pipeline, the asset, and the acquisition lever
#
# **A program is a structured token.** Each asset carries identity across its whole lifecycle: its
# `phase` (Discovery → Phase1 → Phase2 → Phase3 → Filed → Market) is an *attribute*, not a separate
# species, so one `ProjectToken` kind represents every program and we advance it in place. The other
# fields are what the valuation reads — peak sales value, remaining probability of success, and whether
# it arrived organically or via the deal. The kind is host Julia (never serialized): we define it into
# the `ReactiveDynamics` module with the `@register` / `@aagent` idiom, so the engine's bind/advance
# machinery can see the type, and refer to it afterwards as `RD.ProjectToken`.

@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
        phase::Symbol
        npv_peak::Float64        # peak sales value if it reaches market
        pos_remaining::Float64   # cumulative probability of success from here to market
        cost_to_date::Float64    # capital sunk so far
        therapeutic_area::Symbol
        acquired::Bool           # did this program enter via the acquisition?
        acq_time::Float64        # NaN if organic
    end

    using Random: randstring
    function ProjectToken(;
            phase = :Discovery,
            npv_peak = 1000.0,
            pos_remaining = 0.1,
            therapeutic_area = :onc,
            acquired = false,
            acq_time = NaN,
        )
        return ProjectToken(
            "Proj" * randstring(8),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            npv_peak,
            pos_remaining,
            0.0,
            therapeutic_area,
            acquired,
            acq_time,
        )
    end
end

# The acquisition lever injects new programs *by name*: the model document references a constructor in a
# per-network registry, so the file itself carries no code. Acquired programs enter at Phase 2, on the
# same risk scale as our organic ones.

const PROJECT_REGISTRY = Dict{Symbol, Any}(
    :ProjectToken => (state, fields) -> RD.ProjectToken(;
        phase = get(fields, :phase, :Phase2),
        npv_peak = get(fields, :npv_peak, 1200.0),
        pos_remaining = get(fields, :pos_remaining, 0.5),
        therapeutic_area = get(fields, :therapeutic_area, :onc),
        acquired = true,
        acq_time = get(fields, :acq_time, state.t),
    ),
)

# The probability that a program eventually launches, *from* a given phase, is the product of the
# per-phase success probabilities still ahead of it. A Phase-3 program only has to clear Phase3 → Filed
# → Market, so it carries far more remaining value than a Discovery program — the risk scale both the
# organic portfolio and the acquired assets are placed on.

const PHASES = [:Discovery, :Phase1, :Phase2, :Phase3, :Filed, :Market]
const PHASE_POS = Dict(
    :Discovery => 0.45, :Phase1 => 0.6, :Phase2 => 0.4, :Phase3 => 0.65, :Filed => 0.9,
)
function pos_from_phase(ph)
    idx = findfirst(==(ph), PHASES)
    idx === nothing && return 1.0
    return prod(get(PHASE_POS, p, 1.0) for p in PHASES[idx:end])
end

# **The pipeline itself.** One `@select`/`@advance` transition per phase boundary: select an in-phase
# program, hold `scientist` headcount (`@conserved`, returned on completion) and burn `budget`
# (`@rate`), and on a `Binomial` success draw advance the program to the next phase; failures soft-retire
# it. A steady financing inflow tops up the cash reserve each tick. Two synergy parameters — a
# capability lift on late-phase success probability and an operational-efficiency cut to late-phase cycle
# time — sit in the transition expressions at zero, waiting for the deal to switch them on. The pools are
# calibrated so cash genuinely binds: the organic pipeline runs `budget` down to its floor, so the
# priority allocator is rationing a scarce resource rather than idling on a slack model. That is what
# makes the deal's value a contention effect and not a bookkeeping sum.

function build_pipeline_model(; synergy_pos = 0, synergy_eff = 0)
    net = @reaction_network begin
        @deterministic(2.0),
            @select(Project, phase == :Discovery) + 2 * @conserved(scientist) + 2 * @rate(budget) -->
            @advance(phase, :Phase1),
            name => adv_discovery, cycletime => 1.0, probability => 0.45, priority => 1.0
        @deterministic(2.0),
            @select(Project, phase == :Phase1) + 3 * @conserved(scientist) + 3 * @rate(budget) -->
            @advance(phase, :Phase2),
            name => adv_phase1, cycletime => 1.5, probability => 0.6, priority => 1.5
        @deterministic(2.0),
            @select(Project, phase == :Phase2) + 4 * @conserved(scientist) + 5 * @rate(budget) -->
            @advance(phase, :Phase3),
            name => adv_phase2, cycletime => 2.0 - 0.5 * synergy_eff,
            probability => 0.4 + 0.2 * synergy_pos, priority => 2.0
        @deterministic(2.0),
            @select(Project, phase == :Phase3) + 5 * @conserved(scientist) + 8 * @rate(budget) -->
            @advance(phase, :Filed),
            name => adv_phase3, cycletime => 3.0 - 1.0 * synergy_eff,
            probability => 0.65 + 0.15 * synergy_pos, priority => 3.0
        @deterministic(2.0),
            @select(Project, phase == :Filed) + 1 * @conserved(scientist) + 2 * @rate(budget) -->
            @advance(phase, :Market),
            name => adv_filed, cycletime => 1.0, probability => 0.9, priority => 4.0
        @deterministic(16.0), ∅ --> budget, name => financing
    end

    register_structured_species!(net, :Project)
    @prob_init net scientist = 40 budget = 150
    @prob_params net synergy_pos = 0 synergy_eff = 0
    RD.set_params!(net, Dict(:synergy_pos => synergy_pos, :synergy_eff => synergy_eff))

    ## Price the budget burn at 1 currency/unit so the engine's per-program ledger attributes each
    ## program's capital spend during the run. This touches only the ledger rows, never the dynamics —
    ## the trajectory, the pools, and the Δ-rNPV are exactly as they would be without it.
    bi = findfirst(==(:budget), net[:, :specName])
    net[bi, :specCost] = 1.0
    return net
end

# **The starting portfolio** is declarative: an explicit list of `ProjectToken`s across phases, passed
# to the constructor and instantiated before t = 0, so a structured run is reproducible from
# `(model, population, seed)`.

const ORGANIC_PORTFOLIO = [
    (:Discovery, 800.0, :onc),
    (:Discovery, 600.0, :immuno),
    (:Phase1, 1000.0, :onc),
    (:Phase1, 900.0, :cns),
    (:Phase2, 1500.0, :onc),
    (:Phase2, 1200.0, :immuno),
    (:Phase3, 2000.0, :onc),
]
initial_population() = [
    RD.ProjectToken(;
            phase = ph, npv_peak = npv, pos_remaining = pos_from_phase(ph),
            therapeutic_area = area, acquired = false,
        ) for (ph, npv, area) in ORGANIC_PORTFOLIO
]

# **The acquisition is a rule, in the model.** It fires *once*, at the first tick past `T_acq`, and its
# action is a sequence: inject the acquired Phase-2 programs, optionally bump the resource pools (the
# target brings headcount and capital), and optionally flip the capability / efficiency synergy params.
# A "scenario" is simply which of those levers the rule arms — nothing about the pipeline changes. This
# is the endogenous decision channel: the deal is serializable model data, not host patch code.

function acquisition_rule(;
        T_acq = 8.0, n_programs = 3, extra_scientists = 0, extra_budget = 0,
        synergy_pos = false, synergy_eff = false
    )
    actions = RD.ActionStmt[]
    for _ in 1:n_programs
        push!(
            actions,
            AddToken(
                :ProjectToken, [
                    :phase => QuoteNode(:Phase2),
                    :npv_peak => 1400.0,
                    :pos_remaining => pos_from_phase(:Phase2),
                ]
            ),
        )
    end
    extra_scientists > 0 && push!(actions, SetSpecies(:scientist, extra_scientists, :inc))
    extra_budget > 0 && push!(actions, SetSpecies(:budget, extra_budget, :inc))
    (synergy_pos || synergy_eff) && push!(
        actions,
        SetParams([:synergy_pos => (synergy_pos ? 1 : 0), :synergy_eff => (synergy_eff ? 1 : 0)]),
    )
    return Rule(:acquisition, :(@t() > $T_acq), Seq(actions); fire_mode = :once)
end

# One scenario, one seed: build the pipeline, seed the run, and — unless this is the no-deal baseline —
# arm the acquisition rule with the synergies that scenario turns on. The six scenarios are the standard
# BD grid: no deal (`S0`); the deal with programs only (`S1`); plus resource (`S2`), plus capability/PoS
# (`S3`), plus operational-efficiency (`S4`), and the fully-integrated deal with all four (`S5`).

function run_scenario(scenario::Symbol, seed; tspan = 40.0, T_acq = 8.0)
    net = build_pipeline_model()
    prob = ReactionNetworkProblem(
        net; tspan = tspan, dt = 1.0, seed = seed,
        registry = PROJECT_REGISTRY, population = initial_population(),
    )
    if scenario != :S0
        res = scenario in (:S2, :S5)
        pos = scenario in (:S3, :S5)
        eff = scenario in (:S4, :S5)
        push!(
            prob.rules,
            acquisition_rule(;
                T_acq = T_acq, n_programs = 3,
                extra_scientists = res ? 15 : 0, extra_budget = res ? 250 : 0,
                synergy_pos = pos, synergy_eff = eff,
            ),
        )
    end
    simulate(prob)
    return prob
end

# ## 2. Reading portfolio value off a finished run
#
# The engine does no discounting — that is a modeling choice we make in post, over the final token
# population. A program that reached market has realized its peak value; an in-flight program is worth
# its remaining probability of success times its peak, discounted for the expected years still to
# launch. Portfolio rNPV is the sum, net of what we paid for the deal.

const YEARS_TO_MARKET = Dict(
    :Discovery => 9.0, :Phase1 => 7.0, :Phase2 => 5.0, :Phase3 => 3.0, :Filed => 1.0, :Market => 0.0,
)
tokens(prob) = collect(values(inners(getagent(prob, "structured"))))
is_active(t) = get_species(t) != :removed
reached_market(t) = t.phase == :Market

function portfolio_rnpv(prob; discount = 0.1, acq_price = 0.0)
    gross = 0.0
    for t in tokens(prob)
        is_active(t) || continue
        if reached_market(t)
            gross += t.npv_peak
        else
            ttm = get(YEARS_TO_MARKET, t.phase, 8.0)
            gross += t.pos_remaining * t.npv_peak / (1 + discount)^ttm
        end
    end
    return gross - acq_price
end

# The per-run scalars we compare across scenarios: rNPV (netting the price), launches, whether at least
# one program reached market, and — the tell for contention — the *low-water mark* of the cash reserve,
# the lowest `budget` ever fell to. A value near its floor means cash was binding.

m_rnpv(acq_price) = prob -> portfolio_rnpv(prob; acq_price = acq_price)
m_launches(prob) = count(reached_market, tokens(prob))
m_p_launch(prob) = m_launches(prob) >= 1 ? 1.0 : 0.0
m_cash_trough(prob) = minimum(prob.sol[!, "budget"])

# ## 3. With and without the deal: the attributable Δ-rNPV
#
# The counterfactual is the whole method. We hold the pipeline, the seeds, and the starting portfolio
# fixed, and run each scenario as a seeded ensemble: member `k` derives its seed from a single root seed
# so the runs are reproducible and independent across scenarios. `treatment_effect` then reports the
# difference of ensemble means with its unpaired standard error. (We keep the ensemble modest — 40 seeds
# — so the page builds in about a minute; the full study runs 160, which is what tightens the SE enough
# to rank the individual synergies.)

const ACQ_PRICE = 400.0
const NSEED = 40
const ROOT_SEED = 2026

scenario_ensemble(scenario) = ensemble(
    seed -> run_scenario(scenario, seed); nseed = NSEED, root_seed = ROOT_SEED,
)

scenarios = [:S0, :S1, :S2, :S3, :S4, :S5]
ens = Dict(s => scenario_ensemble(s) for s in scenarios)
base = ens[:S0]

# The deal scenarios carry the acquisition price; the baseline does not. Since the price is a
# deterministic constant offset, we net it out of the rNPV difference and the standard error is
# unchanged. The headline is the fully-integrated deal, `S5`, against the no-deal baseline, `S0`:

price(s) = s == :S0 ? 0.0 : ACQ_PRICE
function delta_rnpv(deal_ens; deal_price = ACQ_PRICE)
    te = treatment_effect(base, deal_ens, m_rnpv(0.0))   # gross Δ + unpaired SE
    return (delta = te.delta - deal_price, se = te.se, baseline = te.baseline, deal = te.deal - deal_price)
end

headline = delta_rnpv(ens[:S5])
launch_te = treatment_effect(base, ens[:S5], m_launches)

@printf("Baseline portfolio rNPV (no deal)   : %.0f value units\n", headline.baseline)
@printf("With the fully-integrated deal      : %.0f  (net of the %.0f price)\n", headline.deal, ACQ_PRICE)
@printf("Attributable Δ-rNPV (S5 − S0)        : %+.0f ± %.0f  (1 SE, %d seeds)\n", headline.delta, headline.se, NSEED)
@printf("Δ-launches (extra programs to market): %+.2f\n", launch_te.delta)

# The whole deal-scenario rNPV distribution sits to the right of the no-deal one — the deal shifts the
# *distribution*, not just its mean:

rnpv_samples(e; acq_price = 0.0) = Float64[m_rnpv(acq_price)(m) for m in e.members]
s0_rnpv = rnpv_samples(base)
s5_rnpv = rnpv_samples(ens[:S5]; acq_price = ACQ_PRICE)

histogram(
    s0_rnpv; bins = 15, alpha = 0.5, label = "no deal (S0)",
    xlabel = "portfolio rNPV (value units)", ylabel = "ensemble members",
    title = "In-licensing deal: rNPV with vs without ($NSEED seeds)",
)
histogram!(s5_rnpv; bins = 15, alpha = 0.5, label = "full deal (S5)")
vline!([mean(s0_rnpv), mean(s5_rnpv)]; label = "means", lw = 2, color = :black, ls = :dash)

# ## 4. Why the value is not additive
#
# A spreadsheet prices a deal as a standalone rNPV and adds it to the portfolio. The living model says
# that is wrong on two counts, and both are visible in the scenario grid. First the grid itself — the
# ensemble mean rNPV, the launches, and the cash low-water mark for every scenario:

@printf("%-26s %10s %9s %8s   %s\n", "scenario", "mean rNPV", "launches", "cash⌄", "Δ-rNPV vs S0 (±SE)")
println("-"^82)
labels = Dict(
    :S0 => "baseline (no deal)", :S1 => "deal, programs only", :S2 => "+ resource synergy",
    :S3 => "+ capability/PoS", :S4 => "+ op-efficiency", :S5 => "full (all synergies)",
)
for s in scenarios
    e = ens[s]
    d = delta_rnpv(e; deal_price = price(s))
    @printf(
        "%-26s %10.0f %9.2f %8.1f   %s\n",
        labels[s], summarize(e, m_rnpv(price(s))).mean, summarize(e, m_launches).mean,
        summarize(e, m_cash_trough).mean,
        s == :S0 ? "—" : @sprintf("%+8.0f ± %5.0f", d.delta, d.se),
    )
end

# **First: cash is the binding constraint, in every scenario.** The `cash⌄` column pins to the same
# floor whether or not we do the deal — the organic pipeline alone already runs the reserve down to the
# bottom. The starting reserve is 150; the low-water mark sits far below it throughout. That is the
# system fact a standalone rNPV cannot see: the acquired programs do not run in a vacuum, they draw on a
# reserve that is *already* exhausted, so the deal's value is entirely about what it does to a contended
# system — including whether it eases or worsens the squeeze.
#
# **Second: the synergies are super-additive.** Price each synergy on its own (its marginal Δ over the
# programs-only deal `S1`) and sum them; then compare to the fully-integrated deal `S5`, which turns all
# of them on together. They do not match — and the gap is real, not rounding:

s1 = delta_rnpv(ens[:S1]).delta
marginals = Dict(
    :S2 => delta_rnpv(ens[:S2]).delta - s1,
    :S3 => delta_rnpv(ens[:S3]).delta - s1,
    :S4 => delta_rnpv(ens[:S4]).delta - s1,
)
sum_of_parts = s1 + sum(values(marginals))
full_deal = delta_rnpv(ens[:S5]).delta
interaction = full_deal - sum_of_parts

@printf("programs-only deal (S1)              : %+.0f\n", s1)
@printf("  + resource synergy, marginal       : %+.0f\n", marginals[:S2])
@printf("  + capability/PoS synergy, marginal : %+.0f\n", marginals[:S3])
@printf("  + op-efficiency synergy, marginal  : %+.0f\n", marginals[:S4])
println("-"^52)
@printf("sum of the parts priced standalone   : %+.0f\n", sum_of_parts)
@printf("fully-integrated deal (S5)           : %+.0f\n", full_deal)
@printf("interaction (S5 − sum of parts)      : %+.0f\n", interaction)

# The interaction term is positive: integrating the synergies is worth more than buying them one at a
# time. The mechanism is exactly the contention we just saw — capital and capability compound. More
# capital lets more programs run in parallel; a capability lift pushes more of them through to launch and
# *out* of the resource pool, which frees capacity for the rest. Neither lever, priced alone against a
# cash-starved pipeline, captures what they do together. The synergy stack, with the interaction bar
# that a sum-of-parts valuation would miss:

bar(
    ["programs\n(S1)", "resource", "capability", "op-effic.", "interaction"],
    [s1, marginals[:S2], marginals[:S3], marginals[:S4], interaction];
    legend = false, ylabel = "marginal Δ-rNPV (value units)",
    title = "Deal value is super-additive under contention",
    color = [:steelblue, :steelblue, :steelblue, :steelblue, :goldenrod],
)

# On which synergy is *largest*, this reduced ensemble is honest about its limits: with 40 seeds the
# per-synergy standard error (~±300) is wider than the gaps between the individual synergies, so their
# ranking here is within noise. The full 160-seed study resolves it — and the ranking is the headline
# read for a BD partner: the single largest synergy is **capability/PoS, the target's platform lifting
# the success probability of the programs the company already owns** (marginal ≈ +550), ahead of the
# acquired programs themselves. The value is in the platform, not the pipeline.

# ## 5. Where the capital went — and the decision
#
# The engine attributes its cost ledger *per program* during the run, so we can see which programs the
# scarce capital actually flowed to. On one representative fully-integrated run, sorted by capital
# burned:

demo_prob = run_scenario(:S5, hash((ROOT_SEED, 1)))
led = program_ledger(demo_prob)
tok_by_name = Dict(RD.AlgebraicAgents.getname(t) => t for t in tokens(demo_prob))
led.phase = [haskey(tok_by_name, n) ? tok_by_name[n].phase : :removed for n in led.program]
led.acquired = [haskey(tok_by_name, n) ? tok_by_name[n].acquired : false for n in led.program]
led.npv_peak = [haskey(tok_by_name, n) ? tok_by_name[n].npv_peak : NaN for n in led.program]
sort!(led, :cost_incurred; rev = true)

println("  phase       acquired   npv_peak   capital burned")
for r in eachrow(first(led, min(6, nrow(led))))
    @printf(
        "  %-10s  %-8s   %8.0f   %14.1f\n",
        string(r.phase), string(r.acquired), isnan(r.npv_peak) ? 0.0 : r.npv_peak, r.cost_incurred,
    )
end

launched_capital = sum(led.cost_incurred[led.phase .== :Market]; init = 0.0)
retired_capital = sum(led.cost_incurred[led.species .== :removed]; init = 0.0)
@printf("\n  capital that reached a launched program : %.1f\n", launched_capital)
@printf("  capital sunk into failed/retired programs: %.1f\n", retired_capital)

# ### Reading the result
#
# **This asset is worth about +1,900 value units of risk-adjusted NPV to *this* pipeline, net of a 400
# price** — a clearly value-accretive deal at that price, with a standard error several times below the
# effect. But the number that should change how the deal is priced is not the headline; it is the
# structure behind it. The deal's value is **not** a standalone rNPV we add on: cash is already binding
# before the first acquired program starts, so the asset's worth is entirely what it does to a contended
# system we already own — and the synergies *compound*, so the integrated deal is worth more than the
# sum of its parts. The largest single lever is the platform raising the odds on our organic pipeline,
# not the programs we buy.
#
# That converts directly into BD guidance. The price to beat is the *attributable* Δ-rNPV on the living
# portfolio, not a bolt-on rNPV of the target's programs — and it should be credited to the platform, so
# a deal that is marginal as a pipeline purchase can clear easily once the capability synergy on the
# existing portfolio is counted. The same discipline runs the reverse case: at a higher price, sweep it
# until the Δ crosses zero to find the walk-away number, and read the per-program ledger to see how much
# of the capital is recoverable if a program is killed at its next gate. Every one of those is a
# re-run of this same counterfactual with a different lever — which is the point of pricing a deal on a
# dynamic model rather than a spreadsheet.
#
# The neighboring [flagship case study](marginal_scientist.md) turns the same contention machinery on a
# capacity question — the shadow price of one more scientist — and [the kill-threshold study](kill_a_program.md)
# puts the decision rule itself inside the model.
