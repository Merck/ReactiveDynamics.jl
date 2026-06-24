# BD acquisition-impact demo — HOST Julia (never serialized; the ADR-0006 §B/§C boundary).
#
# Defines the `ProjectToken` structured-token kind (a pharma program, identity preserved across
# phases via phase-as-attribute, ADR 0008 §D) and the per-network registry entries the model
# references BY NAME (the acquisition lever's AddToken constructor). Also the coarse pipeline
# model builder and the scenario grid. This is the "one host file" half of the MVP (§6/§7).

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, register_structured_species!,
    add_structured_token!, Rule, Seq, SetSpecies, SetParams, AddToken, get_species
using Random, Distributions, DataFrames

# ── The ProjectToken kind (host Julia, ADR 0006 §B) ─────────────────────────────────────
# A pharma program: phase is the lifecycle ATTRIBUTE (ADR 0008 §D canonical), the rest are the
# descriptor fields the rNPV roll-up reads. Defined inside ReactiveDynamics scope so the engine's
# bind/advance machinery (which lives there) can see the type — the @register/@aagent idiom.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
        phase::Symbol
        npv_peak::Float64        # peak sales value if it reaches market
        pos_remaining::Float64   # cumulative probability of success from here to market
        cost_to_date::Float64    # capital sunk so far (reconstructed in post; see MVP finding D)
        therapeutic_area::Symbol
        acquired::Bool           # did this program enter via the acquisition?
        acq_time::Float64        # NaN if organic
    end

    using Random: randstring
    # Host constructor (positional fields after the @aagent-injected name/species/bound/past_bonds).
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
            Tuple{Symbol,Float64,ReactiveDynamics.Transition}[],
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

# The registry constructor the model's AddToken lever references BY NAME (ADR 0006 §C / 0010).
# Signature is (state, fields::Dict) — fields are the evaluated AddToken field exprs.
const PROJECT_REGISTRY = Dict{Symbol,Any}(
    :ProjectToken => (state, fields) -> ReactiveDynamics.ProjectToken(;
        phase = get(fields, :phase, :Phase2),
        npv_peak = get(fields, :npv_peak, 1200.0),
        pos_remaining = get(fields, :pos_remaining, 0.5),
        therapeutic_area = get(fields, :therapeutic_area, :onc),
        acquired = true,
        acq_time = get(fields, :acq_time, state.t),
    ),
)

# The phase ladder (phase-as-attribute, one Project kind). Per-phase PoS and cycle time are the
# pipeline parameters; the product of per-phase PoS is the program's overall probability of launch.
const PHASES = [:Discovery, :Phase1, :Phase2, :Phase3, :Filed, :Market]
const PHASE_POS = Dict(            # per-phase probability of success (advance to next)
    :Discovery => 0.45,
    :Phase1 => 0.6,
    :Phase2 => 0.4,
    :Phase3 => 0.65,
    :Filed => 0.9,
)
const PHASE_CT = Dict(             # per-phase cycle time (ticks/years)
    :Discovery => 1.0,
    :Phase1 => 1.5,
    :Phase2 => 2.0,
    :Phase3 => 3.0,
    :Filed => 1.0,
)

# ── The coarse pipeline model (one @select/@advance transition per phase boundary) ──────
# Each advance: select a Project in phase N (with enough remaining PoS), consume scientists
# (@conserved, returned at finish) and burn budget (@rate), and on Binomial(q, PoS) success
# @advance(phase, N+1); failures soft-retire to :removed (finish! default). Authored directly
# here (the @pipeline sugar is ADR 0009 / a later stage); this is the expanded form.
function build_pipeline_model(; synergy_pos = 0, synergy_eff = 0)
    # Each advance fires at most a few instances/tick (genesis is token-gated by @select, so a
    # modest rate suffices); resource demands and financing are calibrated so the organic pipeline
    # flows to launches under contention without gridlocking (later-phase priority wins scientists).
    acs = @ReactionNetworkSchema begin
        # Discovery -> Phase1
        @deterministic(2.0),
        @select(Project, phase == :Discovery) + 2 * @conserved(scientist) + 2 * @rate(budget) -->
        @advance(phase, :Phase1),
        name => adv_discovery, cycletime => 1.0, probability => 0.45, priority => 1.0
        # Phase1 -> Phase2
        @deterministic(2.0),
        @select(Project, phase == :Phase1) + 3 * @conserved(scientist) + 3 * @rate(budget) -->
        @advance(phase, :Phase2),
        name => adv_phase1, cycletime => 1.5, probability => 0.6, priority => 1.5
        # Phase2 -> Phase3  (capability/PoS synergy raises PoS; op-efficiency shortens cycletime —
        # both param-mediated, MVP §2.1: the rule flips synergy_pos/synergy_eff at acquisition)
        @deterministic(2.0),
        @select(Project, phase == :Phase2) + 4 * @conserved(scientist) + 5 * @rate(budget) -->
        @advance(phase, :Phase3),
        name => adv_phase2, cycletime => 2.0 - 0.5 * synergy_eff,
        probability => 0.4 + 0.2 * synergy_pos, priority => 2.0
        # Phase3 -> Filed
        @deterministic(2.0),
        @select(Project, phase == :Phase3) + 5 * @conserved(scientist) + 8 * @rate(budget) -->
        @advance(phase, :Filed),
        name => adv_phase3, cycletime => 3.0 - 1.0 * synergy_eff,
        probability => 0.65 + 0.15 * synergy_pos, priority => 3.0
        # Filed -> Market
        @deterministic(2.0),
        @select(Project, phase == :Filed) + 1 * @conserved(scientist) + 2 * @rate(budget) -->
        @advance(phase, :Market),
        name => adv_filed, cycletime => 1.0, probability => 0.9, priority => 4.0
        # budget replenishment (financing): a steady inflow each tick
        @deterministic(30.0), ∅ --> budget, name => financing
    end

    register_structured_species!(acs, :Project)
    # Declare the resource pools and the synergy params. (@prob_params/@prob_meta eval their RHS
    # in module scope, so synergy values are set via set_params! with the function args, and
    # tspan/dt are passed to the constructor as kwargs.)
    @prob_init acs scientist = 40 budget = 200
    @prob_params acs synergy_pos = 0 synergy_eff = 0
    ReactiveDynamics.set_params!(acs, Dict(:synergy_pos => synergy_pos, :synergy_eff => synergy_eff))
    return acs
end

# ── Initial portfolio — DECLARATIVE initial marking (ADR 0007 §B, the "list of structures" form) ──
# The starting pipeline is a declarative `population[]` value passed to the constructor, so a
# structured run is reproducible from (model, population, seed) (MVP finding H, §8.2 S2) — not
# built by imperative post-construction host code. Here it is an explicit list of host token
# structs (the maintainer's "instantiate as a list of structures"); the count+attribute-dist
# PopulationEntry form is the alternative for "N programs with sampled NPV".
const ORGANIC_PORTFOLIO = [
    (:Discovery, 800.0, :onc),
    (:Discovery, 600.0, :immuno),
    (:Phase1, 1000.0, :onc),
    (:Phase1, 900.0, :cns),
    (:Phase2, 1500.0, :onc),
    (:Phase2, 1200.0, :immuno),
    (:Phase3, 2000.0, :onc),
]

function initial_population()
    pos_from_here() = prod(get(PHASE_POS, p, 1.0) for p in PHASES if p != :Market)
    return [
        ReactiveDynamics.ProjectToken(;
            phase = ph,
            npv_peak = npv,
            pos_remaining = pos_from_here(),
            therapeutic_area = area,
            acquired = false,
        ) for (ph, npv, area) in ORGANIC_PORTFOLIO
    ]
end

# ── The acquisition lever (endogenous Rule, ADR 0010) ───────────────────────────────────
# Fires once at t > T_acq: injects M acquired Phase-2 programs and (synergy 2) bumps scientists,
# (synergies 3/4) flips the synergy_pos/synergy_eff params. Built as a typed Rule the driver arms.
function acquisition_rule(; T_acq = 4.0, n_programs = 3, extra_scientists = 12,
        synergy_pos = false, synergy_eff = false)
    actions = ReactiveDynamics.ActionStmt[]
    for _ = 1:n_programs
        push!(
            actions,
            AddToken(:ProjectToken, [
                :phase => QuoteNode(:Phase2),
                :npv_peak => 1400.0,
                :pos_remaining => 0.4 * 0.65 * 0.9,
            ]),
        )
    end
    extra_scientists > 0 && push!(actions, SetSpecies(:scientist, extra_scientists, :inc))
    (synergy_pos || synergy_eff) && push!(
        actions,
        SetParams([:synergy_pos => (synergy_pos ? 1 : 0), :synergy_eff => (synergy_eff ? 1 : 0)]),
    )
    return Rule(:acquisition, :(@t() > $T_acq), Seq(actions); fire_mode = :once)
end
