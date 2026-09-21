#
# host.jl — the ONE code seam of this study (AGENT_WORKFLOW_CHARTER §3.2, ADR 0006 §B/§C).
#
# Everything about this study that is MODEL lives in the `.rdj.json` documents: places, transitions,
# arcs, modalities, rates, the starting portfolio, the horizon. This file carries the one thing a
# document cannot carry — the Julia TYPE of a structured token — and maps its kind name to a
# constructor. The document then references the kind BY NAME (`"kind": "Project"`) and stays
# eval-free: no Julia source crosses the JSON boundary in either direction.
#
# Read it as the exception that proves the rule. `rdx` prints a notice whenever it includes a host
# file, because this is the only point in the workflow where authored code executes.

using ReactiveDynamics
const RD = ReactiveDynamics

# The token kind must be defined in ReactiveDynamics' own scope: the selection/advancement machinery
# dispatches on the type, so `@register` evaluates the `@aagent` declaration there (the idiom the
# marginal-scientist case study uses, `docs/literate/case_studies/marginal_scientist.jl`).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct StudyProjectToken
        program::Symbol      # :Onc or :Immuno — which franchise this candidate belongs to
        phase::Symbol        # :Active while in trial, :Launched once the trial succeeds
        value::Float64       # a descriptor field, logged per tick
    end

    # A deterministic name counter, NOT `rand`: token names must not consume a random stream, or the
    # run stops being a function of `(model, seed)` alone (CONTRACT §4).
    const STUDY_TOKEN_SEQ = Ref(0)

    function StudyProjectToken(program, phase, value)
        STUDY_TOKEN_SEQ[] += 1
        return StudyProjectToken(
            "P" * string(STUDY_TOKEN_SEQ[]),                        # unique, reproducible token name
            :Project,                                              # kind tag
            nothing,                                               # bound_firing (engine-set)
            Tuple{Symbol, Float64, ReactiveDynamics.Firing}[],      # past_bonds history
            program, phase, value,
        )
    end
end

# Opt this kind into the per-tick trajectory log (ADR 0013): these fields are recorded each tick,
# RNG-free, so `token_trajectory` can reconstruct a candidate's path through the portfolio.
RD.log_token_fields(t::RD.StudyProjectToken) =
    (; program = t.program, phase = t.phase, value = t.value)

# The registry `rdx` looks for. Keys are the kind names the JSON documents reference; values are
# `(state, fields::Dict) -> token` constructors. Field defaults make a `@structured(...)` genesis
# node that omits a field still well-defined.
const RDX_REGISTRY = Dict{Symbol, Any}(
    :Project => (state, f) -> RD.StudyProjectToken(
        get(f, :program, :Onc), get(f, :phase, :Active), get(f, :value, 100.0)
    ),
)
