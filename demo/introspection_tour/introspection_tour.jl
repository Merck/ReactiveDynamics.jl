# =============================================================================
# ReactiveDynamics.jl — INTROSPECTION & EXEC-MAP TOUR (Phase-0.6 analysis/viz)
# =============================================================================
#
# The core_engine_tour demo shows how to BUILD and RUN models; the bd_acquisition
# demo shows a full worked application. This tour is about what you do AFTER a run
# has finished: the READ-ONLY analysis, ensemble, export, and visualization layer
# (ADR 0013 / 0014, CONTRACT §14 / §15). None of it touches the §1–§9 dynamics —
# every piece here is pure post-processing over a finished (or in-progress) run.
#
# Run it:   julia --project=demo/introspection_tour demo/introspection_tour/introspection_tour.jl
# First time: julia --project=demo/introspection_tour -e 'using Pkg; Pkg.instantiate()'
#
# Every construct below is copied from the package's passing semantic tests
# (test/semantic/analysis_observability.jl, visualization.jl, program_ledger.jl)
# — this file invents no API. The tour uses ONE small self-contained structured-
# token model (a two-phase project "advance") with an explicit seed=, and each
# section is preceded by the analysis IDEA and followed by narrated `println`
# output, so running it tells a story.
#
# The headline is §7: the three-layer network "exec map" — a Petri-net diagram of
# the model, decorated with run statistics (starved pools) and a @select cohort's
# path through the net, rendered to an image via Graphviz. That is the maintainer's
# "a system exec map readable for (in)efficiencies, decorated with results" ask.
#
# We need Plots (to render / for the recipe types) and Arrow (for the export
# siblings). The main project demotes both to weakdeps, so this demo carries its
# own Project.toml adding them (which also triggers RDPlotsExt / RDArrowExt).

using ReactiveDynamics
using DataFrames
using Statistics
using Printf
using Plots       # triggers RDPlotsExt — the §15.1 recipes render
using Arrow       # triggers RDArrowExt — the §14.3 export bundle writes .arrow siblings

RD = ReactiveDynamics
const HERE = @__DIR__
const OUTDIR = joinpath(HERE, "output")

banner(title) = (println(); println("="^74); println(title); println("="^74))


# =============================================================================
# §0. One small self-contained model: a two-phase "project advance"
# =============================================================================
#
# The whole tour runs on this single model — small on purpose, so it runs in well
# under a minute. It is the STRUCTURED-token regime (agentic tokens with identity
# and attributes), because the trajectory log, the per-program ledger, and the
# token-path highlighting all key off structured tokens.
#
# A `Project` is a structured token with a `phase` (a Symbol lifecycle attribute)
# and a numeric `value` (so the trajectory log records BOTH a Symbol path and a
# Real path — the Real one drives the medoid / envelope numeric metric). One
# transition `adv` @selects a Phase1 project, meters `budget` at @rate over its
# cycletime, and @advances the project to Phase2 with `probability => 0.6` (a
# Binomial success draw, so the run is genuinely stochastic — the ensemble has
# real spread). `budget` carries a placeCost so the burn is a real ledger cost, and
# the produced :Project carries a placeReward so a successful advance realizes
# reward. `budget` starts DELIBERATELY SCARCE (8) so the pool runs to a trough of
# 0 — which is exactly what the exec map paints as starvation. So this one model
# makes the trajectory log, the exec map's starvation coloring, AND both sides of
# the per-program ledger all non-trivial.
#
# This mirrors test/semantic/analysis_observability.jl::traj_model +
# program_ledger.jl::advance_cost_model exactly — no new API. (Macro arguments are
# LITERAL, so `probability => 0.6` is a schema literal; cost / reward / budget are
# set on the ACSet by index assignment, which is what the kwargs control.)

# The token kind is defined INTO the ReactiveDynamics module via @register (so it
# is referenced as RD.TrajProjectToken). A BaseStructuredToken carries the
# protocol fields (name, species, bound_transition, past_bonds); this kind adds
# `phase` and `value`.
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct TrajProjectToken
        phase::Symbol
        value::Float64
    end
    function TrajProjectToken(phase, value)
        return TrajProjectToken(
            "TP" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            value,
        )
    end
end

# The trajectory-log OPT-IN for this kind (ADR 0013 §A1): record phase + value each
# tick. A kind that does not define this contributes NO rows (opt-in is bounded).
RD.log_token_fields(t::RD.TrajProjectToken) = (; phase = t.phase, value = t.value)

# The model. Macro arguments are literal (evaluated in module scope), so cost /
# reward / budget are set on the ACSet directly via index assignment.
function project_model(; budget0 = 8, cost = 1.0, reward = 10.0)
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase1) + 2 * @rate(budget) --> @advance(phase, :Phase2),
            name => adv, cycletime => 1.0, probability => 0.6
    end
    RD.register_structured_species!(net, :Project)
    bi = findfirst(==(:budget), net[:, :placeName])
    net[bi, :placeInitVal] = Float64(budget0)
    net[bi, :placeCost] = cost
    pi = findfirst(==(:Project), net[:, :placeName])
    net[pi, :placeReward] = reward
    @prob_meta net tspan = 5 dt = 1.0
    return net
end

# Build + simulate one member under a seed. The starting portfolio is a handful of
# Phase1 projects with different values, passed declaratively as `population`.
build_prob(
    seed; pop = [
        RD.TrajProjectToken(:Phase1, 1.0),
        RD.TrajProjectToken(:Phase1, 2.0),
        RD.TrajProjectToken(:Phase1, 3.0),
    ]
) =
    ReactionNetworkProblem(project_model(); seed = seed, population = pop)

const SEED = 20260709
prob = build_prob(SEED)
simulate(prob)

println("Built + simulated the project-advance model (seed = $SEED).")
println("  solution columns : ", names(prob.sol), "  (construction order ⇒ index by name)")
println("  horizon          : tspan = ", prob.tspan, ", dt = ", prob.dt)


# =============================================================================
# §1. The raw run artifacts: prob.sol, prob.log, and the per-program ledger
# =============================================================================
#
# Three artifacts hang off a finished ReactionNetworkProblem and are the substrate
# everything else reads:
#
#   * prob.sol — the solution DataFrame: a "t" column plus one column per species,
#       the pool marking at each tick.
#   * prob.log — the heterogeneous event stream: a vector of tuples
#       (tag::Symbol, t, payload…). Read it by tag (r[1]). Financial rows are
#       :valuation_cost / :valuation_reward / :valuation; a per-tick per-program
#       snapshot is :program_ledger.
#   * program_ledger(prob) — the per-structured-token ledger DataFrame: one row per
#       program with cost_incurred / reward_realized / valuation / net, plus the
#       append-only audit trail via program_ledger_entries(prob, name).

banner("§1. Raw run artifacts: prob.sol, prob.log, the per-program ledger")

# prob.sol — read the budget pool's trajectory (it burns down as projects advance).
budget = prob.sol[!, "budget"]
println("prob.sol — budget pool over time : ", round.(budget; digits = 1))
println(
    "  budget start → end             : ", round(budget[1]; digits = 1),
    " → ", round(budget[end]; digits = 1), "  (burned by the @rate metering)"
)

# prob.log — reduce the ledger rows by tag (the core_engine_tour idiom).
agg_cost = sum(r[3] for r in prob.log if r[1] == :valuation_cost; init = 0.0)
agg_reward = sum(r[3] for r in prob.log if r[1] == :valuation_reward; init = 0.0)
log_tags = unique(r[1] for r in prob.log)
println("prob.log — ", length(prob.log), " event rows; distinct tags: ", log_tags)
println("  aggregate cost (Σ :valuation_cost)   : ", round(agg_cost; digits = 2))
println("  aggregate reward (Σ :valuation_reward): ", round(agg_reward; digits = 2))

# program_ledger — the per-program attribution DataFrame.
led = program_ledger(prob)
println("program_ledger(prob) — one row per program (", nrow(led), " rows):")
for row in eachrow(led)
    @printf(
        "  %-14s  species=%-8s  cost=%5.2f  reward=%5.2f  net=%6.2f\n",
        row.program, row.species, row.cost_incurred, row.reward_realized, row.net
    )
end
# The append-only audit trail for the first program: (t, kind, amount, transition).
if nrow(led) > 0
    first_prog = led[1, :program]
    entries = program_ledger_entries(prob, first_prog)
    println("program_ledger_entries(prob, \"", first_prog, "\") — the audit trail:")
    for e in entries
        @printf("    t=%.1f  %-7s  amount=%.2f  via=%s\n", e[1], e[2], e[3], e[4])
    end
end


# =============================================================================
# §2. The per-token trajectory log (ADR 0013 §14.1)
# =============================================================================
#
# When a token KIND opts in via log_token_fields, the engine records that token's
# opted-in fields at every tick — an 𝓕ₜ-measurable, RNG-free, read-only snapshot
# pushed at the same deterministic seam (and token order) as the per-program
# ledger. Read it in long form with token_trajectory:
#
#   * token_trajectory(prob)         — all rows: columns t, program, species, then
#       each opted-in field (here phase + value); one row per (tick, opted-in token).
#   * token_trajectory(prob, name)   — one token's life.
#   * token_trajectory(prob, pred)   — the rows of tokens CURRENTLY matching a
#       @select predicate (the same selection machinery the dynamics use).

banner("§2. Per-token trajectory log: token_trajectory + predicate scope")

traj = token_trajectory(prob)
println("token_trajectory(prob) — ", nrow(traj), " rows; columns: ", names(traj))
println("  first few rows:")
show(stdout, first(traj, min(6, nrow(traj))); summary = false)
println()

# One token's life (the first program in the log).
one_name = traj[1, :program]
one = token_trajectory(prob, one_name)
println(
    "token_trajectory(prob, \"", one_name, "\") — that one token's ", nrow(one),
    " logged ticks; phases: ", one.phase
)

# Predicate scope: the tokens that ended in :Phase2 (they advanced). TokenPredicate
# reuses the @select machinery (a conjunctive clause list); a QuoteNode wraps the
# literal symbol being compared against.
adv_pred = RD.TokenPredicate(:Project, [RD.Clause(:phase, :(==), QuoteNode(:Phase2))])
adv_rows = token_trajectory(prob, adv_pred)
println(
    "token_trajectory(prob, @select Phase2) — the advanced cohort: ",
    length(unique(adv_rows.program)), " token(s), ", nrow(adv_rows), " rows"
)


# =============================================================================
# §3. "Typical" helpers: representative_token (medoid) + trajectory_envelope
# =============================================================================
#
# Two precise senses of "a typical program", both over the trajectory store:
#
#   * representative_token(prob) — the MEDOID: the token whose numeric logged path
#       is closest to the per-cohort mean path (a documented heuristic). Returns a
#       token name you can feed back into token_trajectory(prob, name).
#   * trajectory_envelope(prob) — the ENVELOPE: for each numeric logged field, the
#       per-tick median and inter-quartile band (q25/q75) across the cohort.
#       Columns: align, field, median, q25, q75, n. This is the single-run input to
#       the ensemble band (ADR 0014 recipe 6).

banner("§3. Typical-program helpers: representative_token + trajectory_envelope")

rep = representative_token(prob)
println("representative_token(prob) — the medoid program : ", rep)
rep_life = token_trajectory(prob, rep)
println("  its `value` path : ", rep_life.value)

env = trajectory_envelope(prob)
println("trajectory_envelope(prob) — columns: ", names(env))
val_band = env[env.field .== :value, :]
println("  the numeric `value` field's median ± IQR band by tick:")
for row in eachrow(val_band)
    @printf(
        "    t=%.1f  median=%.2f  [q25=%.2f, q75=%.2f]  n=%d\n",
        row.align, row.median, row.q25, row.q75, row.n
    )
end


# =============================================================================
# §4. Ensemble analysis: ensemble / summarize / treatment_effect (§14.2)
# =============================================================================
#
# A single run is one sample; the honest signal comes from an ENSEMBLE. `ensemble`
# runs `nseed` independent members, member k seeded deterministically from
# hash((root_seed, k)) (§4 D8) — so member k is reproducible from (root_seed, k)
# regardless of how many members you run or in what order (§4 D9). Each member is
# rebuilt fresh per seed (mode :rebuild). The result is an EnsembleProblem, itself
# a readable AlgebraicAgents hierarchy node.
#
#   * summarize(ens, metric) — reduce a per-run scalar metric across members to
#       (; mean, sem, q25, median, q75, n).
#   * treatment_effect(base, deal, metric) — the unpaired A/B difference of means
#       with se = √(var_b/n_b + var_d/n_d).
#
# Keep the ensemble SMALL (a handful of members) so the demo stays fast.

banner("§4. Ensemble analysis: ensemble / summarize / treatment_effect")

# A small ensemble: each member builds + simulates one project-advance run.
ens = ensemble(s -> (p = build_prob(s); simulate(p); p); nseed = 6, root_seed = 2026)
println("ensemble(build; nseed=6, root_seed=2026)")
println("  members            : ", length(ens.members))
println("  per-member seeds   : ", ens.seeds, "  (= hash((2026, k)) for k in 1:6)")
println("  run mode           : ", ens.mode, "  (:rebuild — mode a ships today)")

# summarize a per-run metric: the reward realized across a member's programs (a
# Binomial success draw makes this genuinely vary run to run — the spread `sem`
# exists to report). The metric is any `member -> Real`; here it reads the ledger.
realized_reward(p) = sum(program_ledger(p).reward_realized)
s = summarize(ens, realized_reward)
@printf(
    "summarize(ens, realized reward) : mean=%.2f ± sem=%.2f   median=%.2f  [q25=%.2f, q75=%.2f]  n=%d\n",
    s.mean, s.sem, s.median, s.q25, s.q75, s.n
)

# EnsembleProblem is an AA-readable node: it exports cross-run observables (the
# across-member mean of its members' observables) and holds each member as a child.
obs = AlgebraicAgents.observables(ens)
println("AlgebraicAgents.observables(ens) : ", obs)
if :budget in obs
    println(
        "  getobservable(ens, :budget)    : ", round(AlgebraicAgents.getobservable(ens, :budget); digits = 2),
        "  (the across-member mean)"
    )
end
println("  AA children (members)          : ", length(AlgebraicAgents.inners(ens)))

# treatment_effect: an A/B lever. The "deal" arm injects more starting capital
# (budget0 20 vs the scarce baseline 8), so it leaves more budget at the horizon.
base = ensemble(
    s -> (
        p = build_prob(s; pop = [RD.TrajProjectToken(:Phase1, 1.0)]);
        simulate(p); p
    ); nseed = 6, root_seed = 9
)
deal = ensemble(
    s -> (
        p = ReactionNetworkProblem(
            project_model(; budget0 = 20); seed = s,
            population = [RD.TrajProjectToken(:Phase1, 1.0)]
        ); simulate(p); p
    ); nseed = 6, root_seed = 9
)
te = treatment_effect(base, deal, p -> last(p.sol.budget))
@printf("treatment_effect(baseline budget0=8, deal budget0=20; final budget):\n")
@printf(
    "  baseline=%.2f  deal=%.2f  Δ=%.2f ± se=%.2f  (n_b=%d, n_d=%d)\n",
    te.baseline, te.deal, te.delta, te.se, te.n_baseline, te.n_deal
)
println("  ⇒ the better-capitalized deal arm leaves more budget (Δ ≥ 0), the A/B signal.")


# =============================================================================
# §5. Export bundles: export_run / export_ensemble (§14.3)
# =============================================================================
#
# A finished run (or ensemble) serializes to a directory bundle — format per
# artifact (one format does not fit all): rectangular artifacts (the trajectory,
# the ledger, the token log) as CSV, heterogeneous ones (the event stream, per-
# token histories, the manifest) as JSON. Because THIS demo loads Arrow, each
# rectangular artifact also gets a byte-faithful `.arrow` sibling (RDArrowExt).
# The run.json manifest pins model_hash + seed, so the bundle is traceable to a
# replayable (model, seed) pair (Invariant 5).

banner("§5. Export bundles: export_run / export_ensemble")

run_dir = joinpath(OUTDIR, "run")
ispath(run_dir) && rm(run_dir; recursive = true)
export_run(prob, run_dir)
println("export_run(prob, \"", relpath(run_dir, HERE), "\") wrote:")
for f in sort(readdir(run_dir))
    println("  ", f)
end
println(
    "  (.arrow siblings present ⇒ Arrow is loaded: RD._arrow_available() = ",
    RD._arrow_available(), ")"
)

ens_dir = joinpath(OUTDIR, "ensemble")
ispath(ens_dir) && rm(ens_dir; recursive = true)
export_ensemble(ens, ens_dir; metric = realized_reward)
println("export_ensemble(ens, \"", relpath(ens_dir, HERE), "\"; metric=realized reward) wrote:")
for f in sort(readdir(ens_dir))
    tag = isdir(joinpath(ens_dir, f)) ? "/  (a per-member export_run bundle)" : ""
    println("  ", f, tag)
end


# =============================================================================
# §6. Result-plot recipes (§15.1) — model-agnostic wrappers over raw artifacts
# =============================================================================
#
# The seven result plots are thin wrapper TYPES (constructible in the core, with
# Plots absent) that Plots.jl @recipes dispatch on when Plots is loaded (RDPlotsExt).
# Each wraps a raw artifact (prob.sol / prob.log / program_ledger / the trajectory
# log / an ensemble), so `plot(MarkingPlot(prob))` works on ANY model. Because this
# demo loads Plots, we actually render a couple to PNG.

banner("§6. Result-plot recipes: MarkingPlot / LedgerPlot / TokenTrajectoryPlot …")

# The wrapper types (constructible regardless of Plots).
println(
    "recipe wrapper types constructed: ",
    (
        RD.MarkingPlot(prob), RD.SaturationPlot(prob), RD.ValuationPlot(prob),
        RD.LedgerPlot(prob), RD.ThroughputPlot(prob),
    ) .|> typeof .|> nameof
)

# Render two to PNG (RDPlotsExt is live since Plots is loaded).
marking_png = joinpath(OUTDIR, "marking.png")
Plots.savefig(Plots.plot(RD.MarkingPlot(prob)), marking_png)
println("  MarkingPlot rendered → ", relpath(marking_png, HERE))

traj_png = joinpath(OUTDIR, "token_value_trajectory.png")
Plots.savefig(Plots.plot(RD.TokenTrajectoryPlot(prob, :value)), traj_png)
println("  TokenTrajectoryPlot(:value) rendered → ", relpath(traj_png, HERE))


# =============================================================================
# §7. THE EXEC MAP — the three-layer network diagram (ADR 0014 §15.2)  ★ headline
# =============================================================================
#
# The maintainer's headline ask: "a system exec map readable for (in)efficiencies,
# decorated with simulation results." It comes in three SEPARABLE layers, each
# usable alone:
#
#   Layer A — network_graph(prob) → NetworkGraph. A pure, dependency-free Petri-net
#       view: species (place) nodes, transition nodes, arcs with stoichiometry +
#       modality. Runs on a deepcopy, so it does NOT perturb the caller's RNG — pure.
#   Layer B — to_graphviz(g) emits DOT; draw_network(prob) renders it through AA's
#       run_graphviz (Graphviz_jll or a system `dot`). No run needed — this is the
#       structure diagram (authoring-time documentation).
#   Layer C — exec_map(prob; highlight) DECORATES Layer A with run statistics:
#       species nodes filled gold where their pool ran to a trough (starvation), and
#       a @select cohort's past_bonds path through the net drawn as thick arcs
#       ("where did these programs go?"). Read-only — never mutates state.
#
# Rendering is best-effort: if no Graphviz backend is installed we still emit the
# DOT source and say so.

banner("§7. The exec map — network_graph → draw_network → exec_map")

# Layer A — the pure structure. Confirm it did not perturb the RNG.
g = network_graph(prob)
println("Layer A network_graph(prob):")
println("  species (places)   : ", [s.name for s in g.species])
println("  transitions        : ", [t.name for t in g.transitions])
println("  arcs               : ", length(g.arcs), " (:in LHS→T and :out T→RHS)")

# Layer B — DOT source (always available, no backend needed) + a render attempt.
dot = to_graphviz(g)
dot_path = joinpath(OUTDIR, "network.dot")
write(dot_path, dot)
println(
    "Layer B to_graphviz(g): wrote DOT source → ", relpath(dot_path, HERE),
    " (", length(dot), " bytes)"
)

# Is a Graphviz backend available? Try to render Layer B to a file. draw_network
# with a `path` writes the rendered output there and returns the path.
network_svg = joinpath(OUTDIR, "network.svg")
graphviz_ok = false
try
    global graphviz_ok
    draw_network(prob; format = "svg", path = network_svg)
    graphviz_ok = isfile(network_svg) && filesize(network_svg) > 0
catch err
    @warn "draw_network: no Graphviz backend available — falling back to DOT" exception = err
end
if graphviz_ok
    println(
        "  draw_network rendered → ", relpath(network_svg, HERE),
        " (", filesize(network_svg), " bytes)"
    )
else
    println(
        "  no Graphviz backend — the structure lives in ", relpath(dot_path, HERE),
        " (render it with `dot -Tsvg`)"
    )
end

# Layer C — the decorated exec map. Highlight the advanced (:Phase2) cohort's path.
# First inspect the overlay DOT directly (always available): starvation fill +
# highlighted arcs. Then render if a backend exists.
starved = [s for (s, v) in RD._pool_troughs(prob) if v <= 0.0]
println("Layer C exec_map — decoration inputs:")
println("  starved species (pool trough ≤ 0) : ", isempty(starved) ? "none" : starved)

# Build the highlight arc set the way exec_map does: each matching token's
# past_bonds map a (species, transition-index) to the SAME node id the graph uses.
hi_arcs = Tuple{Symbol, Symbol}[]
for tok in RD.select_tokens(prob, adv_pred)
    for (sp, _t, tr) in tok.past_bonds
        push!(hi_arcs, (sp, RD._transition_node_name(prob.network, tr.i)))
    end
end
println("  @select(Phase2) token-path arcs   : ", isempty(hi_arcs) ? "none (no cohort bonds)" : unique(hi_arcs))

# Confirm exec_map is read-only, then render (or fall back to DOT).
sol_before = copy(prob.sol)
exec_svg = joinpath(OUTDIR, "exec_map.svg")
exec_ok = false
try
    global exec_ok
    exec_map(prob; highlight = adv_pred, format = "svg", path = exec_svg)
    exec_ok = isfile(exec_svg) && filesize(exec_svg) > 0
catch err
    @warn "exec_map: no Graphviz backend available — falling back to DOT" exception = err
end
# The overlay DOT is always obtainable, backend or not (it is what exec_map renders).
overlay_dot = to_graphviz(g; highlight_species = starved, highlight_arcs = hi_arcs)
exec_dot_path = joinpath(OUTDIR, "exec_map.dot")
write(exec_dot_path, overlay_dot)
@assert prob.sol == sol_before "exec_map must be read-only (Invariant 3)"
println("  exec_map read-only check (prob.sol unchanged) : ", prob.sol == sol_before)
if exec_ok
    println("  exec_map rendered → ", relpath(exec_svg, HERE), " (", filesize(exec_svg), " bytes)")
    println("  ★ Open it: the gold-filled place is a starved pool; the thick arcs are the")
    println("    advanced cohort's path through the net — the (in)efficiency read.")
else
    println(
        "  no Graphviz backend — the decorated exec map DOT is in ",
        relpath(exec_dot_path, HERE), " (render it with `dot -Tsvg`)"
    )
end


# =============================================================================
# §8. Recap — what this tour exercised
# =============================================================================
banner("§8. Recap — the Phase-0.6 analysis & visualization surface")
println(
    """
      §1  Raw artifacts: prob.sol (marking DataFrame), prob.log (tagged event stream,
          reduced by tag to aggregate cost/reward), and program_ledger(prob) +
          program_ledger_entries (the per-program attribution DataFrame + audit trail).
      §2  The per-token trajectory log: token_trajectory(prob) (long form), one token's
          life by name, and predicate-scoped rows via a @select TokenPredicate.
      §3  "Typical" helpers: representative_token (the medoid program) and
          trajectory_envelope (per-tick median + IQR band over a numeric field).
      §4  Ensembles: ensemble (nseed members, hash((root,k)) seeding), summarize
          (mean/sem/quantiles), treatment_effect (unpaired A/B Δ), and the
          EnsembleProblem as an AA-readable node (observables / getobservable / inners).
      §5  Export bundles: export_run / export_ensemble — CSV + JSON core, plus the
          Arrow siblings (RDArrowExt, because this demo loads Arrow), manifest-pinned.
      §6  Result-plot recipes: the model-agnostic wrapper types (MarkingPlot, …) and
          two rendered to PNG via RDPlotsExt.
      §7  ★ The exec map, three layers: network_graph (pure structure) → draw_network
          (Graphviz render of the Petri net) → exec_map (decorated with starvation
          coloring + a @select cohort's token-path highlighting) — the (in)efficiency view.

      All of the above is READ-ONLY post-processing over a finished run — it never
      touches the dynamics. Artifacts were written under demo/introspection_tour/output/.
    """
)
