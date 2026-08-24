# =============================================================================
# ReactiveDynamics.jl — REFINEMENT & OPEN-PORT COMPOSITION TOUR (ADR 0009 / §11)
# =============================================================================
#
# The core_engine_tour demo shows the FLAT modeling vocabulary (§1–§9: transitions,
# resource modalities, the allocator) and bd_acquisition shows a full structured-
# token application. This tour is about the VERTICAL axis the flat contract was
# otherwise silent on: hierarchical REFINEMENT & open-port COMPOSITION (ADR 0009 /
# CONTRACT §11). It answers the maintainer's framing requirement directly — "a
# modeling framework suitable for business processes: compact/expressive
# definition, compositionality, and various levels of granularity with more
# refined dynamics possibly substituted."
#
# Run it:   julia --project=demo/refinement_tour demo/refinement_tour/refinement_tour.jl
# First time: julia --project=demo/refinement_tour -e 'using Pkg; Pkg.instantiate()'
#
# Every construct below is copied from the package's passing semantic tests
# (test/semantic/refinement_composition.jl) and the operators/refine.jl docstrings
# — this file invents no API. The setting is the same pharma R&D pipeline the other
# demos use (Discovery → Phase1 → Phase2 → Phase3 → Filed → Market), but taught
# through the refinement lens, on PLAIN counted places (no structured-token
# machinery) so the refinement mechanics are the star and the demo runs fast.
#
# The whole layer is AUTHORING-time and additive: every operation here produces a
# plain ReactionNetwork that constructs / serializes / simulates exactly like
# a hand-written flat model (it is FORBIDDEN on a live/stepping model — it
# reindexes). The enabling mechanism is the ADR-0003 Phase-2 ArcSpec FK-repoint:
# place identification is repointing an integer FK, not string surgery.

using ReactiveDynamics
using ReactiveDynamics: nrows, row_ids, placename, find_index, arcs, port_role
using Printf

const RD = ReactiveDynamics
banner(title) = (println(); println("="^74); println(title); println("="^74))

# A small helper used throughout §3: the STRUCTURAL SIGNATURE of a named transition
# — its cycletime, prob-of-success, and its arc rows read off the promoted
# ArcSpec table as (place NAME, side, multiplicity). Two transitions with equal
# signatures are structurally identical. Keying by transition NAME (not index) makes
# the comparison robust to the row-reordering that refinement performs.
function trans_signature(m, tname)
    ti = findfirst(i -> m[i, :transName] === tname, collect(row_ids(m, :T)))
    ti === nothing && return nothing
    rows = sort(
        [
            (string(placename(m, r.place)), r.side, r.multiplicity)
                for r in arcs(m) if r.trans == ti && r.place > 0
        ]
    )
    return (ct = m[ti, :transCycleTime], pos = m[ti, :transProbOfSuccess], arcs = rows)
end


# =============================================================================
banner("§1. The coarse portfolio in ONE block: @pipeline (§D)")
# =============================================================================
#
# The dominant business-process shape is a CHAIN of phases. `@pipeline` authors the
# whole chain in one block: each `From => To : (ct, pos)` edge expands to a `flow`-
# genesis routing transition (CONTRACT §2.8) that consumes the upstream phase as an
# upfront LHS (so it fires only once a token exists there — token-flow, not an
# independent Poisson clock) and produces the downstream phase, carrying that edge's
# per-phase cycletime and probability-of-success. Six phases, five edges, one block
# — instead of five hand-wired reactions with repeated @move boilerplate. The
# transition names come out as `flow_<From>_<To>`.

build_portfolio() = @pipeline Project begin
    Discovery => Phase1:(ct = 1.0, pos = 0.45)
    Phase1 => Phase2:(ct = 1.5, pos = 0.6)
    Phase2 => Phase3:(ct = 2.0, pos = 0.4)
    Phase3 => Filed:(ct = 3.0, pos = 0.65)
    Filed => Market:(ct = 1.0, pos = 0.9)
end

portfolio = build_portfolio()
RD.populate_arcs!(portfolio)   # promote the incidence table so we can read it

println("@pipeline expanded the phase chain into a flat ReactionNetwork:")
println("  place (phases)   : ", portfolio[:, :placeName])
println("  transitions        : ", [portfolio[i, :transName] for i in row_ids(portfolio, :T)])
println(
    "  parts              : ", nrows(portfolio, :S), " place, ",
    nrows(portfolio, :T), " transitions"
)
println("  per-edge (ct, pos) :")
for i in row_ids(portfolio, :T)
    @printf(
        "    %-24s ct=%.1f  pos=%.2f\n",
        portfolio[i, :transName], portfolio[i, :transCycleTime],
        portfolio[i, :transProbOfSuccess]
    )
end
# A flow transition consumes its upstream phase on the LHS (the §2.8 flow idiom).
p2ix = find_index(:Phase2, portfolio)
consumes_phase2 = any(r -> r.place == p2ix && r.side === :lhs, arcs(portfolio))
println(
    "  flow_Phase2_Phase3 consumes :Phase2 on its LHS? ", consumes_phase2,
    "  (token-gated genesis, §2.8)"
)


# =============================================================================
banner("§2. Reusable fragments + open ports: @process / @port / @compose (§A,§D,§E)")
# =============================================================================
#
# The other axis of compactness: a NAMED, PARAMETERIZED fragment you instantiate
# many times. `@process name(params…) = begin <reaction lines> end` is a fragment
# factory — inside the body you write ordinary reaction lines with the parameter
# names BARE (e.g. `inp`, `outp`, `ct`); at call time each parameter is substituted
# structurally into the reaction AST BEFORE parsing (eval-free, no `$`-interpolation).
#
# Fragments compose by DECLARED PORTS instead of "remember which names to @equalize".
# A port is a boundary place tagged with a role (§A): :input (consumed-from
# boundary), :output (produced-into boundary), :shared (identified by bare name),
# or the default :private (auto-namespaced). `@port net A => input  B => output`
# tags them (note the `=>` pairs). `@compose f1 f2 …` is `@join` PLUS automatic port
# matching: an :output port of one fragment is identified with a same-named :input
# port of another by the FK-repoint (not string surgery), :private places are
# namespaced per fragment, :shared stay bare.

@process phase_gate(inp, outp; ct, pos) = begin
    1.0, inp --> outp, name => gate, cycletime => ct, probability => pos
end

# Two instances of the SAME fragment, wired head-to-tail on the shared place `Lead`.
screening = phase_gate(:Screen, :Lead; ct = 0.5, pos = 0.85)
lead_opt = phase_gate(:Lead, :Candidate; ct = 0.7, pos = 0.8)
println("phase_gate(:Screen, :Lead; …) — one instance of the reusable fragment:")
println("  place   : ", screening[:, :placeName], "   transitions: ", nrows(screening, :T))
println("  (ct, pos) : ", (screening[1, :transCycleTime], screening[1, :transProbOfSuccess]))

# Tag the boundary: `Lead` is the output of screening and the input of lead_opt — the
# same-named port that @compose will identify. `Screen`/`Candidate` stay dangling
# ports (single-fragment), Lead is the join seam.
@port screening  Screen => input   Lead => output
@port lead_opt   Lead => input      Candidate => output
println(
    "port roles — screening: Screen=", port_role(screening, :Screen),
    " Lead=", port_role(screening, :Lead),
    " | lead_opt: Lead=", port_role(lead_opt, :Lead),
    " Candidate=", port_role(lead_opt, :Candidate)
)

chain = @compose screening lead_opt
RD.populate_arcs!(chain)
names_chain = chain[:, :placeName]
println("@compose screening lead_opt:")
println("  merged places : ", names_chain)
println(
    "  shared port `Lead` collapsed to ONE place? ",
    count(==(:Lead), names_chain) == 1, "  (FK-repoint, not two pools)"
)
println(
    "  private places namespaced per fragment (f1__Screen, f2__Candidate)? ",
    (:f1__Screen in names_chain) && (:f2__Candidate in names_chain)
)
println("  transitions preserved : ", nrows(chain, :T), " (1 + 1, none lost)")
# The promoted incidence table is FK-EXACT: every static FK resolves, and both
# transitions route through the single shared `Lead` index.
leadix = find_index(:Lead, chain)
through_lead = count(r -> r.place == leadix, arcs(chain))
all_fk_ok = all(r -> r.place == 0 || 1 <= r.place <= nrows(chain, :S), arcs(chain))
println("  every arc FK in range?  ", all_fk_ok)
println("  rows routed through `Lead`:  ", through_lead, "  (produced by screening, consumed by lead_opt)")


# =============================================================================
banner("§3. REFINE one transition — the multifidelity payoff (§B)  ★ headline")
# =============================================================================
#
# The core workflow: model the portfolio coarsely for fast what-if, then ZOOM into
# ONE bottleneck when you need fidelity there — and have the rest of the model not
# notice. Here we substitute a detailed Phase-2 sub-model (screening → lead-opt →
# tox → filing, each its own step with its own cycletime and PoS) for the single
# coarse `flow_Phase2_Phase3` transition.
#
# `refine(spec, transition, submodel; ports)` splices the sub-model into the named
# coarse transition in four authoring-time structural moves: (1) namespace the sub's
# :private places; (2) identify the sub's open ports with the parent boundary
# place via `ports` by FK-repoint; (3) append the sub's transitions + remaining
# places/params/obs/EVENTS; (4) drop the coarse transition. It is NON-mutating
# (refine = refine! on a deepcopy). Because move (2) leaves the BOUNDARY places
# (Phase2, Phase3) at their same indices/names, coarse and refined are PLUG-
# COMPATIBLE (Invariant 1): every OTHER transition is structurally untouched.

phase2_detail = @reaction_network begin
    1.0, p2_in  --> screen, name => screening, cycletime => 0.5, probability => 0.85
    1.0, screen --> leadopt, name => lead_opt, cycletime => 0.7, probability => 0.8
    1.0, leadopt --> tox, name => tox_study, cycletime => 0.5, probability => 0.85
    1.0, tox    --> p2_out, name => filing_prep, cycletime => 0.3, probability => 0.7
end
# The sub-model's boundary: `p2_in` is where it plugs onto the parent's Phase2,
# `p2_out` onto Phase3. Everything else (screen/leadopt/tox) is :private → namespaced.
RD.set_port_role!(phase2_detail, :p2_in => :input, :p2_out => :output)

# Record the coarse boundary indices BEFORE the splice (for the invariance check).
phase2_ix_before = find_index(:Phase2, portfolio)
phase3_ix_before = find_index(:Phase3, portfolio)

refined = refine(
    portfolio, :flow_Phase2_Phase3, phase2_detail;
    ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out)
)

tnames_before = [portfolio[i, :transName] for i in row_ids(portfolio, :T)]
tnames_after = [refined[i, :transName]   for i in row_ids(refined, :T)]
println("transitions BEFORE refine : ", tnames_before)
println("transitions AFTER  refine : ", tnames_after)
println("  coarse `flow_Phase2_Phase3` removed?      ", !(:flow_Phase2_Phase3 in tnames_after))
println("  sub-steps spliced in (namespaced …__sub__…)?")
for n in tnames_after
    occursin("__sub__", string(n)) && println("      ", n)
end
println(
    "  refine was NON-mutating (coarse still intact)? ",
    :flow_Phase2_Phase3 in tnames_before
)

# ── The plug-compatibility headline: the boundary is unchanged and every OTHER
#    transition is byte-for-byte structurally identical before and after. ──
println()
println("PLUG-COMPATIBILITY (Invariant 1):")
println("  boundary places keep their indices:")
println(
    "    Phase2 : ", phase2_ix_before, " → ", find_index(:Phase2, refined),
    "   Phase3 : ", phase3_ix_before, " → ", find_index(:Phase3, refined)
)
println(
    "    (and their names — Phase2, Phase3 still present: ",
    (:Phase2 in refined[:, :placeName]) && (:Phase3 in refined[:, :placeName]), ")"
)

untouched = [:flow_Discovery_Phase1, :flow_Phase1_Phase2, :flow_Phase3_Filed, :flow_Filed_Market]
println("  every OTHER transition is structurally IDENTICAL before/after:")
identical = [tn => (trans_signature(portfolio, tn) == trans_signature(refined, tn)) for tn in untouched]
for (tn, same) in identical
    @printf("    %-24s identical? %s\n", tn, same)
end
println(
    "  ⇒ all untouched transitions identical: ", all(last, identical),
    "  — the rest of the portfolio does not notice the zoom."
)

# The sub's PRIVATE places are namespaced (not leaked as bare names).
println(
    "  sub-private places namespaced (bare `screen` NOT present): ",
    !(:screen in refined[:, :placeName]),
    " ; namespaced form present: ",
    any(n -> occursin("__sub__screen", string(n)), refined[:, :placeName])
)


# =============================================================================
banner("§4. Advisory boundary check: refinement_diagnostics (§C)")
# =============================================================================
#
# Refinement does NOT claim the fine model is behaviorally EQUIVALENT to the coarse
# one — that would need a bisimulation the framework can't check. Instead
# `refinement_diagnostics(submodel, coarse_attrs; ports, tol)` returns a
# Vector{String} of ADVISORY warnings (empty = clean), where computable:
#   * port-balance: an :input port consumed by no sub-transition LHS (or an :output
#       produced by no RHS) is a dangling port — a likely modeling error;
#   * on a linear chain: coarse.cycletime ≈ Σ sub cycletimes, coarse.pos ≈ Π sub PoS.
# These are warnings the author can OVERRIDE — the refinement may legitimately change
# the dynamics (that is the point of zooming in). They make the granularity ladder
# auditable without overclaiming equivalence.

# The coarse `flow_Phase2_Phase3` claimed (ct=2.0, pos=0.40). Check the detailed
# sub-model we actually spliced in against those coarse attributes.
coarse_attrs = Dict(:transCycleTime => 2.0, :transProbOfSuccess => 0.4)
warns_ok = RD.refinement_diagnostics(phase2_detail, coarse_attrs)
println("refinement_diagnostics(phase2_detail, coarse (ct=2.0, pos=0.40)):")
Σct = 0.5 + 0.7 + 0.5 + 0.3           # = 2.0  ≈ coarse ct
Πpos = 0.85 * 0.8 * 0.85 * 0.7       # ≈ 0.40 ≈ coarse pos
@printf("  Σ sub cycletimes = %.2f (≈ coarse 2.0) ; Π sub PoS = %.3f (≈ coarse 0.40)\n", Σct, Πpos)
println("  warnings: ", isempty(warns_ok) ? "none — well-matched, SILENT" : warns_ok)

# Now a DELIBERATELY DRIFTED refinement: cycletimes summing to 6.0 and PoS 0.81,
# both far from the coarse (2.0, 0.40) → both aggregate checks fire.
drifted = @reaction_network begin
    1.0, p2_in --> mid, name => slow_a, cycletime => 3.0, probability => 0.9
    1.0, mid   --> p2_out, name => slow_b, cycletime => 3.0, probability => 0.9
end
RD.set_port_role!(drifted, :p2_in => :input, :p2_out => :output)
warns_drift = RD.refinement_diagnostics(drifted, coarse_attrs)
println("refinement_diagnostics(drifted sub (Σct=6.0, ΠPoS=0.81), same coarse):")
for w in warns_drift
    println("  ⚠ ", w)
end

# And a DANGLING PORT: a place declared :input but only ever PRODUCED (RHS) —
# a common wiring mistake the port-balance check catches.
dangling = @reaction_network begin
    1.0, feed --> shelf, name => stock
end
RD.set_port_role!(dangling, :shelf => :input)   # declared :input but only ever produced
println("refinement_diagnostics(sub with an :input port never consumed):")
for w in RD.refinement_diagnostics(dangling, Dict())
    println("  ⚠ ", w)
end
println("  (these are OVERRIDABLE advisories, not equivalence proofs — Invariant 6.)")


# =============================================================================
banner("§5. Round-tripping the ladder: abstract (§B) and JSON (Invariant 5)")
# =============================================================================
#
# `abstract_transitions(spec, [t1,…], :into; lhs, rhs, attrs)` is the INVERSE of
# refine: collapse a connected set of sub-transitions back into a single coarse
# transition whose boundary reaction line is `lhs --> rhs`, carrying summarized
# attrs. It is a structural convenience for climbing back UP the granularity ladder.
# (Honest scope: it drops the sub-transition rows and adds the coarse one; it does
# NOT garbage-collect the now-orphaned internal place — those rows remain, inert.)

sub_transitions = [n for n in tnames_after if occursin("__sub__", string(n))]
collapsed = RD.abstract_transitions(
    refined, sub_transitions, :flow_Phase2_Phase3;
    lhs = [:Phase2], rhs = [:Phase3],
    attrs = Dict(:transCycleTime => 2.0, :transProbOfSuccess => 0.4)
)
println(
    "abstract_transitions collapsed the ", length(sub_transitions),
    " sub-steps back into one coarse `flow_Phase2_Phase3`:"
)
println("  transitions : ", [collapsed[i, :transName] for i in row_ids(collapsed, :T)])
println(
    "  coarse `flow_Phase2_Phase3` restored? ",
    :flow_Phase2_Phase3 in [collapsed[i, :transName] for i in row_ids(collapsed, :T)],
    " ; T count back to ", nrows(collapsed, :T),
    " (coarse was ", nrows(portfolio, :T), ")"
)

# Invariant 5 — a refined spec serializes / reloads as a FLAT model: refinement left
# no runtime trace, it is a plain ModelSpec.
@prob_params refined
json = RD.to_json_model(refined; meta = Dict{String, Any}("tspan" => 5.0))
reloaded = RD.build_network_from_dict(RD.JSON.parse(json))
println("JSON round-trip of the refined model (Invariant 5 — no runtime trace):")
println(
    "  place : ", nrows(refined, :S), " → reload ", nrows(reloaded, :S),
    "   transitions : ", nrows(refined, :T), " → reload ", nrows(reloaded, :T)
)
println(
    "  same place set after reload? ",
    Set(reloaded[:, :placeName]) == Set(refined[:, :placeName])
)


# =============================================================================
banner("§6. It's just a ModelSpec: construct + simulate the refined pipeline")
# =============================================================================
#
# The whole point of Invariant 4 (closure) + additivity: the refined model is an
# ordinary ReactionNetwork, so it constructs and simulates exactly like a
# hand-authored flat model — no special runtime for the refined structure. We seed a
# batch of Discovery projects and run the multifidelity pipeline to the horizon,
# reproducibly from (model, seed). (Macro args are LITERAL — evaluated in module
# scope — so @prob_init / @prob_meta take literal counts, the tests' idiom.)

const SEED = 20260710
function build_sim_model()
    m = refine(
        build_portfolio(), :flow_Phase2_Phase3, phase2_detail;
        ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out)
    )
    @prob_init m Discovery = 60
    @prob_params m
    @prob_meta m tspan = 30 dt = 1.0
    return m
end

prob = ReactionNetworkProblem(build_sim_model(); seed = SEED)
simulate(prob)
println("Constructed + simulated the REFINED pipeline (seed = $SEED):")
println("  solution columns : ", names(prob.sol))
launched = prob.sol[end, "Market"]
horizon = prob.sol[end, "t"]
@printf("  60 Discovery projects → %.0f reached :Market by t=%.0f\n", launched, horizon)
println("  final marking (pools that ran through the detailed Phase-2 sub-net):")
for c in ["Discovery", "Phase2", "Phase3", "Filed", "Market"]
    @printf("    %-10s end = %5.1f\n", c, prob.sol[end, c])
end
# The detailed sub-net's private pools appear in the solution too — the zoom is live.
subcols = [c for c in names(prob.sol) if occursin("__sub__", c)]
println("  detailed Phase-2 sub-net pools present in prob.sol: ", subcols)


# =============================================================================
banner("§7. Recap — the granularity ladder ADR 0009 / CONTRACT §11 gives you")
# =============================================================================
println(
    """
      §1  @pipeline — a chain of phases authored in ONE block, each edge a `flow`
          routing transition carrying (ct, pos): the coarse portfolio, compactly.
      §2  @process / @port / @compose — a reusable parameterized fragment, instantiated
          twice and composed by DECLARED open ports (output↔input identified by FK-
          repoint, private namespaced, shared bare) — compositionality without
          remembering which names to @equalize.
      §3  ★ refine — substitute a finer sub-model for one coarse transition, PLUG-
          COMPATIBLY: the boundary places keep their indices/names, so every OTHER
          transition is structurally untouched. The portfolio doesn't notice Phase-2
          became four sub-steps. (Non-mutating; refine! is the in-place form.)
      §4  refinement_diagnostics — ADVISORY §C checks (dangling ports; Σct / ΠPoS drift
          on a linear chain) that make the ladder auditable, WITHOUT claiming behavioral
          equivalence (Invariant 6 — warnings you may override).
      §5  abstract — climb back UP the ladder (collapse sub-steps into one coarse
          transition); and JSON round-trip proving a refined spec reloads as a FLAT model
          (Invariant 5 — refinement leaves no runtime trace).
      §6  …and it's just a ModelSpec: the refined pipeline constructs + simulates like any
          hand-written flat model, reproducibly from (model, seed) (Invariant 4 closure).

      This is the maintainer's ask — "various levels of granularity with more refined
      dynamics possibly substituted" — realized as cheap, collision-safe, authoring-time
      structural operations over the ADR-0003 ArcSpec FK table. All of §11 is
      FORBIDDEN on a live/stepping model (it reindexes); it operates on a static
      ReactionNetwork only.
    """
)
