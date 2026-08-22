# # Composition & granularity — the granularity ladder
#
# **What this covers.** How a model is assembled from smaller models, and how one coarse
# transition is *substituted* by a finer sub-model without the rest of the network noticing.
# This is the *vertical* axis the flat authoring surface (transitions, modalities, the
# allocator) is silent on: **compositionality** and **levels of granularity**. We author a
# pharma R&D pipeline coarsely for fast what-if, then zoom into one bottleneck when we need
# fidelity there — and we finish on the guarantee that makes the zoom trustworthy.
#
# **Who this is for.** Readers comfortable with the [introductory tutorial](../tutorials/introductory.md)
# who want the structural operators. We stay in the *classical* regime — every phase is a
# plain counted pool — so the refinement mechanics are the star and every block runs fast.
#
# **The one invariant to hold onto.** Every operator here is **authoring-time**: it rewrites
# a static `ReactionNetwork` before it is constructed into a runnable problem, and it
# **reindexes** the store as it does so. That makes the whole layer *forbidden on a live,
# stepping model* — you compose, refine, and abstract a network, and only then hand the
# result to `ReactionNetworkProblem`. The theory (open-port semantics, the FK-repoint that
# makes place-identification cheap, the plug-compatibility invariants) lives in the normative
# [operational-semantics contract §11](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/CONTRACT_DRAFT.md)
# and [ADR 0009](https://github.com/Merck/ReactiveDynamics.jl/blob/rework/spec/adr/0009-refinement-and-composition.md); here we exercise the operators.

using ReactiveDynamics
using Printf
using Plots                     # inline figures

# A handful of *internal* store accessors let us look inside a network to check our work.
# They are not part of the modeling surface — they read the typed struct-of-columns store
# (row ids, the promoted arc-incidence table) so we can assert an operator did what it
# claims. We import them explicitly to keep that boundary visible.
using ReactiveDynamics: nrows, row_ids, find_index, arcs, placename,
    populate_arcs!, port_role
const RD = ReactiveDynamics

# The **structural signature** of a named transition — its cycletime, probability-of-success,
# and the arc rows read off the promoted incidence table as `(place name, side,
# stoich)`. Two transitions with equal signatures are structurally identical; keying by
# transition *name* (not row index) makes the comparison robust to the row-reordering that
# composition and refinement perform. We use it to *prove* plug-compatibility in §4.
function trans_signature(m, tname)
    ti = findfirst(i -> m[i, :transName] === tname, collect(row_ids(m, :T)))
    ti === nothing && return nothing
    rows = sort(
        [
            (string(placename(m, r.place)), r.side, r.stoich)
                for r in arcs(m) if r.trans == ti && r.place > 0
        ]
    )
    return (ct = m[ti, :transCycleTime], pos = m[ti, :transProbOfSuccess], arcs = rows)
end

# ## 1. Manual composition: `@join` and `@equalize`
#
# The lowest rung. `@join` takes the **union** of two networks' place, transitions, and
# parameters (and their events and observables), optionally *identifying* shared places
# across the two via equations. `@equalize` collapses two places *within* one network into a
# single pool and rewrites every reference. Both are the **manual, no-declared-ports** path:
# you name the place to identify by hand. Both operate on a static network, before
# construction.
#
# Two reaction sub-systems each consume a shared resource `A`; we join them, identifying the
# two `A`s as one pool, and count the parts of the merged network.

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

# `@equalize` collapses two conceptually-identical place `A` and `A2` into one pool.

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

# `@join`/`@equalize` work, but they make you *remember which names to identify*. The next
# rung declares the boundary once, on the fragment, and lets composition match it.

# ## 2. Declared open ports: `@process` / `@port` / `@compose`
#
# A **fragment** is a named, parameterized model factory. `@process name(params…) = begin …
# end` writes ordinary reaction lines with the parameters *bare*; at call time each parameter
# is substituted structurally into the reaction AST before parsing (eval-free — there is no
# `$`-interpolation in the DSL).

@process phase_gate(inp, outp; ct, pos) = begin
    1.0, inp --> outp, name => gate, cycletime => ct, probability => pos
end

# We instantiate the *same* fragment twice, wired head-to-tail on the shared place `Lead`.

screening = phase_gate(:Screen, :Lead; ct = 0.5, pos = 0.85)
lead_opt = phase_gate(:Lead, :Candidate; ct = 0.7, pos = 0.8)
println("phase_gate(:Screen, :Lead; …) — one instance of the reusable fragment:")
println("  place   : ", screening[:, :placeName], "   transitions: ", nrows(screening, :T))
println("  (ct, pos) : ", (screening[1, :transCycleTime], screening[1, :transProbOfSuccess]))

# A **port** is a boundary place tagged with a role: `:input` (consumed-from boundary),
# `:output` (produced-into boundary), `:shared` (identified by bare name), or the default
# `:private` (auto-namespaced on compose). `@port` tags them via `place => role` pairs
# (written with `=>`, not `=`). Here `Lead` is the *output* of screening and the *input* of
# lead_opt — the same-named port `@compose` will identify; `Screen`/`Candidate` stay dangling.

@port screening  Screen => input   Lead => output
@port lead_opt   Lead => input      Candidate => output
println(
    "port roles — screening: Screen=", port_role(screening, :Screen),
    " Lead=", port_role(screening, :Lead),
    " | lead_opt: Lead=", port_role(lead_opt, :Lead),
    " Candidate=", port_role(lead_opt, :Candidate)
)

# `@compose` is `@join` **plus automatic port matching**: an `:output` port of one fragment
# is identified with a same-named `:input` port of another by repointing an integer foreign
# key (not string surgery), `:private` places are namespaced per fragment, and `:shared`
# places stay bare.

chain = @compose screening lead_opt
populate_arcs!(chain)   # promote the incidence table so we can read it
names_chain = chain[:, :placeName]
println("@compose screening lead_opt:")
println("  merged places : ", names_chain)
println(
    "  shared port `Lead` collapsed to ONE place? ",
    count(==(:Lead), names_chain) == 1, "  (FK-repoint, not two pools)"
)
leadix = find_index(:Lead, chain)
through_lead = count(r -> r.place == leadix, arcs(chain))
println(
    "  rows routed through `Lead` : ", through_lead,
    "  (produced by screening, consumed by lead_opt ⇒ one seam, not two)"
)
println("  transitions preserved      : ", nrows(chain, :T), " (1 + 1, none lost)")

# ## 3. The coarse portfolio in ONE block: `@pipeline`
#
# The dominant business-process shape is a **chain of phases**. `@pipeline` authors the whole
# chain in one block: each `From => To : (ct, pos)` edge expands to a `flow` routing
# transition that consumes the upstream phase as an upfront left-hand side — so it fires only
# once a token exists there (token-flow, not an independent Poisson clock) — and produces the
# downstream phase, carrying that edge's cycletime and probability-of-success. Six phases,
# five edges, one block; the transition names come out as `flow_<From>_<To>`.

build_portfolio() = @pipeline Project begin
    Discovery => Phase1:(ct = 1.0, pos = 0.45)
    Phase1 => Phase2:(ct = 1.5, pos = 0.6)
    Phase2 => Phase3:(ct = 2.0, pos = 0.4)
    Phase3 => Filed:(ct = 3.0, pos = 0.65)
    Filed => Market:(ct = 1.0, pos = 0.9)
end

portfolio = build_portfolio()
populate_arcs!(portfolio)
println("@pipeline expanded the phase chain into a flat ReactionNetwork:")
println("  place (phases) : ", portfolio[:, :placeName])
println(
    "  parts            : ", nrows(portfolio, :S), " place, ",
    nrows(portfolio, :T), " transitions"
)
println("  per-edge (ct, pos):")
for i in row_ids(portfolio, :T)
    @printf(
        "    %-24s ct=%.1f  pos=%.2f\n",
        portfolio[i, :transName], portfolio[i, :transCycleTime],
        portfolio[i, :transProbOfSuccess]
    )
end

# ## 4. Zooming in: `refine` one transition, plug-compatibly
#
# The headline workflow: keep the portfolio coarse, but substitute a **detailed Phase-2
# sub-model** (screening → lead-opt → tox → filing, each its own timed, probabilistic step)
# for the single coarse `flow_Phase2_Phase3` transition. The sub-model declares its boundary
# as ports: `p2_in` plugs onto the parent's `Phase2`, `p2_out` onto `Phase3`; everything else
# is `:private` and gets namespaced.

phase2_detail = @reaction_network begin
    1.0, p2_in  --> screen, name => screening, cycletime => 0.5, probability => 0.85
    1.0, screen --> leadopt, name => lead_opt, cycletime => 0.7, probability => 0.8
    1.0, leadopt --> tox, name => tox_study, cycletime => 0.5, probability => 0.85
    1.0, tox    --> p2_out, name => filing_prep, cycletime => 0.3, probability => 0.7
end
set_port_role!(phase2_detail, :p2_in => :input, :p2_out => :output)

# Record the coarse boundary indices *before* the splice, so we can show they survive it.
phase2_ix_before = find_index(:Phase2, portfolio)
phase3_ix_before = find_index(:Phase3, portfolio)

# `refine(spec, transition, submodel; ports)` splices the sub-model into the named coarse
# transition: it namespaces the sub's private places, identifies the sub's ports with the
# parent boundary places by FK-repoint, appends the sub's transitions, and drops the coarse
# transition. It is **non-mutating** (`refine` = `refine!` on a `deepcopy`); `portfolio` is
# left intact.

refined = refine(
    portfolio, :flow_Phase2_Phase3, phase2_detail;
    ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out)
)

tnames_before = [portfolio[i, :transName] for i in row_ids(portfolio, :T)]
tnames_after = [refined[i, :transName] for i in row_ids(refined, :T)]
println("transitions BEFORE refine : ", tnames_before)
println("transitions AFTER  refine : ", tnames_after)
println("  coarse `flow_Phase2_Phase3` removed? ", !(:flow_Phase2_Phase3 in tnames_after))

# **Plug-compatibility (the point of the whole rung).** Because the boundary places keep
# their indices *and* their names, every transition *other* than the one we refined is
# structurally byte-for-byte identical before and after. The rest of the portfolio does not
# notice that Phase-2 became four sub-steps.

println()
println("PLUG-COMPATIBILITY:")
println(
    "  boundary places keep their indices — Phase2: ", phase2_ix_before, " → ",
    find_index(:Phase2, refined), "   Phase3: ", phase3_ix_before, " → ", find_index(:Phase3, refined)
)
untouched = [:flow_Discovery_Phase1, :flow_Phase1_Phase2, :flow_Phase3_Filed, :flow_Filed_Market]
identical = [tn => (trans_signature(portfolio, tn) == trans_signature(refined, tn)) for tn in untouched]
for (tn, same) in identical
    @printf("    %-24s identical before/after? %s\n", tn, same)
end
println(
    "  ⇒ all untouched transitions identical: ", all(last, identical),
    "  — the zoom is local to Phase 2."
)
sub_transitions = [n for n in tnames_after if occursin("__sub__", string(n))]
println("  sub-steps spliced in (namespaced): ", sub_transitions)

# ## 5. It's just a `ReactionNetwork`: construct and simulate the refined pipeline
#
# Refinement leaves no runtime trace — the refined model is an ordinary `ReactionNetwork`, so
# it constructs and simulates exactly like a hand-authored flat model. We attach a marking
# and a horizon and run both the coarse and the refined portfolio, reproducibly from
# `(model, seed)`. (Macro arguments are literal, so `@prob_init` takes a literal count; the
# sub-model's private pools default to zero.)

const SEED = 20260718
function attach_and_build(net; seed)
    @prob_init net Discovery = 60
    @prob_params net
    @prob_meta net tspan = 30 dt = 1.0
    return ReactionNetworkProblem(net; seed = seed)
end

coarse_prob = attach_and_build(build_portfolio(); seed = SEED)
refined_prob = attach_and_build(
    refine(
        build_portfolio(), :flow_Phase2_Phase3, phase2_detail;
        ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out)
    ); seed = SEED
)
simulate(refined_prob)

launched = refined_prob.sol[end, "Market"]
horizon = refined_prob.sol[end, "t"]
@printf(
    "Refined pipeline (seed=%d): 60 Discovery projects → %.0f reached :Market by t=%.0f\n",
    SEED, launched, horizon
)
subcols = [c for c in names(refined_prob.sol) if occursin("__sub__", c)]
println("  the detailed Phase-2 sub-net's pools appear in prob.sol: ", subcols)
println("  ⇒ the zoom is live: the refined model runs like any flat network.")

# A structural view of the splice — the coarse Phase-2 edge versus the four-step sub-net —
# renders through Graphviz. Rendering needs a backend, so we wrap it: the DOT source is
# always obtainable from `to_graphviz` even when no backend is present, so a hiccup can never
# fail the build.

structure = try
    coarse_svg = draw_network(coarse_prob; format = "svg")
    refined_svg = draw_network(refined_prob; format = "svg")
    HTML(
        "<div style=\"display:flex;flex-wrap:wrap;gap:1rem;justify-content:center\">" *
            "<figure style=\"margin:0\"><figcaption>coarse portfolio</figcaption>" * coarse_svg * "</figure>" *
            "<figure style=\"margin:0\"><figcaption>Phase 2 refined</figcaption>" * refined_svg * "</figure></div>"
    )
catch err
    @warn "draw_network: no Graphviz backend — showing the DOT source instead" exception = err
    Text(to_graphviz(network_graph(coarse_prob)))
end

# ### `abstract` — climbing back up the ladder
#
# `abstract_transitions` is the inverse of `refine`: it collapses a connected set of
# sub-transitions back into one coarse transition whose boundary reaction line is `lhs -->
# rhs`, carrying summarized attributes. It is a structural convenience for moving *up* the
# granularity ladder (it drops the sub rows and adds the coarse one; it does not garbage-collect
# the now-inert internal place).

collapsed = abstract_transitions(
    refined, sub_transitions, :flow_Phase2_Phase3;
    lhs = [:Phase2], rhs = [:Phase3],
    attrs = Dict(:transCycleTime => 2.0, :transProbOfSuccess => 0.4)
)
tnames_collapsed = [collapsed[i, :transName] for i in row_ids(collapsed, :T)]
println(
    "abstract_transitions collapsed the ", length(sub_transitions),
    " sub-steps back into one coarse `flow_Phase2_Phase3`:"
)
println(
    "  coarse transition restored? ", :flow_Phase2_Phase3 in tnames_collapsed,
    " ; T count back to ", nrows(collapsed, :T), " (coarse was ", nrows(portfolio, :T), ")"
)

# ## 6. The granularity-substitution guarantee
#
# Refinement does **not** claim the fine model is behaviorally *equivalent* to the coarse one
# — that would need a bisimulation the framework cannot check, and zooming in legitimately
# changes the dynamics. What it *does* offer is **plug-compatibility** (§4, structural) plus an
# **advisory boundary check** that makes the ladder auditable: on a linear chain, a coarse
# transition and its refined sub-model *agree in aggregate* when the coarse cycletime matches
# the **sum** of the sub cycletimes and the coarse probability matches the **product** of the
# sub probabilities. `refinement_diagnostics(submodel, coarse_attrs)` returns exactly those
# warnings (empty ⇒ clean); the author may override them.
#
# The coarse `flow_Phase2_Phase3` claimed `(ct = 2.0, pos = 0.40)`. We check the detailed
# sub-model we actually spliced in against those coarse attributes — and, for contrast, a
# deliberately *drifted* sub-model whose aggregates do not match.

coarse_attrs = Dict(:transCycleTime => 2.0, :transProbOfSuccess => 0.4)

Σct = 0.5 + 0.7 + 0.5 + 0.3           # sum of the four sub cycletimes
Πpos = 0.85 * 0.8 * 0.85 * 0.7        # product of the four sub PoS
warns_ok = refinement_diagnostics(phase2_detail, coarse_attrs)

drifted = @reaction_network begin
    1.0, p2_in --> mid, name => slow_a, cycletime => 3.0, probability => 0.9
    1.0, mid   --> p2_out, name => slow_b, cycletime => 3.0, probability => 0.9
end
set_port_role!(drifted, :p2_in => :input, :p2_out => :output)
warns_drift = refinement_diagnostics(drifted, coarse_attrs)

println("Plug-compatible sub-model (the one we spliced in):")
@printf("  Σ sub cycletimes = %.2f  (coarse ct = 2.00)\n", Σct)
@printf("  Π sub PoS        = %.3f  (coarse pos = 0.400)\n", Πpos)
println("  refinement_diagnostics : ", isempty(warns_ok) ? "clean — no warnings" : warns_ok)
println()
println("Drifted sub-model (Σct = 6.0, ΠPoS = 0.81):")
for w in warns_drift
    println("  ⚠ ", w)
end

# The aggregate agreement is easiest to *see*: the four sub-step cycletimes stack up to the
# coarse cycletime, and the four sub-step success probabilities multiply back to the coarse
# probability — for the plug-compatible refinement, not the drifted one.

p_ct = bar(
    ["coarse\n(1 step)", "refined\n(Σ 4 steps)", "drifted\n(Σ 2 steps)"],
    [2.0, Σct, 6.0]; legend = false, title = "cycletime", ylabel = "ticks",
    color = [:steelblue :seagreen :firebrick],
)
p_pos = bar(
    ["coarse\n(1 step)", "refined\n(Π 4 steps)", "drifted\n(Π 2 steps)"],
    [0.4, Πpos, 0.81]; legend = false, title = "prob. of success", ylabel = "probability",
    color = [:steelblue :seagreen :firebrick],
)
plot(p_ct, p_pos; layout = (1, 2), size = (760, 320), plot_title = "Granularity substitution: refined ≈ coarse, drifted ⇏ coarse")

# ### Reading the result
#
# **The plug-compatible Phase-2 sub-model agrees with the coarse transition on both boundary
# aggregates** — Σ cycletimes = 2.0 against the coarse 2.0, and Π success = 0.40 against the
# coarse 0.40 — so `refinement_diagnostics` returns *clean*. The drifted sub-model, whose
# steps sum to 6.0 ticks and multiply to 0.81, trips both checks. That is the
# granularity-substitution guarantee in operational form: you may model Phase 2 at whatever
# resolution the question demands, and a cheap, deterministic check tells you when the finer
# model still *stands in* for the coarser one at the boundary — without ever claiming a
# behavioral equivalence the framework cannot honestly prove.
#
# What matters is the *kind* of guarantee, not the specific numbers: a coarse portfolio for
# fast what-if, a finer sub-model spliced in plug-compatibly where fidelity is needed, and an
# advisory boundary check that keeps the two rungs honest — all as authoring-time structural
# rewrites over a static network, never touched once a run is live.

# ## Recap
#
# The granularity ladder, bottom to top:
#
# 1. `@join` / `@equalize` — manual, no-declared-ports union and collapse;
# 2. `@process` / `@port` / `@compose` — reusable fragments composed by *declared open ports*,
#    matched automatically by FK-repoint;
# 3. `@pipeline` — a whole phase chain authored in one block as `flow` routing transitions;
# 4. `refine` — substitute a finer sub-model for one coarse transition, **plug-compatibly**
#    (boundary places keep their indices/names; every other transition is untouched);
# 5. `abstract` — the inverse, collapsing sub-steps back into one coarse transition;
# 6. `refinement_diagnostics` — the advisory Σ-ct / Π-PoS boundary check that makes the ladder
#    auditable, the granularity-substitution guarantee in computable form.
#
# Every rung is an authoring-time rewrite of a static `ReactionNetwork` — forbidden on a
# live, stepping model, because it reindexes the store. The [serialization deep-dive](serialization.md)
# shows the flip side: a composed or refined model is still just data, and round-trips through
# JSON with no runtime trace.
