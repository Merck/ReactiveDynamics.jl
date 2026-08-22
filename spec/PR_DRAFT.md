# DRAFT PR — ReactiveDynamics.jl: native discrete-event engine rework

> Draft while the two companion papers and the refined case-study HTML are finished; the engine rework, the semantic suite, and the documentation site are complete and land together in this PR.

## What this PR does

This rework pivots ReactiveDynamics.jl from a SciML and Catlab base to a **native, dependency-light discrete-event engine** for timed, stochastic, resource-constrained modeling of business and R&D processes — budgeting, ledgers, what-if analysis, rNPV — with first-class support for agentic, structured resources and for composition. The model is now accommodated within AlgebraicAgents.jl, which provides the modeling interface and unlocks these goals: agentic resources and composability with third-party models.

The motivation: years of unresolved technical debt were limiting the framework's applicability and publication prospects, and burying its conceptual novelty under technical-level constraints. The rework revisits how the core pillars are implemented — the DSL, the data store, and the simulation engine, focusing on correctness and intended semantics — introduces proper semantic tests, and delivers a documentation site.

Everything is specified contract-first: a normative operational-semantics contract (`spec/CONTRACT_DRAFT.md` §1–§15) and 15 Architecture Decision Records (`spec/adr/`), with the engine built to satisfy them. Implementation state is tracked in `spec/STATUS.md`.

The classical modeling surface stays familiar. A plain-place SIR, end to end:

```julia
using ReactiveDynamics

sir = @reaction_network begin
    α * S * I, S + I --> 2I, name => I2R    # a bare numeric rate is a stochastic (Poisson) intensity
    β * I,     I     --> R,  name => R2S
end
@prob_init  sir S = 999 I = 10 R = 0
@prob_params sir α = 0.0001 β = 0.01
@prob_meta  sir tspan = 250 dt = 0.1

prob = ReactionNetworkProblem(sir; seed = 1)   # seed= owns the per-run RNG — the only route to reproducibility
simulate(prob)
prob.sol[!, "I"]                                # read solution columns BY NAME (order is construction order)
```

## Headline changes

**Modeling capabilities**

This is the consequential part of the rework: the framework now natively expresses what previously required host-side workarounds — structured tokens with identity, value-qualified selection, in-model decision rules, resource contention semantics, and hierarchical composition.

- Structured, agentic tokens with live instantiation and query; a host-function registry replaces `@register`-time `eval` (ADR 0006). A token is a first-class entity with attributes and a stable identity, not an anonymous count:

  ```julia
  @register begin
      @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
          phase::Symbol      # the canonical "phase-as-attribute" (ADR 0008 §D)
          npv::Float64
      end
  end
  ```
- Token filtration: `TokenPredicate` and `@select` pick tokens by an 𝓕ₜ-measurable predicate; phase-as-attribute is canonical; `@advance` and `SetField` write fields while preserving identity (ADR 0008). A pipeline step selects a value-qualified subset and advances the same object in place:

  ```julia
  @deterministic(1.0),
      @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(phase, :Phase3),
      name => fasttrack, cycletime => 1.0, probability => 1.0
  ```
- Resource modality as a per-participation tag on the LHS — `@conserved` (returned at finish), `@rate` (drawn per in-flight tick), `@nonblock` (claimed, not held) — so contention is properly modeled, not hand-coded:

  ```julia
  pipeline = @reaction_network begin
      @deterministic(2.0),
          @select(Project, phase == :Phase2) + 4 * @conserved(scientist) + 5 * @rate(budget) -->
          @advance(phase, :Phase3),
          name => adv_phase2, cycletime => 2.0, probability => 0.4, priority => 2.0
  end
  ```
- Genesis as a first-class transition product: `@structured(:Kind, field = …)` mints a fresh token on the RHS (the agentic `∅ --> place`), registry-resolved so it too serializes eval-free (ADR 0006, ADR 0005). Field expressions read live state (`@t()`) and the seeded RNG:

  ```julia
  @deterministic(1.0),
      ∅ --> @structured(:Project, phase = :Phase1, npv = rand(state.rng, Normal(120.0, 20.0)), born = @t()),
      name => genesis
  ```
- Rules, triggers, and conditional transitions form the endogenous decision channel (ADR 0010), with action callbacks `SetTokens` and `Invoke` (ADR 0011). A management lever lives in the model as a typed `Rule`, not as host patch code:

  ```julia
  Rule(:series_b, :(@t() > 2.0),
       Seq([SetMarking(:cash, 500, :inc),          # inject capital
            SetParams([:synergy => 1]),            # flip a parameter the rates read
            AddToken(:Project, [:phase => QuoteNode(:Phase2), :npv => 175.0])]);  # add a program
       fire_mode = :once)
  ```
- Hierarchical refinement and open-port composition: `@pipeline`, `@process`, and `@compose` author coarsely; `refine` substitutes a finer sub-process for one step via FK splice, non-mutating and leaving boundary places in place (ADR 0009):

  ```julia
  portfolio = @pipeline Project begin
      Discovery => Phase1:(ct = 1.0, pos = 0.45)
      Phase1    => Phase2:(ct = 1.5, pos = 0.6)
      Phase2    => Phase3:(ct = 2.0, pos = 0.4)
  end
  refined = refine(portfolio, :flow_Phase2_Phase3, phase2_detail;
                   ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out))
  ```
- Interface lifecycle (authoring → construction → live), a declarative, serializable `population[]` initial marking, and `dump_state` with `restore` (ADR 0007). The starting portfolio is reproducible input, not imperative post-construction code:

  ```julia
  ReactionNetworkProblem(pipeline_model(); seed = 1, registry = REGISTRY,
      population = [ProjectToken(:Phase1, 120.0), ProjectToken(:Phase2, 200.0), ProjectToken(:Phase3, 300.0)])
  ```
- AlgebraicAgents integration: ReactiveDynamics as an AA hierarchy node — outbound via `getobservable` and parameters, inbound via wired external coupling (`inputs[]`, `ExternalRef`, `_prestep!` with a one-tick Jacobi lag) (ADR 0012). The foreign-agent topology lives host-side in `add_wire!`, never in the model document:

  ```julia
  add_wire!(root; from = market, to = rd,      from_var_name = "sentiment", to_var_name = "sentiment")
  add_wire!(root; from = rd,     to = finance,  from_var_name = "cash",      to_var_name = "rd_cash")
  ```
  AA moved to the published registry 0.4 release. `@agentize` provides thin auto-naming sugar over the `ReactionNetworkProblem` constructor (ADR 0001, ADR 0012) with no second construction path.

**Engine and semantics**

- Native discrete-event engine (`ReactionNetworkProblem` stepped via AA's `_step!`); the SciML stack (`DifferentialEquations`, `OrdinaryDiffEq`, `DiffEqBase`) and the old `DiscreteProblem` transform were removed entirely — no dependencies, no code (ADR 0001). A SciML interop adapter is left as a possible future package extension, but none ships here.
- Priority-weighted progressive-fill (water-filling) resource allocator — work-conserving, deterministic, dependency-free (ADR 0002).
- Append-only mutation with soft deactivation, so transitions, places, and parameters can be added — and transitions retired — mid-simulation without breaking position-indexed compiled closures (ADR 0004).
- `AbstractRNG` and `seed=` threaded through every draw; a run is fully determined by `(model, seed)` (CONTRACT §4).
- Construction-time modality validation (`validate_modalities`, CONTRACT §1.4): rejects the three illegal modality configurations (`{:nonblock,:conserved}`; `:rate` with concrete `cycletime==0`; `:rate` on a structured place) with a clear `ArgumentError` before any tick, replacing late or silent failures.

**Data store and serialization**

- ACSets and Catlab **dropped** in favor of a dependency-free, typed struct-of-columns IR; the transition–place relation (the net's arcs) promoted to a first-class typed `ArcSpec` incidence table (ADR 0003).
- A single eval-free JSON serialization with a typed `ExprNode` IR: `from_json_model` and `to_json_model` round-trip, plus `validate`. This closes the import-time RCE and retires the TOML, CSV, and JLD2 format zoo (ADR 0005).
- Post-ACSets naming pass (ADR 0015): `@reaction_network` (was `@ReactionNetworkSchema`), `net` (was `acs`), store type `ReactionNetwork`; store verbs renamed to store vocabulary (`nrows`, `row_ids`, `column`, `cell`, `find_rows`, …) **and unexported**. Old names survive one release as `@deprecate` shims; `GeneratedExpressions` dropped.

A model IS data: the same pipeline as an eval-free JSON document, which `validate` checks and `from_json_model` loads to a run bit-identical to the DSL build (host token kinds referenced by name through a registry — the document carries no Julia):

```json
{ "rd_format":"reactive-dynamics-model", "version":"1.0",
  "meta":{ "tspan":6.0, "dt":1.0 },
  "places":[ {"name":"Project","structured":true} ],
  "transitions":[
    {"id":"adv12","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":1.0} ],
  "arcs":[
    {"transition":"adv12","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase1"]]}},
    {"transition":"adv12","side":"rhs","advance":{"field":"phase","value":"Phase2"}} ] }
```

```julia
diags = validate(JSON.parse(json); registry = REGISTRY)   # [] ⇒ clean; a dangling FK is a Diagnostic, never an eval
prob  = from_json_model(json; seed = 7, registry = REGISTRY, population = pop)
to_json_model(prob)                                        # the inverse — a live model back to a document, loss-free
```

**Analysis and visualization**

- Per-token trajectory log with `representative_token` and `trajectory_envelope`; an ensemble runner (`ensemble`, `summarize`, `treatment_effect`, with `EnsembleProblem` as an AA node); a results export bundle — CSV and JSON in core, Arrow as a weak dependency (ADR 0013). The ensemble runner supports both rebuild-per-seed (mode a) and reinit-reseed member reuse (`ensemble(...; mode=:reinit)`, mode b via `_reinit!(state; seed=…)`), with member-for-member equivalence of the two modes under test.
- Model-agnostic Plots recipes behind `RDPlotsExt`; the three-layer network "exec map" (`network_graph`, `draw_network`, `exec_map`) with bottleneck and starvation coloring and `@select` token highlighting (ADR 0014).
- `Plots` and `Arrow` demoted to weak dependencies with the `RDPlotsExt` and `RDArrowExt` package extensions; `Pluto`, `PlutoUI`, `IJulia`, and `DifferentialEquations` dropped from dependencies.

**Documentation and records**

- Contract-first: the normative operational-semantics specification (`spec/CONTRACT_DRAFT.md` §1–§15) and 15 ADRs (`spec/adr/`), all statuses truthed-up to Implemented with commit citations. The durable engineering artifacts live under top-level `spec/`, kept separate from `docs/` (the Documenter static-pages site).
- Repo-root `CLAUDE.md` agent guide; `spec/STATUS.md` as the single state index; completed handoff plans archived under `spec/history/`.

## Tests

`test/semantic/*.jl`, two-tier — T1 characterization and T2 acceptance — with real assertions over operational semantics, not smoke tests. The suite is fully green (801 pass / 0 broken, no skips). `test/Project.toml` declares the test-only dependencies (Plots, Arrow, DataFrames, Distributions) so the weak-dependency extension paths are exercised, not skipped. A GitHub Actions scaffold (tests matrix, docs, Runic, TagBot, CompatHelper) has been added; the suite and the formatter ([Runic](https://github.com/fredrikekre/Runic.jl)) are also run locally.

## Documentation and tutorials

Documentation and tutorials are now part of this rework (previously slated for a standalone follow-up; that work — branch `docs-tutorials` — has been consolidated into `rework`). The seven runnable `demo/` tours are migrated **in place** as the single source of truth: their `.jl` sources are the Literate inputs the site ingests, so there is no duplicated model code to drift. The site is a [Documenter.jl](https://documenter.juliadocs.org) + [Literate.jl](https://fredrikekre.github.io/Literate.jl) build (HTML, GitHub Pages); every code block executes against the current engine at build time, pinned by an explicit `seed=`. It is chartered and tracked in `spec/DOCS_CHARTER.md`.

The structure is Diátaxis, routed on the landing page by reader intent:

- **Tutorials** (learning) — three tiered, cumulative onboarding tutorials (**introductory** → **advanced** → **expert**) plus two technical **deep-dives** (serialization; composition & granularity). Each tutorial ends by *computing* a decision-relevant quantity — a marginal effect, a treatment effect with its standard error — not a feature recap.
- **Applied case studies** (understanding) — three question-titled, headline-number-first decision case studies: *"What is the marginal eNPV of the Nth scientist?"* (the shadow price of the binding resource; flagship), *"What is this in-licensing asset worth to this pipeline?"* (BD/M&A — rNPV is not additive under contention), and *"When should you kill a program?"* (the value-maximizing kill threshold). The flagship and BD studies get a refined, self-contained HTML presentation.
- **API reference** (information) — capability-organized autodoc pages off the real (post-ADR-0015) export surface: authoring, structured tokens, rules & actions, construction & simulation, composition, serialization, analysis & visualization, and AlgebraicAgents coupling — plus a JSON model-schema page.
- **The "why"** — rather than re-hosting the semantics as web pages, the site links to the normative contract (`spec/CONTRACT_DRAFT.md` §1–§15) and the ADRs at their `spec/` home, and to two companion papers: a **methodology paper** (the operational semantics and backing concepts, arXiv-first) and an **HBR-style adoption/value paper** for non-technical executives (motivation → solution → examples).

The tiered tutorials, deep-dives, case studies, and the full API reference are delivered and render green. Still tracked as follow-on (see `spec/DOCS_CHARTER.md` §10): the two companion papers (D3a/D3b) and the refined-HTML presentations for the flagship and BD case studies (E3).

## Genuinely deferred (gates recorded; not part of this PR)

- Entity-level refinement — a structured token hosting its own sub-network (ADR 0009 §F, a future ADR).
- Threaded ensemble backend (`ensemble(...; parallel=true)` is accepted but runs sequentially).
- AA `Opera`-level implicit fixed-point coupling (current coupling is an explicit one-tick-lag Jacobi scheme).
- `dump_state` in-flight limitation — implemented, but it refuses to dump mid-cycle transitions (a scoped Milestone-1 constraint, not a gap).
