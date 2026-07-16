# DRAFT PR — ReactiveDynamics.jl: native discrete-event engine rework

This intentionally stays draft until the docs/tutorial refinement lands.

## What this is

I decided to address the years of accummulated/unresolved technical debt that would limit the applicability and publication prospects of the framework, with the intent to revisit/solidify the core foundations of the framework (DSL, data store, simulation engine), and proper semantic tests, and improve documentation. We pivot from a SciML/Catlab embedding to a **native, dependency-light discrete-event engine** for timed, stochastic, resource-constrained business/R&D process modeling (enabling budgeting, ledgers, what-if, rNPV). The whole rework is specified contract-first: a normative operational-semantics contract (`docs/CONTRACT_DRAFT.md` §1–§15) and 15 Architecture Decision Records (`docs/adr/`), with the engine built to satisfy them. State of the implementation is tracked in `docs/STATUS.md`.

The whole modeling surface is one small script. A classical (plain-species) SIR, end to end:

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

**Engine & semantics**
- Native discrete-event engine (`ReactionNetworkProblem` stepped via AA's `_step!`); SciML demoted to optional (ADR 0001).
- Priority-weighted progressive-fill (water-filling) resource allocator — work-conserving, deterministic, dependency-free (ADR 0002).
- Append-only + soft-deactivate runtime mutation, so transitions/species/params can be added and transitions retired mid-simulation without breaking position-indexed compiled closures (ADR 0004).
- `AbstractRNG`/`seed=` threaded through every draw; a run is fully determined by `(model, seed)` (CONTRACT §4).
- Construction-time modality validation (`validate_modalities`, CONTRACT §1.4): rejects the three illegal modality configs (`{:nonblock,:conserved}`; `:rate` with concrete `cycletime==0`; `:rate` on a structured species) with a clear `ArgumentError` before any tick, replacing late/silent failures.

Resource modality is a per-participation tag on the LHS — `@conserved` (returned at finish), `@rate` (drawn per in-flight tick), `@nonblock` (claimed, not held) — so contention is modeled, not hand-coded:

```julia
pipeline = @reaction_network begin
    @deterministic(2.0),
        @select(Project, phase == :Phase2) + 4 * @conserved(scientist) + 5 * @rate(budget) -->
        @advance(phase, :Phase3),
        name => adv_phase2, cycletime => 2.0, probability => 0.4, priority => 2.0
end
```

**Data store & serialization**
- ACSets/Catlab **dropped** for a dependency-free typed struct-of-columns IR; the transition↔reactant relation promoted to a first-class typed `ReactantSpec` incidence table (ADR 0003).
- Single eval-free JSON serialization + typed `ExprNode` IR with `from_json_model`/`to_json_model` round-trip + `validate`; closes the import-time RCE; drops the TOML/CSV/JLD2 zoo (ADR 0005).
- Post-ACSets naming rename (ADR 0015): `@reaction_network` (was `@ReactionNetworkSchema`), `net` (was `acs`), store type `ReactionNetwork`, store verbs renamed to store vocabulary (`nrows`/`row_ids`/`column`/`cell`/`find_rows`/…) **and unexported**; old names survive one release as `@deprecate` shims; `GeneratedExpressions` dropped.

A model IS data: the same pipeline as an eval-free JSON document, which `validate` checks and `from_json_model` loads to a run bit-identical to the DSL build (host token kinds referenced BY NAME through a registry — the document carries no Julia):

```json
{ "rd_format":"reactive-dynamics-model", "version":"1.0",
  "meta":{ "tspan":6.0, "dt":1.0 },
  "species":[ {"name":"Project","structured":true} ],
  "transitions":[
    {"id":"adv12","rate":1.0,"rate_mode":"deterministic","cycletime":1.0,"prob_of_success":1.0} ],
  "reactants":[
    {"transition":"adv12","side":"lhs","predicate":{"kind":"Project","clauses":[["phase","==","Phase1"]]}},
    {"transition":"adv12","side":"rhs","advance":{"field":"phase","value":"Phase2"}} ] }
```

```julia
diags = validate(JSON.parse(json); registry = REGISTRY)   # [] ⇒ clean; a dangling FK is a Diagnostic, never an eval
prob  = from_json_model(json; seed = 7, registry = REGISTRY, population = pop)
to_json_model(prob)                                        # the inverse — a live model back to a document, loss-free
```

**Modeling language**
- Structured/agentic tokens with live instantiation/query; host-function registry replaces `@register` (eval-free) (ADR 0006). A token is a first-class entity with attributes and a stable identity, not an anonymous count:

  ```julia
  @register begin
      @aagent BaseStructuredToken AbstractStructuredToken struct ProjectToken
          phase::Symbol      # the canonical "phase-as-attribute" (ADR 0008 §D)
          npv::Float64
      end
  end
  ```
- Interface lifecycle (authoring → construction → live) + declarative serializable `population[]` initial marking + `dump_state`/`restore` (ADR 0007). The starting portfolio is reproducible INPUT, not imperative post-construction code:

  ```julia
  ReactionNetworkProblem(pipeline_model(); seed = 1, registry = REGISTRY,
      population = [ProjectToken(:Phase1, 120.0), ProjectToken(:Phase2, 200.0), ProjectToken(:Phase3, 300.0)])
  ```
- Token filtration: `TokenPredicate`/`@select` selects tokens by 𝓕ₜ-measurable predicate; phase-as-attribute canonical; `@advance`/`SetField` field-writes preserving identity (ADR 0008). A pipeline step selects a value-qualified subset and advances the same object in place:

  ```julia
  @deterministic(1.0),
      @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(phase, :Phase3),
      name => fasttrack, cycletime => 1.0, probability => 1.0
  ```
- Genesis as a first-class transition product: `@structured(:Kind, field = …)` mints a fresh token on the RHS (the agentic `∅ --> species`), registry-resolved so it too serializes eval-free (ADR 0006/0005). Field exprs read live state (`@t()`) and the seeded RNG:

  ```julia
  @deterministic(1.0),
      ∅ --> @structured(:Project, phase = :Phase1, npv = rand(state.rng, Normal(120.0, 20.0)), born = @t()),
      name => genesis
  ```
- Hierarchical refinement & open-port composition: `@pipeline`/`@process`/`@compose` author coarsely; `refine` substitutes a finer sub-process for one step via FK-splice, non-mutating and leaving boundary species in place (ADR 0009):

  ```julia
  portfolio = @pipeline Project begin
      Discovery => Phase1:(ct = 1.0, pos = 0.45)
      Phase1    => Phase2:(ct = 1.5, pos = 0.6)
      Phase2    => Phase3:(ct = 2.0, pos = 0.4)
  end
  refined = refine(portfolio, :flow_Phase2_Phase3, phase2_detail;
                   ports = Dict(:Phase2 => :p2_in, :Phase3 => :p2_out))
  ```
- Rules/triggers & conditional transitions — the endogenous decision channel (ADR 0010); action callbacks `SetTokens`/`Invoke` (ADR 0011). A management lever lives IN the model as a typed `Rule`, not host patch code:

  ```julia
  Rule(:series_b, :(@t() > 2.0),
       Seq([SetSpecies(:cash, 500, :inc),          # inject capital
            SetParams([:synergy => 1]),            # flip a parameter the rates read
            AddToken(:Project, [:phase => QuoteNode(:Phase2), :npv => 175.0])]);  # add a program
       fire_mode = :once)
  ```
- AlgebraicAgents integration: RD as an AA hierarchy node (outbound `getobservable`/params) + inbound wired external coupling (`inputs[]`/`ExternalRef`/`_prestep!` one-tick Jacobi lag) (ADR 0012). The foreign-agent topology lives host-side in `add_wire!`, never in the model document:

  ```julia
  add_wire!(root; from = market, to = rd,      from_var_name = "sentiment", to_var_name = "sentiment")
  add_wire!(root; from = rd,     to = finance,  from_var_name = "cash",      to_var_name = "rd_cash")
  ```
  AA moved to the published registry 0.4 release. `@agentize` provides thin auto-naming sugar over the `ReactionNetworkProblem` constructor (ADR 0001/0012) with no second construction path.

**Analysis & visualization**
- Per-token trajectory log + `representative_token`/`trajectory_envelope`; ensemble runner (`ensemble`/`summarize`/`treatment_effect`, `EnsembleProblem` as an AA node); results export bundle (CSV/JSON core, Arrow weakdep) (ADR 0013). The ensemble runner supports both rebuild-per-seed (mode a) and reinit-reseed member reuse (`ensemble(...; mode=:reinit)`, mode b via `_reinit!(state; seed=…)`), with mode (a)≡mode (b) member-for-member equivalence tested.
- Model-agnostic Plots recipes behind `RDPlotsExt`; the three-layer network "exec map" (`network_graph`/`draw_network`/`exec_map`) with bottleneck/starvation coloring and `@select` token highlighting (ADR 0014).
- `Plots`/`Arrow` demoted to weakdeps with `RDPlotsExt`/`RDArrowExt` package extensions; `Pluto`/`PlutoUI`/`IJulia`/`DifferentialEquations` dropped from deps.

**Documentation & records**
- Contract-first: normative operational-semantics spec (`docs/CONTRACT_DRAFT.md` §1–§15) and 15 ADRs (`docs/adr/`), all statuses truthed-up to Implemented with commit citations.
- Repo-root `CLAUDE.md` agent guide; `docs/STATUS.md` as the single state index; the completed handoff plans archived under `docs/history/`.

## Tests

`test/semantic/*.jl`, two-tier (T1-characterization / T2-acceptance) — real assertions over operational semantics, not smoke tests. `test/Project.toml` declares the test-only deps (Plots, Arrow, DataFrames, Distributions) so the weakdep extension paths are exercised, not skipped. No CI is configured in the repo — the suite and the formatter ([Runic](https://github.com/fredrikekre/Runic.jl)) are run locally.

## Documentation and Tutorials

Docs/tutorial refinement is deliberately NOT in this PR. This PR should not be merged until the rework until the demos/tutorials are polished and coherent against the final (post-ADR-0015) surface. This will be addressed through a standalone PR into `rework`. There are currently seven demos illustrating the various facets of the framework along with applied examples. The rough plan is to improve this further and
- have onboarding tutorials for tiered user expertise level (introductory, advanced, expert); these should be in literate form and clearly didactive, on sample problems; consider the end-to-end modeling workflow in varying depth of detail;
- in addition, worked out examples demonstrating the modeling value on examples such as "business development"/MA valuation, optimization/decision making in R&D.
- have usual docs for the API, also interfacing the tutorials.

These should come in literate form and be interfaced from the official github pages documentation, with the applied examples deserving more refined presentations (HTML, and beyond).

## Genuinely deferred (gates recorded, NOT part of this or the follow-up)

- Entity-level refinement — a structured token hosting its own sub-network (ADR 0009 §F, a future ADR).
- Threaded ensemble backend (`ensemble(...; parallel=true)` is accepted but runs sequentially).
- AA `Opera`-level implicit/fixed-point coupling (current coupling is explicit one-tick-lag Jacobi).
- `dump_state` in-flight limitation — implemented but refuses to dump mid-cycle transitions (a scoped Milestone-1 constraint, not a gap).
