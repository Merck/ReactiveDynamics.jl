# Documentation & Tutorials Charter — ReactiveDynamics.jl

> **The tracking spine for the documentation rework.** A standalone PR into `rework` (branch `docs-tutorials`). This is the durable "what are we building, in what order, and how do we know a piece is done" index for the docs effort, the counterpart to [STATUS.md](STATUS.md) for the engine. Last updated 2026-07-17. Re-verify any `file:line` before acting — line numbers drift.

## 1. Why this exists and what "done" means

The engine rework (§1–§15 of [CONTRACT_DRAFT.md](CONTRACT_DRAFT.md), ADRs 0001–0015) is implemented and green (801 pass / 0 broken). What is missing is the reader-facing surface: onboarding tutorials, applied case studies, and an API reference that match the post-ADR-0015 API. The current [docs/src/index.md](../docs/src/index.md) is a single stale API page referencing removed macros (`@ReactionNetwork`, `@optimize`, `@fit`, `@problematize`, `@plot`, `@import_network`, …), and [docs/make.jl](../docs/make.jl) still uses the deprecated `DocumenterMarkdown` backend (not even declared in `docs/Project.toml`). The seven runnable tours under [demo/](../demo) are the current source of truth for working code but are not wired into a published site.

This PR closes that gap. **Done** for the whole PR means: a Documenter site, publishable to GitHub Pages, comprising three tiered tutorials, two focused deep-dives, three applied case studies, a complete API reference, and an explanation layer built on the contract — every code block executed against the current engine, every tutorial ending in a manager-actionable number, and the `demo/` sources migrated in as the single source of truth (no duplicated model code). The per-facet acceptance criteria are in §9; the progress tracker is in §10.

## 2. Design principles (non-negotiable)

These are the brand and quality commitments the whole effort is measured against. They come from the maintainer framing in the PR draft ([PR_DRAFT.md](PR_DRAFT.md) §"Documentation and tutorials") and are binding on every page.

- **Every applied example is a decision case study with a headline number, not a feature tour.** Title each case study as the *question it answers*; put the number and the verdict in the first paragraph; let the framework mechanics appear as "how we got it." This converts docs effort into artifacts that share the skeleton of a decision memo (question → number → implication) so a case study can *seed* a memo rather than compete with it. A feature-organized gallery ("modalities demo", "composition demo") signals tool-builder; a question-organized one signals operator. Same content, different brand transfer.
- **Diátaxis, with the contract promoted to the explanation layer.** Tutorials are learning-oriented; case studies are understanding-oriented; API docs are reference; and `CONTRACT_DRAFT.md` becomes the explanation layer — the differentiator versus every ad-hoc discrete-event-simulation package and the strongest argument in a JOSS / software-paper review. See §3.
- **The decision-number discipline applies to *every* tutorial, even the introductory one.** Each tutorial ends by extracting a number a manager would act on — not a recap of features. This is what keeps the learning-oriented tier honest to the operator brand.
- **Single source of truth: migrate the demos in place.** The `demo/*.jl` literate sources become the Literate.jl inputs the site ingests. Model-building code is written once; there is no parallel copy in `docs/` to drift. See §4 for the demo→docs mapping.
- **No invented API.** Every construct on every page must trace to a passing semantic test (`test/semantic/*.jl`) or an engine docstring, exactly as the demos already assert. When docs and code disagree, the code is truth (per [CLAUDE.md](../CLAUDE.md)).
- **Executable and reproducible.** Every code block runs during the site build (Literate `@example`/executed markdown); every run is pinned by an explicit `seed=`, so the rendered numbers are reproducible from `(model, seed)` (CONTRACT §4).

## 3. Diátaxis structure

| Quadrant | Orientation | What lives here | Source |
|---|---|---|---|
| **Tutorials** | learning | Three tiered onboarding tutorials (introductory / advanced / expert) + two focused technical deep-dives (serialization, composition). End-to-end modeling workflows at increasing depth. | Workstream A (§5) |
| **Case studies** (how-to / understanding) | understanding | Three applied decision case studies, question-titled, each with a headline number. The flagship gets a refined HTML presentation. | Workstream B (§6) |
| **Reference** | information | The complete public API, autodoc-generated from docstrings, organized by capability; the attribute/shorthand tables; the JSON model schema. | Workstream C (§7) |
| **Explanation** | understanding | The normative operational-semantics contract (§1–§15) and the ADRs, framed for readers as *why the engine behaves this way* — the differentiator. | Workstream D (§8) |

The site's landing page routes a reader by intent: "new here" → introductory tutorial; "what can it do for my decision" → case studies; "how do I call X" → reference; "why does it behave this way" → explanation/contract.

## 4. Toolchain, mechanics, and the migrate-in-place rule

**Build stack.** [Documenter.jl](https://documenter.juliadocs.org) for the site (HTML backend, GitHub Pages deploy) + [Literate.jl](https://fredrikekre.github.io/Literate.jl) for the executable tutorial/case-study sources. This replaces the deprecated `DocumenterMarkdown` backend currently in `docs/make.jl`. Literate sources under `docs/literate/**` are processed to `docs/src/**` markdown with `@example` blocks that Documenter executes at build time, so code runs end-to-end and its output is captured into the page. The one demo already on Literate ([demo/wires_viz_tour](../demo/wires_viz_tour), with its `build.jl` render recipe) is the working precedent for the executed-Literate → HTML path and the SVG-inlining trick for self-contained pages.

**Target docs tree** (the structure the charter builds toward; `[x]` = exists after this PR's exemplar, `[ ]` = tracked):

```
docs/
  make.jl                        # [ ] Documenter + Literate site build (replaces DocumenterMarkdown)
  Project.toml                   # [ ] + Documenter, Literate (and Plots/Arrow for case-study renders)
  literate/                      # single source of truth — migrated from demo/
    tutorials/
      introductory.jl            # [x] EXEMPLAR (this PR) — migrated + reshaped from core_engine_tour
      advanced.jl                # [ ] from agentic_pipeline (+ modalities/ensemble from core_engine_tour)
      expert.jl                  # [ ] from agentic_pipeline §6-7 + aa_integration + refinement_tour
    deep_dives/
      serialization.jl           # [ ] model-is-data JSON round-trip, registry, validate, RCE boundary
      composition.jl             # [ ] @join/@equalize/@compose/@pipeline/refine granularity ladder
    case_studies/
      marginal_scientist.jl      # [ ] FLAGSHIP — shadow price of the binding resource
      inlicensing_value.jl       # [ ] BD/M&A — rNPV is not additive under contention
      kill_a_program.jl          # [ ] kill-threshold Rule variants via treatment_effect
  src/
    index.md                     # [ ] landing / intent router (rewrite of the stale API page)
    tutorials/                   # generated from literate/tutorials
    case_studies/                # generated from literate/case_studies
    reference/                   # autodocs by capability (Workstream C)
    explanation/                 # contract + ADRs framed for readers (Workstream D)
    assets/                      # existing diagram*.png etc. (audit for staleness)
```

**Migrate-in-place rule (the anti-drift contract).** For each of the seven demos, the demo's `.jl` becomes (or is refactored into) the Literate source the site ingests; the demo directory keeps its README and its `--project=demo/x` runnability by pointing at the same source, so a demo stays independently runnable AND is the site's input. Where a demo has a demo-local `Project.toml` for weakdeps (Plots/Arrow), the site build reuses that dependency set. No page re-authors model code that a demo already contains; a tutorial that needs a subset of a demo extracts it into a shared included file rather than copying. The demo→tutorial mapping:

| Demo | Format today | Migrates into | Notes |
|---|---|---|---|
| `core_engine_tour` | script, `--project=.`, Statistics only | **introductory** tutorial (§1,§2,§8 subset) + **advanced** (modalities §3, allocator §4, genesis §5, ledger §6) + **composition deep-dive** (§7) | Plain-species; no weakdeps. The intro tier is a strict, reshaped subset. |
| `agentic_pipeline` | script, `--project=.`, self-contained token kind | **advanced** tutorial (structured tokens, `@select`/`@advance`, rules) + **serialization deep-dive** (§6 JSON) + **expert** (checkpoint §7) | The light companion to `bd_acquisition`. |
| `introspection_tour` | Literate-ready script, demo-local env (Plots/Arrow) | **advanced/expert** tutorial (trajectory log, ensemble, export, exec map) + feeds case-study visuals | The exec map is a hero visual for the flagship case study. |
| `refinement_tour` | script, thin demo-local env | **composition deep-dive** (granularity ladder) + **expert** | Structural; `@pipeline`/`@compose`/`refine`/`abstract`. |
| `aa_integration` | script, `--project=.` | **expert** tutorial (AA hierarchy node, wires, `_prestep!` Jacobi lag) | Prose/`println` twin of `wires_viz_tour`. |
| `wires_viz_tour` | **Literate + HTML** (`build.jl`) | **expert** tutorial (drawn) + toolchain precedent | Already the target format; its `build.jl` is the render recipe to generalize. |
| `bd_acquisition` | scripts + JSON model + HTML presentation | **flagship + in-licensing case studies** (refined HTML) | Its `MVP_BD_DEMO.md` design doc defines the decision-case-study DNA; its `presentation.html` is the refined-HTML precedent. |

**Refined HTML for case studies.** The flagship case study (and the in-licensing one) get a refined, self-contained HTML presentation beyond the Documenter page, following the existing `demo/bd_acquisition/presentation.html` precedent (client-side SVG charts from a real ensemble dump, HBR-style narrative). The [HTML→PDF recipe](../demo/bd_acquisition) (headless Chrome, print stylesheet) remains available for a PDF artifact.

## 5. Workstream A — Tiered onboarding tutorials

Learning-oriented, literate, didactic, on sample problems, framed as end-to-end modeling workflows at increasing depth. **Each ends in a manager-actionable number** (the discipline in §2). Tier boundaries follow the maintainer's mapping in the PR draft.

### A1 — Introductory tutorial *(EXEMPLAR delivered this PR)*

- **Sample problem.** A plain-species timed pipeline. Opens on the classic **SIR** epidemic to teach the metalanguage and a conservation invariant, then pivots to a **3-phase R&D pipeline** (Discovery → Preclinical → Clinical → Approved) as plain counted species — a source feeds Discovery; each phase is a timed, probabilistic transition.
- **Concepts.** `@reaction_network` DSL; `@prob_init`/`@prob_params`/`@prob_meta`; mass-action (Poisson) vs `@deterministic` rates; the timed lifecycle (`cycletime`, `probability`, `capacity`); `ReactionNetworkProblem(...; seed=)` + `simulate`; reading `prob.sol` **by column name** (construction order, not authoring order) and `prob.u`/`find_index`; seed reproducibility and a deterministically-seeded ensemble (`hash((root,k))`) for mean ± spread.
- **APIs.** `@reaction_network`, `@prob_init`, `@prob_params`, `@prob_meta`, `@deterministic`, `ReactionNetworkProblem`, `simulate`, `prob.sol`, `prob.u`, `find_index`, `mean`/`std`.
- **The number it ends on.** *"Over a 5-year horizon this pipeline delivers a mean of N approvals (p10 = X); a fourth clinical slot would/would not move it."* A capacity/throughput read a portfolio manager acts on — computed from a seeded ensemble, not a point estimate.
- **Source.** `core_engine_tour` §1 (SIR), §2 (lifecycle), §8 (determinism/ensemble), reshaped; the decision-number closer is new.
- **Target length.** ~250–350 rendered lines; the tightest, most narrated tier.
- **Status:** ✅ delivered (`docs/literate/tutorials/introductory.jl`, rendered).

### A2 — Advanced tutorial

- **Sample problem.** A structured-token R&D portfolio: projects as first-class entities with `phase`/`npv`, advanced through a lifecycle, under contended resources.
- **Concepts.** Structured/agentic tokens (`@structured_token`/`@register`, `BaseStructuredToken`); `@select`/`@advance` (phase-as-attribute, value-qualified selection); the **three resource modalities** (`@conserved`/`@rate`/`@nonblock`) and the priority allocator under contention; `population[]` declarative initial marking; ensembles and `treatment_effect`; reading portfolio KPIs. The modality truth-table and allocator (currently `core_engine_tour` §3–§4) land here, not in the intro tier.
- **APIs.** `@structured_token`, `register_structured_species!`, `@select`, `@advance`, `@conserved`/`@rate`/`@nonblock`, `priority`, `Rule`/`SetSpecies`/`SetParams`/`AddToken`, `population=`, `ensemble`, `summarize`, `treatment_effect`.
- **The number it ends on.** A treatment-effect Δ between two policies (e.g. a Series-B raise vs not) with a standard error — "the raise buys +Δ expected launches ± se."
- **Source.** `agentic_pipeline` §0–§5 + `core_engine_tour` §3–§4; ensemble/`treatment_effect` from `introspection_tour` §4.
- **Target length.** ~400–550 lines.
- **Status:** ⬜ not started.

### A3 — Expert tutorial

- **Sample problem.** The portfolio as a node in a larger heterogeneous system: composed from fragments, coupled to sibling agents, checkpointed and replayed.
- **Concepts.** Model-is-data JSON round-trip and the registry (points to the serialization deep-dive); `@compose`/`@pipeline`/`refine`/`abstract` (points to the composition deep-dive); rules and the endogenous decision channel in depth; AA wiring (`inputs[]`, `ExternalRef`, `_prestep!` one-tick Jacobi lag, `add_wire!`); `dump_state`/`restore`; the analysis/exec-map layer.
- **APIs.** `from_json_model`/`to_json_model`/`@import_model`/`@export_model`, `validate`, `refine`/`@compose`/`@pipeline`, `entangle!`/`add_wire!`/`getobservable`, `dump_state`/`restore`/`reinit!`, `exec_map`/`draw_network`.
- **The number it ends on.** A coupled-system read — e.g. the cash trajectory the finance sibling reconstructs off a wire, or a coarse-vs-refined aggregate-rNPV agreement within ensemble CI (the granularity-substitution guarantee).
- **Source.** `aa_integration` + `wires_viz_tour` (drawn) + `refinement_tour` + `agentic_pipeline` §6–§7 + `introspection_tour` §7.
- **Target length.** ~500–700 lines; may split across sub-pages.
- **Status:** ⬜ not started.

### A4 — Deep-dive: serialization *(technical focus)*

- **Scope.** A model IS data: the single eval-free JSON document; `from_json_model`/`to_json_model` round-trip; `validate` (clean vs a deliberately broken model → a `Diagnostic`, never an eval); the host-function registry (host Julia referenced by name, never carried in the document); the RCE boundary the design closes; the typed `ExprNode` IR. Structured genesis (`@structured(:Kind, …)`) as the eval-free RHS twin of `AddToken`.
- **Source.** `agentic_pipeline` §4, §6; `serialization_ir.jl` tests; ADR 0005/0006.
- **Status:** ⬜ not started.

### A5 — Deep-dive: composition & granularity *(technical focus)*

- **Scope.** The granularity ladder: `@join`/`@equalize` (manual, no-port), `@compose`/`@pipeline`/`@process` (declared open ports, matched by FK-repoint), `refine`/`abstract` (boundary-matched FK-splice, plug-compatibility), `refinement_diagnostics` (advisory Σ-ct/Π-PoS checks). The authoring-time-only invariant (forbidden on a live model).
- **Source.** `refinement_tour` + `core_engine_tour` §7; `refinement_composition.jl` tests; ADR 0009.
- **Status:** ⬜ not started.

## 6. Workstream B — Applied decision case studies

Understanding-oriented. **Question-titled, headline-number-first** (§2). In this order (flagship first, because it is the framework's raison d'être).

### B1 — "What is the marginal eNPV of the Nth scientist?" *(flagship, refined HTML)*

- **The number.** The **shadow price of the binding resource** — the marginal expected-NPV of one more scientist (or one more unit of budget), computed via ensemble runs at capacity N vs N+1. First paragraph states the number and the verdict ("hire the Nth scientist iff their fully-loaded cost < shadow price of $X").
- **The insight it carries.** Under contention the value of a resource is not its accounting cost; it is what the *system* does with it at the margin — a number a spreadsheet cannot produce.
- **Mechanics exercised (the "how we got it").** A portfolio under `@conserved` FTE and `@rate` budget contention; the ADR-0002 priority allocator; `priority`; the three modalities; `treatment_effect` across the N-vs-N+1 arms; the **exec map as the hero visual** (starvation coloring shows which resource actually binds).
- **Hero visual.** The three-layer `exec_map` (`network_graph` → `draw_network` → `exec_map`) with the binding pool painted gold.
- **Memo candidate.** Directly seeds a capacity/hiring memo.
- **Source.** New model built on `introspection_tour` (exec map, ensemble) + `core_engine_tour` (allocator, modalities); calibrated so a resource genuinely binds (per the `bd_acquisition` calibration discipline).
- **Status:** ⬜ not started.

### B2 — "What is this in-licensing asset worth to *this* pipeline?" *(BD/M&A, refined HTML)*

- **The number.** The deal's attributable **Δ-rNPV** on the *living* portfolio — with the non-obvious verdict that **rNPV is not additive under contention**: the asset's value depends on whose resources it cannibalizes. The existing 160-seed run gives the shape (S0 → S5, capability/PoS the largest synergy, cash the binding constraint); the case study leads with the headline Δ and the "value is in the platform, not the pipeline" read.
- **The insight.** A deal's worth is a *system* effect on programs you already own, not a standalone rNPV you add on.
- **Mechanics exercised.** Structured tokens; `population[]`; the endogenous acquisition `Rule` (`Seq[AddToken, SetSpecies, SetParams]`, `fire_mode=:once`); `@compose` to splice the asset in for the "combined-company-from-t=0" view; a with/without ensemble comparison (`treatment_effect`); the per-program ledger.
- **Memo candidate.** Directly reusable as a BD memo (gate-6 / salvage-value framing).
- **Source.** `demo/bd_acquisition` (migrated), reframed question-first; its `MVP_BD_DEMO.md` is the design record.
- **Status:** ⬜ scaffolded (demo exists); needs question-first reframe + HTML.

### B3 — "When should you kill a program?"

- **The number.** The kill-threshold that maximizes expected portfolio value — comparing kill-threshold `Rule` variants via `treatment_effect`, reporting the Δ between "kill below PoS θ₁" and "kill below θ₂."
- **The insight.** Demonstrates the endogenous decision channel — the capability hardest to fake in competing DES tools: the kill rule is *in the model*, state-contingent, and serializable, not host patch code.
- **Mechanics exercised.** `Rule` with a state-contingent guard; `SetTokens`/`Deactivate`/soft-retire; `treatment_effect` across threshold variants; the per-program ledger for attributing the saved capital.
- **Source.** New model built on `agentic_pipeline` (rules, structured tokens) + `introspection_tour` (`treatment_effect`).
- **Status:** ⬜ not started.

## 7. Workstream C — API / reference documentation

- **Rebuild the stale reference.** `docs/src/index.md` currently `@docs`-references removed macros (`@ReactionNetwork`, `@import_network`, `@export_network`, `@import_solution`, `@export_solution`, `@problematize`, `@plot`, `@optimize`, `@fit`, `@fit_and_plot`, `@build_solver`) and misses the entire post-rework surface. Replace with capability-organized autodoc pages driven by the real export surface.
- **Organize by capability**, mirroring the module map in [INVENTORY.md](../INVENTORY.md): *Authoring* (`@reaction_network`, `@push`, `@add_species`, `@mode`, `@aka`, `@name_transition`, cost/reward/valuation, `@prob_*`); *Structured tokens* (`@structured_token`, `register_structured_species!`, `@select`/`@advance`, `PopulationEntry`); *Rules & actions* (`Rule`, the `ActionStmt` family); *Construction & simulation* (`ReactionNetworkProblem`, `@agentize`, `simulate`, `reinit!`); *Composition* (`@join`, `@equalize`, `@compose`, `@pipeline`, `refine`, `abstract`); *Serialization* (`from_json_model`/`to_json_model`, `@import_model`/`@export_model`, `validate`, the `ExprNode` IR); *Analysis & viz* (`token_trajectory`, `ensemble`, `summarize`, `treatment_effect`, `export_run`, `network_graph`/`draw_network`/`exec_map`, the plot recipe types); *AA coupling* (`getobservable`, `add_wire!`, `ExternalRef`). Keep the attribute/shorthand tables (they are still accurate) and the rate-semantics note.
- **The JSON model schema** documented as a reference page (the §8 serialization schema, with the ADR-0005 document shape).
- **Docstring coverage.** src carries ~146 docstring blocks already; the reference is largely `@autodocs`/`@docs` assembly plus filling gaps flagged during assembly.
- **Status:** ⬜ not started.

## 8. Workstream D — Explanation layer

- **Promote the contract.** `CONTRACT_DRAFT.md` §1–§15 becomes the explanation quadrant — framed for a reader as *why the engine behaves this way* (the modality truth table, the single-clock time model, the determinism/seeding obligations, the object model, composition/serialization semantics). This is the JOSS/software-paper differentiator.
- **Surface the ADRs** as linked "decisions and rejected alternatives" reading, with the status table from `adr/README.md`.
- **Editorial note.** The contract stays normative in `spec/`; the explanation pages *link to and excerpt* it rather than forking it, to avoid a second drifting copy (same anti-drift rule as the demos).
- **Status:** ⬜ not started.

## 9. Workstream E — Build, render, and deploy

- **`docs/make.jl`** rewritten for Documenter HTML + a Literate pre-pass over `docs/literate/**`. Replaces `DocumenterMarkdown`.
- **`docs/Project.toml`** gains `Documenter`, `Literate` (and `Plots`/`Arrow`/`DataFrames`/`Distributions` for the case-study renders, mirroring the demo-local envs). `ReactiveDynamics` path-dev'd.
- **Self-contained HTML render recipe** for the case studies generalized from `demo/wires_viz_tour/build.jl` (executed Literate → markdown → inlined-SVG HTML) and `demo/bd_acquisition/build_presentation.jl`.
- **GitHub Pages deploy** via `deploydocs` (repo already targets `github.com/Merck/ReactiveDynamics.jl.git`). No CI exists in-repo (no `.github/workflows/`); document the local build+deploy command, and note CI as an optional follow-up.
- **Formatter.** All Literate `.jl` sources pass Runic (`julia -m Runic --check .`), same as the rest of the tree.
- **Status:** ⬜ not started (exemplar renders via a scoped path this PR; full site build tracked).

## 10. Progress tracker

Legend: ✅ done · 🟡 in progress · ⬜ not started.

| Facet | ID | Deliverable | Source | Status |
|---|---|---|---|---|
| Charter | — | This document | — | ✅ |
| Tutorials | A1 | Introductory (exemplar, rendered) | core_engine_tour | ✅ |
| Tutorials | A2 | Advanced | agentic_pipeline + core §3-4 | ⬜ |
| Tutorials | A3 | Expert | aa_integration + wires + refinement | ⬜ |
| Tutorials | A4 | Deep-dive: serialization | agentic_pipeline §4,6 | ⬜ |
| Tutorials | A5 | Deep-dive: composition | refinement_tour + core §7 | ⬜ |
| Case study | B1 | Marginal eNPV of the Nth scientist (flagship, HTML) | new + introspection_tour | ⬜ |
| Case study | B2 | In-licensing asset value (BD/M&A, HTML) | bd_acquisition | 🟡 (demo exists) |
| Case study | B3 | When to kill a program | new + agentic_pipeline | ⬜ |
| Reference | C1 | Capability-organized API pages | src docstrings | ⬜ |
| Reference | C2 | JSON model schema page | CONTRACT §8 / ADR 0005 | ⬜ |
| Explanation | D1 | Contract-as-explanation pages | CONTRACT_DRAFT.md | ⬜ |
| Explanation | D2 | ADR reading surface | spec/adr | ⬜ |
| Build | E1 | make.jl (Documenter + Literate) | wires_viz_tour/build.jl | ⬜ |
| Build | E2 | docs/Project.toml + deploy | — | ⬜ |
| Build | E3 | Refined-HTML render recipe (case studies) | bd_acquisition/build_presentation.jl | ⬜ |

## 11. Sequencing

1. **Charter + exemplar (this PR increment).** This document + the introductory tutorial (A1) authored in the target Literate style and rendered end-to-end, as the quality bar. *(done)*
2. **Toolchain (E1/E2).** Wire `make.jl` + `docs/Project.toml` so the full Documenter+Literate site builds with A1 in place. Land the reference (C1/C2) and explanation (D1/D2) scaffolds so the site is navigable.
3. **Tutorial tiers (A2/A3) + deep-dives (A4/A5).** Migrate the remaining demos in.
4. **Case studies (B1/B2/B3).** Flagship first; refined HTML for B1/B2 (E3).
5. **Polish + deploy.** Runic pass, link audit, GitHub Pages deploy, README refresh to point at the published site.

Each numbered step is a reviewable batch; the charter's tracker (§10) is updated as facets land.

## 12. Acceptance criteria (per page)

A tutorial or case-study page is *done* when: every code block executes clean during the site build; it ends in (tutorial) or leads with (case study) a manager-actionable number; the page's constructs each trace to a passing semantic test or a docstring; the run is `seed=`-pinned and reproducible; the source passes Runic; and — for case studies — the title is the question, and the mechanics read as "how we got it," not as the subject. A reference page is done when every documented symbol resolves against the current export surface (the `exports_resolve.jl` invariant) and no removed symbol is referenced.

## 13. Decisions log

- **2026-07-17 — Migrate demos in place** (single source of truth), not keep-alongside or absorb-and-retire. Rationale: no duplicated model code to drift; demos stay independently runnable.
- **2026-07-17 — This-PR scope: charter + one rendered exemplar** (the introductory tutorial), not full scaffolding. The exemplar sets the quality bar; the rest is tracked in §10 and sequenced in §11.
- **2026-07-17 — Introductory tier is a strict subset** (plain-species timed pipeline only); the modality truth-table and the allocator move to the advanced tier, per the PR-draft tier mapping.
