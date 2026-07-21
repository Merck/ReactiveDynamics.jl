# Documentation & Tutorials Charter — ReactiveDynamics.jl

> **The tracking spine for the documentation rework.** A standalone PR into `rework` (branch `docs-tutorials`). This is the durable "what are we building, in what order, and how do we know a piece is done" index for the docs effort, the counterpart to [STATUS.md](STATUS.md) for the engine. Last updated 2026-07-21 (Workstreams A2–A5, B1–B3, C1/C2 landed; the full site builds green. Workstream D re-scoped 2026-07-21: the web explanation quadrant D1/D2 — contract-as-explanation pages and an ADR reading surface — is **dropped**; the reader-facing "why" is carried instead by two peer papers, and the site links directly to the normative `spec/` artifacts. The single arXiv paper is now split into **two peer papers** — D3a methodology, D3b an HBR-style adoption/value paper for non-technical executives.). Re-verify any `file:line` before acting — line numbers drift.

## 1. Why this exists and what "done" means

The engine rework (§1–§15 of [CONTRACT_DRAFT.md](CONTRACT_DRAFT.md), ADRs 0001–0015) is implemented and green (801 pass / 0 broken). What was missing was the reader-facing surface: onboarding tutorials, applied case studies, and an API reference matching the post-ADR-0015 API. As of 2026-07-18 that surface is largely built: the landing page is now an intent-router, `make.jl` runs a Documenter HTML + Literate build (replacing the deprecated `DocumenterMarkdown` backend), the seven runnable tours under [demo/](../demo) have been migrated in place as the site's Literate sources, and the reference is capability-organized off the real export surface. (The historical starting point: `docs/src/index.md` was a single stale API page referencing removed macros — `@ReactionNetwork`, `@optimize`, `@fit`, `@problematize`, `@plot`, `@import_network`, … — and the demos were not wired into a published site.) What remains is Workstream D (the two papers) and E3/E4 (refined case-study HTML + the paper builds).

This PR closes that gap. **Done** for the whole PR means: a Documenter site, publishable to GitHub Pages, comprising three tiered tutorials, two focused deep-dives, three applied case studies, and a complete API reference — with the reader-facing "why" carried by two peer papers (a methodology paper and an executive adoption/value paper) that the site links to, and the normative contract and ADRs linked directly at their `spec/` home. Every code block is executed against the current engine, every tutorial ends in a decision-relevant computed quantity, and the `demo/` sources are migrated in as the single source of truth (no duplicated model code). The per-facet acceptance criteria are in §12; the progress tracker is in §10.

## 2. Design principles (non-negotiable)

These are the brand and quality commitments the whole effort is measured against. They come from the maintainer framing in the PR draft ([PR_DRAFT.md](PR_DRAFT.md) §"Documentation and tutorials") and are binding on every page.

- **Every applied example is a decision case study with a headline number, not a feature tour.** Title each case study as the *question it answers*; put the number and the verdict in the first paragraph; let the framework mechanics appear as "how we got it." This converts docs effort into artifacts that share the skeleton of a decision memo (question → number → implication) so a case study can *seed* a memo rather than compete with it. A feature-organized gallery ("modalities demo", "composition demo") signals tool-builder; a question-organized one signals operator. Same content, different brand transfer.
- **Diátaxis, with the "why" carried by the papers, not re-hosted on the web.** Tutorials are learning-oriented; case studies are understanding-oriented; API docs are reference. The explanation quadrant — the *why the engine behaves this way* — is **not** rebuilt as a set of web pages (that was D1/D2, dropped 2026-07-21); it is carried by the normative `CONTRACT_DRAFT.md` + ADRs (linked at their `spec/` home) and, in scholarly and executive registers, by the two peer papers (§8). This is the differentiator versus every ad-hoc discrete-event-simulation package and the strongest argument in a JOSS / software-paper review. See §3.
- **The decision-relevant-quantity discipline applies to *every* tutorial, even the introductory one.** Each tutorial ends by *computing* a quantity that supports an informed decision — a marginal effect, a shadow price, a distribution with its spread — not a recap of features. Frame it as information a decision rests on (what the number *is* and what it implies), in a computational-framework register — **not** as an aspirational "answer for a manager." The payoff is a computed result the framework produces and a spreadsheet cannot, stated plainly. This keeps the learning-oriented tier honest to the decision-modeling purpose without leaning on role-flattering language.
- **Single source of truth: migrate the demos in place.** The `demo/*.jl` literate sources become the Literate.jl inputs the site ingests. Model-building code is written once; there is no parallel copy in `docs/` to drift. See §4 for the demo→docs mapping.
- **No invented API — an internal acceptance rule, never reader-facing.** Every construct on every page must trace to a passing semantic test (`test/semantic/*.jl`) or an engine docstring, and this is verified in review/CI (see §12). But this is an *internal quality gate*: do **not** surface it to the reader. Tutorial prose must not contain internal-QA meta ("every construct here is covered by the test suite", "this API is stable", "the tutorial invents no API") — that is the maintainer's assurance to check in the background, not part of the reader's learning path. When docs and code disagree, the code is truth (per [CLAUDE.md](../CLAUDE.md)).
- **Executable and reproducible.** Every code block runs during the site build (Literate `@example`/executed markdown); every run is pinned by an explicit `seed=`, so the rendered numbers are reproducible from `(model, seed)` (CONTRACT §4).

## 3. Diátaxis structure

| Quadrant | Orientation | What lives here | Source |
|---|---|---|---|
| **Tutorials** | learning | Three tiered onboarding tutorials (introductory / advanced / expert) + two focused technical deep-dives (serialization, composition). End-to-end modeling workflows at increasing depth. | Workstream A (§5) |
| **Case studies** (how-to / understanding) | understanding | Three applied decision case studies, question-titled, each with a headline number. The flagship gets a refined HTML presentation. | Workstream B (§6) |
| **Reference** | information | The complete public API, autodoc-generated from docstrings, organized by capability; the attribute/shorthand tables; the JSON model schema. | Workstream C (§7) |
| **Explanation** | understanding | The *why the engine behaves this way* is carried by the normative operational-semantics contract (§1–§15) and the ADRs — linked at their `spec/` home, not re-hosted as web pages — **and by two peer papers**: a methodology paper (scholarly narrative of the semantics and backing concepts, arXiv-first) and an HBR-style adoption/value paper for non-technical executives (motivation → solution → examples). | Workstream D (§8) |

The site's landing page routes a reader by intent: "new here" → introductory tutorial; "what can it do for my decision" → case studies; "how do I call X" → reference; "why does it behave this way" → the normative contract/ADRs and the two papers.

## 4. Toolchain, mechanics, and the migrate-in-place rule

**Build stack.** [Documenter.jl](https://documenter.juliadocs.org) for the site (HTML backend, GitHub Pages deploy) + [Literate.jl](https://fredrikekre.github.io/Literate.jl) for the executable tutorial/case-study sources. This replaces the deprecated `DocumenterMarkdown` backend currently in `docs/make.jl`. Literate sources under `docs/literate/**` are processed to `docs/src/**` markdown with `@example` blocks that Documenter executes at build time, so code runs end-to-end and its output is captured into the page. The one demo already on Literate ([demo/wires_viz_tour](../demo/wires_viz_tour), with its `build.jl` render recipe) is the working precedent for the executed-Literate → HTML path and the SVG-inlining trick for self-contained pages.

**Julia documentation conventions to follow (from a survey of Catalyst, SciML/DifferentialEquations, Agents.jl, JuMP, ModelingToolkit, and the Documenter/Literate guides).** These are how prominent Julia computational-framework packages actually write docs; the tutorials adopt them so RD reads as a native Julia framework:

- **Register.** Technical, first-person-plural instructional voice ("we author…", "let's run…"). Application framing is fine as a *hook* and a *closer*, but the payoff is a **computed quantity or observed behavior the reader just produced**, never a marketing/aspirational claim (this is why "answer for a manager" is out — §2). JuMP is the austere end (numbers only); Catalyst/SciML close on a rendered result.
- **Plots as the payoff.** For a trajectory-producing engine, an inline figure is the idiomatic close (Catalyst/SciML render `plot(sol)`); when the natural output is a scalar, showing the number is legitimate (JuMP). RD tutorials render inline figures *and* report the decision quantity.
- **Intro shape.** Catalyst (the closest analogue — a reaction-network DSL) builds cumulatively around one running model; JuMP/SciML/Agents lead with a complete copy-pasteable example then deconstruct. RD's introductory tier follows the **cumulative buildup** (decision 2026-07-17), consistent with Catalyst.
- **Theory goes in its own section, never inline in a tutorial.** The "why" (Petri-net / discrete-event semantics, rNPV rationale) lives in the normative contract/ADRs and the two papers (Workstream D), not inline — matching JuMP's *Background* and ModelingToolkit's *Internals* split. Tutorials link out to those; they don't embed the theory.
- **No internal-QA meta in reader prose** (confirmed universal across all six packages surveyed) — see §2. *User-facing* reproducibility/gotcha notes ("thread the seed; never touch the global RNG") ARE idiomatic and welcome. The "traces to a test" gate is enforced via Literate `#src` hidden self-test lines (which never reach the reader) and review/CI, not prose.
- **Diátaxis without naming it.** None of the packages advertise "we use Diátaxis" (JuMP's style guide cites its ancestor, Divio); they just implement the four-way split. RD does the same — §3 structures by Diátaxis, but reader-facing nav is by intent, not by quadrant label.
- **Literate mechanics.** `# ## Heading` for prose headings, `##` for real code comments, last-expression-shows-output, `savefig`+`![]()` or a returned plot object for figures, `@setup`/`#src` to hide boilerplate and self-tests. Standard `docs/src/` + `docs/make.jl` layout, `index.md` as Home.

**Target docs tree** (the structure the charter builds toward; `[x]` = exists after this PR's exemplar, `[ ]` = tracked):

```
docs/
  make.jl                        # [x] Documenter + Literate site build (replaces DocumenterMarkdown)
  build_literate.jl              # [x] standalone executed-Literate → self-contained HTML (fast authoring/preview)
  Project.toml                   # [x] Documenter + Literate + Plots (inline figures; also triggers RDPlotsExt); Arrow/DataFrames later
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
    assets/                      # existing diagram*.png etc. (audit for staleness)
                                 # NB: no explanation/ web pages — D1/D2 dropped 2026-07-21; the site links
                                 #     directly to spec/CONTRACT_DRAFT.md, spec/adr/, and the two papers below
paper/                           # [ ] the two peer papers (Workstream D) — format/home TBD at authoring
  methodology/                   # [ ] D3a — scholarly semantics/computational-concepts paper (arXiv-first)
  adoption/                      # [ ] D3b — HBR-style adoption/value paper for non-technical executives
                                 #     (LaTeX+arXiv vs Markdown+pandoc; see §8). NORMATIVE source stays spec/CONTRACT_DRAFT.md
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
| `bd_acquisition` | scripts + JSON model | **flagship + in-licensing case studies** (refined HTML) | Its `MVP_BD_DEMO.md` design doc defines the decision-case-study DNA. |

**Refined HTML for case studies.** The flagship case study (and the in-licensing one) get a refined, self-contained HTML presentation beyond the Documenter page: client-side SVG charts from a real ensemble dump, narrative-led. A headless-Chrome + print-stylesheet pass remains available for a PDF artifact.

## 5. Workstream A — Tiered onboarding tutorials

Learning-oriented, literate, didactic, on sample problems, framed as end-to-end modeling workflows at increasing depth. **Each ends in a manager-actionable number** (the discipline in §2). Tier boundaries follow the maintainer's mapping in the PR draft.

### A1 — Introductory tutorial *(EXEMPLAR delivered this PR)*

- **Sample problem.** A plain-species timed pipeline. Opens on the classic **SIR** epidemic to teach the metalanguage and a conservation invariant, then pivots to a **3-phase R&D pipeline** (Discovery → Preclinical → Clinical → Approved) as plain counted species — a source feeds Discovery; each phase is a timed, probabilistic transition.
- **Concepts.** `@reaction_network` DSL; `@prob_init`/`@prob_params`/`@prob_meta`; mass-action (Poisson) vs `@deterministic` rates; the timed lifecycle (`cycletime`, `probability`, `capacity`); `ReactionNetworkProblem(...; seed=)` + `simulate`; reading `prob.sol` **by column name** (construction order, not authoring order) and `prob.u`/`find_index`; seed reproducibility and a deterministically-seeded ensemble (`hash((root,k))`) for mean ± spread.
- **APIs.** `@reaction_network`, `@prob_init`, `@prob_params`, `@prob_meta`, `@deterministic`, `ReactionNetworkProblem`, `simulate`, `prob.sol`, `prob.u`, `find_index`, `mean`/`std`.
- **The quantity it ends on.** The marginal effect of a fourth clinical trial slot on approvals over the horizon (≈ +4, computed from a seeded ensemble with its standard error) — a marginal, system-level quantity that informs a capacity decision, not a point estimate. Framed in a computational register (what the number is and what it implies), not as an "answer for a manager".
- **Hierarchy.** The two model sections (SIR, then the timed pipeline) are the spine; seeding/ensembles is a *supporting* subsection (§3, "running the model many times"), not a peer of the models; the decision-quantity is the payoff (§4).
- **Figures.** Inline plots as the Julia-idiomatic payoff: the SIR S/I/R curves after §1, and the 3-slot-vs-4-slot approval-distribution histogram at §4 (alongside the reported number). Rendered via `Plots` (a docs-env dep; also triggers `RDPlotsExt`).
- **Source.** `core_engine_tour` §1 (SIR), §2 (lifecycle), §8 (determinism/ensemble), reshaped; the decision-quantity closer is new.
- **Target length.** ~250–350 rendered lines; the tightest, most narrated tier.
- **Status:** ✅ delivered (`docs/literate/tutorials/introductory.jl`, rendered).

### A2 — Advanced tutorial

- **Sample problem.** A structured-token R&D portfolio: projects as first-class entities with `phase`/`npv`, advanced through a lifecycle, under contended resources.
- **Concepts.** Structured/agentic tokens (`@structured_token`/`@register`, `BaseStructuredToken`); `@select`/`@advance` (phase-as-attribute, value-qualified selection); the **three resource modalities** (`@conserved`/`@rate`/`@nonblock`) and the priority allocator under contention; `population[]` declarative initial marking; ensembles and `treatment_effect`; reading portfolio KPIs. The modality truth-table and allocator (currently `core_engine_tour` §3–§4) land here, not in the intro tier.
- **APIs.** `@structured_token`, `register_structured_species!`, `@select`, `@advance`, `@conserved`/`@rate`/`@nonblock`, `priority`, `Rule`/`SetSpecies`/`SetParams`/`AddToken`, `population=`, `ensemble`, `summarize`, `treatment_effect`.
- **The number it ends on.** A treatment-effect Δ between two policies (e.g. a Series-B raise vs not) with a standard error — "the raise buys +Δ expected launches ± se."
- **Source.** `agentic_pipeline` §0–§5 + `core_engine_tour` §3–§4; ensemble/`treatment_effect` from `introspection_tour` §4.
- **Target length.** ~400–550 lines.
- **Status:** ✅ delivered (`docs/literate/tutorials/advanced.jl`, rendered; closes on a Series-B `treatment_effect` Δ +0.52 launches ± 0.24 SE). NB the migration surfaced a stale demo assumption: the `@rate(cycletime=0)` combination is now REJECTED at construction (the validator hard-errors), so the "silent free input" prose was corrected.

### A3 — Expert tutorial

- **Sample problem.** The portfolio as a node in a larger heterogeneous system: composed from fragments, coupled to sibling agents, checkpointed and replayed.
- **Concepts.** Model-is-data JSON round-trip and the registry (points to the serialization deep-dive); `@compose`/`@pipeline`/`refine`/`abstract` (points to the composition deep-dive); rules and the endogenous decision channel in depth; AA wiring (`inputs[]`, `ExternalRef`, `_prestep!` one-tick Jacobi lag, `add_wire!`); `dump_state`/`restore`; the analysis/exec-map layer.
- **APIs.** `from_json_model`/`to_json_model`/`@import_model`/`@export_model`, `validate`, `refine`/`@compose`/`@pipeline`, `entangle!`/`add_wire!`/`getobservable`, `dump_state`/`restore`/`reinit!`, `exec_map`/`draw_network`.
- **The number it ends on.** A coupled-system read — e.g. the cash trajectory the finance sibling reconstructs off a wire, or a coarse-vs-refined aggregate-rNPV agreement within ensemble CI (the granularity-substitution guarantee).
- **Source.** `aa_integration` + `wires_viz_tour` (drawn) + `refinement_tour` + `agentic_pipeline` §6–§7 + `introspection_tour` §7.
- **Target length.** ~500–700 lines; may split across sub-pages.
- **Status:** ✅ delivered (`docs/literate/tutorials/expert.jl`, rendered; closes on the finance sibling reconstructing RD's cash off a wire, verifying the one-tick Jacobi lag exactly). Delivered as one page (~436 lines), not split.

### A4 — Deep-dive: serialization *(technical focus)*

- **Scope.** A model IS data: the single eval-free JSON document; `from_json_model`/`to_json_model` round-trip; `validate` (clean vs a deliberately broken model → a `Diagnostic`, never an eval); the host-function registry (host Julia referenced by name, never carried in the document); the RCE boundary the design closes; the typed `ExprNode` IR. Structured genesis (`@structured(:Kind, …)`) as the eval-free RHS twin of `AddToken`.
- **Source.** `agentic_pipeline` §4, §6; `serialization_ir.jl` tests; ADR 0005/0006.
- **Status:** ✅ delivered (`docs/literate/deep_dives/serialization.jl`, rendered; closes on byte-identical DSL-vs-JSON-reload trajectories + the broken model yielding a `Diagnostic` returned as data, never an eval).

### A5 — Deep-dive: composition & granularity *(technical focus)*

- **Scope.** The granularity ladder: `@join`/`@equalize` (manual, no-port), `@compose`/`@pipeline`/`@process` (declared open ports, matched by FK-repoint), `refine`/`abstract` (boundary-matched FK-splice, plug-compatibility), `refinement_diagnostics` (advisory Σ-ct/Π-PoS checks). The authoring-time-only invariant (forbidden on a live model).
- **Source.** `refinement_tour` + `core_engine_tour` §7; `refinement_composition.jl` tests; ADR 0009.
- **Status:** ✅ delivered (`docs/literate/deep_dives/composition.jl`, rendered; closes on the granularity-substitution guarantee — a plug-compatible sub-model's Σ-ct/Π-PoS matches the coarse transition and `refinement_diagnostics` passes clean, warns on a drifted sub-model).

## 6. Workstream B — Applied decision case studies

Understanding-oriented. **Question-titled, headline-number-first** (§2). In this order (flagship first, because it is the framework's raison d'être).

### B1 — "What is the marginal eNPV of the Nth scientist?" *(flagship, refined HTML)*

- **The number.** The **shadow price of the binding resource** — the marginal expected-NPV of one more scientist (or one more unit of budget), computed via ensemble runs at capacity N vs N+1. First paragraph states the number and the verdict ("hire the Nth scientist iff their fully-loaded cost < shadow price of $X").
- **The insight it carries.** Under contention the value of a resource is not its accounting cost; it is what the *system* does with it at the margin — a number a spreadsheet cannot produce.
- **Mechanics exercised (the "how we got it").** A portfolio under `@conserved` FTE and `@rate` budget contention; the ADR-0002 priority allocator; `priority`; the three modalities; `treatment_effect` across the N-vs-N+1 arms; the **exec map as the hero visual** (starvation coloring shows which resource actually binds).
- **Hero visual.** The three-layer `exec_map` (`network_graph` → `draw_network` → `exec_map`) with the binding pool painted gold.
- **Memo candidate.** Directly seeds a capacity/hiring memo.
- **Source.** New model built on `introspection_tour` (exec map, ensemble) + `core_engine_tour` (allocator, modalities); calibrated so a resource genuinely binds (per the `bd_acquisition` calibration discipline).
- **Status:** ✅ page delivered (`docs/literate/case_studies/marginal_scientist.jl`, rendered; leads with the shadow price of the 5th scientist ≈ +$19M ± $2.5M, the scientist bench genuinely binds, exec-map hero visual paints the binding pool). The refined-HTML presentation is tracked separately as E3.

### B2 — "What is this in-licensing asset worth to *this* pipeline?" *(BD/M&A, refined HTML)*

- **The number.** The deal's attributable **Δ-rNPV** on the *living* portfolio — with the non-obvious verdict that **rNPV is not additive under contention**: the asset's value depends on whose resources it cannibalizes. The existing 160-seed run gives the shape (S0 → S5, capability/PoS the largest synergy, cash the binding constraint); the case study leads with the headline Δ and the "value is in the platform, not the pipeline" read.
- **The insight.** A deal's worth is a *system* effect on programs you already own, not a standalone rNPV you add on.
- **Mechanics exercised.** Structured tokens; `population[]`; the endogenous acquisition `Rule` (`Seq[AddToken, SetSpecies, SetParams]`, `fire_mode=:once`); `@compose` to splice the asset in for the "combined-company-from-t=0" view; a with/without ensemble comparison (`treatment_effect`); the per-program ledger.
- **Memo candidate.** Directly reusable as a BD memo (gate-6 / salvage-value framing).
- **Source.** `demo/bd_acquisition` (migrated), reframed question-first; its `MVP_BD_DEMO.md` is the design record.
- **Status:** ✅ page delivered (`docs/literate/case_studies/inlicensing_value.jl`, rendered; leads with the attributable Δ-rNPV +1947 ± 367 on a reduced 40-seed ensemble and the not-additive-under-contention verdict — cash binds, synergies are super-additive; honestly notes the reduced seed count vs the 160-seed demo). The refined-HTML presentation is tracked separately as E3.

### B3 — "When should you kill a program?"

- **The number.** The kill-threshold that maximizes expected portfolio value — comparing kill-threshold `Rule` variants via `treatment_effect`, reporting the Δ between "kill below PoS θ₁" and "kill below θ₂."
- **The insight.** Demonstrates the endogenous decision channel — the capability hardest to fake in competing DES tools: the kill rule is *in the model*, state-contingent, and serializable, not host patch code.
- **Mechanics exercised.** `Rule` with a state-contingent guard; `SetTokens`/`Deactivate`/soft-retire; `treatment_effect` across threshold variants; the per-program ledger for attributing the saved capital.
- **Source.** New model built on `agentic_pipeline` (rules, structured tokens) + `introspection_tour` (`treatment_effect`).
- **Status:** ✅ delivered (`docs/literate/case_studies/kill_a_program.jl`, rendered; leads with the value-maximizing kill threshold θ\*=0.5 worth +48.5 ± 12.6 vs never-killing, with a genuine interior optimum — the endogenous kill Rule lives in the model, state-contingent and serializable).

## 7. Workstream C — API / reference documentation

- **Rebuild the stale reference.** `docs/src/index.md` currently `@docs`-references removed macros (`@ReactionNetwork`, `@import_network`, `@export_network`, `@import_solution`, `@export_solution`, `@problematize`, `@plot`, `@optimize`, `@fit`, `@fit_and_plot`, `@build_solver`) and misses the entire post-rework surface. Replace with capability-organized autodoc pages driven by the real export surface.
- **Organize by capability**, mirroring the module map in [INVENTORY.md](INVENTORY.md): *Authoring* (`@reaction_network`, `@push`, `@add_species`, `@mode`, `@aka`, `@name_transition`, cost/reward/valuation, `@prob_*`); *Structured tokens* (`@structured_token`, `register_structured_species!`, `@select`/`@advance`, `PopulationEntry`); *Rules & actions* (`Rule`, the `ActionStmt` family); *Construction & simulation* (`ReactionNetworkProblem`, `@agentize`, `simulate`, `reinit!`); *Composition* (`@join`, `@equalize`, `@compose`, `@pipeline`, `refine`, `abstract`); *Serialization* (`from_json_model`/`to_json_model`, `@import_model`/`@export_model`, `validate`, the `ExprNode` IR); *Analysis & viz* (`token_trajectory`, `ensemble`, `summarize`, `treatment_effect`, `export_run`, `network_graph`/`draw_network`/`exec_map`, the plot recipe types); *AA coupling* (`getobservable`, `add_wire!`, `ExternalRef`). Keep the attribute/shorthand tables (they are still accurate) and the rate-semantics note.
- **The JSON model schema** documented as a reference page (the §8 serialization schema, with the ADR-0005 document shape).
- **Docstring coverage.** src carries ~146 docstring blocks already; the reference is largely `@autodocs`/`@docs` assembly plus filling gaps flagged during assembly.
- **Status:** ✅ delivered (C1 + C2). Nine capability-organized `@docs` pages under `docs/src/reference/` (authoring, structured tokens, rules & actions, construction & simulation, composition, serialization, analysis & viz, AA coupling) plus the JSON model-schema page; AA verbs (`simulate`/`reinit!`/`entangle!`/`add_wire!`/`getobservable`) documented in prose as reexported AlgebraicAgents functions; deprecated symbols excluded. Every `@docs` symbol resolves against the current export surface and the package precompiles clean. The gap-fill added docstrings to the closed action family + `Rule`, the `ReactionNetworkProblem` live-state type + constructor, the eval-free `ExprNode` IR + JSON (de)serializer, and `@append_transitions` (committed `0761cf3`; a concurrent src-review pass `dd8ce5d` backfilled the remaining legacy-file docstrings).

## 8. Workstream D — The two papers

**Re-scoped 2026-07-21.** The earlier plan rebuilt the "why" as a web explanation quadrant (D1 contract-as-explanation pages, D2 an ADR reading surface) anchored by a single arXiv paper. That web layer is **dropped**: re-hosting the contract and ADRs as Documenter pages is a third copy that drifts from the normative `spec/` source for no reader payoff the papers don't cover better. Instead the reader-facing "why" is carried by **two peer papers**, and the site links directly to the normative artifacts (`spec/CONTRACT_DRAFT.md`, `spec/adr/`) at their canonical home on GitHub.

The two papers are peers, not a long-and-short of one text: they lead differently, target different readers, and can be socialized independently.

- **D3a — the methodology paper.** Technical / computational-concepts focused: the operational semantics, the *why* behind them, and the backing concepts, as scholarly narrative. Leads with the applied use cases as evidence (the three Workstream-B case studies), in the register of a methods/software paper. arXiv-first; downstream feeds a thin JOSS `paper.md` and/or a peer-reviewed methods paper.
- **D3b — the adoption / value paper (HBR-style).** For technically non-savvy executives: it does **not** lead with the use cases. It builds the motivation (why R&D/business-process decisions are dynamic, contended, and decision-laden — not a spreadsheet), then the solution (what the framework brings), up to *some* illustrative examples drawn from the same applied work, in the register of a Harvard Business Review article. The artifact to hand an exec to explain "what we bring."

Both are the same substance seen from opposite ends: D3a leads with mechanism and lands on value; D3b leads with the problem and value and lands on *some* mechanism. Neither forks the normative contract (§8.4).

### D3a — Methodology paper (arXiv-first, technical)

- **Purpose.** A scholarly narrative arguing the operational semantics and backing concepts — the artifact to *share for socialization* with a technical audience, and to cite. Longer than a JOSS paper (which is ~250–1000 words of Summary + Statement of need); here the semantics themselves are the substance.
- **Fixed outline** (the substance exists as `CONTRACT_DRAFT.md` §1–§15 + the ADRs; the work is reframing it as scholarly narrative — the *why* and the backing concepts, not the normative *what*):
  1. **Abstract.**
  2. **Introduction** — modeling R&D/business processes as timed, stochastic, resource-contended, *decision-laden* systems; rNPV / what-if; why the question is dynamic, not a spreadsheet.
  3. **Statement of need & positioning** — vs ad-hoc DES (SimPy-style), vs chemical reaction networks (Gillespie/Catalyst), vs system dynamics, vs spreadsheet rNPV; the DyVE + AlgebraicAgents lineage. *(Also seeds a later JOSS statement-of-need.)*
  4. **The conceptual model** — transitions as stateful recipes, species as resources, the ontology; why it is *not* a CRN despite the DSL surface.
  5. **Operational semantics** — the formal core: single discrete clock, the ordered per-tick step, instance lifecycle, the modality truth table, the invariants (CONTRACT §1–§3).
  6. **Determinism & reproducibility** — the `(model, seed)` contract, RNG threading, ensembles (§4). A genuine differentiator — most DES tools are informal here.
  7. **Resource allocation under contention** — the priority-weighted progressive-fill allocator (ADR 0002 / §1.5).
  8. **Models as data** — the eval-free typed `ExprNode` IR + single-JSON serialization, the RCE boundary, registry-by-name; why it matters for agentic authoring/exchange (§8, ADR 0005/0006).
  9. **Structured/agentic tokens** — first-class entities, filtration, phase-as-attribute (§9, ADR 0006/0008).
  10. **The endogenous decision channel** — rules, guards, conditional transitions; decisions *in* the model (§12, ADR 0010/0011).
  11. **Composition & refinement** — open ports, boundary-matched FK-splice refinement, granularity substitution, the algebraic properties (§7/§11, ADR 0009).
  12. **Integration** — RD as an AlgebraicAgents node; coupled heterogeneous simulation, the one-tick Jacobi lag (§13, ADR 0012).
  13. **Worked case studies** — the three Workstream-B decision case studies as evidence the semantics deliver decision value (shared skeleton: question → number → implication).
  14. **Related work.**
  15. **Discussion, limitations & deferred work** — the honest deferrals (entity-level refinement, threaded ensemble backend, `Opera` implicit coupling, `dump_state` in-flight constraint).
  16. **Availability & reproducibility** — license (MIT), the semantic test suite, `(model, seed)` reproducibility.
  17. **References.**
- **Status:** ⬜ not started (scoped, format TBD).

### D3b — Adoption / value paper (HBR-style, executive)

- **Purpose.** A persuasive, largely non-technical narrative for executives and decision-makers who will *not* read the methodology paper — the artifact that explains "what we bring" and why it matters, shareable inside an organization. Register: a Harvard Business Review article (concrete, narrative, one core idea per section, minimal notation), not a tutorial or a feature tour.
- **Shape (leads with the problem, not the use cases).** Unlike the case studies (question-first) and D3a (mechanism-first), D3b earns the examples:
  1. **The hook** — a familiar decision an organization gets wrong with a spreadsheet: a resource everyone treats as its accounting cost, a deal valued standalone, a program killed too late or too early.
  2. **Why the usual tools fail** — rNPV-on-a-spreadsheet and static models ignore time, contention, and the decisions embedded in the process; the value of a resource is what the *system* does with it at the margin.
  3. **The idea** — model the organization as a living, timed, resource-contended, decision-laden process; simulate it; read decision-relevant quantities (shadow prices, attributable Δ-value, value-maximizing thresholds) a spreadsheet cannot produce.
  4. **What the framework brings** — in plain language: structured entities with identity, resource modalities/contention, decisions in the model, composition, reproducibility — each tied to a decision it unlocks, not to an API.
  5. **Some examples** — *illustrative* distillations of the three case studies (the Nth-scientist shadow price; in-licensing value that is not additive under contention; the kill-threshold), each as a short "the number and what it changed" vignette — a taste, not the full worked study.
  6. **What it takes to adopt** — how this fits an organization's decision process (seeds a decision memo), what is and isn't in scope, honest limits.
  7. **Close** — the one-sentence "what we bring."
- **Source.** The same applied work as D3a §13 and Workstream B, re-narrated executive-first; the `bd_acquisition` presentation's HBR-style narrative is a register precedent, and its refined-HTML presentation (E3) is a natural companion artifact for this paper.
- **Status:** ⬜ not started (scoped, format TBD).

### D3 — Shared format decision

- **Format/home — DEFERRED to authoring time** (decision 2026-07-17, unchanged by the 2026-07-21 split). Two candidates per paper, chosen when authoring starts: **(a) LaTeX under `paper/<name>/`, arXiv-primary** — best typesetting/math, direct `.tar.gz` submission, the site links to the PDF; **(b) Markdown** — single source, PDF generated via pandoc→LaTeX (more brittle for math/refs). D3a leans (a) (math-heavy); D3b may take (b) (prose-heavy, few equations). The outlines above and the artifact roles (§8.4) hold regardless of format.

### D-future — downstream venues (fed by the methodology paper, not scheduled)

- **JOSS `paper.md`** — a thin (~750-word) Summary + Statement-of-need vehicle for a citable DOI, distilled from D3a §2–§3. Prerequisites JOSS will check: CI (a `.github/workflows/` scaffold now exists on this branch — see `6255579`), a tagged release + archival DOI (Zenodo), and community guidelines (`CONTRIBUTING` — now present). Tracked as gates, not committed work.
- **Peer-reviewed methods paper** (JSS / SoftwareX / a decision-sciences or DES domain venue) — where the operational-semantics contract *is* the reviewed substance, expanded from D3a with a comparison against ad-hoc DES packages. A larger effort; parked.

### D4 — Artifact roles (the anti-fork rule)

Three roles, non-overlapping, so nothing drifts:

- **`spec/CONTRACT_DRAFT.md` (+ `spec/adr/`) stays NORMATIVE** — the source of truth for engine *behavior*; both papers and the docs site cite/link it, never restate it normatively. There is no web copy of it (D1/D2 dropped).
- **The two papers are the scholarly and executive narratives built on the contract** — D3a argues the semantics for a technical audience; D3b argues the value for executives. Both summarize and motivate; neither is normative.
- **The docs site links out** — to the normative contract/ADRs on GitHub and to the two papers — rather than hosting a third copy of the "why."

## 9. Workstream E — Build, render, and deploy

- **`docs/make.jl`** rewritten for Documenter HTML + a Literate pre-pass over `docs/literate/**`. Replaces `DocumenterMarkdown`.
- **`docs/Project.toml`** gains `Documenter`, `Literate` (and `Plots`/`Arrow`/`DataFrames`/`Distributions` for the case-study renders, mirroring the demo-local envs). `ReactiveDynamics` path-dev'd.
- **Self-contained HTML render recipe** for the case studies generalized from `demo/wires_viz_tour/build.jl` (executed Literate → markdown → inlined-SVG HTML).
- **GitHub Pages deploy** via `deploydocs` (repo already targets `github.com/Merck/ReactiveDynamics.jl.git`). No CI exists in-repo (no `.github/workflows/`); document the local build+deploy command, and note CI as an optional follow-up (also a JOSS gate, §8 D-future).
- **Paper builds (E4)** — when each paper's format lands: either `latexmk` over `paper/<name>/*.tex` (LaTeX route) or a pandoc→LaTeX pass over the Markdown route, per paper (D3a and D3b build independently). Deferred with the format decision (§8 D3).
- **Formatter.** All Literate `.jl` sources pass Runic (`julia -m Runic --check .`), same as the rest of the tree.
- **Status:** 🟡 E1/E2 fully landed — `make.jl` runs the Literate pre-pass over all eight tutorial/deep-dive/case-study sources and wires the complete nav (Tutorials, Case studies, Reference); `docs/Project.toml` carries the demo-union deps; the **full site builds green**. The in-text links that previously pointed at the not-yet-authored `explanation/*` pages are retargeted (2026-07-21) to the normative `spec/` artifacts on GitHub now that D1/D2 are dropped, so `warnonly` no longer masks a permanently-dangling internal link. E3 (refined case-study HTML for B1/B2) and E4 (the two paper builds) tracked.

## 10. Progress tracker

Legend: ✅ done · 🟡 in progress · ⬜ not started.

| Facet | ID | Deliverable | Source | Status |
|---|---|---|---|---|
| Charter | — | This document | — | ✅ |
| Tutorials | A1 | Introductory (exemplar, rendered) | core_engine_tour | ✅ |
| Tutorials | A2 | Advanced | agentic_pipeline + core §3-4 | ✅ (closes on Series-B Δ +0.52 launches ± 0.24) |
| Tutorials | A3 | Expert | aa_integration + wires + refinement | ✅ (closes on the off-wire cash reconstruction, Jacobi lag) |
| Tutorials | A4 | Deep-dive: serialization | agentic_pipeline §4,6 | ✅ (byte-identical DSL-vs-reload trajectory) |
| Tutorials | A5 | Deep-dive: composition | refinement_tour + core §7 | ✅ (granularity-substitution Σct/ΠPoS match) |
| Case study | B1 | Marginal eNPV of the Nth scientist (flagship, HTML) | new + introspection_tour | ✅ page (shadow price +$19M ± $2.5M); refined HTML = E3 |
| Case study | B2 | In-licensing asset value (BD/M&A, HTML) | bd_acquisition | ✅ page (Δ-rNPV +1947 ± 367, not additive); refined HTML = E3 |
| Case study | B3 | When to kill a program | new + agentic_pipeline | ✅ (θ\*=0.5 worth +48.5 ± 12.6 vs never-kill) |
| Reference | C1 | Capability-organized API pages | src docstrings | ✅ (9 pages; every @docs symbol resolves) |
| Reference | C2 | JSON model schema page | CONTRACT §8 / ADR 0005 | ✅ |
| Explanation | ~~D1~~ | ~~Contract-as-explanation pages~~ | — | ❌ dropped 2026-07-21 (site links to normative `spec/` instead) |
| Explanation | ~~D2~~ | ~~ADR reading surface~~ | — | ❌ dropped 2026-07-21 (site links to `spec/adr/` instead) |
| Explanation | **D3a** | **Methodology paper (arXiv-first, technical; outline §8 fixed, format TBD)** | CONTRACT §1–§15 + ADRs + case studies | ⬜ |
| Explanation | **D3b** | **Adoption/value paper (HBR-style, executive; leads with the problem, not the use cases)** | motivation → solution → some examples (Workstream B) | ⬜ |
| Explanation | D-future | JOSS `paper.md` + methods paper (gates: CI scaffold ✅, release/DOI ⬜, CONTRIBUTING ✅) | fed by D3a | ⬜ (parked) |
| Build | E1 | make.jl (Documenter + Literate) | wires_viz_tour/build.jl | ✅ (full nav wired: all 8 tutorials + 9 reference pages) |
| Build | E2 | docs/Project.toml (site build green) | — | ✅ (demo-union deps; full site builds green; explanation links retargeted to `spec/`) |
| Build | E3 | Refined-HTML render recipe (case studies) | wires_viz_tour/build.jl | ⬜ (B1/B2 Literate pages done; refined-HTML presentation pending — companion to D3b) |
| Build | E4 | Paper builds (latexmk or pandoc), one per paper | — | ⬜ (with D3a/D3b formats) |

## 11. Sequencing

1. **Charter + exemplar + toolchain (this PR increment).** This document + the introductory tutorial (A1) authored in the target Literate style and rendered end-to-end, plus the `make.jl` + `docs/Project.toml` Documenter+Literate build (E1/E2), which builds the site green with A1 in place. *(done)*
2. **Reference (C1/C2).** Land the capability-organized reference so the site is navigable end-to-end. *(done. The former "explanation scaffolds (D1/D2)" step is dropped — the site links to the normative `spec/` artifacts instead of re-hosting them.)*
3. **Tutorial tiers (A2/A3) + deep-dives (A4/A5).** Migrate the remaining demos in. *(done — all four rendered green.)*
4. **Case studies (B1/B2/B3).** Flagship first; refined HTML for B1/B2 (E3). These become the papers' worked evidence (D3a §13; D3b's illustrative examples). *(Literate pages done — all three rendered green with their headline numbers; the E3 refined-HTML presentation for B1/B2 remains.)*
5. **The two papers (D3a/D3b).** Once the case studies exist as evidence, decide each paper's format (§8 D3), then author them against the fixed outlines — D3a mechanism-first (arXiv, technical), D3b problem-first (HBR-style, executive) — and wire their builds (E4). Post the methodology preprint for socialization. *(D3a's semantics sections may run in parallel with 3–4, as they do not depend on the case studies; D3b's examples do.)*
6. **Polish + deploy.** Runic pass, link audit, GitHub Pages deploy, README refresh to point at the published site + the papers. Optionally open the remaining D-future gates (tagged release + Zenodo DOI) toward a JOSS submission (CI scaffold and `CONTRIBUTING` now in place).

Each numbered step is a reviewable batch; the charter's tracker (§10) is updated as facets land.

## 12. Acceptance criteria (per page)

A tutorial or case-study page is *done* when: every code block executes clean during the site build; it ends in (tutorial) or leads with (case study) a computed, decision-relevant quantity, framed in a computational register (not as an "answer for a manager"); the page's constructs each trace to a passing semantic test or a docstring — **verified in review/CI, and NOT stated to the reader** (no internal-QA meta in reader prose: no "covered by the test suite", "invents no API", "this API is stable"); the run is `seed=`-pinned and reproducible; the source passes Runic; and — for case studies — the title is the question, and the mechanics read as "how we got it," not as the subject. A reference page is done when every documented symbol resolves against the current export surface (the `exports_resolve.jl` invariant) and no removed symbol is referenced.

## 13. Decisions log

- **2026-07-17 — Migrate demos in place** (single source of truth), not keep-alongside or absorb-and-retire. Rationale: no duplicated model code to drift; demos stay independently runnable.
- **2026-07-17 — This-PR scope: charter + one rendered exemplar** (the introductory tutorial), not full scaffolding. The exemplar sets the quality bar; the rest is tracked in §10 and sequenced in §11.
- **2026-07-17 — Introductory tier is a strict subset** (plain-species timed pipeline only); the modality truth-table and the allocator move to the advanced tier, per the PR-draft tier mapping.
- **2026-07-17 — Voice & hierarchy calibration (maintainer).** (1) Drop the "answer for a manager" register — it reads as role-flattering ("little-guy aspirational"); every tutorial instead *computes a decision-relevant quantity* framed as information a decision rests on, in a computational-framework register. (2) The "no invented API / traces to the test suite" rule is an **internal** acceptance gate verified in the background, **never** surfaced in reader prose. (3) Fix the tutorial hierarchy: the model sections are the spine, seeding/ensembles is a *supporting* subsection (not a peer of the first/second model), the decision-quantity is the payoff. (4) Calibrate tutorials to typical Julia computational-framework conventions (§4 build stack + the conventions notes), keeping scope as a *computational framework* tutorial. Applied to the A1 exemplar; binding on all tiers.
- **2026-07-17 — Workstream D is anchored by an academic paper, arXiv-first.** *(Partly superseded 2026-07-21 — see the next entry: the anchor is now two peer papers, and the D1/D2 web explanation layer is dropped.)* The explanation layer is framed as a scholarly paper (semantics + why + backing concepts) posted to arXiv for socialization and citation, downstream feeding a thin JOSS `paper.md` and/or a peer-reviewed methods paper. Rationale: arXiv has no infra prerequisites (unlike JOSS, which needs CI + a release DOI first), suits long-form, and doesn't foreclose the deeper venues. The §8 outline is fixed; the paper's format/home (LaTeX+arXiv vs Markdown+pandoc) is **deferred to authoring time**. The normative contract stays in `spec/CONTRACT_DRAFT.md`; the paper cites it, never forks it (§8.4 anti-fork rule).
- **2026-07-21 — Drop the web explanation quadrant (D1/D2); split the paper in two (maintainer).** (1) **Drop D1 (contract-as-explanation pages) and D2 (ADR reading surface).** Rebuilding the contract and ADRs as Documenter pages is a third copy that drifts from the normative `spec/` source for no reader payoff the papers don't cover; the site links directly to `spec/CONTRACT_DRAFT.md` and `spec/adr/` on GitHub instead. In-text `explanation/*` links across the tutorials/deep-dives/case-studies/reference are retargeted accordingly. (2) **Split the single arXiv paper into two peer papers.** D3a — a technical/computational-concepts methodology paper (mechanism-first, arXiv). D3b — an HBR-style adoption/value paper for technically non-savvy executives that does *not* lead with the use cases but builds motivation → solution → *some* examples, to explain "what we bring." Same substance, opposite lead-ins and audiences; both non-normative (§8.4 unchanged). Format/home for each still deferred to authoring time.
