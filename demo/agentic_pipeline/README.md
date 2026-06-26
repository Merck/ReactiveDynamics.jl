# Agentic R&D-Pipeline Demo

A focused, didactic tour of ReactiveDynamics.jl's **production & agentic** capabilities on one small, self-contained scenario: a tiny R&D project portfolio advanced through a lifecycle pipeline, steered by an in-model decision rule, then serialized to an eval-free JSON model, validated, checkpointed, and replayed. It is the *light* companion to the heavier [`demo/bd_acquisition/`](../bd_acquisition/) case study — same engine machinery, much simpler numbers, narrated step by step.

## Run it

```bash
julia --project=. demo/agentic_pipeline/agentic_pipeline.jl
```

One command runs the whole tour end to end. The structured-token kind and its registry are defined inside the script, so it is fully self-contained; every run uses an explicit `seed=` and is reproducible.

## The scenario

Each project is a **structured token** of kind `:Project` carrying a lifecycle `phase` (Phase1 → Phase2 → Phase3 → Launched) and a net-present value `npv`. Projects advance through the pipeline via predicate-selected transitions; risky stage gates fail a fraction of attempts (the program soft-retires). A management lever — a Series-B raise that injects capital, flips a synergy parameter, and adds a fresh program — lives in the model as a typed `Rule`. The same pipeline is then written as JSON, validated, loaded (matching the DSL build bit-for-bit), checkpointed mid-run, and replayed.

## What each section shows

| § | Section | Engine capability | ADR |
|---|---|---|---|
| 0 | The structured-token KIND | A project is a first-class entity carrying attributes (`phase`, `npv`) and a stable identity — not an anonymous scalar count; defined in RD scope via the `@register`/`@aagent` idiom, resolved by name through a registry | [0006](../../docs/adr/0006-structured-tokens.md) |
| 1 | Phase-as-attribute pipeline + declarative portfolio | One `:Project` kind with a `phase` **attribute** (not a species per phase); steps are `@select(Project, …) --> @advance(phase, :Next)` preserving identity; a `probability` gate soft-retires failures; the starting portfolio is the declarative `population[]` initial marking (both authoring forms) so the run is reproducible input | [0008](../../docs/adr/0008-token-filtration.md), [0007](../../docs/adr/0007-interface-and-initial-state.md) |
| 2 | Predicate selection | `@select` with a continuous clause (`npv > θ`) binds only the qualifying subset; equal-priority ties break deterministically by `creation_index` (first added wins), so selection reproduces under `(model, seed)` | [0008](../../docs/adr/0008-token-filtration.md) |
| 3 | In-model decision rule | A `fire_mode=:once` `Rule` (the management lever) whose `Seq` action injects capital (`SetSpecies`), flips a parameter (`SetParams`), and adds a program (`AddToken`); `set_guard!` withholds a line until funded. The decision lives in the model, not in host patch code | [0010](../../docs/adr/0010-rules-and-conditional-transitions.md) |
| 4 | Population write | `SetTokens` with `@field` revalues a `@select`-ed sub-population (write down every Phase-2 valuation 10%), reading each token's own current field; `@field` is legal only in a write value, never in a predicate | [0011](../../docs/adr/0011-action-callbacks-and-general-code.md) |
| 5 | Model-as-data (eval-free JSON) | The same pipeline as a JSON document: `validate` (clean pass + a deliberately broken model yielding a diagnostic), `from_json_model` matching the DSL build bit-for-bit under the same seed, the security point (a malicious string param loads as inert data, never executed), and loading the document from a file with `@import_model`. (The canonical direction is authored-document → load, exactly how `demo/bd_acquisition/model.rdj.json` is consumed; emitting a *live* model back to a full JSON document is not yet complete, so the run's reproducible artifacts are the authored input plus the solution trajectory via `@export_solution_as_table`.) | [0005](../../docs/adr/0005-serialization-json-ir.md) |
| 6 | Checkpoint & replay | `dump_state` at a clean tick boundary (model uses `cycletime => 0.0` so nothing is in-flight) → `restore` → resume reproduces the continuation; `_reinit!` resets state, rebuilds the t=0 population, and re-arms once-rules so a re-run reproduces the first trajectory | [0007](../../docs/adr/0007-interface-and-initial-state.md) |
| 7 | Recap | The through-line: a model — species, pipeline, levers, and starting portfolio — is reproducible DATA, fully determined by `(model, population, rules, seed)` | — |

## Why these choices

- **Structured tokens over plain species** — a plain Float64 species is an anonymous count; a structured token carries attributes and a stable identity, so a project can be *selected by value* and *followed* through every phase to launch (the same object, mutated in place by `@advance`).
- **Phase-as-attribute over a species-per-phase** — one `:Project` kind with a `phase` field keeps identity intact across advances and avoids an explosion of species; a pipeline step is a predicate-filtered `@select`/`@advance` pair.
- **The decision rule lives in the model** — a management lever is a typed `Rule`, part of the model and the `(model, rules, seed)` triple, not external host code reaching in mid-run; this is what makes a scenario reproducible rather than an imperative script.
- **Eval-free JSON matters for security and agentic authoring** — `from_json_model` never `eval`s or `Meta.parse`s a model field, so an untrusted party (or an LLM) can author a model as JSON and the worst a malicious field can do is fail validation; host code is referenced by name through a registry, never carried in the document.

## Relationship to the test suite

Every construct in the script is copied from the engine's passing semantic tests: [`test/semantic/token_filtration.jl`](../../test/semantic/token_filtration.jl) (structured tokens, `@select`/`@advance`), [`test/semantic/rules_decisions.jl`](../../test/semantic/rules_decisions.jl) (Rules, the action family, `SetTokens`), [`test/semantic/initial_state.jl`](../../test/semantic/initial_state.jl) (`population[]`, `reinit`, `dump_state`/`restore`), and [`test/semantic/serialization_ir.jl`](../../test/semantic/serialization_ir.jl) (JSON model, `validate`, `@import_model`).
