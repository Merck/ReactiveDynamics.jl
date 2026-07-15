# Introspection & exec-map tour (analysis / observability / visualization)

A single, runnable, literate walkthrough of what you do AFTER a run finishes: the READ-ONLY analysis, ensemble, export, and visualization layer of ReactiveDynamics (Phase-0.6 — ADR 0013 / 0014, CONTRACT §14 / §15). The [core_engine_tour](../core_engine_tour) demo shows how to BUILD and RUN models and [bd_acquisition](../bd_acquisition) shows a full worked application; this one is about reading a finished run — the trajectory log, the per-program ledger, ensembles, export bundles, and the headline network "exec map". None of it touches the §1–§9 dynamics: every construct is pure post-processing over a finished (or in-progress) `ReactionNetworkProblem`, lifted from the passing semantic tests (`test/semantic/analysis_observability.jl`, `visualization.jl`, `program_ledger.jl`) — the tour invents no API.

## Run it

```bash
julia --project=demo/introspection_tour demo/introspection_tour/introspection_tour.jl
```

First time (resolve the demo-local env — takes a few minutes):

```bash
julia --project=demo/introspection_tour -e 'using Pkg; Pkg.instantiate()'
```

The script is literate: every section opens with a block comment explaining the analysis idea, then runs the code, then `println`-narrates the result, so running it once tells the whole story top to bottom. It writes its artifacts (the rendered exec map, the recipe PNGs, the DOT sources, and the export bundles) under `demo/introspection_tour/output/` (git-ignored — regenerated each run).

### Why a demo-local `Project.toml`

The main ReactiveDynamics project is intentionally LEAN: `Plots` and `Arrow` are weakdeps (they back the `RDPlotsExt` / `RDArrowExt` extensions), so `julia --project=.` does not resolve them. This tour needs `Plots` (to render the exec map through AlgebraicAgents' Graphviz backend and for the §15.1 recipe types) and benefits from `Arrow` (the §14.3 export bundle then writes byte-faithful `.arrow` siblings). So this environment adds them on top of a path-dev'd ReactiveDynamics; loading them here also triggers the two extensions. `AlgebraicAgents` (0.4, the ADR 0012/0013/0014 integration target) is resolved as an ordinary registered dependency.

## The model in one paragraph

The whole tour runs on ONE small self-contained model — a two-phase "project advance". A `Project` is a structured token with a `phase` (a Symbol lifecycle attribute) and a numeric `value`, so the trajectory log records both a Symbol path and a Real path. One transition `adv` `@select`s a Phase1 project, meters `budget` at `@rate` over its cycletime, and `@advance`s it to Phase2 with `probability => 0.6` (a Binomial success draw, so the run is genuinely stochastic and the ensemble has real spread). `budget` starts deliberately SCARCE (8) and carries a `specCost`, so the pool burns to a trough of 0 — which is exactly what the exec map paints as starvation. The produced `:Project` carries a `specReward`, so a successful advance realizes reward. That single model makes the trajectory log, the exec map's starvation coloring, and both sides of the per-program ledger all non-trivial. It mirrors the tests' `traj_model` / `advance_cost_model` exactly.

## What each section exercises

| § | Section | Capability exercised |
|---|---|---|
| 0 | The model | A structured-token `@register`ed kind with `log_token_fields` opt-in; `@reaction_network` + index-assigned cost / reward / budget; `ReactionNetworkProblem(...; seed=, population=)` + `simulate` |
| 1 | Raw run artifacts | `prob.sol` (the marking DataFrame, read by column name), `prob.log` (the tagged event stream, reduced by tag to aggregate cost / reward), and `program_ledger(prob)` + `program_ledger_entries` (the per-program attribution DataFrame + append-only audit trail) |
| 2 | Per-token trajectory log | `token_trajectory(prob)` (long form: `t, program, species, <field>…`), one token's life by name, and predicate-scoped rows via a `@select` `TokenPredicate` / `Clause` (the same selection machinery the dynamics use) |
| 3 | "Typical" helpers | `representative_token(prob)` (the MEDOID program — closest to the cohort mean path) and `trajectory_envelope(prob)` (per-tick median + IQR band over each numeric logged field) |
| 4 | Ensemble analysis | `ensemble(build; nseed, root_seed)` (member `k` seeded `hash((root_seed, k))`), `summarize` (mean / sem / quantiles), `treatment_effect` (the unpaired A/B Δ with `se = √(var_b/n_b + var_d/n_d)`), and the `EnsembleProblem` as an AA-readable node (`observables` / `getobservable` / `inners`) |
| 5 | Export bundles | `export_run(prob, dir)` and `export_ensemble(ens, dir; metric)` — the CSV + JSON core (`trajectory` / `ledger` / `tokens` / `events.json` / `tokens.json` / `run.json`), plus the `.arrow` siblings (because this demo loads Arrow), with a manifest pinning `model_hash` + `seed` |
| 6 | Result-plot recipes | The model-agnostic wrapper types (`MarkingPlot`, `SaturationPlot`, `ValuationPlot`, `LedgerPlot`, `ThroughputPlot`, `TokenTrajectoryPlot`), two of them rendered to PNG via `RDPlotsExt` |
| 7 | ★ The exec map | The three separable layers: `network_graph` (pure Petri-net structure, does not perturb the RNG) → `to_graphviz` / `draw_network` (DOT, rendered via AA's `run_graphviz`) → `exec_map(prob; highlight)` (Layer A decorated with starvation coloring and a `@select` cohort's `past_bonds` token-path highlighting) — the maintainer's "system exec map readable for (in)efficiencies, decorated with results" |
| 8 | Recap | A closing summary of the whole surface |

## The exec map (§7) — the headline

The exec map is the maintainer's headline ask: a system diagram you can read for inefficiencies, decorated with the run's results. It is built in three layers, each usable alone:

- **Layer A — `network_graph(prob)`** returns a plain `NetworkGraph` (species/place nodes, transition nodes, arcs with stoichiometry + modality). It is a pure function of the model — no plotting dependency, no simulation, and it runs on a `deepcopy` so it does NOT perturb the caller's RNG.
- **Layer B — `to_graphviz(g)` / `draw_network(prob)`** emits Graphviz DOT (species as circles, transitions as boxes, arcs colored by resource modality) and renders it through AlgebraicAgents' `run_graphviz`. The DOT string is always obtainable even with no Graphviz backend; rendering is the only step that needs one.
- **Layer C — `exec_map(prob; highlight)`** decorates Layer A with finished-run statistics: species nodes filled gold where their pool ran to a trough (starvation), and — given a `highlight::TokenPredicate` — the matching cohort's `past_bonds` path through the net drawn as thickened arcs. It is read-only: it never mutates state or re-runs dynamics.

On this model the exec map fills `budget` gold (it starves to 0) and thickens the `Project → adv` arc (the advanced cohort's path) — the (in)efficiency read at a glance.

Rendering is best-effort. If a Graphviz backend (`Graphviz_jll` or a system `dot`) is present, `draw_network` / `exec_map` write an SVG; if not, the demo still emits the DOT source and prints a clear note, so you can render it later with `dot -Tsvg`. (On the authoring machine Graphviz 12.2.1 was present, so the SVGs render.)

## The load-bearing points (narrated in-line)

These are surfaced in the script so a reader does not trip over them: macro arguments are LITERAL (so `probability => 0.6` is a schema literal and cost / reward / budget are set by index assignment, not kwargs into the DSL — the tests' idiom); the trajectory-log opt-in is per-KIND via `log_token_fields`, so a kind that does not define it contributes no rows; `TokenPredicate` clauses compare against a `QuoteNode`-wrapped literal symbol; `ensemble` seeds members `hash((root_seed, k))` and always runs mode `:rebuild` (the reinit-reseed mode b is gated on ADR 0007 §D); `treatment_effect` is the UNPAIRED estimator (the two arms desync the shared-RNG-free streams); and the exec map's token-path highlighting maps a bond's transition INDEX to the same node id the graph uses (not the per-instance `"<name>_@<t>"`), so a highlighted arc actually matches a graph arc.
