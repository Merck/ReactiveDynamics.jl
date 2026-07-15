# CLAUDE.md — ReactiveDynamics.jl agent guide

Agent-facing orientation for the `ref-agents` branch. Re-verify any `file:line` before acting — line numbers drift; grep to confirm.

## What this project is

ReactiveDynamics.jl (RD) is a timed, stochastic, resource-constrained Petri net / discrete-event engine for system-dynamics-style modeling of business / R&D processes (budgeting, ledgers, what-if, rNPV) — NOT a chemical reaction network, despite the Catalyst-derived DSL surface. A **transition** is a stateful recipe that spawns in-flight instances at a Poisson rate, occupies shared finite resources (**species**) over a cycle time, and completes with a terminal probability-of-success that emits RHS products. Resources carry a **modality** governing allocation (`:nonblock`/`:conserved`/`:rate`), and a cost/reward/valuation **ledger** accumulates into a per-step log. Structured/agentic **tokens** are first-class entities (a "project" is an agent with attributes, custom behavior, and history) that can be instantiated, selected by predicate, advanced through lifecycle phases, and audited per-program. RD is part of the DyVE family and sits ON TOP of [AlgebraicAgents.jl](https://github.com/Merck/AlgebraicAgents.jl) (AA): a `ReactionNetworkProblem` IS an AA `@aagent`, so a network is a node in a larger heterogeneous AA hierarchy.

## Where the durable spec lives

- **`docs/STATUS.md`** — the single "what is the state, what is left" index. Start here.
- **`docs/CONTRACT_DRAFT.md`** — the normative operational-semantics spec (§1–§15). The durable spine.
- **`docs/adr/`** — Architecture Decision Records 0001–0015 (append-only), with the status table in `docs/adr/README.md`. To change a decision, add a new ADR; do not rewrite history.
- **`INVENTORY.md`** — the current-source map (module map, public-API audit, static store, stepping trace, AA touchpoints).
- **`readme.md`** — the user-facing package README (about + four worked sketches).
- `docs/PHASE0_REVIEW.md` and `docs/MVP_BD_DEMO.md` are records/design docs (Phase-0 sign-off; the BD acquisition demo that drove §12). Superseded planning docs are archived under `docs/history/`.

## Current architecture (post-ADR-0015)

The static authoring/IR store is a dependency-free typed struct-of-columns (ACSets and Catlab were dropped, ADR 0003). ADR 0015 retired the last ACSets-lineage vocabulary, so the CURRENT names are:

- Store type **`ReactionNetwork`** (was `ReactionNetworkSchema`); the type-level object model is `const SCHEMA`.
- Authoring macro **`@reaction_network`** (was `@ReactionNetworkSchema`).
- Conventional variable **`net`** (was `acs`); the `ReactionNetworkProblem` field **`.network`** (was `.acs`); store fields **`.counts`**/**`.columns`** (were `.parts`/`.subparts`).
- Store verbs renamed to store vocabulary AND UNEXPORTED (internal, `RD.`-qualified): `nrows`, `row_ids`, `col_row_ids`, `column`/`cell`, `set_cell!`, `add_row!`/`add_rows!`, `rem_rows!`, `find_rows` (a `findall`, NOT an FK follow), `merge_networks!`, `build_network_from_dict`.
- Old names survive ONE release as `@deprecate`/`@deprecate_binding` shims (`src/ReactiveDynamics.jl:427-435`; the deprecated `@ReactionNetworkSchema` macro alias at `src/interface/create.jl:83`). GeneratedExpressions was dropped.

Module layout (`src/`, include order orchestrated in `src/ReactiveDynamics.jl`): `state.jl` (live sim state, the `@aagent` structs), `exprnode.jl` (the closed eval-free `ExprNode` IR), `compilers.jl` (expr→closure), `interface/*` (the DSL: `create.jl`/`update.jl` macros, `agents.jl` structured tokens, `aa_coupling.jl` the AA read/coupling surface, `checkpoint.jl` dump/restore, `solve.jl` `@agentize`), `operators/*` (`joins.jl`/`equalize.jl`/`refine.jl`), `solvers.jl` (THE step loop + constructor + `_step!`/`_reinit!`/`_projected_to`), `predicates.jl` (`TokenPredicate`/`@select`), `actions.jl` (the closed action family + `Rule`), `ledger.jl` (per-program ledger), `serialize.jl` (single-JSON eval-free (de)serialization), `analysis.jl`/`export.jl`/`visualize.jl` (Phase-0.6 result-inspection layer). Extensions: `ext/RDPlotsExt.jl` (Plots recipes + `_draw`, weakdep), `ext/RDArrowExt.jl` (Arrow export, weakdep).

Everything specified in the contract (§1–§15) and ADRs 0001–0015 is IMPLEMENTED and green. See `docs/STATUS.md` for the feature→commit table and the small remaining-work list (3 modality-validation `@test_skip` pins; `dump_state`'s clean-tick-boundary constraint; deferred entity-level refinement, threaded ensemble backend, and Opera implicit coupling).

## Build / test / dev

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install deps
julia --project=. -e 'using Pkg; Pkg.test()'          # run the semantic suite
julia --project=. -e 'using ReactiveDynamics'         # load the package
```

Tests live in `test/semantic/*.jl` (entry `test/runtests.jl` → `test/semantic/runtests.jl`) under a two-tier scheme: **T1-characterization** (runs against the engine) and **T2-acceptance** (target behavior). The suite is green; the only non-passing tests are deliberate `@test_skip`s — 3 CONTRACT §1.4 modality-validation placeholders (`test/semantic/modality_genesis.jl:172,193,214`) plus a few weakdep environment guards (Plots/Arrow/Graphviz not loaded). There are NO live `@test_broken` pins. `test/Project.toml` declares the test-only deps (Plots, Arrow, DataFrames, Distributions, …) that the main project keeps as weakdeps. AlgebraicAgents is the registered 0.4 release (no `[sources]` pin). Julia ≥ 1.11.

Run `julia -e 'using JuliaFormatter; format(".")'` before committing (config in `.JuliaFormatter.toml`). CI: `.github/workflows/` (Tests, Documenter, Formatter, CompatHelper, TagBot).

## Dev-loop gotchas

- **Determinism is contractual.** All stochasticity draws from a state-owned `rng`; a run is determined by `(model, seed)` (CONTRACT §4). Thread the seed; never reach for the global RNG. `_reinit!` restores the exact stream (and takes a `seed` kwarg for ensemble mode (b) reseed).
- **Append-only store.** Compiled attribute closures hard-code each species/param position; runtime mutation is append-only + soft-deactivate (ADR 0004). Never reorder or delete rows on a live/stepping model; reindexing macros (`equalize!`, refinement) are authoring-time only.
- **Eval-free everywhere.** The serializer and IR never `Meta.parse`/`eval` a model field — the closed `ExprNode`/action whitelist + `validate` is the trust boundary (ADR 0005/0006). Custom host functions go through the per-network registry, not `@eval`. Do not add runtime `eval`.
- **Adding a new `src/` file** requires clearing the compiled cache; the `dev/` Revise test server (`dev/run.sh`) otherwise avoids recompilation across edits.
- **`NodeRef`, not `Ref`**, in the IR; JSON alloc-strategy fields are Symbols. Structured-token iteration uses a `(species, creation_index)` total order for determinism.
- **Worktrees** go under `~/worktrees/<repo-short>/<branch>` (never repo-siblings or nested inside the tree).

## Working rules

Source of truth for decisions is `docs/adr/`; the normative spec is `docs/CONTRACT_DRAFT.md`. When docs contradict code, the CODE is truth. Commit at the end of a coherent logical batch (not per-file). Do not hard-wrap prose in Markdown — one continuous line per paragraph/list-item.
