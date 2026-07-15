# ADR 0015: Post-ACSets naming — retire the categorical vocabulary and the GeneratedExpressions dependency

Status: Proposed — 2026-07-15

Date: 2026-07-15

Supersedes/relates to: [ADR 0003](0003-data-store.md) (dropped ACSets for a dep-free typed store, but deliberately KEPT the ACSets-lineage names for signature compatibility — this ADR pays down that named debt). Relates to [ADR 0005](0005-serialization-json-ir.md) (the JSON IR key strings are a SEPARATE surface, out of scope here — see §Open questions). Does not touch the runtime engine semantics ([ADR 0001](0001-discrete-event-engine.md)) or the allocator ([ADR 0002](0002-priority-weighted-allocation.md)); this is a pure naming + dependency change.

## Context

ADR 0003 removed ACSets.jl and Catlab entirely: the static authoring/IR store is now a dependency-free typed struct-of-columns (`src/ReactiveDynamics.jl:39-330`), with `const SCHEMA` (`ReactiveDynamics.jl:47`) as the type-level object model and zero homs/morphisms/colimits anywhere in `src/` (ADR 0003 §Rationale: "its entire categorical surface … is unused here"). But ADR 0003 explicitly kept the old vocabulary — "The store type keeps the name `ReactionNetworkSchema` … so every existing signature … compiles unchanged" (`ReactiveDynamics.jl:141`). The retention reason on record is compile-time compatibility, NOT semantic aptness.

The result is a public API whose central names now misdescribe the implementation on two axes:

1. **Instance-vs-schema inversion.** `@ReactionNetworkSchema` produces a `ReactionNetworkSchema` struct that is a *populated network instance* holding data columns (`parts::Dict`, `subparts::NamedTuple`, `reactants::Vector{ReactantSpec}`). That is a model value, not a schema. The thing that genuinely IS the schema — the type-level object model — already exists separately as `const SCHEMA`. Naming the instance `…Schema` is like naming a `DataFrame` value `TableSchema`.

2. **Dead categorical vocabulary.** The name is a fossil of `@present TheoryReactionNetwork(FreeSchema)` → `const ReactionNetworkSchema = @acset_type(...)`. The exported store verbs `nparts, parts, dom_parts, subpart, set_subpart!, add_part!, add_parts!, rem_parts!, incident` are re-implemented ACSets verbs (`ReactiveDynamics.jl:15`). Two actively mislead: `incident` is just `findall` (it implies a foreign-key/hom follow that does not exist), and `subpart`/`part` are category-theoretic terms for a store with no morphisms. The variable name `acs` (≈640 sites across `src/` + `test/`) abbreviates "ACSet" — a data structure the code no longer contains. Docstrings and comments still say "build … an acset" (20 `acset` mentions in `src/`).

Separately, the package still depends on **GeneratedExpressions.jl** (`Project.toml:12`, reexported at `ReactiveDynamics.jl:7`) purely for its `generate(...)` brace-comprehension expander, called at exactly two sites (`create.jl:73,79`). No current model, test, or demo uses the `{$a*r[$l], l=1:$r, dlm=+}` comprehension syntax; `generate` runs as an identity pass in every real invocation. The feature is advertised only in `readme.md:185-239`, whose examples are AlgebraicJulia-era (they still call the model an "attributed C-set") and are effectively dead documentation. The dependency also complicates the export-resolution invariant: `test/semantic/exports_resolve.jl:27-28` must special-case `getproperty(RD, :GeneratedExpressions)` as an unpoliced reexport surface.

Two stale-only lineage tokens exist in `src/` and should die in the same sweep: `TheoryReactionNetwork` and `DiscreteProblem`, both in the dead docstring at `create.jl:40` (neither type exists anymore).

This is the moment to pay the debt: the store is stable, the ADR-0003 migration is complete, and the names are the last thing still pointing at the abandoned lineage.

## Decision

Rename the ACSets-lineage public and internal surface to store/network vocabulary, remove the GeneratedExpressions dependency, and correct the stale documentation — in one coordinated, mechanical pass, staged behind deprecation shims where the surface is public.

### Tier 1 — the network type, macro, variable, and field (headline)

| Now | New | Notes |
|---|---|---|
| `struct ReactionNetworkSchema` | `struct ReactionNetwork` | The instance type. `ReactionNetwork` is free (the old `const ReactionNetwork` alias was renamed away pre-ACSets). Now reads correctly against `const SCHEMA`. |
| `@ReactionNetworkSchema` | `@reaction_network` | Lowercase, matching Julia macro idiom and the Catalyst.jl precedent (`@reaction_network` / `ReactionSystem`). |
| variable `acs` (~640 sites) | `net` | Reads as `net::ReactionNetwork`. Pure word-boundary find-replace across `src/`, `test/`, `demo/`. |
| field `.acs` on `ReactionNetworkProblem` (34 sites) | `.network` | `state.acs::ReactionNetworkSchema` → `state.network::ReactionNetwork`. Field of an `@aagent` struct (`state.jl:64`). |
| `SampleableValues`-doc "an acset" / "the acset" (comments, docstrings) | "a network" / "the network" | 20 `acset` mentions in `src/` comments + docstrings. |

Compatibility: add `Base.@deprecate_binding ReactionNetworkSchema ReactionNetwork` and keep a thin deprecated `@ReactionNetworkSchema` that expands to `@reaction_network` (with a `Base.depwarn`). Old caller code keeps working for one release with a warning, then the shims are removed.

### Tier 2 — the store API shim (unexport AND rename now)

These are an internal store implementation detail that leaked into the export list only because pre-migration callers used the ACSets verbs bare. Rename to store vocabulary AND remove from the export surface (they become `RD.`-qualified internals). The `const SCHEMA` name stays (it genuinely is the schema).

| Now (exported) | New (internal, `RD.`-qualified) | Semantics |
|---|---|---|
| `nparts(net, obj)` | `nrows(net, obj)` | row count for an object |
| `parts(net, obj)` | `row_ids(net, obj)` | `Base.OneTo(count)` |
| `dom_parts(net, attr)` | `col_row_ids(net, attr)` | rows of the attr's owning object |
| `subpart(net, attr)` / `subpart(net, i, attr)` | `column(net, attr)` / `cell(net, i, attr)` | whole-column COPY / single cell |
| `set_subpart!(net, i, attr, v)` | `set_cell!(net, i, attr, v)` | write a cell |
| `add_part!(net, obj; kwargs...)` | `add_row!(net, obj; kwargs...)` | append a row |
| `add_parts!(net, obj, m)` | `add_rows!(net, obj, m)` | append m rows |
| `rem_parts!(net, obj, idxs)` | `rem_rows!(net, obj, idxs)` | swap-and-pop delete (semantics unchanged) |
| `incident(net, val, attr)` | `find_rows(net, val, attr)` | `findall` over a column (NOT an FK follow) |
| field `.parts` | `.counts` | `Dict{Symbol,Int}` row counts |
| field `.subparts` | `.columns` | `NamedTuple` of `AttrColumn`s |

The eight `propertynames(net.columns)` reflection loops (compilers/solvers/joins/equalize) are updated to the new field name; iteration order is preserved (the `ALLATTRS` order is unchanged). `columns(SCHEMA)` / `columns(SCHEMA, obj)` accessors keep their names (they already describe columns).

Two more `acs`-in-the-name exports move to network vocabulary:

| Now (exported) | New |
|---|---|
| `union_acs!` (21 sites) | `merge_networks!` |
| `build_acs_from_dict` (13 sites) | `build_network_from_dict` |
| internal `merge_acs!` (5 sites) | `merge_network!` |

Because Tier 2 removes names from the export surface, the deprecation is a re-export of the OLD name bound to the new function plus a `depwarn`, held for one release. `union_acs!` and `build_acs_from_dict` were exported and get `@deprecate` (function-level, which warns AND forwards).

### Tier 3 — dependency removal, consistency, dead text

- **Remove GeneratedExpressions.** Delete the dep (`Project.toml:12`, `[compat]:60`), drop `@reexport using GeneratedExpressions` (`ReactiveDynamics.jl:7`), and inline the two `create.jl` call sites (`:73,79`) to skip the comprehension pass — `make_ReactionNetwork` operates on the block expression directly. Delete the readme comprehension section (`readme.md:35, 49, 185-239`) and the `@generate`/`@fileval` mentions. Update `exports_resolve.jl:27-28` to drop the `GeneratedExpressions` reexport-module reference (AlgebraicAgents remains). The `{…, a=1:n, dlm=+}` brace-comprehension feature is RETIRED — it is used by no model/test/demo and its only documentation was stale.
- **Standardize the `ReactionNetworkProblem` variable.** It is `state` in the engine core (`state.jl`, `solvers.jl`) and `prob` in the viz/serialize/export/plot layers. Standardize the PUBLIC/boundary API on `prob` (matches the `ReactionNetworkProblem` constructor and SciML convention); keep `state` only inside the hot `_step!`/`_reinit!` loop where it is entrenched and local. (Lower priority; do not let it balloon the diff.)
- **Delete the dead docstring** at `create.jl:40` ("outputs an instance of `TheoryReactionNetwork` that can be converted to a `DiscreteProblem`") — both types are gone. Replace with an accurate one-line description referencing `ReactionNetwork`.
- **Correct stale docs** that still describe ACSets/`BasicSchema`/`TheoryReactionNetwork`/`@present FreeSchema` as PRESENT: `INVENTORY.md:5,230-233`, `docs/HANDOFF_PLAN.md:40`, `docs/adr/0014-visualization.md:4`, `REVIEW.md:26,119`, and the `readme.md:31` "attributed C-set" line. These describe an intermediate or pre-migration state and now contradict `src/`.
- **`dt` / `tstep`.** `CONTRACT_DRAFT.md:108` already flags `tstep` as an internal alias and rename candidate. Out of scope for THIS ADR (it is not lineage debt), but noted so a future pass folds it in; not blocking.

## Rationale

- **Names should describe the implementation, not its history.** The store has zero categorical structure; `schema`/`subpart`/`incident`/`part`/`acs` all point at a paradigm that was removed. `ReactionNetwork` + `column`/`cell`/`row`/`find_rows` describe exactly what the struct-of-columns is.
- **Resolve the `SCHEMA` collision.** With the instance renamed to `ReactionNetwork`, `const SCHEMA` (the type-level model) and `ReactionNetwork` (an instance built from it) finally read correctly against each other.
- **Ecosystem legibility.** `@reaction_network` / `ReactionNetwork` matches Catalyst.jl, the reference Julia reaction-network package, so the API is immediately legible to that audience.
- **Shrink the public surface.** The Tier-2 verbs never should have been public; unexporting them narrows the API RD must keep stable and removes nine misleading names from tab-completion.
- **One fewer dependency, one less stale feature.** GeneratedExpressions buys a comprehension metalanguage nothing uses; removing it deletes a dep, a reexport special-case in the export test, and a block of dead readme documentation.
- **Non-breaking rollout.** Deprecation shims mean external code keeps working for a release; the change is mechanical and verifiable by the existing suite (500 pass / 7 broken).

## Consequences

- Large but mechanical diff: type ≈89 sites, macro ≈141 sites, `acs` variable ≈640 sites, `.acs` field 34 sites, Tier-2 verbs across `src/` + `test/` + `demo/`. Deprecation shims localize risk.
- The JSON IR (ADR 0005) key strings and the `const SCHEMA` object symbols (`:S,:T,:E,:obs,:P,:M`, the `spec*`/`trans*` column names) are UNCHANGED — serialized models still round-trip. Only Julia-level identifiers move.
- One release carries deprecation warnings for `ReactionNetworkSchema`, `@ReactionNetworkSchema`, `union_acs!`, `build_acs_from_dict`; then a follow-up removes the shims.
- The brace-comprehension feature is gone. If a future model needs template expansion, reintroduce it as an owned, tested internal utility rather than a dependency (recorded in Open questions).
- Docs across `INVENTORY.md`, `HANDOFF_PLAN.md`, `REVIEW.md`, `readme.md`, ADR 0014 stop contradicting the source.

## Migration / handoff plan

Ordered so the suite stays green at each step:

1. **Tier 1 defs + shims.** Rename the struct and its constructors/accessors to `ReactionNetwork`; add `@deprecate_binding`. Rename the macro to `@reaction_network`; add the deprecated `@ReactionNetworkSchema` alias. Rename the `.acs` field to `.network` on `ReactionNetworkProblem` and all field readers.
2. **Tier 1 sweep.** Word-boundary replace `acs` → `net` across `src/`, `test/`, `demo/` (skip `docs/build/`), and the 20 `acset` comment/docstring mentions → "network".
3. **Tier 2.** Rename the store verbs + `.parts`/`.subparts` fields; update the eight reflection loops; unexport the verbs; add re-export deprecation shims for the previously-exported ones. Rename `union_acs!`/`build_acs_from_dict`/`merge_acs!`.
4. **Tier 3 dep removal.** Drop GeneratedExpressions from `Project.toml` + reexport; inline `create.jl:73,79`; update `exports_resolve.jl:27-28`; delete the readme comprehension section; delete the `create.jl:40` dead docstring.
5. **Tier 3 docs + consistency.** Fix the stale ACSets references listed above; standardize boundary-API `ReactionNetworkProblem` locals on `prob`.
6. **Verify.** `Pkg.test()` green (target 500 pass / 7 broken, unchanged); grep confirms zero remaining `ReactionNetworkSchema`/`acs`/`subpart`/`incident`/`GeneratedExpressions` tokens outside the deprecation shims and `docs/build/`; a smoke `using ReactiveDynamics; @reaction_network begin … end` builds and simulates.

Update `docs/adr/README.md` with the 0015 row and this file's summary.

## Open questions

- **Deprecation window.** One minor release, or remove the shims immediately given RD is pre-1.0 with no external users on record? (Recommendation: one release — cheap insurance, and it documents the rename in-code.)
- **`incident` → `find_rows` vs `find_incident`.** `find_rows` drops the misleading FK connotation entirely; if any future FK-following accessor is planned on `ReactantSpec`, reserve a distinct name for it (`follow_fk`?) so `find_rows` stays a pure column search.
- **Template expansion, if ever needed.** If model-template comprehensions return as a requirement, spec them as an owned internal utility (a ~40-line brace expander over the block AST) in a new ADR — not a dependency.
- **`dt`/`tstep`** rename (`CONTRACT_DRAFT.md:108`) — fold into a later naming pass; not part of this ADR.
