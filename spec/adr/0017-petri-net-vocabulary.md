# ADR 0017: Should the engine retire "species" for standard Petri-net vocabulary?

Status: Accepted + Implemented — 2026-08-22 (branch `rename/petri-vocabulary`). Deciders: maintainer + package author. The rename landed as one commit per tier: Tier 2 store column symbols `spec*`→`place*` plus the two `occursin` reflection filters (`adaa0ed`); Tier 1 exported names `@add_species`→`@add_place`, `SetSpecies`→`SetMarking`, `ReactantSpec`→`ArcSpec`, `reactant_specs`→`arcs`, `specname`→`placename`, `register_structured_species!`→`register_token_kind!` (`ebf648c`); Tier 5 local identifiers, struct fields (`ArcSpec.species`→`.place`, `ReactionNetwork.reactants`→`.arcs`) and code prose (`fdbe572`, plus stragglers in `aeb7a61`); Tier 4 the spec, the docs site and the demo prose, with the term dictionary promoted into a CONTRACT **Glossary** and a `docs/src/glossary.md` page (`aeb7a61`); Tier 3 the serialized keys `"species"`→`"places"` and `"set_species"`→`"set_marking"` last and on its own. Every previously-exported name survives ONE release as a `@deprecate`/`@deprecate_binding` shim in the single block at `src/ReactiveDynamics.jl:488-524`, and the loader accepts the legacy JSON keys for ONE release with a deprecation warning — the open question about saved models outside the repository was resolved in the negative by the maintainer (2026-08-22), so there is no permanent alias and both shim sets are dropped together at the next minor release. `SCHEMA`'s `:S` object symbol and the `token`/`transition` words are unchanged. Suite fully green at the final commit — 821 pass / 0 fail / 0 broken / 0 skipped, against 801 at the branch base `c671875` (the +20 are new regression pins: the `_validate_node!` place-pool branch the rename exposed as uncovered, and the retired-key deprecation test).

Date: 2026-08-21

Supersedes/relates to: [ADR 0015](0015-post-acsets-naming.md) (retired the ACSets-lineage store names but deliberately scoped `species`/`reactant` and the JSON keys OUT — this ADR is the second half of that debt). Relates to [ADR 0003](0003-data-store.md) (dropped ACSets; its Phase-3 interop view was rejected, so no consumer of the AlgebraicPetri/Catalyst vocabulary remains), [ADR 0005](0005-serialization-json-ir.md) (the JSON key strings, in scope here), [ADR 0008](0008-token-filtration.md) (fixed *token* as the word for a structured entity), [ADR 0014](0014-visualization.md) (`MarkingPlot`, and the drawing layer already speaks of arcs). Pure naming change: no runtime semantics, no allocator behaviour, no wire-format semantics beyond key spelling.

## Problem

The engine's central noun is wrong twice over, and half the codebase already knows it. A resource pool is called a **species** — chemistry's word — in a tool whose own front page opens by denying that it is a chemical reaction network. The newer layers already use the standard words: the initial state is a *marking* (ADR 0007), the counts-over-time plot is `MarkingPlot`, the drawing code talks about *arcs*. The old core still says *species* and *reactant*: 448 occurrences in `src/`, 438 in `spec/`, 154 in `demo/`, 163 across the docs sources. A reader therefore meets two vocabularies for one concept, with no glossary reconciling them.

Doing nothing costs a translation tax on every tutorial, paper, and review conversation — "species means a resource pool, which is a place" — and it leaves the one word we do get right, *token*, reading to newcomers as the language-model sense.

## Options

**Option A — Do nothing.** Keep `species` and `reactant` everywhere; at most add a glossary line.

- **Gain:** zero effort, zero churn, no deprecation window; the demos and the drafted ADRs keep their current wording.
- **Cost:** the split vocabulary stays in the public documentation, and newcomers are taught chemistry words for a business-process engine.
- **Breaks when:** the two companion papers publish or an outside user adopts the API — after that the wording is expensive to change, because published examples pin it.

**Option B — Prose only.** Documentation, spec, and figures adopt place/marking/arc; the code keeps `species`.

- **Gain:** cheap (about a day), no shims, nothing can break for a caller.
- **Cost:** docs and code then disagree by design — the tutorial says *place* while the function is `@add_species` — and generated docstrings leak the old word back onto the site.
- **Breaks when:** a reader moves from the tutorial to the code, which is the normal path.

**Option C — Rename the whole surface, one release of shims.** Identifiers, store column names, and serialized keys all move to standard vocabulary; old names warn for one release.

- **Gain:** one vocabulary across code, spec, docs, and saved models; a saved model becomes readable by anyone who knows Petri nets.
- **Cost:** about 1 400 text sites, 7 exported names, 15 serialized-key sites, plus a permanent read alias so old saved models still parse.
- **Breaks when:** saved models exist outside this repository that nobody can regenerate — then the key rename buys nothing and the alias is pure overhead.

**Option D — Rename the code, freeze the file format.** As C, but saved-model keys and the internal column symbols keep the old spelling.

- **Gain:** no compatibility question at all; saved models round-trip untouched. This is the line ADR 0015 drew, for that reason.
- **Cost:** the exported artifact keeps the chemistry word, so anyone reading or hand-writing a saved model meets it again — and the field name no longer matches the function that fills it.
- **Breaks when:** the saved-model format becomes the primary way models are exchanged (its stated purpose: a machine-authored artifact).

## Recommendation

Take **Option C**, in two stages: prose and identifiers first, serialized keys second. The price is one release of the loader reading both spellings and warning on the old — the identifier shim discipline, dropped in the same release — for one vocabulary across code, spec, and format. Keep **token** — it is the canonical Petri word, no synonym is better, and one gloss at first use settles the language-model confusion. Confidence: moderate — the contested half is the file format, not the rename. Reverse to Option D on finding one saved model outside this repository that cannot be re-emitted.

## What changes

A model author writes `@add_place` instead of `@add_species`, and reads `SetMarking` where an action sets a pool count. Attribute names in authored models change spelling — `specName` becomes `placeName` — so any model that names an attribute directly needs a one-word edit. Old spellings work for one release with a warning, exactly as in ADR 0015.

Documentation keeps two registers deliberately: the spec and the API say *place*, *marking*, *arc*; the tutorials say *resource pool*, glossed once as "a place, in Petri-net terms". A tutorial reader needs no Petri theory; a spec reader meets no private dialect.

One-time work: about a day of renaming plus half a day re-reading the docs, one person in a single pass, verified by `julia --project=. -e 'using Pkg; Pkg.test()'` staying green and a grep gate confirming the old spellings survive only inside the shims. CI is unaffected.

## Open questions

- Does any saved model exist outside this repository — the single fact that decides C versus D? (Maintainer, before stage 2.) **RESOLVED 2026-08-22 (maintainer):** no saved models exist outside the repository, so the key rename lands with the rest of the change behind one-release shims.
- Should the shims be dropped at the next minor release, in step with the ADR 0015 shims, or kept a cycle longer? (Package author, at release time.)
- Is *place* or *pool* the better spec-level word, given that the domain audience reads "place" as a location? (Recommendation: `place` in the spec and API, `pool` in tutorial prose — resolve in review.)
- Does the rename wait for the two companion papers to fix their terminology first, so all three land consistent? (Maintainer.)
- Do the structured-token registry names move too (`register_structured_species!` → a token-kind verb), or is that a separate coloured-net naming pass? (Package author.) **RESOLVED in implementation:** it moved, to `register_token_kind!` — a registered kind is a *colour set*, not a place, so leaving it spelled `species` would have been the one remaining CRN word on the exported surface.

Technical detail, evidence, the full term dictionary, and the ordered migration: see Appendix.

## Appendix

### Evidence for Problem

Occurrences of the word `species` (case-insensitive, word boundary), re-measured at implementation time on the branch base `c671875` (`grep -rhoiE '\bspecies\b' <tree> | wc -l`):

| tree | count | as drafted (2026-08-21, on `rework`) |
|---|---|---|
| `src/` | 448 | 448 |
| `spec/` | 411 | 438 |
| `test/` | 177 | 177 |
| `demo/` | 154 | 154 |
| `docs/src/` | 31 | 97 |
| `docs/literate/` | 66 | 66 |

Two draft figures were wrong and are corrected above. `docs/src/` was never 97 — that is the `docs/src/` + `docs/literate/` total (31 + 66), double-counted. `spec/` reads 411 on the branch base rather than 438 because ADR 0016 (21 hits) and the presentation deck are not on `main` yet. The same word appears 131 times as `reactant`/`reactants` in `src/`, 90 in `spec/`, 53 in `test/`, 21 in `demo/`, 23 across `docs/`.

The split is not hypothetical — the two vocabularies are already interleaved in shipped code:

| standard term already in use | where | old term in the same engine |
|---|---|---|
| *marking* | `src/interface/checkpoint.jl` (8), `src/state.jl` (4), `src/solvers.jl` (3), `src/interface/agents.jl` (6), CONTRACT §7/§10 (12), ADR 0007 (27) | `specInitVal` |
| `MarkingPlot` (exported) | `src/analysis.jl:482,490` | its own docstring reads "species/token COUNTS" |
| *arc* | `src/visualize.jl` (22) | `ReactantSpec`, `reactant_specs` |
| *place* | — | `:S` object, `specName` |

Lineage: the vocabulary is inherited from `@present TheoryReactionNetwork(FreeSchema)` and the AlgebraicPetri/Catalyst surface. ADR 0003 removed the dependency; ADR 0003 Phase 3 (an ACSets interop view) was rejected outright. Nothing downstream now requires the names to match Catalyst, which is the only reason they were kept in ADR 0015.

The category error is worth stating precisely: `specInitVal` is a component of the **marking** — a property of the *net's state*, not of a substance. Chemistry conflates the two because a chemical species has intrinsic properties (mass, charge) that a place does not. Under RD's actual ontology, `specCost`/`specReward`/`specValuation` are per-token ledger rates attached to a place, and `specModality` is an arc discipline, not a property of a substance at all.

### Term dictionary (the durable artifact of this ADR)

The mapping from RD's concepts to published Petri-net vocabulary, with the extension family each borrows from. Adopting the right-hand column is what Option C means concretely; it also tells a reader which literature answers a question about the engine.

| RD concept | standard term | family / note |
|---|---|---|
| fungible species | **place**; its count is that place's **marking** | classical P/T net. `M₀` is the initial marking |
| structured species | place with a **colour set**; its tokens are distinguishable | Coloured Petri nets (Jensen). Token attributes are colours |
| reactant / product entry | input **arc** / output arc; stoichiometry is the **arc weight** (inscription) | classical |
| LHS / RHS multiset | **preset** `•t` / **postset** `t•` | classical |
| transition | **transition** | already canonical — do not rename |
| rate expression | **firing rate**, marking-dependent | Stochastic PN / GSPN. RD's Poisson intensity with unbounded concurrency is **infinite-server** semantics |
| `capacity` on a transition | **k-server** semantics | GSPN. Distinct from **place capacity**, an unrelated classical notion — the docs must disambiguate |
| `cycletime` | **firing duration** | Timed PN (Ramchandani). NOT Time PN, whose parameter is a firing *interval* |
| `probability` on completion | random switch / probabilistic firing outcome | GSPN |
| `priority` + the allocator | **priority** and random switches, extended by RD's rationing rule | GSPN + ADR 0002 |
| `:nonblock` modality | **read arc** / test arc | classical extension |
| `:conserved` modality | **self-loop** / side condition | classical |
| `:rate` modality | **continuous transition** | hybrid / continuous PN |
| conservation invariant (`S+I+R`) | **P-invariant** (S-invariant) | classical structural analysis |
| resources held over a duration, contended | closest published relative: **Queueing Petri nets** (QPN) | Bause |
| structured token selected by predicate | ~ **binding element** (a transition plus a variable binding) | Coloured PN |

*Token* is canonical and stays. The disambiguation belongs in prose at first use: "a token is a discrete unit of resource sitting in a place — the Petri-net sense, unrelated to language-model tokens".

### Option C detail — tiers

**Tier 1 — exported names (7 in scope).** Each keeps a one-release shim, per the `@deprecate` / `@deprecate_binding` pattern at `src/ReactiveDynamics.jl:172,445-453` (line numbers on the branch base; the draft cited `444-452`, which is the preceding comment line).

| now | new | note |
|---|---|---|
| `@add_species` | `@add_place` | `src/interface/update.jl:3,243` — the draft cited `create.jl`, which is wrong: the macro is defined and exported in `update.jl` |
| `SetSpecies` | `SetMarking` | `src/actions.jl:15,28` — the action sets a place's count; pairs with `MarkingPlot` |
| `ReactantSpec` | `ArcSpec` | `src/ReactiveDynamics.jl:135` — it is an arc record (place, weight, modality), one per LHS/RHS entry |
| `reactant_specs` | `arcs` | `src/ReactiveDynamics.jl:183` |
| `specname` | `placename` | `src/ReactiveDynamics.jl:186` |
| `register_structured_species!` | `register_token_kind!` | `src/interface/agents.jl:3,28` — it registers a *colour set*, not a place (see Open questions) |
| `MarkingPlot` | unchanged | already correct; its docstring loses "species" |

**Tier 1b — names the draft missed, forced by the grep gate.** Three unexported functions are named in demo code, in docstrings, or across `src/` often enough that renaming them silently would break a reader's muscle memory, so they got shims too (`export_old = false`): `get_species` → `get_place`, `set_species!` → `set_place!`, `populate_reactant_specs!` → `populate_arcs!`. Five struct fields moved with NO shim, exactly as ADR 0015 did for its field renames: `ReactionNetwork.reactants` → `.arcs`, `ArcSpec.species` → `.place`, `BaseStructuredToken.species` → `.place`, `PopulationEntry.species` → `.place`, `NetworkGraph.species` → `.places`, plus the `StateDump.tokens` NamedTuple field `species` → `place`. The draft did not mention `@aka`'s `alias_default` entry `:S => :species` either; it survives as `_AKA_LEGACY_NAMES` (`src/interface/update.jl:520`), a one-release authoring alias so `@aka net species = resource` keeps working.

**Tier 2 — store column symbols (206 sites in `src/` + `test/`).** `specName` (80) → `placeName`, `specModality` (27) → `placeModality`, `specCost` (**21**, not 20 as drafted — the total of 206 is right), `specValuation` (18), `specInitVal` (18), `specReward` (14), `specStructured` (13), `specRole` (8), `specInitUncertainty` (7) likewise; the `:S` object symbol may stay (`:S` reads as *places* under either vocabulary) or become `:P` — but `:P` is taken by the parameter object, so keep `:S`. One site is a column symbol CONSTRUCTED at runtime and therefore invisible to a grep for the literal names — `src/interface/update.jl:191`, `Symbol(:spec, uppercasefirst(string(valuation_type)))`, which reaches `specCost`/`specReward`/`specValuation` from `@cost`/`@reward`/`@valuation`. The draft omitted it; missing it would have broken those three macros with no compile error.

**Load-bearing constraint:** two reflection loops filter columns by the literal substring `"spec"` — `src/operators/joins.jl:33` and `src/operators/equalize.jl:58` (`!occursin("spec", string(attr)) && continue`). These must move in the same commit as the column rename, or `@join`/`equalize!` silently stop seeing any place column. This is the one place where the rename is not mechanical, and the reason Tier 2 must be a single atomic change rather than an incremental sweep.

**Tier 3 — serialized keys (stage 2 of the recommendation).** The draft scoped this at "16 sites in `src/`" and listed only the top-level `"species"` array plus the `"set_species"` verb. That undercounts the wire surface, and the omissions are not optional: the draft's own mandate to change the `validate` diagnostic path to `arcs[$i].place` is incoherent unless the array those paths index is itself renamed. The complete set, as implemented:

| wire name | new | where |
|---|---|---|
| top-level `"species"` array | `"places"` | `src/serialize.jl` writer + reader + `validate` |
| top-level `"reactants"` array | `"arcs"` | *omitted by the draft*; renaming it is what makes the `arcs[$i]` diagnostic paths real |
| per-arc `"species"` field | `"place"` | *omitted by the draft* (`src/serialize.jl` arc reader/writer) |
| population entry `"species"` | `"place"` | *omitted by the draft* (`population[]`, ADR 0007 §B) |
| action verb `"set_species"` | `"set_marking"` | plus the `ACTION_VERBS` symbol `:set_species` → `:set_marking` |
| `NodeRef` kind `"species"` | `"place"` | *omitted by the draft*: `REF_KINDS = (:species, :param, :obs)` → `(:place, :param, :obs)` |
| export-bundle manifest `"species"` | `"places"` | *omitted by the draft* (`src/export.jl:133`) |
| per-token trajectory record `"species"` | `"place"` | *omitted by the draft* (`src/export.jl:80`) |
| `validate` diagnostic paths | `places[$i]`, `arcs[$i].place`, `population[$i].place` | user-visible error text, not a format change |
| result-frame column `:species` | `:place` | *omitted by the draft*: `token_trajectory` (`src/analysis.jl`) and the program ledger (`src/ledger.jl`) — a `DataFrame` column name, so user analysis code sees it |
| `@select`/`@advance` field `:species` | `:place` | *omitted by the draft*: a DSL field name (`src/predicates.jl`, `src/solvers.jl`), accepted under both spellings for one release |

**Compatibility, per the maintainer's 2026-08-22 resolution.** The writer emits only the new keys. The reader accepts the legacy `"species"` / `"reactants"` / `"set_species"` spellings for ONE release and emits a `depwarn` when it meets one, and those shims retire together with the Tier-1 name shims — no permanent alias, since a permanent second accepted spelling would rebuild the two-vocabulary problem this ADR exists to remove. Every fixture in the repository is regenerated to the new keys; the legacy-key test stays, reframed as a deprecation test (old keys still load, and warn). The one exception to the warning rule is the `@select`/`@advance` field alias: a `depwarn` there would fire inside the step loop, so that alias is silent for its one release.

**Tier 4 — prose (`spec/` 411, `docs/` 97, `demo/` 154).** ADRs are append-only, so every earlier ADR keeps its wording; CONTRACT §1 (modality truth table), §2, §5.5, §7, §10 are amended in place, since the contract is a living normative document — in practice §8–§15 needed amending too. Add a **Glossary** section to the CONTRACT carrying the dictionary above, and a `docs/src/glossary.md` page to the docs site.

**Tier 5 — local variable names.** `species` as a loop variable, `sp`, `r.species` on the arc record. Word-boundary replace; no shim needed.

**Tier 5 trap, worth recording.** A `\bspecies\b` sweep is blind to underscore-joined identifiers, because `_` is a word character: `species_ixs`, `species_from`, `species_to`, `plain_species`, `free_blocked_species`, `species_modalities`, `extract_reactants` all survive it silently. A second grep family (`_species|species_|_reactant|reactant_`) is required, and it belongs in the verification gate next to the word-boundary one. Two half-completed `sp` → `pl` renames in uncovered branches (`src/ledger.jl`, `src/interface/agents.jl`) also got through the suite as latent `UndefVarError`s; the discipline that catches them is to grep the OLD name repo-wide immediately after each rename, not at the end.

### Migration plan (suite green at every step)

1. Tier 2 in one commit: column symbols + the two `occursin("spec", …)` filters + `ALLATTRS` order preserved. Run the suite.
2. Tier 1: rename the definitions, add the shims, keep the exports. Run the suite plus `test/semantic/exports_resolve.jl`.
3. Tier 5 sweep across `src/`, `test/`, `demo/`; skip `docs/build/`.
4. Tier 4: CONTRACT amendments + glossary, docs sources, demo prose. Re-run the suite — `test/semantic/serialization_ir.jl:375-390` greps source *text* for eval-free violations, and docstring prose has tripped it before.
5. Flip this ADR's `Status:` and its `spec/adr/README.md` row, correcting whatever the draft got wrong along the way.
6. Stage 2 (Tier 3), last and on its own commit: writer emits new keys, reader takes both for one release with a `depwarn`, `validate` path strings, in-repo fixtures regenerated, and a deprecation test that loads a document written under the old keys.
7. `julia -m Runic --inplace src test ext docs/make.jl docs/build_literate.jl demo` (there is no `dev/` on `main`), then a grep gate over `src/`, `test/`, `ext/`, `demo/`, `docs/src/`, `docs/literate/`: no `species`/`reactant`/`spec[A-Z]` token and no `_species`/`species_`/`_reactant`/`reactant_` identifier outside the shim block, the Tier-3 read shim, and the deliberate deprecation tables.

### Consequences

- Large, mostly mechanical diff (~1 400 text sites), localized behind shims for one release exactly as ADR 0015 did; the existing suite is the verification surface.
- The forward/backward asymmetry the draft worried about is not a real cost. With no saved documents outside this repository (Open questions, resolved), the only consumer of the format is RD itself, and every in-repo fixture is regenerated in the same commit. A document written by the new writer will not load on an older release — but nobody holds one, and the older release is one `git tag` away.
- One observable output string changes, and it is not a format key: a `@join` block with no explicit `:alias` now generates the shared place name `shared_place_N` instead of `shared_species_N` (`src/operators/joins.jl:175`). Nothing in the suite, the demos, or any fixture pins the old spelling, and no serialized document carries it; but because `observables` is name-sorted, such a model can see that list reorder. This is the single behavioural consequence of the rename and it is worth stating rather than burying.
- The units ADR (drafted separately) and the retired docs keep the old vocabulary, since ADRs are append-only. The README status table gains a pointer so a reader knows which vocabulary a given ADR predates.
- The two companion papers and the presentation deck must be resynced before publication, or they will teach a vocabulary the code no longer uses.
- The demo tours are the most-read runnable code; their prose rename is the highest-value part of Tier 4 and should not be deferred.

### Rejected micro-alternatives

- **`pool` instead of `place` at the spec level.** Warmer for the business audience, but it is not a term any Petri-net paper uses, so it re-creates the private dialect this ADR removes. Kept as the *tutorial* register instead (Open questions).
- **Keep `species` as a permanent alias.** Two names forever is the current problem with extra steps.
- **Rename `token` to `entity` or `unit`.** Loses the canonical term and the whole literature attached to it, to dodge a confusion that one sentence of prose fixes.
- **Rename `transition`.** Already canonical.
- **`@add_place` → `@place`.** Shorter, but breaks the `@add_*` family symmetry in the authoring macros.
