# ADR 0017: Should the engine retire "species" for standard Petri-net vocabulary?

Status: Draft — 2026-08-21. Proposed; no code has moved. Deciders: maintainer + package author.

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

Take **Option C**, in two stages: prose and identifiers first, serialized keys second behind a read alias. The trade-off accepted is a permanent one-line alias in the loader in exchange for a single vocabulary across code, spec, and file format. Keep **token** — it is the canonical Petri word, no synonym is better, and one gloss at first use settles the language-model confusion. Confidence: moderate — the contested half is the file format, not the rename. Reverse to Option D on finding one saved model outside this repository that cannot be re-emitted.

## What changes

A model author writes `@add_place` instead of `@add_species`, and reads `SetMarking` where an action sets a pool count. Attribute names in authored models change spelling — `specName` becomes `placeName` — so any model that names an attribute directly needs a one-word edit. Old spellings work for one release with a warning, exactly as in ADR 0015.

Documentation keeps two registers deliberately: the spec and the API say *place*, *marking*, *arc*; the tutorials say *resource pool*, glossed once as "a place, in Petri-net terms". A tutorial reader needs no Petri theory; a spec reader meets no private dialect.

One-time work: about a day of renaming plus half a day re-reading the docs, one person in a single pass, verified by `julia --project=. -e 'using Pkg; Pkg.test()'` staying green and a grep gate confirming the old spellings survive only inside the shims. CI is unaffected.

## Open questions

- Does any saved model exist outside this repository — the single fact that decides C versus D? (Maintainer, before stage 2.)
- Should the shims be dropped at the next minor release, in step with the ADR 0015 shims, or kept a cycle longer? (Package author, at release time.)
- Is *place* or *pool* the better spec-level word, given that the domain audience reads "place" as a location? (Recommendation: `place` in the spec and API, `pool` in tutorial prose — resolve in review.)
- Does the rename wait for the two companion papers to fix their terminology first, so all three land consistent? (Maintainer.)
- Do the structured-token registry names move too (`register_structured_species!` → a token-kind verb), or is that a separate coloured-net naming pass? (Package author.)

Technical detail, evidence, the full term dictionary, and the ordered migration: see Appendix.

## Appendix

### Evidence for Problem

Occurrences of the word `species` (case-insensitive, word boundary), 2026-08-21 on `rework`:

| tree | count |
|---|---|
| `src/` | 448 |
| `spec/` | 438 |
| `test/` | 177 |
| `demo/` | 154 |
| `docs/src/` | 97 |
| `docs/literate/` | 66 |

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

**Tier 1 — exported names (7 in scope).** Each keeps a one-release shim, per the `@deprecate` / `@deprecate_binding` pattern at `src/ReactiveDynamics.jl:172,444-452`.

| now | new | note |
|---|---|---|
| `@add_species` | `@add_place` | `src/interface/create.jl` |
| `SetSpecies` | `SetMarking` | `src/actions.jl:15` — the action sets a place's count; pairs with `MarkingPlot` |
| `ReactantSpec` | `ArcSpec` | `src/ReactiveDynamics.jl:131` — it is an arc record (place, weight, modality), one per LHS/RHS entry |
| `reactant_specs` | `arcs` | `src/ReactiveDynamics.jl:179` |
| `specname` | `placename` | `src/ReactiveDynamics.jl:186` |
| `register_structured_species!` | `register_token_kind!` | `src/interface/agents.jl:3` — it registers a *colour set*, not a place (see Open questions) |
| `MarkingPlot` | unchanged | already correct; its docstring loses "species" |

**Tier 2 — store column symbols (206 sites in `src/` + `test/`).** `specName` (80) → `placeName`, `specModality` (27) → `placeModality`, `specCost` (20), `specValuation` (18), `specInitVal` (18), `specReward` (14), `specStructured` (13), `specRole` (8), `specInitUncertainty` (7) likewise; the `:S` object symbol may stay (`:S` reads as *places* under either vocabulary) or become `:P` — but `:P` is taken by the parameter object, so keep `:S`.

**Load-bearing constraint:** two reflection loops filter columns by the literal substring `"spec"` — `src/operators/joins.jl:33` and `src/operators/equalize.jl:58` (`!occursin("spec", string(attr)) && continue`). These must move in the same commit as the column rename, or `@join`/`equalize!` silently stop seeing any place column. This is the one place where the rename is not mechanical, and the reason Tier 2 must be a single atomic change rather than an incremental sweep.

**Tier 3 — serialized keys (16 sites in `src/`, stage 2 of the recommendation).** `"species"` → `"places"` in `src/serialize.jl:102,250,380,687,688,696,740,741,742,780,973` and `src/export.jl:80,133`; `"set_species"` → `"set_marking"`. The reader accepts either key permanently (`get(d, "places", get(d, "species", []))`); the writer emits only the new one. `validate` diagnostics change their path strings (`reactants[$i].species` → `arcs[$i].place`), which is user-visible text in error messages, not a format change.

**Tier 4 — prose (`spec/` 438, `docs/` 163, `demo/` 154).** ADRs are append-only, so every earlier ADR keeps its wording; CONTRACT §1 (modality truth table), §2, §5.5, §7, §10 are amended in place, since the contract is a living normative document. Add a **Glossary** section to the CONTRACT carrying the dictionary above, and one to the docs site.

**Tier 5 — local variable names.** `species` as a loop variable, `sp`, `r.species` on the arc record. Word-boundary replace; no shim needed.

### Migration plan (suite green at every step)

1. Tier 2 in one commit: column symbols + the two `occursin("spec", …)` filters + `ALLATTRS` order preserved. Run the suite.
2. Tier 1: rename the definitions, add the shims, keep the exports. Run the suite plus `test/semantic/exports_resolve.jl`.
3. Tier 5 sweep across `src/`, `test/`, `demo/`; skip `docs/build/`.
4. Tier 4: CONTRACT amendments + glossary, docs sources, demo prose. Re-run the suite — `test/semantic/serialization_ir.jl:375-390` greps source *text* for eval-free violations, and docstring prose has tripped it before.
5. Stage 2 (Tier 3), gated on the Open question about outside models: writer + reader alias + `validate` path strings, with a round-trip test that loads a fixture written under the old keys.
6. `julia -m Runic --inplace src test ext dev docs/make.jl demo`, then a grep gate: no `species`/`reactant`/`spec[A-Z]` token outside the shim block and `docs/build/`.

### Consequences

- Large, mostly mechanical diff (~1 400 text sites), localized behind shims for one release exactly as ADR 0015 did; the existing suite is the verification surface.
- Model documents written before stage 2 keep loading forever; documents written after it will not load on an older release. That asymmetry is the price of the key rename and the reason it is staged separately.
- The units ADR (drafted separately) and the retired docs keep the old vocabulary, since ADRs are append-only. The README status table gains a pointer so a reader knows which vocabulary a given ADR predates.
- The two companion papers and the presentation deck must be resynced before publication, or they will teach a vocabulary the code no longer uses.
- The demo tours are the most-read runnable code; their prose rename is the highest-value part of Tier 4 and should not be deferred.

### Rejected micro-alternatives

- **`pool` instead of `place` at the spec level.** Warmer for the business audience, but it is not a term any Petri-net paper uses, so it re-creates the private dialect this ADR removes. Kept as the *tutorial* register instead (Open questions).
- **Keep `species` as a permanent alias.** Two names forever is the current problem with extra steps.
- **Rename `token` to `entity` or `unit`.** Loses the canonical term and the whole literature attached to it, to dodge a confusion that one sentence of prose fixes.
- **Rename `transition`.** Already canonical.
- **`@add_place` → `@place`.** Shorter, but breaks the `@add_*` family symmetry in the authoring macros.
