# ADR 0018: Should the in-flight instance be called a *firing*, and the static recipes *transitions*?

Status: **Accepted + Implemented** — 2026-08-24. Drafted from [ADR 0017](0017-petri-net-vocabulary.md)'s amendment §A7, which identified this as the one category error the Petri-net rename left in place and deferred it here because it is not a chemistry word. Accepted as **Option C** and landed the same day, in one commit on branch `rename/petri-vocabulary`, with the maintainer's four open questions resolved in §Resolutions below. Deciders: maintainer + package author.

Date: 2026-08-24

Relates to: [ADR 0017](0017-petri-net-vocabulary.md) (retired `species`/`reactant`, kept `transition` unchanged on the grounds that it was already canonical — true of the word, not of what it labels), [ADR 0015](0015-post-acsets-naming.md) (the previous naming window; its one-release shims and this rename should land in one migration, not two), [ADR 0012](0012-algebraicagents-integration.md) (the in-flight instance is an AlgebraicAgents node, so its *name string* is part of the agent path — unchanged here), [ADR 0007](0007-interface-and-initial-state.md) (`dump_state` refuses to run while instances are in flight, and says so in its error text). Pure naming change: no runtime semantics, no allocator behaviour, no saved-document change.

## Problem

Two names are swapped. The type that models one in-flight instance is called `Transition`, and the static transitions it is spawned from live in a field called `transition_recipes` — "recipes" being a word invented to escape the collision. Standard usage is the reverse: a transition is the static object, and one execution of it is a **firing**.

The plain word `transitions` is not free either: it is taken by a third table, the per-tick snapshot of every transition's attributes evaluated for the current tick. So the swap is a three-way one, and freeing `transitions` for the static rows means naming that table too.

The documentation already knows this. "Firing" appears 117 times across the repository — 54 in `spec/`, 31 in `docs/` — while no name in the code spells it. A reader who learns *firing* from the tutorial cannot grep for it.

Doing nothing keeps the exact translation tax ADR 0017 was written to remove, one layer down: every review conversation and every docstring pays "a `Transition` is a firing; the transitions are the recipes." It also strands the amendment, which committed this to a follow-up ADR.

## Options

**Option A — Do nothing.** Keep `Transition` for the instance and `transition_recipes` for the static rows.

- **Gain:** zero effort and zero churn, one release after a rename the maintainer has just absorbed.
- **Cost:** the codebase keeps a name that means the opposite of what the literature and our own prose mean by it.
- **Breaks when:** an outside reader works from the papers or the glossary into the source, which is the normal path.

**Option B — Prose only.** A glossary row records that the instance type is a firing; no code moves.

- **Gain:** an afternoon's work; nothing can break for a caller.
- **Cost:** documents the defect instead of fixing it, and the generated docstrings keep contradicting the page that explains them.
- **Breaks when:** someone writes a new subsystem and reasonably follows the code's vocabulary rather than the glossary's.

**Option C — Rename the whole family.** The instance type becomes `Firing`, the live vector `ongoing_firings`, the static field `transitions`, the per-tick snapshot it displaces `sampled_transitions`, and a token's back-pointer `bound_firing`.

- **Gain:** one vocabulary at both layers; *firing* becomes greppable; the invented word "recipes" disappears.
- **Cost:** about 250 text sites; one more deprecation binding, because every structured-token example names the type; and a field rename for anyone reading the live vector.
- **Breaks when:** a caller reads a renamed field and cannot be updated, since a struct field cannot carry the cheap alias a type can.

**Option D — Rename the type only.** `Transition`→`Firing`; both field names stay.

- **Gain:** fixes the type that appears in printed output and in error messages, with no field migration at all.
- **Cost:** leaves `ongoing_transitions` holding a vector of firings, so the swap survives in the field a reader actually touches.
- **Breaks when:** immediately, for the reader in the Problem — the field is what demos and tutorials read.

## Recommendation

Adopt **Option C**, and land it in the same release window as ADR 0017's shim removal. The trade-off: model authors absorb one more field rename now, in exchange for a single migration rather than two in consecutive releases. Confidence is high — the blast radius is measured, the saved-document format does not change, and the renamed type takes the same one-line alias ADR 0017 gave its own renamed names.

One observation reverses this: if the maintainer wants ADR 0017's deprecations to expire before any further renaming, hold C for the release after next and take nothing in the meantime.

## What changes

In `src/state.jl` the instance type becomes `Firing`, the live vector `ongoing_firings`, and the static attribute columns `transitions` — which displaces the per-tick snapshot of evaluated attributes to `sampled_transitions`, the name the function that fills it already used. In `src/interface/agents.jl` a token's back-pointer becomes `bound_firing`. Roughly 250 text sites in all.

Three things deliberately do not move: the store's transition columns, because renaming them breaks the substring filter at `src/operators/joins.jl:40` — the reflection trap ADR 0017 hit at `place`; the saved-document format, whose key is already `"transitions"`; and agent paths, built from a name string rather than the type.

Two things change for a model author. Anyone reading the in-flight vector renames that field access — four demo sites and one tutorial do so today. Anyone naming the type when declaring a token's bond history writes `RD.Firing`, and the old spelling warns for one release through an alias in its own block at `src/ReactiveDynamics.jl:531-542`, directly below ADR 0017's shim block (`:486-529`).

One-time work: roughly half a day, one commit, one row in the glossary's renamed table. Verify with `julia --project=. -e 'using Pkg; Pkg.test()'` (≈4 minutes).

## Resolutions

The four questions the brief left open, decided on acceptance so nothing dangles:

- **Ship now, in ADR 0017's window.** The rename landed in the same branch and release as the Petri-net vocabulary, so model authors absorb one migration rather than two.
- **No read alias for the renamed fields.** A glossary row and the renamed table carry them. Overloading `getproperty` on `ReactionNetworkProblem` would put a branch in the step loop's hot field access, and the fields were never exported names.
- **`Firing` stays unexported**, reached as `ReactiveDynamics.Firing` — the spelling every structured-token example already writes. The old name resolves for one release through `@deprecate_binding Transition Firing false`.
- **Column-family reflection keeps its substring test**, and stays a known trap rather than a blocker: the three filters that read a column prefix (`src/operators/joins.jl`, `src/operators/equalize.jl`) are documented in place. Replacing them with a declared per-column role is a follow-up whenever the store schema is next revised — no naming decision waits on it.

Technical detail, evidence, and agent-executable steps: see Appendix.

## Appendix

### Evidence for Problem

Current state, verified 2026-08-24 on branch `rename/petri-vocabulary` at the ADR 0017 amendment commit.

The two swapped names, in source:

```julia
# src/state.jl:48 — the RUNTIME instance, holding the canonical STATIC word
@aagent struct Transition
    i::Int                                            # originating :T row
    trans::Dict{Symbol, Any}                          # per-instance attributes
    bound_tokens::Vector{AbstractAlgebraicAgent}      # blocking binds   (ADR 0017 amendment)
    nonblock_tokens::Vector{AbstractAlgebraicAgent}   # read-only binds
    binding::Vector                                   # kind => token assignment
    t::Float64                                        # spawn time
    q::Float64                                        # allocated quantity
    state::Float64                                    # progress through the cycle
end

# src/state.jl:89  — the STATIC transitions, under an invented word
transition_recipes::Dict{Symbol, Vector}
# src/state.jl:101 — the live instances
ongoing_transitions::Vector{Transition}
```

`transition_recipes` is materialized at construction from the `:T` rows: `compile_attrs` returns it as its second value under the local name `transitions`, which `src/solvers.jl:978` then aliased to `transition_recipes`.

**Correction to this ADR's first draft.** That draft read the alias at `solvers.jl:978` as evidence that `transitions` was a free name. It is not. The same local is REBOUND thirty lines later (`solvers.jl:1001-1006`) to a fresh `Dict{Symbol, Vector}(a => [] for a in transitions_attrs)`, and *that* dict is what the constructor passes to the `transitions` field. So the state holds two distinct transition tables:

| Field | What it holds | Lifetime |
|:--- |:--- |:--- |
| `transition_recipes` | the compiled static columns, one entry per `:T` row (`:trans` recipe body, `transRate`, `transActivated`, `transGuard`, `transHash`, `transToSpawn`, …) | built once at construction, mutated only by the append-only runtime mutators and the soft gates |
| `transitions` | the same column family with each value EVALUATED for the current tick, plus `transLHS` (resolved arcs), `transRHS` and `transFiring` | emptied and refilled every tick by `sample_transitions!` (`src/state.jl`); what `state[i, :trans*]` reads |

Renaming `transition_recipes` → `transitions` therefore requires renaming the per-tick table first. It becomes **`sampled_transitions`** — the word its producer (`sample_transitions!`) and its reader (`get_sampled_transition`) already use. Collapsing both onto one name is not merely ambiguous, it is a silent corruption: the sweep that did so left `push!(state.transitions[k], state.transitions[k][i])` in `sample_transitions!` and dropped the static dict entirely at `solvers.jl:978`, and it fails at load time on `duplicate field name: "transitions" is not unique`, which is how it was caught.

Occurrence counts (`git grep -o`, tracked files only):

| Name | src | test | spec | docs | demo | total |
|:--- |---:|---:|---:|---:|---:|---:|
| `Transition` (as a word) | 24 | — | — | — | — | 92 all-dirs |
| `ongoing_transitions` | 28 | 9 | 10 | 1 | 4 | 52 |
| `transition_recipes` | — | — | — | — | — | 27 |
| `bound_transition` (+ its accessor pair) | 18 | 0 | 20 | 6 | 3 | 47 |
| `firing` / `fires` (prose, case-insensitive) | 14 | 6 | 54 | 31 | 12 | 117 / 79 |

The last row is the argument: the concept is already named *firing* everywhere except in the code that implements it.

Reachability of the type — it is **not** exported (every `export` line in `src/` was checked; none names `Transition`), but it is nonetheless part of the public surface in practice, which is the fact that decides the shim question:

- It is **documented** as an API name: `docs/src/reference/structured_tokens.md:21,24` explains it and lists it in the reference's `@docs` block.
- It is **written by hand in every structured-token example**, qualified, because `BaseStructuredToken`'s `past_bonds` field is typed by it — `Tuple{Symbol, Float64, ReactiveDynamics.Transition}[]` at `demo/agentic_pipeline/agentic_pipeline.jl:75,382`, `demo/bd_acquisition/host.jl:42`, `demo/introspection_tour/introspection_tour.jl:85`, `docs/literate/case_studies/{inlicensing_value.jl:66, kill_a_program.jl:49, marginal_scientist.jl:48}`, `docs/literate/deep_dives/serialization.jl:31`, `docs/literate/tutorials/advanced.jl:59`. Ten sites, all in authoring position.
- It reaches users two more ways: as the element type of `state.ongoing_transitions`, and by AlgebraicAgents traversal (`getagent`, `entangle!`) of a live hierarchy, where the type name shows in printed output.

So "unexported" is not an argument for skipping a shim here. A renamed type takes `@deprecate_binding Transition Firing` — one line in the kind of block ADR 0017 already established at `src/ReactiveDynamics.jl:486-529`, which is why Option C's cost over Option D is a field rename and not a type rename.

De-facto public field reads, which a rename must update or grandfather:

- `demo/agentic_pipeline/agentic_pipeline.jl:660,692`
- `demo/core_engine_tour/core_engine_tour.jl:142`
- `docs/literate/tutorials/expert.jl:300`
- `test/semantic/{allocation_conservation_lifecycle.jl:246, initial_state.jl:191,219, modality_genesis.jl:348,370,373, reference_models.jl:401,415, determinism_composition_bugs.jl:108}`

User-visible strings that would be reworded: `src/interface/checkpoint.jl:51-52` (`"dump_state: N in-flight transition(s) — dump is …"`), and the `:rate`-modality construction error at `src/solvers.jl:115`.

### Option C detail

Exact mechanical mapping:

| Now | Proposed | Sites |
|:--- |:--- |:--- |
| `struct Transition` | `struct Firing` | `src/state.jl:49` + 91 references |
| `ongoing_transitions` | `ongoing_firings` | 52 |
| `transition_recipes` | `transitions` | 27 |
| `transitions` (the per-tick evaluated snapshot) | `sampled_transitions` | 21 in `src/`, 0 elsewhere |
| `transitions_attrs` (its key set, a construction local) | `sampled_attrs` | `src/solvers.jl:1001` |
| `Transition[]` (constructor call) | `Firing[]` | `src/solvers.jl:998` |
| `bound_transition` (token field) | `bound_firing` | `src/interface/agents.jl:19` + 46 references |
| `get_bound_transition` / `set_bound_transition!` | `get_bound_firing` / `set_bound_firing!` | `src/interface/agents.jl:166,167`; both unexported |
| — | `@deprecate_binding Transition Firing false` | new block in `src/ReactiveDynamics.jl:531-542`, below ADR 0017's shim block (`:486-529`) |
| `TransitionNode` (viz) | unchanged | `src/visualize.jl:31` — it draws the *static* node, correctly named |
| the compiled closure's second parameter | `transition` → `firing`, and the in-model marker `@transition` gains the canonical spelling `@firing` (the old one still lowers to the same identifier) | `src/compilers.jl`, `src/actions.jl`, `src/exprnode.jl` |
| `trans_` (the `Firing`-valued step-loop local, trailing-underscored to dodge the `trans` column family) | `firing` | 30 occurrences in `src/solvers.jl` and nowhere else, now 0 |

Held back on purpose, with the reason each one is load-bearing:

1. **`SCHEMA` `:T` and the `trans*` columns** (`src/ReactiveDynamics.jl:62-72`, `transPriority`, `transRate`, `transCycleTime`, `transProbOfSuccess`, `transCapacity`, `transMaxLifeTime`, `transPreAction`, `transPostAction`, `transMultiplier`, `transName`, …). These name the static transition and are already right. Renaming them would also break `src/operators/joins.jl:40`, `!occursin("trans", string(attr)) && continue`, which selects the transition column family by substring — the mirror of the `"place"` filter at `joins.jl:33` and `equalize.jl:58` that ADR 0017 had to fix in its own Tier 2 commit (`adaa0ed`). `SCHEMA` declaration order is load-bearing too (`ALLATTRS`).
2. **The serialized document.** `"transitions"` is already the top-level key (`src/serialize.jl:156,289,765`), and per-arc/per-statement `"transition"` keys refer to the static object (`:493,494,528,530,786,1014`). No wire-format change, therefore no `_legacy_key` read alias and no second deprecation set — the sharpest contrast with ADR 0017, whose Tier 3 needed both.
3. **AA agent names and paths.** A live instance is entangled under the string `"$(state[i, :transName])_@$(state.t)"` (`src/solvers.jl:374`), built from the store column, not from the Julia type. Saved agent paths, `getagent` calls and wiring diagrams are unaffected.
4. **The state dump.** `bound_transition` is on the `_PROTOCOL_FIELDS` **exclusion** list (`src/interface/checkpoint.jl:37-40`), so it is never written into a `StateDump` — only a token's modeling attributes are. It is `nothing` at dump time anyway, because `dump_state` refuses to run unless the in-flight set is empty (ADR 0007 §C). Renaming the field therefore cannot invalidate a saved dump; the tuple in the list gets the new symbol and nothing else moves.

### Rejected micro-alternatives

- **`Firing` vs `FiringInstance` / `TransitionInstance` / `Instance`.** `Firing` is the literature's noun for exactly this object (Murata; Jensen's *binding element* is the transition-plus-binding pair, and its execution is a firing). The `*Instance` forms restate the type's role in its own name, and `Instance` alone is unusably generic.
- **`Firing` vs `Event`.** Rejected: `state.log` already carries an event channel, and ADR 0010's rules layer uses "event" for the decision channel. Reusing it would rebuild the ambiguity this ADR removes.
- **Keeping `transition_recipes` while renaming the type.** This is Option D in the brief. Kept as a foil, because the field is what callers touch.
- **`transition_defs` / `transition_specs` for the static field.** Both are placeholder-flavoured, and `*Spec` is now spoken for by `ArcSpec` (ADR 0017 Tier 1). The plain word is available once `Transition` stops occupying it, which is the point of the pairing.
- **`bound_firing` vs `firing` for the token back-pointer.** The `bound_` prefix earns its place: it pairs with `isblocked` (`src/interface/agents.jl:158`) and with the ADR 0017 amendment's `bound_tokens`/`nonblock_tokens` on the other end of the same relation, so the two sides of a binding read as a pair. A bare `firing` field would also collide conceptually with the `past_bonds` history, which holds firings the token is no longer bound to.
- **A `getproperty` alias for the old field names.** Deferred to the first open question. Cost: `Transition`/`Firing` and `ReactionNetworkProblem` are `@aagent` structs, so overloading `getproperty` on the state would add a branch to a hot field access in the step loop. A glossary row plus the renamed-table entry is the cheaper default, given the field is not an exported name and ADR 0017's shim precedent covers exported names only.

### Migration steps (agent-executable)

Single commit; do the three renames in one pass so the suite is never red between them. Two hard constraints on the sweep, both learned from ADR 0017's execution:

1. **Drive it from `git ls-files`.** Untracked generated artifacts under `demo/*/output/` match these patterns and must not be rewritten. Also note that `zsh` does not word-split an unquoted variable, so a `perl -pi -e '…' $FILES` form passes the whole list as one filename and silently does nothing.
2. **Exclude `spec/adr/0001`–`0017` entirely.** ADRs are append-only, and six of them name the type in prose (`0001:11`, `0003:48,71`, `0006:13`, `0007:83,89,132`). They keep the old spelling by the same rule that left the retired CRN vocabulary in place below the ADR 0017 boundary.

```bash
cd ~/worktrees/rd-review/petri-vocabulary   # or a fresh worktree off the merge base

# 0. FIRST free the target name: the per-tick snapshot vacates `transitions`. src/ only — no test,
#    doc or demo reads it. Do this before step 1 or the two tables collapse onto one field.
perl -pi -e 's/\bstate\.transitions\b/state.sampled_transitions/g;
             s/\bwork\.transitions\b/work.sampled_transitions/g;
             s/\btransitions_attrs\b/sampled_attrs/g' \
  src/state.jl src/solvers.jl src/visualize.jl
#    then by hand: the struct field (`src/state.jl`), the construction local and the constructor
#    argument (`src/solvers.jl:1006` + the call), and — inside `sample_transitions!` — every read of
#    the STATIC table, which must stay `state.transitions` while every push target becomes
#    `state.sampled_transitions`. `src/actions.jl`'s four `transActivated`/`transGuard`/`transHash`
#    writes address the static table and keep the plain name.

# 1. the three field names + the accessor pair — mechanical, no ambiguity, safe to sweep
git ls-files -- src test docs demo spec | grep -v '^spec/adr/00' | while IFS= read -r f; do
  perl -pi -e 's/\bongoing_transitions\b/ongoing_firings/g;
               s/\btransition_recipes\b/transitions/g;
               s/\b(get_|set_)?bound_transition\b/${1}bound_firing/g' "$f"
done

# 2. the type name — sweep CODE only; \bTransition\b is exact there
git ls-files -- src test demo docs/literate | grep -E '\.jl$' | while IFS= read -r f; do
  perl -pi -e 's/\bTransition\b/Firing/g' "$f"
done
git grep -nw 'Firing' -- src test demo docs/literate   # review every hit before committing

# 3. hand-edit the four prose/reference sites where BOTH senses appear
#    spec/CONTRACT_DRAFT.md   — :127,176,182,247,500 mean the INSTANCE (rename);
#                               :246,375,447,472,480,553 mean the static TABLE (keep)
#    spec/INVENTORY.md        — :22,150,172,199 mean the instance (rename); :17 names the type in
#                               a load-order note (rename)
#    docs/src/reference/structured_tokens.md:21,24 — the documented API name (rename + shim note)
#    src/state.jl:41-46       — the docstring, which already says "in-flight transition instance"

# 4. the shim and the two user-facing strings
#    src/ReactiveDynamics.jl:531-542   add `@deprecate_binding Transition Firing false` — its own
#                                      block, below ADR 0017's shim block at :486-529
#    src/interface/checkpoint.jl:51-52 "in-flight transition(s)" → "in-flight firing(s)"
#    src/solvers.jl:115                KEEP "transition $(…[:transName])" — it names the static row

# 5. the Firing-VALUED parameter name, in every signature that takes one: `transition` → `firing`
#    in src/{solvers,ledger,predicates,visualize,state}.jl, src/interface/agents.jl, and the
#    compiled-closure/action path src/{compilers,actions,exprnode}.jl. Leave alone the parameters
#    that receive a transition NAME (`priority(a, transition)`, `Activate.transition::Symbol`).
#    Same step, same new name: the `Firing`-valued LOCAL `trans_` → `firing`, 30 occurrences in
#    src/solvers.jl and nowhere else. Record it here — it is not covered by the `\bTransition\b`
#    sweep of step 2, and a prose sweep that misses it leaves `trans_` cited in the CONTRACT.

julia -m Runic --inplace src test ext docs/make.jl demo
julia --project=. -e 'using Pkg; Pkg.test()'
```

Step 3 is the real work and cannot be automated: `Transition` in `spec/` and in the CONTRACT means the static model table about as often as it means the runtime instance (§6.6 defines a `Transition` table with a `Transition.id` FK — that one is correctly named and stays). A blanket prose sweep would invert the very distinction this ADR draws. `\bTransition\b` is exact in `.jl` sources, though: it does not match `TransitionNode` (`src/visualize.jl:31`, the static draw node), `@append_transitions` or `abstract_transitions`, so step 2 needs no repair pass — only a review.

### Verification

- Full suite green at the count current when this lands. The exact count is environment-dependent, because three `@test_skip` environment guards survive in the suite: **835 pass / 0 fail / 0 broken / 0 skipped** locally when a Graphviz backend is installed, and **833 pass / 0 fail / 1 broken / 1 skipped** on CI, where the Graphviz guard at `test/semantic/visualization.jl:132` is taken (Graphviz is a system binary the runners do not install; the Plots and Arrow guards at `visualization.jl:203` and `analysis_observability.jl:333` never trip, since both packages are declared in `test/Project.toml`). Against 828 at the ADR 0017 amendment commit, this ADR's own gate accounts for +6. No new tests are required for the rename itself — it is covered by every existing lifecycle test that reads the live vector — but the retired-vocabulary gate in `test/semantic/exports_resolve.jl` gains `transition_recipes`, `ongoing_transitions` and `bound_transition`, so neither the invented word nor a swapped name can reappear in `src/`. The `Transition` → `Firing` alias is asserted resolvable in the same file.
- `git grep -now 'firing' -- src` should be non-zero afterwards; it is 14 today and all of it is prose.
- Round-trip a saved document written before the rename (`demo/bd_acquisition/model.rdj.json`) and diff the re-export: it must be byte-identical modulo key ordering, since no serialized key moves.
- Load the package before running anything else. The two-transition-table collision is a *load-time* error (`duplicate field name`), so `julia --project=. -e 'using ReactiveDynamics'` is the cheapest check that step 0 was done correctly.
