# Glossary

ReactiveDynamics is a **Petri net** engine, and since v0.3 its API uses standard Petri-net vocabulary — *place*, *marking*, *arc* — in place of the chemical-reaction-network words it inherited from its Catalyst-derived DSL surface (see [ADR 0017](https://github.com/Merck/ReactiveDynamics.jl/blob/main/spec/adr/0017-petri-net-vocabulary.md)). This page is the term dictionary: what each concept is called here, and which published family of Petri nets to read if you want the theory behind it.

**A token is a discrete unit of resource sitting in a place** — the Petri-net sense of the word, unrelated to language-model tokens. When this documentation says "token" it always means the former.

## Two registers, deliberately

The reference pages and the [normative contract](https://github.com/Merck/ReactiveDynamics.jl/blob/main/spec/CONTRACT_DRAFT.md) say **place**, **marking** and **arc**, because those are the names in the code and in the literature. The tutorials and case studies say **resource pool** for the same thing — it is the word a portfolio or capacity discussion actually uses — glossed once at first use as "a place, in Petri-net terms". Both registers describe one object; nothing in the engine distinguishes them.

If your readers want a third register, you do not have to argue with ours: [`@aka`](@ref) renames the objects per model, and its own example is `@aka net place = resource transition = reaction`. Drift where your domain wants it, in the model file, and leave the API canonical.

## Core vocabulary

| Term | What it is in RD |
|:--- |:--- |
| **place** | A resource pool: the `:S` rows of a network, one per named resource. Authored with `@add_place` or by naming it in a transition line. |
| **marking** | The quantity of tokens in each place — the net's state. `state.u` is the current marking; `placeInitVal` (`@prob_init`) declares the initial marking `M₀`. |
| **token** | One discrete unit of resource in a place. *Fungible* tokens are pure quantity (cash, headcount); *structured* tokens are agents with attributes, identity and history. |
| **arc** | One participation of a place in a transition — an [`ArcSpec`](@ref) row carrying `(transition, place, side, multiplicity, modality)`. Input arcs form the left-hand side, output arcs the right-hand side. |
| **arc multiplicity** | How many tokens an arc moves per firing (`2I` is an arc of multiplicity 2). May be a time-varying expression. The literature also says *arc weight*; the field and JSON key are `multiplicity`, because `weight` already names `@choose` alternative weights and the allocator's fill-rate weights. |
| **transition** | A stateful *recipe* that spawns in-flight instances, occupies its input places for a `cycletime`, then completes with probability `probability` and emits its output places. Already the canonical Petri-net word; unchanged. |
| **preset / postset** | The input places of a transition (`•t`) and its output places (`t•`) — the left- and right-hand sides of a reaction line. |
| **firing** | One in-flight instance of a transition running to completion — the engine type [`Firing`](@ref), heaped in `state.ongoing_firings` while the static recipes stay in `state.transitions`. |
| **binding** | Which specific tokens a firing holds — `bound_tokens` (the blocking binds), `nonblock_tokens` (the read-only ones) and `binding` (the `kind => token` assignment) on a live firing. A transition plus a binding is a *binding element*. Distinct from an **arc**, which is static topology and has no runtime instance. |
| **colour set** | The attribute schema of a structured place's tokens. Declared with `@structured_token` and registered by [`register_token_kind!`](@ref). |
| **reaction network** | A Petri net presented in *reaction notation* — the arrow-form lines you author (`3*@conserved(scientist) + @rate(budget) --> compound`). Under the standard correspondence a reaction network and a Petri net are the same object: species↔places, reactions↔transitions. RD keeps the chemistry register for the **notation** (`@reaction_network`, "reaction line") and Petri vocabulary for the **object model** — you author in reaction lines and get a net of places, transitions and arcs. It is not a *chemical* reaction network: kinetics is just the archetypal instance of the ontology. |

## Which literature answers a question

The mapping from RD's features onto the published Petri-net extensions, so you know where to look.

| RD concept | Standard term | Family / note |
|:--- |:--- |:--- |
| fungible resource pool | **place**; its count is that place's **marking** | classical place/transition net; `M₀` is the initial marking |
| structured-token pool | place with a **colour set**; its tokens are distinguishable | Coloured Petri nets (Jensen) — token attributes are *colours* |
| left-/right-hand-side entry | input **arc** / output arc; multiplicity is the **arc weight** (inscription) | classical |
| LHS / RHS multiset | **preset** `•t` / **postset** `t•` | classical |
| transition | **transition** | already canonical |
| `rate` expression | **firing rate**, marking-dependent. RD's Poisson intensity with unbounded concurrency is **infinite-server** semantics | Stochastic PN / GSPN |
| `capacity` on a transition | **k-server** semantics — distinct from *place capacity*, an unrelated classical notion RD does not implement | GSPN |
| `cycletime` | **firing duration** (not a firing *interval*: this is Timed PN, not Time PN) | Timed PN (Ramchandani) |
| `probability` on completion | random switch / probabilistic firing outcome | GSPN |
| `priority` + the allocator | **priority** and random switches, extended by RD's rationing rule | GSPN + [ADR 0002](https://github.com/Merck/ReactiveDynamics.jl/blob/main/spec/adr/0002-priority-weighted-allocation.md) |
| `:nonblock` modality | **read arc** / test arc | classical extension |
| `:conserved` modality | **self-loop** / side condition | classical |
| `:rate` modality | **continuous transition** | hybrid / continuous PN |
| a conservation invariant (`S+I+R` constant) | **P-invariant** (S-invariant) | classical structural analysis |
| resources held over a duration under contention | closest published relative: **Queueing Petri nets** (QPN) | Bause |
| a structured token selected by predicate | ≈ **binding element** (a transition plus a variable binding) | Coloured PN |

## Renamed in v0.3

Each retired name still resolves for one release and warns; see [ADR 0017](https://github.com/Merck/ReactiveDynamics.jl/blob/main/spec/adr/0017-petri-net-vocabulary.md) for the full migration.

| Retired | Use instead |
|:--- |:--- |
| `@add_species` | [`@add_place`](@ref) |
| `SetSpecies` | [`SetMarking`](@ref) |
| `ReactantSpec` | [`ArcSpec`](@ref) |
| `reactant_specs` | [`arcs`](@ref) |
| `specname` | [`placename`](@ref) |
| `register_structured_species!` | [`register_token_kind!`](@ref) |
| the `specName`/`specInitVal`/`specModality`/… store columns | `placeName`/`placeInitVal`/`placeDefaultModality`/… |
| `stoich`, on an arc | `multiplicity` |

The serialized [JSON document](reference/json_schema.md) renamed with them. The exporter emits only the new keys, so re-exporting a document written before v0.3 migrates it; the loader reads the retired keys for one release and warns.

| Retired key | Use instead |
|:--- |:--- |
| top-level `"species"` array | `"places"` |
| top-level `"reactants"` array | `"arcs"` |
| an arc's `"species"` | `"place"` |
| a `population[]` entry's `"species"` | `"place"` |
| the action verb `"set_species"` | `"set_marking"` |
| a `ref` node's `"kind": "species"` | `"kind": "place"` |
| the result-frame / export column `:species` | `:place` |
| the `@select`/`@advance` field `:species` | `:place` (accepted silently — a warning would fire once per tick) |
| an arc's `"stoich"` | `"multiplicity"` |

[ADR 0018](https://github.com/Merck/ReactiveDynamics.jl/blob/main/spec/adr/0018-firing-vs-transition.md) landed in the same release and fixed one further inversion: the type for an in-flight instance was called `Transition`, while the static transitions sat in a field called `transition_recipes`. The type keeps a forwarding alias for one release; the renamed *fields* cannot (a struct field carries no forwarding binding), so a model that reads them updates at once. No serialized key moved — the saved document already said `"transitions"`.

| Retired | Use instead |
|:--- |:--- |
| the type `ReactiveDynamics.Transition` | [`Firing`](@ref) — one in-flight execution |
| `state.transition_recipes` | `state.transitions` — the static transition table, compiled once |
| `state.transitions` (the per-tick evaluated snapshot) | `state.sampled_transitions` — the same table realized for the current tick, which is what `state[i, :trans…]` reads |
| `state.ongoing_transitions` | `state.ongoing_firings` |
| a token's `bound_transition` | `bound_firing` |
| the in-model marker `@transition` in an attribute expression | `@firing` (the old spelling still compiles) |
