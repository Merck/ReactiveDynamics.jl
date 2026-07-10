# Refinement & open-port composition tour (the granularity ladder)

A single, runnable, literate walkthrough of ReactiveDynamics' hierarchical **refinement** and open-port **composition** layer (ADR 0009 / CONTRACT §11) — the vertical granularity axis the flat §7 composition was otherwise silent on. The [core_engine_tour](../core_engine_tour) demo shows the flat modeling vocabulary and [bd_acquisition](../bd_acquisition) shows a full structured-token application; this one shows how to build a process COARSELY, then substitute a FINER sub-process for one step without disturbing the rest of the model. It answers the maintainer's framing requirement directly: "a modeling framework suitable for business processes — compact/expressive definition, compositionality, and various levels of granularity with more refined dynamics possibly substituted." Every construct is lifted from the passing semantic test `test/semantic/refinement_composition.jl` and the `src/operators/refine.jl` docstrings — the tour invents no API.

## Run it

```bash
julia --project=demo/refinement_tour demo/refinement_tour/refinement_tour.jl
```

First time (resolve the demo-local env — takes a few minutes):

```bash
julia --project=demo/refinement_tour -e 'using Pkg; Pkg.instantiate()'
```

The script is literate: every section opens with a block comment explaining the modeling idea, then runs the code, then `println`-narrates the result, so running it once tells the whole story top to bottom. It is fast (well under a minute after compilation) and writes nothing to disk — it is a structural tour, printing the transition lists and the promoted `ReactantSpec` incidence table rather than rendering artifacts.

### Why a demo-local `Project.toml`

This tour is structural authoring, so its environment is deliberately thin: a path-dev'd `ReactiveDynamics` plus `AlgebraicAgents` (from which `simulate` is re-exported) and `Printf` for the narrated output. It needs no `Plots` / `Arrow` — nothing here renders or exports. `AlgebraicAgents` (0.4) is resolved as an ordinary registered dependency; there is no GitHub `[sources]` pin.

## The setting in one paragraph

The same pharma R&D pipeline the other demos use — `Discovery → Phase1 → Phase2 → Phase3 → Filed → Market`, one routing transition per phase boundary, each carrying a per-phase cycletime and probability-of-success — but taught through the refinement lens. To keep the refinement mechanics the star (and the demo fast), the phases are PLAIN counted species (a pool per phase), not the full structured-token `@select`/`@advance` machinery of `bd_acquisition`. The headline move: the coarse `flow_Phase2_Phase3` transition (a single atomic Phase-2 → Phase-3 advance) is `refine`d into a detailed four-step sub-model (screening → lead-opt → tox → filing, each with its own cycletime and PoS) — plug-compatibly, so the rest of the portfolio does not notice.

## What each section exercises

| § | Section | Capability exercised |
|---|---------|----------------------|
| 1 | The coarse portfolio | `@pipeline Name begin From => To : (ct=…, pos=…) … end` — expands a phase chain into N `flow`-genesis routing transitions (`flow_<From>_<To>`), each consuming its upstream phase as an upfront LHS (§2.8 token-flow) and carrying the per-edge (ct, pos) |
| 2 | Reusable fragments + ports | `@process name(params…) = begin … end` (a parameterized fragment factory, eval-free param substitution), `@port acs A => input B => output` (tag open-port roles via `=>` pairs), `@compose f1 f2` (`@join` + automatic output↔input port matching by FK-repoint: shared port collapses to ONE species, private species namespaced) |
| 3 | ★ refine | `refine(spec, transition, submodel; ports=Dict(boundary => sub_port, …))` — non-mutating splice of a finer sub-model into a coarse transition. Demonstrates PLUG-COMPATIBILITY (Invariant 1): boundary species keep their indices/names and every OTHER transition is structurally identical before/after |
| 4 | Advisory boundary check | `refinement_diagnostics(submodel, coarse_attrs; ports, tol)` — a `Vector{String}` of ADVISORY warnings: silent on a well-matched refinement, firing on Σ-cycletime / Π-PoS drift and on a dangling port (an `:input` never consumed) |
| 5 | Round-tripping the ladder | `abstract_transitions(spec, [subs…], :into; lhs, rhs, attrs)` — the inverse collapse back to one coarse transition; and `to_json_model` / `build_acs_from_dict` proving a refined spec reloads as a FLAT model (Invariant 5) |
| 6 | It's just a ModelSpec | The refined pipeline constructs (`ReactionNetworkProblem(...; seed=)`) + `simulate`s exactly like a hand-written flat model, reproducibly from (model, seed) (Invariant 4 closure) — 60 Discovery projects flow through the detailed Phase-2 sub-net to `:Market` |
| 7 | Recap | A closing summary tying the ladder to the maintainer's ask and to ADR 0009 / CONTRACT §11 |

## The refine payoff (§3) — the headline

`refine` splices a sub-model into a named coarse transition in four authoring-time structural moves: (1) namespace the sub's `:private` species; (2) identify the sub's open `:input`/`:output` ports with the parent's boundary species (per the `ports` map) by the ADR-0003 `ReactantSpec` FK-repoint — repoint an integer FK, no string surgery; (3) append the sub's transitions + remaining species/params/observables/events; (4) drop the coarse transition. Because move (2) leaves the BOUNDARY species (here `Phase2`, `Phase3`) at their same indices, names, and attributes, the coarse and refined models are **plug-compatible**: every transition NOT in the refined set is byte-for-byte structurally unchanged. The demo checks this explicitly — it prints the boundary indices (`Phase2 : 3 → 3`, `Phase3 : 4 → 4`) and confirms all four untouched transitions have identical structural signatures before and after. That is the multifidelity payoff: zoom the bottleneck, and the rest of the portfolio does not notice.

## The load-bearing points (narrated in-line)

These are surfaced in the script so a reader does not trip over them:

- **`@port` uses `=>` pairs, space-separated** (`@port acs A => input B => output`), NOT `=` (which macro-call syntax would parse as a keyword argument). The default role is `:private` (auto-namespaced); `:input`/`:output` are open ports matched by `@compose`; `:shared` is identified by bare name.
- **`@compose` identifies open ports by FK-repoint**, so a shared output→input port collapses to ONE species (the demo asserts `count(==(:Lead), names) == 1`), not two pools that happen to share a name — and it cannot corrupt a species name colliding inside a subexpression the way the old `recursively_substitute_vars!` path could.
- **`refine` is non-mutating** (`refine = refine!` on a `deepcopy`); the demo shows the coarse model is still intact after the splice. Sub-transitions come out namespaced `<coarse>__sub__<name>`; the sub's private species come out `<coarse>__sub__<species>`, so a bare private name never leaks.
- **`refinement_diagnostics` are ADVISORY (Invariant 6)** — warnings the author may override, not an equivalence proof. The linear-chain aggregate checks (Σct, ΠPoS) are best-effort and skipped if any needed attribute is non-numeric.
- **Macro arguments are LITERAL** — `@prob_init` / `@prob_meta` evaluate their right-hand sides in module scope, so a simulation config uses literal counts (a loop-local variable there fails with an `UndefVarError`); this matches the other demos' idiom. Reproducibility comes from the `seed=` construction kwarg alone.
- **The whole layer is AUTHORING-time and FORBIDDEN on a live/stepping model** — `refine`/`abstract`/`@compose` reindex (they drop the coarse `:T` and port `:S` rows), so they operate on a static `ReactionNetworkSchema`, never a stepping `ReactionNetworkProblem` (Invariant 3).

## Honest scope

- The demo uses PLAIN counted species (a pool per phase), not the structured-token `@select`/`@advance` idiom of `bd_acquisition`. This is a deliberate teaching choice — it keeps the refinement mechanics unobscured and the demo fast; the refinement layer itself is agnostic to which regime the species are in. (Refining a *structured* transition whose LHS selects tokens by an ADR-0008 predicate is an open question in ADR 0009 — the sub-model's input port would then be a predicated structured species.)
- `abstract_transitions` collapses the sub-transition rows and adds the coarse one; it does NOT garbage-collect the now-orphaned internal species rows (they remain, inert). The demo notes this in §5.
- Entity-level refinement — a structured token hosting its OWN sub-network — is explicitly scoped out of ADR-0009 v1 (a separate future ADR); this tour is process-structural refinement only.
- The §4 aggregate consistency checks are well-defined for a linear chain; for a sub-model with branches/loops "the" aggregate is not yet pinned down (an ADR-0009 open question).
