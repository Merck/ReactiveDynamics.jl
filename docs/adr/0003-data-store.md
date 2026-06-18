# ADR 0003: Data store for the ReactiveDynamics authoring/IR layer

Status: Proposed — pending maintainer confirmation

Date: 2026-06-18

Supersedes/relates to: ADR 0001 (keep the native discrete-event engine ReactionNetworkProblem), ADR 0002 (priority-weighted water-filling allocation). This ADR decides only the STATIC authoring/IR layer; the runtime engine and allocator are unaffected.

## Context

ReactiveDynamics is a timed, stochastic, resource-constrained generalized stochastic Petri net (GSPN) / discrete-event system. It currently uses ACSets.jl as the static authoring/IR layer: `using ACSets` (src/ReactiveDynamics.jl:3), a `BasicSchema` (ReactiveDynamics.jl:30-76), `@acset_type` (line 78), and the concrete alias `const ReactionNetworkSchema` (lines 80-88). We must decide whether to keep ACSets, extend it, or replace it, weighing the maintainer's priorities: production-readiness, dependency-minimalism, and agentic authoring (a typed, serializable, validatable spec an LLM can emit with NO eval of generated source — Workstream F).

Ground-truth facts, verified against source on branch ref-agents (not REVIEW.md):

- The schema has SIX objects (:S,:T,:E,:obs,:P,:M) and ZERO homs — the homs list is literally `[]` (src/ReactiveDynamics.jl:32). There is no categorical structure in the schema.
- There are ZERO categorical operations anywhere in src/ (no pushout/colimit/limit/migrate/oapply/ACSetTransformation/FinFunction). All ACSets usage is typed-columnar CRUD plus name-based linear reverse-lookup; `incident` is never a foreign-key follow but always equivalent to findall over an attribute column.
- The runtime hot path is already ACSets-free. ReactionNetworkProblem.getindex (src/state.jl:73-83) routes state[i,:attr] to plain Dict{Symbol,Vector} columns (state.attrs/state.transitions); only transPreAction/transPostAction fall through to state.acs. ACSets is solely the static authoring/IR layer; the compiled runtime store is already Dict-of-Vectors (state.jl:43,55).
- The transition↔reactant relation (which species each transition consumes/produces, with stoich + modality) is NOT stored as structure. It is re-parsed from the `trans` Expr at runtime by extract_reactants, called from src/state.jl:194 and src/solvers.jl:418.
- Composition is manual name-matching, not categorical. union_acs! (src/operators/joins.jl:14-26) loops incident(acs1, name, :specName) + add_part!; equalize! (src/operators/equalize.jl) merges species by name, with the only structural deletion at equalize.jl:52 and a fragile fallback at equalize.jl:55-63 that rewrites species names via recursively_substitute_vars! over EVERY attribute Expr — silently corrupting incidence if a species name collides inside a stoich/rate subexpression.
- ACSets is a LIGHT, stable dependency. Project.toml declares only `ACSets = \"0.2\"` (Project.toml:6,36); Catlab/GATlab/AlgebraicRewriting are not deps. ACSets has only ever shipped 0.1.x/0.2.x, so `compat \"0.2\"` cannot drift into a breaking 0.3.0. The \"Catlab version-churn bit-rot\" premise does NOT transfer to this dependency.
- There is an active RCE surface on import: load_network runs eval(Meta.parseall(...)) on prmVal (loadsave.jl:65) and on a `registered` source section (loadsave.jl:72), and the Base.convert hooks at ReactiveDynamics.jl:101-102 eval string attributes during assignment.
- Modality is an UNvalidated Set{Symbol} (ReactiveDynamics.jl:144); the legal set is {:nonblock,:conserved,:rate}.
- There is no Manifest.toml; dependency resolution floats.

## Options considered

- LENS A — Drop ACSets; hand-written dependency-free typed IR (struct of columnar object-tables) AND promote the transition↔reactant relation to a first-class typed ReactantSpec incidence table (integer FKs into S and T).
- LENS B — Drop ACSets; struct-of-Dict-columns IR with a signature-preserving shim; structure UNCHANGED (reactants stay parsed from the `trans` Expr).
- LENS C — Keep ACSets; add a :R Reactant object with real homs rspec:R→S, rtrans:R→T, turning the schema into a bipartite Petri net.
- LENS D — Dependency-free typed-struct IR as the canonical core; ship an optional to_acset/from_acset view behind a Julia package extension (weakdep on ACSets).

## Decision (proposed)

Replace ACSets with a dependency-free typed-struct-of-columns IR, executed as a phased migration, and reach LENS A's destination via LENS D's mechanism and LENS B's phasing:

- Phase 1: behavior-preserving store swap (LENS B core). A struct of typed columns per object, a single `const SCHEMA` source of truth replacing the propertynames(acs.subparts) reflection sites, and a signature-preserving shim implementing getindex/setindex!/nparts/parts/add_part!/add_parts!/incident/subpart/set_subpart!/rem_parts!. Runtime untouched (state.jl:73-83 already bypasses the store). Guarded by Phase-0 characterization tests.
- Phase 2: promote the transition↔reactant relation to a first-class typed ReactantSpec incidence table with integer FKs (LENS A's structural upgrade), with an Expr escape-hatch for the legitimately dynamic @choose/@move/@structured/expression-valued-stoich reactants. This makes specs FK-checkable for agentic authoring, makes the relation round-trippable in TOML/JSON, and makes species-merge in equalize! structurally exact instead of fragile string substitution.
- Phase 3 (optional): ship a tested to_acset/from_acset adapter behind a weakdep package extension (LENS D's interop hedge), preserving an optional AlgebraicPetri/Catlab view of the static spec at near-zero core cost.

If the maintainer chooses to stop after Phase 1, the terminal state is LENS B/D — a clean, conservative, dependency-free swap — which is itself fully defensible.

## Rationale

All three judge panels ranked the three drop-ACSets options (A, B, D) above the keep-ACSets option (C) on the priority axes, and independently grafted the same shape: phase the swap separately from the structural promotion, add the weakdep adapter as the interop hedge, and land the store-independent hardening fixes regardless. The panels split only on which lens to name first (two LENS A, one LENS D; B and D near-tied seconds) — a packaging disagreement, not a directional one. Phasing resolves it: Phase 1 is LENS B/D, Phase 2 reaches LENS A, Phase 3 adds LENS D's hedge.

The justification is FITNESS, not dependency weight. ACSets is light and stable and its entire categorical surface (homs, colimits, migrations) is unused here (zero homs at ReactiveDynamics.jl:32; zero categorical ops in src/). The wins are: an eval-free, JSON-Schema-emittable, validatable typed IR for agentic authoring; removal of an unused abstraction; and (Phase 2) promotion of the model's defining relation from a re-parsed Expr to inspectable, FK-checkable structure that also makes composition correct. Interop is impossible for the running GSPN engine (stateful in-flight Transition instances at state.jl:14-26, expression-valued rates) and off the BD/rNPV roadmap, so it is worth only an optional static-spec view — exactly Phase 3.

LENS C is not chosen: it keeps ACSets (no dependency-minimalism gain), does not close the eval-on-import RCE, and introduces a dual `trans`-Expr-vs-R-rows representation that can drift, while its strengths (compositionFit, interop) land on the least-weighted axes. Its one genuine insight — repointing a structural reference rather than string-substituting during species merge — is captured by Phase 2 without keeping ACSets or adding homs.

## Consequences

Positive: removes a dependency with zero functional loss; collapses store and IR into one inspectable typed struct an LLM can emit and a validator can check; enables eval-free import and modality enum validation; (Phase 2) makes the bipartite relation explicit, validatable, and round-trippable, and makes equalize!/union_acs! structurally exact; (Phase 3) keeps interop as a tested optional view at near-zero core cost.

Negative / risks: the migration touches a broad surface (~70 indexing sites plus nparts/parts/incident/add_part! sites) — uniform and mechanical, but large. The shim must faithfully reproduce three subtle invariants: subpart returning the LIVE column (broadcast-assign in update.jl), the acs[rowvec,:attr] slice (solvers.jl:557), and String→typed coercion (relocated to the loader, made eval-free). Phase 2 introduces genuinely new FK-remap/offset logic in union_acs!/equalize! and a parser split — the highest-regression area — which is why it is isolated as a separate PR behind Phase-0 tests. Phase 3 adds a second representation kept in sync by a round-trip test.

## Fixes required regardless of the chosen option

These are store-independent and should ride along with the store work:

1. Close the eval-on-import RCE: replace eval(Meta.parseall(...)) at loadsave.jl:65 and loadsave.jl:72 with parse-not-eval, and gate the `registered` source-injection path behind an explicit opt-in unsafe API (off by default).
2. Remove the eval-based Base.convert hooks at ReactiveDynamics.jl:101-102; parse modality and FoldedObservable structurally.
3. Validate modality ⊆ {:nonblock,:conserved,:rate} at construction/load time (currently unvalidated, ReactiveDynamics.jl:144), and reject the :conserved+:nonblock clash at validation time rather than deep in finish!.
4. Fix the broken free_blocked_species! (undefined `q`, solvers.jl ~512).
5. Commit a Manifest.toml for reproducible builds (currently absent).
6. Add Phase-0 characterization tests over the tutorial models before any store change.

## Open questions

See the accompanying review notes: scope/timing of Phase 2; whether to ship the Phase 3 weakdep adapter at all; backward compatibility for already-serialized models; whether the `registered` source-injection feature can be removed outright; stoichiometry typing; and whether to auto-derive a published JSON-Schema from `const SCHEMA` in this workstream.