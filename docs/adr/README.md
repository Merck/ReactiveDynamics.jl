# Architecture Decision Records

One ADR per real decision in the ReactiveDynamics.jl rework. Each records the context, the decision, consequences, and the options considered and rejected. ADRs are append-only: to change a decision, add a new ADR that supersedes the old one rather than editing history.

| # | Title | Status |
|---|---|---|
| [0001](0001-discrete-event-engine.md) | Native discrete-event engine (keep, harden; not SciML) | Accepted 2026-06-17 |
| [0002](0002-priority-weighted-allocation.md) | Priority-weighted resource allocation under contention | Accepted 2026-06-17 |
| [0003](0003-data-store.md) | Data store: drop ACSets for a dep-free typed IR (phased) | Accepted 2026-06-18 |
| [0004](0004-runtime-mutation.md) | Runtime mutation contract (append-only, mutate-during-simulation) | Accepted 2026-06-18 |
| [0005](0005-serialization-json-ir.md) | Single JSON serialization + typed ExprNode IR (agentic artifact) | Accepted 2026-06-18 |
| [0006](0006-structured-tokens.md) | Structured/agentic tokens (live instantiation, queries) + the custom-function registry replacing `@register` | Accepted 2026-06-20 |
| [0007](0007-interface-and-initial-state.md) | Interface & initial-state contract: lifecycle phases, declarative initial marking, state dump/restore | Proposed 2026-06-21 |
| [0008](0008-token-filtration.md) | Agentic species under a filtration: predicate selection + state advance (`TokenPredicate`/`@select`, `SetField`/`@advance`); phase-as-attribute canonical | Proposed 2026-06-21 |
| [0009](0009-refinement-and-composition.md) | Hierarchical refinement, open-port composition, and compact process authoring | Proposed 2026-06-21 |
| [0010](0010-rules-and-conditional-transitions.md) | Rules/triggers & conditional transitions: the endogenous decision channel (repairs the no-op event channel) | Accepted 2026-06-21 |
| [0011](0011-action-callbacks-and-general-code.md) | Action callbacks: population-level token writes (`SetTokens`) + the general-code escape hatch (`Invoke`), eval-free | Accepted 2026-06-21 |

ADRs 0007–0011 are the Phase-0.5 modeling-language extension increment (2026-06-21), answering the maintainer's business-process modeling asks: a refined/compact/compositional metalanguage with substitutable granularity (0009), a contracted interface with declarative + dumpable initial state (0007), filtration-based selection of agentic-species tokens (0008), an endogenous, state-contingent decision channel — in-model rules/callbacks/conditional transitions — that makes an expressive BD demo possible (0010), and full action expressivity — population-level token-field writes plus a general-code escape hatch through the eval-free registry — so callbacks can run arbitrary queries/code without reopening the RCE (0011). They map to CONTRACT §10, §11, §9.5, §12, and the §12.3/§8.4 amendment respectively. 0007–0009 are PROPOSED pending sign-off — every `file:line` is read-verified against `ref-agents`, but unlike 0001–0006 they have not yet had an adversarial running-engine verification pass. 0010 and 0011 are maintainer-ratified (2026-06-21) and were driven by the [MVP_BD_DEMO.md](../MVP_BD_DEMO.md) exercise; 0010's repair of the event channel (`solvers.jl:323`) and 0011's reuse of the ADR-0006 §C registry boundary are verified against source.

Pending ADRs flagged by the rework brief (not yet written): the AlgebraicAgents integration surface (Workstream E — `@agentize`, `getobservable`/wire bridge for bidirectional observables); entity-level refinement (a structured token hosting its own sub-network — deferred from ADR 0009 §F). The typed `ExprNode` IR + modality re-modeling are specified across [ADR 0003](0003-data-store.md) (store + reactant promotion), [ADR 0005](0005-serialization-json-ir.md) (ExprNode), and [CONTRACT_DRAFT.md](../CONTRACT_DRAFT.md) §1 (modality truth table). A genesis-semantics decision (source vs routing transitions; the `{poisson,scheduled,flow,capacity}` mode tag) is recorded in CONTRACT_DRAFT.md §2/§3 and may be promoted to its own ADR. See [../../INVENTORY.md](../../INVENTORY.md) for the current-source map.
