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

Pending ADRs flagged by the rework brief (not yet written): the AlgebraicAgents integration surface (Workstream E — `@agentize`, `getobservable`/wire bridge for bidirectional observables). The typed `ExprNode` IR + modality re-modeling are now specified across [ADR 0003](0003-data-store.md) (store + reactant promotion), [ADR 0005](0005-serialization-json-ir.md) (ExprNode), and [CONTRACT_DRAFT.md](../CONTRACT_DRAFT.md) §1 (modality truth table). A genesis-semantics decision (source vs routing transitions; the `{poisson,scheduled,flow,capacity}` mode tag) is recorded in CONTRACT_DRAFT.md §2/§3 and may be promoted to its own ADR. See [../../INVENTORY.md](../../INVENTORY.md) for the current-source map.
