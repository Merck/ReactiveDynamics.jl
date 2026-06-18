# Architecture Decision Records

One ADR per real decision in the ReactiveDynamics.jl rework. Each records the context, the decision, consequences, and the options considered and rejected. ADRs are append-only: to change a decision, add a new ADR that supersedes the old one rather than editing history.

| # | Title | Status |
|---|---|---|
| [0001](0001-discrete-event-engine.md) | Native discrete-event engine (keep, harden; not SciML) | Accepted 2026-06-17 |
| [0002](0002-priority-weighted-allocation.md) | Priority-weighted resource allocation under contention | Accepted 2026-06-17 |
| [0003](0003-data-store.md) | Data store: drop ACSets for a dep-free typed IR (phased) | **Proposed** — pending maintainer confirmation |

Pending ADRs flagged by the rework brief (not yet written): the typed `AttrType`/`ExprNode` IR + modality re-modeling (REVIEW.md #10/#13; partly pre-decided by [ADR 0003](0003-data-store.md) Phase 2 and the modality truth table in [CONTRACT_DRAFT.md](../CONTRACT_DRAFT.md)), and the AlgebraicAgents integration surface (Workstream E — `@agentize`, `getobservable`/wire bridge). See [../../INVENTORY.md](../../INVENTORY.md) for the current-source map.
