# Architecture Decision Records

One ADR per real decision in the ReactiveDynamics.jl rework. Each records the context, the decision, consequences, and the options considered and rejected. ADRs are append-only: to change a decision, add a new ADR that supersedes the old one rather than editing history.

| # | Title | Status |
|---|---|---|
| [0001](0001-discrete-event-engine.md) | Native discrete-event engine (keep, harden; not SciML) | Accepted 2026-06-17 |
| [0002](0002-priority-weighted-allocation.md) | Priority-weighted resource allocation under contention | Accepted 2026-06-17 |

Pending ADRs flagged by the rework brief (not yet written): the ACSets/Catlab data-store fork (Workstream C), the typed `AttrType`/IR + modality re-modeling (REVIEW.md #10/#13), and serialization safety (close the eval-on-import RCE, Workstream F). See [../../INVENTORY.md](../../INVENTORY.md) for the current-source map and the [rework brief] for the full fork list.
