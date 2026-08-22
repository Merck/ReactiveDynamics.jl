```@meta
CurrentModule = ReactiveDynamics
```

# Composition

The composition surface spans a granularity ladder: manual, no-port merges (`merge_networks!`/`@join`, `equalize!`/`@equalize`); declared open-port composition matched by foreign-key repoint (`@compose`/`@pipeline`/`@process`, `compose`, the `@port`/`set_port_role!`/`port_role` role surface); and boundary-matched refinement/abstraction (`refine`/`refine!`, `abstract!`/`abstract_transitions`), with `refinement_diagnostics` supplying advisory plug-compatibility checks. Composition and refinement are authoring-time-only operations — never apply them to a live, stepping model.

```@docs
merge_networks!
@join
equalize!
@equalize
refine
refine!
abstract_transitions
abstract!
set_port_role!
port_role
@port
@compose
@pipeline
@process
refinement_diagnostics
compose
```

## Namespacing internals

The merge operators namespace each fragment's private places before copying rows, so two fragments cannot conflate their state. These helpers are internal (unexported), documented here because the composition docstrings reference them.

```@docs
prepend!
prepend_obs!
normalize_name
```
