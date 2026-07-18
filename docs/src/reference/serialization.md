```@meta
CurrentModule = ReactiveDynamics
```

# Serialization

A model *is* data: a network round-trips through a single eval-free JSON document. Rate/attribute expressions are carried by the closed, typed `ExprNode` IR (never a captured host closure), host functions are referenced by name through a per-network registry, and load never `eval`s a model field — the typed IR plus validation is the trust boundary that closes the remote-code-execution surface. `from_json_model`/`to_json_model` are the file-level round-trip; the `*_to_dict`/`*_from_dict` and `build_network_from_dict` helpers expose the intermediate dict form; `@import_model`/`@export_model` are the authoring-macro wrappers.

The concrete JSON document shape is documented on the [JSON model schema](json_schema.md) page.

Loading routes through the (unexported) `validate` pass, which checks a candidate document against the closed IR and returns diagnostics rather than ever evaluating a field; call it as `ReactiveDynamics.validate(dict)` when inspecting a document directly.

```@docs
from_json_model
to_json_model
@import_model
@export_model
node_to_dict
node_from_dict
model_to_dict
build_network_from_dict
```

## The `ExprNode` IR

The closed, eval-free expression IR: every rate/attribute expression is one of these typed node kinds, with `to_expr`/`from_expr` converting between a node tree and a Julia `Expr`. `OP_WHITELIST`, `DIST_WHITELIST`, and `REF_KINDS` are the closed vocabularies the IR (and validation) admit.

```@docs
ExprNode
Const
NodeRef
Call
Sample
TimeRef
Choose
Field
ExternalRef
to_expr
from_expr
OP_WHITELIST
DIST_WHITELIST
REF_KINDS
```
