```@meta
CurrentModule = ReactiveDynamics
```

# Structured tokens

Structured/agentic tokens are first-class entities: a token (a "project") carries typed attributes (its `phase`, `npv`, cost-to-date), custom behavior, and its own history, and can be instantiated, selected by predicate, and advanced through lifecycle phases. A structured species is declared with `@structured_token` and registered against a network before tokens are added to the initial marking.

```@docs
@structured_token
register_structured_species!
add_structured_token!
AbstractStructuredToken
BaseStructuredToken
PopulationEntry
log_token_fields
```
