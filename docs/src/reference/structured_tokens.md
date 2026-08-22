```@meta
CurrentModule = ReactiveDynamics
```

# Structured tokens

Structured/agentic tokens are first-class entities: a token (a "project") carries typed attributes (its `phase`, `npv`, cost-to-date), custom behavior, and its own history, and can be instantiated, selected by predicate, and advanced through lifecycle phases. A structured place — one whose tokens are distinguishable, a *colour set* in Coloured-Petri-net terms — is declared with `@structured_token` and registered against a network before tokens are added to the initial marking.

```@docs
@structured_token
register_token_kind!
add_structured_token!
AbstractStructuredToken
BaseStructuredToken
PopulationEntry
log_token_fields
```

## Runtime token instances

A `Transition` is a live in-flight instance of a transition recipe; the tokens it occupies are its bound structured agents. `ArcSpec` is the promoted transition↔place incidence row (an *arc*) — the typed, foreign-key-exact form of a reaction line that structured-token binding and place merges repoint against.

```@docs
Transition
ArcSpec
arcs
placename
```
