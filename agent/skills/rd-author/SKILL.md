---
name: rd-author
description: Phase 2 of an rd-study — compile a decision brief into a ReactiveDynamics model document (.rdj.json), one file per experimental arm, and iterate to zero validation diagnostics. Covers places, transitions, arcs, the three modality axes, structured tokens, expression nodes, and the host.jl registry seam. Use after rd-frame, before rd-doctor.
---

# rd-author — compile the brief into a model document

The brief's slots table is your input; a validate-clean `.rdj.json` per arm is your output. **A model is data** — one JSON document, no Julia, closed node vocabulary. That is not a stylistic preference: the document is what makes the model reviewable, diffable, replayable, and safe to accept from an agent.

**Read `agent/reference/model_schema.md` now.** It is generated from the engine, so it is the authority on every key, column, operator and verb. This skill is the judgement layer on top of it; when the two disagree, the generated reference wins.

## The shape

```json
{
  "rd_format": "reactive-dynamics-model", "version": 1.0,
  "meta":        { "tspan": 48.0, "dt": 1.0 },
  "params":      [ { "name": "hit_rate", "value": 0.55 } ],
  "places":      [ { "name": "scientist", "init": 4.0 } ],
  "transitions": [ { "id": "trial", "rate": 1.5, "cycletime": 4.0 } ],
  "arcs":        [ { "transition": "trial", "place": "scientist", "side": "lhs" } ],
  "population":  [ ]
}
```

`meta.tspan` is required. `transitions[].id` is the FK every arc points at. Arcs are a flat list, not nested under transitions — which is what lets you diff one arc without touching anything else.

## Work in this order

1. **Places** — one per slot in the brief's resource and outcome rows. `init`, and `cost`/`reward` where value flows.
2. **Transitions** — one per process. `rate`, `cycletime`, `prob_of_success`, `priority`.
3. **Arcs** — the wiring, with modality. This is where the modeling lives and where the errors are.
4. **Population** — starting tokens, if the study has structured items.
5. **Validate, repair, repeat** until zero diagnostics.
6. **Arms** — copy the baseline, change the lever, nothing else.

## Modality: the axis that decides whether your model is right

Every lhs arc carries three independent axes. **Omitting `modality` means `consumed`.**

```json
{ "transition": "trial", "place": "scientist", "side": "lhs",
  "modality": { "allocation": "upfront", "return": "conserved", "blocking": "block" } }
```

| The resource is… | Write | Because |
|---|---|---|
| people, equipment, bench, a licence — **occupied then returned** | `"return": "conserved"` | the pool recovers when the cycle ends |
| cash, materials, a consumable — **spent** | `"return": "consumed"` (default) | it is gone |
| drawn continuously while work runs (budget/month, power) | `"allocation": "perstep"` | metered each tick over the cycle, not taken up front |
| merely *tested* — a gate that neither holds nor spends | `"blocking": "nonblock"` | a read arc |

Ask yourself of every single lhs arc: *when this finishes, does the pool go back up?* Getting `conserved` wrong on headcount produces a model that drains its own capacity to zero and never recovers — it will still validate, still run, and be wrong for the entire horizon. `rd-doctor`'s `modality_intent` check exists for exactly this, and it only works if the brief pre-registered `--expect-reusable`.

Illegal at construction: `{nonblock, conserved}`; `perstep` with `cycletime == 0`; `perstep` on a structured place.

## Expression nodes

Anywhere a number can vary, write a node instead of a literal. A bare number is accepted as shorthand for `const`.

```json
{ "node": "const", "value": 0.55, "symbol": false }
{ "node": "ref",   "kind": "param", "name": "hit_rate" }
{ "node": "call",  "op": "*", "args": [ {"node":"ref","kind":"place","name":"backlog"}, 0.1 ] }
{ "node": "sample","dist": "Normal", "args": [ 100.0, 15.0 ] }
```

`"symbol": true` on a `const` makes the value a **Symbol**, not a string — this is how you write `:Onc` or `:Active`. Getting it wrong means a predicate silently never matches, because `"Onc" != :Onc`. It is the second most common authoring error after modality.

Operators, distributions and reference kinds come from closed whitelists (see the generated reference). There is **no escape to arbitrary Julia** — anything outside the whitelist is a validation error, not a fallback. If a model seems to need one, the answer is a param, an observable, or a registry function, not a source string.

## Structured tokens

Use a structured place when each item carries state that drives decisions (a program, a phase, a value). Use a plain place when items are interchangeable — plain places are cheaper and simpler, so do not reach for tokens to model a queue of identical things.

```json
{ "name": "Project", "structured": true }
```

Three arc forms operate on tokens. None of them carries a `place` key — the `kind` is the place:

```json
// lhs: select a token to work on
{ "transition": "trial", "side": "lhs",
  "predicate": { "kind": "Project", "clauses": [
      ["program", "==", {"node":"const","value":"Onc","symbol":true}],
      ["phase",   "==", {"node":"const","value":"Active","symbol":true}] ] } }

// rhs: advance the bound token's lifecycle
{ "transition": "trial", "side": "rhs",
  "advance": { "field": "phase", "value": {"node":"const","value":"Launched","symbol":true} } }

// rhs: create a new token
{ "transition": "intake", "side": "rhs",
  "structured": { "kind": "Project", "fields": [
      { "name": "program", "value": {"node":"const","value":"Onc","symbol":true} },
      { "name": "value",   "value": {"node":"const","value":100.0} } ] } }
```

Starting tokens go in `population[]`, whose `attributes` take literals only — `":Onc"` denotes a Symbol:

```json
{ "place": "Project", "kind": "Project", "count": 3,
  "attributes": { "program": ":Onc", "phase": ":Active", "value": 100.0 } }
```

### `host.jl` — the one code seam

A token's Julia *type* is the only thing a document cannot carry. A study that uses tokens ships a small `host.jl` next to the model: the type, and a registry mapping kind names to constructors. Copy `agent/studies/marginal_scientist/host.jl` and change the fields. Two rules: keep it **tiny** — it declares a type, it does not contain modeling logic, which belongs in the document where it can be reviewed; and keep it **deterministic** — never call `rand` against the global RNG in it, because that breaks the `(model, seed)` replay guarantee the whole workflow rests on.

Pass it with `--host=path/to/host.jl` on every `rdx` command. Models without tokens need no host file.

## The validate loop

```bash
julia --project=. agent/scripts/rdx.jl validate agent/studies/<slug>/model.rdj.json --host=.../host.jl
```

Exit `0` clean, `1` diagnostics, `2` the document could not even be read. Every diagnostic carries a **JSON path** (`transitions[2].rate`, `arcs[7].modality`) — go to that path and fix it. Do not guess, and do not work around a diagnostic by deleting the construct that triggered it.

Reach zero diagnostics before moving on. `validate` is a static pass: it proves the document is *well-formed*, never that the model is *sensible*. That is phase 3's job, and phase 3 is not optional.

## Arms — the diff is the lever

One file per arm, in the same directory, differing in **as few lines as possible**:

```diff
  { "name": "scientist",
-   "init": 4.0 }
+   "init": 5.0 }
```

Each arm is a complete, self-contained, validate-clean document — not a patch. That costs some duplication and buys two things worth far more: the diff between two files *is* the intervention, reviewable at a glance; and each arm replays on its own forever. Name arms for the lever (`arm_hire_fifth.rdj.json`), and `diff` them before reporting to confirm nothing else drifted.

## Review with the user

Show the model back in the brief's vocabulary, not as JSON: a table of places with their roles, transitions with what they consume and produce, and the modality of each resource in words ("a trial *occupies* one scientist for four months and gives them back"). Then offer the `rdx draw` picture:

```bash
julia --project=. agent/scripts/rdx.jl draw <model> --out=figures/net.svg
```

The user is the only one who can catch a model that is correct but not about their system. Make that easy for them.
