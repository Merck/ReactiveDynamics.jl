#!/usr/bin/env julia
#
# gen_reference.jl — generate `agent/reference/model_schema.md` FROM THE SOURCE (charter W6, Tier 2).
#
#   julia --project=. agent/scripts/gen_reference.jl            # write the extract
#   julia --project=. agent/scripts/gen_reference.jl --check    # fail if it is stale
#
# The skills need to know the authoring vocabulary exactly: which columns exist, which operators and
# distributions are on the whitelist, which action verbs and node tags a document may carry. That
# information already exists in `src/` as the object model and the whitelist tuples, so restating it
# by hand would create a second source of truth that silently rots — the failure the docs charter's
# anti-fork rule (§D4) exists to prevent. This script reads it out of the loaded module instead, and
# `--check` is the CI gate that the committed extract still matches.
#
# Hand-written PROSE about how to author a model belongs in the skills (Tier 3). This file carries
# only what can be derived.

using ReactiveDynamics
const RD = ReactiveDynamics
using JSON

const OUT = joinpath(@__DIR__, "..", "reference", "model_schema.md")

# The store's column families, in the object-model order a reader cares about.
const OBJECTS = [
    (:S, "places[]", "one resource pool or structured-token container per row"),
    (:T, "transitions[]", "one static transition recipe per row"),
    (:P, "params[]", "one model parameter per row"),
    (:obs, "observables[]", "one named computed value per row"),
    (:M, "meta{}", "solver settings, as an object not an array — every key is passed through as a constructor keyword"),
]

# JSON key ⇒ (store column, meaning). The mapping itself lives in `src/serialize.jl`'s read/emit
# pairs (`build_network_from_dict`'s jsonkey→col tuples and `_place_to_dict`/`_transition_to_dict`),
# which are function-local and so cannot be read out of the module — this table is therefore the one
# HAND-MAINTAINED part of the extract. `check_facts` below guards it: every column named here must
# exist in `SCHEMA`, and every JSON key must appear as a string literal in the serializer. A rename
# on either side breaks the generator rather than silently producing a lying reference.
const NO_COLUMN = "—"   # a document-only key, with no store column behind it
const JSON_KEYS = [
    ("meta.tspan", "metaVal", "REQUIRED — the horizon, in ticks"),
    ("meta.dt", "metaVal", "tick length (default `tunit`)"),
    ("meta.tunit", "metaVal", "the time unit `tspan`/`dt` are expressed in"),
    ("params[].name", "prmName", "parameter name; target of a `ref` node with `kind: \"param\"`"),
    ("params[].value", "prmVal", "a JSON number, never an expression"),
    ("observables[].name", "obsName", "observable name; target of a `ref` node with `kind: \"obs\"`"),
    ("places[].name", "placeName", "the place's name; referenced by arcs and by `ref` nodes"),
    ("places[].init", "placeInitVal", "initial marking (omit for 0)"),
    ("places[].structured", "placeStructured", "`true` for a token container"),
    ("places[].cost", "placeCost", "per-unit ledger cost charged when the place is drawn"),
    ("places[].reward", "placeReward", "per-unit ledger reward credited when the place is PRODUCED"),
    ("places[].valuation", "placeValuation", "per-unit carrying valuation"),
    ("places[].modality", "placeDefaultModality", "default modality for arcs on this place that omit one"),
    ("transitions[].id", "transName", "REQUIRED — the FK arcs point at"),
    ("transitions[].name", "transName", "display name (defaults to `id`)"),
    ("transitions[].rate", "transRate", "spawn rate — a node or a bare number"),
    ("transitions[].rate_mode", NO_COLUMN, "`\"poisson\"` (default) or `\"deterministic\"` for a fixed count"),
    ("transitions[].cycletime", "transCycleTime", "ticks an instance stays in flight"),
    ("transitions[].prob_of_success", "transProbOfSuccess", "terminal success probability, in [0,1]"),
    ("transitions[].priority", "transPriority", "allocator weight under contention (ADR 0002)"),
    ("transitions[].capacity", "transCapacity", "max simultaneous in-flight instances"),
    ("transitions[].multiplier", "transMultiplier", "scales every arc's multiplicity"),
    ("transitions[].max_lifetime", "transMaxLifeTime", "abandon an instance after this long"),
    ("transitions[].pre_action", "transPreAction", "action fired when an instance spawns — **NOT read on import today**: export writes it, the loader drops it without a diagnostic (charter W30); use `rules[]`"),
    ("transitions[].post_action", "transPostAction", "action fired when an instance completes — **NOT read on import today**, as `pre_action` (charter W30)"),
    ("arcs[].transition", "ArcSpec.trans", "REQUIRED — the owning transition's `id`"),
    ("arcs[].side", "ArcSpec.side", "REQUIRED — `\"lhs\"` (input) or `\"rhs\"` (output)"),
    ("arcs[].place", "ArcSpec.place", "the place, for a static arc"),
    ("arcs[].multiplicity", "ArcSpec.multiplicity", "coefficient (default 1)"),
    ("arcs[].modality", "ArcSpec.modality", "`{allocation, return, blocking}` — see below"),
    ("arcs[].predicate", NO_COLUMN, "lhs token selection: `{kind, clauses:[[field,op,node]]}` (a `@select`) — no `place` key; the `kind` IS the place"),
    ("arcs[].advance", NO_COLUMN, "rhs field write on the bound token: `{field, value:node}` (an `@advance`) — no `place` key"),
    ("arcs[].structured", NO_COLUMN, "rhs token genesis: `{kind, fields:[{name,value:node}]}` (a `@structured`) — no `place` key"),
    ("population[].place", NO_COLUMN, "structured place the starting tokens go on"),
    ("population[].kind", NO_COLUMN, "registry key resolved to a host constructor"),
    ("population[].count", NO_COLUMN, "how many to instantiate (default 1) — applied by `rdx`, see below"),
    ("population[].attributes", NO_COLUMN, "field ⇒ literal; `\":Sym\"` denotes a Symbol — applied by `rdx`, see below"),
]

# Store columns a JSON document CANNOT carry. Worth naming explicitly: an authoring agent that sees
# them in the column dump above would otherwise invent a key for them.
const NOT_IN_DOCUMENT = [
    ("placeRole", "authoring-time role tag, set by the DSL"),
    ("placeInitUncertainty", "initial-marking uncertainty, set by the DSL"),
    ("trans", "the raw reaction line — DSL provenance, reconstructed on import"),
]

# The three modality axes and their legal values, with the DSL marker each triple corresponds to.
const MODALITIES = [
    ("allocation", "`upfront` (default) · `perstep`", "`perstep` is the `@rate` marker — metered each tick over the cycle"),
    ("return", "`consumed` (DEFAULT) · `conserved`", "`conserved` is the `@conserved` marker — held for the cycle, returned intact"),
    ("blocking", "`block` (default) · `nonblock`", "`nonblock` is the `@nonblock` marker — a read/test arc that neither holds nor spends"),
]

# The JSON shape of each action verb — the `action` of a `rules[]` entry — as `stmt_from_dict` reads
# it. Hand-maintained like JSON_KEYS (the shapes live in `src/serialize.jl`'s `stmt_from_dict`
# branches), and guarded the same way: `check_facts` asserts the verbs here are exactly
# `RD.ACTION_VERBS`, and that every example parses with `stmt_from_dict` and serializes back through
# `stmt_to_dict` to the same keys. So an example can be pasted into a document as is.
const ACTION_SHAPES = [
    ("set_marking", """{"verb": "set_marking", "name": "A", "value": {"node": "const", "value": 0}, "mode": "set"}""", "write place `name`'s marking; `\"mode\": \"inc\"` adds `value` instead (default `set`)"),
    ("set_params", """{"verb": "set_params", "assigns": [{"name": "beta", "value": {"node": "const", "value": 0.5}}]}""", "write one or more params"),
    ("set_field", """{"verb": "set_field", "field": "phase", "value": {"node": "const", "value": 2}}""", "write a field of the BOUND token — illegal in a rule, which has none; a document writes token fields with `set_tokens` or an `arcs[].advance`"),
    ("set_tokens", """{"verb": "set_tokens", "predicate": {"kind": "Project", "clauses": [["phase", "==", {"node": "const", "value": 1}]]}, "assigns": [{"name": "phase", "value": {"node": "const", "value": 2}}]}""", "write fields on every token the predicate selects"),
    ("add_token", """{"verb": "add_token", "kind": "Project", "fields": [{"name": "phase", "value": {"node": "const", "value": 1}}]}""", "create a token; `kind` resolves through the registry"),
    ("activate", """{"verb": "activate", "transition": "t1"}""", "re-open a transition's spawn gate"),
    ("deactivate", """{"verb": "deactivate", "transition": "t1"}""", "close it; instances already in flight run to completion"),
    ("invoke", """{"verb": "invoke", "fn": "notify", "args": [{"node": "const", "value": 1}]}""", "call a registry host function by name"),
    ("log", """{"verb": "log", "msg": "budget exhausted"}""", "append a message (a string or a node) to the run log"),
    ("seq", """{"verb": "seq", "stmts": [{"verb": "deactivate", "transition": "t1"}, {"verb": "log", "msg": "t1 closed"}]}""", "run statements in order"),
]

# Guard the hand-maintained table against drift on either side (see JSON_KEYS). Throws rather than
# emitting a reference the engine would contradict — a wrong reference is worse than a missing one,
# because the agent reading it has no way to notice.
function check_facts()
    all_cols = Set(String(c) for obj in values(RD.SCHEMA) for c in keys(obj))
    # A key must be READ by something on the verified path. `src/serialize.jl` is the engine's reader;
    # `rdx.jl` is the reader for the `population[]` keys the engine does not yet consume itself
    # (charter W14). Both count — what must not happen is documenting a key nothing reads.
    readers = join(
        read(joinpath(@__DIR__, "..", "..", p), String)
            for p in ("src/serialize.jl", "agent/scripts/rdx.jl")
    )
    problems = String[]

    for (json, col, _) in JSON_KEYS
        # `ArcSpec.x` names a struct field, not a store column; check those against the struct.
        if startswith(col, "ArcSpec.")
            f = Symbol(col[(length("ArcSpec.") + 1):end])
            f in fieldnames(RD.ArcSpec) || push!(problems, "$json: ArcSpec has no field `$f`")
        elseif col != NO_COLUMN
            col in all_cols || push!(problems, "$json: `$col` is not a column in SCHEMA")
        end
        # `meta` is the one open key set: `_meta_kwargs` passes EVERY key straight through as a
        # constructor keyword, so a meta key need not appear as a literal anywhere. The check's
        # closed-vocabulary assumption does not hold there — exempt it rather than weaken the rule
        # for everything else. (`meta.tspan` is still enforced: `from_json_model` requires it.)
        startswith(json, "meta.") && continue
        leaf = last(split(json, '.'))
        occursin("\"$leaf\"", readers) ||
            push!(problems, "$json: key `\"$leaf\"` is read by neither src/serialize.jl nor rdx.jl")
    end
    for (col, _) in NOT_IN_DOCUMENT
        col in all_cols || push!(problems, "NOT_IN_DOCUMENT: `$col` is not a column in SCHEMA")
    end
    for (sym, _, _) in OBJECTS
        haskey(RD.SCHEMA, sym) || push!(problems, "OBJECTS: SCHEMA has no object `:$sym`")
    end
    verbs = [Symbol(v) for (v, _, _) in ACTION_SHAPES]
    Set(verbs) == Set(RD.ACTION_VERBS) && length(verbs) == length(RD.ACTION_VERBS) || push!(
        problems,
        "ACTION_SHAPES: verbs $(sort(verbs)) are not exactly RD.ACTION_VERBS $(sort(collect(RD.ACTION_VERBS)))",
    )
    for (verb, example, _) in ACTION_SHAPES
        d = JSON.parse(example)
        try
            back = RD.stmt_to_dict(RD.stmt_from_dict(d))
            Set(keys(back)) == Set(keys(d)) ||
                push!(problems, "ACTION_SHAPES: `$verb` keys $(sort(collect(keys(d)))) ≠ serialized $(sort(collect(keys(back))))")
        catch e
            push!(problems, "ACTION_SHAPES: the `$verb` example does not parse: $(sprint(showerror, e))")
        end
    end

    isempty(problems) && return nothing
    error(
        "gen_reference: the hand-maintained JSON↔column table is stale:\n  " *
            join(problems, "\n  ") *
            "\nFix agent/scripts/gen_reference.jl to match src/, then regenerate."
    )
end

function render()
    check_facts()
    io = IOBuffer()
    println(io, "<!-- GENERATED by agent/scripts/gen_reference.jl — do not edit by hand. -->")
    println(io, "<!-- Regenerate: julia --project=. agent/scripts/gen_reference.jl -->")
    println(io)
    println(io, "# Model-document reference (generated)")
    println(io)
    println(io, "The authoring vocabulary of a `.rdj.json` model document, extracted from the loaded `ReactiveDynamics` module so it cannot drift from the engine. The normative semantics live in `spec/CONTRACT_DRAFT.md`; this is the surface.")
    println(io)
    println(io, "A model is ONE JSON object. `meta.tspan` is required; everything else is optional and defaulted.")
    println(io)

    println(io, "## Top-level keys")
    println(io)
    println(io, "`rd_format` · `version` · `meta` · `params` · `places` · `transitions` · `arcs` · `observables` · `population` · `rules` · `inputs`")
    println(io)

    println(io, "## Object model")
    println(io)
    println(io, "| Object | JSON | Rows are |")
    println(io, "|---|---|---|")
    for (sym, arr, what) in OBJECTS
        println(io, "| `:$sym` | `$arr` | $what |")
    end
    println(io)
    println(io, "Column names per object, as the store declares them (`RD.SCHEMA`):")
    println(io)
    for (sym, arr, _) in OBJECTS
        cols = sort(collect(keys(RD.SCHEMA[sym])))
        isempty(cols) && continue
        println(io, "- **`:$sym`** — ", join(("`$c`" for c in cols), ", "))
    end
    println(io)

    println(io, "## JSON keys")
    println(io)
    println(io, "| JSON path | Store column | Meaning |")
    println(io, "|---|---|---|")
    for (json, col, meaning) in JSON_KEYS
        println(io, "| `$json` | `$col` | $meaning |")
    end
    println(io)
    println(io, "`meta` is the one OPEN key set: every key is passed through as a constructor keyword, so it accepts more than the rows above (and silently ignores what the constructor does not use — `alloc_strategy` is accepted-and-ignored since ADR 0002 made the priority-weighted allocator unconditional). Everything else is a closed vocabulary.")
    println(io)
    println(io, "Columns with NO document key — set only by the Julia DSL, so do not invent a JSON spelling for them:")
    println(io)
    for (col, why) in NOT_IN_DOCUMENT
        println(io, "- `$col` — $why")
    end
    println(io)

    println(io, "## Arc modality — the three axes")
    println(io)
    println(io, "An input arc's `modality` is an object with three independent axes. **Omitting it means `consumed`**, which is the single most common modeling error: a reusable resource (headcount, equipment, bench capacity) declared by default is *spent* rather than occupied, so its pool drains and never recovers.")
    println(io)
    println(io, "| Axis | Values | Note |")
    println(io, "|---|---|---|")
    for (axis, values, note) in MODALITIES
        println(io, "| `$axis` | $values | $note |")
    end
    println(io)
    println(io, "Illegal combinations are rejected at construction (CONTRACT §1.4): `{nonblock, conserved}`; `perstep` with `cycletime == 0`; `perstep` on a structured place.")
    println(io)

    println(io, "## Expression nodes (the closed IR)")
    println(io)
    println(io, "Every time-varying attribute is a node tree, never a source string — this is the eval-free trust boundary (ADR 0005). A bare number is accepted as shorthand for `const`.")
    println(io)
    println(io, "| Tag | Shape |")
    println(io, "|---|---|")
    println(io, "| `const` | `{\"node\":\"const\",\"value\":v,\"symbol\":false}` — `symbol: true` makes the value a Symbol name |")
    println(io, "| `ref` | `{\"node\":\"ref\",\"kind\":\"place\"\\|\"param\"\\|\"obs\",\"name\":\"n\"}` |")
    println(io, "| `call` | `{\"node\":\"call\",\"op\":\"+\",\"args\":[…]}` |")
    println(io, "| `sample` | `{\"node\":\"sample\",\"dist\":\"Poisson\",\"args\":[…]}` |")
    println(io, "| `timeref` | `{\"node\":\"timeref\"}` — the current tick |")
    println(io, "| `choose` | `{\"node\":\"choose\",\"alts\":[{\"weight\":w,\"value\":node}]}` |")
    println(io, "| `field` | `{\"node\":\"field\",\"name\":\"f\"}` — a bound token's field |")
    println(io, "| `externalref` | `{\"node\":\"externalref\",\"port\":\"p\"}` — a declared AA input port |")
    println(io)
    println(io, "**Operator whitelist** (`RD.OP_WHITELIST`): ", join(("`$op`" for op in RD.OP_WHITELIST), " · "))
    println(io)
    println(io, "**Distribution whitelist** (`RD.DIST_WHITELIST`): ", join(("`$d`" for d in RD.DIST_WHITELIST), " · "))
    println(io)
    println(io, "**Reference kinds** (`RD.REF_KINDS`): ", join(("`$k`" for k in RD.REF_KINDS), " · "))
    println(io)
    println(io, "Anything outside these tuples is a validation error, not a fallback. There is no escape to arbitrary Julia.")
    println(io)

    println(io, "## Token selection and lifecycle")
    println(io)
    println(io, "**Predicate clause operators** (`RD.PRED_OP_WHITELIST`): ", join(("`$op`" for op in RD.PRED_OP_WHITELIST), " · "))
    println(io)
    println(io, "A clause is a 3-element array `[field, op, node]`. A predicate is `{\"kind\":\"K\",\"clauses\":[…]}` and selects tokens of kind `K` satisfying every clause.")
    println(io)
    println(io, "### Starting population — read by `rdx`, not by the engine")
    println(io)
    println(io, "`population[]` declares the t=0 tokens. `validate` checks it (the `place` must be a declared structured place; the `kind` must be in the registry), but `from_json_model` does NOT instantiate it — it constructs the network and drops the array. `rdx` closes the gap: it reads `population[]` itself and passes `PopulationEntry(place, kind; count, attributes)` to the constructor. The practical consequence: **a document's starting tokens materialize under `rdx`, and silently do not under a bare `from_json_model` call.** (Charter W14 — the fix belongs upstream in the engine.)")
    println(io)
    println(io, "`attributes` values are literals only, on a closed whitelist: numbers, strings, booleans, and `\":Name\"` for a Symbol. Anything else is rejected — an attribute is not an expression.")
    println(io)

    println(io, "## Rule action verbs")
    println(io)
    println(io, "The endogenous decision channel (ADR 0010): a `rules[]` entry is `{id, guard, action, fire_mode}`. Action verbs (`RD.ACTION_VERBS`):")
    println(io)
    println(io, join(("`$v`" for v in RD.ACTION_VERBS), " · "))
    println(io)
    println(io, "`guard` is a node evaluated once per tick: a boolean fires the action once when true, a number `v` fires it `Poisson(v)` times. `fire_mode` is `\"every_tick\"` (default) or `\"once\"`, which disables the rule after it first fires. A `value` is a node or a bare literal. Each example below parses as is (the generator checks it).")
    println(io)
    println(io, "| Verb | JSON | Does |")
    println(io, "|---|---|---|")
    for (verb, example, meaning) in ACTION_SHAPES
        println(io, "| `$verb` | `$example` | $meaning |")
    end
    println(io)
    println(io, "`invoke` and a `add_token` kind are resolved through the per-network registry BY NAME — the host-code seam, never inline source.")
    println(io)

    println(io, "## Validation")
    println(io)
    println(io, "`ReactiveDynamics.validate(dict; registry) -> Vector{Diagnostic}` is a pure, eval-free static pass, and `from_json_model` refuses to construct unless it returns empty. Each `Diagnostic` carries a `severity`, a JSON `path` and a message — the path is what makes a failure machine-repairable. Run it (via `rdx validate`) before every simulate.")
    return String(take!(io))
end

function main(argv)
    text = render()
    if "--check" in argv
        isfile(OUT) || (println(stderr, "gen_reference: $OUT is missing — run without --check"); return 1)
        if read(OUT, String) == text
            println("gen_reference: $(relpath(OUT)) is up to date")
            return 0
        end
        println(stderr, "gen_reference: $(relpath(OUT)) is STALE — regenerate it (the engine's schema/whitelists moved)")
        return 1
    end
    mkpath(dirname(OUT))
    write(OUT, text)
    println("gen_reference: wrote $(relpath(OUT))")
    return 0
end

exit(main(ARGS))
