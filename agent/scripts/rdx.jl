#!/usr/bin/env julia
#
# rdx — the ReactiveDynamics agent CLI (AGENT_WORKFLOW_CHARTER W7).
#
# One entry point for every phase of the agentic workflow, designed to be driven by an agent over
# Bash: JSON on stdout, human notes on stderr, and an exit code that means something. It is a THIN
# transport over the engine's own functions — `validate`, `from_json_model`, `simulate`, `ensemble`,
# `summarize`, `treatment_effect`, `draw_network`/`exec_map` — and deliberately implements no
# semantics of its own. In particular the schema check CALLS `ReactiveDynamics.validate`; it never
# reimplements it (charter §5.4 architectural rule 1).
#
#   julia --project=. agent/scripts/rdx.jl <command> [args...] [--flag=value ...]
#
# Commands
#   validate MODEL                          static, eval-free schema check → diagnostics
#   run      MODEL --seed=N                  one seeded run → per-place/per-transition summary
#   ensemble MODEL --seeds=N --metric=M      seeded ensemble → mean ± sem for each metric
#   effect   BASE VARIANT --metric=M --seeds=N   unpaired treatment effect Δ ± se
#   diagnose MODEL --seed=N                  the §3.3 plausibility gate
#   draw     MODEL --out=PATH                the net (Layer A), or --exec-map over a run (Layer C)
#
# Exit codes: 0 = ok · 1 = findings (invalid model, failing checks, insignificant effect) · 2 = error.
# The 0/1 split is what lets an agent branch without parsing prose; `ok` in the envelope mirrors it.

using ReactiveDynamics
const RD = ReactiveDynamics
using AlgebraicAgents: simulate
using JSON
using Statistics

# ════════════════════════════════════════════════════════════════════════════════════════
# Argument parsing and the output envelope
# ════════════════════════════════════════════════════════════════════════════════════════

# Options that legitimately take several values. Repeating one of these accumulates; a comma list is
# equivalent, so `--tally=a --tally=b` and `--tally=a,b` mean the same thing. Getting this wrong is
# not a cosmetic bug: every one of these keys drives a CHECK, so a dropped value is a check that
# silently did not run while the envelope still reports `pass`.
const REPEATABLE = Set(["metric", "tally", "expect_binds", "expect_reusable", "trace"])

# Split ARGS into positionals and `--key=value` / `--flag` options. Keys keep their dashes turned
# into underscores so `--root-seed=1` reads as `opts["root_seed"]`.
function parse_args(argv)
    pos, opts = String[], Dict{String, String}()
    for a in argv
        if startswith(a, "--")
            body = a[3:end]
            k, v = occursin('=', body) ? split(body, '='; limit = 2) : (body, "true")
            key, val = replace(String(k), '-' => '_'), String(v)
            if !haskey(opts, key)
                opts[key] = val
            elseif key in REPEATABLE
                opts[key] *= "," * val
            elseif opts[key] != val
                # Last-wins on a scalar option is how you run with a seed you did not intend. Refuse
                # rather than silently pick one.
                error("rdx: --$(k) given twice with different values ($(opts[key]) and $(val)); it takes one value")
            end
        else
            push!(pos, a)
        end
    end
    return pos, opts
end

optf(opts, k, default) = haskey(opts, k) ? parse(Float64, opts[k]) : default
opti(opts, k, default) = haskey(opts, k) ? parse(Int, opts[k]) : default
# A repeatable option, as a list. Either spelling reaches here as one comma-joined string.
optlist(opts, k) = haskey(opts, k) ? String.(split(opts[k], ','; keepempty = false)) : String[]

# The uniform envelope every command emits. `findings` is the machine-readable "what is wrong" list
# an agent repairs against; `ok` is false exactly when the exit code is nonzero.
function emit(command, ok, result; findings = Any[], error = nothing)
    JSON.print(
        stdout,
        Dict(
            "command" => command, "ok" => ok, "result" => result,
            "findings" => findings, "error" => error,
        ),
        2,
    )
    println(stdout)
    return ok ? 0 : 1
end

note(msg) = println(stderr, msg)

# ════════════════════════════════════════════════════════════════════════════════════════
# Model loading — the host-code seam, made explicit
# ════════════════════════════════════════════════════════════════════════════════════════

# A study's optional `host.jl` is the ONE place the workflow executes authored code (ADR 0006 §B/§C;
# charter §3.2). It must define `RDX_REGISTRY::Dict{Symbol,Any}` mapping a token-kind name to a
# constructor; `PROJECT_REGISTRY` is accepted as an alias so the existing bd_acquisition host works
# unchanged. Absent a host file the registry is empty and the model must be registry-free — which is
# the data-only path the workflow prefers.
#
# Called exactly ONCE, from `main` before the command runs: including twice would redefine the token
# type, and the registry's constructor closures are born in a world newer than this script's methods,
# so the command itself is dispatched through `invokelatest` (see `main`).
function load_registry(opts)
    haskey(opts, "host") || return Dict{Symbol, Any}()
    path = abspath(opts["host"])
    isfile(path) || error("--host: no such file: $path")
    note("rdx: including host code from $path (executes authored Julia)")
    Base.include(Main, path)
    for name in (:RDX_REGISTRY, :PROJECT_REGISTRY)
        # `invokelatest`: the binding was created by the `include` above, i.e. in a world newer than
        # this method's. Julia 1.12 warns on the direct access and will error on it in future.
        isdefined(Main, name) && return Base.invokelatest(getfield, Main, name)
    end
    error("--host: $path defines neither RDX_REGISTRY nor PROJECT_REGISTRY")
end

read_model(path) = read(abspath(path), String)

# One declarative initial-marking attribute, JSON → the value `PopulationEntry` expects. Closed by
# design (ADR 0005): numbers, strings, and booleans pass through, a `":Sym"` string becomes the
# Symbol it names (the engine unwraps a `QuoteNode` without evaluating it), and anything else is
# rejected rather than guessed at. No JSON attribute can become executable Julia this way.
function population_attr(v, path)
    v isa Bool && return v
    v isa Real && return v
    if v isa AbstractString
        return startswith(v, ':') ? QuoteNode(Symbol(v[2:end])) : String(v)
    end
    return error(
        "$path: initial-marking attributes must be a number, string, \":Symbol\" or boolean; got $(typeof(v))"
    )
end

# `population[]` (ADR 0007 §B) read out of the DOCUMENT. `validate` already checks this key's
# well-formedness (`src/serialize.jl` rule 6), but `from_json_model` only accepts a `population=`
# kwarg and never reads `d["population"]` — so a model file cannot today declare its own starting
# tokens. Reading it here keeps a study's `model.rdj.json` self-contained (and hence its ensemble
# members and `reinit!` reproducible from the file alone) without changing the engine; closing the
# gap upstream is charter W14.
function document_population(doc)
    out = Any[]
    for (i, pe) in enumerate(get(doc, "population", []))
        place = Symbol(get(pe, "place", get(pe, "species", "")))
        attrs = Dict{Symbol, Any}(
            Symbol(f) => population_attr(v, "population[$i].attributes.$f")
                for (f, v) in get(pe, "attributes", Dict{String, Any}())
        )
        push!(
            out,
            RD.PopulationEntry(
                place, Symbol(pe["kind"]);
                count = Int(get(pe, "count", 1)), attributes = attrs,
            ),
        )
    end
    return out
end

# Parse → validate → build. Returns `(problem, diagnostics)`; the problem is `nothing` exactly when
# the document is invalid. `validate` gates construction inside `from_json_model` already — we call
# it first only to report diagnostics in our own envelope shape rather than as an exception.
function build_problem(path, opts, registry; seed)
    json = read_model(path)
    doc = JSON.parse(json)
    diags = safe_validate(doc, registry)
    isempty(diags) || return (nothing, diags)
    prob = model_builder(json, registry)(seed)
    return (prob, RD.Diagnostic[])
end

# The `seed -> problem` closure. One definition, so a single run and an ensemble member are built the
# same way — including the document-declared population.
function model_builder(json, registry)
    population = document_population(JSON.parse(json))
    return seed -> RD.from_json_model(
        json; seed = seed, registry = registry, population = population
    )
end

# The builder `ensemble` needs: build AND simulate. `ensemble` does not run its members itself unless
# given `max_t` ("`build` may already simulate", `src/analysis.jl`), and a build-only closure yields
# an ensemble of un-run models whose every metric reads zero. Simulating inside the closure also
# means the document's own `meta.tspan` is the default horizon, with `--tspan` as the override.
function simulating_builder(json, registry, tspan)
    build = model_builder(json, registry)
    return function (seed)
        prob = build(seed)
        tspan === nothing ? simulate(prob) : simulate(prob, tspan)
        return prob
    end
end

diag_dicts(diags) = [
    Dict("severity" => String(d.severity), "path" => d.path, "message" => d.msg) for d in diags
]

# `RD.validate` is a pure static pass, but it is not TOTAL over arbitrary JSON: a few reads are
# unguarded (e.g. `tr["id"]`, `src/serialize.jl:766`), so a first draft missing a required key throws
# instead of returning a diagnostic. An exception here would leave the agent an exit-2 stack trace
# and nothing to repair against, so a throw is converted into a diagnostic of its own. Making the
# validator total upstream is charter W15; this wrapper stays regardless, since a repair loop must
# never be able to dead-end on malformed input.
function safe_validate(doc, registry)
    try
        return RD.validate(doc; registry = registry)
    catch e
        return [
            RD.Diagnostic(
                :error, "<document>",
                "the validator could not read this document: $(sprint(showerror, e)). A required " *
                    "key is probably missing or has the wrong type — check that every transition has " *
                    "an `id`, every place a `name`, and every arc a `place`/`side`.",
            ),
        ]
    end
end

# ════════════════════════════════════════════════════════════════════════════════════════
# Reading a finished run
# ════════════════════════════════════════════════════════════════════════════════════════

place_names(prob) = String.(prob.network[:, :placeName])

# Map a transition's `transHash` — the key the log rows carry — to its authored name, so
# per-transition statistics come back in the vocabulary the model file uses. The hash lives on the
# COMPILED static table (`src/compilers.jl`: `coalesce(transName, gensym())`), not on the network, and
# the two are row-aligned; an unnamed transition falls back to its gensym so it is still countable.
function hash_to_name(prob)
    out = Dict{Symbol, String}()
    for (h, n) in zip(prob.transitions[:transHash], prob.network[:, :transName])
        out[Symbol(h)] = (n === nothing || ismissing(n)) ? String(Symbol(h)) : String(n)
    end
    return out
end

# Total spawns per transition, from the `(:new_transitions, t, (hash, q)...)` rows.
function spawn_counts(prob)
    names = hash_to_name(prob)
    out = Dict{String, Float64}(n => 0.0 for n in values(names))
    for row in prob.log
        row[1] === :new_transitions || continue
        for el in row[3:end]
            el isa Tuple || continue
            h, q = Symbol(el[1]), Float64(el[2])
            out[get(names, h, String(h))] = get(out, get(names, h, String(h)), 0.0) + q
        end
    end
    return out
end

# Total completions per transition (all, and the successful subset), from the `terminated_*` rows,
# whose payload is splatted `hash => count` pairs.
function completion_counts(prob, key::Symbol)
    names = hash_to_name(prob)
    out = Dict{String, Float64}(n => 0.0 for n in values(names))
    for row in prob.log
        row[1] === key || continue
        for el in row[3:end]
            el isa Pair || continue
            n = get(names, Symbol(el.first), String(Symbol(el.first)))
            out[n] = get(out, n, 0.0) + Float64(el.second)
        end
    end
    return out
end

# The per-tick fill fraction of in-flight instances (`:saturation`). A value below 1 means an
# instance did not get its full per-tick claim — i.e. a resource is genuinely contended. This is a
# sharper contention signal than a pool trough, because it fires even when a pool never hits zero.
function min_saturation(prob)
    m = 1.0
    for row in prob.log
        row[1] === :saturation || continue
        for el in row[3:end]
            el isa Tuple && (m = min(m, Float64(el[2])))
        end
    end
    return m
end

# Pool trough / peak / final per place, straight off `prob.sol` (the same signal `exec_map` paints
# with, `src/visualize.jl` `_pool_troughs`).
function pool_stats(prob)
    out = Dict{String, Any}()
    for s in place_names(prob)
        s in names(prob.sol) || continue
        col = prob.sol[!, s]
        isempty(col) && continue
        out[s] = Dict(
            "initial" => first(col), "final" => last(col),
            "min" => minimum(col), "max" => maximum(col),
        )
    end
    return out
end

# Summed ledger flows off the aggregate log rows.
function ledger_totals(prob)
    reward = cost = 0.0
    for row in prob.log
        row[1] === :valuation_reward && (reward += Float64(row[3]))
        row[1] === :valuation_cost && (cost += Float64(row[3]))
    end
    return (reward = reward, cost = cost)
end

function run_summary(prob)
    lt = ledger_totals(prob)
    return Dict(
        "t_end" => isempty(prob.sol) ? nothing : last(prob.sol.t),
        "ticks" => nrow_sol(prob),
        "places" => pool_stats(prob),
        "spawns" => spawn_counts(prob),
        "completions" => completion_counts(prob, :terminated_all),
        "successes" => completion_counts(prob, :terminated_success),
        "min_saturation" => min_saturation(prob),
        "ongoing_at_end" => length(prob.ongoing_firings),
        "ledger" => Dict("reward" => lt.reward, "cost" => lt.cost, "net" => lt.reward - lt.cost),
    )
end

nrow_sol(prob) = isempty(prob.sol) ? 0 : size(prob.sol, 1)

# The whole per-tick series for named places, not just the trough and the peak. `pool_stats` answers
# "did it run out"; this answers "for how long, and when" — which is the difference between a
# contention claim and a chart of one. Opt-in, because a 120-tick run over 11 places is a wall of
# numbers in an envelope meant to be read.
function pool_traces(prob, wanted)
    out = Dict{String, Any}("t" => collect(prob.sol.t))
    for s in wanted
        s in names(prob.sol) || error(
            "rdx run --trace: no place named $(repr(s)); this run has " *
                join(sort(place_names(prob)), ", ")
        )
        out[s] = collect(prob.sol[!, s])
    end
    return out
end

# ════════════════════════════════════════════════════════════════════════════════════════
# Named metrics — an ensemble reduction without user code
# ════════════════════════════════════════════════════════════════════════════════════════

# `summarize`/`treatment_effect` take a `metric(prob) -> Real` function. Letting the agent supply
# Julia here would reopen the code path the JSON artifact exists to close, so metrics are NAMED and
# resolved from a closed table instead. Unknown names are an error, not a silent zero.
#
#   final:<place> · min:<place> · max:<place> · sum:<place>   read prob.sol
#   spawns[:<transition>] · completions[:<t>] · successes[:<t>]   counts (total when unqualified)
#   reward · cost · net                                        summed ledger flows
#   enpv                                                       net discounted at --annual, × --scale
function resolve_metric(spec::AbstractString, opts)
    parts = split(spec, ':'; limit = 2)
    head = String(parts[1])
    arg = length(parts) == 2 ? String(parts[2]) : ""

    head == "final" && return p -> last(p.sol[!, arg])
    head == "min" && return p -> minimum(p.sol[!, arg])
    head == "max" && return p -> maximum(p.sol[!, arg])
    head == "sum" && return p -> sum(p.sol[!, arg])
    head == "reward" && return p -> ledger_totals(p).reward
    head == "cost" && return p -> ledger_totals(p).cost
    head == "net" && return p -> (lt = ledger_totals(p); lt.reward - lt.cost)

    if head in ("spawns", "completions", "successes")
        getter = head == "spawns" ? (p -> spawn_counts(p)) :
            head == "completions" ? (p -> completion_counts(p, :terminated_all)) :
            (p -> completion_counts(p, :terminated_success))
        return isempty(arg) ? (p -> sum(values(getter(p)))) : (p -> get(getter(p), arg, 0.0))
    end

    if head == "enpv"
        # Discounted net value, the case-study reduction (docs/literate/case_studies): per-tick
        # reward minus cost, discounted at an annual rate over ticks-as-months, then scaled into
        # the unit the study reports in. Both knobs are flags so the metric name stays data.
        annual = optf(opts, "annual", 0.08)
        scale = optf(opts, "scale", 1.0)
        periods = optf(opts, "periods_per_year", 12.0)
        return function (p)
            acc = 0.0
            for row in p.log
                flow = row[1] === :valuation_reward ? Float64(row[3]) :
                    row[1] === :valuation_cost ? -Float64(row[3]) : 0.0
                flow == 0.0 && continue
                acc += flow / (1 + annual)^(Float64(row[2]) / periods)
            end
            return acc * scale
        end
    end

    return error(
        "unknown metric $(repr(spec)); expected one of final:<place>, min:<place>, max:<place>, " *
            "sum:<place>, spawns[:<transition>], completions[:<t>], successes[:<t>], reward, cost, net, enpv"
    )
end

summary_dict(s) = Dict{String, Any}(
    "mean" => s.mean, "sem" => s.sem, "q25" => s.q25,
    "median" => s.median, "q75" => s.q75, "n" => s.n,
)

# The plausibility checks live in their own file because they are the part of this CLI that carries
# modeling judgement rather than transport. They read the run through the helpers above.
include(joinpath(@__DIR__, "diagnose.jl"))

# ════════════════════════════════════════════════════════════════════════════════════════
# Commands
# ════════════════════════════════════════════════════════════════════════════════════════

function cmd_validate(pos, opts, registry)
    length(pos) == 1 || error("usage: rdx validate MODEL")
    diags = safe_validate(JSON.parse(read_model(pos[1])), registry)
    errs = count(d -> d.severity === :error, diags)
    return emit(
        "validate", isempty(diags),
        Dict("diagnostic_count" => length(diags), "error_count" => errs);
        findings = diag_dicts(diags),
    )
end

function cmd_run(pos, opts, registry)
    length(pos) == 1 || error("usage: rdx run MODEL --seed=N")
    seed = opti(opts, "seed", 1)
    prob, diags = build_problem(pos[1], opts, registry; seed = seed)
    prob === nothing && return emit(
        "run", false, Dict("stage" => "validate"); findings = diag_dicts(diags)
    )

    haskey(opts, "tspan") ? simulate(prob, optf(opts, "tspan", 0.0)) : simulate(prob)
    result = run_summary(prob)
    result["seed"] = seed

    wanted = optlist(opts, "trace")
    isempty(wanted) || (result["trace"] = pool_traces(prob, wanted))

    if haskey(opts, "out")
        dir = abspath(opts["out"])
        RD.export_run(prob, dir)
        result["bundle"] = dir
        note("rdx: wrote run bundle to $dir")
    end
    return emit("run", true, result)
end

function cmd_ensemble(pos, opts, registry)
    length(pos) == 1 || error("usage: rdx ensemble MODEL --seeds=N --metric=M[,M...]")
    metrics = optlist(opts, "metric")
    isempty(metrics) && error("rdx ensemble: at least one --metric is required")
    nseed = opti(opts, "seeds", 24)
    root = opti(opts, "root_seed", 2026)
    json = read_model(pos[1])

    diags = safe_validate(JSON.parse(json), registry)
    isempty(diags) || return emit(
        "ensemble", false, Dict("stage" => "validate"); findings = diag_dicts(diags)
    )

    tspan = haskey(opts, "tspan") ? optf(opts, "tspan", 0.0) : nothing
    ens = RD.ensemble(
        simulating_builder(json, registry, tspan); nseed = nseed, root_seed = root
    )

    # `--per-seed` adds the raw member values next to the summary. A spread is the one thing a reader
    # has to see rather than be told: it is what makes an error bar credible instead of decorative.
    per_seed = haskey(opts, "per_seed")

    out = Dict{String, Any}()
    for m in metrics
        f = resolve_metric(m, opts)
        out[m] = summary_dict(RD.summarize(ens, f))
        per_seed && (out[m]["samples"] = [Float64(f(mem)) for mem in ens.members])
    end
    return emit(
        "ensemble", true,
        Dict("nseed" => nseed, "root_seed" => root, "metrics" => out),
    )
end

function cmd_effect(pos, opts, registry)
    length(pos) == 2 || error("usage: rdx effect BASELINE VARIANT --metric=M --seeds=N")
    metrics = optlist(opts, "metric")
    isempty(metrics) && error("rdx effect: at least one --metric is required")
    nseed = opti(opts, "seeds", 24)
    root = opti(opts, "root_seed", 2026)
    tspan = haskey(opts, "tspan") ? optf(opts, "tspan", 0.0) : nothing

    arms = Dict{String, Any}()
    for (label, path) in (("baseline", pos[1]), ("variant", pos[2]))
        json = read_model(path)
        diags = safe_validate(JSON.parse(json), registry)
        isempty(diags) || return emit(
            "effect", false, Dict("stage" => "validate", "arm" => label);
            findings = diag_dicts(diags)
        )
        arms[label] = RD.ensemble(
            simulating_builder(json, registry, tspan);
            nseed = nseed, root_seed = root, name = label,
        )
    end

    out, findings = Dict{String, Any}(), Any[]
    for m in metrics
        te = RD.treatment_effect(arms["baseline"], arms["variant"], resolve_metric(m, opts))
        # The charter's spread discipline (§3.3): an effect inside ±2 se is not a finding, and the
        # CLI says so in the envelope rather than leaving the agent to eyeball it.
        significant = te.se > 0 && abs(te.delta) > 2 * te.se
        out[m] = Dict(
            "delta" => te.delta, "se" => te.se, "baseline" => te.baseline,
            "variant" => te.deal, "n_baseline" => te.n_baseline, "n_variant" => te.n_deal,
            "significant_2se" => significant,
        )
        significant || push!(
            findings,
            Dict(
                "check" => "effect_vs_noise", "verdict" => "fail", "metric" => m,
                "message" => "Δ = $(round(te.delta; sigdigits = 4)) is within ±2 se " *
                    "($(round(2 * te.se; sigdigits = 4))) — not distinguishable from noise at n = $nseed. " *
                    "Raise the seed count or strengthen the lever before reporting this number.",
            )
        )
    end
    return emit(
        "effect", isempty(findings),
        Dict("nseed" => nseed, "root_seed" => root, "metrics" => out); findings = findings
    )
end

function cmd_draw(pos, opts, registry)
    length(pos) == 1 || error("usage: rdx draw MODEL --out=PATH [--exec-map --seed=N]")
    haskey(opts, "out") || error("rdx draw: --out=PATH is required")
    out = abspath(opts["out"])
    seed = opti(opts, "seed", 1)
    prob, diags = build_problem(pos[1], opts, registry; seed = seed)
    prob === nothing && return emit(
        "draw", false, Dict("stage" => "validate"); findings = diag_dicts(diags)
    )

    execmap = get(opts, "exec_map", "false") == "true"
    execmap && (haskey(opts, "tspan") ? simulate(prob, optf(opts, "tspan", 0.0)) : simulate(prob))

    # Charter §5.4 architectural rule 2: Graphviz is a SYSTEM binary, not a Julia package (the same
    # dependency behind the one @test_skip CI takes). A missing renderer degrades to a finding and a
    # DOT fallback — it never fails the workflow for a cosmetic reason.
    try
        execmap ? RD.exec_map(prob; path = out) : RD.draw_network(prob; path = out)
        return emit("draw", true, Dict("artifact" => out, "layer" => execmap ? "exec_map" : "structure"))
    catch e
        dot = out * ".dot"
        write(dot, RD.to_graphviz(RD.network_graph(prob)))
        return emit(
            "draw", true, Dict("artifact" => dot, "layer" => "dot_source", "rendered" => false);
            findings = [
                Dict(
                    "check" => "graphviz_available", "verdict" => "warn",
                    "message" => "Graphviz could not render ($(sprint(showerror, e))). Wrote DOT " *
                        "source to $dot instead; install a Graphviz backend for the diagram.",
                ),
            ]
        )
    end
end

function cmd_diagnose(pos, opts, registry)
    length(pos) == 1 || error("usage: rdx diagnose MODEL --seed=N")
    seed = opti(opts, "seed", 1)
    doc = JSON.parse(read_model(pos[1]))
    prob, diags = build_problem(pos[1], opts, registry; seed = seed)
    prob === nothing && return emit(
        "diagnose", false, Dict("stage" => "validate"); findings = diag_dicts(diags)
    )
    haskey(opts, "tspan") ? simulate(prob, optf(opts, "tspan", 0.0)) : simulate(prob)

    checks = RDXDiagnose.run_checks(prob, doc, opts)
    failed = filter(c -> c["verdict"] == "fail", checks)
    return emit(
        "diagnose", isempty(failed),
        Dict(
            "seed" => seed, "checks" => checks,
            "passed" => count(c -> c["verdict"] == "pass", checks),
            "warned" => count(c -> c["verdict"] == "warn", checks),
            "failed" => length(failed),
            "run" => run_summary(prob),
        );
        findings = failed,
    )
end

const COMMANDS = Dict(
    "validate" => cmd_validate, "run" => cmd_run, "ensemble" => cmd_ensemble,
    "effect" => cmd_effect, "diagnose" => cmd_diagnose, "draw" => cmd_draw,
)

function main(argv)
    isempty(argv) && return emit("", false, nothing; error = "usage: rdx <$(join(sort(collect(keys(COMMANDS))), '|'))> ...")
    cmd = argv[1]
    haskey(COMMANDS, cmd) || return emit(cmd, false, nothing; error = "unknown command $(repr(cmd))")
    try
        # Inside the try: a bad option list is a caller error like any other and must come back as the
        # envelope with exit 2, not as a Julia stack trace the calling agent cannot parse.
        pos, opts = parse_args(argv[2:end])
        # Host code loads ONCE, here, before any command runs — and the command is then dispatched
        # through `invokelatest`, because the registry's constructors were defined in a world newer
        # than these methods and a direct call would hit a `MethodError: method too new`.
        registry = load_registry(opts)
        return Base.invokelatest(COMMANDS[cmd], pos, opts, registry)
    catch e
        emit(cmd, false, nothing; error = sprint(showerror, e))
        return 2
    end
end

exit(main(ARGS))
