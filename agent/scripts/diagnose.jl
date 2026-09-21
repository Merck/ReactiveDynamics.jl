#
# diagnose.jl — the plausibility gate (AGENT_WORKFLOW_CHARTER §3.3, tracker W4).
#
# `ReactiveDynamics.validate` proves a model document is WELL-FORMED. It cannot prove the model means
# anything: a schema-clean net whose transitions never fire, whose resources are never contended, or
# whose horizon truncates the answer is exactly the plausible-but-wrong artifact an LLM produces by
# default. These checks close that gap, and they are the reason the workflow may quote a number at
# all.
#
# Every check is derived from a finished seeded run plus the model document — no new engine code:
#   · `prob.sol`             per-place trajectory (the pool troughs `exec_map` paints with)
#   · `prob.log`             `:new_transitions` spawns, `:saturation` fill fractions,
#                            `:terminated_*` completions, `:valuation_*` ledger flows
#   · `prob.ongoing_firings` instances still in flight when the horizon ended
#   · the JSON document      static intent: declared modalities, success probabilities
#
# Verdicts are `pass` / `warn` / `fail`. Only `fail` blocks: it means a number read off this run would
# answer a different question than the one asked. `warn` is a caveat the briefing must carry.
module RDXDiagnose

using Statistics

# The check record every function below returns. `evidence` is deliberately structured rather than
# prose, so an agent can repair against it and a report can quote it.
mk(check, verdict, message, evidence = Dict{String, Any}()) = Dict{String, Any}(
    "check" => check, "verdict" => verdict, "message" => message, "evidence" => evidence
)

# ── helpers over the finished run ───────────────────────────────────────────────────────

# `Main` holds the reader helpers defined by rdx.jl; diagnose consumes them rather than duplicating
# the log-row decoding in two places.
spawns(prob) = Main.spawn_counts(prob)
successes(prob) = Main.completion_counts(prob, :terminated_success)
pools(prob) = Main.pool_stats(prob)
saturation(prob) = Main.min_saturation(prob)

# Places whose pool ran to (or below) zero at its trough — the classic starvation signal. A place
# that is empty for the WHOLE run is not starved, it is unused: requiring `max > 0` keeps an
# never-populated structured place from reading as the binding constraint.
starved_places(prob) = String[
    s for (s, st) in pools(prob) if st["min"] <= 0.0 && st["max"] > 0.0
]

# Mean of the last `frac` of a column, and of the `frac` window before it — the convergence probe.
function tail_windows(col, frac = 0.25)
    n = length(col)
    n < 8 && return (nothing, nothing)
    w = max(2, floor(Int, n * frac))
    return (mean(col[(n - w + 1):n]), mean(col[(n - 2w + 1):(n - w)]))
end

# ── the checks ──────────────────────────────────────────────────────────────────────────

# 1. A transition that never spawned is a dead branch of the model — its rate is too low for the
#    horizon, its inputs are never satisfiable, or it is gated off. Either way a whole mechanism the
#    author believed was in the model is absent from the result.
function check_transitions_fire(prob, doc, opts)
    sp = spawns(prob)
    dead = sort(String[t for (t, q) in sp if q <= 0.0])
    isempty(dead) && return mk(
        "transitions_fire", "pass",
        "All $(length(sp)) transitions fired at least once.", Dict("spawns" => sp)
    )
    return mk(
        "transitions_fire", "fail",
        "Dead transition(s): $(join(dead, ", ")). Each never spawned over the horizon — check the " *
            "rate, whether the input places are ever simultaneously available, and whether the " *
            "transition is deactivated. A dead transition silently removes a mechanism from the model.",
        Dict("dead" => dead, "spawns" => sp)
    )
end

# 2. Spawning without completing means nothing ever reached an outcome: cycle times longer than the
#    horizon, or a success probability of zero.
function check_completions_exist(prob, doc, opts)
    su = successes(prob)
    total = sum(values(su); init = 0.0)
    total > 0 && return mk(
        "completions_exist", "pass",
        "$(Int(round(total))) successful completions over the horizon.", Dict("successes" => su)
    )
    return mk(
        "completions_exist", "fail",
        "No transition completed successfully. Nothing reached an outcome, so any headline number " *
            "is zero by construction — check cycle times against the horizon and that success " *
            "probabilities are above zero.",
        Dict("successes" => su, "ongoing_at_end" => length(prob.ongoing_firings))
    )
end

# 3. THE GATE. If no resource is ever contended, the model is a deterministic flow calculation with
#    extra steps and a spreadsheet answers it more cheaply. Contention is what makes a marginal
#    quantity — a shadow price, an attributable delta — mean anything at all.
function check_contention_exists(prob, doc, opts)
    sat, starved = saturation(prob), starved_places(prob)
    contended = sat < 1.0 || !isempty(starved)
    ev = Dict("min_saturation" => sat, "starved_places" => starved)
    contended && return mk(
        "contention_exists", "pass",
        "Resources are contended (minimum in-flight fill fraction $(round(sat; sigdigits = 3))" *
            (isempty(starved) ? "" : "; pools run to zero: $(join(starved, ", "))") * ").",
        ev
    )
    return mk(
        "contention_exists", "fail",
        "Nothing binds: every in-flight instance received its full claim every tick and no pool ran " *
            "to zero. Under no contention the marginal value of a resource is just its accounting " *
            "cost, which a spreadsheet computes — so this model cannot support a shadow-price or " *
            "attributable-value conclusion. Either scarcify the intended resource (lower its initial " *
            "marking, raise demand or the cycle time) or stop and say a simulation is not needed.",
        ev
    )
end

# 4. Contention in the WRONG place answers a different question than the one asked. If the brief says
#    people are scarce and cash turns out to bind, the headline number is a cash number.
function check_intended_resource_binds(prob, doc, opts)
    want = get(opts, "expect_binds", "")
    isempty(want) && return mk(
        "intended_resource_binds", "warn",
        "No --expect-binds given, so which resource binds was not checked against intent. Pass the " *
            "place the brief names as scarce to verify the model answers the question asked.",
        Dict("starved_places" => starved_places(prob))
    )
    wanted = String.(split(want, ','))
    starved = starved_places(prob)
    hit = intersect(wanted, starved)
    ev = Dict("expected" => wanted, "starved_places" => starved)
    !isempty(hit) && return mk(
        "intended_resource_binds", "pass",
        "The intended resource binds: $(join(hit, ", ")) runs to zero.", ev
    )
    return mk(
        "intended_resource_binds", "fail",
        "Expected $(join(wanted, ", ")) to bind, but the pools that run to zero are " *
            "$(isempty(starved) ? "none" : join(starved, ", ")). A marginal number read off this run " *
            "prices the wrong constraint.",
        ev
    )
end

# 5. The highest-yield static check. The DEFAULT modality is `consumed`, so a resource the author
#    thinks of as reusable headcount is silently spent unless it was declared `@conserved`. This
#    compares declared intent against the document, before any dynamics are involved.
function check_modality_intent(prob, doc, opts)
    want = get(opts, "expect_reusable", "")
    isempty(want) && return mk(
        "modality_intent", "warn",
        "No --expect-reusable given. Any resource that is held and returned (headcount, equipment, " *
            "bench capacity) must be declared conserved; the default is consumed.",
        Dict()
    )
    wanted = String.(split(want, ','))
    arcs = get(doc, "arcs", get(doc, "reactants", []))
    bad, ok = Dict{String, Any}(), String[]
    for place in wanted
        lhs = [
            a for a in arcs if get(a, "side", "") == "lhs" &&
                String(get(a, "place", get(a, "species", ""))) == place
        ]
        if isempty(lhs)
            bad[place] = "no input arc references this place"
            continue
        end
        offenders = [
            a for a in lhs if String(get(get(a, "modality", Dict()), "return", "consumed")) != "conserved"
        ]
        isempty(offenders) ? push!(ok, place) :
            (bad[place] = "$(length(offenders)) of $(length(lhs)) input arcs are not conserved")
    end
    isempty(bad) && return mk(
        "modality_intent", "pass",
        "Every resource declared reusable is conserved on all its input arcs: $(join(ok, ", ")).",
        Dict("conserved" => ok)
    )
    return mk(
        "modality_intent", "fail",
        "Reusable-resource mismatch: " * join(["$k ($v)" for (k, v) in bad], "; ") *
            ". A conserved arc occupies its place for the cycle time and returns the tokens intact; " *
            "the default consumed arc spends them. Declaring headcount consumed drains the pool and " *
            "understates capacity for the whole run.",
        Dict("mismatched" => bad, "conserved" => ok)
    )
end

# 6. In-flight instances at the horizon, or trajectories still moving in the final window, mean the
#    reported number is a snapshot of a transient rather than the system's behaviour.
function check_horizon_sufficient(prob, doc, opts)
    ongoing = length(prob.ongoing_firings)
    moving = Dict{String, Any}()
    for s in Main.place_names(prob)
        s in names(prob.sol) || continue
        late, early = tail_windows(prob.sol[!, s])
        late === nothing && continue
        scale = max(abs(early), 1.0)
        rel = (late - early) / scale
        abs(rel) > 0.2 && (moving[s] = round(rel; sigdigits = 3))
    end
    ev = Dict("ongoing_at_end" => ongoing, "still_moving" => moving)
    (ongoing == 0 && isempty(moving)) && return mk(
        "horizon_sufficient", "pass",
        "The horizon covers the process: nothing in flight at the end and no pool still trending.", ev
    )
    parts = String[]
    ongoing > 0 && push!(parts, "$ongoing instance(s) still in flight at the horizon")
    isempty(moving) || push!(parts, "pools still trending in the final quarter ($(join(keys(moving), ", ")))")
    return mk(
        "horizon_sufficient", "warn",
        join(parts, "; ") * ". The reported number is horizon-truncated — extend tspan until the " *
            "quantity of interest settles, or state the truncation as a limitation.",
        ev
    )
end

# 7. A pool that only ever grows usually means a missing sink or a mis-scaled intake: work arrives
#    faster than any capacity can clear it, and the queue — not the decision — dominates the result.
function check_bounded_pools(prob, doc, opts)
    # `--tally=a,b`: places the author declares to be cumulative counters (a launch tally, a units-
    # shipped total). Those SHOULD only grow — a product place with no sink is how the ledger realizes
    # reward at completion — so flagging them would train the reader to ignore this check.
    tallies = Set(String.(split(get(opts, "tally", ""), ','; keepempty = false)))
    growing = Dict{String, Any}()
    for (s, st) in pools(prob)
        s in tallies && continue
        base = max(abs(st["initial"]), 1.0)
        if st["final"] >= st["max"] && st["final"] > 3 * base
            late, early = tail_windows(prob.sol[!, s])
            late !== nothing && late > early && (growing[s] = Dict("final" => st["final"], "initial" => st["initial"]))
        end
    end
    isempty(growing) && return mk(
        "bounded_pools", "pass",
        "No pool grows without bound" *
            (isempty(tallies) ? "." : " (declared tallies exempt: $(join(sort(collect(tallies)), ", ")))."),
        Dict("tallies" => sort(collect(tallies)))
    )
    return mk(
        "bounded_pools", "warn",
        "Unbounded growth in $(join(keys(growing), ", ")): each ends at its maximum, well above its " *
            "initial level, and is still rising. Check for a missing downstream transition or an " *
            "intake rate that outruns any plausible capacity — an ever-growing queue can swamp the " *
            "effect being measured.",
        Dict("growing" => growing)
    )
end

# 8. A model with no ledger flows cannot produce a value conclusion. Also catches the specific error
#    of attaching reward to a stock rather than to a product place, which credits value at intake.
function check_ledger_wired(prob, doc, opts)
    lt = Main.ledger_totals(prob)
    ev = Dict("reward" => lt.reward, "cost" => lt.cost)
    (lt.reward != 0.0 || lt.cost != 0.0) && return mk(
        "ledger_wired", "pass",
        "Ledger flows are non-trivial (reward $(round(lt.reward; sigdigits = 4)), cost $(round(lt.cost; sigdigits = 4))).",
        ev
    )
    return mk(
        "ledger_wired", "warn",
        "No cost or reward accrued. Any value-denominated conclusion will be zero — set the per-unit " *
            "cost/reward/valuation on the places that carry money, and put reward on the PRODUCT " *
            "place so it realizes at completion rather than at intake.",
        ev
    )
end

# 9. Success probabilities pinned at 1 make the model deterministic in outcome: the ensemble spread
#    collapses and the risk-adjustment the framework exists to compute is absent.
function check_stochasticity(prob, doc, opts)
    ps = Any[]
    for t in get(doc, "transitions", [])
        haskey(t, "prob_of_success") && push!(ps, t["prob_of_success"])
    end
    risky = count(p -> !(p isa Real) || p < 1.0, ps)
    ev = Dict("transitions_with_pos" => length(ps), "below_one" => risky)
    (isempty(ps) || risky > 0) && return mk(
        "stochasticity", "pass",
        isempty(ps) ? "No explicit success probabilities (all transitions complete)." :
            "$risky of $(length(ps)) transitions carry a success probability below 1.",
        ev
    )
    return mk(
        "stochasticity", "warn",
        "Every declared success probability is 1.0 — nothing can fail, so the model is " *
            "outcome-deterministic and an ensemble will show no spread. If the process has attrition, " *
            "that is the number the decision is most sensitive to.",
        ev
    )
end

const CHECKS = [
    check_transitions_fire, check_completions_exist, check_contention_exists,
    check_intended_resource_binds, check_modality_intent, check_horizon_sufficient,
    check_bounded_pools, check_ledger_wired, check_stochasticity,
]

"""
    run_checks(prob, doc, opts) -> Vector{Dict}

Run every plausibility check against a finished run `prob` and its source document `doc`. A check
that throws is reported as a `warn` naming the failure rather than aborting the gate — a diagnostic
tool that dies on an edge case is worse than one that admits it could not decide.
"""
function run_checks(prob, doc, opts)
    out = Dict{String, Any}[]
    for f in CHECKS
        try
            push!(out, f(prob, doc, opts))
        catch e
            push!(
                out, mk(
                    string(nameof(f)), "warn",
                    "Check could not run: $(sprint(showerror, e))", Dict()
                )
            )
        end
    end
    return out
end

end # module
