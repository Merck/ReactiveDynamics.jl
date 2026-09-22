#!/usr/bin/env julia
#
# acceptance.jl — the agent layer's known-answer test (charter W1 §11.3).
#
#   julia --project=. agent/scripts/acceptance.jl
#
# The workflow's only real test is whether it reproduces an answer someone already computed by hand.
# `agent/studies/marginal_scientist/` is that study: the published case study
# (`docs/literate/case_studies/marginal_scientist.jl`) builds the same model in the Julia DSL and
# reports a shadow price of ≈ +$19M ± $2.5M. This script drives the same question through the
# AGENT path — two JSON documents, the `rdx` CLI, `host.jl` as the only Julia — and checks the answer
# survives the trip. If it does not, the agent-facing surface is no longer the same engine as the
# hand-authored one, which is the failure this whole layer exists to preclude.
#
# It asserts at two strengths, deliberately:
#
#   HARD  — the scientific claim: the delta lands inside the published band and clears its own noise.
#           Robust to anything that reshuffles the RNG stream without changing the semantics.
#   SOFT  — bit-reproducibility against the recorded value. A drift here is reported, not fatal:
#           Julia's `hash` is not stable across versions, so the seed stream can move for reasons
#           that are not regressions. Worth knowing about; not worth blocking a merge over.
#
# Exit 0 = the study reproduces. Exit 1 = it does not.

using JSON

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const STUDY = joinpath("agent", "studies", "marginal_scientist")
const RDX = joinpath("agent", "scripts", "rdx.jl")
const HOST = joinpath(STUDY, "host.jl")

# The published answer, and the exact value this workflow produced when the study landed (f8379b0).
const PUBLISHED_DELTA = 19.0
const PUBLISHED_SE = 2.5
const RECORDED_DELTA = 18.99102052953539
const RECORDED_SE = 2.5279131809320163

failures = String[]
notes = String[]

fail(msg) = (push!(failures, msg); println("  FAIL  $msg"))
pass(msg) = println("  ok    $msg")
note(msg) = (push!(notes, msg); println("  note  $msg"))

# Run an `rdx` subcommand, returning (exit_code, parsed_json_or_nothing). stderr is passed through so
# the host-include notice and any warning stay visible in the log.
function rdx(args...)
    cmd = Cmd(`julia --project=$ROOT $RDX $args`; dir = ROOT)
    out = IOBuffer()
    code = try
        run(pipeline(cmd; stdout = out, stderr = stderr))
        0
    catch e
        e isa ProcessFailedException ? e.procs[1].exitcode : rethrow()
    end
    text = String(take!(out))
    parsed = try
        JSON.parse(text)
    catch
        isempty(text) || println(text)
        nothing
    end
    return code, parsed
end

println("acceptance: the marginal-scientist known-answer study")
println()

# ── 1. both arms are well-formed documents ────────────────────────────────────────────────────────
println("1. validate")
for arm in ("model.rdj.json", "arm_hire_fifth.rdj.json")
    code, _ = rdx("validate", joinpath(STUDY, arm), "--host=$HOST")
    code == 0 ? pass("$arm: zero diagnostics") : fail("$arm: validate exited $code")
end

# ── 2. the arms differ ONLY in the lever ──────────────────────────────────────────────────────────
# The two-arm convention (charter W1) is that the diff IS the intervention. A drifted second arm
# would still validate and still produce a delta — just not the delta of the stated lever.
println("\n2. the diff is the lever")
let
    a = JSON.parsefile(joinpath(ROOT, STUDY, "model.rdj.json"))
    b = JSON.parsefile(joinpath(ROOT, STUDY, "arm_hire_fifth.rdj.json"))
    init(d, name) = only(p for p in d["places"] if p["name"] == name)["init"]
    for k in ("meta", "transitions", "arcs", "params", "population")
        get(a, k, nothing) == get(b, k, nothing) || fail("arms differ in `$k` — the lever is not the only change")
    end
    ia, ib = init(a, "scientist"), init(b, "scientist")
    (ia == 4.0 && ib == 5.0) ? pass("scientist.init 4.0 → 5.0, nothing else") :
        fail("scientist.init is $ia → $ib, expected 4.0 → 5.0")
end

# ── 3. the plausibility gate still passes on the brief's PRE-REGISTERED intent ─────────────────────
println("\n3. diagnose")
let
    code, d = rdx(
        "diagnose", joinpath(STUDY, "model.rdj.json"), "--host=$HOST", "--seed=1",
        "--expect-binds=scientist", "--expect-reusable=scientist", "--tally=launch"
    )
    if d === nothing
        fail("diagnose produced no JSON")
    else
        r = get(d, "result", Dict())
        nfail = get(r, "failed", -1)
        code == 0 && nfail == 0 ?
            pass("0 fail ($(get(r, "passed", "?")) pass, $(get(r, "warned", "?")) warn)") :
            fail("diagnose exited $code with $nfail failing check(s)")
    end
end

# ── 4. the number ─────────────────────────────────────────────────────────────────────────────────
println("\n4. the shadow price")
let
    code, d = rdx(
        "effect", joinpath(STUDY, "model.rdj.json"), joinpath(STUDY, "arm_hire_fifth.rdj.json"),
        "--host=$HOST", "--metric=enpv", "--scale=0.001", "--annual=0.08",
        "--seeds=48", "--root-seed=2026"
    )
    if d === nothing
        fail("effect produced no JSON")
    else
        e = get(get(get(d, "result", Dict()), "metrics", Dict()), "enpv", Dict())
        delta = get(e, "delta", NaN)
        se = get(e, "se", NaN)
        sig = get(e, "significant_2se", false)
        println("        delta $(round(delta, digits = 3))  se $(round(se, digits = 3))  exit $code")

        # HARD: the scientific claim.
        abs(delta - PUBLISHED_DELTA) <= PUBLISHED_SE ?
            pass("delta \$$(round(delta, digits = 2))M is inside the published \$$(PUBLISHED_DELTA)M ± \$$(PUBLISHED_SE)M band") :
            fail("delta \$$(round(delta, digits = 2))M is OUTSIDE the published \$$(PUBLISHED_DELTA)M ± \$$(PUBLISHED_SE)M band")
        sig ? pass("effect clears 2·se — a real marginal gain, not ensemble noise") :
            fail("effect does not clear 2·se at 48 members per arm")

        # SOFT: bit-reproducibility.
        if isapprox(delta, RECORDED_DELTA; atol = 1.0e-9) && isapprox(se, RECORDED_SE; atol = 1.0e-9)
            pass("bit-identical to the recorded run")
        else
            note(
                "the seed stream moved: delta $delta / se $se vs recorded $RECORDED_DELTA / " *
                    "$RECORDED_SE. Not necessarily a regression (Julia's `hash` is not stable across " *
                    "versions), but if the engine did not change, find out why."
            )
        end
    end
end

println()
if isempty(failures)
    println("acceptance: PASS — the workflow reproduces the published answer through the JSON path.")
    isempty(notes) || println("  ($(length(notes)) note(s) above.)")
    exit(0)
else
    println("acceptance: FAIL — $(length(failures)) check(s) did not hold:")
    for f in failures
        println("  · $f")
    end
    exit(1)
end
