# twin.jl — the NocoBase → DyVE bridge for the Pharma R&D digital twin.
#
# Three things live here, and nothing else:
#
#   1. `compile_model(snapshot, arm)` — the MIGRATION: it reads a NocoBase-shaped snapshot (the
#      collections of the prototype deck at tick 42) and emits an eval-free ReactiveDynamics model
#      document (`.rdj.json`). Rows become places, transitions, arcs, params and rules; the NocoBase
#      workflows become rules. The document is plain JSON, so it can be validated, diffed, stored back
#      in NocoBase, or handed to the agent layer's `rdx` CLI.
#   2. The host seam — the `ProgramToken` kind, the `rebalance` portfolio controller (the stand-in for
#      the deck's MPC / NocoBase WF5 "external solver"), and `ClinOpsAgent`, an AlgebraicAgents agent
#      that simulates patient accrual (Poisson, site_count × rate) and feeds the net through a wire.
#   3. Metrics — risk-adjusted eNPV with patent-life erosion, spend, and the per-tick series the
#      dashboards and the NocoBase write-back read.
#
# No ReactiveDynamics or AlgebraicAgents source is modified; everything is public API plus the
# documented host seam (a token kind and a registry function).

using ReactiveDynamics
using ReactiveDynamics: ReactionNetworkProblem, from_json_model, validate
using AlgebraicAgents
using Random, Distributions, Statistics, Dates, DataFrames
import JSON

const RD = ReactiveDynamics

# ── 1. Structured token: a NocoBase `therapeutics_programs` row becomes a first-class agent ──────────

@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct ProgramToken
        code::Symbol
        target::String
        modality::String
        area::String
        stage::String
        next_stage::String
        pos::Float64        # cumulative probability of technical success (pos_estimate)
        stage_pos::Float64  # probability of passing the readout that ends the current stage
        enpv::Float64       # unadjusted eNPV, USD M (enpv_usd_m)
        readout_t::Float64  # month of a passed readout; NaN while pending (written by @advance)
    end
end

ReactiveDynamics.log_token_fields(t::ReactiveDynamics.ProgramToken) = (; readout_t = t.readout_t)

function program_token(p)
    return ReactiveDynamics.ProgramToken(
        p.code, :Program, nothing, Tuple{Symbol, Float64, ReactiveDynamics.Firing}[],
        Symbol(p.code), p.target, p.modality, p.area, p.stage, p.next_stage,
        p.pos, p.stage_pos, p.enpv, NaN,
    )
end

# ── 2. Read the NocoBase snapshot and join it into one row per program ───────────────────────────────

const SNAPSHOT_PATH = joinpath(@__DIR__, "data", "nocobase_snapshot.json")
const SHARED = ("ms", "gmp", "viv", "gpu", "clin")   # shared instrument / department pools
const GMP_SLOTS_PER_BATCH_MONTH = 3.0
const GMP_BATCH_CRO = 0.3                              # $M per month while a GMP batch runs

load_snapshot(path = SNAPSHOT_PATH) = JSON.parsefile(path)

pid(code) = "p" * split(code, "-")[2]                  # "PRG-104" → "p104" (a legal place/param name)
_num(x) = x === nothing ? NaN : Float64(x)

"""
    program_table(snap) -> Vector{NamedTuple}

Join `therapeutics_programs` with `compounds_candidates`, `clinical_trials`, `studies` and
`resource_allocations` on `program_code`, one row per program. This is the relational read
NocoBase's `?appends=` would do.
"""
function program_table(snap)
    byc(coll) = Dict(r["program_code"] => r for r in snap[coll])
    cmp, trl, std, alc = byc("compounds_candidates"), byc("clinical_trials"), byc("studies"),
        byc("resource_allocations")
    return map(snap["therapeutics_programs"]) do r
        c = r["program_code"]
        a, s = alc[c], std[c]
        cp, tr = get(cmp, c, nothing), get(trl, c, nothing)
        months = Float64(s["months_remaining_full_speed"])
        clinical = tr !== nothing
        (;
            code = c, id = pid(c), target = r["target"], modality = r["modality"],
            area = r["therapeutic_area"], stage = r["current_stage"], next_stage = r["next_stage"],
            pos = Float64(r["pos_estimate"]), stage_pos = Float64(r["stage_pos"]),
            enpv = Float64(r["enpv_usd_m"]), risk_adj = r["pos_estimate"] * r["enpv_usd_m"],
            months, study_code = s["study_code"],
            wet_alloc = Float64(a["wet_fte"]), wet_full = Float64(a["wet_fte_full_speed"]),
            shared = [(k, Float64(a[k])) for k in SHARED if a[k] > 0],
            cro = Float64(a["cro_usd_m_month"]),
            herg = cp === nothing ? NaN : _num(cp["herg_ic50_um"]),
            clinical,
            trial_id = clinical ? tr["trial_id"] : "",
            need = clinical ? (tr["target_enrollment"] - tr["current_enrollment"]) / months : 0.0,
            sites = clinical ? Float64(tr["site_count"]) : 0.0,
            accrual_rate = clinical ? Float64(tr["accrual_per_site_month"]) : 0.0,
            kits = cp === nothing ? NaN : _num(cp["kits_on_hand"]),
            reorder = cp === nothing ? NaN : _num(cp["reorder_point_kits"]),
            batch = cp === nothing ? NaN : _num(cp["batch_kits"]),
        )
    end
end

resource_table(snap) = Dict(r["resource_id"] => r for r in snap["resources_capacity"])

# ── 3. The scenario arms (the decision the deck's headline claims to answer) ─────────────────────────

const ARMS = (
    static = (kill_at = 12.0, replan = :annual, label = "Static annual plan"),
    kill_only = (kill_at = 0.0, replan = :annual, label = "Tox kill switch only"),
    replan_only = (kill_at = 12.0, replan = :monthly, label = "Monthly re-plan only"),
    closed_loop = (kill_at = 0.0, replan = :monthly, label = "Closed loop"),
)
const ANNUAL_REVIEW = 12.0   # month of the annual portfolio review in the static arms

# ── 4. ExprNode helpers (the closed, eval-free expression IR of the model document) ──────────────────

C(v) = Dict("node" => "const", "value" => v)
SYM(s) = Dict("node" => "const", "value" => string(s), "symbol" => true)
PAR(n) = Dict("node" => "ref", "kind" => "param", "name" => string(n))
PLC(n) = Dict("node" => "ref", "kind" => "place", "name" => string(n))
CALL(op, args...) = Dict("node" => "call", "op" => op, "args" => collect(args))
SAMPLE(d, args...) = Dict("node" => "sample", "dist" => d, "args" => collect(args))
PORT(p) = Dict("node" => "externalref", "port" => p)
const NOW = Dict("node" => "timeref")
const PERSTEP = Dict("allocation" => "perstep", "return" => "consumed", "blocking" => "block")

setmark(place, value; mode = "set") =
    Dict("verb" => "set_marking", "name" => place, "value" => value, "mode" => mode)
rule(id, guard, action; once = false) =
    Dict("id" => id, "guard" => guard, "action" => action, "fire_mode" => once ? "once" : "every_tick")

# ── 5. The migration: NocoBase snapshot → ReactiveDynamics model document ────────────────────────────

"""
    compile_model(snap, arm; horizon = 24.0) -> Dict

Compile the NocoBase snapshot into a ReactiveDynamics model document for scenario `arm`.

Mapping (see README.md for the full table):

  - `therapeutics_programs` row → a `ProgramToken` in the structured place `Program`;
  - `studies` row (the in-flight study of the current stage) → transition `study_<p>`, spawned once
    (enabler place `start_<p>`), bound to its program by `@select(Program, code == …)`; cycle time =
    remaining months × Gamma noise; terminal success = `stage_pos` (a failed readout retires the token,
    which is the NocoBase `Terminated` status);
  - `resources_capacity` row → a shared place refreshed to its monthly capacity by a rule, consumed
    per month (`perstep`) by every study that needs it; cost per unit from the row;
  - `resource_allocations` junction (wet-lab FTE) → a per-program quota place `wet_<p>` refreshed to the
    param `q_<p>`; the controller re-plans by writing `q_<p>`;
  - `clinical_trials` row → a patient place `pts_<p>` fed by the `ClinOpsAgent` through port
    `accrual_<p>`; the study consumes the planned accrual per month, so a lag slows it (WF4);
  - spawn gating: every transition's rate reads its own enabler place (`start_<p>`, `slot_<p>`), so the
    rate is 0 while an instance is in flight. A rate above `capacity` is otherwise queued in the
    engine's spawn backlog and fires later, after the condition that raised it has gone (README, F5);
  - `compounds_candidates.kits_on_hand` → place `kits_<p>`; when stock drops below the reorder point,
    a GMP batch transition occupies Suite B for the 90-day lead time and restocks (WF3);
  - WF2 tox branch → the `kill_<p>` rule: hERG IC50 below the threshold halts the study (`act_<p> = 0`)
    from month `kill_at` (0 = the automatic kill switch, 12 = acted on at the annual review);
  - WF5 portfolio rebalancing → the `replan` rule, which invokes the host `rebalance` function;
  - `Studies.status` → param `alive_<p>` (study in flight and not halted), kept in sync each month by
    the host `sync_status` function, so a GMP batch is only reordered for a live study.
"""
function compile_model(snap, arm::Symbol; horizon = 24.0)
    cfg = ARMS[arm]
    progs = program_table(snap)
    res = resource_table(snap)
    params, transitions, arcs, rules, inputs = Any[], Any[], Any[], Any[], Any[]
    places = Any[Dict("name" => "Program", "structured" => true)]
    addp(n, v) = push!(params, Dict("name" => n, "value" => Float64(v)))
    addplace(n, init; cost = nothing) =
        push!(places, cost === nothing ? Dict("name" => n, "init" => init) : Dict("name" => n, "init" => init, "cost" => cost))

    addp("kill_at", cfg.kill_at)
    addp("herg_threshold", snap["kill_switch"]["herg_ic50_um_threshold"])
    addp("cap_wet", res["wet_fte"]["total_capacity_units"])
    addp("cro_envelope", snap["budget"]["quarterly_budget_usd_m"] / 3)

    refresh = Any[]
    for k in SHARED
        addp("cap_$k", res[k]["total_capacity_units"])
        addplace(k, res[k]["total_capacity_units"]; cost = res[k]["cost_usd_m_per_unit_month"])
        push!(refresh, rule("refresh_$k", C(true), setmark(k, PAR("cap_$k"))))
    end
    addplace("cro", snap["budget"]["quarterly_budget_usd_m"] / 3; cost = 1.0)
    push!(refresh, rule("refresh_cro", C(true), setmark("cro", PAR("cro_envelope"))))
    addplace("readouts_passed", 0)

    kills, accrual = Any[], Any[]
    wet_cost = res["wet_fte"]["cost_usd_m_per_unit_month"]
    for p in progs
        id, t = p.id, "study_$(p.id)"
        addp("act_$id", 1)
        addp("q_$id", p.wet_alloc)
        addp("prio_$id", 1)
        addp("ps_$id", p.stage_pos)
        addp("alive_$id", 1)
        addplace("start_$id", 1)
        addplace("go_$id", 1)
        addplace("wet_$id", p.wet_alloc; cost = wet_cost)
        addplace("killed_$id", 0)

        push!(
            transitions, Dict(
                "id" => t, "name" => t, "rate_mode" => "deterministic", "rate" => PLC("start_$id"), "capacity" => 1,
                "cycletime" => CALL("*", C(p.months), SAMPLE("Gamma", C(25.0), C(0.04))),
                "prob_of_success" => PAR("ps_$id"), "priority" => PAR("prio_$id"),
            )
        )
        lhs(tr, place, m; perstep = true) = push!(
            arcs, perstep ?
                Dict("transition" => tr, "side" => "lhs", "place" => place, "multiplicity" => m, "modality" => PERSTEP) :
                Dict("transition" => tr, "side" => "lhs", "place" => place, "multiplicity" => m)
        )
        lhs(t, "start_$id", 1; perstep = false)
        push!(
            arcs, Dict(
                "transition" => t, "side" => "lhs",
                "predicate" => Dict("kind" => "Program", "clauses" => [["code", "==", SYM(p.code)]]),
            )
        )
        lhs(t, "go_$id", 1)
        lhs(t, "wet_$id", p.wet_full)
        for (k, v) in p.shared
            lhs(t, k, v)
        end
        lhs(t, "cro", p.cro)
        push!(arcs, Dict("transition" => t, "side" => "rhs", "advance" => Dict("field" => "readout_t", "value" => NOW)))
        push!(arcs, Dict("transition" => t, "side" => "rhs", "place" => "readouts_passed", "multiplicity" => 1))

        push!(refresh, rule("refresh_wet_$id", C(true), setmark("wet_$id", PAR("q_$id"))))
        push!(refresh, rule("refresh_go_$id", C(true), setmark("go_$id", PAR("act_$id"))))

        if p.clinical
            addplace("pts_$id", round(p.need; digits = 2))
            lhs(t, "pts_$id", round(p.need; digits = 2))
            push!(inputs, Dict("port" => "accrual_$id", "default" => C(0.0)))
            push!(accrual, rule("accrue_$id", C(true), setmark("pts_$id", PORT("accrual_$id"); mode = "inc")))
        end
        if p.clinical && !isnan(p.kits)
            addplace("kits_$id", p.kits)
            lhs(t, "kits_$id", round(p.need; digits = 2))
            b = "gmp_batch_$id"
            push!(
                transitions, Dict(
                    "id" => b, "name" => b, "rate_mode" => "deterministic", "capacity" => 1,
                    "rate" => CALL(
                        "&&", CALL("&&", CALL("<", PLC("kits_$id"), C(p.reorder)), CALL(">", PAR("alive_$id"), C(0))),
                        CALL(">", PLC("slot_$id"), C(0))
                    ),
                    "cycletime" => C(Float64(snap["gmp_batch"]["lead_time_months"])),
                    "prob_of_success" => C(1.0), "priority" => PAR("prio_$id"),
                )
            )
            addplace("slot_$id", 1)
            lhs(b, "slot_$id", 1; perstep = false)
            push!(arcs, Dict("transition" => b, "side" => "rhs", "place" => "slot_$id", "multiplicity" => 1))
            lhs(b, "gmp", Float64(get(snap["gmp_batch"], "slots_per_month", GMP_SLOTS_PER_BATCH_MONTH)))
            lhs(b, "cro", GMP_BATCH_CRO)
            push!(arcs, Dict("transition" => b, "side" => "rhs", "place" => "kits_$id", "multiplicity" => p.batch))
        end
        if !isnan(p.herg)
            addp("herg_$id", p.herg)
            guard = CALL("&&", CALL(">=", NOW, PAR("kill_at")), CALL("<", PAR("herg_$id"), PAR("herg_threshold")))
            act = Dict(
                "verb" => "seq", "stmts" => [
                    Dict("verb" => "set_params", "assigns" => [Dict("name" => "act_$id", "value" => C(0.0))]),
                    setmark("killed_$id", C(1.0)),
                ]
            )
            push!(kills, rule("kill_$id", guard, act; once = true))
        end
    end

    replan_guard = cfg.replan === :monthly ? C(true) : CALL(">=", NOW, C(ANNUAL_REVIEW))
    replan = rule("replan", replan_guard, Dict("verb" => "invoke", "fn" => "rebalance", "args" => []); once = cfg.replan !== :monthly)

    # Rule order matters: decisions (kill, re-plan) first, then the capacity refreshes that apply them
    # to next month's pools, then the accrual that arrives over the AlgebraicAgents wire.
    push!(rules, rule("sync_status", C(true), Dict("verb" => "invoke", "fn" => "sync_status", "args" => [])))
    append!(rules, kills)
    push!(rules, replan)
    append!(rules, refresh)
    append!(rules, accrual)

    return Dict(
        "rd_format" => "reactive-dynamics-model", "version" => "1.0",
        "meta" => Dict("tspan" => horizon, "dt" => 1.0),
        "params" => params, "places" => places, "transitions" => transitions, "arcs" => arcs,
        "rules" => rules, "inputs" => inputs,
    )
end

# ── 6. The portfolio controller (host seam; NocoBase WF5 / the deck's MPC step) ──────────────────────

"""
    make_rebalance(progs, replans) -> Function

The `rebalance` registry function the `replan` rule invokes. It re-plans the wet-lab FTE pool across
the studies still in flight: programs are ranked by risk-adjusted eNPV per wet-lab FTE of full-speed
demand, and each gets `min(demand, what is left)`. Under the erosion objective used here (value lost
is proportional to slip, slip to the unfunded share of demand) this greedy ranking is the exact
optimum of the fractional knapsack, so it is an upper bound on what a smarter MPC could do in this
model. It also sets every study's allocator priority to its risk-adjusted eNPV, so contended shared
pools (GMP Suite B, mass spec) go to value first. Deterministic and RNG-free (ADR 0011 O1–O4); each
call appends the new plan to `replans` for the write-back.
"""
function make_rebalance(progs, replans)
    order = sort(collect(progs); by = p -> (-p.risk_adj / p.wet_full, p.code))
    return function (state, firing)
        inflight = Set(string(state.network[f.i, :transName]) for f in state.ongoing_firings)
        remaining = Float64(state.p[:cap_wet])
        plan = Dict{String, Float64}()
        for p in order
            active = ("study_$(p.id)" in inflight) && state.p[Symbol("act_$(p.id)")] > 0
            q = active ? min(p.wet_full, remaining) : 0.0
            remaining -= q
            state.p[Symbol("q_$(p.id)")] = q
            state.p[Symbol("prio_$(p.id)")] = p.risk_adj / 100
            plan[p.code] = q
        end
        (isempty(replans) || last(replans)[2] != plan) && push!(replans, (state.t, plan))
        return nothing
    end
end

"The `sync_status` registry function: `alive_<p>` = 1 while program p's study is in flight and not halted."
function make_sync_status(progs)
    return function (state, firing)
        inflight = Set(string(state.network[f.i, :transName]) for f in state.ongoing_firings)
        for p in progs
            live = ("study_$(p.id)" in inflight) && state.p[Symbol("act_$(p.id)")] > 0
            state.p[Symbol("alive_$(p.id)")] = live ? 1.0 : 0.0
        end
        return nothing
    end
end

# ── 7. Clinical operations as an AlgebraicAgents agent (the Corvane "contract bus" idea, minimal) ────

"""
    ClinOpsAgent

Simulates monthly patient accrual for every recruiting trial: `Poisson(site_count × rate)`, the deck's
`np.random.poisson(prog.trial.site_count * 1.8)`. It is a separate AlgebraicAgents agent with its own
seeded RNG, wired into the portfolio net through one port per trial (`accrual_<p>`). Before
publishing, each value is checked against a producer-side contract (range [0, 40] patients per trial
per month, the Corvane `enroll_to_supply` shape); a violation is rejected (the last good value is
kept) and recorded as an event.
"""
@aagent struct ClinOpsAgent
    trials::Vector{Any}
    rng::Random.Xoshiro
    accrual::Dict{String, Float64}
    enrolled::Dict{String, Float64}
    history::Vector{Any}
    violations::Vector{Any}
    t::Float64
    dt::Float64
    horizon::Float64
end

const ACCRUAL_CONTRACT = (lo = 0.0, hi = 40.0)

function ClinOpsAgent(name::AbstractString, progs, seed, horizon)
    trials = Any[(; port = "accrual_$(p.id)", code = p.code, trial_id = p.trial_id, sites = p.sites, rate = p.accrual_rate, need = p.need) for p in progs if p.clinical]
    a = ClinOpsAgent(name, trials, Random.Xoshiro(hash((seed, :clinops))), Dict{String, Float64}(), Dict{String, Float64}(), Any[], Any[], 0.0, 1.0, Float64(horizon))
    _draw_accrual!(a)
    return a
end

function _draw_accrual!(a::ClinOpsAgent)
    for tr in a.trials
        n = Float64(rand(a.rng, Poisson(tr.sites * tr.rate)))
        if !(ACCRUAL_CONTRACT.lo <= n <= ACCRUAL_CONTRACT.hi)
            push!(a.violations, (a.t, tr.code, n))
            n = get(a.accrual, tr.port, 0.0)
        end
        a.accrual[tr.port] = n
        a.enrolled[tr.port] = get(a.enrolled, tr.port, 0.0) + n
        push!(a.history, (a.t, tr.code, n))
    end
    return a
end

AlgebraicAgents.observables(a::ClinOpsAgent) = Symbol.(tr.port for tr in a.trials)
AlgebraicAgents.getobservable(a::ClinOpsAgent, name::Union{Symbol, AbstractString}) = a.accrual[String(name)]
AlgebraicAgents._step!(a::ClinOpsAgent) = (a.t += a.dt; _draw_accrual!(a); a.t)
AlgebraicAgents._projected_to(a::ClinOpsAgent) = a.t > a.horizon ? true : a.t

# ── 7b. A host-side speed-up, no engine change ────────────────────────────────────────────────────────

"""
    cache_compiles!(rd)

The engine re-compiles a few per-tick expressions (arc multiplicities, `@advance` values, rule
values) through `rd.wrap_fun` every tick. They are the same expressions every tick, so memoizing the
compiler on the problem object is safe: the closures are pure functions of `(state, firing)`. It
makes a run about 10x faster and leaves the trajectory byte-identical (checked in `run_twin.jl`).
"""
function cache_compiles!(rd)
    inner, cache = rd.wrap_fun, Dict{Any, Any}()
    rd.wrap_fun = ex -> ex isa Union{Expr, Symbol} ? get!(() -> inner(ex), cache, ex) : inner(ex)
    return rd
end

# ── 8. Build and run one twin: the AlgebraicAgents hierarchy {portfolio net, clinical ops} ───────────

"""
    build_twin(seed; arm, snap) -> (; root, rd, clin, doc, progs, replans)

Compile the snapshot for `arm`, construct the ReactiveDynamics net from the JSON document (validated,
eval-free), seed the program tokens, and wire the `ClinOpsAgent` into it under one AlgebraicAgents
root. `simulate(root)` runs the whole hierarchy.
"""
function build_twin(seed; arm::Symbol = :closed_loop, snap = load_snapshot(), horizon = 24.0, fast = true)
    progs = program_table(snap)
    doc = compile_model(snap, arm; horizon)
    replans = Any[]
    registry = Dict{Symbol, Any}(:rebalance => make_rebalance(progs, replans), :sync_status => make_sync_status(progs))
    rd = from_json_model(JSON.json(doc); seed = seed, registry = registry, population = [program_token(p) for p in progs])
    fast && cache_compiles!(rd)
    clin = ClinOpsAgent("clinical_ops", progs, seed, horizon)
    root = FreeAgent("rd_twin")
    entangle!(root, rd)
    entangle!(root, clin)
    for tr in clin.trials
        add_wire!(root; from = clin, to = rd, from_var_name = tr.port, to_var_name = tr.port)
    end
    return (; root, rd, clin, doc, progs, replans)
end

run_twin(seed; kwargs...) = (sys = build_twin(seed; kwargs...); simulate(sys.root); sys)

# ── 9. Metrics: slip, patent-life erosion, risk-adjusted eNPV, spend ─────────────────────────────────

const EROSION = 0.015   # eNPV lost per month of slip (patent-life erosion, NocoBase WF4); a sensitivity

trans_names(prob) = Dict(prob.transitions[:transHash][i] => string(prob.network[i, :transName]) for i in eachindex(prob.transitions[:transHash]))

"Per-study (t, fill) pairs from the allocator's `:saturation` log: 1.0 = full speed, 0 = stalled."
function study_saturation(prob)
    names = trans_names(prob)
    out = Dict{String, Vector{Tuple{Float64, Float64}}}()
    for r in prob.log
        r[1] === :saturation || continue
        for (h, f) in r[3:end]
            n = names[h]
            startswith(n, "study_") && push!(get!(out, n[7:end], Tuple{Float64, Float64}[]), (Float64(r[2]), Float64(f)))
        end
    end
    return out
end

"Months of slip per program up to month `upto`: Σ (1 − fill) over the months its study was in flight."
function slips(prob; upto = Inf)
    sat = study_saturation(prob)
    return Dict(id => sum((1 - f for (t, f) in v if t <= upto); init = 0.0) for (id, v) in sat)
end

"Month the kill rule halted program `id`, or `Inf` if it never did."
killed_at(prob, id) = (k = findfirst(>(0), prob.sol[!, Symbol("killed_$id")]); k === nothing ? Inf : prob.sol.t[k])
killed_by(prob, id, upto) = (tk = killed_at(prob, id); isfinite(tk) && tk <= upto)

"""
    riskadj_value(prob, progs; ε = EROSION, upto = Inf) -> Float64

Portfolio risk-adjusted eNPV (USD M) at month `upto`: Σ over programs not killed by then of
`enpv × pos × (1 − ε × slip)`. Readout outcomes enter only through their expectation (the pass/fail
re-mark is a martingale), which removes the binary-outcome noise and leaves the scheduling effect
the decision is about. See `realized_value` for the outcome-dependent version.
"""
function riskadj_value(prob, progs; ε = EROSION, upto = Inf)
    s = slips(prob; upto)
    return sum(progs) do p
        killed_by(prob, p.id, upto) ? 0.0 : p.enpv * p.pos * max(0.0, 1 - ε * get(s, p.id, 0.0))
    end
end

"Outcome-dependent eNPV at the horizon: passed readouts re-marked to pos/stage_pos, failed ones and kills at zero."
function realized_value(prob, progs; ε = EROSION)
    s = slips(prob)
    toks = Dict(String(t.code) => t for t in values(AlgebraicAgents.inners(AlgebraicAgents.getagent(prob, "structured"))))
    return sum(progs) do p
        tok = toks[p.code]
        isfinite(killed_at(prob, p.id)) && return 0.0
        tok.place === :removed && return 0.0
        pos = isnan(tok.readout_t) ? p.pos : p.pos / p.stage_pos
        p.enpv * pos * max(0.0, 1 - ε * get(s, p.id, 0.0))
    end
end

spend_series(prob) = RD.log_scalar_series(prob, :valuation_cost)
spend(prob) = sum(last(spend_series(prob)))
net_value(prob, progs; ε = EROSION) = riskadj_value(prob, progs; ε) - spend(prob)

"Per-tick utilization of each pool: allocated ÷ capacity, from the allocator's `:allocation` log."
function utilization(prob, snap)
    res = resource_table(snap)
    idx(n) = findfirst(==(Symbol(n)), prob.network[:, :placeName])
    ts = Float64[]
    u = Dict(k => Float64[] for k in (SHARED..., "wet_fte"))
    wet = [idx(n) for n in prob.network[:, :placeName] if startswith(string(n), "wet_")]
    for r in prob.log
        r[1] === :allocation || continue
        push!(ts, r[2])
        a = r[3]
        for k in SHARED
            push!(u[k], a[idx(k)] / res[k]["total_capacity_units"])
        end
        push!(u["wet_fte"], sum(a[wet]) / res["wet_fte"]["total_capacity_units"])
    end
    return ts, u
end
