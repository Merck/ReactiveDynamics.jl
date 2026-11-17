#!/usr/bin/env julia
#
# run_twin.jl — run the Pharma R&D digital twin end to end.
#
#   julia --project=demo/nocobase_pharma_twin -e 'using Pkg; Pkg.instantiate()'   # first time only
#   julia --project=demo/nocobase_pharma_twin demo/nocobase_pharma_twin/run_twin.jl [--nseed=40]
#       [--snapshot=path/to/snapshot.json] [--out=path/to/outdir]
#
# --snapshot runs a what-if on an edited copy of the snapshot; --out keeps its results apart from the
# baseline in out/. Both default to this folder's data/nocobase_snapshot.json and out/.
#
# What it does, in order (each step prints a short section):
#
#   §1  read the NocoBase snapshot (data/nocobase_snapshot.json) and check it against the deck
#   §2  compile it into one ReactiveDynamics model document per scenario arm, validate, save
#   §3  run one coupled twin (RD portfolio net + AA clinical-ops agent) and print a month-by-month log
#   §4  run a seeded ensemble per arm and report the effect of each lever with error bars
#   §5  write the results back in NocoBase's own shape (simulation_ticks, simulation_events, plan)
#   §6  draw the model (Petri net, AA wiring) and the dashboard charts into out/
#
# Everything under out/ is regenerated on every run. Same seed, same numbers (on Julia 1.12; seeds come from `hash`).

include(joinpath(@__DIR__, "twin.jl"))
# chartkit.jl lives in the repository's agent layer; a standalone copy of this folder carries it in vendor/.
let repo_copy = joinpath(@__DIR__, "..", "..", "agent", "scripts", "chartkit.jl")
    include(isfile(repo_copy) ? repo_copy : joinpath(@__DIR__, "vendor", "chartkit.jl"))
end
using Printf

argstr(name, default) = (a = findfirst(startswith("--$name="), ARGS); a === nothing ? default : split(ARGS[a], "="; limit = 2)[2])
argval(name, default) = parse(Int, string(argstr(name, default)))
const NSEED = argval("nseed", 40)
const SNAPSHOT = abspath(argstr("snapshot", SNAPSHOT_PATH))
const OUT = abspath(argstr("out", joinpath(@__DIR__, "out")))
mkpath(OUT)
const ROOT_SEED = 2026

banner(s) = (println(); println("="^96); println(s); println("="^96))
usd(x) = @sprintf("\$%.0fM", x)
usd1(x) = @sprintf("\$%.1fM", x)

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§1  NocoBase snapshot → one joined row per program (tick $(load_snapshot(SNAPSHOT)["tick_index"]))")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

snap = load_snapshot(SNAPSHOT)
# eNPV lost per month of slip; a valuation input only (the simulation itself does not read it)
const EPS = Float64(get(get(snap, "valuation", Dict()), "erosion_per_month", EROSION))
SNAPSHOT == abspath(SNAPSHOT_PATH) || println("what-if snapshot: ", SNAPSHOT)
progs = program_table(snap)
res = resource_table(snap)

@printf("%-8s %-9s %-15s %-13s %6s %7s %9s %6s %8s %s\n", "program", "target", "modality", "stage", "PoS", "eNPV", "risk-adj", "months", "wet FTE", "shared pools / month")
for p in progs
    @printf(
        "%-8s %-9s %-15s %-13s %5.1f%% %7.0f %9.1f %6.0f %4.0f/%-3.0f %s\n", p.code, p.target, p.modality, p.stage,
        100p.pos, p.enpv, p.risk_adj, p.months, p.wet_alloc, p.wet_full, join(["$k $(Int(v))" for (k, v) in p.shared], ", ")
    )
end
total_riskadj = sum(p.risk_adj for p in progs)
@printf("portfolio risk-adjusted eNPV at tick 42: %s (sum of the deck's rows: \$1,955.2M)\n", usd1(total_riskadj))
println("\nPool load at tick 42, recomputed from the allocations junction (deck values in brackets):")
for k in ("wet_fte", SHARED...)
    alloc = sum(r[k] for r in snap["resource_allocations"])
    @printf("  %-26s %5.0f / %-5.0f  %5.1f%%   [%s / %s]\n", res[k]["resource_name"], alloc, res[k]["total_capacity_units"], 100alloc / res[k]["total_capacity_units"], res[k]["allocated_units"], res[k]["total_capacity_units"])
end
burn = sum(r["cro_usd_m_month"] for r in snap["resource_allocations"]) +
    sum(sum(r[k] for r in snap["resource_allocations"]) * res[k]["cost_usd_m_per_unit_month"] for k in ("wet_fte", SHARED...))
@printf("monthly burn at the tick-42 allocation: %s → %s per quarter [deck: QTD burn \$142M of \$150M]\n", usd1(burn), usd1(3burn))

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§2  Compile to ReactiveDynamics model documents (eval-free JSON), validate, save")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

probe_registry = Dict{Symbol, Any}(:rebalance => (state, firing) -> nothing, :sync_status => (state, firing) -> nothing)
for arm in keys(ARMS)
    doc = compile_model(snap, arm)
    diags = validate(JSON.parse(JSON.json(doc)); registry = probe_registry)
    isempty(diags) || error("model for arm $arm failed validation:\n" * join(string.(diags), "\n"))
    path = joinpath(OUT, "model_$(arm).rdj.json")
    open(io -> JSON.print(io, doc, 2), path, "w")
    @printf(
        "  %-12s %-38s %2d places  %2d transitions  %2d arcs  %2d rules  %d AA input ports  validate: OK → %s\n", arm, ARMS[arm].label,
        length(doc["places"]), length(doc["transitions"]), length(doc["arcs"]), length(doc["rules"]), length(doc["inputs"]), relpath(path, @__DIR__)
    )
end
println("The four documents differ only in two params (kill_at) and the replan rule's guard: the arms are the decision, not four models.")

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§3  One coupled run: {portfolio net (RD), clinical ops (AA agent)} under one AlgebraicAgents root")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

demo_seed = UInt64(hash((ROOT_SEED, 1)))
sys = run_twin(demo_seed; arm = :closed_loop, snap)
rd = sys.rd
println("hierarchy: ", getname(sys.root), " ⊇ {", join(sort([getname(a) for a in values(inners(sys.root))]), ", "), "}  wires: ", length(getopera(sys.root).wires))

slow = build_twin(demo_seed; arm = :closed_loop, snap, fast = false)
simulate(slow.root)
println("compile cache leaves the trajectory identical: ", slow.rd.sol == rd.sol)

traj = RD.token_trajectory(rd)
removed_at(code) = (r = findfirst(i -> traj.program[i] == code && traj.place[i] === :removed, 1:nrow(traj)); r === nothing ? NaN : traj.t[r])
println("\nmonth-by-month (closed loop, seed 1):")
sat = study_saturation(rd)
for t in 0.0:3.0:24.0
    row = String[]
    for p in progs
        v = get(sat, p.id, Tuple{Float64, Float64}[])
        k = findfirst(x -> x[1] == t, v)
        push!(row, k === nothing ? "   ·  " : @sprintf("%5.0f%%", 100v[k][2]))
    end
    @printf("  m%-3.0f %s\n", t, join(row, " "))
end
println("        ", join([rpad(p.code, 6) for p in progs], " "), "   (study speed = allocator fill; · = no study in flight)")
println("\noutcomes:")
for p in progs
    tok = only(t for t in values(inners(getagent(rd, "structured"))) if t.code == Symbol(p.code))
    status = isfinite(killed_at(rd, p.id)) ? @sprintf("halted by kill switch at m%.0f", killed_at(rd, p.id)) :
        !isnan(tok.readout_t) ? @sprintf("readout PASSED at m%.0f → %s", tok.readout_t, p.next_stage) :
        tok.place === :removed ? @sprintf("readout FAILED at m%.0f → Terminated", removed_at(p.code)) : "study still running at m24"
    @printf("  %-8s %-13s %s\n", p.code, p.stage, status)
end
println("\nre-plans issued by the controller (wet-lab FTE per program):")
for (t, plan) in sys.replans
    @printf("  m%-3.0f %s\n", t, join(["$(c) $(Int(plan[c]))" for c in sort(collect(keys(plan)))], "  "))
end
println("clinical-ops contract violations rejected at the wire: ", length(sys.clin.violations))
@printf(
    "per-program ledger vs. aggregate: Σ program cost %.2f + unattributed %.2f = %.2f; aggregate spend %.2f\n",
    sum(RD.program_ledger(rd).cost_incurred), rd.unattributed_cost, sum(RD.program_ledger(rd).cost_incurred) + rd.unattributed_cost, spend(rd)
)

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§4  Ensembles: $(NSEED) seeded runs per arm, common seeds across arms (paired comparison)")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

ens = Dict{Symbol, Any}()
systems = Dict{Symbol, Vector{Any}}()
for arm in keys(ARMS)
    t0 = time()
    systems[arm] = Any[]
    ens[arm] = RD.ensemble(; nseed = NSEED, root_seed = ROOT_SEED, name = "ensemble_$(arm)") do seed
        s = run_twin(seed; arm, snap)
        push!(systems[arm], s)
        disentangle!(s.rd)   # RD.ensemble re-parents each member; detach it from the finished twin first
        s.rd
    end
    @printf("  %-12s %d members in %.1fs\n", arm, NSEED, time() - t0)
end

metric_value(m) = riskadj_value(m, progs; ε = EPS)
metric_spend(m) = spend(m)
metric_net(m) = net_value(m, progs; ε = EPS)
vals(arm, f) = Float64[f(m) for m in ens[arm].members]
paired(arm, f; base = :static) = (d = vals(arm, f) .- vals(base, f); (delta = mean(d), se = std(d) / sqrt(length(d))))

println("\n  arm                                     risk-adj eNPV        spend (24 mo)       net value          Δ net vs static (paired)")
for arm in keys(ARMS)
    v, s, n = RD.summarize(ens[arm], metric_value), RD.summarize(ens[arm], metric_spend), RD.summarize(ens[arm], metric_net)
    d = paired(arm, metric_net)
    @printf(
        "  %-38s %8.1f ± %-6.1f  %7.1f ± %-6.1f  %8.1f ± %-6.1f  %s\n", ARMS[arm].label, v.mean, v.sem, s.mean, s.sem, n.mean, n.sem,
        arm === :static ? "—" : @sprintf("%+6.1f ± %.1f  (%+.1f%%)", d.delta, d.se, 100d.delta / n.mean)
    )
end
te = RD.treatment_effect(ens[:static], ens[:closed_loop], metric_net)
head = paired(:closed_loop, metric_net)
@printf("\n  headline: closed loop vs static annual plan = %+.1f ± %.1f \$M net over 24 months (paired); unpaired SE from RD.treatment_effect: %.1f\n", head.delta, head.se, te.se)
@printf("  as a lift on the static arm's net value: %+.1f%%   [deck claims +14.8%%]\n", 100head.delta / mean(vals(:static, metric_net)))

kill_d, replan_d = paired(:kill_only, metric_net), paired(:replan_only, metric_net)
inter = head.delta - kill_d.delta - replan_d.delta
@printf("  decomposition: kill switch %+.1f, monthly re-plan %+.1f, interaction %+.1f\n", kill_d.delta, replan_d.delta, inter)
spend_d, value_d = paired(:closed_loop, metric_spend), paired(:closed_loop, metric_value)
@printf("  where it comes from: spend avoided %+.1f ± %.1f, risk-adjusted value protected %+.1f ± %.1f\n", -spend_d.delta, spend_d.se, value_d.delta, value_d.se)

println("\n  sensitivity to the patent-life erosion rate (eNPV lost per month of slip):")
sens = Any[]
for ε in (0.0, 0.0075, 0.015, 0.03, 0.05)
    d = paired(:closed_loop, m -> net_value(m, progs; ε))
    push!(sens, (ε, d.delta, d.se))
    @printf("    ε = %4.1f%%/month   Δ net %+6.1f ± %.1f\n", 100ε, d.delta, d.se)
end

realized_sd = std(vals(:closed_loop, m -> realized_value(m, progs; ε = EPS)))
@printf("\n  outcome noise: one run's realized eNPV has SD %.0f \$M across seeds, %.0fx the headline effect —\n", realized_sd, realized_sd / head.delta)
println("  a single simulated future cannot show the effect; only the ensemble expectation can.")

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§5  Write-back in NocoBase's shape: simulation_ticks, simulation_events, resource_allocations")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

anchor = Date(snap["anchor_date"])
function nocobase_rows(s, arm, scenario_id)
    rd, progs_ = s.rd, s.progs
    ts, u = utilization(rd, snap)
    cost_t, cost_v = spend_series(rd)
    ticks, events = Any[], Any[]
    names_ = trans_names(rd)
    for (i, t) in enumerate(ts)
        push!(
            ticks, Dict(
                "scenario_id" => scenario_id, "tick_index" => snap["tick_index"] + Int(t), "sim_date" => string(anchor + Month(Int(t))),
                "portfolio_enpv" => round(riskadj_value(rd, progs_; upto = t, ε = EPS); digits = 1),
                "snapshot_metrics" => Dict(
                    "spend_month_usd_m" => round(cost_v[i]; digits = 2), "spend_cum_usd_m" => round(sum(cost_v[1:i]); digits = 1),
                    "utilization" => Dict(k => round(v[i]; digits = 3) for (k, v) in u),
                    "readouts_passed" => rd.sol.readouts_passed[min(i + 1, nrow(rd.sol))],
                ), "is_synthetic" => true,
            )
        )
    end
    ev(t, asset, type, sev; payload = Dict()) = push!(events, Dict("scenario_id" => scenario_id, "tick" => snap["tick_index"] + Int(t), "asset" => asset, "event_type" => type, "severity" => sev, "payload" => payload, "is_synthetic" => true))
    traj_ = RD.token_trajectory(rd)
    for p in progs_
        tk = killed_at(rd, p.id)
        isfinite(tk) && ev(tk, p.code, "Kill switch: hERG IC50 $(p.herg) uM < 1.0 — study halted, capacity released", "CRITICAL")
        tok = only(t for t in values(inners(getagent(rd, "structured"))) if t.code == Symbol(p.code))
        !isnan(tok.readout_t) && ev(tok.readout_t, p.code, "$(p.stage) readout passed → $(p.next_stage)", "INFO")
        r = findfirst(i -> traj_.program[i] == p.code && traj_.place[i] === :removed, 1:nrow(traj_))
        r === nothing || ev(traj_.t[r], p.code, "$(p.stage) readout failed → Terminated", "WARNING")
    end
    for r in rd.log
        r[1] === :new_transitions || continue
        for (h, q) in r[3:end]
            q > 0 && startswith(names_[h], "gmp_batch_") && ev(r[2], "PRG-" * names_[h][12:end], "GMP synthesis batch started (Suite B, 3-month lead time)", "INFO")
        end
    end
    for p in progs_
        p.clinical || continue
        enrolled = cumsum(Float64[n for (t, c, n) in s.clin.history if c == p.code])
        for (i, e) in enumerate(enrolled)
            ratio = e / (p.need * i)
            if ratio < 0.85
                ev(i - 1, p.code, @sprintf("Accrual lag: %.0f%% of plan", 100ratio), "WARNING"; payload = Dict("trial_id" => p.trial_id))
                break
            end
        end
    end
    for p in progs_
        sl = get(slips(rd), p.id, 0.0)
        sl >= 1 && !isfinite(killed_at(rd, p.id)) && ev(24, p.code, @sprintf("Slip %.1f months vs. full-speed plan (patent-life erosion %.1f%% of eNPV)", sl, 100EPS * sl), sl >= 3 ? "WARNING" : "INFO")
    end
    sort!(events; by = e -> (e["tick"], e["asset"]))
    plans = [Dict("scenario_id" => scenario_id, "tick" => snap["tick_index"] + Int(t), "program_code" => c, "resource_id" => "wet_fte", "allocated_units" => q) for (t, plan) in s.replans for (c, q) in sort(collect(plan))]
    return ticks, events, plans
end

scenarios = Any[]
all_ticks, all_events, all_plans = Any[], Any[], Any[]
for (k, arm) in enumerate(keys(ARMS))
    sid = "SCN-$(lpad(k, 2, '0'))"
    push!(scenarios, Dict("scenario_id" => sid, "name" => ARMS[arm].label, "parent_id" => k == 1 ? nothing : "SCN-01", "params" => Dict("kill_at_month" => ARMS[arm].kill_at, "replan" => string(ARMS[arm].replan), "root_seed" => ROOT_SEED, "nseed" => NSEED), "is_synthetic" => true))
    tk, evs, pl = nocobase_rows(systems[arm][1], arm, sid)
    append!(all_ticks, tk); append!(all_events, evs); append!(all_plans, pl)
end
summary = Dict(
    "nseed" => NSEED, "root_seed" => ROOT_SEED, "erosion_per_month" => EPS,
    "arms" => Dict(
        string(arm) => Dict(
            "label" => ARMS[arm].label,
            "riskadj_enpv" => RD.summarize(ens[arm], metric_value), "spend" => RD.summarize(ens[arm], metric_spend), "net" => RD.summarize(ens[arm], metric_net),
            "delta_net_vs_static" => paired(arm, metric_net),
        ) for arm in keys(ARMS)
    ),
    "sensitivity" => [Dict("erosion" => e, "delta" => d, "se" => s) for (e, d, s) in sens],
)
# One row per (arm, seed): the numbers behind every ensemble statistic, so a reader can recompute
# them, chart the spread, or slice by program without re-running. Seeds are shared across arms.
per_seed = Any[]
for arm in keys(ARMS), (k, m) in enumerate(ens[arm].members)
    sl = slips(m)
    toks = Dict(String(t.code) => t for t in values(inners(getagent(m, "structured"))))
    push!(
        per_seed, Dict(
            "arm" => string(arm), "seed_index" => k, "seed" => string(ens[arm].seeds[k]),
            "riskadj_enpv" => round(metric_value(m); digits = 3), "spend" => round(metric_spend(m); digits = 3),
            "net" => round(metric_net(m); digits = 3), "realized_enpv" => round(realized_value(m, progs; ε = EPS); digits = 3),
            "programs" => [
                Dict(
                    "program_code" => p.code, "slip_months" => round(get(sl, p.id, 0.0); digits = 3),
                    "killed_at" => isfinite(killed_at(m, p.id)) ? killed_at(m, p.id) : nothing,
                    "readout" => isfinite(killed_at(m, p.id)) ? "halted" : !isnan(toks[p.code].readout_t) ? "passed" : toks[p.code].place === :removed ? "failed" : "pending",
                    "readout_t" => isnan(toks[p.code].readout_t) ? nothing : toks[p.code].readout_t,
                ) for p in progs
            ],
        )
    )
end
writejson(name, obj) = open(io -> JSON.print(io, obj, 2), joinpath(OUT, name), "w")
writejson("scenarios.json", scenarios)
writejson("simulation_ticks.json", all_ticks)
writejson("simulation_events.json", all_events)
writejson("resource_allocations_plan.json", all_plans)
writejson("ensemble_summary.json", summary)
writejson("ensemble_per_seed.json", per_seed)
println("  scenarios.json                 $(length(scenarios)) rows (scenario_id scopes every row below; S01 is the parent)")
println("  simulation_ticks.json          $(length(all_ticks)) rows (tick_index, sim_date, portfolio_enpv, snapshot_metrics)")
println("  simulation_events.json         $(length(all_events)) rows (kills, readouts, GMP batches, accrual lags, slips)")
println("  resource_allocations_plan.json $(length(all_plans)) rows (the controller's wet-lab FTE plan per re-plan)")
println("  ensemble_summary.json          the §4 numbers, for the dashboard")
println("  ensemble_per_seed.json         $(length(per_seed)) rows (one per arm × seed: value, spend, net, per-program slip and readout)")
println("\n  closed-loop event feed (member 1):")
for e in filter(e -> e["scenario_id"] == "SCN-04", all_events)
    @printf("    T-%03d  %-8s %-9s %s\n", e["tick"], e["asset"], e["severity"], e["event_type"])
end

# ═════════════════════════════════════════════════════════════════════════════════════════════════
banner("§6  Figures → out/")
# ═════════════════════════════════════════════════════════════════════════════════════════════════

# (a) the model itself: one program's slice of the Petri net (the full net is 6 copies of it plus shared pools)
# These four drawings need Graphviz's `dot` on the PATH (https://graphviz.org). Without it they are
# skipped, and every number and chart below is still produced.
if Sys.which("dot") === nothing
    println("  Graphviz `dot` not found: skipping net_one_program.svg, net_full.svg, net_exec_map.svg, aa_wiring.svg")
else
    one = deepcopy(snap)
    one["therapeutics_programs"] = filter(r -> r["program_code"] == "PRG-518", one["therapeutics_programs"])
    one_sys = build_twin(demo_seed; arm = :closed_loop, snap = one)
    RD.draw_network(one_sys.rd; path = joinpath(OUT, "net_one_program.svg"))
    RD.draw_network(rd; path = joinpath(OUT, "net_full.svg"))
    RD.exec_map(rd; path = joinpath(OUT, "net_exec_map.svg"))
    AlgebraicAgents.run_graphviz(joinpath(OUT, "aa_wiring.svg"), AlgebraicAgents.wiring_diagram(sys.root))
    println("  net_one_program.svg, net_full.svg, net_exec_map.svg (Graphviz via RD.draw_network / RD.exec_map)")
    println("  aa_wiring.svg (AlgebraicAgents.wiring_diagram: the hierarchy and its $(length(getopera(sys.root).wires)) wires)")
end

# (b) dashboard charts (agent/scripts/chartkit.jl: dependency-free SVG)
bars = [(label = ARMS[arm].label, value = paired(arm, metric_net).delta, se = paired(arm, metric_net).se) for arm in (:kill_only, :replan_only, :closed_loop)]
write(joinpath(OUT, "chart_lever_effects.svg"), ChartKit.bars_error(bars; unit = "\$M", title = "What each lever is worth", subtitle = "Δ net value vs. static annual plan, 24 months, $(NSEED) paired seeds, ±2 SE", x_label = "Δ net value (\$M)", sortby = "none"))

write(
    joinpath(OUT, "chart_value_bridge.svg"), ChartKit.waterfall(
        (label = "Static plan", value = mean(vals(:static, metric_net))),
        [(label = "Kill switch", value = kill_d.delta, se = kill_d.se), (label = "Monthly re-plan", value = replan_d.delta, se = replan_d.se), (label = "Interaction", value = inter)],
        (label = "Closed loop", value = mean(vals(:closed_loop, metric_net)));
        unit = "\$M", title = "Net value bridge: static plan → closed loop", subtitle = "Risk-adjusted eNPV minus 24-month spend, mean of $(NSEED) seeds", y_label = "\$M"
    )
)

write(
    joinpath(OUT, "chart_paired.svg"), ChartKit.paired_dots(
        [(label = "Static annual plan", samples = vals(:static, metric_net)), (label = "Closed loop", samples = vals(:closed_loop, metric_net))];
        delta = head.delta, delta_se = head.se, unit = "\$M", title = "Same futures, two policies", subtitle = "Net value per seed; each seed is one simulated future run under both policies", x_label = "net value (\$M)"
    )
)

ts_s, u_s = utilization(systems[:static][1].rd, snap)
ts_c, u_c = utilization(systems[:closed_loop][1].rd, snap)
write(
    joinpath(OUT, "chart_wetlab.svg"), ChartKit.trace_ceiling(
        ts_s, [(label = "Static plan", y = 100 .* u_s["wet_fte"]), (label = "Closed loop", y = 100 .* u_c["wet_fte"], emphasis = true)];
        ceiling = 100.0, ceiling_label = "52 FTE", x_label = "month", y_label = "wet-lab FTE in use (%)", unit = "%",
        title = "Wet-lab FTE: freed capacity is re-used", subtitle = "Share of the 52 research FTE working on a live study, member 1"
    )
)
write(
    joinpath(OUT, "chart_gmp.svg"), ChartKit.trace_ceiling(
        ts_c, [(label = "GMP Suite B", y = 100 .* u_c["gmp"], emphasis = true), (label = "Mass spec", y = 100 .* u_c["ms"]), (label = "Vivarium", y = 100 .* u_c["viv"])];
        ceiling = 100.0, ceiling_label = "capacity", x_label = "month", y_label = "utilization (%)", unit = "%",
        title = "Shared pools under the closed loop", subtitle = "Monthly utilization, member 1; GMP batches compete with studies for Suite B"
    )
)
for f in ("chart_lever_effects", "chart_value_bridge", "chart_paired", "chart_wetlab", "chart_gmp")
    println("  $f.svg")
end
println("\nDone. Re-run with --nseed=N to change the ensemble size.")
