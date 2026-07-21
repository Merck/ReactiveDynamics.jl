# BD acquisition-impact demo — export the canonical ensemble to JSON.
#
# Runs the SAME S0–S5 grid as run_demo.jl / figures.jl (160 seeds, root 2026) and writes the
# per-seed metric vectors + summary stats to presentation_data.json — the real numbers behind the
# result figures, suitable for downstream client-side rendering (no faked distributions).
# Dependency-free JSON emission (manual string building) to avoid a JSON dep.
#
# Run:  julia --project=. demo/bd_acquisition/export_data.jl

using Printf, Statistics

const HERE = @__DIR__
include(joinpath(HERE, "host.jl"))
include(joinpath(HERE, "analysis.jl"))
include(joinpath(HERE, "run_demo.jl"))   # main() is PROGRAM_FILE-guarded, so it does not auto-run

jarr(v) = "[" * join((@sprintf("%.4f", x) for x in v), ",") * "]"

function export_data(; root_seed = 2026, nseed = 160)
    scenarios = [:S0, :S1, :S2, :S3, :S4, :S5]
    @info "Running S0–S5 grid for export ($nseed seeds)…"
    # Each scenario is an engine EnsembleProblem (RD.ensemble via scenario_ensemble); per-member
    # vectors and summary stats come from the engine-backed analysis helpers. The acquisition price
    # is netted in the metric closures (deals carry it, S0 does not).
    R = Dict(
        s => scenario_ensemble(seed -> run_scenario(s, seed); root_seed = root_seed, nseed = nseed)
            for s in scenarios
    )
    base = R[:S0]
    price(s) = s == :S0 ? 0.0 : ACQ_PRICE
    out = joinpath(HERE, "presentation_data.json")
    open(out, "w") do io
        println(io, "{")
        println(io, "  \"nseed\": $nseed, \"root_seed\": $root_seed, \"acq_price\": $(ACQ_PRICE),")
        println(io, "  \"scenarios\": {")
        for (i, s) in enumerate(scenarios)
            ens = R[s]
            te = treatment_effect(base, ens; deal_price = price(s))
            comma = i < length(scenarios) ? "," : ""
            println(io, "    \"$s\": {")
            @printf(
                io, "      \"mean_rnpv\": %.4f, \"sem_rnpv\": %.4f,\n",
                mean_rnpv(ens; acq_price = price(s)), sem_rnpv(ens; acq_price = price(s))
            )
            @printf(io, "      \"mean_launches\": %.4f, \"p_launch\": %.4f,\n", mean_launches(ens), p_launch(ens))
            @printf(
                io, "      \"mean_cash_trough\": %.4f, \"mean_sci_trough\": %.4f,\n",
                mean_cash_trough(ens), mean_sci_trough(ens)
            )
            @printf(
                io, "      \"delta_rnpv\": %.4f, \"se_delta_rnpv\": %.4f, \"delta_launches\": %.4f,\n",
                te.delta_rnpv, te.se_delta_rnpv, te.delta_launches
            )
            println(io, "      \"rnpv\": ", jarr(rnpv_samples(ens; acq_price = price(s))), ",")
            println(io, "      \"launches\": ", jarr(launch_samples(ens)), ",")
            println(io, "      \"cash_trough\": ", jarr(cash_trough_samples(ens)), ",")
            println(io, "      \"sci_trough\": ", jarr(sci_trough_samples(ens)))
            println(io, "    }$comma")
        end
        println(io, "  }")
        println(io, "}")
    end
    println("wrote ", out)
    return R
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_data()
end
