# BD acquisition-impact demo — export the canonical ensemble to JSON for the HTML brief.
#
# Runs the SAME S0–S5 grid as run_demo.jl / figures.jl (160 seeds, root 2026) and writes the
# per-seed metric vectors + summary stats to presentation_data.json, which is inlined into the
# static HTML presentation (figures rendered client-side as SVG from the real numbers — no faked
# distributions). Dependency-free JSON emission (manual string building) to avoid a JSON dep.
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
    R = Dict(
        s => ensemble(seed -> run_scenario(s, seed); root_seed = root_seed, nseed = nseed,
                      acq_price = (s == :S0 ? 0.0 : ACQ_PRICE)) for s in scenarios
    )
    base = R[:S0]
    out = joinpath(HERE, "presentation_data.json")
    open(out, "w") do io
        println(io, "{")
        println(io, "  \"nseed\": $nseed, \"root_seed\": $root_seed, \"acq_price\": $(ACQ_PRICE),")
        println(io, "  \"scenarios\": {")
        for (i, s) in enumerate(scenarios)
            ms = R[s]
            te = treatment_effect(base, ms)
            comma = i < length(scenarios) ? "," : ""
            println(io, "    \"$s\": {")
            @printf(io, "      \"mean_rnpv\": %.4f, \"sem_rnpv\": %.4f,\n", mean_rnpv(ms), sem_rnpv(ms))
            @printf(io, "      \"mean_launches\": %.4f, \"p_launch\": %.4f,\n", mean_launches(ms), p_launch(ms))
            @printf(io, "      \"mean_cash_trough\": %.4f, \"mean_sci_trough\": %.4f,\n",
                    mean_cash_trough(ms), mean_sci_trough(ms))
            @printf(io, "      \"delta_rnpv\": %.4f, \"se_delta_rnpv\": %.4f, \"delta_launches\": %.4f,\n",
                    te.delta_rnpv, te.se_delta_rnpv, te.delta_launches)
            println(io, "      \"rnpv\": ", jarr([m.rnpv for m in ms]), ",")
            println(io, "      \"launches\": ", jarr([m.launches for m in ms]), ",")
            println(io, "      \"cash_trough\": ", jarr([m.cash_trough for m in ms]), ",")
            println(io, "      \"sci_trough\": ", jarr([m.sci_trough for m in ms]))
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
