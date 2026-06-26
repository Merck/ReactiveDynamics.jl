using SafeTestsets, BenchmarkTools

@time begin
    # Semantic suite: real assertions over operational semantics (contract + ADRs).
    # See test/semantic/runtests.jl for the tier scheme (T1 characterization / T2 acceptance).
    # (The former assertion-free tutorial smoke tests were retired; the worked examples now live
    # under demo/ — see demo/core_engine_tour/ and demo/agentic_pipeline/.)
    @time @safetestset "Semantic tests" begin
        include("semantic/runtests.jl")
    end
end
