using SafeTestsets, BenchmarkTools

@time begin
    # Phase-0 semantic suite: real assertions over operational semantics (contract + ADRs).
    # See test/semantic/runtests.jl for the tier scheme (T1 characterization / T2 acceptance).
    @time @safetestset "Semantic tests" begin
        include("semantic/runtests.jl")
    end

    # Legacy tutorial smoke tests: assert only "did not throw" (kept for regression coverage).
    @time @safetestset "Tutorial tests" begin
        include("tutorial_tests.jl")
    end
end
