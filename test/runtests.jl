using SafeTestsets, BenchmarkTools

@time begin
    @time @safetestset "Tutorial tests" begin include("tutorial_tests.jl") end
end
