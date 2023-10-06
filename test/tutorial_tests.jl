using ReactiveDynamics

@safetestset "example" "../tutorial/example.jl"
@safetestset "joins" "../tutorial/joins/joins.jl"
@safetestset "loadsave" "../tutorial/loadsave/loadsave.jl"
# @safeinclude "optimize" "../tutorial/optimize/optimize.jl"
# @safeinclude "solution wrap" "../tutorial/optimize/optimize_custom.jl"
@safetestset "toy pharma model" "../tutorial/toy_pharma_model.jl"
