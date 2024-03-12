using ReactiveDynamics

include("safeinclude.jl")

@safeinclude "example" "../tutorial/example.jl"
@safeinclude "joins" "../tutorial/joins/joins.jl"
@safeinclude "loadsave" "../tutorial/loadsave/loadsave.jl"
# @safeinclude "optimize" "../tutorial/optimize/optimize.jl"
# @safeinclude "solution wrap" "../tutorial/optimize/optimize_custom.jl"
@safeinclude "toy pharma model" "../tutorial/toy_pharma_model.jl"
