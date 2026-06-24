using ReactiveDynamics

include("safeinclude.jl")

@safeinclude "example" "../tutorial/example.jl"
@safeinclude "joins" "../tutorial/joins/joins.jl"
# NOTE: the "loadsave" tutorial is RETIRED (ADR 0005): it exercised the legacy TOML/CSV loader,
# whose import-time eval of attribute strings + `registered` function bodies was an RCE. The model
# format is now the single eval-free model.rdj.json — see demo/bd_acquisition/model.rdj.json and
# @import_model/@export_model. The old tutorial/loadsave files remain for reference but are no
# longer loadable.
# @safeinclude "loadsave" "../tutorial/loadsave/loadsave.jl"
# @safeinclude "optimize" "../tutorial/optimize/optimize.jl"
# @safeinclude "solution wrap" "../tutorial/optimize/optimize_custom.jl"
@safeinclude "toy pharma model" "../tutorial/toy_pharma_model.jl"
