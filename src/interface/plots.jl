# Result plotting (ADR 0014 / CONTRACT §15.1).
#
# All result-plot code now lives in the `RDPlotsExt` package extension (ext/RDPlotsExt.jl), loaded
# only when the user has `Plots` available — `Plots` is a weakdep (ADR 0014 §A2 / open question:
# author recipes in the extension from the start so the only choice is `[deps]` vs `[weakdeps]`).
# That extension carries the seven model-agnostic `@recipe`s AND the generic `_draw` reduction
# (`AlgebraicAgents._draw(::ReactionNetworkProblem, vars)`, place trajectories from `prob.sol`),
# which was the one live plot here. The never-called `plot_df` helper was removed (ADR 0014 §A1).
#
# This file is intentionally (almost) empty: it keeps no `using Plots` in the core so the package
# precompiles and loads with `Plots` absent. See ext/RDPlotsExt.jl for the recipes and `_draw`.
