# NOTE (WS-4 housekeeping): the `export @agentize` that used to head this file was DELETED — the
# macro was never defined anywhere (a dangling export; the export was the only `agentize` token in
# the repo). Agentization is implicit in the `ReactionNetworkProblem` constructor — the real public
# contract is `prob = ReactionNetworkProblem(net[, u0, p]; …)` then `simulate(prob[, n])` (ADR 0001).
# AA `@agentize` sugar over the constructor remains explicitly-deferred future work per ADR 0012.
#
# NOTE (ADR 0014 §A1 / CONTRACT §15.1): this file previously carried ~195 lines of dead SciML
# plotting — `plot_summary`/`plot_ensemble_sol`/`first_sol`/`plot_from_log` and the `@plot` macro —
# all written against SciML `EnsembleSummary`/`EnsembleSolution` objects that ADR 0001 removed when
# it demoted SciML. Those types were never imported, so the macro was dead-on-arrival on
# `ref-agents`. The result-plotting story is now the model-agnostic `@recipe` set in the `RDPlotsExt`
# package extension (ext/RDPlotsExt.jl); the live generic `_draw` reduction moved there too (it
# needs `Plots`, which is now a weakdep). Nothing here imports `Plots` anymore.
