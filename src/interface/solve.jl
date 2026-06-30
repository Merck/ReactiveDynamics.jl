export @agentize

# NOTE (ADR 0014 §A1 / CONTRACT §15.1): this file previously carried ~195 lines of dead SciML
# plotting — `plot_summary`/`plot_ensemble_sol`/`first_sol`/`plot_from_log` and the `@plot` macro —
# all written against SciML `EnsembleSummary`/`EnsembleSolution` objects that ADR 0001 removed when
# it demoted SciML. Those types were never imported, so the macro was dead-on-arrival on
# `ref-agents`. The result-plotting story is now the model-agnostic `@recipe` set in the `RDPlotsExt`
# package extension (ext/RDPlotsExt.jl); the live generic `_draw` reduction moved there too (it
# needs `Plots`, which is now a weakdep). Nothing here imports `Plots` anymore.
