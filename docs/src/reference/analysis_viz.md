```@meta
CurrentModule = ReactiveDynamics
```

# Analysis & visualization

The read-only inspection layer over a finished run: per-token trajectory helpers, the ensemble runner and its summaries, treatment-effect comparisons across policy arms, the results-export bundle, and the three-layer network exec map. All operate over run state produced by [`ReactionNetworkProblem`](@ref); none mutate the model.

```@docs
token_trajectory
representative_token
trajectory_envelope
ensemble
summarize
treatment_effect
EnsembleProblem
export_run
export_ensemble
network_graph
to_graphviz
draw_network
exec_map
NetworkGraph
```

## Plot specs

Typed plot specifications consumed by the `RDPlotsExt` recipes (loaded when `Plots` is present). Each wraps a run (or ensemble) and the variables to render.

```@docs
MarkingPlot
SaturationPlot
ValuationPlot
LedgerPlot
TokenTrajectoryPlot
EnsembleBar
TreatmentEffectPlot
ThroughputPlot
```

The `@export_solution_as_table` / `@export_solution_as_csv` macros are the older, solution-output path predating the [`export_run`](@ref) bundle; prefer `export_run`/`export_ensemble` for new work.
