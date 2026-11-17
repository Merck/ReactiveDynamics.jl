# Handoff prompt

Paste the block below into your coding agent (Claude Code, or any agent that can run shell commands), started in this folder.

---

You are working in a folder that contains a runnable simulation of a synthetic pharma R&D portfolio (a fictional company; all data synthetic, and no result is a forecast for a real portfolio), a "digital twin" built from a NocoBase snapshot with ReactiveDynamics.jl and AlgebraicAgents.jl. The material was built for a talk on using simulation to answer a strategic portfolio question: is a closed-loop policy (stop a study as soon as a tox flag is raised, and re-plan wet-lab staff monthly) worth it against a static annual plan? The talk page is `talk.html`; read it first to see the story and the numbers.

Then read `AGENTS.md`. It describes the files, how to run the model, which inputs you can change, and the rules for reporting a number. Follow those rules.

Do this first:

1. Install and check the setup: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`, then a quick run with `julia --project=. run_twin.jl --nseed=5 --out=scratch/smoke`.
2. Reproduce the baseline into a separate folder: `julia --project=. run_twin.jl --out=scratch/baseline`. Confirm the closed loop's delta in `scratch/baseline/ensemble_summary.json` is +53.3 ± 1.7. If it is not, stop and tell me, with the output of `julia --version` and `Pkg.status()`.

Then help me build my own views. Start with these, and ask me before going further:

- At what erosion rate (eNPV lost per month of slip) would the closed loop reach the prototype's +14.8%?
- What is the kill switch worth if the static plan acts on a tox flag at month 3 or 6 instead of month 12?
- How large a kit batch, or how short a lead time, would stop manufacturing from limiting PRG-518?

For each view, give me: the question, the inputs you changed (field and value), the result as a paired delta with its standard error, one sentence on what it means for the decision, and the assumption it depends on most. Make charts as standalone SVG or HTML files in `views/`. Keep every what-if in its own folder under `runs/`, and never write into `out/`, which is the reference the talk page quotes.

Write for a business audience: say what a number means for the decision before saying how it was computed, and avoid model jargon (places, transitions, arcs, tokens) unless I ask for it.
