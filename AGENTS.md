# AGENTS.md

Orientation for a coding agent working in this repository: Codex, Jules, Cursor, Gemini CLI, or anything else that reads an `AGENTS.md`. Claude Code reads [`CLAUDE.md`](CLAUDE.md) instead, which carries the same facts plus the engine's internals.

## What this is

ReactiveDynamics.jl is a Julia engine for timed, stochastic, resource-constrained process models: business and R&D processes, not chemistry, despite the reaction-network notation. [`agent/`](agent) is a skill layer on top of it that carries a user from a business question to a defensible number. Models are JSON documents (`*.rdj.json`), never Julia source; the eval-free `validate` pass is the trust boundary.

## Start here

1. Read [`agent/README.md`](agent/README.md), then [`agent/skills/rd-study/SKILL.md`](agent/skills/rd-study/SKILL.md), the orchestrator. It runs four phases (FRAME → AUTHOR → DIAGNOSE → REPORT) and reads one phase skill at a time.
2. A skill named in the text (`rd-frame`, `rd-author`, `rd-doctor`, `rd-report`) lives at `agent/skills/<name>/SKILL.md`. Skills are plain Markdown playbooks: the YAML header is metadata, and everything they ask you to do is a shell command.
3. The authoring vocabulary is [`agent/reference/model_schema.md`](agent/reference/model_schema.md), generated from the source and never edited by hand. [`agent/studies/portfolio_bid/model.rdj.json`](agent/studies/portfolio_bid/model.rdj.json) is a pure-data model to copy as a starting point.

## Setup

- Julia ≥ 1.12 ([`juliaup`](https://github.com/JuliaLang/juliaup) is the easiest install), then once: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`. Expect several minutes on a fresh machine, mostly precompilation. No `Manifest.toml` is committed, so dependency versions resolve fresh; CI tests both the latest compatible versions and the declared lower bounds.
- Optional: Graphviz (`dot`) for `rdx draw`, which writes a `.dot` file without it; Python ≥ 3.8, standard library only, for the `launch_supply` helper scripts.
- The commands below are bash. On Windows use WSL, or run the same `julia` invocations from PowerShell; Julia itself is cross-platform.

## The CLI

Every action goes through one command, run from the repository root: `julia --project=. agent/scripts/rdx.jl <validate|run|diagnose|ensemble|effect|draw> MODEL [--flag=value ...]`. It prints one JSON envelope on stdout and notes on stderr. Exit `0` = ok, `1` = findings (an invalid model, a failing check, an effect inside the noise), `2` = error. Branch on the exit code, not on prose.

The repair loop: edit the JSON, run `rdx validate`, fix whatever sits at each `findings[].path`, repeat until exit 0; then `rdx run MODEL --seed=1` and `rdx diagnose MODEL --seed=1` with the expectations from the brief. Finding paths are Julia-style and **1-based**: `transitions[2]` is the second transition, `[1]` in Python or jq.

Each call pays Julia's startup, roughly 10–30 s. To amortize it, `include("agent/scripts/rdx.jl")` once in a `julia --project=.` session and call `main(["validate", "path/to/model.rdj.json"])`: it prints the same envelope and returns the exit code. `ensemble` and `effect` recompile the model for every seed (charter W24), so keep seed counts small while prototyping, and run one Julia process at a time — concurrent processes on one package depot block each other.

## Rules

- Author models as JSON documents only. A study's `host.jl` (a token type plus a registry, passed with `--host=`) is the one sanctioned place for Julia code.
- Do not edit `src/` to make a study work. Known engine gaps are worked around in `agent/` and tracked in [`spec/AGENT_WORKFLOW_CHARTER.md`](spec/AGENT_WORKFLOW_CHARTER.md) §12.
- Quote no number without its error bar: an effect inside ±2·se is noise, and `rdx effect` exits 1 on it. Before a write-up ships, have a second, fresh-context session critique it — the usual failure is in the narrative around the numbers, not in the numbers.
- After changing the schema, the whitelists or the action verbs, regenerate the reference with `julia --project=. agent/scripts/gen_reference.jl`; CI runs it with `--check`.
