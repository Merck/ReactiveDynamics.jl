# `agent/` — the agent-facing surface of ReactiveDynamics.jl

This directory is the portable skill/plugin layer that carries a user from a business question to a defensible number: a decision frame, a validate-clean model document, a plausibility gate, and a briefing whose every figure is replayable from `(model, seed)`.

It is chartered in **[`spec/AGENT_WORKFLOW_CHARTER.md`](../spec/AGENT_WORKFLOW_CHARTER.md)** — read that for the why, the four phases, the layer division (skill vs MCP vs subagent vs hook), and the tracker. This file is the operational entry point.

Why it lives in-repo rather than in a separate plugin repository: the skills' reference material is *generated from this source tree* and guarded by the same drift gate as the docs. A skill that drifts from the engine is worse than no skill.

## Layout

| Path | What it is |
|---|---|
| `scripts/rdx.jl` | the CLI — one entry point per phase, JSON on stdout, meaningful exit codes |
| `scripts/diagnose.jl` | the plausibility checks (`rdx diagnose`), the part that carries modeling judgement |
| `skills/` | the skills themselves (`rd-study` orchestrator + one per phase) |
| `reference/` | generated extracts of the authoring schema — never hand-edited |
| `studies/` | worked studies; `marginal_scientist/` is the known-answer acceptance test |

## The CLI

Everything runs through `rdx`, from the repository root:

```bash
julia --project=. agent/scripts/rdx.jl <command> [args] [--flag=value ...]
```

| Command | Does |
|---|---|
| `validate MODEL` | static, eval-free schema check → diagnostics with JSON paths |
| `run MODEL --seed=N` | one seeded run → per-place troughs, per-transition spawns/completions, ledger totals; `--out=DIR` writes an export bundle |
| `diagnose MODEL --seed=N` | the plausibility gate — nine checks against the brief's declared intent |
| `ensemble MODEL --seeds=N --metric=M` | seeded ensemble → mean ± sem, quartiles |
| `effect BASE VARIANT --metric=M --seeds=N` | two-arm treatment effect Δ ± se, with a significance verdict |
| `draw MODEL --out=PATH` | the net; `--exec-map --seed=N` decorates it with run statistics |

**Exit codes carry the meaning**: `0` ok · `1` findings (invalid model, failing checks, an effect inside the noise) · `2` error. An agent branches on the code and repairs against `findings[]`; it never has to parse prose.

**Metrics are named, not code** — `final:<place>`, `min:<place>`, `max:<place>`, `sum:<place>`, `spawns[:<transition>]`, `completions[:<t>]`, `successes[:<t>]`, `reward`, `cost`, `net`, `enpv`. This is deliberate: letting the caller pass a Julia reduction would reopen the code path the JSON artifact exists to close. `enpv` takes `--annual`, `--periods-per-year` and `--scale`.

**`--host=FILE` is the one place authored Julia executes.** A structured token's Julia *type* is the only thing a model document cannot carry, so a study that uses tokens ships a `host.jl` defining the type and a registry mapping kind names to constructors. `rdx` prints a notice to stderr whenever it includes one. Models without structured tokens need no host file at all.

## Worked example — the known-answer study

`studies/marginal_scientist/` reproduces the published case study's shadow price (≈ +\$19M ± \$2.5M) entirely through the JSON path. It is the acceptance test for the whole sequence, and the reference for what each phase's artifact should look like: `BRIEF.md` (frame) → `model.rdj.json` + `arm_hire_fifth.rdj.json` (author) → `DOCTOR.md` (diagnose) → `FINDINGS.md` (report).

```bash
S=agent/studies/marginal_scientist
julia --project=. agent/scripts/rdx.jl diagnose $S/model.rdj.json --host=$S/host.jl \
  --seed=1 --expect-binds=scientist --expect-reusable=scientist --tally=launch
julia --project=. agent/scripts/rdx.jl effect $S/model.rdj.json $S/arm_hire_fifth.rdj.json \
  --host=$S/host.jl --metric=enpv --scale=0.001 --seeds=48 --root-seed=2026
```

The second command prints `delta: 18.991, se: 2.528`. If it ever prints something else, the workflow has regressed against the engine — that is what makes it a test rather than a demo.

## Portability

The CLI is plain Julia over Bash, and the skills are Markdown. Nothing here depends on Claude Code: any agent that can run a shell command and read a file can drive this workflow. Claude Code's skill/plugin packaging is the *distribution* mechanism, not the substrate.
