---
name: rd-study
description: Run a decision study on a resource-constrained, timed process using ReactiveDynamics.jl — answer questions like "what is the fifth scientist worth", "where is the bottleneck", "should we take this deal", "how much capacity do we need". Use when someone needs a defensible number out of a system where something binds, something takes time, and something is decided mid-flight, and a spreadsheet cannot produce it. Orchestrates four phases: frame the decision, author a model document, prove it plausible, report the answer.
---

# rd-study — a decision study, end to end

You are running a **study**, not building software. The deliverable is a number with an error bar, a decision rule that uses it, and enough evidence that someone can act on it. The model is instrumentation.

## The one qualifying test

Before anything else, check the question against three conditions. A simulation earns its cost only when **all three** hold:

1. **Something binds** — a finite resource that runs out (headcount, bench, line capacity, cash, reviewers).
2. **Something takes time** — work occupies that resource for a while before yielding anything.
3. **Something is decided mid-flight** — when demands collide, something has to choose who waits.

Two out of three is a spreadsheet. Say so plainly and offer the arithmetic instead — you will have saved the user a week. Three out of three and the answer is a *marginal response to relieving a constraint*, which has no closed form; that is the study.

Run this test out loud in your first response. It is the highest-value thing in this skill.

## The four phases

Each has its own skill. Read the phase skill when you enter the phase, not before — they are long and only one applies at a time.

| Phase | Skill | Produces | Done when |
|---|---|---|---|
| 1. Frame | `rd-frame` | `BRIEF.md` | every slot is filled and the failure modes are pre-registered |
| 2. Author | `rd-author` | `model.rdj.json` (+ one file per arm) | `rdx validate` returns zero diagnostics |
| 3. Diagnose | `rd-doctor` | `DOCTOR.md` | `rdx diagnose` exits 0, or every failure is understood and stated |
| 4. Report | `rd-report` | `FINDINGS.md` | the headline, the decision rule, and the caveats are all in it |

**The gates are real.** Do not author before the brief's slots are filled — you will interview the user twice. Do not quote a number before DIAGNOSE passes — a schema-clean model can be confident nonsense, and the failure mode of this whole workflow is a beautifully formatted wrong answer. If DIAGNOSE fails, go back to AUTHOR; the brief's pre-registered failure modes exist precisely so that this bounce is a short loop and not a redesign.

## Layout

One directory per study, under `agent/studies/<slug>/` in this repo or anywhere the user prefers:

```
<slug>/
  BRIEF.md              phase 1
  host.jl               only if the study uses structured tokens
  model.rdj.json        phase 2 — the baseline arm
  arm_<lever>.rdj.json  phase 2 — one file per counterfactual, differing in as few lines as possible
  DOCTOR.md             phase 3
  FINDINGS.md           phase 4
  figures/              exec map, result plots
```

`agent/studies/marginal_scientist/` is the worked exemplar and the workflow's known-answer test. When unsure what an artifact should look like, read its counterpart there.

## The tool

Everything mechanical goes through one CLI, from the repository root:

```bash
julia --project=. agent/scripts/rdx.jl <command> [args] [--flag=value ...]
```

`validate` · `run` · `diagnose` · `ensemble` · `effect` · `draw`. Exit code `0` ok · `1` findings · `2` error. **Branch on the exit code and repair against `findings[]`** — every finding carries a JSON path into the document, which is what makes a failure mechanically fixable. Never parse the prose. `agent/README.md` has the full command and metric tables.

Two standing rules about the tool:

- **Never hand-edit a number into a report.** Every figure in FINDINGS must come from a command you actually ran, and the command goes in the report. If you cannot produce a number by running something, it does not go in.
- **Determinism is contractual.** A run is fully determined by `(model, seed)`. Always pass `--seed` / `--root-seed` explicitly and record it, so anyone can replay the exact number rather than one near it.

## Working style

Interview in the user's vocabulary — headcount, trials, queues, backlog — never in Petri-net vocabulary. The translation into places and arcs is your job, and it happens in phase 2 where the user reviews the result. A user who has to learn the formalism to answer your questions has been failed by the workflow.

Tell the user where you are (`phase 2 of 4`) and what you still need. Keep the artifacts on disk as you go; they are the interface, and a study that dies mid-conversation should be resumable from its directory alone.

State uncertainty where it exists. A warning carried into the report as a scope caveat is worth more than a warning suppressed, and far more than a warning you did not look for.

## When to stop

The study is done when the user has a number they can defend to someone who was not in the room. Not when the model is elegant, not when every warning is cleared — when the decision is supported. If the honest answer is "the effect is inside the noise" or "four scientists are plenty, the fifth is worth their cost", that IS the finding. Report it.
