# `agent/` — the agent-facing surface of ReactiveDynamics.jl

This directory is the portable skill/plugin layer that carries a user from a business question to a defensible number: a decision frame, a validate-clean model document, a plausibility gate, and a briefing whose every figure is replayable from `(model, seed)`.

It is chartered in **[`spec/AGENT_WORKFLOW_CHARTER.md`](../spec/AGENT_WORKFLOW_CHARTER.md)** — read that for the why, the four phases, the layer division (skill vs MCP vs subagent vs hook), and the tracker. This file is the operational entry point.

Why it lives in-repo rather than in a separate plugin repository: the skills' reference material is *generated from this source tree* and guarded by the same drift gate as the docs. A skill that drifts from the engine is worse than no skill.

**Evaluating rather than building? Open [`overview.html`](overview.html).** One self-contained page: what the layer does and why the gate in the middle is the differentiator, how it is built (skills, subagents, MCP, the hook, the knowledge tiers), the integration surface, how to drive it from Claude Code today — and then the engagement-scale `portfolio_bid` study as the demonstration, with the briefing it produced attached in full. Regenerate with `julia --project=. agent/scripts/render_session.jl agent/overview.json`.

**[`pitch.html`](pitch.html)** is the same capability for a reader who will not evaluate an architecture: one page, no CLI and no component inventory, the money case as the proof, and an explicit list of what the layer does *not* claim. Both pages quote the same study, so a number that changes in `studies/portfolio_bid/` changes in both — grep before editing either.

## Layout

| Path | What it is |
|---|---|
| `scripts/rdx.jl` | the CLI — one entry point per phase, JSON on stdout, meaningful exit codes |
| `scripts/diagnose.jl` | the plausibility checks (`rdx diagnose`), the part that carries modeling judgement |
| `scripts/gen_reference.jl` | regenerates `reference/model_schema.md` from the loaded module; `--check` is the drift gate |
| `scripts/acceptance.jl` | the known-answer test — drives the published shadow price through the JSON path |
| `scripts/render_session.jl` | renders a page document — authored sections, a captured transcript, or both — as one self-contained HTML file |
| `scripts/chartkit.jl` | data-driven SVG charts (bars with error bars, paired dots, waterfall, traced ceiling) — no plotting dep, so a study's figures build anywhere the CLI runs |
| `skills/` | the skills themselves (`rd-study` orchestrator + one per phase) |
| `reference/` | generated extracts of the authoring schema — never hand-edited |
| `studies/` | worked studies; `marginal_scientist/` is the known-answer acceptance test, `portfolio_bid/` the one at engagement scale |
| `sessions/` | captured end-to-end sessions — what the workflow is like to *use*, not just what it produces |
| `overview.json` → `overview.html` | the capability overview page: authored sections, then the `portfolio_bid` transcript imported as its demo |
| `pitch.json` → `pitch.html` | the same capability as a one-page pitch — no tooling surface, the case as the argument, the non-claims stated |
| `figures.jl` → `figures/` | the three structural diagrams on those pages (phase pipeline, the same sequence in plain language, integration surface), generated from tables that cite the charter's tracker ids |

## The skills

`rd-study` is the entry point; it runs the four phases and reads one phase skill at a time.

| Skill | Phase | Produces | Gate to clear |
|---|---|---|---|
| `rd-study` | orchestrator | the study directory | the three-part qualifying test — something binds, something takes time, something is decided mid-flight |
| `rd-frame` | 1. frame | `BRIEF.md` | every slot filled, failure modes pre-registered |
| `rd-author` | 2. author | `model.rdj.json` + one file per arm | `rdx validate` returns zero diagnostics |
| `rd-doctor` | 3. diagnose | `DOCTOR.md` | `rdx diagnose` exits 0, or every failure is understood and stated |
| `rd-report` | 4. report | `FINDINGS.md` | headline, decision rule and caveats all present |

The gates are the point. DIAGNOSE standing between AUTHOR and REPORT is what keeps the workflow's most likely bad output — a confident wrong number — from reaching a memo.

## The CLI

Everything runs through `rdx`, from the repository root:

```bash
julia --project=. agent/scripts/rdx.jl <command> [args] [--flag=value ...]
```

| Command | Does |
|---|---|
| `validate MODEL` | static, eval-free schema check → diagnostics with JSON paths |
| `run MODEL --seed=N` | one seeded run → per-place troughs, per-transition spawns/completions, ledger totals; `--trace=P[,P...]` adds the full per-tick series for those places, `--out=DIR` writes an export bundle |
| `diagnose MODEL --seed=N` | the plausibility gate — nine checks against the brief's declared intent |
| `ensemble MODEL --seeds=N --metric=M` | seeded ensemble → mean ± sem, quartiles; `--per-seed` adds the raw member values |
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

Both checks are wrapped as runnable gates, and both run in the `Agent layer` CI workflow on any change under `agent/` or `src/`:

```bash
julia --project=. agent/scripts/acceptance.jl              # the known-answer study, end to end
julia --project=. agent/scripts/gen_reference.jl --check    # the generated reference is not stale
```

`acceptance.jl` asserts at two strengths: the *scientific* claim is hard (the delta lands inside the published \$19.0M ± \$2.5M band and clears 2·se), while bit-reproducibility against the recorded run is a soft note — Julia's `hash` is not stable across versions and the member seeds derive from it, so the stream can move without any semantics changing.

## Worked example — a study at engagement scale

`studies/portfolio_bid/` is the same four artifacts on a decision the size of real work: a VP of corporate development pricing three in-licensed clinical assets against a 40-FTE development group that the company's own pipeline already fills. Eleven places, eight transitions, twenty-nine arcs, **nine arms**, 120 seeds per arm — and no structured tokens and no `host.jl`, so the entire study is data. The headline is a number an NPV model structurally cannot produce: the same three assets are worth **+\$744M ± 144** with the FTE ceiling lifted and **+\$241M ± 90** against the real group, so **\$503M ± 169 of the seller's valuation is stranded in a queue** — and twelve more development FTE are worth **+\$370M ± 106 given the deal** against **−\$8M ± 91 without it**.

It is also where the layer's numeric discipline is visible, in three files that are worth copying into any study of this size: `README.md` is the **provenance file** — every command, every number it produced, and the rule that the report may quote nothing absent from it; `figures.jl` draws through `chartkit.jl` straight out of `results/*.json`, so a figure cannot drift from the number it illustrates; and `numbers.jl` **re-derives every published table** from those same envelopes, simulating nothing. The last one exists because `validate` and `diagnose` both inspect the *model*, and neither can see a number mistyped into prose — which is a mistake this study made and caught that way.

## Captured sessions — what it is like to use

`studies/` shows you what the workflow *produces*. `sessions/` shows you what it is like to *use*: which questions get asked, which defaults get proposed, what happens when the user's numbers contradict the assistant's arithmetic. Each session is a JSON transcript rendered to one self-contained HTML page — no external CSS, fonts or JS, figures inlined.

```bash
julia --project=. agent/scripts/render_session.jl agent/sessions/field_service.json
```

A session's appendix publishes the study's own deliverable rather than a copy of it: an appendix entry takes `"file": "…/FINDINGS.md"` in place of `"text"`, imported at render time, so the page cannot quietly disagree with the report the study shipped.

Two kinds of transcript live here and each page states which it is. A **blind capture** is recorded by two agents that cannot see each other's reasoning — evidence rather than marketing, and so far the layer's most effective bug-finder. A **scripted walkthrough** composes the dialogue around a finished study, where every command shown was really run and every number is a real measurement traceable to the study's `README.md`. The protocol and the limits of each are in [`sessions/README.md`](sessions/README.md).

## Keeping the reference honest

`reference/model_schema.md` is generated, so the skills cannot describe an engine that no longer exists. One table inside the generator is hand-maintained — the JSON-key ⇒ store-column mapping, which lives function-local inside `src/serialize.jl` and cannot be read out of the module. It is guarded rather than trusted: the generator asserts every column it names exists in `SCHEMA` and every JSON key it documents is read by `src/serialize.jl` or `rdx.jl`, and refuses to emit otherwise. Regenerate with `gen_reference.jl` after any change to the schema, the whitelists or the action verbs.

## Portability

The CLI is plain Julia over Bash, and the skills are Markdown. Nothing here depends on Claude Code: any agent that can run a shell command and read a file can drive this workflow. Claude Code's skill/plugin packaging is the *distribution* mechanism, not the substrate.
