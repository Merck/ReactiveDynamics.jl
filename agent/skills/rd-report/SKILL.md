---
name: rd-report
description: Phase 4 of an rd-study — turn simulation output into a decision briefing. Runs the two-arm treatment effect or ensemble summary, checks the effect against the noise, and writes FINDINGS.md to the case-study signature: headline number with error bar, decision rule, why a spreadsheet cannot produce it, collapsible formalism, hero visual, replay provenance. Use after rd-doctor passes.
---

# rd-report — the number, and the decision it supports

DOCTOR passed, so you are now allowed to quote a number. Your output is `FINDINGS.md`, written to the same signature the published case studies use — `agent/studies/marginal_scientist/FINDINGS.md` is the exemplar, and it is worth reading before you write a line.

The audience is someone who has to act and was not in the room. They need the answer, the rule that uses it, and a reason to trust it. In that order.

## Get the numbers

Two arms — the usual case, a lever:

```bash
julia --project=. agent/scripts/rdx.jl effect <baseline> <arm> --host=<host.jl> \
  --metric=enpv --scale=0.001 --annual=0.08 --seeds=48 --root-seed=2026
```

One arm — a level, a bottleneck, a distribution:

```bash
julia --project=. agent/scripts/rdx.jl ensemble <model> --host=<host.jl> \
  --metric=final:launch --seeds=48 --root-seed=2026
```

Member *k* of each arm is seeded from `hash((root_seed, k))`, so both arms are compared on the **same seed stream** — the paired comparison is what makes a modest effect detectable at all. Record the root seed in the report; it is what turns a number into a replayable number.

**Is the effect real?** `effect` reports `significant_2se`. If it is false, the honest headline is *"we cannot distinguish this from noise at N members"* — either raise `--seeds` and rerun, or report that. Never report a point estimate whose sign the data does not support. Sample-size honesty is cheap here: members are independent, so if 48 doesn't separate it, try 200 before concluding.

## The signature

### 1. Headline — the number, in the title

`## **+$19.0M ± $2.5M** — hire them`

The number and the verdict, in units the reader thinks in, above everything else. Then one paragraph: the two levels being compared, the horizon, and how many standard errors from zero. If the reader stops here, they should already be able to act.

### 2. The decision rule

Not "the fifth scientist is worth \$19M" but **"hire iff the fully loaded cost over 48 months is below \$19M."** Turn the number into the comparison the user actually makes — that sentence is the deliverable and everything else is support. Where it is natural, extend it: a shadow price that falls as the constraint relaxes traces a diminishing-returns curve, so say where hiring should *stop*, not just that the next hire clears.

Then the caveats, in the same section, as scope rather than apology: what the number counts and what it doesn't, carried down from DOCTOR's warnings. A truncated horizon makes the headline a **lower bound** — state that as a direction, not a disclaimer.

### 3. Why a spreadsheet cannot produce this

One short section, and do not skip it — it is the section that justifies the method to a skeptical reader and is the strongest pitch the framework has. Be specific to *this* study: name the three things that had to be true at once (the resource genuinely ran out; work *held* it for a while; something had to decide who waits) and note that changing any one changes the answer. Then name the mechanism: the marginal unit doesn't produce the outcome, it lets work *start* that would otherwise have queued, some fraction of which pays off later, probabilistically, wherever the allocator sent it. No closed form gets there.

### 4. The model — collapsible

```markdown
<details>
<summary>Formalism — places, transitions, modalities (click to expand)</summary>
...
</details>
```

Collapsed because the decision-maker does not need it and the reviewer does. Inside: a table of the elements in operator language with their formal role, the **units** stated explicitly, and the arm diff as a fenced `diff` block — that diff *is* the lever, so show it rather than describing it. Note that the model is data and name the only Julia in the study (`host.jl`, if any).

### 5. Hero visual — collapsible

```bash
julia --project=. agent/scripts/rdx.jl draw <model> --host=<host.jl> --exec-map --seed=1 --out=figures/exec_map.svg
```

The execution map decorates the structure with run statistics, so the starved pool is visible as the thing the system is short of. One figure that carries the argument beats four that decorate it. Say in the caption what the reader should see — and note what *didn't* starve, since that is how the reader knows which constraint is the real one. If Graphviz is unavailable `rdx` writes `.dot` instead and warns; report a `.dot` rather than fabricating a figure.

### 6. How we got it, and provenance

Number the four phases with links to `BRIEF.md`, the model documents, `DOCTOR.md`. Then the **exact command** that produced the headline, verbatim and runnable from the repository root, and the statement that `(model, seed)` determines the result — rerun it and you get \$18.991M again, not a number near it.

## Rules that hold without exception

- **Every number traces to a command you ran.** No hand-computed figures, no numbers from an earlier draft of the model, no rounding a `18.991` into a `20` because it reads better. Round for presentation in the prose, keep the precise value in provenance.
- **Report the effect against its noise, always.** A delta without a standard error is not a finding.
- **Carry DOCTOR's warnings into the caveats.** If a warning doesn't appear in FINDINGS, you suppressed it.
- **A null result is a result.** "The constraint doesn't bind; the next hire is worth their cost" is a finding that saves the user money. Deliver it with the same confidence as a positive one.
- **Units in the report.** Pools in \$k, results in \$M, ticks in months — state the conversion. This is the error that survives every check and lands in the memo off by a factor of 1000.

## Offer next, don't do it unasked

The infrastructure is now built, so the second question is cheap — a third arm is one file and one command. Offer the obvious follow-ups: the next step on the curve (5→6, 6→7) to find where returns diminish; a different lever (the priority rule instead of the headcount); a sensitivity sweep on the input the user was least sure of in the brief. Name them, let the user pick.
