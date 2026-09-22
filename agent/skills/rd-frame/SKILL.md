---
name: rd-frame
description: Phase 1 of an rd-study — turn a vague business question into a decision brief through Socratic interview. Elicits the decision, the lever, the scarce resource, the process that consumes it, the contention rule, the horizon, and the headline metric; pre-registers what could go wrong. Produces BRIEF.md. Use before authoring any ReactiveDynamics model.
---

# rd-frame — from a question to a brief

Your output is `BRIEF.md`: a document in the *operator's* vocabulary where every slot nonetheless maps 1:1 onto the object model, so AUTHOR can compile it without a second interview. Read `agent/studies/marginal_scientist/BRIEF.md` first — it is what "already a model in all but syntax" looks like.

**Interview, don't interrogate.** Ask two or three questions at a time, propose a default for each, and let the user correct you. "I'll assume a trial holds one scientist for four months — is that right?" gets an answer; "how long does a trial take?" gets a meeting. Users know their system and do not know which details matter; deciding what matters is your job.

## Start with the decision, not the system

The most common failure of this phase is modeling a system nobody has to decide anything about. So open with:

1. **What decision hangs on this?** If nothing does, stop and ask what the user actually wants — a number nobody acts on is not worth a study.
2. **What is the lever?** Something the user can actually change: a headcount, a budget, a priority rule, a deal, a policy. A study with no lever is a description, not a decision.
3. **What would make you act differently?** This is the answer's *shape* — a threshold, a ranking, a go/no-go. Write it down before computing anything; it is what keeps you from producing a number the user cannot use.

Then run the qualifying test from `rd-study` explicitly and record the verdict in the brief. If it fails, say so and offer the arithmetic.

## The slots

The brief's centerpiece is a table with one row per slot. Fill every row; an unfilled row becomes an invented default in phase 2, which is how models quietly answer the wrong question.

| Slot | What to ask | Where it lands |
|---|---|---|
| **Scarce resource** | "What runs out? What do people queue for?" | a plain place, `init` = how many |
| **Reusable or spent?** | "When work finishes, do you get it back?" — people/equipment yes, cash/materials no | `conserved` vs the default `consumed` arc modality |
| **Metered resource** | "Anything drawn continuously while work runs?" (budget, power) | a `perstep` arc out of a deep pool |
| **Work in progress** | "What flows through — and does each item have attributes that matter?" | a plain place if items are interchangeable; a structured place if each carries state (program, phase, value) |
| **Starting state** | "What is already in flight on day one?" | `init` for plain places, `population[]` for tokens |
| **Intake** | "Where does new work come from, and how fast?" | a transition, `rate` per tick |
| **The process** | "What does one unit of work consume, for how long, and how often does it succeed?" | a transition: lhs arcs, `cycletime`, `prob_of_success` |
| **Contention rule** | "When two things want the resource at once, who wins?" | `priority` weights |
| **Outcome that carries value** | "What produces the payoff, and when is it realized?" | a place with `reward`, produced on the rhs |
| **Horizon and tick** | "Over what period is this decided? What is the natural unit?" | `meta.tspan`, `meta.dt` |
| **The lever** | "What exactly changes between the arms?" | the diff between arm documents |
| **Headline metric** | "What number do you want?" | a named `rdx` metric |

### The questions that earn their keep

- **"When the work finishes, do you get the resource back?"** The default arc modality is *consumed*. Headcount declared by default drains to zero and never returns — the single most common silent error in this formalism, and it makes a model that understates capacity for the whole run and overstates the value of a hire. Ask it about every resource, explicitly.
- **"When both want it at once, who waits?"** If the user says "it works itself out", there is a real rule and nobody has written it down. Getting it on paper is often worth more to the user than the simulation.
- **"What is the unit?"** \$k vs \$M, months vs weeks. Unit confusion is the error that survives every check and lands in the memo off by 1000. Pin it here and restate it in the report.
- **"How would you know this model is wrong?"** The answers become the pre-registered failure modes below.

## Pre-register what could go wrong

This is what separates a brief from a summary, and it is where the workflow gets its credibility. Before any model exists, write down the ways the answer could be an artifact — then hand DIAGNOSE the flags that test each one:

| Failure mode | Pre-registered as |
|---|---|
| the resource doesn't actually bind → there is no study | `--expect-binds=<place>` |
| a reusable resource declared as consumed | `--expect-reusable=<place>` |
| the *wrong* thing binds → the headline is about the wrong lever | `--expect-binds=<place>` |
| a legitimate cumulative counter looks like runaway growth | `--tally=<place>` |
| horizon too short to realize the outcome | expect the `horizon_sufficient` warn; decide now whether it is scope or defect |
| the effect is inside the ensemble noise | the `effect_vs_noise` check on `rdx effect` |

Pre-registration is the point: a check the brief demanded *before* the run existed is evidence, while a check invented afterwards to explain a number is decoration.

## Finish

Write `BRIEF.md` with: the decision · why not a spreadsheet · the qualifying test verdict · the slots table · what could go wrong · what "done" looks like. Show the user the slots table and ask what is wrong with it. Then hand to `rd-author`.

Do not write JSON in this phase. If you catch yourself reaching for it, you have stopped listening.
