# `agent/sessions/` — captured rd-study sessions

A study's four artifacts (`BRIEF.md`, `model.rdj.json`, `DOCTOR.md`, `FINDINGS.md`) show you what the workflow *produced*. They do not show you what it is actually like to use, which is the thing a prospective user needs to see: which questions get asked, which defaults get proposed, what happens when the user's numbers contradict the assistant's arithmetic, and how a business question turns into a document with a validator behind it. This directory holds transcripts of whole sessions, rendered as single-file HTML so they can be mailed, committed, or opened in five years without a build step.

```bash
julia --project=. agent/scripts/render_session.jl agent/sessions/<name>.json
# → agent/sessions/<name>.html — self-contained: no external CSS, fonts, or JS; SVGs inlined
```

The transcript JSON is the source of truth and is meant to be read as diffable text; the HTML is generated and safe to regenerate. The schema is documented in the header of `agent/scripts/render_session.jl`: a `title`/`subtitle`/`headline` block, a `meta` map rendered as the provenance strip, and `turns[]`, each with a `role` (`user` / `assistant` / `note`), the turn text, an optional `phase` label that starts a new band, an optional `tools[]` of `{cmd, out}` pairs rendered as a collapsed evidence block, and an optional `figure` path that gets inlined.

## How these are captured, and why you can trust them

A transcript is only worth showing if it was not written to flatter the thing it demonstrates. The sessions here are **role-played by two separate agents that cannot see each other's reasoning**, under a protocol designed so that the interesting behaviour cannot be faked:

1. The **operator** agent is given a business situation and, *before speaking a single word*, writes its facts to a file — call volumes, durations, resource counts, contract terms, penalty figures, today's observed backlog. That file is a pre-commitment: the operator must answer from it, and cannot invent a convenient number later to rescue a conclusion.
2. The **assistant** agent starts from an empty context, is pointed at `agent/skills/rd-study/SKILL.md`, and is forbidden from reading the directory the fact sheet lives in. Everything it knows about the business, it elicited.
3. A relay passes messages verbatim in both directions. Neither agent sees the other's tool calls, scratch work, or instructions.

The consequence is the property that makes the transcript evidence rather than marketing: **every number in the model was extracted by asking**, and every number in the result was produced by running the engine. Where the assistant proposed a placeholder and the operator corrected it, you can watch the correction land.

The fact sheet is committed next to the transcript, with a line-by-line audit of what the interview recovered and what it missed: [`field_service_facts.md`](field_service_facts.md). **Read the audit before the transcript if you are evaluating rather than browsing** — it is where the workflow's weaknesses are, and it names three: the interview took a *range* for a quantity the operator had a point estimate for (which is how the tick-grid defect survived to phase 4), it never modelled the value of a completed repair, and it carried none of the three quantities the operator had explicitly flagged as guesses into the study as sensitivities. The pattern is that the interview was strong on contractual facts and weak on exactly the numbers its informant had labelled soft.

What this does *not* prove: that a real human would answer as fluently as a role-played one, or that a real operation matches the shape the persona described. A captured session demonstrates the workflow, not the state of anyone's field-service department. The header of each rendered page says so. In particular this persona was *cooperative* — it answered in good faith and volunteered its own uncertainty markers — so the capture says little about the harder case of an informant who misremembers or contradicts themselves between sessions.

## Sessions

| Session | Question | Outcome |
|---|---|---|
| [`field_service`](field_service.json) | A lab-equipment vendor's field-service manager has a budget meeting in six weeks and three candidate levers: hire 2 more technicians, buy more diagnostic equipment, or change the dispatch priority rule. | The lever the manager arrived asking about ranks worst of the three: 2 technicians avoid \$57,504 ± 21,741 of annual penalties for \$190k/yr, while 2 diagnostic laptops avoid \$228,050 ± 18,712 for \$38k once. The free policy change cannot be separated from zero on any of five metrics at 48 seeds. **Then the operator pushed back** — their own 4.3 h rule of thumb turned out to be unrepresentable on the model's 1-hour tick grid — and the assistant had to find, report and repair a defect in its own model on a finer grid, which moved the laptop estimate *up* to \$288,644 ± 31,502 and retired three published figures. The repair is the most useful part of the transcript. |
