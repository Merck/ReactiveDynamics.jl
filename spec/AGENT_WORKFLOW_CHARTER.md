# Agentic Modeling Workflow Charter — ReactiveDynamics.jl

> **The tracking spine for the agentic-workflow effort**: a portable skill/plugin layer that takes a user from a business question to a calibrated ReactiveDynamics model, a reproducible run, and a decision briefing. The counterpart to [DOCS_CHARTER.md](DOCS_CHARTER.md) (reader-facing surface) and [STATUS.md](STATUS.md) (engine state) for the *agent-facing* surface. Nothing here is normative about engine behavior — [CONTRACT_DRAFT.md](CONTRACT_DRAFT.md) and [adr/](adr/) remain the source of truth, and this layer **cites** them rather than restating them (the §D4 anti-fork rule of the docs charter, extended to a fourth consumer). Drafted 2026-09-22. Re-verify any `file:line` before acting.

## 1. Why this exists

The engine is built, green and documented (835 pass, §1–§15 implemented, a full Documenter site with three tiered tutorials and three decision case studies). What gates adoption is not a missing feature — it is **the distance between "I have a business question" and "I have a model whose scarce resource genuinely binds."** Closing that distance today requires tacit expert knowledge that is spread across the contract, eighteen ADRs, and the calibration discipline recorded in [demo/bd_acquisition/MVP_BD_DEMO.md](../demo/bd_acquisition/MVP_BD_DEMO.md): that the default modality is *consumed* and headcount must be declared `@conserved`; that `@rate` with a concrete `cycletime == 0` is rejected at construction; that macro arguments are literal so a swept parameter needs a builder function; that a reward on a product place realizes at completion and not at intake; that a marginal number is meaningless unless something is actually contended.

That body of knowledge is exactly what an agent carrying the spec can hold. **The agentic workflow is the framework's distribution channel: it turns a framework that requires an expert into one that requires a question.**

**Done** for this effort means: a user in Claude Code (or any skill-capable harness) says what decision they face, is interviewed into a written decision frame, watches a ReactiveDynamics model be authored as data and statically validated, gets a plausibility verdict *before* any headline number is quoted, and ends with a reproducible `(model, seed)` briefing in the house signature — question, number, verdict, then how we got it. Every artifact lands on disk, reviewable; nothing important lives only in chat scrollback.

## 2. The value thesis — where the leverage actually is

Three claims, in descending order of how much they justify the effort.

**2.1 RD is the discrete-event engine an agent can safely author.** This is the differentiated claim and it is a property of the design, not a hope. A model IS data: one eval-free JSON document ([reference/json_schema.md](../docs/src/reference/json_schema.md), ADR 0005) whose every expression is a node-tagged tree from a closed whitelist, whose every host construct is referenced *by name* through a per-network registry, and which is never `Meta.parse`d or `eval`d on load. `ReactiveDynamics.validate(dict)` is a **pure static pass returning `Vector{Diagnostic}`** — severity, JSON path, message — that `from_json_model` gates construction on. The source comment at `src/serialize.jl:652-654` says it outright: *"An LLM runs this to self-check a model before from_json."* The engine was built for this.

Compare the alternative. A SimPy or AnyLogic model is code; an LLM authoring one is doing code generation, with an arbitrary-execution surface and no critic short of running it. Here the agent emits structured data, gets machine-readable diagnostics with paths it can repair against, and — because `(model, seed)` determines a run contractually (CONTRACT §4) — any number it quotes replays exactly. **Closed IR + pure validator + contractual determinism is the combination that makes agentic authoring trustworthy, and no competing DES tool has it.** This belongs in the D3a methodology paper's statement of need, not just here.

**2.2 The expert knowledge that makes a model *bind* is spec-encoded and transferable.** The framework's central insight — under contention the value of a resource is not its accounting cost but what the system does with it at the margin — only pays off in a model where something is actually scarce. The three case studies each open by explaining the modeling choices that make that true ("three modeling choices make the scientist bench the true bottleneck"). An agent with the modality truth table (CONTRACT §1.1), the allocator semantics (ADR 0002) and the calibration discipline can make those choices and, more importantly, **verify them** — which is §3.3 below and the highest-value new artifact in this charter.

**2.3 The output is already shaped like a decision memo.** The docs charter's §2 discipline (question-titled, headline number first, verdict, why-a-spreadsheet-cannot, then "how we got it") is a signature the workflow can emit directly, and it composes with the user's existing `decision-brief` and `business-writing` skills. A simulation that ends in a memo is adopted; one that ends in a trajectory dump is not.

**2.4 The honest limit — when to say no.** If nothing is contended, nothing has duration, and no decision is taken *inside* the horizon, a spreadsheet is the right tool and the workflow must say so and stop. A skill that qualifies out builds more trust than one that always produces a model, and it protects the framework's reputation from a corpus of agent-built models that answer nothing. The qualifying test is three questions: **does something bind, does something take time, is something decided mid-flight?** One "yes" is thin; two is a real study; three is the framework's home ground.

## 3. The workflow — four phases, each with an artifact and a gate

The phases are sequential, resumable, and each ends in a file. A phase may not be skipped, but a returning user re-enters at whichever phase the study directory's contents imply.

### 3.1 Phase FRAME — decision framing (formalism-independent)

**Register: the operator's, not the engine's.** The interview never says *place*, *arc* or *modality*; it says resource pool, step, how long, what's scarce. The glossary's "two registers, deliberately" section already licenses this ([docs/src/glossary.md](../docs/src/glossary.md) §"Two registers") — it is the same object, named for the room you are in.

**Start at the end.** The first questions are the decision, the number that would settle it, and what would change the user's mind. A model built before the question is known answers nothing; the case studies are titled with the question for the same reason. This also yields the study's acceptance test: *the briefing must contain the number named here, with a spread.*

**Interview formalism-independently, record in slots that compile.** This is the phase's central design constraint. A genuinely unstructured interview produces a beautiful narrative nobody can turn into a net. So the *questions* stay in business language while the *record* is a slot table that maps 1:1 onto the object model: decision + deciding number + horizon and tick; entities that flow and whether they have identity (fungible count vs structured token with attributes and history); resources, and for each one whether it is **reusable** (held then returned → `@conserved`), **consumed** (spent → default), or **metered per unit time** (burn rate → `@rate`), plus whether it merely gates without being held (`@nonblock`); steps, with duration, success probability and concurrency limit; what is uncertain and how; which decisions are taken *inside* the run and on what trigger (the endogenous channel, ADR 0010); what money or value accrues where; and the contention story — who competes for what, and who should win (priority).

**Ask where it could go wrong, explicitly.** A dedicated slot for *the modeling risks the user already suspects* — what is being assumed away, which number is a guess, where the data is weak. This is carried through to the briefing's limitations section rather than discovered late, and it is where the Socratic register earns its keep: the productive questions are the ones that surface a contradiction between two of the user's own answers (an intake rate that outruns any plausible capacity, a horizon shorter than one cycle of the process, a "scarce" resource nothing competes for).

**Output:** `BRIEF.md` — the decision frame, human-reviewable, with the slot table and the risks. **Gate:** the qualifying test (§2.4) answered, the deciding number named, and the frame read back to the user in their own words for confirmation. This is the one phase that should *not* be automated past the user.

### 3.2 Phase AUTHOR — frame → ReactiveDynamics model

**The artifact is `model.rdj.json`, not Julia.** Authoring as data is what makes this phase safe and checkable. The agent emits the document, runs `RD.validate` (unexported — call it qualified), and repairs against the returned paths until the diagnostic vector is empty. That loop is fast, deterministic, eval-free and needs no simulation — it is a compiler in the loop, and it should run to a clean pass *before* anything is executed.

**The one code seam, named honestly.** A structured-token *kind* with new attribute fields needs a host Julia struct and a registry entry (ADR 0006 §B/§C; the pattern is `demo/bd_acquisition/host.jl`). Pure-JSON authoring cannot introduce one. Two ways out, and this is [Open decision 2](#11-open-decisions): ship a *generic* token kind with a free-form attribute bag in the plugin's support code (or upstream in RD) so the common case stays data-only; or have the agent emit a small `host.jl` that is presented for review as code. The distinction matters to the trust story and must not be blurred — JSON is inert, `host.jl` is executed.

**Author the model, then explain it back.** Before running anything, the agent renders the net (`draw_network`, Layer A — structure only, no run needed) and walks the user through it against the frame: this pool is your bench, it is `@conserved` because you said people come back, these two steps compete and oncology wins because you gave it priority. A user who cannot follow the picture will not trust the number. This is the review surface the user asked about, and today it is a diagram plus a slot-to-arc table — see §7 for where it goes next.

**Output:** `model.rdj.json` (+ `host.jl` only if unavoidable), `study.jl` (the runnable driver), and a `MODEL.md` narrating the frame→model mapping. **Gate:** `validate` clean, the net drawn, every `BRIEF.md` slot traceable to a row in the document, and no slot silently dropped.

### 3.3 Phase DIAGNOSE — the plausibility gate (the highest-value new artifact)

**`validate` catches schema errors; it cannot catch a plausible-but-wrong model.** That gap is where LLM-authored simulations fail, and nothing in the engine closes it today. This phase is a new checklist that runs a short seeded smoke simulation and reports, as structured data, whether the model is *dynamically* worth believing. Every check below is computable from `prob.sol`, `prob.log`, the per-program ledger and the pool troughs that `exec_map` already derives (`_pool_troughs`, `src/visualize.jl:240` — `minimum(prob.sol[!, col])`, so no new engine code is needed for the core signals).

- **Did anything fire?** A transition with zero firings over the horizon is dead — rate too low, LHS never satisfiable, or gated off. A dead transition silently removes a whole branch of the model.
- **Does anything bind?** No pool ever reaching its trough means no contention, which means the marginal question is trivial and the framework added nothing over arithmetic. **This is the gate that decides whether the study is worth reporting.**
- **Does the *intended* thing bind?** The starved pool must be the one `BRIEF.md` names as scarce. If budget binds when the frame says people are scarce, the headline number answers a different question than the one asked.
- **Modality-vs-intent mismatch.** A resource the frame calls reusable but the document declares consumed (the default, hence the most common silent error) — a direct slot-to-arc comparison, and the single highest-yield check in the list.
- **Is the horizon long enough?** In-flight firings at `tspan`, outputs still rising, or zero completions of the terminal step mean the number is horizon-truncated rather than converged.
- **Unbounded growth.** A pool diverging monotonically indicates a missing sink or a mis-scaled intake.
- **Ledger sanity.** Value realizing at intake rather than at completion; costs identically zero; rewards attached to a stock rather than a product.
- **Degenerate stochasticity.** Success probabilities pinned at 1.0, or a single-seed conclusion that flips across seeds — the model is not actually risk-adjusting.
- **Is the effect bigger than the noise?** Once an ensemble runs, `treatment_effect` returns a standard error; an effect inside ±2 se is not a finding. Enforcing this is what keeps the workflow from manufacturing confident nonsense.

**Output:** `DOCTOR.md` — per-check verdict with the evidence, plus the exec map as the visual. **Gate:** no failing check, or each failure explicitly acknowledged and carried into the briefing's limitations. A red "nothing binds" should *block* the headline number, not footnote it.

This phase is also where calibration iteration lives: making something bind at a sensible scale over a long-enough horizon is a search, and it is the natural place for an autonomous loop (§10).

### 3.4 Phase REPORT — the briefing

Run the real thing: a seeded ensemble, the comparison the frame calls for (`ensemble` / `summarize` / `treatment_effect`), the export bundle (`export_run` / `export_ensemble`, which writes a `run.json` manifest carrying model hash, seed, tspan and the embedded model — so the bundle is self-describing and replayable).

Then emit the **house signature**, which is the docs charter's case-study shape and needs no invention: the question as the title; the number and the verdict in the first paragraph; the decision rule that falls out of it; why a spreadsheet cannot produce it; then *how we got it* — the model (collapsible), the mechanics, the figures; then limitations (from `BRIEF.md`'s risk slots and `DOCTOR.md`'s findings); then reproduction (`(model, seed)`, RD version, one command). Prose register delegates to the `business-writing` / `decision-brief` skills; charts delegate to `dataviz`; the network visuals come from `exec_map`.

**Output:** `FINDINGS.md`, the run bundles, figures, and optionally `card.html`. **Gate:** every number in the briefing traces to a bundle artifact, carries a spread, and replays from the recorded seed.

## 4. The artifact contract — a study directory

Not invented: this is [demo/bd_acquisition](../demo/bd_acquisition) generalized. That demo already separates `model.rdj.json` (data) from `host.jl` (the ADR-0006 boundary, in its own words) from `run_demo.jl`/`analysis.jl`/`figures.jl` (drivers) from `MVP_BD_DEMO.md` (the design record), and even carries an `RD_BD_FROM_JSON=1` toggle proving the JSON and DSL paths agree. The workflow's job is to produce that shape from an interview.

```
studies/<slug>/
  BRIEF.md            # FRAME    — the decision frame + slot table + risks
  model.rdj.json      # AUTHOR   — the model AS DATA (eval-free, validate-clean)
  host.jl             # AUTHOR   — ONLY if a custom token kind is unavoidable (reviewed as code)
  study.jl            # AUTHOR   — the driver: build, run, ensemble, export
  MODEL.md            # AUTHOR   — the frame→model mapping, narrated
  DOCTOR.md           # DIAGNOSE — per-check verdict + evidence + exec map
  runs/<tag>/         # REPORT   — export_run / export_ensemble bundles (run.json manifest)
  figures/            # REPORT   — SVG/PNG
  FINDINGS.md         # REPORT   — the signature briefing
  card.html           # optional — self-contained model + results card
  study.lock.json     # provenance: RD version/commit, seeds, model hash, gates passed
```

Phase is inferred from which files exist plus `study.lock.json`, so the workflow is resumable across sessions with no chat state. `study.lock.json` is what makes a quoted number auditable months later.

## 5. The capability layer — what is a skill, what is MCP, and why

The division follows the frontier pattern: **skills carry knowledge and procedure; MCP carries capability the model does not have; a plugin is the distribution unit; subagents buy context isolation and parallelism; hooks buy deterministic enforcement.** Applied here:

### 5.1 Plugin `reactive-dynamics` — the distribution unit

Bundles the skills, agents, commands, `.mcp.json` and hooks below. Slash commands as the entry points: `/rd:study` (start or resume), `/rd:frame`, `/rd:author`, `/rd:doctor`, `/rd:report`, `/rd:explain` (look something up in the spec).

### 5.2 Skills — five, with progressive disclosure

Deliberately not ten. The split is by *mode*, because the modes need disjoint knowledge and disjoint registers, and you do not want the JSON schema loaded during the interview.

| Skill | Role | Knowledge it loads |
|---|---|---|
| `rd-study` | Orchestrator/router. Owns the study directory, phase inference, and the gates. Thin. | The workflow, the artifact contract, the gate table. |
| `rd-frame` | The interview → `BRIEF.md`. Operator register; never names the formalism. | The slot table, the qualifying test, the question bank, the contradiction patterns. |
| `rd-author` | `BRIEF.md` → `model.rdj.json`. Owns the validate/repair loop and the registry seam. | The JSON schema, the DSL surface, the modality truth table, the attribute contract, the authoring foot-guns. |
| `rd-doctor` | The §3.3 plausibility gate → `DOCTOR.md`. | The check list, how each signal is computed, what each failure implies. |
| `rd-report` | Results → `FINDINGS.md`. Delegates register and charts. | The signature template, the analysis surface, the spread discipline. |

### 5.3 Subagents — three, each for a real reason

- **`rd-modeler`** — runs the author/validate/repair loop in isolation. Justified by context hygiene: the loop generates noisy tool output (diagnostic vectors, partial documents) that the main thread does not need; it returns the finished document and a summary.
- **`rd-critic`** — adversarially reviews `model.rdj.json` against `BRIEF.md` in **fresh context**. Justified on the same grounds as `code-reviewer`: a model's author is the worst reviewer of its own assumptions, and the failure mode here (a plausible model answering the wrong question) is precisely what an independent reader catches. This is the analogue of the existing `pr-review-toolkit` pattern, pointed at models instead of diffs.
- **`rd-calibrator`** — parallel parameter sweeps to find the regime where the intended resource binds. Expensive and embarrassingly parallel; isolate it. *Phase 2 — only once `rd-doctor` can score a candidate.*

### 5.4 MCP server `rd` — and an honest account of why it is needed

In Claude Code, running Julia needs no MCP: `Bash` plus a checked-in script is the idiomatic default, and the skill ships those scripts regardless (they are also what makes the layer portable to harnesses with no MCP). MCP earns its place here for two specific reasons, neither of which is "capability":

1. **A persistent Julia session.** `using ReactiveDynamics` plus Plots is tens of seconds of load and JIT per cold process. An interactive modeling loop that pays that on every validate/run/redraw is unusable. A long-lived worker turns a ~45 s round trip into sub-second, and *session* is exactly the thing a stateless `Bash` call cannot hold — it can also retain the constructed `ReactionNetworkProblem` between calls, so redrawing the exec map under a different `@select` highlight is a styling pass over a finished run rather than a re-simulation (which is precisely what §15.2 says the overlay is). A precompiled sysimage is the complementary optimization, not a substitute: it cuts load time but still cannot hold a live `prob`.
2. **Cross-client portability.** Claude Desktop and other MCP hosts have no Bash. The same tool surface reaches them.

**What it connects to — nothing external.** Worth stating plainly, because "MCP" usually connotes a third-party service: this is a **local stdio server acting as thin transport over the local Julia runtime**, calling RD's own functions — `validate`, `from_json_model`/`to_json_model`, `ReactionNetworkProblem`/`simulate`, `ensemble`/`summarize`/`treatment_effect`, `export_run`/`export_ensemble`, `draw_network`/`exec_map`, and `prob.sol`/`prob.log`/the ledger for the §3.3 checks. No network egress, no credentials, no remote dependency.

**Tool surface — keep it small** (each definition costs context): `rd_validate`, `rd_simulate`, `rd_ensemble`, `rd_treatment_effect`, `rd_diagnose`, `rd_draw`. Plus the spec as MCP **resources** (contract sections, ADRs, glossary, JSON schema) so any host can attach the normative text — the portable answer to "leverage the ADRs/docs" for hosts that cannot simply read the repo.

**Architectural rule 1: `rd_validate` must call `ReactiveDynamics.validate`, never reimplement it.** A TypeScript or Python re-implementation of the schema check would fork the trust boundary and drift from `src/serialize.jl`, which is authoritative by the json-schema reference's own admission. The server is a thin transport over the engine, and thin is the whole design.

**Architectural rule 2: `rd_draw` degrades, never fails, without Graphviz.** The renderer needs a Graphviz **system binary**, not a Julia package — the same dependency that makes `test/semantic/visualization.jl:132` the one `@test_skip` CI actually takes (835 local vs 834 on CI). A modeling session on a machine without it must still frame, author, diagnose and report; it simply loses the diagram. Any tool that hard-errors on a missing system binary would make the whole workflow environment-fragile for a cosmetic reason.

### 5.5 Hooks — one, doing what only a hook can

A `PostToolUse` hook on writes to `*.rdj.json` that runs `validate` and surfaces the diagnostics. Hooks are for deterministic enforcement rather than persuasion: no amount of skill prose reliably prevents an invalid model being left on disk, and a hook simply makes it impossible to do so silently.

## 6. The knowledge layer — three tiers, and the anti-drift rule

The spec and the docs are not background reading for this layer; they are its working material. But a skill that *copies* the contract becomes a fourth divergent copy of the semantics — exactly what the docs charter's §D4 anti-fork rule exists to prevent (and why its D1/D2 web explanation quadrant was dropped). So the knowledge tiers are explicit, and only one of them is new prose.

**Tier 1 — pointers (never copied).** The skills cite normative artifacts by section: CONTRACT §1.1 for the modality truth table, §4 for determinism, §5 for the attribute contract, §8 for serialization, §9.5 for token filtration, §12 for rules, §14/§15 for analysis and visualization; ADR 0002 for the allocator, 0005/0006 for the IR and the registry boundary, 0008 for predicates, 0010 for the endogenous channel, 0017/0018 for vocabulary. Exposed as MCP resources and read on demand — progressive disclosure over the real documents, not a summary of them. *This is what "actively leverage the ADRs/docs" means operationally: the agent reads the normative text at the moment of need.*

**Tier 2 — generated extracts (cannot drift, by construction).** A script regenerates, from source: the public export surface, the JSON document's key set, the modality truth table, the op/dist/ref whitelists (`OP_WHITELIST`, `DIST_WHITELIST`, `REF_KINDS`, `PRED_OP_WHITELIST`), the attribute defaults. The suite already carries an `exports_resolve.jl` invariant asserting the docs' symbols resolve against the live surface and a retired-vocabulary grep gate; the same mechanism extends to the skill's extract, so **the skill's API reference is validated by the test suite that validates the engine's.** A drifted extract fails CI. This is the "no invented API" gate of the docs charter, inherited.

*Delivered as `agent/scripts/gen_reference.jl` → `agent/reference/model_schema.md`, gated by `--check` in the `Agent layer` workflow.* One honest limitation surfaced in the build: the JSON-key ⇒ store-column mapping is **function-local** inside `src/serialize.jl`'s read/emit pairs, not an exported table, so it cannot be read out of the module and that one table is hand-maintained. It is guarded rather than trusted — the generator asserts every column it names exists in `SCHEMA` (or is a real `ArcSpec` field) and every JSON key it documents is actually read by `src/serialize.jl` or `rdx.jl`, and *refuses to emit* otherwise. A rename on either side breaks the generator instead of silently producing a reference that lies. That assertion immediately earned its cost: it is what surfaced W14, by proving `population[].count`/`attributes` are read by nothing in `src/`.

**Tier 3 — genuinely new procedural prose (the only hand-written part).** The interview and its question bank; the qualifying test; the slot table and its mapping to the object model; the plausibility check list; the calibration discipline (currently tacit in `MVP_BD_DEMO.md` and the case-study openings); the briefing template. None of this exists anywhere yet, and all of it is procedure rather than semantics — so it cannot contradict the contract.

## 7. The review surface, and the road to an editor pane

The user's target is a model editor beside the chat, Claude-Science-like. The honest sequencing:

**Today: the study directory *is* the editor surface.** The model is a file; the user's editor shows it; Claude Code renders every edit as a diff. `draw_network` gives structure before any run, `exec_map` gives the annotated result after.

**Next increment: `card.html`, a self-contained model card.** Because the model round-trips to JSON, a static HTML+JS page can render the net, the attribute tables, the slot→arc mapping and the results panel **with no Julia in the loop** — a genuinely nice property of model-as-data. The recipe is already in the repo twice over: `demo/wires_viz_tour/build.jl` (executed Literate → inlined-SVG, self-contained HTML) and the finished brand board at `spec/design_system.html`, plus the headless-Chrome print path for a PDF. Editing in the card writes back JSON, which the hook validates. That is a live editor pane in the only sense that matters, reusing three assets that already exist.

**Later: bidirectional editing** (card mutates the document, chat and card share one file as state) is a real build, and worth doing only once the workflow's phases have proved themselves.

## 8. Reuse — existing skills and plugins this composes with

| Existing | Used for |
|---|---|
| `decision-brief`, `business-writing` | The briefing's register and structure. `rd-report` should *delegate*, not re-derive prose craft. |
| `dataviz` | Every chart. Its palette/validator discipline replaces ad-hoc plot styling and keeps figures consistent with the brand board. |
| `skill-creator` | Scaffolding the five skills and the plugin manifest correctly the first time. |
| `pr-review-toolkit` / `code-review` | The pattern `rd-critic` copies; also reviews the plugin's own code. |
| `ralph-loop` | The autonomous calibrate→diagnose iteration (§10). |
| `claude-code-guide` | Authoritative answers on hook/MCP/plugin mechanics during the build. |

## 9. Portability

Three layers, decreasing portability, deliberately ordered so the valuable part is the portable part:

1. **`skills/` + `scripts/`** — markdown plus Julia CLI scripts. Portable to any harness that reads skills, including the Agent SDK; usable by hand with no agent at all.
2. **`.mcp.json` + the `rd` server** — any MCP host.
3. **`plugin.json`, commands, hooks** — Claude Code specific.

The knowledge and the procedure live in layer 1. A non-Claude harness loses the slash commands and the validate hook, not the workflow.

## 10. Sequencing

1. **Spine first, thin.** All four phases minimally working end to end on one worked example, over `Bash` + scripts, no MCP, no subagents. The value is in the *sequence*; a vertical slice of the interview alone produces a `BRIEF.md` nobody can compile. The worked example should be a study whose answer is already known — re-derive the marginal-scientist shadow price or the in-licensing Δ-rNPV from an interview — because that is the only way to tell whether the workflow found the right answer or merely a confident one.
2. **`rd-doctor` for real.** The plausibility checks are the differentiator and the hardest part; build them against the three case-study models (known-good) and against deliberately broken variants (known-bad).
3. **MCP server** with the persistent session, once the script surface has stabilized enough to wrap.
4. **Subagents** `rd-modeler` and `rd-critic`.
5. **`card.html`.**
6. **Package as a plugin**, publish to a marketplace, and dogfood on a study nobody knows the answer to.

**Where `ralph-loop` fits.** Calibration is a search — adjust intake, capacity, horizon and priority until the intended resource binds, the effect clears the noise floor, and the doctor is green. That is a textbook autonomous loop: a concrete goal state, a machine-checkable score, bounded iterations. The rule is that the loop must be **gated on the numeric acceptance criteria from `DOCTOR.md`**, never on a subjective "looks good", or it will iterate forever. Use it for step 2's calibration and for sweeping the check list against broken models; do not use it for the FRAME phase, which needs the user in the room.

## 11. Decisions taken, and what is still open

**Resolved 2026-09-22 (maintainer).**

1. **The plugin lives in this repo, under a top-level `agent/`**, published to a marketplace from the subdirectory. The decisive argument is Tier-2 anti-drift: in-repo, CI regenerates the extract and fails on drift against the same export surface `exports_resolve.jl` already guards. Accepting the coupling to RD's release cadence is the price, and it is the right one — a separate repo is cleaner for distribution but re-opens the drift problem the docs charter closed.
2. **The structured-token seam: generic kind first, `host.jl` as the documented escape hatch.** A generic attribute-bag token kind keeps the common case data-only and inert; a typed `host.jl` remains available for studies that need real fields, presented for review *as code*. This mirrors how the engine treats `Invoke` (ADR 0011) — a sanctioned, named, visible exception rather than a blurred boundary. It implies a small upstream addition to RD (tracker W13), which must ship with semantic tests and a `log_token_fields` default like any other engine change.
3. **First increment: spine-thin across all four phases** (§10 step 1), over `Bash` + scripts, no MCP and no subagents, on a study whose answer is already known — re-deriving the marginal-scientist shadow price from an interview. That is the only setup in which a passing run distinguishes "the workflow found the right answer" from "the workflow found a confident one."

**Still open.**

4. **Whether `rd-doctor`'s checks belong upstream.** Several are general-purpose model hygiene, not agent-specific — a `diagnose(prob)` in `src/analysis.jl` would serve human modelers equally and would be covered by the semantic suite. Deferred until the check list has been exercised against real models; if it lands upstream it becomes CONTRACT §14 material and a new ADR.

## 12. Tracker

Legend: ✅ done · 🟡 in progress · ⬜ not started.

| Facet | ID | Deliverable | Status |
|---|---|---|---|
| Charter | — | This document | ✅ |
| Spine | W1 | Four phases minimally end-to-end on a known-answer study | ✅ `f8379b0` |
| Frame | W2 | `rd-frame` skill: slot table, question bank, qualifying test | ✅ |
| Author | W3 | `rd-author` skill: JSON authoring + validate/repair loop | ✅ |
| Diagnose | W4 | `rd-doctor`: the plausibility check list + `DOCTOR.md` | ✅ `f8379b0` |
| Report | W5 | `rd-report`: the signature briefing, delegating register and charts | ✅ |
| Knowledge | W6 | Tier-2 generated extract + the CI drift gate | ✅ `gen_reference.jl` + `.github/workflows/AgentLayer.yml` |
| Runtime | W7 | Julia CLI scripts (`validate`/`run`/`ensemble`/`draw`) | ✅ `f8379b0` |
| Runtime | W8 | MCP server `rd` with a persistent session + spec resources | ⬜ |
| Agents | W9 | `rd-modeler`, `rd-critic` (`rd-calibrator` phase 2) | ⬜ |
| Enforcement | W10 | `PostToolUse` validate hook on `*.rdj.json` | ⬜ |
| Surface | W11 | `card.html` self-contained model card | ⬜ |
| Packaging | W12 | `plugin.json`, commands, marketplace publication | ⬜ |
| Upstream | W13 | Generic attribute-bag token kind in `src/` (+ semantic tests, `log_token_fields` default) — the §11.2 data-only path | ⬜ |
| Upstream | W14 | `from_json_model` should instantiate `d["population"]` — today it validates the array and drops it, so a document's t=0 tokens materialize under `rdx` (which reads it and builds `PopulationEntry`s) and silently do not under a bare engine call | ⬜ |
| Upstream | W15 | Make `validate` total — a transition missing `"id"` throws out of `src/serialize.jl` instead of returning a `Diagnostic`, which costs the repair loop the JSON path it repairs against. `rdx` wraps the throw as a synthetic `:error` diagnostic; the engine should not need the wrapper | ⬜ |
| Upstream | W16 | **`validate` accepts a bare scalar as a node's CHILD; the constructor rejects it.** `_validate_node!` returns early for a non-dict (`src/serialize.jl:671`, *"a bare literal scalar — fine"*) but `node_from_dict` recurses unconditionally into an `AbstractDict` method (`:97`, `:99`, `:103`), so a number inside a `call` arg, a `sample` arg or a `choose` alt value passes validate with ZERO diagnostics and then throws `MethodError: no method matching node_from_dict(::Float64)`. The documented shorthand "a bare number is accepted for `const`" is therefore true at the top level and false one level down. One-line fix in the direction the reference already promises: `node_from_dict(x::Union{Real,Bool}) = Const(x)` | ⬜ |
| Upstream | W17 | **A `choose` node in an ATTRIBUTE position passes validate and then dies at compile.** `Choose` lowers to a `@choose(…)` marker that only the *reaction-line* resolver consumes (`src/interface/reaction_parser.jl`) — it is not a real macro, so a `choose` in `transitions[].rate` reaches `compile_attrs` and throws `UndefVarError: @choose not defined`. `choose` is legal in arc positions only, and `validate` does not know the difference. Either give `_validate_node!` the position context to reject it, or make `Choose` compile in an attribute | ⬜ |
| Surface | W18 | Captured-session renderer (`render_session.jl`) + the blind-capture protocol + one full capture | ✅ |

### W1 — what the spine actually delivered (2026-09-22, `f8379b0`)

`agent/studies/marginal_scientist/` reproduces the published shadow price **to the digit** — `rdx effect` returns `delta 18.991, se 2.528` against the case study's ≈ +\$19M ± \$2.5M — entirely through the JSON-document path, with `host.jl` (a token type and a registry) as the only Julia. That is the §11.3 acceptance criterion met: the agent-facing surface and the hand-authored surface are demonstrably the same engine, not an approximation of it.

Three findings worth carrying forward:

- **The gate discriminates.** On the `demo/bd_acquisition` fixture `rdx diagnose` reports 2 fails (five dead transitions; the wrong resource binding) and exits 1; on the study it reports 8 pass / 1 warn / 0 fail and exits 0. The checks were built to be capable of failing, and they do.
- **One warning was carried, not cleared.** The study's horizon truncates (two trials in flight at t=48), so `FINDINGS.md` states the headline as a **lower bound** rather than suppressing the warn. That is the intended discipline (§3.3) exercised once for real.
- **Two arms, two complete documents.** The convention is one self-contained validate-clean document per arm in the same directory, differing in as few lines as possible — here exactly one (`"init": 4.0` → `5.0`). Accepting the duplication buys a diff that *is* the intervention and an arm that replays forever on its own; a patch-and-overlay scheme would forfeit both.

### Captured sessions — the demo surface, and the best bug-finder so far (2026-09-22)

The four artifacts show what the workflow *produced*; they do not show what using it is like, which is what a prospective user needs. `agent/scripts/render_session.jl` renders a transcript JSON into one self-contained HTML page (phase bands, collapsed evidence blocks holding the exact commands, inlined SVGs, print stylesheet), and `agent/sessions/` holds the captures. The protocol is what makes a transcript evidence rather than marketing: two agents role-play, the operator **pre-commits its facts to a file before speaking** and the assistant starts from empty context and cannot read that file, so every number in the model was demonstrably elicited and every number in the result was produced by running the engine. Documented in `agent/sessions/README.md`, including what it does not prove.

Running one unscripted session on a NEW problem found more defects than any amount of re-reading the skills, because it exercised paths the known-answer study never touches:

- **W16 and W17, both upstream**, both instances of one theme already visible in W15: `validate` says yes and the constructor throws. The assistant hit W16 on its first `rdx run`, worked around it in its own document, and carried on — correct behaviour (`src/` is not the agent's to edit mid-study), but it means the class is only visible if someone reads the transcript. Reproduced and narrowed independently before recording: W16 is not `choose`-specific, it fires for a bare scalar under `call` and `sample` too.
- **A silent under-check in the gate itself.** `rdx`'s option parser stored options in a `Dict`, so `--tally=a --tally=b` — the natural spelling, and the one the skills' prose invites — kept only `b` and ran ONE check while the envelope still reported `pass`. A gate that quietly checks less than it was asked to is worse than no gate. Fixed: the four list-valued options (`metric`, `tally`, `expect_binds`, `expect_reusable`) accumulate, so both spellings are equivalent, and a *scalar* option repeated with two different values is now a hard envelope error rather than last-wins — last-wins on `--seed` is how you run with a seed you did not intend.
- **One check does not discriminate.** `--expect-binds=laptop` passed and so did `--expect-binds=tech`, the counter-hypothesis, on the same document. The assistant caught this itself and recorded it in `DOCTOR.md` rather than citing the pass as proof. `intended_resource_binds` currently tests *"does the named place run to zero"*, which several places can satisfy at once; it wants a relative measure (which pool is starved *most*, or for the largest fraction of the horizon) before it can rank two candidate constraints.

- **The failure mode is in the narrative, not the numbers — and no validator catches it.** The assistant honoured every `effect_vs_noise` finding the tool returned, then made four claims the tool had no opinion about: a return ratio formed by dividing two point estimates where one could not be separated from zero; a causal mechanism inferred from a single seed; an improvement in calibration between two figures 1.31 σ apart; and a ledger decomposition offered as evidence when it is an **identity** (`reward` is *defined* as a fixed linear combination of the two breach counts, so the decomposition closes for any two runs, signal or noise). All four were caught by a second reader applying the same ±2·se rule to the prose that `rdx` applies to the metrics, and all four were retracted in-transcript before the operator acted. **Implication for W9:** `rd-critic` is not optional polish — the checks it must run are on the *report*, not the model, and they are exactly the ones a schema validator cannot express. Candidate rules: no ratio whose denominator is a non-significant effect; no mechanism from n = 1; no comparison of two estimates without their difference's σ; no arithmetic identity cited as evidence.
- **A divergent baseline has no annual anything.** At the stressed duration the 4-laptop arm never settles, so every dollar figure measured against it scales with the horizon. The gate should flag this: if a run's terminal queue is still climbing, per-year framings of any difference against it are horizon-dependent and must be reported as such, or replaced by the settled/not-settled statement, which is horizon-free. Detectable mechanically — regress the queue on the last quartile of ticks and check the slope against zero.

## 13. Decisions log

- **2026-09-22 — Author models as JSON, not Julia.** The eval-free document plus the pure `validate` pass is the trust boundary and the repair loop's critic; generating DSL source would forfeit both. Host Julia is the documented exception, presented as code.
- **2026-09-22 — The knowledge layer cites and generates; it does not copy.** Extends the docs charter's §D4 anti-fork rule to this fourth consumer. Only procedural prose (the interview, the checks, the calibration discipline) is hand-written here.
- **2026-09-22 — Skill-and-script first, MCP second.** MCP is justified by session persistence and cross-client reach, not by capability; the portable layer must remain usable without it.
- **2026-09-22 — DIAGNOSE needs a skeptical reader, not only a check-runner.** A blind capture produced four claims that every automated gate passed and a second reader rejected, all of them ratios, mechanisms or comparisons asserted at a confidence the sample size did not support. `rdx` can only test what it can name; the statistical discipline the CLI applies to `metrics` has to be applied to the prose by something that reads the prose. This makes `rd-critic` (W9) load-bearing rather than a refinement, and it is why captured sessions are kept as a standing practice rather than a one-off demo.
- **2026-09-22 — A plausibility gate stands between the model and the headline number.** `validate` checks schema, not sense. Without §3.3 the workflow's most likely output is a confident wrong number, which would damage the framework more than no workflow at all.
- **2026-09-22 — The plugin lives in-repo under `agent/`; the token seam is generic-kind-first with `host.jl` as a named escape hatch; the first increment is spine-thin across all four phases on a known-answer study** (§11.1–§11.3, maintainer-ratified). The through-line of all three: keep the knowledge layer un-drifted by construction, keep the authored artifact inert by default, and validate the *sequence* before deepening any single phase.
- **2026-09-22 — One self-contained document per experimental arm, not a patch.** Arms live side by side in the study directory and differ in as few lines as possible. The duplication is the price of two properties worth more than it: the `diff` between two files *is* the intervention, reviewable at a glance by someone who knows nothing about the formalism, and every arm replays from `(model, seed)` on its own forever without needing the baseline it was derived from.
- **2026-09-22 — The known-answer test asserts at two strengths.** `agent/scripts/acceptance.jl` treats the *scientific* claim as hard (the delta lands inside the published band and clears 2·se) and *bit-reproducibility* against the recorded run as soft (reported, not fatal). Julia's `hash` is explicitly not stable across versions and the member seeds derive from it, so the seed stream can move without any semantics changing; blocking a merge on the last digits would train people to ignore the gate. The soft note still prints, because if the engine did not change, a moved stream is worth investigating.
- **2026-09-22 — The agent layer gets its own CI workflow, not a step in `CI.yml`.** `Agent layer` runs the reference drift gate and the known-answer study on one configuration, path-filtered to `agent/**` and `src/**`. It guards the agent-facing surface rather than the engine, so it does not belong on the version × OS matrix — and keeping it separate means the acceptance study's 40 seconds are paid once per relevant change rather than six times per push.
- **2026-09-22 — A captured session is a first-class artifact of this layer, and it is captured blind.** `agent/sessions/` holds rendered transcripts of whole sessions, produced by two agents role-playing under a pre-commitment protocol (the operator writes its facts down before speaking; the assistant cannot read them). Two reasons this earns its keep beyond the demo value: a transcript is the only artifact that shows the *interview*, which is where most of the workflow's value is created and none of it is currently documented; and an unscripted session on an unfamiliar problem is the layer's most productive test — the first one found two upstream engine defects (W16, W17), a silent under-check in `rdx`'s own option parsing, and a plausibility check that cannot discriminate between two candidate constraints.
- **2026-09-22 — A list-valued option accumulates; a scalar option repeated is an error.** `rdx` used a plain `Dict` for options, so repeating a flag silently dropped every value but the last — and every one of the list-valued flags drives a CHECK, so the failure mode was a gate reporting `pass` for a check that never ran. `metric`/`tally`/`expect_binds`/`expect_reusable` now accumulate (comma and repeated spellings are equivalent), and a scalar option given two different values fails loudly with exit 2 instead of picking one.
