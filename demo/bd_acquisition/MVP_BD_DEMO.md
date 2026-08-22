# MVP — Acquisition Impact on a Pharma Pipeline Portfolio (BD demo)

> **Status — Phase-0.5 design artifact, 2026-06-21.** A demo specification written against the contract ([CONTRACT_DRAFT.md](../../spec/CONTRACT_DRAFT.md) §1–§12) and ADRs 0001–0010. Its second purpose is adversarial: to surface limitations of the draft contract *before* Phase-1 implementation. Findings are catalogued in §8; the endogenous-decision-channel finding (B) has been promoted to [ADR 0010](../../spec/adr/0010-rules-and-conditional-transitions.md) + CONTRACT §12, and the concurrent Phase-0.5 increment (ADRs 0007/0008/0009 → CONTRACT §10/§9.5/§11) resolves or softens findings C, F, G and the declarative-initial-portfolio gap (cross-referenced inline). Domain is pharma throughout (the engine's demonstrators are pharma; rNPV is the canonical BD metric), but the structure is domain-general.

## 1. The thesis — what the demo sells

The artifact is **a reproducible counterfactual on a living portfolio**: run the same pipeline under the same conditions, once *without* an acquisition and once *with* it, and read off the change in portfolio value, launches, and capital need *attributable to the deal* — including the deal's effect on the programs the company already owns. BD partners do not need another spreadsheet rNPV; they need *"given my pipeline and my finite scientists and budget, what does this acquisition do to the whole system over time?"* That is a dynamic, resource-contended, stochastic, decision-laden question — exactly what this engine is for. The demo proves the framework is a reusable artifact: one host file + one JSON model + a scenario grid yields a defensible deal analysis.

This is the ADR-0006 Business-Development north-star made concrete, extended with the ADR-0010 endogenous decision channel so that *acquisition decisions themselves can be state-contingent rules inside the model*, not just exogenous host interventions.

## 2. What to model

| Element | RD representation | Notes |
|---|---|---|
| **Program / asset** | a **structured token** of one `ProjectToken <: AbstractStructuredToken` KIND, identity preserved across phases (§9.5 phase-as-attribute) | descriptor fields: `phase`, `npv_peak`, `pos_remaining`, `cost_to_date`, `therapeutic_area`, `acquired::Bool`, `acq_time` |
| **Pipeline phase** | a `phase` ATTRIBUTE value `∈ {:Discovery, :Phase1, :Phase2, :Phase3, :Filed, :Market}` on the single `:Project` KIND (§9.5/ADR 0008 §D, the maintainer-canonical encoding) | not a separate species kind; portfolio counts per phase come from `@select`/`ntokens` reductions (§9.3) over the `phase` field |
| **Phase advance** | a **transition**, `genesis: flow` (≥1 upfront LHS, §2.8): selects an in-phase token via `@select(Project, phase==:PhaseN)` (`TokenPredicate`, §9.5) + `scientists` (`@conserved`) + `budget` (`@rate`) for a `cycletime`; on `Binomial(q, PoS)` success (solvers.jl:413) advances it via `@advance(phase, :PhaseN1)` (lowers to `SetField{:phase,…}`, §9.5) | failure / `maxlifetime` timeout ⇒ program **soft-retired** `:removed` (shelved/killed); `past_bonds` keeps phase history for the audit trail |
| **Scientists / capacity** | plain `@conserved` species (held during a phase, returned on completion) | the scarce resource the ADR-0002 allocator rations by `priority` |
| **Budget / cash** | plain `@rate` (burn) + `@conserved` reserve species; replenished by a `scheduled` financing transition or a `every_tick` rule | drives the capital ledger (`:valuation_cost`, solvers.jl:308-311) |
| **Initial portfolio** | declarative `population[]` array in the model JSON (CONTRACT §10.3): `{species, kind, instances:[…]}` of `ProjectToken`s across phases | makes the starting pipeline part of the reproducible input (§8.2 S2), not imperative host code — closes the gap noted at finding **H** |
| **The acquisition** | **endogenous**: an ADR-0010 **rule** whose action is `Seq[AddToken(ProjectToken, …), SetParams(synergy…), SetSpecies(scientists,+…)]`, `fire_mode: once`. The lever lives *in the model JSON* and may be scheduled (`@t() > T_acq`) **or** state-contingent (`portfolio_proxy < X`) | NOT the legacy host patch and NOT the old no-op event channel; `AddToken` lowers to `add_structured_token!` (index-safe, ADR 0006 §9.1) |

### 2.1 The four synergy mechanisms (each independently toggle-able)

This is what makes it a BD *tool* rather than a "+N programs" toy. Each synergy is gated by a single `synergy_*` param the acquisition rule flips, so a scenario is a choice of which synergies are on.

1. **Pipeline synergy** — more shots on goal (the injected programs). Baseline of any deal.
2. **Resource synergy** — the target brings headcount/capital: the rule's `SetSpecies(scientists, +Δ)` / `SetParams(budget_rate, ·)` bumps the pools at acquisition.
3. **Capability / PoS synergy** — the target's platform raises `prob_of_success` of *existing* same-area programs: the relevant transition's `prob_of_success` ExprNode reads `Ref{param} synergy_pos`, which the rule flips.
4. **Operational-efficiency synergy** — shared CRO/process shortens `cycletime`: the transition's `cycletime` ExprNode reads `Ref{param} synergy_eff`.

**Milestone-1 rule: synergies are param-mediated** (read a `Ref{param}` a rule sets). Token-pool-mediated synergy (a `prob_of_success` reading sibling counts via `TokenAgg`) is now *deterministic* — CONTRACT §9.5 (ADR 0008) pins the measurable observation point that ADR 0006 left open (finding C downgraded) — but param-mediation stays the Milestone-1 default for simplicity; token-mediated synergy is a clean Milestone-2 extension. Note finding **D-bis**: PoS/cycletime synergy affects only programs that *enter* a phase after the deal, not those mid-flight, because those attributes are frozen at spawn — a deliberate modeling choice the demo must state.

**Phase representation — phase-as-attribute (canonical, maintainer ruling).** Per CONTRACT §9.5 / ADR 0008 §D, a token's lifecycle state is a TOKEN ATTRIBUTE, not a species KIND: one `:Project` KIND carries a `phase` field, selected by `@select(Project, phase==:Phase2 && npv>θ)` and advanced by `@advance(phase, :Phase3)` (which lowers to `SetField{:phase, …}`, the field-write generalization of `@move`). The per-KIND-per-phase encoding survives only as the degenerate `species`-field predicate (`@move` back-compat), not a co-equal alternative — no kind explosion, continuous predicates like `npv>θ` work, uniform `past_bonds` history on one identity. The MVP uses phase-as-attribute throughout; the table above names phases for readability, but they are `phase`-field values of one `ProjectToken` kind, advance transitions select their input via the §9.5 `TokenPredicate` (`@select`) and advance it via `@advance`/`SetField`.

## 3. Granularity, and "refined dynamics substituted"

- **Coarse (the MVP you show):** 6 phases × 1 advance-transition each, ~4 plain resource species, tokens with ~6 fields, ~3 rules. Authored compactly with the §11.4 `@pipeline` sugar (the chain expands to `flow`-genesis guarded transitions with the phase-as-attribute idiom). The whole model is **one `model.rdj.json` + one host file** (the `ProjectToken` type/constructor/`priority` override). Compact, expressive, LLM-authorable.
- **Refined (the granularity-substitution demo):** use the §11.2 `refine(spec, :t_p2_advance, phase2_submodel; ports)` operator to replace the coarse Phase-2 advance with a sub-pipeline `{enroll → dose → readout}` carrying its own resource draws and interim kill points — boundary-matched at the in/out port species by FK-repoint (§7.4/J7), so every other transition is untouched (§11.5 plug-compatibility). Then check coarse vs refined agree on aggregate rNPV within ensemble CI, with §11.3 advisory diagnostics (`coarse.cycletime ≈ Σ critical-path`, `coarse.pos ≈ Π sub-PoS`) flagging where they shouldn't. This is *granularity substitution* as a first-class operator — **finding F is resolved by §11/ADR 0009** (it was OPEN when this demo was first drafted). Caveat (§11.4): refinement is process-structural; a token hosting its OWN sub-network (entity-level refinement) is deferred to a future ADR.

## 4. Scenarios — the experiment grid

| Scenario | Acquisition | Synergies |
|---|---|---|
| **S0** Baseline | none | — |
| **S1** Deal, pipeline-only | inject M programs | (1) |
| **S2** + resource | " | (1,2) |
| **S3** + capability/PoS | " | (1,3) |
| **S4** + op-efficiency | " | (1,4) |
| **S5** Full | " | (1–4) |

Plus two decision sweeps: **timing** (`T_acq` early vs late, or the state-contingent trigger threshold) and **price** (acquisition cost → breakeven). Every cell is run as an **ensemble of N seeds** (§4 D8 per-member seeding from a root seed) → distributions, not point estimates. Because a run is fully determined by `(model.rdj.json, seed)` (§8.2 S2), the grid is a clean reproducible matrix: each cell is one `(model_hash, scenario_params, seed)` triple.

### 4.1 Counterfactual semantics (decision locked 2026-06-21)

The headline number is Δ-rNPV = rNPV(Sₖ) − rNPV(S0). **We compute it as an ensemble-averaged difference of means, not a per-seed paired difference**, and we keep the engine's single state-owned RNG stream as-is (§4 D2/D5). Rationale and the known cost are in §8 finding A: a per-seed paired Δ is contaminated because S1 makes *more* draws than S0 after the deal, desynchronising the shared RNG so every organic program's draws shift too. The ensemble *mean* Δ is unbiased; we accept the wider CI (more seeds) in exchange for not extending §4 now. **Future work (documented, not scheduled):** per-entity RNG substreams keyed by `(root_seed, entity_id)` to recover common-random-numbers variance reduction and a meaningful per-seed paired Δ — now recorded as CONTRACT **§4.6 (D-future)** and the ADR-0010 open question. **Partial mitigation available now:** CONTRACT §10.5 `dump_state`/`restore` lets you fork a SINGLE run at the lever tick and apply the lever to only one copy — the two branches share an identical pre-lever history and RNG state up to the fork, so the divergence is strictly the lever's causal cone. This gives a cleaner paired comparison than two independent runs (though draws still desync post-fork); the demo uses it for the illustrative single-seed "before/after" narrative while reporting the ensemble Δ as the headline.

## 5. Outcome metrics

- **Portfolio rNPV** (headline) — end-of-run and as a time series: `Σ_active pos_remaining · npv_peak · discount(t) − cost_to_date − acq_price`.
- **Δ-rNPV (treatment effect)** = rNPV(Sₖ) − rNPV(S0), ensemble-averaged (§4.1).
- **# launches** (tokens reaching `:Market`), **time-to-first-launch**, **P(≥1 launch)**.
- **Peak capital requirement** (max cumulative budget draw) — the financing ask.
- **Scientist utilization over time** — does the deal *starve* the organic pipeline? An honest **diseconomy** the demo surfaces: acquired late-stage tokens whose `priority` override (ADR 0006) wins scarce scientists crowd out organic Phase-1 starts.
- **Synergy decomposition** — rNPV(S5) vs Σ of individual synergy contributions → super/sub-additivity.
- **Breakeven price** — from the price sweep: "don't pay more than X."

Discounting and rNPV roll-ups are **post-processing** over the ledger Arrow streams (§8.5) plus token-query snapshots (§9.3 `sumattr`/`active`). The engine itself does no discounting (§8 finding D).

## 6. Implementation against the contract

| MVP element | Contract / ADR mechanism |
|---|---|
| Model artifact | one `model.rdj.json` — §8.1, ADR 0005 document shape |
| `ProjectToken` type/constructor/`priority` override | **host Julia**, never serialized; `from_json(…; registry=Dict(:ProjectToken=>…))` — ADR 0006 §B/§C |
| Compact authoring | `@pipeline` chain sugar + `@process` fragment — §11.4, ADR 0009 |
| Phase representation | one `:Project` KIND, `phase` attribute, phase-as-attribute (canonical) — §9.5, ADR 0008 §D |
| Declarative initial portfolio | `population[]` array (`{species,kind,instances}`) — §10.3, ADR 0007; makes the start pipeline part of the reproducible input (§8.2 S2) |
| Phase advance | `genesis: flow` (≥1 upfront LHS, validated) — §2.8; input selected by `@select(Project, phase==…)` `TokenPredicate` — §9.5; advanced by `@advance(phase,…)` → `SetField{:phase,…}` preserving identity + `past_bonds` — §9.5/§9.1 |
| Success / kill / shelve | `Binomial(q,PoS)` §3.2; timeout ⇒ q=0 §2.5; soft-retire `:removed` §9.1 |
| **Acquisition lever (endogenous)** | **ADR 0010 rule**: guard ExprNode + `Seq[AddToken, SetParams, SetSpecies]`, `fire_mode: once` — CONTRACT §12 |
| **Conditional pipeline lines** | **ADR 0010** `guard::ExprNode` on the transition (e.g. Phase-3 start fires iff `cash > phase3_cost`) — §12 |
| Scarce-resource rationing | ADR 0002 weighted progressive filling; acquired-token `priority` override |
| Determinism / counterfactual | §4 D1/D6 (seed), D8/D9 (ensemble); ensemble-averaged Δ (§4.1); fork-at-tick via §10.5 `dump_state`/`restore` |
| Metrics | ledger Arrow streams §8.5 + token reductions §9.3 |
| Granularity substitution | `refine`/`abstract` boundary-matched FK-splice — §11.2, ADR 0009 |
| "Combined company from t=0" view | authoring-time `@compose`/`@join` + `equalize` to pool shared `budget`/`scientists` — §7/§11.4 — distinct from the runtime lever (§8 finding E) |

## 7. Recommended Milestone-1 scope cut

Build, in order:
1. Host `ProjectToken` (type, constructor, `priority` override) + the per-network registry plumbing (ADR 0006 §C).
2. Coarse `model.rdj.json`: one `:Project` KIND (phase-as-attribute §9.5), declarative `population[]` initial portfolio (§10.3), `@pipeline` chain of guarded `flow` advances (`@select` input / `@advance` output), **param-mediated** synergies, and the three rules (financing `every_tick`; acquisition `once`; one conditional Phase-3 guard).
3. Host driver (`from_json` → `simulate`); the lever is *inside* the model, so the driver only sets scenario params + seed. For the single-seed before/after narrative, `dump_state` at the lever tick → `restore` twice → arm the rule on one branch (§10.5).
4. Post-processor: ledger Arrow + token snapshots → rNPV / Δ / launches / capital / utilization.

**Defer to Milestone 2:** the §11.2 `refine` granularity demo (build it once the coarse model is validated), TokenAgg-mediated synergy (now deterministic via §9.5 but not needed for M1), per-tick structured `getobservable` portfolio time-series (§9.4), the `@compose`-two-companies path.

**Hard Phase-1 prerequisites (bug list + Phase-0.5 ADRs):** Inv-6 prune fix (solvers.jl:501 — else shelved programs re-emit RHS and corrupt the ledger every tick); token-bind bugs (solvers.jl:288/452/455/476); §4 D5 RNG threading; **the ADR-0010 event-eval repair (solvers.jl:323) and guard skip-gate (state.jl:179)** — without these the endogenous lever does not fire; the ADR-0007 `reinit!` completeness fix (rebuild the §10.3 token population, reset the RNG and the `once`-rule latches) — without it the ensemble re-runs are not reproducible; the §9.5 `TokenPredicate`/`@advance`(`SetField`) surface — without it phase-as-attribute advance is not expressible.

## 8. Contract weak-points this exercise surfaces

The exercise's primary deliverable. Status tags: **RESOLVED** (folded into a contract change this round), **DEFERRED** (decision taken to not act now, documented), **OPEN** (flagged for a future contract increment).

**A. Counterfactual at fixed seed has no clean home in §4. — DEFERRED.** Δ-rNPV(S1−S0) per-seed is the money shot, but the two runs diverge in *number of RNG draws* after the deal (S1 has more programs ⇒ more Poisson/Binomial), so the single state-owned stream (§4 D2/D5) desynchronises and every organic program's draws shift — the per-seed Δ conflates the deal effect with mean-zero RNG reshuffling. §4 gives independent ensemble *members* (D8/D9) but **no per-entity substreams within a trajectory**, so common-random-numbers variance reduction is unavailable. *Decision (2026-06-21):* keep the single RNG; use ensemble-averaged Δ (§4.1); record per-entity substreams keyed by `(root_seed, entity_id)` as §4 future work.

**B. No endogenous, state-contingent decision channel — only exogenous host patches. — RESOLVED (ADR 0010 + §12).** Events were a no-op (Inv-7, solvers.jl:323) and §8.6 made the lever a host-applied patch — scheduled only, never "acquire *if* portfolio rNPV drops below X." ADR 0010 repairs the event channel into a **Rule** (guard ExprNode + closed action set), adds a stateless **`guard::ExprNode`** on transitions (AND-ed with the latching `transActivated`, state.jl:179), extends the action set to the full endogenous set `{SetSpecies, SetParams, AddToken, Activate, Deactivate, Log, Seq}`, and adds `fire_mode {every_tick, once}`. The acquisition lever is now in-model, serializable, reproducible, and LLM-authorable.

**C. The most expressive synergy (cross-entity reads) had an unresolved determinism hole. — RESOLVED (§9.5 / ADR 0008).** "My Phase-3 PoS rises because the target gave me a regulatory team" wants a `prob_of_success` reading sibling token counts via `TokenAgg`; ADR 0006's open question was the `TokenAgg` observation point. §9.5 (ADR 0008) pins it: token-reading subexpressions (in a `TokenPredicate` AND a `TokenAgg`) observe the pre-`evolve!` reflected counts (`update_u_structured!`, solvers.jl:649/653) and MUST be 𝓕ₜ-measurable (no `Sample` node), so D1 reproducibility is now well-defined. *Residual:* the demo still defaults to param-mediated synergy for M1 simplicity (§2.1); token-mediated synergy is a clean, now-deterministic M2 extension.

**D. Two non-composing accounting systems; no per-tick hook to accrue economics onto an in-flight token. — OPEN (recommend doc clarification + evaluate `on_step`).** The plain-species ledger (§8.5 `:valuation*`) is pool-level and per-tick; per-program economics (`npv_estimate`, `cost_to_date`, `pos_remaining`) live per-token. They do not auto-reconcile, and rNPV/discounting is entirely **post-hoc** (the engine never discounts, though §8.5 speaks of "the discounted streams"). And pre/post actions fire only at **spawn and finish** (§3.2 stages 3, 9) — there is **no per-tick hook on an ongoing instance** — so "accrue burn onto the bound `ProjectToken` each tick" is not expressible; you accrue at the plain `budget` ledger and reattribute to programs in post via `past_bonds` + phase durations. *Recommendation:* state in the contract that per-program economics are reconstructed in post-processing, and evaluate whether the action model needs a per-tick `on_step` hook — the ADR-0010 action infrastructure is the natural home (its open question #4 already flags this).

**D-bis. TVE re-evaluation semantics are inconsistent, and it changes what "operational-efficiency synergy" means. — OPEN (recommend §5 split).** §5 marks `cycletime`/`prob_of_success`/`max_lifetime` "TVE? = yes," but `finish!` reads them from the instance's **frozen spawn-time snapshot** (`trans_[:transCycleTime]`/`[:transProbOfSuccess]`, solvers.jl:409-413), whereas `rate`/`priority` are genuinely re-read per tick. So a mid-run synergy lowering `cycletime` or raising `PoS` affects **only programs entering the phase after the deal, not those mid-flight** — a real, surprising consequence the demo must choose deliberately. *Recommendation:* split §5's TVE column into `{recipe-per-tick: rate, priority, multiplier}` vs `{instance-frozen-at-spawn: cycletime, prob_of_success, max_lifetime}` and pin the freezing rule.

**E. "Acquisition" has two incompatible readings. — DEFERRED (use injection at runtime).** Structural merge of two pipelines is **authoring-time** `@compose`/`@join`/`equalize`, and is **forbidden on a live model** (J8 — `rem_parts!` reindexes, breaking frozen closures); the mid-run lever is **append-only token injection** (§9.1, ADR 0010 `AddToken`). So you cannot truly merge two companies' models (pooling a shared `budget` via `equalize`) *mid-run* — the merge must precede `simulate`. Also `@join` currently **drops `:E` and `:obs`** (J4) — now scheduled to be fixed en route by §11.4 `@compose` (closes J4 and the J9 `include_model` bug). *Decision:* the runtime lever is `AddToken`; `@compose`/`@join` is reserved for the "combined company from t=0" view.

**F. Compositionality without refinement. — RESOLVED (§11 / ADR 0009).** §7 was flat name-merge, not hierarchical/morphic; there was no "this transition *is* a sub-network" with a coarse/fine consistency guarantee (REVIEW.md's "structure is implicit, not morphic"). §11 adds open ports (§11.1) + `refine`/`abstract` as a boundary-matched FK-splice (§11.2, reusing the §7.4/J7 reactant-promotion repoint), with advisory coarse/fine consistency diagnostics (§11.3) and `@pipeline`/`@process`/`@compose` compact authoring (§11.4). The granularity demo (§3) now exhibits a refinement *operator*, not just two hand-built models. *Residual:* entity-level refinement (a token hosting its own sub-network) is deferred to a future ADR (§11.4).

**G. Per-tick portfolio KPIs are not first-class. — PARTIALLY RESOLVED (§9.4 names the path).** A time series of token-aggregated portfolio value needs `getobservable` over the structured container (§9.4 specifies implementing it as the coupled-agent read path); the observation point is now pinned (§9.5, ex-finding C). The remaining work is implementing `getobservable` on the container (still flagged unimplemented in §9.4). Until then Milestone-1 reconstructs portfolio time-series in post; end-of-run token reductions are first-class today.

**H. Structured-token initial state was not declarative. — RESOLVED (§10.3 / ADR 0007).** The starting portfolio was built by imperative host code (`add_structured_token!` in a loop), so it lived outside the model document and broke §8.2 S2 ("(model.json, seed) determines the run") for structured models. §10.3 adds a declarative serializable `population[]` array (the structured analogue of `placeInitVal`), instantiated at construction before t=0 with seeded attribute draws and per-species creation indices. The MVP's initial pipeline is now reproducible input.

**I. A Rule action could not write a token field over a `@select`-ed set, and actions could not run general code. — RESOLVED (ADR 0011 + §12.3 amendment).** §9.5/ADR 0008 added `SetField`/`@advance` to write the *firing instance's bound token* (a transition-post-action idiom), but a standalone Rule could not perturb a *population* of tokens, and the closed verb set could not express arbitrary callback logic. ADR 0011 closes both, eval-free: **(1) `SetTokens{predicate, assigns}`** — a population-level field write over a §9.5 `TokenPredicate`, legal in a Rule (it carries its own predicate, needs no bound token), iterated in the §9.2 `(species, creation_index)` order — "write down all Phase-2 oncology `pos_remaining` by 10% on a competitor readout" is one declarative action; **(2) `Invoke{fn, args}`** — the general-code escape hatch, lowering to `registry[fn](state, transition, …)` (the ADR-0006 §C registry-by-name boundary in statement position), so a callback can run arbitrary host queries/mutations/logic with the file still carrying only a name (§8.4 S4 holds verbatim; author obligations O1–O4 cover determinism/append-only/tick-boundary/purity). The action family is now `{SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}`. *Trade made explicit:* declarative verbs stay statically validatable + LLM-authorable (the default); `Invoke` is the trusted-but-unverified tier for logic that exceeds them.

## 9. Decisions taken (2026-06-21)

1. **Counterfactual:** ensemble-averaged Δ; keep single RNG; per-entity substreams = documented future work (finding A → CONTRACT §4.6).
2. **Decision channel:** endogenous rules/callbacks + conditional transitions are **core** — promoted to ADR 0010 + CONTRACT §12 (finding B). Ratified forks: guard field **and** rule-driven Activate/Deactivate; `fire_mode {every_tick, once}`.
3. **Action expressivity:** actions/callbacks may write a `@select`-ed token population (`SetTokens`) and run general host code (`Invoke`, via the eval-free ADR-0006 registry) — ADR 0011 + §12.3 (findings I + the general-code ask). Full family `{SetSpecies, SetParams, SetField, SetTokens, AddToken, Activate, Deactivate, Invoke, Log, Seq}`.
4. **Phase representation:** phase-as-attribute is canonical (§9.5/ADR 0008 §D); `@select` selects, `@advance`/`SetField` advances.
5. **Domain:** pharma throughout.

## 10. Next steps

- (this round) ADR 0010 (§12) + ADR 0011 (§12.3/§8.4 amendment) + CONTRACT §4.6 written; ADR index updated; findings re-statused against the Phase-0.5 increment (B/C/F/H/I resolved, G partial, A/D/D-bis open).
- (Phase-1, gated) implement the ADR-0010/0011 surface + the §7 prerequisite bug fixes; then scaffold `host/ProjectToken.jl` + `models/pipeline.rdj.json` + `analysis/rnpv.jl` so each §8 finding becomes an executable test against the contract. The `Invoke` escape hatch needs a registered statement-callback calling-convention tag (ADR 0006 §C) and a `validate` arity check; `SetTokens` reuses the ADR-0008 predicate evaluator.
