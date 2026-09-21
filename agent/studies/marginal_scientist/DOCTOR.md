# DOCTOR — plausibility record

> The DIAGNOSE artifact (AGENT_WORKFLOW_CHARTER §3.3). A schema-clean model can still be nonsense; this is the record that it isn't, and it is the precondition for quoting any number. Every line below is machine-produced evidence, not an opinion.

**Command** (from the repository root):

```bash
julia --project=. agent/scripts/rdx.jl diagnose \
  agent/studies/marginal_scientist/model.rdj.json \
  --host=agent/studies/marginal_scientist/host.jl \
  --seed=1 --expect-binds=scientist --expect-reusable=scientist --tally=launch
```

**Verdict: 8 pass · 1 warn · 0 fail** → exit 0, clear to report. The `--expect-*` flags are the brief's pre-registered intent (`BRIEF.md`, "What could go wrong"), so these are not self-graded checks: the tool compares the run against what the brief claimed *before* the run existed.

| Check | Verdict | Evidence |
|---|---|---|
| `transitions_fire` | pass | all 4 transitions spawned — `intake_onc` 49, `intake_imm` 49, `trial_onc` 39, `trial_imm` 13 |
| `completions_exist` | pass | 126 successful completions over 48 months |
| `contention_exists` | pass | `scientist` runs to a trough of **0** |
| `intended_resource_binds` | pass | the starved pool is `scientist` — the resource the brief named, so the headline is a hiring number |
| `modality_intent` | pass | every `scientist` input arc carries `return: conserved` — the bench is occupied and returned, not spent |
| `horizon_sufficient` | **warn** | 2 trials still in flight at t=48; `Project` +35% and `launch` +35% across the final quarter |
| `bounded_pools` | pass | `launch` exempt as a declared tally (`--tally=launch`); nothing else grows unbounded |
| `ledger_wired` | pass | reward 112,000 \$k, cost 3,920 \$k |
| `stochasticity` | pass | both trials carry `prob_of_success` 0.55 — outcomes are genuinely uncertain |

## The one warning, and why it stands

The horizon truncates: at month 48 two trials are mid-cycle and the candidate queue is still growing at roughly a third per quarter. This is **scope, not error** — the study asks what the fifth scientist is worth *over a 48-month plan*, which is the horizon a hiring decision is actually made over. The consequence is stated rather than fixed: the reported shadow price is the value realized within 48 months and therefore a *lower bound* on the value of relieving the bench, since candidates still queued at the horizon carry option value the number does not count. Extending `meta.tspan` would raise it; the decision does not need it to.

The 35%-per-quarter growth in `Project` is the queue the study is about, and it is the reason the bench binds. The intake rate (1.5/month per franchise) genuinely outruns a four-person bench running four-month trials — that is the modeled world, not a mis-scaling.

## Contention is real but not visible in saturation

`min_saturation` is 1.0: every trial that *started* got its full per-tick claim. The contention shows up one step earlier, as trials that could not start at all — which is why the `scientist` trough at 0 and the growing `Project` queue are the operative evidence here. A `@conserved` bench rations by refusing admission, not by under-serving what is already running; a `@rate` resource is the one that shows up as fractional saturation. Both signals are checked, and it only takes one.
