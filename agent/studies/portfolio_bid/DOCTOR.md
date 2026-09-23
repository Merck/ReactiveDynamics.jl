# DOCTOR — plausibility record

> The DIAGNOSE artifact (AGENT_WORKFLOW_CHARTER §3.3). `validate` reports zero diagnostics on all nine documents, which proves they are well-formed and proves nothing about whether they are about the VP's system. This is the record that they are, and it is the precondition for quoting any number in [`FINDINGS.md`](FINDINGS.md). Every line below is machine-produced evidence; the commands and the raw envelopes are catalogued in [`README.md`](README.md), the brief that pre-registered the expectations is [`BRIEF.md`](BRIEF.md), and the workflow is [`../../README.md`](../../README.md).

**Commands** (from the repository root). Two arms, because an arm that pulls the lever can change which resource binds — and here it changes which transitions exist at all:

```bash
julia --project=. agent/scripts/rdx.jl diagnose agent/studies/portfolio_bid/model.rdj.json \
    --seed=1 --expect-binds=fte --expect-reusable=fte --tally=own_market

julia --project=. agent/scripts/rdx.jl diagnose agent/studies/portfolio_bid/arm_both.rdj.json \
    --seed=1 --expect-binds=fte --expect-reusable=fte --tally=own_market,acq_market
```

**Verdict: `model` 7 pass · 1 warn · 1 fail → exit 1. `arm_both` 8 pass · 1 warn · 0 fail → exit 0.** The `--expect-*` flags are the brief's pre-registered intent, including the expectation that the baseline would fail `transitions_fire`, so none of this is self-graded: the tool compares the run against what the brief claimed before the run existed.

| Check | `model` | `arm_both` | Evidence |
|---|---|---|---|
| `transitions_fire` | **fail** | pass | baseline: `adv_acq_p2`, `adv_acq_p3` and `adv_acq_filed` never spawn — there is nothing acquired to advance. `arm_both`: all 8 transitions spawn. |
| `completions_exist` | pass | pass | work finishes on both sides: over 120 seeds the baseline averages 2.842 own launches per decade, `arm_both` 2.808 own plus 0.975 acquired |
| `contention_exists` | pass | pass | `fte` is exhausted, not merely busy: on seed 1 the baseline leaves ≤ 1 FTE free on **96 of 122 ticks (79%)**, `arm_both` on 52 of 122 (43%) |
| `intended_resource_binds` | pass | pass | the starved pool is `fte` — the resource the brief named. The exec map paints 9 of 11 places as contended; `cash` is one of the two that is not, so the headline is a capacity number and not a financing one |
| `modality_intent` | pass | pass | every `fte` arc is `allocation: upfront`, `return: conserved`, `blocking: block` — people are occupied for the cycle time and returned, and a programme that cannot claim its people does not start. `cash` is `perstep`/`consumed`: money does not come back |
| `horizon_sufficient` | **warn** | **warn** | **7** firings (baseline) and **9** (`arm_both`) still in flight at month 120, with pools still trending through the final quarter |
| `bounded_pools` | pass | pass | `own_market` (and `acq_market` on `arm_both`) exempt as declared tallies; the preclinical backlog peaks at 25.7 (`max:own_p1`) and is not flagged — a queue under a ceiling, not a leak |
| `ledger_wired` | pass | pass | both sides of the ledger move: `cash` metered at 0.6–3.0 \$M/month per in-flight programme against launches credited 900 (own) / 1400 (acquired) \$M on production |
| `stochasticity` | pass | pass | every `adv_*` carries a probability of success strictly inside (0, 1) — 0.38 to 0.90 — and `discovery` is Poisson at 0.35/month |

## The one failure, and why it is explained rather than fixed

`transitions_fire` fails on `model` because the three acquired-pipeline transitions never spawn. That is the counterfactual arm behaving exactly as a counterfactual must: in a world where we do not bid, there are no in-licensed assets, and a transition that advances them has nothing to advance. The same check passes on `arm_both` with all eight transitions firing, which is the evidence that the dead transitions are dead for the modelled reason and not because an arc asks for a pool that is empty or a rate that is zero. The engine's rule is that exit `1` forbids quoting a number; the resolution here is not to suppress the check but to pair it — every conclusion drawn about the acquired pipeline is drawn from arms in which that pipeline is alive and fully exercised.

The general point is worth more than this instance. **A study that reports "diagnose is green" without naming the arm it ran on is hiding something.** Green on the treated arm and green on the control are different claims, and in a counterfactual design they are usually not both available; the honest artifact is a per-arm verdict line, which is why there are two above.

## The warning both arms carry, and its consequence

Neither arm finishes its work inside the window: 7 and 9 firings are still in flight at month 120, and the pools have not flattened. **Every eNPV in this study is therefore horizon-truncated and a lower bound** — a Phase III started in year nine draws its cash inside the window and books its launch outside it, so it enters the ledger as pure cost. This is scope rather than error, because a bid is decided over a planning horizon and not over steady state, but it has one consequence sharp enough that it must travel into the report: it is the reason **+24 FTE measures *less* than +12** (+306 ± 116 against +370 ± 106). Twenty-four extra people genuinely start more work; the extra launches simply fall past month 120 where the ledger cannot book them. Read that pair as a horizon artifact, not as diminishing returns — and in either reading, buy the cheaper dose.

## Contention is real, and saturation is the wrong gauge for it

A `conserved` resource rations by **refusing admission**: a Phase II that cannot claim its six people never starts, so everything that did start was fully served and fractional saturation looks healthy. The operative evidence is therefore the FTE trough and the queue, not a utilization percentage — `fte` sitting at ≤ 1 free for 79% of the baseline's ticks, and a preclinical backlog that peaks near 26 in a pipeline that opened with 4.

There is a trap here that this study walked into and had to route around. The *baseline* is the most saturated arm of the three traced: 96 of 122 ticks near the ceiling, against 81 with the deal and 52 with the deal plus 12 FTE. Utilization falls as the system gets less constrained, which is backwards, and the reason is fragmentation — with no acquired work competing, the group packs tightly into its own pipeline, whereas a 6-FTE Phase II start cannot fit into the gaps a 40-FTE ceiling leaves once bigger programmes are resident. So the argument in the report is made from queue-wait and launch counts, never from a saturation percentage. If a diagnostic number moves the wrong way for a mechanical reason, the fix is to stop quoting that number, not to explain it.

## What the hand pass found

Magnitudes are believable for the world described: ≈4 new preclinical candidates a year from `discovery`, 2.8 own launches per decade out of eleven in-flight assets plus intake, and a \$6,000M budget that is never the thing anyone waits for. The single-seed traces agree qualitatively across arms about which resource binds, and they make the constraint concrete: the three acquired assets land in `acq_p2` at month 8 and begin Phase II in months 25, 31 and 39 at today's capacity — the last of them idle on the shelf for 31 months. Fund twelve more people and two of the three start in month 9. Clear to report.
