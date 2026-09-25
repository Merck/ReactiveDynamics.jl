The option the operator took to ExCo, *reserve now and review at Kestrel's readout*, is one JSON document: `arm_cmo_exit.rdj.json`. Because it is data, not code, every sentence the operator said can be traced to the element it became. The map below shows the whole model at the level of its parts. The rows after it follow one quote at a time: what the operator said, the element as the document holds it, and what that element does. The document stores each expression as a tree of typed nodes; the rows print those trees as formulas so they can be read. In this study the documents were emitted by a generator script (`gen_models.py`) and each was then validated as data, so what is shown here is the validated document, not the script.

!figure agent/studies/launch_supply/figures/model_map.svg | The model behind the recommendation, one box per part. Generated from `arm_cmo_exit.rdj.json`: all 48 places, 31 transitions and 27 rules, each in exactly one box. Solid arrows are the model's arcs (the number is how many), dotted arrows are rates that read the live state, dashed arrows are rules setting a marking, or a parameter that a rate reads.

> *Operator, turn 3:* “Effective is about 13 after changeovers and two 2-week shutdowns a year; the last 12 months averaged 12.8.”

```
place line · starts at 128

transition k_base_own
  takes    line [conserved] + cash×15
  gives    k_stock + own_used
  rate     r_k_base_own
  holds    1 month(s)
  succeeds 0.97
  priority 100
```

The line is a pool of 128 tenths of a batch per month: the measured 12.8, not the nominal 16. Every batch filled on it, here Kestrel's launch stock, takes one tenth for the month it runs and gives it back after (`conserved`). That is why Kestrel, Tern and the marketed products can collide on it. The `cash×15` is the $150k variable cost per batch, charged in $1k units.

> *Operator, turn 3:* “the fill site has to be in the BLA with three PPQ batches before we file”

```
rule t_file  (once)
  when  (((t_ok > 0.5) && (t >= ((14 + t_shift) + 3))) && !(t_filed > 0.5)) && (((cmo_now > 0.5) && (cmo_ppq > 0.5)) || (!(cmo_now > 0.5) && (t_ppq_done >= 30)))
  then  t_filed := 1
        t_file_t := t

rule t_cmo_pas  (once)
  when  ((!(cmo_now > 0.5) && (cmo_ppq > 0.5)) && (t_filed > 0.5)) && (t >= (max(cmo_ready_t, (t_file_t + 12)) + 6))
  then  t_cmo_ok := 1
```

The constraint the operator said people forget became a guard. Tern files only when its readout succeeded, three months have passed, and a fill site has its three PPQ batches, either at the CMO or on the own line (`t_ppq_done >= 30` tenths). If the CMO qualifies after the filing, `t_cmo_pas` adds it six months later, the Prior Approval Supplement from the same answer.

> *Operator, turn 14:* “I'd move Tern's PPQ in-house, file Tern on our own line, and give the CMO notice.”

```
rule cmo_exit_k  (once)
  when  (((cmo_now > 0.5) && (exit_at_k > 0.5)) && !(k_ok > 0.5)) && (t >= (5 + k_slip))
  then  cmo_active := 0
        cmo := 0
        cmo_ppq := 0
        cmo_now := 0
        cmo_exited := 1

param exit_at_k = 1    (the base model has 0: this is the option)

transition cmo_fee_y2
  takes    cash×5000
  gives    nothing (a pure cost)
  rate     1 * ((((cmo_now > 0.5) || (cmo_exited > 0.5)) && (t == 14)) && ((!(!(t_ok > 0.5) && (!(cmo_reassign > 0.5) || !(k_ok > 0.5))) && !(cmo_exited > 0.5)) || (walk2 > 0.5)))
  priority 1
```

The operator's rule, added in turn 15, is one rule and one switch. At Kestrel's readout (month 5 plus the drawn slip), if Kestrel failed, the reservation is dropped and Tern's PPQ falls back to the own line through the `t_file` guard above. The forfeited year's fee is still charged. `exit_at_k` is the only line that differs from the rigid reservation, `arm_cmo.rdj.json`. That is what the inertness check tested: on the same seeds, the two options came out identical in every world where Tern does not succeed alone.

+++~ 7 more: what else the operator said, and where it went

> *Operator, turn 3:* “Kestrel is 60% and Tern is 45%”

```
param p_k = 0.6

param p_t = 0.45

rule draws  (once)
  when  t >= 0
  then  k_ok := 1 * Bernoulli(p_k)
        t_ok := 1 * Bernoulli(p_t)
        k_slip := floor(Uniform(0, 3))
        t_shift := floor(Uniform(-3, 4))
        cmo_lead := floor(Uniform(12, 16))
        build_lead := floor(Uniform(30, 37))
        cmo2_lead := floor(Uniform(12, 16))
        k_success := k_ok
        t_success := t_ok
```

The odds are parameters. The readout outcomes and every lead time the operator gave as a range are drawn once, at the start, from the seeded stream. That is what lets the same seed give every option the same readout dates and lead times, the paired comparison the findings rest on.

> *Operator, turn 3:* “about 3% of batches get rejected”

```
transition t_make_own
  takes    line [conserved] + cash×15
  gives    t_stock + own_used
  rate     r_t_make_own
  holds    1 month(s)
  succeeds 0.97
  priority 1
```

The reject rate is the `succeeds 0.97` on every fill transition: a rejected batch still used its line time and its cost, but gives no stock.

> *Operator, turn 5:* “Model it as competing for the same batches, not as free.”

```
transition t_ppq_own
  takes    line [conserved] + cash×15
  gives    t_ppq_done + own_used
  rate     r_ppq
  holds    1 month(s)
  succeeds 0.97
  priority 1
```

Tern's in-house PPQ is a fill transition like any other. It draws on the same `line` pool as Kestrel's launch stock, so the collision the operator predicted is something the engine produces, not an assumption written in.

> *Operator, turn 9:* “I would never starve a launch to top up a buffer.”

```
transition k_base_own
  takes    line [conserved] + cash×15
  gives    k_stock + own_used
  rate     r_k_base_own
  holds    1 month(s)
  succeeds 0.97
  priority 100

transition k_buf_own
  takes    line [conserved] + cash×15
  gives    k_stock + own_used
  rate     r_k_buf_own
  holds    1 month(s)
  succeeds 0.97
  priority 0.5

param k_floor = 2
```

Kestrel's stock is split into two transitions over the same pool. The part up to the floor (`k_floor`, 2 months of forward stock) runs at priority 100, and the buffer above it at priority 0.5, below Tern's batches at 1. When the line is short, the buffer waits and the launch does not.

> *Operator, turn 3:* “$5M a year for up to 4 batches a month”

```
transition cmo_fee
  takes    cash×5000
  gives    nothing (a pure cost)
  rate     1 * (((cmo_active > 0.5) && ((t >= cmo_start) && ((12 * floor(((t - cmo_start) / 12))) == (t - cmo_start)))) && !(((cmo_now > 0.5) || (cmo_exited > 0.5)) && (t == 14)))
  priority 1

rule cmo_ready  (once)
  when  ((cmo_active > 0.5) && !(cmo_resigned > 0.5)) && (t >= (cmo_start + cmo_lead))
  then  cmo := 40
        cmo_ppq := 1
        cmo_ready_t := t
```

The fee is a transition that takes 5000 ($5M in $1k units) from cash once a year while the contract is active. When the tech transfer finishes, `cmo_ready` fills the CMO pool with 40 tenths, the 4 batches a month.

> *Operator, turn 7:* “Lost. In immunology an unsupplied patient goes to a competitor or stays on their current therapy.”

```
transition k_ship
  takes    k_orders + k_stock
  gives    k_sold
  rate     1000
  priority 10

transition k_lose
  takes    k_orders
  gives    k_lost
  rate     1000
  priority 0
```

An order that finds stock ships at priority 10 and earns the batch's value. One that does not is moved to `k_lost` the same month, never queued. That is the operator's answer, and also the engine's semantics: unmet demand is dropped, so it has to be counted explicitly.

> *Operator, turn 1:* “Pay a CMO a reservation fee now to hold a fill-finish slot.”

```
param cmo_now = 1    (the base model has 0: this is the option)

rule cmo_sign_now  (once)
  when  (cmo_now > 0.5) && (t >= 1)
  then  cmo_active := 1
        cmo_start := 2
```

The option itself is a switch. `cmo_now = 1` against the base model's 0 is the whole difference between reserving and waiting, apart from `exit_at_k` above. Every other element is shared, so a difference in the result can only come from the decision.

+++

*Generated by `julia --project=. agent/scripts/render_mirror.jl agent/studies/launch_supply/mirror.json`. Every quote is checked verbatim against the transcript, and every element is printed from the model document.*
