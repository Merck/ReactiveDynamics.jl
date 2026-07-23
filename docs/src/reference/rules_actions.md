```@meta
CurrentModule = ReactiveDynamics
```

# Rules & actions

Rules are the endogenous decision channel: a `Rule` pairs a state-contingent guard with an ordered action program, so decisions (raise capital, acquire an asset, kill a program) live *in* the model rather than in host patch code. Actions are a closed, serializable family (`ActionStmt` subtypes) applied by [`apply_action!`](@ref); the engine fires enabled rules each step via [`fire_rules!`](@ref).

```@docs
Rule
ActionStmt
SetSpecies
SetParams
SetField
SetTokens
AddToken
Activate
Deactivate
Invoke
Log
Seq
apply_action!
fire_rules!
activate!
deactivate!
set_guard!
```

## Per-program ledger

Cost/reward/valuation actions accumulate into a per-program ledger — the model's accounting side. These read the ledger a finished (or in-progress) run has built up.

```@docs
ProgramLedger
program_ledger
program_ledger_entries
```
