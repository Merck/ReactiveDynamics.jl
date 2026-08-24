# Phase-0 semantic tests — Modality Truth Table (§1) and Genesis Modes (§2.8)
#
# Tiers: T1-characterization runs against the CURRENT engine (locks in behavior or pins a known
# bug via @test_broken); T2-acceptance encodes TARGET behavior and is wrapped (@test_skip + a
# commented reference block) because it names APIs that do not exist yet. Assertions were largely
# verified on Julia 1.12.5 during drafting; re-verify file:line citations before acting.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
using Statistics

@testset "Modality Truth Table (§1) and Genesis Modes (§2.8)" begin

    # [mod-row1-upfront-consumed-block-raw] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.3 truth-table row 1 ({} ⇒ upfront,consumed,block); §3.2 stage 3 + stage 8 ('Plain consumed tokens are not returned')
    # note: Verified live: material 100→98→96→94, widget→3, valuation_cost rows all 10.0. Locks in raw-
    # note: consumption semantics + the (:valuation_cost,t,scalar) ledger shape (solvers.jl:308-311).
    # note: placeInitVal/sol column order is spec order: material is col 1.
    @testset "Row 1 (upfront/consumed/block): empty modality is raw consumption at the arc multiplicity, never returned" begin
        net = @reaction_network begin
            @deterministic(1.0), 2 * material --> widget, name => build
        end
        @prob_init net material = 100 widget = 0
        @prob_params net
        @cost net material = 5.0
        prob = ReactionNetworkProblem(net, Dict(); tspan = 3, dt = 1.0)
        simulate(prob)
        df = prob.sol
        @test df.material[1] == 100.0
        # 2 material burned per tick, monotone non-increasing, never credited back
        @test df.material[2] == 98.0 && df.material[3] == 96.0 && df.material[4] == 94.0
        @test df.widget[4] == 3.0   # one widget emitted per completed (ct=0) tick
        # ledger: each spawn tick charges 2 units * cost 5.0 = 10.0 at placeCost
        costs = [r[3] for r in prob.log if r[1] == :valuation_cost]
        @test all(c -> c == 10.0, costs[1:3])
    end

    # [mod-row2-upfront-conserved-block-return] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.3 row 2 ({:conserved} ⇒ upfront,conserved,block); §3.2 stage 8 (':conserved tokens return q·multiplicity·(cycleTime if :rate else 1)'); finish! solvers.jl:438-442
    # note: Verified live (steady cash=94.0). Distinguishes conserved (returns; pool plateaus above 0) from raw-
    # note: consume (row 1, monotone drain). The exact floor 94 reflects the in-flight backlog at ct=3,rate=1;
    # note: assertion is on the plateau invariant + value.
    @testset "Row 2 (upfront/conserved/block): @conserved holds q·s and credits it back at finish" begin
        # steady-state proof of return: 1 conserved-holder spawned/tick, each holds 3 cash for ct=3 ticks.
        # Backlog of in-flight holders is bounded, so cash settles to a constant (held, not consumed).
        net = @reaction_network begin
            @deterministic(1.0), 3 * @conserved(cash) --> product, name => hold, cycletime => 3.0
        end
        @prob_init net cash = 100 product = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 12, dt = 1.0)
        simulate(prob)
        df = prob.sol
        # cash is HELD (debited at spawn) but RETURNED at finish: it reaches a steady floor, never drains to 0
        tail = df.cash[(end - 3):end]
        @test all(==(tail[1]), tail)            # constant in steady state
        @test tail[1] > 0                        # not consumed away (would hit 0 if raw-consumed)
        @test tail[1] == 94.0                    # observed steady floor: 100 - 2*3 (two cohorts mid-flight)
        @test df.product[end] > 0                # RHS still emitted on success
    end

    # [mod-row3-perstep-consumed-block-metered] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.3 row 3 ({:rate} ⇒ perstep,consumed,block); §1.5 (dt_scale=Δt, gated C>0); get_reqs_ongoing! solvers.jl:35-37
    # note: Verified live: fuel 1000→999→997→994→991→988→985 (draws 1,2,3,3,3,3). Metered/flow consumption is
    # note: dt-scaled and only accrues for ct>0 instances. Contrast mod-perstep-ct0-footgun.
    @testset "Row 3 (perstep/consumed/block): @rate draws q·s·Δt each ongoing tick, gated on cycletime>0" begin
        net = @reaction_network begin
            @deterministic(1.0), @rate(fuel) --> trip, name => drive, cycletime => 3.0
        end
        @prob_init net fuel = 1000 trip = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 6, dt = 1.0)
        simulate(prob)
        df = prob.sol
        # one new in-flight instance per tick; each draws 1*1*dt=1 fuel/tick while alive (ct=3)
        # so per-tick fuel draw ramps 1,2,3,3,... as the in-flight population builds toward 3
        draws = -diff(df.fuel)
        @test draws[1] == 1.0 && draws[2] == 2.0 && draws[3] == 3.0
        @test all(d -> d == 3.0, draws[3:end])   # saturates at 3 concurrent ct=3 instances
        @test df.fuel[1] == 1000.0
    end

    # [mod-row4-perstep-conserved-block-rented] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.3 row 4 ({:rate,:conserved}); §1.5 (return = +q·s·C, the per-tick integral); finish! solvers.jl:438-442 (:rate ⇒ ×transCycleTime)
    # note: Verified live: fuel settles at 998.0 (vs row-3 metered which keeps draining). The conserved return
    # note: uses the (:rate ? transCycleTime : 1) factor at solvers.jl:442. Stacked macro syntax
    # note: @rate(@conserved(fuel)) verified to attach BOTH tags (reaction_parser.jl:69-76 unions macro names
    # note: down the nesting).
    @testset "Row 4 (perstep/conserved/block): @rate+@conserved is a rented hold — drawn per tick, fully returned at finish" begin
        net = @reaction_network begin
            @deterministic(1.0), 2 * @rate(@conserved(fuel)) --> out, name => rc, cycletime => 2.0
        end
        @prob_init net fuel = 1000 out = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 8, dt = 1.0)
        simulate(prob)
        df = prob.sol
        # rented throughput: per-tick draws are exactly offset by q·s·C returns at finish, so the pool plateaus high
        tail = df.fuel[(end - 3):end]
        @test all(==(tail[1]), tail)        # steady state
        @test tail[1] == 998.0              # observed floor (one cohort's in-flight reservation)
        @test tail[1] > 990.0               # net hold is small vs raw consumption — characterizes 'rented'
        @test df.out[end] > 0
    end

    # [mod-row5-perstep-consumed-nonblock] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.3 row 5 ({:nonblock} ⇒ perstep,consumed,nonblock; token freed every step); §3.4 Invariant 1; FIXED solvers.jl:512
    # note: STAGE-A FIX (was a KNOWN BUG): free_blocked_places! previously referenced an undefined bare `q` and
    # note: raised `UndefVarError: q` on the 2nd tick once a :nonblock instance was in-flight. The free path now
    # note: credits `trans.q * tok.multiplicity` back every step, so the run completes. CRUCIAL: cycletime MUST be >0 to
    # note: keep an instance in-flight across a tick boundary and exercise free_blocked_places! at all.
    # note: Verified live: sensor=[10,9,8,8,8,8,8] (non-negative, finite) — the per-step free credit makes the
    # note: :nonblock pool plateau (held token returned each tick) rather than crash. Assertions are robust
    # note: invariants (no throw; non-negative; finite), not exact sensor numerals.
    # action: simulate(prob) — runs to completion
    @testset "Row 5 (perstep/consumed/nonblock): in-flight @nonblock token frees q·s every step and runs to completion" begin
        # cycletime>0 keeps a :nonblock instance in-flight across a tick boundary, so free_blocked_places!
        # (step 3 of _step!) iterates it; the freed :nonblock resource is credited back rather than crashing on `q`.
        net = @reaction_network begin
            @deterministic(1.0),
                @nonblock(sensor) --> reading,
                name => measure,
                cycletime => 3.0
        end
        @prob_init net sensor = 10 reading = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 5, dt = 1.0)
        @test (simulate(prob); true)              # FIXED: free_blocked_places! no longer hits undefined `q`
        df = prob.sol
        # the freed :nonblock resource is credited back every step ⇒ non-negative, finite trajectory
        @test all(>=(-1.0e-9), df.sensor)
        @test all(isfinite, df.sensor) && all(isfinite, df.reading)
    end

    # [mod-perstep-ct0-footgun] tier=T1-characterization expectedStatus=pass-now
    # contract: §1.4 illegal rule 'allocation = perstep requires transCycleTime > 0'; validate_modalities (solvers.jl); §2.4
    # note: PREVIOUSLY this foot-gun constructed and ran SILENTLY (fuel stayed flat at 100 while out grew,
    # note: because the @rate draw is gated on C>0 at build_requirements! ~117 and C defaulted to 0). The
    # note: §1.4 construction-time validator (validate_modalities) now REJECTS it up front, so the silent-
    # note: nothing behavior is unreachable — the foot-gun is caught, not run. This T1 now characterizes the
    # note: post-validator engine (construction throws); the T2 mod-construct-rejects-perstep-ct0 below is
    # note: the acceptance spec for the same model.
    @testset "Illegal: perstep (@rate) with cycletime=0 is rejected at construction (was a silent foot-gun)" begin
        net = @reaction_network begin
            @deterministic(1.0), @rate(fuel) --> out, name => r0
        end                                   # cycletime defaults to 0.0
        @prob_init net fuel = 100 out = 0
        @prob_params net
        # validate_modalities fires in the constructor: `allocation = perstep` requires cycletime > 0.
        @test_throws ArgumentError ReactionNetworkProblem(net, Dict(); tspan = 4, dt = 1.0)
    end

    # [mod-mode-macro-unions-modality] tier=T1-characterization expectedStatus=pass-now
    # contract: §5.4 placeDefaultModality; FIXED update.jl:108 (now uses `:placeDefaultModality` not bare `placeDefaultModality`)
    # note: STAGE-A FIX (was a KNOWN BUG): mode!/@mode previously raised `UndefVarError: placeDefaultModality` because
    # note: update.jl:108 referenced a bare `placeDefaultModality` instead of the column symbol `:placeDefaultModality`. @mode now
    # note: unions the named modality into the place's modality set. Verified live: after `@mode net X conserved`,
    # note: `net[1,:placeDefaultModality] == Set([:conserved])`. Note `net[1,:placeDefaultModality]` indexes the SCHEMA (net),
    # note: not the problem.
    # action: invoke `@mode net X conserved`
    @testset "@mode unions :conserved into the place's modality set (placeDefaultModality)" begin
        net = @reaction_network begin
            1.0, X --> Y, name => t1
        end
        @mode net X conserved
        @test :conserved in net[1, :placeDefaultModality]
    end

    # [mod-construct-rejects-nonblock-conserved] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: §1.4 illegal row ('blocking = nonblock requires return = consumed'); §1.1 typed Modality; replaces error at solvers.jl:461-465 with construction-time validation
    # note: Today this does NOT raise ArgumentError at construction — it constructs fine and then, on simulate,
    # note: throws `UndefVarError: q` from free_blocked_places! (solvers.jl:512) on the SECOND tick (verified),
    # note: i.e. the q-bug masks the intended :conserved+:nonblock error at solvers.jl:461-465, which is itself
    # note: never reached for in-flight tokens. So this @test_throws ArgumentError currently FAILS (wrong type,
    # note: wrong site, wrong time) — it encodes the §1.4 target of a single construction-time validation rule.
    # note: Will pass once the typed Modality validation is added to the constructor.
    # action: construct the problem and expect a construction-time rejection
    @testset "TARGET: nonblock+conserved rejected at CONSTRUCTION (not deep in finish!)" begin
        # A single token cannot be both held-until-finish and freed-every-step.
        net = @reaction_network begin
            @deterministic(1.0), @nonblock(@conserved(x)) --> y, name => bad, cycletime => 2.0
        end
        @prob_init net x = 10 y = 0
        @prob_params net
        # validate_modalities fires in ReactionNetworkProblem(...) (solvers.jl), before any step runs
        @test_throws ArgumentError ReactionNetworkProblem(net, Dict(); tspan = 5, dt = 1.0)
    end

    # [mod-construct-rejects-perstep-ct0] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: §1.4 illegal rule 'allocation = perstep requires transCycleTime > 0'
    # note: Today constructs and runs silently (see mod-perstep-ct0-footgun: fuel stays 100). Encodes §1.4's
    # note: construction-time validation. Pairs with the T1 characterization that pins the current silent-
    # note: nothing behavior.
    # action: construct and expect rejection
    @testset "TARGET: perstep modality with cycletime=0 rejected at construction" begin
        net = @reaction_network begin
            @deterministic(1.0), @rate(fuel) --> out, name => r0
        end                                  # cycletime defaults to 0.0
        @prob_init net fuel = 100 out = 0
        @prob_params net
        @test_throws ArgumentError ReactionNetworkProblem(net, Dict(); tspan = 4, dt = 1.0)
    end

    # [mod-construct-rejects-perstep-structured] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: §1.4 illegal rule 'allocation = perstep requires a non-structured (countable) place'; current deep error get_reqs_ongoing! solvers.jl:38-42
    # note: Today the structured+:rate clash only errors DEEP in get_reqs_ongoing! at simulate time
    # note: (solvers.jl:38-42), and `set_structured!` is a target helper (the engine currently sets structured-
    # note: ness via the `placeStructured` schema column / @structured token authoring, not a one-liner). This
    # note: test assumes the §1.1 typed re-model + a construction validator; it errors-until-implemented on both
    # note: the helper and the validation.
    # action: construct and expect rejection
    @testset "TARGET: perstep (@rate) on a structured/agent place rejected at construction" begin
        # Mark a place structured, then try to draw it per-step. TARGET: construction-time rejection.
        net = @reaction_network begin
            @deterministic(1.0), @rate(robot) --> task, name => run, cycletime => 3.0
        end
        @prob_init net robot = 5 task = 0
        @prob_params net
        # Mark `robot` structured via the ACTUAL engine mechanism: register_token_kind!
        # sets the `placeStructured` schema flag (the reference block's `set_structured!` was a
        # hypothetical target helper — this is the real API, per ADR 0006/0007 authoring).
        register_token_kind!(net, :robot)
        @test_throws ArgumentError ReactionNetworkProblem(net, Dict(); tspan = 4, dt = 1.0)
    end

    # [gen-poisson-source-dt-invariant] tier=T1-characterization expectedStatus=pass-now
    # contract: §2.8 poisson mode (Source, expand_rate Poisson(dt·rate) create.jl:151); §2.3 dt-invariance on the Poisson path; §2.6
    # note: Verified live: m1≈101.7, m2≈100.3. Empty LHS (∅) bypasses the upfront gate (reqs==0 ⇒ Inf,
    # note: solvers.jl:121) so the Poisson proposal is the realized count — a pure source. atol is a sampling
    # note: tolerance for n=200 (Var≈100 ⇒ SE of mean≈0.7; 8.0 is generous). Uses Random.seed! on the GLOBAL rng
    # note: because no per-state seeding exists yet (§4); a later D6 seed= kwarg should replace it.
    # action: ensemble of 200 runs at dt=1.0 and dt=0.5, compare means against rate·tspan = 100
    @testset "Genesis `poisson`: empty-LHS source is dt-invariant in expectation (ensemble mean ≈ rate·tspan)" begin
        function total_spawned(dt; seed, rate = 2.0, tspan = 50.0)
            Random.seed!(seed)            # NOTE: relies on GLOBAL rng today (no seeding API — §4 D2 unmet)
            net = @reaction_network begin
                2.0, ∅ --> arrival, name => inflow
            end
            @prob_init net arrival = 0
            @prob_params net
            prob = ReactionNetworkProblem(net, Dict(); tspan = tspan, dt = dt)
            simulate(prob)
            return prob.u[1]
        end
        m1 = mean(total_spawned(1.0; seed = s) for s in 1:200)
        m2 = mean(total_spawned(0.5; seed = s) for s in 1:200)
        @test isapprox(m1, 100.0; atol = 8.0)     # E = rate*tspan = 2*50
        @test isapprox(m2, 100.0; atol = 8.0)
        @test isapprox(m1, m2; atol = 8.0)        # halving dt preserves the expected total
    end

    # [gen-floor-dt-invariant] tier=T1-characterization expectedStatus=pass-now (WS-3 fix landed)
    # contract: §2.3 'Discretization hazard — ceil on the spawn count'; §2.8 'scheduled/batch counts MUST be integer-valued'; FIX solvers.jl (floor, not ceil)
    # note: Previously `qs .= ceil.(Int, qs)` rounded a fractional @deterministic count 0.3→1 EVERY tick
    # note: regardless of dt, so halving dt (doubling ticks) roughly DOUBLED the total (t1≈11, t2≈21) —
    # note: an upward bias as dt→0. WS-3 replaced it with `floor` (§2.3 'condition ceil to the poisson
    # note: path'): the Poisson path is already integer so floor is a no-op there, while a fractional
    # note: @deterministic count truncates to 0 (0.3→0) identically at every dt. Result: dt-invariant.
    # note: Integer @deterministic counts (the only in-contract form, §2.8) are unaffected — verified by
    # note: the integer sub-case below, which spawns the same whole total at dt=1.0 and dt=0.5.
    # action: at dt=1.0 vs dt=0.5, the spawned total is now equal (both 0 for the fractional count).
    @testset "Genesis `floor` path: @deterministic count is dt-invariant (WS-3 fix)" begin
        # Fractional deterministic count: floor(0.3)=0 at BOTH dt=1.0 and dt=0.5 — no upward bias, invariant.
        function frac_total(dt; tspan = 10.0)
            net = @reaction_network begin
                @deterministic(0.3), ∅ --> a, name => src
            end
            @prob_init net a = 0
            @prob_params net
            prob = ReactionNetworkProblem(net, Dict(); tspan = tspan, dt = dt)
            simulate(prob)
            return prob.u[1]
        end
        t1 = frac_total(1.0)
        t2 = frac_total(0.5)
        @test isapprox(t1, t2; atol = 1.0)   # dt-invariance: ceil→floor fix (was the @test_broken pin)
        @test t1 ≈ 0.0                       # floor(0.3)=0, so nothing spawns (the old ceil forced 1/tick)
        # Integer deterministic count: floor(2)=2, so the fix does NOT suppress a legitimate whole count.
        # (@deterministic is a bare per-tick count with no dt-scaling per §2.3, so its TOTAL is inherently
        # dt-dependent — finer dt ⇒ more ticks ⇒ larger total; that is correct, not the bug. We assert only
        # that integer counts still spawn, i.e. floor left them intact.)
        function int_total(dt; tspan = 10.0)
            net = @reaction_network begin
                @deterministic(2), ∅ --> a, name => src
            end
            @prob_init net a = 0
            @prob_params net
            prob = ReactionNetworkProblem(net, Dict(); tspan = tspan, dt = dt)
            simulate(prob)
            return prob.u[1]
        end
        i1 = int_total(1.0)
        i2 = int_total(0.5)
        @test i1 > 0                         # floor(2)=2: integer counts are unaffected by the fix
        @test i2 > i1                        # finer dt ⇒ more ticks ⇒ larger total for a per-tick count
    end

    # [gen-flow-triggered-zero-then-positive] tier=T1-characterization expectedStatus=pass-now
    # contract: §2.8 flow mode (Routing; EXISTING upfront-LHS gate solvers.jl:110-128); 'flow-triggered genesis works today with zero new mechanism'
    # note: Verified live: product 0,0,0,2,4,6,8,10 while feed plateaus at 2. The upfront-LHS gate (reqs>0 ⇒
    # note: floor(alloc/multiplicity), solvers.jl:121-124) clamps the rate-100 proposal to available feed tokens — the
    # note: §2.8 'flow' idiom (high nominal rate + upstream place as consumed LHS). Demonstrates token-bounded
    # note: firing min(proposal, tokens) = tokens.
    @testset "Genesis `flow`: non-empty upfront LHS spawns ZERO when input empty, fires once upstream deposits tokens" begin
        net = @reaction_network begin
            @deterministic(2.0), ∅ --> feed, name => upstream
            @deterministic(100.0), feed --> product, name => router
        end
        @prob_init net feed = 0 product = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 5, dt = 1.0)
        simulate(prob)
        df = prob.sol
        # t=0: feed=0 ⇒ router (nominal rate 100) is TOKEN-GATED to 0; product stays 0 through the first interval
        @test df.product[1] == 0.0 && df.product[2] == 0.0
        # once feed accumulates, router fires bounded by available tokens (2/tick deposited, 2/tick routed)
        @test df.product[end] > 0
        @test all(>=(0.0), df.feed)         # gate never debits below available tokens
        # realized routing == upstream deposit rate, NOT the nominal rate=100
        @test maximum(-diff(df.product)) <= 2.0 + 1.0e-9
    end

    # [gen-capacity-overflow-deferral] tier=T1-characterization expectedStatus=pass-now
    # contract: §2.8 capacity mode + add_to_spawn! deferral; §3.4 Invariant 3; FIXED state.jl:251-256
    # note: STAGE-A FIX (was a KNOWN BUG with TWO defects in add_to_spawn!): `findfirst(pred, length(vec))` passed
    # note: an Int not a range (raised a `MethodError`), and on match it did `:transHash += n` (Symbol += Float64)
    # note: instead of `:transToSpawn += n`. With both repaired, the over-capacity overflow path runs: the proposal
    # note: (3/tick) exceeds capacity (5) while instances are in-flight (ct=10), and the surplus is DEFERRED via
    # note: add_to_spawn! rather than crashing. Verified live: simulate completes and the live concurrent count is
    # note: exactly 5 (≤ capacity). Characterizes the working Invariant-3 deferral (live count never exceeds
    # note: capacity, overflow carried forward).
    # action: simulate(prob) — runs to completion
    @testset "Genesis `capacity`: over-capacity proposal is deferred via add_to_spawn!, live count bounded by capacity" begin
        # proposal (3/tick) eventually exceeds capacity (5) while instances are in-flight (ct=10),
        # triggering the overflow-deferral path add_to_spawn!.
        net = @reaction_network begin
            @deterministic(3.0),
                1 * @conserved(slot) --> job,
                name => start,
                cycletime => 10.0,
                capacity => 5
        end
        @prob_init net slot = 100 job = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 8, dt = 1.0)
        @test simulate(prob) !== nothing   # FIXED: add_to_spawn! deferral no longer hits MethodError / Symbol +=
        # Invariant 3: live concurrent instances never exceed capacity; overflow is carried forward, not dropped.
        @test count(t -> t[:transHash] == prob[1, :transHash], prob.ongoing_firings) <=
            5
    end

    # [gen-capacity-clamp-no-overflow] tier=T1-characterization expectedStatus=pass-now
    # contract: §2.8 capacity mode (transCapacity gate solvers.jl:147-153); §3.4 Invariant 3 (gate present)
    # note: Verified-adjacent: same gate path as gen-capacity-overflow-deferral-bug but kept under capacity so
    # note: add_to_spawn! is never called. Characterizes the working clamp half of Invariant 3. @rate(fuel) with
    # note: ct=3 keeps instances alive long enough to accumulate live count without consuming the upfront pool
    # note: (fuel=1000 generous).
    @testset "Genesis `capacity`: when proposal ≤ capacity, concurrent instances are bounded with no deferral" begin
        # proposal (1/tick) never exceeds capacity (3); live count rises to 3 (ct=3) and holds — no add_to_spawn! crash.
        net = @reaction_network begin
            @deterministic(1.0),
                @rate(fuel) --> job,
                name => start,
                cycletime => 3.0,
                capacity => 3
        end
        @prob_init net fuel = 1000 job = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 8, dt = 1.0)
        # inspect prob.ongoing_firings and the :new_transitions log
        @test simulate(prob) !== nothing                       # does NOT hit the deferral bug
        h = prob[1, :transHash]
        @test count(t -> t[:transHash] == h, prob.ongoing_firings) <= 3   # capacity bound holds
        # each tick proposes exactly 1 (deterministic), within capacity, so no overflow is ever deferred
        spawncounts = [
            v for r in prob.log if r[1] == :new_transitions for
                (hh, v) in r[3:end] if hh == h
        ]
        @test all(c -> c <= 3, spawncounts)
    end

    # [gen-scheduled-periodic-calendar] tier=T1-characterization expectedStatus=pass-now
    # contract: §2.8 scheduled mode (Source; '@deterministic bare-count + @periodic idiom'); periodic() state.jl:239 + compilers.jl:67; §2.7
    # note: Verified live: cohort 0,0,0,3,3,6,6,9 over t=0..8 — steps of 3 at boundaries (one save-row lag from
    # note: periodic() needing length(sol.t)>1, state.jl:241). IMPORTANT: idiom must use the MACRO
    # note: `@periodic(2.0)` inside the rate (compiles to periodic(state,2.0) via the reserved-name rewrite
    # note: compilers.jl:98); the bare CALL `periodic(2.0)` errors (MethodError: needs (state,period)). This is
    # note: the BD acquisition-lever idiom and sidesteps the event channel by design.
    @testset "Genesis `scheduled`: @deterministic(N*@periodic(p)) fires N spawns at each calendar boundary" begin
        # scheduled idiom (NOT the event channel): rate is 0 except at multiples of period, where it is N.
        net = @reaction_network begin
            @deterministic(3 * @periodic(2.0)), ∅ --> cohort, name => intake
        end
        @prob_init net cohort = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 7, dt = 1.0)
        simulate(prob)
        df = prob.sol
        # spawns occur only at period boundaries (t=2,4,6): cohort steps up by 3 there, flat between
        deltas = diff(df.cohort)
        @test count(d -> d == 3.0, deltas) == 3        # exactly three boundary firings over tspan=7
        @test count(d -> d == 0.0, deltas) == length(deltas) - 3
        @test df.cohort[end] == 9.0                    # 3 boundaries * 3 each
    end

    # [gen-event-action-repaired] tier=T1-characterization expectedStatus=pass-now
    # contract: §3.4 Invariant 7 (repaired); ADR 0010 §A — the scheduled in-model lever via the event/Rule channel
    # note: Stage B repairs the event channel: each :E row is lifted at construction to a Rule
    # note: {guard=trigger, action=RawExpr(action), every_tick}, and fire_rules! (src/actions.jl, _step! step
    # note: 10) evaluates the guard and RUNS the action (no longer the bare no-op fetch). A scheduled
    # note: injection `(@t() > T) && (budget += N)` now takes effect on schedule — the §2.8 acquisition-lever
    # note: idiom is live (and the typed Rule/SetMarking form in rules_decisions.jl is the preferred surface).
    # note: `budget` must be a real place column (it is set by the action), so it appears on an inert
    # note: transition to enter :S.
    @testset "Event action takes effect on schedule — Invariant 7 met (the in-model lever)" begin
        # A scheduled budget injection expressed as an event: from t>2 the action adds 999 each tick.
        net = @reaction_network begin
            0.0, budget --> budget, name => budget_holder   # inert; declares `budget` as a place
            0.0, raw --> product, name => t1
            (@t() > 2.0) && (budget += 999.0)
        end
        @prob_init net raw = 10 product = 0 budget = 0
        @prob_params net
        prob = ReactionNetworkProblem(net, Dict(); tspan = 5, dt = 1.0, seed = 1)
        simulate(prob)
        # Invariant 7: the event fired on schedule, so budget rose from 0 once t crossed 2.
        @test prob.sol[!, "budget"][end] > 0.0
        @test prob.sol[!, "budget"][findfirst(<=(2.0), prob.sol[!, "t"])] == 0.0  # nothing before the trigger
    end
end
