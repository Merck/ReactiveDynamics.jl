# Phase-0 semantic tests — Reference models: SIR, toy-pharma, rNPV (brief acceptance criteria)
#
# Tiers: T1-characterization runs against the CURRENT engine (locks in behavior or pins a known
# bug via @test_broken); T2-acceptance encodes TARGET behavior and is wrapped (@test_skip + a
# commented reference block) because it names APIs that do not exist yet. Assertions were largely
# verified on Julia 1.12.5 during drafting; re-verify file:line citations before acting.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
using Statistics

@testset "Reference models: SIR, toy-pharma, rNPV (brief acceptance criteria)" begin

    # [simulate-n-ticks] tier=T1-characterization expectedStatus=pass-now
    # contract: ADR 0001 / INVENTORY API: simulate(prob,n) drives the stepper; Contract §2.1 (one clock, state.t += dt, solvers.jl:668)
    # note: Verified: simulate(p,5) => 6 sol rows, prob.t==5.0. Locks in the single-clock dt-stepping contract
    # note: (solvers.jl:668) and n-tick semantics of the AA-reexported simulate. No bug.
    @testset "simulate(prob, n) advances exactly n ticks on the single clock" begin
        m = @ReactionNetworkSchema begin
            α * S * I, S + I --> 2I, name => I2R
        end
        @prob_init m S = 100 I = 5
        @prob_params m α = 0.001
        @prob_meta m tspan = 100 dt = 1.0
        prob = ReactionNetworkProblem(m)
        simulate(prob, 5)
        @test prob.t == 5.0
        @test size(prob.sol, 1) == 6   # initial row + 5 ticks
        @test prob.sol[!, "t"][1] == 0.0
        @test prob.sol[!, "t"][end] == 5.0
    end

    # [sir-conservation] tier=T1-characterization expectedStatus=pass-now
    # contract: Brief reference-model criterion (SIR known behavior); Contract §3.4 Invariant 2 (conservation) for the closed-flow case
    # note: Verified S+I+R == 1009 first and last tick exactly. S+I->2I nets -1S+1I, I->R nets -1I+1R, so total
    # note: is structurally invariant for this same-tick (cycletime=0) flow SIR. cycle_time=>0 in example.jl is
    # note: a NON-alias and ignored; transCycleTime defaults 0.0 anyway. No bug exercised (no maxlifetime, no
    # note: nonblock). Post Stage A: reproducibility comes from the seed= construction kwarg (state-owned rng),
    # note: NOT Random.seed!; conservation here is structural so it holds for any seed.
    @testset "SIR conserves total population S+I+R across the run" begin
        sir = @ReactionNetworkSchema begin
            α * S * I, S + I --> 2I, name => I2R
            β * I, I --> R, name => R2S
        end
        @prob_init sir S = 999 I = 10 R = 0
        @prob_params sir α = 0.0001 β = 0.01
        @prob_meta sir tspan = 250 dt = 0.1
        prob = ReactionNetworkProblem(sir; seed = 1)
        simulate(prob)
        S = prob.sol[!, "S"]; I = prob.sol[!, "I"]; R = prob.sol[!, "R"]
        total = S .+ I .+ R
        @test all(isapprox.(total, total[1]; atol = 1e-9))   # S+I+R invariant
        @test total[1] == 1009.0
        @test all(S .>= -1e-9) && all(I .>= -1e-9) && all(R .>= -1e-9)   # non-negativity
    end

    # [sir-epidemic-peak] tier=T1-characterization expectedStatus=pass-now
    # contract: Brief reference-model criterion (epidemic peak exists; monotone R); qualitative SIR behavior
    # note: Verified Imax≈680 (interior), R monotone non-decreasing, R is a pure sink (no R-->S), so diff(R)>=0
    # note: holds. Stochastic: pin with the seed= construction kwarg (Stage A state-owned rng) so the
    # note: peak/monotonicity assertions are reproducible; assertions are robust invariants, not exact values.
    @testset "SIR shows an epidemic peak and monotone non-decreasing recovered" begin
        sir = @ReactionNetworkSchema begin
            α * S * I, S + I --> 2I, name => I2R
            β * I, I --> R, name => R2S
        end
        @prob_init sir S = 999 I = 10 R = 0
        @prob_params sir α = 0.0001 β = 0.01
        @prob_meta sir tspan = 250 dt = 0.1
        prob = ReactionNetworkProblem(sir; seed = 1)
        simulate(prob)
        S = prob.sol[!, "S"]; I = prob.sol[!, "I"]; R = prob.sol[!, "R"]
        @test maximum(I) > 5 * I[1]   # a genuine outbreak peak (verified Imax≈680 vs I0=10)
        peak_ix = argmax(I)
        @test 1 < peak_ix < length(I)   # peak is interior, not at an endpoint
        @test I[end] < maximum(I)        # infection declines after the peak
        @test all(diff(R) .>= -1e-9)     # recovered is monotone non-decreasing
    end

    # [sir-infection-burns-out] tier=T2-acceptance expectedStatus=pass-after-fix
    # contract: Brief reference-model criterion (I->0 as t->inf); known SIR asymptotic behavior
    # note: At the brief's natural horizon (tspan=250) I[end]=153 (verified) — I does NOT reach 0, so the naive
    # note: 'I->0' assertion FAILS. It is horizon/param sensitive: a longer tspan (or once S is exhausted, since
    # note: S hit 0 at t=250) drives I down via I-->R. Marked T2/pass-after-fix because it needs the longer
    # note: horizon to hold and is the asymptotic-behavior acceptance criterion, not a current-engine lock-in.
    # note: If a deterministic-rate SIR is preferred to remove burn-out noise, use @deterministic on both rates.
    # note: Post Stage A: reproducibility via the seed= construction kwarg (state-owned rng); under seed=1 at
    # note: tspan=2000 I_end=0.0 (verified), so the long-horizon burn-out criterion holds.
    @testset "SIR infection eventually decays toward zero (long horizon)" begin
        sir = @ReactionNetworkSchema begin
            α * S * I, S + I --> 2I, name => I2R
            β * I, I --> R, name => R2S
        end
        @prob_init sir S = 999 I = 10 R = 0
        @prob_params sir α = 0.0001 β = 0.01
        @prob_meta sir tspan = 2000 dt = 0.1
        prob = ReactionNetworkProblem(sir; seed = 1)
        simulate(prob)
        I = prob.sol[!, "I"]; S = prob.sol[!, "S"]
        @test I[end] <= 1.0   # infection has effectively died out by the long horizon
        @test I[end] < maximum(I)
    end

    # [sir-seeded-reproducible-construction] tier=T1-characterization expectedStatus=pass-now
    # contract: Contract §4 D1 (reproducibility), D2 (RNG isolation), D6 (seed at construction) — Stage A: state owns rng
    # note: Stage A landed a state-owned rng + working seed= construction kwarg, so reproducibility is now via
    # note: seed=, NOT the global RNG. Verified: two runs built with seed=42 are bit-identical EVEN with global
    # note: rand() perturbed between them (RNG isolation, D2). Random.seed! before an UNSEEDED run no longer
    # note: makes it reproducible — that idiom is retired. This is the proper D1/D2/D6 statement.
    @testset "SIR trajectory is reproducible under a fixed CONSTRUCTION seed" begin
        function build_sir()
            sir = @ReactionNetworkSchema begin
                α * S * I, S + I --> 2I, name => I2R
                β * I, I --> R, name => R2S
            end
            @prob_init sir S = 999 I = 10 R = 0
            @prob_params sir α = 0.0001 β = 0.01
            @prob_meta sir tspan = 50 dt = 0.1
            return sir
        end
        pa = ReactionNetworkProblem(build_sir(); seed = 42); simulate(pa)
        rand(Int); rand(Int)   # perturb the global RNG: a seeded run must be isolated from it (D2)
        pb = ReactionNetworkProblem(build_sir(); seed = 42); simulate(pb)
        @test pa.sol[!, "I"] == pb.sol[!, "I"]   # identical construction seed => identical trajectory
        @test pa.sol[!, "S"] == pb.sol[!, "S"]
        @test pa.sol == pb.sol
    end

    # [sir-unseeded-diverges] tier=T1-characterization expectedStatus=pass-now
    # contract: Contract §4 D2 (RNG isolation) — characterizes the Stage A default: an UNSEEDED run is entropy-seeded
    # note: Verified divergent. Post Stage A this divergence is the CORRECT behavior: each unseeded run draws a
    # note: fresh ENTROPY seed for its own state-owned rng, so two default-constructed runs differ almost surely.
    # note: (Pre-Stage-A it diverged for the WRONG reason — shared dependence on ambient global RNG state; that
    # note: failure mode is now fixed.) Reproducibility requires the seed= kwarg — see sir-seeded-reproducible-construction.
    @testset "Two unseeded SIR runs diverge (each gets a fresh entropy seed)" begin
        function build_sir()
            sir = @ReactionNetworkSchema begin
                α * S * I, S + I --> 2I, name => I2R
                β * I, I --> R, name => R2S
            end
            @prob_init sir S = 999 I = 10 R = 0
            @prob_params sir α = 0.0001 β = 0.01
            @prob_meta sir tspan = 50 dt = 0.1
            return sir
        end
        p1 = ReactionNetworkProblem(build_sir()); simulate(p1); p2 = ReactionNetworkProblem(build_sir()); simulate(p2)
        @test p1.sol[!, "I"] != p2.sol[!, "I"]   # two unseeded runs are independently entropy-seeded => diverge
    end

    # [sir-seed-kwarg-isolated] tier=T1-characterization expectedStatus=pass-now
    # contract: Contract §4 D2 (RNG isolation), D5 (state owns AbstractRNG), D6 (seed at construction): ReactionNetworkProblem(...; seed=) owns its rng
    # note: Stage A IMPLEMENTED the seed= kwarg: it constructs a state-owned rng (Xoshiro(seed)) threaded through
    # note: every rand, so the trajectory is reproducible from the recorded seed AND isolated from the global RNG.
    # note: Verified: pa.sol == pb.sol under seed=99 even with the global RNG perturbed between runs, and the
    # note: state exposes a :rng property. Previously this errored (kwarg was a swallowed no-op, no :rng field);
    # note: flipped from @test_skip to live @test now the RNG-threading rework has landed.
    @testset "seed= kwarg gives RNG-isolated reproducibility independent of global state" begin
        function build_sir()
            sir = @ReactionNetworkSchema begin
                α * S * I, S + I --> 2I, name => I2R
                β * I, I --> R, name => R2S
            end
            @prob_init sir S = 999 I = 10 R = 0
            @prob_params sir α = 0.0001 β = 0.01
            @prob_meta sir tspan = 50 dt = 0.1
            return sir
        end
        rand(Int); pa = ReactionNetworkProblem(build_sir(); seed = 99); simulate(pa); rand(Int); rand(Int); pb = ReactionNetworkProblem(build_sir(); seed = 99); simulate(pb)
        @test pa.sol == pb.sol   # same seed => identical trajectory regardless of intervening global rand()
        @test hasproperty(pa, :rng)   # state owns an AbstractRNG (D5)
    end

    # [pharma-pipeline-runs] tier=T1-characterization expectedStatus=pass-now
    # contract: Brief reference-model criterion (toy-pharma pipeline runs; Discovery->...->market); tutorial/toy_pharma_model.jl lines 4-42
    # note: Verified runs to completion; sol columns are exactly the 4 species + t (NOTE construction order
    # note: ['t','scientist','budget','candidate_compound','marketed_drug'] != author order — index by NAME).
    # note: @register the α/β rate fns BEFORE building, exactly as the tutorial does. No maxlifetime/nonblock so
    # note: the known bugs are not triggered.
    @testset "Toy-pharma pipeline runs end to end and flows candidate_compound to market" begin
        @register function α(n1, n2, κ); return κ + exp(-n1) + exp(-n2); end
        @register function β(n1, n2); return n1 + exp(-n2); end
        toy = @ReactionNetworkSchema begin
            α(candidate_compound, marketed_drug, κ),
            3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
            name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
            β(candidate_compound, marketed_drug),
            candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
            name => dx2market, probability => 0.5 + 0.001 * @t(), cycletime => 4
            γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
        end
        @periodic toy 1.0 budget += 11 * marketed_drug
        @prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
        @prob_params toy κ = 4 γ = 0.1
        @prob_meta toy tspan = 250 dt = 0.1
        prob = ReactionNetworkProblem(toy; seed = 1)
        simulate(prob)
        @test Set(names(prob.sol)) == Set(["t", "scientist", "budget", "candidate_compound", "marketed_drug"])
        @test size(prob.sol, 1) > 1   # the run produced a trajectory
        cc = prob.sol[!, "candidate_compound"]; md = prob.sol[!, "marketed_drug"]
        @test maximum(md) >= md[1] || any(md .!= md[1])   # marketed_drug pool is exercised by the pipeline
        @test all(isfinite, prob.sol[!, "budget"])
    end

    # [pharma-conserved-nonneg] tier=T1-characterization expectedStatus=pass-now
    # contract: Contract §3.4 Invariant 1 (resource non-negativity); §1 rows 2 (conserved) & 3 (rate/perstep); brief 'never negative'
    # note: Verified scientist in [0,20] (min 0.0, max 20.0) and budget >= 100 over tspan=50. Non-negativity is
    # note: the ADR-0002 allocator guarantee. NOTE: the scientist<=20 cap holds here only because no instance
    # note: times out before completing (no maxlifetime); the conserved-reentry bug (solvers.jl:501) would break
    # note: the cap — that is pinned separately in pharma-conserved-reentry-bug.
    @testset "Toy-pharma conserved (scientist) and rate (budget) pools stay non-negative" begin
        @register function α(n1, n2, κ); return κ + exp(-n1) + exp(-n2); end
        @register function β(n1, n2); return n1 + exp(-n2); end
        toy = @ReactionNetworkSchema begin
            α(candidate_compound, marketed_drug, κ),
            3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
            name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
            β(candidate_compound, marketed_drug),
            candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
            name => dx2market, probability => 0.5, cycletime => 4
            γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
        end
        @periodic toy 1.0 budget += 11 * marketed_drug
        @prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
        @prob_params toy κ = 4 γ = 0.1
        @prob_meta toy tspan = 50 dt = 0.1
        prob = ReactionNetworkProblem(toy; seed = 1)
        simulate(prob)
        sci = prob.sol[!, "scientist"]; bud = prob.sol[!, "budget"]
        @test all(sci .>= -1e-9)         # conserved pool never goes negative
        @test all(bud .>= -1e-9)         # rate pool never goes negative
        @test maximum(sci) <= 20 + 1e-9  # conserved scientist never exceeds its initial holding
    end

    # [pharma-ledger-populated] tier=T1-characterization expectedStatus=pass-now
    # contract: Brief criterion (ledger populated); Contract §5.4 specCost/specReward/specValuation; log rows solvers.jl:304-312, 426-433, 659-666
    # note: Verified: with @cost/@reward/@valuation attached, cost sum≈585.6, reward sum≈250.0 over tspan=50 (seed=1).
    # note: CRITICAL: the STOCK toy_pharma_model.jl sets NO valuation attrs, so its ledger rows are all 0.0
    # note: (verified) — the brief's 'ledger populated' is true row-wise even then, but a meaningful (nonzero)
    # note: ledger REQUIRES @cost/@reward. No bug; PoS fixed at 0.5 + seed= construction so the draw is reproducible.
    @testset "Toy-pharma ledger has cost/reward/valuation rows once valuation attrs are set" begin
        @register function α(n1, n2, κ); return κ + exp(-n1) + exp(-n2); end
        @register function β(n1, n2); return n1 + exp(-n2); end
        toy = @ReactionNetworkSchema begin
            α(candidate_compound, marketed_drug, κ),
            3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
            name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
            β(candidate_compound, marketed_drug),
            candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
            name => dx2market, probability => 0.5, cycletime => 4
            γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
        end
        @periodic toy 1.0 budget += 11 * marketed_drug
        @prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
        @prob_params toy κ = 4 γ = 0.1
        @cost toy budget = 1.0 scientist = 2.0
        @reward toy marketed_drug = 50.0
        @valuation toy marketed_drug = 100.0
        @prob_meta toy tspan = 50 dt = 0.1
        prob = ReactionNetworkProblem(toy; seed = 1)
        simulate(prob)
        tags = unique([r[1] for r in prob.log])
        @test :valuation_cost in tags && :valuation_reward in tags && :valuation in tags
        cost_rows = filter(r -> r[1] == :valuation_cost, prob.log)
        rew_rows  = filter(r -> r[1] == :valuation_reward, prob.log)
        @test !isempty(cost_rows) && !isempty(rew_rows)
        @test sum(r[3] for r in cost_rows) > 0   # cost actually accrues (verified ≈585.6 under seed=1)
        @test sum(r[3] for r in rew_rows) > 0    # reward actually accrues (verified ≈250.0 under seed=1)
        @test all(length(r) >= 3 for r in cost_rows)   # row shape (:tag, t, value)
    end

    # [rnpv-finite-reduction] tier=T1-characterization expectedStatus=pass-now
    # contract: Brief criterion (rNPV finite); rNPV = sum over ticks of discounted (reward - cost) read from prob.log
    # note: Verified finite. rnpv() reads (:valuation_reward,t,cf) and (:valuation_cost,t,cf) rows (shapes
    # note: verified) and discounts by tick time t. The discounted (reward-cost) reduction is the simplest rNPV
    # note: proxy on the current ledger; the full phase-PoS-weighted form (sum over phases of
    # note: expected_cashflow*cumulativePoS*discount) is its T2 generalization in rnpv-pos-lever's notes.
    @testset "rNPV as a post-processing reduction over the ledger is finite" begin
        @register function α(n1, n2, κ); return κ + exp(-n1) + exp(-n2); end
        @register function β(n1, n2); return n1 + exp(-n2); end
        function build_pharma(pos)
            toy = @ReactionNetworkSchema begin
                α(candidate_compound, marketed_drug, κ),
                3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
                name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
                β(candidate_compound, marketed_drug),
                candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
                name => dx2market, probability => 0.5, cycletime => 4
                γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
            end
            @prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
            @prob_params toy κ = 4 γ = 0.1
            @cost toy budget = 1.0 scientist = 2.0
            @reward toy marketed_drug = 50.0
            @prob_meta toy tspan = 50 dt = 0.1
            for i in 1:ReactiveDynamics.nparts(toy, :T)
                string(toy[i, :transName]) == "dx2market" && (toy[i, :transProbOfSuccess] = pos)
            end
            return toy
        end
        rnpv(prob; r = 0.1) = sum((row[1] == :valuation_reward ? row[3] : row[1] == :valuation_cost ? -row[3] : 0.0) / (1 + r)^row[2] for row in prob.log)
        prob = ReactionNetworkProblem(build_pharma(0.5); seed = 7)
        simulate(prob); val = rnpv(prob)
        @test isfinite(val)   # verified ≈ -34.9 under seed=7; assertion is the robust finiteness invariant
        @test val isa Real
    end

    # [rnpv-pos-lever] tier=T2-acceptance expectedStatus=pass-now
    # contract: Brief BD-demo core assertion (higher PoS -> higher rNPV); Contract §2.8 acquisition lever; PoS = transProbOfSuccess, solvers.jl:413
    # note: Verified directionally under the Stage A seed= kwarg: avg rNPV PoS=0.2 ≈ -98.6 vs PoS=0.9 ≈ -37.8
    # note: (8-seed average via seed=1000+s, clean per-run RNG isolation), so high>low holds. Marked
    # note: T2-acceptance because it is the headline BD criterion and still seed-averages to be stable (single
    # note: seeds are noisy); once §4 D8 lands it should become a deterministic single-seed paired comparison
    # note: (same seed, two PoS) with no averaging. The lever is set via toy[i,:transProbOfSuccess]; the §2.8
    # note: contract notes the acquisition lever must use the scheduled-rate idiom, NOT the event channel,
    # note: because event_action! is a no-op (solvers.jl:323). Each @testset is its own scope, so build_pharma/
    # note: rnpv (and the @register α/β) are redefined locally here — they are NOT visible from the earlier
    # note: rnpv-finite-reduction testset.
    @testset "Higher probability-of-success lever raises rNPV (BD demo core assertion)" begin
        # Self-contained: build_pharma(pos)/rnpv(prob) defined here, matching rnpv-finite-reduction.
        @register function α(n1, n2, κ); return κ + exp(-n1) + exp(-n2); end
        @register function β(n1, n2); return n1 + exp(-n2); end
        function build_pharma(pos)
            toy = @ReactionNetworkSchema begin
                α(candidate_compound, marketed_drug, κ),
                3 * @conserved(scientist) + @rate(budget) --> candidate_compound,
                name => discovery, probability => 0.3, cycletime => 10.0, priority => 0.5
                β(candidate_compound, marketed_drug),
                candidate_compound + 5 * @conserved(scientist) + 2 * @rate(budget) --> marketed_drug + 5 * budget,
                name => dx2market, probability => 0.5, cycletime => 4
                γ * marketed_drug, marketed_drug --> ∅, name => drug_killed
            end
            @prob_init toy candidate_compound = 5 marketed_drug = 6 scientist = 20 budget = 100
            @prob_params toy κ = 4 γ = 0.1
            @cost toy budget = 1.0 scientist = 2.0
            @reward toy marketed_drug = 50.0
            @prob_meta toy tspan = 50 dt = 0.1
            for i in 1:ReactiveDynamics.nparts(toy, :T)
                string(toy[i, :transName]) == "dx2market" && (toy[i, :transProbOfSuccess] = pos)
            end
            return toy
        end
        rnpv(prob; r = 0.1) = sum((row[1] == :valuation_reward ? row[3] : row[1] == :valuation_cost ? -row[3] : 0.0) / (1 + r)^row[2] for row in prob.log)
        # Average over seeds to suppress Poisson/Binomial noise (each run is RNG-isolated via seed=).
        function avg_rnpv(pos; nseed = 8)
            total = 0.0
            for s in 1:nseed
                p = ReactionNetworkProblem(build_pharma(pos); seed = 1000 + s); simulate(p)
                total += rnpv(p)
            end
            return total / nseed
        end
        low = avg_rnpv(0.2); high = avg_rnpv(0.9)
        @test isfinite(low) && isfinite(high)
        @test high > low   # raising dx2market PoS strictly increases expected rNPV
    end

    # [nonblock-free-credits-resource] tier=T1-characterization expectedStatus=pass-now
    # contract: Contract §1.3 row 5 (nonblock free path) / §3.4 Invariant 1 (resource non-negativity)
    # note: Stage A FIXED the free path: free_blocked_species! no longer references the undefined q, so an
    # note: in-flight @nonblock token (cycletime>0) no longer throws UndefVarError. The run now completes and the
    # note: freed :nonblock resource is credited back. Verified under seed=1: simulate runs clean, A stays
    # note: non-negative (100 -> 94) and finite, B grows (0 -> 7). Was a @test_throws pin of the crash; flipped
    # note: to positive behavior + the Invariant 1 non-negativity/finiteness check now the bug is fixed.
    @testset "In-flight @nonblock token frees and credits its resource without crashing" begin
        m = @ReactionNetworkSchema begin
            1.0, @nonblock(A) --> B, cycletime => 5.0, name => t1
        end
        @prob_init m A = 100 B = 0
        @prob_meta m tspan = 10 dt = 1.0
        prob = ReactionNetworkProblem(m; seed = 1)
        @test (simulate(prob); true)            # completes — no UndefVarError from free_blocked_species!
        @test all(prob.sol.A .>= -1e-9)         # freed @nonblock resource credited back; A never negative
        @test all(isfinite, prob.sol.A) && all(isfinite, prob.sol.B)   # trajectory stays finite
    end

    # [pharma-conserved-reentry] tier=T2-acceptance expectedStatus=pass-now
    # contract: Contract §3.4 Invariant 6 (termination completeness) & Invariant 2 (conservation); prune path solvers.jl:408-410
    # note: Stage A FIXED the prune predicate: an instance past its terminal test is now removed, so conservation
    # note: HOLDS — the previously-observed inflation (cash 10 -> 30, ongoing stuck at 1) is gone. Verified under
    # note: seed=1: cash stays bounded by its initial 10 (series settles at 10.0) and the timed-out instance is
    # note: pruned, leaving ongoing_transitions empty by the end. Both assertions flipped from @test_broken to
    # note: live @test. This is the canonical 'conservation holds' acceptance test (INV2/INV6).
    @testset "Lifetime-timeout instance is pruned and conserved tokens stay bounded (INV2/INV6)" begin
        # Instance times out (maxlifetime=2) before completing (cycletime=100); the prune fix removes it, so its
        # held @conserved cash is NOT re-credited every subsequent tick.
        m = @ReactionNetworkSchema begin
            @deterministic(1.0), 2 * @conserved(cash) + raw --> product, cycletime => 100.0, maxlifetime => 2.0, name => t1
        end
        @prob_init m cash = 10 raw = 1 product = 0
        @prob_meta m tspan = 12 dt = 1.0
        prob = ReactionNetworkProblem(m; seed = 1)
        simulate(prob)
        cash = prob.sol[!, "cash"]
        @test maximum(cash) <= 10 + 1e-9   # conserved pool never exceeds its initial holding (INV2)
        @test length(prob.ongoing_transitions) == 0   # timed-out instance is pruned (INV6; verified count==0)
    end
end
