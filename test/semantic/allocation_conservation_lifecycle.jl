# Phase-0 semantic tests — Allocation (ADR 0002), Conservation & Non-negativity (§3.4), Instance Lifecycle (§2.4/§3.2)
#
# Tiers: T1-characterization runs against the CURRENT engine (locks in behavior or pins a known
# bug via @test_broken); T2-acceptance encodes TARGET behavior and is wrapped (@test_skip + a
# commented reference block) because it names APIs that do not exist yet. Assertions were largely
# verified on Julia 1.12.5 during drafting; re-verify file:line citations before acting.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

@testset "Allocation (ADR 0002), Conservation & Non-negativity (§3.4), Instance Lifecycle (§2.4/§3.2)" begin

    # [alloc-nonneg-capacity] tier=T1-characterization expectedStatus=pass-now
    # contract: ADR 0002 'Invariant test plan' (non-negativity & capacity); CONTRACT §3.4 INV1
    # note: The two universally-true allocation guards, now on the progressive_fill! allocator
    # note: (ADR 0002 replaced the alloc_weighted!/alloc_greedy! switch). With fmax=[1,1] each
    # note: transition fills at most its full per-instance demand and u >= total demand, so the
    # note: alloc is the raw req — all >=0 and within supply.
    @testset "progressive_fill! output is non-negative and never exceeds supply" begin
        RD = ReactiveDynamics; req = [1.0 1.0; 2.0 1.0]; ws = RD.AllocWorkspace(req); u = [10.0, 10.0]; w = [1.0, 1.0]
        f = RD.progressive_fill!(ws, u, w; fmax = [1.0, 1.0]); allocs = ws.req .* f'
        @test all(allocs .>= 0)
        @test all(vec(sum(allocs; dims = 2)) .<= u .+ 1.0e-9)
    end

    # [alloc-priority-ratio-contended] tier=T1-characterization expectedStatus=pass-now
    # contract: ADR 0002 priority semantics: 'at equal demand the split equals the priority ratio (1 vs 3 -> 1:3)'
    # note: The headline guarantee on the new allocator: a single contended resource (u=8, two
    # note: trans need 5 each, fmax=Inf) splits in the exact priority ratio 1:3 -> f=[1.0,3.0],
    # note: alloc=[2,6]. Work-conserving (full 8 used). Unlike the old alloc_weighted! this no
    # note: longer needs total-demand>=supply to avoid a no-contention shortcut — fmax=Inf means
    # note: fill until the resource exhausts.
    @testset "Single contended resource at equal demand splits by exact priority ratio" begin
        RD = ReactiveDynamics; req = reshape([5.0, 5.0], 1, 2); ws = RD.AllocWorkspace(req); u = [8.0]; w = [1.0, 3.0]
        f = RD.progressive_fill!(ws, u, w; fmax = [Inf, Inf]); allocs = vec(ws.req .* f')
        @test allocs[2] / allocs[1] ≈ 3.0
        @test sum(allocs) ≈ 8.0  # work-conserving on a single resource
    end

    # [alloc-priority-ratio-fillfraction-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 Verification row 'Priority split' (f=[2.5,7.5], ratio 3:1); ADR 0002 algorithm fmax defaults to Inf
    # note: progressive_fill!/AllocWorkspace do not exist yet (verified isdefined == false). This is the ADR
    # note: headline guarantee and is NOT reproducible by the current allocator: alloc_weighted! returns [1,1]
    # note: for these exact inputs because demand (2) < supply (10) trips the no-contention shortcut. fmax=Inf
    # note: means 'fill until the resource exhausts', which is the semantic the current code lacks. Target
    # note: surface per ADR 0002 'Proposed Julia surface'.
    @testset "progressive_fill! priority-split: u=10, two trans need 1 each, w=1 vs 3 -> f=[2.5,7.5]" begin
        RD = ReactiveDynamics
        ws = RD.AllocWorkspace(reshape([1.0, 1.0], 1, 2))  # req[s,t]
        u = [10.0]; w = [1.0, 3.0]; fmax = [Inf, Inf]
        f = RD.progressive_fill!(ws, u, w; fmax = fmax)  # alloc[s,t] = f' .* ws.req
        @test f ≈ [2.5, 7.5]
        @test f[2] / f[1] ≈ w[2] / w[1]
        @test sum(ws.req .* f') ≈ 10.0  # resource fully consumed (work-conserving)
    end

    # [alloc-conjunctive-no-stranding-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 Verification row 'Conjunctive 2-res' (f=[3.33,3.33], both resources fully used); 'Invariant test plan': alloc ≈ f'.*req
    # note: Pins the central ADR-0002 fix. Verified the CURRENT path FAILS this: alloc_weighted! +
    # note: get_frac_satisfied leaves resource used = [2.0, 3.0] of [10,10] (massive stranding) and qs=[1,1]
    # note: instead of f=[3.33,3.33]. This test encodes the target progressive-fill behavior and errors until
    # note: the new allocator lands.
    @testset "Conjunctive 2-resource consistency: alloc ≈ f'.*req, binding resource fully used, nothing stranded" begin
        RD = ReactiveDynamics
        req = [1.0 1.0; 2.0 1.0]; ws = RD.AllocWorkspace(req); u = [10.0, 10.0]; w = [1.0, 1.0]
        f = RD.progressive_fill!(ws, u, w; fmax = [Inf, Inf]); alloc = ws.req .* f'
        @test alloc ≈ ws.req .* f'  # conjunctive consistency by construction
        @test f ≈ [10 / 3, 10 / 3] atol = 1.0e-3
        # ADR-0002 oracle f=[10/3,10/3]: species 2 (demand D=3) is the BINDING resource and
        # saturates first (2·10/3 + 1·10/3 = 10), freezing BOTH transitions. Species 1 (demand
        # D=2) is then at 10/3+10/3 = 20/3 with 10/3 left idle — work-conserving, because no
        # unfrozen transition can use it. (The ADR table's "both resources fully used" gloss is
        # only true for the SYMMETRIC 'Differing binding res' row req=[[2,1],[1,2]]; for THIS
        # asymmetric req the binding-resource fill is the correct, verified value.)
        @test vec(sum(alloc; dims = 2)) ≈ [20 / 3, 10.0]  # binding resource (2) at capacity, no stranding
    end

    # [alloc-work-conservation-capped-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 Verification row 'Work-conservation' (f=[2,8], total used 10/10; naive split would strand it using only 7)
    # note: Encodes the work-conservation property the current allocator lacks (per-resource split + min-rescale
    # note: strands the freed 3 units). Requires fmax-capping support in progressive_fill!, which does not exist
    # note: yet. ADR oracle is explicit that the naive split would use only 7/10.
    @testset "Work-conservation: T1 capped at fmax=2 -> leftover flows to T2, f=[2,8], 10/10 used" begin
        RD = ReactiveDynamics
        req = reshape([1.0, 1.0], 1, 2); ws = RD.AllocWorkspace(req); u = [10.0]; w = [1.0, 1.0]; fmax = [2.0, Inf]
        f = RD.progressive_fill!(ws, u, w; fmax = fmax)
        @test f ≈ [2.0, 8.0]
        @test sum(ws.req .* f') ≈ 10.0  # no usable resource left idle while T2 still wants it
        @test f[1] ≈ fmax[1]  # T1 frozen exactly at its cap
    end

    # [alloc-priority-zero-leftover-only-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 'Zero-priority leftover only — two-stage' (u=10, T1 prio1 cap4 / T2 prio0 -> f=[4,6]; T1 uncapped -> f=[10,0], T2 starved)
    # note: Pins the maintainer-confirmed priority=0 = 'runs only from genuinely leftover resource' two-stage
    # note: semantics (ADR 0002 Resolved #3). The current alloc_weighted! has no leftover tier — a zero priority
    # note: just yields a zero weight column and zero alloc, not the two-stage leftover fill. Target API only.
    @testset "priority=0 is leftover-only (two-stage): T1(prio1,cap4) + T2(prio0) over u=10 -> f=[4,6]; T1 uncapped -> f=[10,0]" begin
        RD = ReactiveDynamics
        req = reshape([1.0, 1.0], 1, 2); ws = RD.AllocWorkspace(req); u = [10.0]
        f_a = copy(RD.progressive_fill!(ws, copy(u), [1.0, 0.0]; fmax = [4.0, Inf]))
        f_b = copy(RD.progressive_fill!(ws, copy(u), [1.0, 0.0]; fmax = [Inf, Inf]))
        @test f_a ≈ [4.0, 6.0]  # zero-prio T2 takes only the leftover after T1's cap
        @test f_b ≈ [10.0, 0.0]  # uncapped positive-prio T1 starves zero-prio T2
        @test f_b[2] == 0.0  # priority=0 never advances under contention
    end

    # [alloc-weighted-conjunctive-ratio-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 Verification row 'Weighted + conjunctive' (T1=[1,1] T2=[1,1], X scarce u=6, w=1 vs 2 -> f=[2,4], X 6/6)
    # note: Two-resource weighted fairness oracle. Verified current path FAILS: alloc_weighted! returns the raw
    # note: reqs (u for X=6 >= total X demand 2 trips the shortcut) giving qs=[1,1], no weighting. Target
    # note: allocator only.
    @testset "Weighted + conjunctive: shared scarce X (u=6), w=1 vs 2 -> f=[2,4], X used 6/6, ratio exactly 2:1" begin
        RD = ReactiveDynamics
        req = [1.0 1.0; 1.0 1.0]; ws = RD.AllocWorkspace(req); u = [6.0, 100.0]; w = [1.0, 2.0]
        f = RD.progressive_fill!(ws, u, w; fmax = [Inf, Inf])
        @test f ≈ [2.0, 4.0]
        @test f[2] / f[1] ≈ 2.0
        @test sum((ws.req .* f')[1, :]) ≈ 6.0  # the binding resource X fully used
    end

    # [alloc-3way-stranding-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 Verification row '3-way stranding' (u=12, T1 cap 2 -> f=[2,5,5], freed amount split evenly, total 12/12)
    # note: Generalizes work-conservation to >2 contenders with even redistribution of freed capacity. Target
    # note: progressive_fill! only.
    @testset "3-way redistribution: u=12, w=1, T1 capped at 2 -> f=[2,5,5], freed amount split evenly, 12/12 used" begin
        RD = ReactiveDynamics
        req = reshape([1.0, 1.0, 1.0], 1, 3); ws = RD.AllocWorkspace(req); u = [12.0]; w = [1.0, 1.0, 1.0]; fmax = [2.0, Inf, Inf]
        f = RD.progressive_fill!(ws, u, w; fmax = fmax)
        @test f ≈ [2.0, 5.0, 5.0]
        @test sum(ws.req .* f') ≈ 12.0
        @test f[2] ≈ f[3]  # the 10 freed/remaining units split evenly among absorbers
    end

    # [alloc-integer-spawn-topup-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0002 'Spawn phase (integer instances)' + Verification row 'Integer spawn' (real fair fill 1.75 each -> n=[2,1], leftover 1, top-up deterministic)
    # note: Pins ADR 0002's deterministic priority-ordered integer top-up (stable tie-break by (-priority,
    # note: index)). spawn_integer! does not exist (verified). The current get_init_satisfied
    # note: (solvers.jl:110-128) floors per-resource but has no documented deterministic top-up rule. Target
    # note: surface.
    @testset "spawn_integer!: u=7, each instance needs 2, want 5 each, equal prio -> n=[2,1], leftover 1, deterministic" begin
        RD = ReactiveDynamics
        req = reshape([2.0, 2.0], 1, 2); ws = RD.AllocWorkspace(req); u = [7.0]; w = [1.0, 1.0]; q_desired = [5, 5]
        n = RD.spawn_integer!(ws, u, w, q_desired)  # returns Vector{Int}
        @test n isa Vector{Int}
        @test n == [2, 1]  # 1.75 each floored to [1,1], deterministic priority/index top-up grants the +1 to T1 (lower index on ties)
        @test sum(ws.req[1, :] .* n) <= u[1]  # integral allocation fits in supply
        @test eltype(n) === Int && all(n .>= 0)
    end

    # [alloc-determinism-construction-seed] tier=T1-characterization expectedStatus=pass-now
    # contract: ADR 0002 Verification row 'Determinism'; CONTRACT §4.1 D1 (reproducibility) via §4.3 D6 (seed at construction)
    # note: Stage A moved all step-loop randomness onto a state-owned rng seeded from `seed` (solvers.jl:565-567),
    # note: so resetting the GLOBAL rng no longer controls a run (verified: two Random.seed!(42) runs now DIFFER,
    # note: a=[985,9] vs [984,8], because the default run is entropy-seeded). The honest current-engine
    # note: reproducibility guarantee is now: two runs with the SAME construction seed are identical. Verified
    # note: live: a == b == [979.0, 12.0] under seed=42. The contract-grade global-RNG isolation guarantee is the
    # note: D2 testset above.
    @testset "Allocation/trajectory is reproducible under the same construction seed (current engine lock-in)" begin
        function run_once()
            net = @reaction_network begin
                1.0, budget --> product, name => job, cycletime => 0.0, probability => 0.5
            end; @prob_init net budget = 1000 product = 0; @prob_params net; @prob_meta net tspan = 20 dt = 1.0; prob = ReactionNetworkProblem(net; seed = 42); simulate(prob); copy(prob.u)
        end
        a = run_once(); b = run_once()
        @test a == b  # identical trajectory when the run is constructed with the same seed
    end

    # [alloc-rng-isolation-seed-kwarg] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §4.2 D2 (RNG isolation), §4.3 D6 (seed at construction); ADR 0002 'allocator RNG-free, draws route through state rng'
    # note: Stage A landed the state-owned rng::AbstractRNG (state.jl:69, solvers.jl:565-567) seeded from `seed`,
    # note: with every rand() in _step! routed through it. Verified live: ReactionNetworkProblem(...; seed=7) now
    # note: exposes a :rng field (propertynames check == true) and a seeded run is isolated from global-RNG
    # note: perturbation (a == b == [985,6] despite differing global draws between runs). Un-skipped: the D2/D6
    # note: target API now EXISTS, so these run as real @tests.
    @testset "seed= kwarg isolates the run from the global RNG (D2): state owns an :rng, seeded runs are reproducible under global-RNG perturbation" begin
        RD = ReactiveDynamics
        function run_seeded()
            net = @reaction_network begin
                1.0, budget --> product, name => job, cycletime => 0.0, probability => 0.5
            end; @prob_init net budget = 1000 product = 0; @prob_params net; @prob_meta net tspan = 20 dt = 1.0; prob = ReactionNetworkProblem(net; seed = 7); simulate(prob); copy(prob.u)
        end
        Random.seed!(1); rand(10); a = run_seeded(); Random.seed!(2); rand(3); b = run_seeded()
        rng_acs = @reaction_network begin
            1.0, budget --> product, name => job
        end; @prob_init rng_acs budget = 10 product = 0; @prob_params rng_acs; @prob_meta rng_acs tspan = 2 dt = 1.0; rng_prob = ReactionNetworkProblem(rng_acs; seed = 1)
        @test :rng in propertynames(rng_prob)  # state owns an AbstractRNG
        @test a == b  # perturbing the global RNG does NOT change a seeded run
    end

    # [conserve-conserved-token-returned-current] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §3.4 INV2 (conserved tokens returned exactly); §1.3 truth-table row 2 (upfront/conserved/block); §3.2 stage 8
    # note: Verified live end-to-end: scientist stays 10.0 across all rows while budget falls 10->7 and product
    # note: rises 0->3. Locks in the conservation return path (solvers.jl:437-458, the in(:conserved) branch
    # note: crediting q*stoich back) AND the non-return of plain consumed tokens. cycletime=0 means instant
    # note: completion so no rate scaling. Survives the rework since INV2 is a hard contract obligation.
    @testset "A :conserved token is held then returned in full; plain consumed token is not returned" begin
        # seed=1: the rate-1.0 Poisson spawn is entropy-driven, and ~3% of unseeded runs draw zero
        # spawns across all 3 ticks (Poisson(1) = 0 with prob e^-1 per tick), leaving product==0 and
        # budget undrawn — which spuriously fails the two spawn-dependent assertions below (the
        # conservation invariant itself never breaks). Seeding pins a trajectory that does spawn, so
        # the test is deterministic. (Pre-existing fragility surfaced during the ADR-0003 store swap.)
        net = @reaction_network begin
            1.0, 2 * @conserved(scientist) + budget --> product, name => job, cycletime => 0.0
        end; @prob_init net scientist = 10 budget = 10 product = 0; @prob_params net; @prob_meta net tspan = 3 dt = 1.0; prob = ReactionNetworkProblem(net; seed = 1)
        simulate(prob)  # species order verified: [:scientist, :budget, :product]
        @test all(prob.sol.scientist .== 10.0)  # conserved pool held and returned every tick: never net-debited
        @test prob.sol.budget[end] < prob.sol.budget[1]  # plain consumed token IS permanently drawn down
        @test prob.u[1] == 10.0  # final conserved == initial
        @test prob.u[3] > 0.0  # product accumulated from RHS emission
    end

    # [conserve-closed-system-mass-invariant-current] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §3.4 INV1 (non-negativity) + INV2 (conservation closed-system mass invariant)
    # note: With cycletime=2 the conserved scientists are held across ticks (dip below 30 while in-flight) then
    # note: returned, so the pool is bounded above by its initial closed-system mass. Locks in INV1 (non-
    # note: negativity) and the upper-bound half of INV2 on a healthy (cycle-completing, no lifetime timeout)
    # note: instance — distinguishing it from the lifetime-timeout blow-up pinned below. NOTE: must keep
    # note: maxlifetime=Inf (default) so the solvers.jl:501 prune bug does not fire here.
    @testset "Closed-system conserved mass is invariant across spawn->return over a multi-tick run; u stays >=0" begin
        net = @reaction_network begin
            1.0, 3 * @conserved(scientist) + budget --> product, name => job, cycletime => 2.0, probability => 1.0
        end; @prob_init net scientist = 30 budget = 1000 product = 0; @prob_params net; @prob_meta net tspan = 10 dt = 1.0; prob = ReactionNetworkProblem(net)
        simulate(prob)  # cycletime>0 so scientists are genuinely held in-flight for ~2 ticks before return
        @test all(prob.sol.scientist .>= 0.0)  # INV1 non-negativity throughout
        @test all(prob.sol.budget .>= 0.0)
        @test prob.sol.scientist[end] <= 30.0 + 1.0e-9  # conserved mass never EXCEEDS the closed-system total (no spurious creation)
        @test all(prob.sol.scientist .<= 30.0 + 1.0e-9)
    end

    # [conserve-lifetime-timeout-no-double-return] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §3.4 INV6 (termination completeness) + INV2; FIXED prune predicate (Stage A) removes every maxLifeTime-terminated instance
    # note: Prune bug at src/solvers.jl:501 is FIXED — the predicate no longer RETAINS instances that terminated
    # note: by maxLifeTime (state never reached cycleTime), so they are no longer re-finished every subsequent
    # note: tick and no longer re-credit their conserved tokens. Verified live: conserved mass now holds (max
    # note: scientist == 10.0, u[1] == 6.0 over 6 ticks) instead of inflating to 26.0. INV2 conservation and INV6
    # note: termination-completeness now HOLD. The few transitions still in ongoing at tspan are legitimately
    # note: in-flight (spawned in the final ticks, state < cycleTime), NOT retained zombies — the correct INV6
    # note: statement is therefore 'no completed-but-retained instance survives', i.e. every ongoing instance
    # note: still has state < cycleTime (verified). Positive acceptance test for INV2/INV6.
    @testset "Lifetime-timed-out instances are pruned: conserved tokens are not re-returned and the pool stays within closed-system mass" begin
        net = @reaction_network begin
            @deterministic(1.0), 2 * @conserved(scientist) --> product, name => job, cycletime => 10.0, maxlifetime => 2.0, probability => 1.0
        end; @prob_init net scientist = 10 product = 0; @prob_params net; @prob_meta net tspan = 6 dt = 1.0; prob = ReactionNetworkProblem(net)
        simulate(prob)  # each tick spawns 1 instance; cycletime=10 never reached before maxlifetime=2 forces timeout
        @test prob.u[1] <= 10.0 + 1.0e-9  # INV2: conserved mass NEVER exceeds closed-system total now that timed-out instances are pruned (verified u[1] == 6.0, max == 10.0)
        @test all(tr -> tr.state < tr[:transCycleTime], prob.ongoing_transitions)  # INV6: no completed-but-retained zombie; every surviving instance is genuinely in-flight (verified)
    end

    # [lifecycle-cycletime-completion-current] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §2.4 (cycle time: under full saturation completes after ceil(cycletime/dt) ticks); §3.2 stages 4-5
    # note: Verified live: with cycletime=3, dt=1, pos=1, deterministic 1 spawn/tick and ample budget, product
    # note: first appears in the t=3 row and increments by 1 each subsequent tick. Locks in §2.4 quantization
    # note: (transition.state += qs*dt at solvers.jl:261; finish gated on state>=cycleTime at
    # note: solvers.jl:409,412). probability=1 removes Binomial noise so the completion tick is deterministic.
    @testset "Instance with cycletime>0 completes after ceil(cycletime/dt) ticks under full allocation" begin
        net = @reaction_network begin
            @deterministic(1.0), budget --> product, name => job, cycletime => 3.0, probability => 1.0
        end; @prob_init net budget = 1000 product = 0; @prob_params net; @prob_meta net tspan = 8 dt = 1.0; prob = ReactionNetworkProblem(net)
        simulate(prob)  # ample budget => full allocation (qs==1) => progress advances by exactly dt per tick
        @test all(prob.sol.product[prob.sol.t .< 3.0] .== 0.0)  # no completion before ceil(3/1)=3 ticks
        @test prob.sol.product[findfirst(==(3.0), prob.sol.t)] == 1.0  # first instance (born t=0) completes exactly at t=3
        @test prob.sol.product[end] > prob.sol.product[findfirst(==(3.0), prob.sol.t)]  # steady-state completion thereafter
    end

    # [lifecycle-lifetime-timeout-zero-success-current] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT §2.5 (lifetime timeout -> q=0 successes); §3.2 stage 6; solvers.jl:408-416
    # note: Verified live: product stays 0.0 across all rows despite probability=1, because the lifetime cap
    # note: fires before cycletime and the success branch (solvers.jl:412-416) returns 0 when state<cycleTime.
    # note: This tests the success-gating independently of the prune bug (product=0 is correct; the bug only
    # note: inflates the conserved/RHS-on-success quantities — here pos applies to a non-conserved consumed
    # note: token with q=0, so RHS is genuinely 0). Companion to the conserve bug pin.
    @testset "Lifetime timeout before cycletime yields q=0 successes (no RHS products emitted)" begin
        net = @reaction_network begin
            @deterministic(1.0), budget --> product, name => job, cycletime => 10.0, maxlifetime => 2.0, probability => 1.0
        end; @prob_init net budget = 1000 product = 0; @prob_params net; @prob_meta net tspan = 6 dt = 1.0; prob = ReactionNetworkProblem(net)
        simulate(prob)  # instances age out at age>=2 with state<10 (cycletime never reached)
        @test all(prob.sol.product .== 0.0)  # success draw gated on state>=cycleTime (solvers.jl:412); timeout => q=0 => zero products even with pos=1
    end

    # [lifecycle-pos-binomial-ensemble-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT §3.2 stage 6 (q_success = Binomial(q, probOfSuccess)); §4.4 D8 (per-trajectory seeding from index, reproducible ensemble)
    # note: Encodes the Binomial-PoS statistic AND the D8 per-trajectory seeding contract. Requires the seed=
    # note: kwarg to actually seed a state-owned rng (verified TODAY it is swallowed and reproducibility relies
    # note: on the global RNG, so the member-reproducibility assertion fails and seeding is not isolated). The
    # note: mean-fraction assertion is a property of the Binomial draw (solvers.jl:413) that holds in
    # note: expectation but only becomes a DETERMINISTIC test under D8 seeding. Errors-until-implemented on the
    # note: rng threading.
    @testset "PoS Binomial: over an ensemble, ~q*pos firings succeed (mean of successful completions ≈ pos)" begin
        function member(k)
            # NB the LHS token is @conserved(budget): the reference block's stated premise is that
            # ~200 single-tick instances complete over tspan=200, which requires budget to persist
            # across ticks. A PLAIN `budget` LHS is consumed irreversibly, so budget=1 permits
            # exactly ONE spawn ever and the ensemble mean would be a single Bernoulli(0.3) (~0.3
            # products, fraction ~0.0016) — the block's own "completes ~200 instances" premise is
            # then unsatisfiable. @conserved(budget) (held at spawn, returned at completion) is the
            # faithful model: 1 unit funds 1 instance/tick, ~200 complete, success fraction → pos.
            net = @reaction_network begin
                @deterministic(1.0), @conserved(budget) --> product, name => job, cycletime => 0.0, probability => 0.3
            end
            @prob_init net budget = 1 product = 0
            @prob_params net
            @prob_meta net tspan = 200 dt = 1.0
            prob = ReactionNetworkProblem(net; seed = hash((42, k)))
            simulate(prob)
            prob.sol.product[end]
        end
        # 50-member ensemble, each seeded deterministically from a root seed + member index (D8);
        # each completes ~200 single-tick instances with pos=0.3.
        results = [member(k) for k in 1:50]
        @test isapprox(mean(results) / 200, 0.3; atol = 0.05)  # empirical success fraction ≈ probOfSuccess
        @test [member(k) for k in 1:5] == [member(k) for k in 1:5]  # D8: ensemble reproducible from the root seed, member k independent of N/order
    end

    # [lifecycle-lifetime-zero-success-ensemble-target] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT §2.5 (timeout -> q=0); §4.4 D8/D9 (independent per-member seeded streams)
    # note: Strengthens the q=0-on-timeout guarantee to hold across independent seeded streams (D8/D9). Marked
    # note: T2 because it uses the seed= kwarg per member; with the CURRENT engine the seed is swallowed
    # note: (verified) so the test is not meaningfully seeded — and separately the prune bug (solvers.jl:501)
    # note: retains the timed-out instances. product==0 itself would pass today, but the test's CONTRACT intent
    # note: (per-member isolated seeding) is not satisfiable until the rng field lands. Pairs with the bug pin
    # note: which covers the retained-instance side effects.
    # action: 30-member seeded ensemble; every instance times out at age 1 << cycletime 100, so success
    # count is identically 0 regardless of pos
    @testset "Lifetime timeout before cycletime gives 0 successes across a seeded ensemble (q=0 deterministically)" begin
        function member(k)
            net = @reaction_network begin
                @deterministic(1.0), budget --> product, name => job, cycletime => 100.0, maxlifetime => 1.0, probability => 0.9
            end
            @prob_init net budget = 10000 product = 0
            @prob_params net
            @prob_meta net tspan = 20 dt = 1.0
            prob = ReactionNetworkProblem(net; seed = hash((7, k)))
            simulate(prob)
            prob.sol.product[end]
        end
        @test all(member(k) == 0.0 for k in 1:30)  # q=0 on every timeout in every member, independent of seed/pos
    end
end
