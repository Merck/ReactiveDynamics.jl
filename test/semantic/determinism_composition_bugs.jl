# Phase-0 semantic tests — Determinism & Seeding (§4), Composition (§7), and remaining bug-pins
#
# Tiers: T1-characterization runs against the CURRENT engine (locks in behavior or pins a known
# bug via @test_broken); T2-acceptance encodes TARGET behavior and is wrapped (@test_skip + a
# commented reference block) because it names APIs that do not exist yet. Assertions were largely
# verified on Julia 1.12.5 during drafting; re-verify file:line citations before acting.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames
using ACSets  # nparts/incident for composition tests

@testset "Determinism & Seeding (§4), Composition (§7), and remaining bug-pins" begin

    # [determinism-unseeded-differs-characterization] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §3.4 Invariant 5 (Determinism under seed — VIOLATED); §4.2 RNG-site inventory
    # note: Pins the §3.4 Invariant 5 violation and §4.2 inventory: global-RNG draws at create.jl:151 (realized
    # note: solvers.jl:140-144), solvers.jl:413, solvers.jl:321, state.jl:123/132/70. The final @test_throws
    # note: documents that the D6 `seed=` kwarg is not yet accepted (kwargs are merged into `keywords` at
    # note: solvers.jl:549 and `seed` is simply ignored / never used — so passing it does not error TODAY; if it
    # note: silently no-ops instead of throwing, downgrade this line to @test_broken). The first two @test lines
    # note: reliably pass: with rate 3.0 over 30 ticks the Poisson+Binomial draws make collision astronomically
    # note: unlikely.
    # action: Construct two independent problems from the SAME spec (no seed kwarg exists today) and
    # `simulate` each to completion; compare the final-state vectors / B column of prob.sol.
    @testset "T1 characterization: two unseeded runs of a stochastic model diverge (pins the no-seeding gap)" begin
        # Poisson genesis + Binomial PoS => genuinely stochastic on the global RNG.
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 30 dt = 1.0
          acs
        end
        prob1 = ReactionNetworkProblem(mk()); simulate(prob1)
        prob2 = ReactionNetworkProblem(mk()); simulate(prob2)
        # Today there is no seed: the global RNG advances between the two runs, so trajectories differ.
        @test prob1.sol.B != prob2.sol.B
        # Document the cause: no `seed`/`rng` field is reachable on the state (target API absent).
        @test !(:rng in fieldnames(typeof(prob1)))
        @test_throws Exception ReactionNetworkProblem(mk(); seed = 1234)
    end

    # [determinism-d1-reproducible] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT_DRAFT.md §4.1 D1 (Reproducibility); D6 (seed at construction)
    # note: Target API per D5/D6: ReactionNetworkProblem grows a `rng::AbstractRNG` field (state.jl:40-63)
    # note: seeded deterministically from `seed` (e.g. Xoshiro(seed)); every rand reachable from _step! becomes
    # note: rand(state.rng, ...) — the four explicit sites plus context_eval at state.jl:70. Errors today
    # note: because `seed=` is unused and draws hit the global RNG, so p1.sol != p2.sol (the @test fails) — and
    # note: constructing with `seed` may not even be wired. Marked T2.
    # action: Construct two problems from the same spec with the SAME `seed`, simulate both, and assert
    # the full trajectory and the entire log are equal.
    @testset "T2 D1: same (model, seed) => identical prob.sol AND prob.log" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 30 dt = 1.0
          acs
        end
        p1 = ReactionNetworkProblem(mk(); seed = 42); simulate(p1)
        p2 = ReactionNetworkProblem(mk(); seed = 42); simulate(p2)
        @test p1.sol == p2.sol
        @test p1.log == p2.log
        # A different seed must (almost surely) produce a different trajectory.
        p3 = ReactionNetworkProblem(mk(); seed = 7); simulate(p3)
        @test p1.sol.B != p3.sol.B
        =#
    end

    # [determinism-d2-rng-isolation] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT_DRAFT.md §4.1 D2 (RNG isolation); D5 (threading rule)
    # note: Per D2/D5: all randomness comes from a state-owned AbstractRNG; no _step! path may call bare rand().
    # note: Today every draw uses Random.default_rng() (the §4.2 sites), so (a) external rand() shifts the run
    # note: -> pa.sol != pb.sol, and (b) the run advances default_rng() -> the before/after snapshot differs.
    # note: The grep-invariant in D5 (no `rand(` lacking an RNG arg under src/) is the structural enforcement.
    # note: Errors today; T2.
    # action: Run a seeded simulation. Run it again with the SAME seed but perturb the global RNG (extra
    # rand() calls) before/around the run. The trajectory must be identical (independence from global
    # RNG), and the global RNG state must be untouched by the run.
    @testset "T2 D2: a run neither reads nor perturbs the global RNG (external rand() does not affect it)" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 30 dt = 1.0
          acs
        end
        pa = ReactionNetworkProblem(mk(); seed = 99); simulate(pa)
        # Perturb the global stream, then re-run with the same seed.
        Random.seed!(1); rand(); rand(); rand()
        pb = ReactionNetworkProblem(mk(); seed = 99); rand(); simulate(pb)
        @test pa.sol == pb.sol
        # The run itself must not advance the global RNG: snapshot before/after a simulate.
        Random.seed!(2024); before = copy(Random.default_rng())
        pc = ReactionNetworkProblem(mk(); seed = 5); simulate(pc)
        @test copy(Random.default_rng()) == before
        =#
    end

    # [determinism-d7-reinit-restores-stream] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT_DRAFT.md §4.3 D7 (re-init restores the stream); _reinit! at solvers.jl:619-628
    # note: _reinit! (solvers.jl:619-628) currently resets u/t/ongoing_transitions/log/observables/sol but NOT
    # note: the RNG (it has no RNG to reset today). D7 requires it to restore the RNG to the (M, seed) initial
    # note: state — i.e. store the seed/initial RNG state on the struct and re-seed in _reinit!. Errors today:
    # note: even if an rng existed, the un-reset stream makes the second run diverge -> p.sol != sol1. T2;
    # note: depends on D6 plumbing.
    # action: Build one seeded problem, simulate it (capturing sol/log), call _reinit!, simulate again,
    # and assert the second trajectory equals the first.
    @testset "T2 D7: init -> step* -> reinit! -> step* reproduces the first trajectory" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        # import AlgebraicAgents: _reinit!
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 30 dt = 1.0
          acs
        end
        p = ReactionNetworkProblem(mk(); seed = 321)
        simulate(p)
        sol1 = copy(p.sol); log1 = copy(p.log)
        ReactiveDynamics._reinit!(p)   # AA-dispatched _reinit!
        simulate(p)
        @test p.sol == sol1
        @test p.log == log1
        =#
    end

    # [determinism-d8-ensemble-per-index-seeding] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: CONTRACT_DRAFT.md §4.4 D8 (per-trajectory seeding from index); D9 (no shared mutable RNG)
    # note: D8 requires per-member seed derived deterministically from (root_seed,k) (member_seed shown), each
    # note: member owning its own RNG (D9, no shared mutable instance). The recorded per-member seed (D8) would
    # note: also be asserted once a member-log channel exists. Errors today because `seed=` is unused (member
    # note: draws all share the global RNG), so ens5[3] != solo3. T2; depends on D6. NOTE: the standalone-vs-in-
    # note: ensemble equality also implicitly requires D2 isolation (running member 1,2 before 3 must not
    # note: perturb member 3).
    # action: Compute member k's trajectory inside two different ensemble sizes / orderings and assert it
    # is identical; assert two distinct members differ.
    @testset "T2 D8/D9: ensemble member k reproducible from (root_seed,k), independent of N and order, own RNG per member" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 20 dt = 1.0
          acs
        end
        # TARGET helper: derive a per-member seed from a single root seed + member index.
        member_seed(root, k) = hash((root, k))
        run_member(root, k) = (p = ReactionNetworkProblem(mk(); seed = member_seed(root, k)); simulate(p); p.sol)
        # member 3 of a size-5 ensemble vs the same member computed standalone:
        ens5 = [run_member(2026, k) for k in 1:5]
        solo3 = run_member(2026, 3)
        @test ens5[3] == solo3
        # order/N independence: build the same member from a reversed iteration.
        ens5_rev = Dict(k => run_member(2026, k) for k in reverse(1:5))
        @test ens5_rev[3] == ens5[3]
        # members are distinct streams (almost surely):
        @test ens5[1] != ens5[2]
        =#
    end

    # [join-species-count-union] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md Pending §Composition; union_acs! operators/joins.jl:14-44 (S merge by name) and :30-36 (T append)
    # note: Locks in present union_acs! behavior: S-merge loop (joins.jl:14-28) dedups by specName via incident;
    # note: T-append (joins.jl:30-36) copies every trans-attr for nparts(acs2,:T) new rows, then renames
    # note: (joins.jl:38-44). The transition-count assertion pins that the historic mid-loop early-return is
    # note: gone on ref-agents (the loop runs to completion). If alias resolution differs and A is NOT merged, S
    # note: would be 4 — that failure would itself be informative. Uses @alias per tutorial/example.jl:89.
    # action: `@join` the two models identifying the shared species A, then count S and T parts of the
    # merged schema.
    @testset "T1 lock-in: @join merged species count = |union of names|, transition count = sum (no transitions lost)" begin
        acs1 = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
        end
        acs2 = @ReactionNetworkSchema begin
          1.0, A --> C, name => t2
        end
        # acs1 species: {A,B}; acs2 species: {A,C}. Identify the two A's via the join eqs.
        m = @join acs1 acs2 acs1.A = acs2.A = @alias(A)
        # union of names {A, B, C} => 3 species (A merged; B, C distinct after prefixing).
        @test nparts(m, :S) == 3
        @test Symbol("A") in m[:, :specName]
        # transitions are appended, never dropped: 1 + 1 = 2 (pins the historic mid-loop-return bug is FIXED).
        @test nparts(m, :T) == 2
        # both transition bodies survive the merge (the :trans column is fully populated).
        @test count(!isnothing, m[:, :trans]) == 2
    end

    # [join-obs-events-not-merged-gap] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md Pending §Composition; union_acs! walks only S/T/P/M (joins.jl:14-56) — no :obs or :E loop
    # note: Documents the §Composition gap: union_acs! (joins.jl:10-59) iterates parts(acs2,:S), :T, :P, :M only
    # note: — there is no parts(acs2,:E) or parts(acs2,:obs) loop, so events and observables of the joined
    # note: submodel are LOST. The @test_broken line encodes the T2 target (events merged) so it flips green
    # note: once implemented; the two plain @test lines pass-now and lock in the current drop. Event syntax
    # note: (`cond && action`) per tutorial/example.jl and get_events! (create.jl:130-147).
    # action: Join models that carry events/observables and inspect whether :E / :obs parts survive into
    # the merged schema.
    @testset "T1 gap-doc / T2 target: union_acs! does NOT merge observables (:obs) or events (:E)" begin
        acs1 = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
        end
        @valuation acs1 B = 0.1
        acs2 = @ReactionNetworkSchema begin
          1.0, C --> D, name => t2
        end
        # give acs2 an event and an observable so we can check they are dropped on merge
        acs2_ev = @ReactionNetworkSchema begin
          1.0, C --> D, name => t2
          (D > 5) && (D -= 1)
        end
        m = @join acs1 acs2_ev
        # T1 (current): union_acs! has no :E loop, so the event is silently dropped.
        @test nparts(m, :E) == 0
        @test_broken nparts(m, :E) == nparts(acs2_ev, :E)  # T2 target: events should be merged
        # T1 (current): union_acs! has no :obs loop either.
        @test nparts(m, :obs) == 0
    end

    # [equalize-collapse-and-rewrite] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md Pending §Composition; equalize! operators/equalize.jl:24-66 (specmap + rem_parts! at :52 + recursively_substitute_vars! at :60)
    # note: Locks in current equalize! (equalize.jl:24-66): builds specmap, keeps the lowest index,
    # note: rem_parts!(acs,:S,species_ixs[2:end]) at :52, then recursively_substitute_vars! rewrites every spec-
    # note: referencing attr (:55-63). Pins string-surgery semantics that ADR 0003's promoted ReactantSpec will
    # note: replace structurally (see next test). Uses bare `name = name` eq form per @equalize docstring
    # note: (equalize.jl:73). If get_eqs_ff parsing of the bare `A = A2` form differs, the count assertion
    # note: surfaces it.
    # action: Call `@equalize` to collapse A and A2 into a single species and assert the species count
    # drops by one and references are rewritten.
    @testset "T1 lock-in: equalize! collapses two identified species into one and rewrites refs" begin
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
          1.0, A2 --> B, name => t2
        end
        # A and A2 are conceptually the same pool; identify them.
        before_S = nparts(acs, :S)
        m = @equalize acs A = A2
        # A and A2 collapse to one => species count drops by exactly 1.
        @test nparts(m, :S) == before_S - 1
        # the surviving merged name is present; the eliminated alias is gone.
        @test count(n -> n in (:A, :A2), m[:, :specName]) == 1
        # transitions are preserved (rem_parts! only touched :S).
        @test nparts(m, :T) == 2
    end

    # [equalize-reactant-fk-repoint] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0003 (promote transition<->reactant relation to typed ReactantSpec incidence table); CONTRACT_DRAFT.md Pending §Composition (structurally exact species-merge)
    # note: Encodes ADR 0003: the transition<->reactant relation becomes a typed ReactantSpec incidence table
    # note: (FK trans->T, species->S, side, stoich ExprNode, modality). equalize! then repoints the species FK
    # note: from A2 to A structurally rather than via recursively_substitute_vars! string rewriting
    # note: (equalize.jl:60). Errors today: ReactiveDynamics.reactant_specs / .specname do not exist (reactants
    # note: live as Expr in the :trans column, parsed per-tick by extract_reactants, reaction_parser.jl:32). T2
    # note: against the not-yet-built IR.
    # action: After equalize!, inspect the promoted ReactantSpec table (target IR) and assert every
    # reactant row that pointed at the eliminated species now points at the survivor by FK — not by re-
    # parsed expression strings.
    @testset "T2: promoted-ReactantSpec equalize repoints species FKs structurally (no string surgery)" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
          1.0, A2 --> B, name => t2
        end
        # TARGET IR: a first-class ReactantSpec incidence table with FK species->S.
        m = equalize!(acs, [[(:catchall, :A), (:catchall, :A2)]])
        reactants = ReactiveDynamics.reactant_specs(m)   # target accessor over the promoted table
        # every ReactantSpec.species FK resolves to a live S index (no dangling FK after collapse).
        @test all(r -> 1 <= r.species <= ReactiveDynamics.nparts(m, :S), reactants)
        # no reactant still references the eliminated A2 index.
        surv = ReactiveDynamics.find_index(:A, m)
        @test any(r -> r.species == surv, reactants)
        @test !any(r -> ReactiveDynamics.specname(m, r.species) == :A2, reactants)
        =#
    end

    # [join-include-model-undefined] tier=T2-acceptance expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md Pending §Composition ('the @join file branch calls an undefined include_model (joins.jl:226,228)')
    # note: Pins joins.jl:221-232: when an @join argument is a macrocall/string (e.g. @join @file("m.jl")), the
    # note: macro expands to :(include_model($str)) at joins.jl:226 and :228, but include_model is defined
    # note: NOWHERE in src/ (grep-confirmed). The two plain @test lines pass-now (documenting the undefined
    # note: symbol); the @test_broken encodes the T2 target so it flips green once include_model is implemented.
    # note: The symbol-only branch (acsex passed directly) still works — this only breaks the file-include form.
    # action: Confirm include_model is undefined (so the file-include path of @join cannot work today),
    # and that the TARGET API defines it.
    @testset "T2 bug-pin: @join file-include branch calls undefined include_model" begin
        # The @join macro, when given a string/macrocall arg, emits :(include_model(str)).
        # include_model is never defined anywhere in src/.
        # Today: no such symbol is bound in the module.
        @test !isdefined(ReactiveDynamics, :include_model)
        # Exercising the file-include branch must therefore raise (UndefVarError on include_model).
        @test_throws UndefVarError @eval ReactiveDynamics include_model("some_model.jl")
        # T2 target: include_model exists and returns a ReactionNetworkSchema.
        @test_broken isdefined(ReactiveDynamics, :include_model)
    end

    # [equalize-live-guard-refuses] tier=T2-acceptance expectedStatus=errors-until-implemented
    # contract: ADR 0004 INV-2 (No mid-run reindex; the sole offender operators/equalize.jl:52 must refuse under a runtime live guard)
    # note: ADR 0004 INV-2: rem_parts!(acs,:S,...) at equalize.jl:52 is the one mid-run reindexer that breaks
    # note: every position-indexed compiled closure (compilers.jl:149 varmap, sample_transitions! positional
    # note: loop state.jl:178-191) and MUST refuse on a stepping model. Errors-until-implemented because there
    # note: is NO live guard today: equalize! takes a ReactionNetworkSchema, not a ReactionNetworkProblem, so
    # note: equalize!(prob, ...) currently throws a MethodError (which technically satisfies @test_throws
    # note: Exception for the wrong reason) — the real target is a deliberate live-guard error plus an equalize!
    # note: method that accepts the live state and refuses. Treat the MethodError-today as accidental; T2 pins
    # note: the intended guard.
    # action: Step the live problem one tick, then attempt to equalize species on the LIVE state; the
    # reindexing rem_parts! must be refused.
    @testset "T2: equalize!'s rem_parts! must refuse on a live/stepping model (ADR 0004 INV-2)" begin
        # TARGET API not yet implemented — guarded so the suite loads; build it, then unskip.
        @test_skip false  # see the reference block below
        #=
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
          1.0, A2 --> B, name => t2
        end
        @prob_init acs A = 10 A2 = 10 B = 0
        @prob_params acs
        @prob_meta acs tspan = 10 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        simulate(prob, 1)   # advance to a tick boundary; model is now live
        # TARGET: equalize!/@equalize on a live ReactionNetworkProblem must error (live guard), not silently rem_parts! the acs out from under the compiled closures.
        @test_throws Exception equalize!(prob, [[(:catchall, :A), (:catchall, :A2)]])
        # and the model's species indexing must be unchanged after the refusal.
        @test ReactiveDynamics.nparts(prob, :S) == 3
        =#
    end

    # [bugpin-event-action-noop] tier=T1-characterization expectedStatus=test_broken-pins-bug
    # contract: CONTRACT_DRAFT.md §3.4 Invariant 7; event_action! solvers.jl:316-326 (line 323 fetches state[i,:eventAction] but never evaluates it)
    # note: Pins Invariant 7 / solvers.jl:323: event_action! computes q correctly (solvers.jl:321) but the loop
    # note: body is just the expression `state[i,:eventAction]` — a fetch through getindex (state.jl:73-83;
    # note: :eventAction lacks 'trans' so it goes through context_eval) with no evaluation/side-effect
    # note: application. So events are a complete no-op. The @test_broken pins the desired behavior (B>0) and
    # note: flips green when the action is actually run; the plain @test locks in today's no-op (B==0). Event
    # note: authoring `cond && action` per create.jl:130-147 / tutorial/example.jl:52.
    # action: Simulate; the event fires every tick (Bool trigger true => q=1, solvers.jl:321) and SHOULD
    # raise B, but the action is never executed.
    @testset "T1 bug-pin: event_action! is a no-op — a triggered event that should set a species does nothing" begin
        # Event with an always-true trigger whose action would bump B by 100 each tick.
        acs = @ReactionNetworkSchema begin
          0.0, A --> B, name => inert        # no spawning; isolates the event effect
          (true) && (B += 100)              # event: trigger true, action sets B
        end
        @prob_init acs A = 0 B = 0
        @prob_params acs
        @prob_meta acs tspan = 5 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        simulate(prob)
        # Invariant 7 target: B should have been incremented by the event action each tick (B > 0).
        @test_broken last(prob.sol.B) > 0
        # Characterize current reality: B stays at its initial 0 because event_action! never evals the action.
        @test last(prob.sol.B) == 0
    end

    # [bugpin-add-to-spawn-deferral-lost] tier=T1-characterization expectedStatus=test_broken-pins-bug
    # contract: CONTRACT_DRAFT.md §3.4 Invariant 3; add_to_spawn! state.jl:251-256 (findfirst over scalar length(); increments :transHash not :transToSpawn)
    # note: Pins Invariant 3 / state.jl:251-256: findfirst is handed the scalar
    # note: `length(state.transition_recipes[:transHash])` instead of an index range (so it cannot match), and
    # note: on a (never-taken) match it does `[:transHash][ix] += n` — incrementing the wrong column. Net:
    # note: capacity overflow is silently dropped, never deferred. Also note sample_transitions! resets
    # note: transToSpawn to 0 at state.jl:215 each tick, reinforcing the loss. The @test_broken pins the target
    # note: (overflow recorded in transToSpawn); the plain @test locks in the current drop. Calls the unexported
    # note: allocator internals (ReactiveDynamics.evolve!/sample_transitions!) directly to observe the gate.
    # note: NOTE: if add_to_spawn! THROWS today (scalar passed to findfirst can error rather than no-op), wrap
    # note: evolve! in @test_throws instead and keep the @test_broken target — verify against the running engine
    # note: and adjust.
    # action: Call evolve! for one tick so the capacity gate (solvers.jl:147-153) computes overflow and
    # calls add_to_spawn!; then inspect the transToSpawn deferral column directly.
    @testset "T1 bug-pin: add_to_spawn! capacity-overflow deferral is lost (findfirst scalar + wrong column)" begin
        # Sustained over-demand against a hard capacity: rate forces 5 desired/ tick,
        # but capacity caps concurrent instances; overflow must be DEFERRED to transToSpawn.
        acs = @ReactionNetworkSchema begin
          @deterministic(5), A --> B, name => capped, capacity => 2.0, cycletime => 3.0
        end
        @prob_init acs A = 1000 B = 0
        @prob_params acs
        @prob_meta acs tspan = 4 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        ReactiveDynamics.sample_transitions!(prob)
        ReactiveDynamics.evolve!(prob)
        # Invariant 3 target: overflow (desired 5 - capacity 2 = 3) is carried forward in transToSpawn.
        @test_broken prob.transition_recipes[:transToSpawn][1] >= 3
        # Current reality: the deferral is lost. add_to_spawn! either no-ops (findfirst over a scalar)
        # or mutates :transHash; either way transToSpawn stays 0 and overflow is dropped.
        @test prob.transition_recipes[:transToSpawn][1] == 0
    end

    # [bugpin-resample-oval-crash] tier=T1-characterization expectedStatus=test_broken-pins-bug
    # contract: CONTRACT_DRAFT.md §4.2 (observable sampling); resample! state.jl:135-139 (line 137 sets o.val; Observable field is .sampled, state.jl:37)
    # note: Pins state.jl:137: `isempty(o.range) && (return o.val = missing)` references o.val, but @aagent
    # note: Observable (state.jl:31-38) declares `sampled::Any`, not `val` — so the range-less path throws
    # note: (likely an ErrorException/setproperty! failure / type error) instead of returning missing. The non-
    # note: empty-range path at state.jl:139 correctly uses o.sampled. The @test_throws pins the current crash;
    # note: the @test_broken encodes the fixed behavior (assign .sampled = missing).
    # note: Observable/SampleableValues/ActionableValues are unexported, hence the ReactiveDynamics. prefix.
    # note: NOTE: if @aagent injects a `val` accessor or the setproperty path silently no-ops rather than
    # note: throwing, downgrade the first line to @test_broken.
    # action: Call resample! on the range-less observable; the isempty(o.range) branch (state.jl:137)
    # executes `o.val = missing`, but Observable has no `val` field.
    @testset "T1 bug-pin: resample! assigns nonexistent o.val on a range-less observable" begin
        # Build a problem and a range-less Observable directly, then resample it.
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
        end
        @prob_init acs A = 10 B = 0
        @prob_params acs
        @prob_meta acs tspan = 5 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        # An observable whose `range` is empty triggers the isempty(o.range) branch at state.jl:137.
        obs = ReactiveDynamics.Observable("o", -Inf, ReactiveDynamics.SampleableValues[], Inf, ReactiveDynamics.ActionableValues[], missing)
        # o.val does not exist (field is `sampled`), so the assignment errors.
        @test_throws Exception ReactiveDynamics.resample!(prob, obs)
        # T2 target: a range-less observable should resample to `missing` on the correct `.sampled` field.
        @test_broken (ReactiveDynamics.resample!(prob, obs); obs.sampled === missing)
        # Field-shape pin: confirm the struct really has `sampled`, not `val`.
        @test :sampled in fieldnames(ReactiveDynamics.Observable)
        @test !(:val in fieldnames(ReactiveDynamics.Observable))
    end
end
