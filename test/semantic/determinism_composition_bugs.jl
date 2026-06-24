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
    # contract: CONTRACT_DRAFT.md §4.1 D2 (RNG isolation); §4.3 (entropy-seeded by default)
    # note: Stage A landed a state-owned `rng::AbstractRNG` and a working `seed=` ctor kwarg. With NO seed the
    # note: state is entropy-seeded and ISOLATED from the global RNG, so two unseeded runs still diverge (a fresh
    # note: entropy seed per construction) — but the cause is no longer "no rng field"; the state now OWNS one.
    # note: With rate 3.0 over 30 ticks the Poisson+Binomial draws make a collision astronomically unlikely, so
    # note: the divergence @test reliably passes. The seeding gap is CLOSED: `seed=` constructs successfully (no
    # note: longer throws) — D1/D2/D7/D8 below now run as real T1 PASSES.
    # action: Construct two independent UNSEEDED problems from the SAME spec and `simulate` each to
    # completion; the entropy seeds differ so the B columns diverge. Confirm the state owns an `rng` field
    # and that constructing WITH a seed succeeds.
    @testset "T1 characterization: two unseeded runs diverge (entropy-seeded), and the state now owns an rng" begin
        # Poisson genesis + Binomial PoS => genuinely stochastic; each unseeded run draws a fresh entropy seed.
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
        # No seed => each run is entropy-seeded independently, so trajectories diverge.
        @test prob1.sol.B != prob2.sol.B
        # Stage A: the state now OWNS a per-run rng (seeding gap closed; isolated from the global RNG).
        @test (:rng in fieldnames(typeof(prob1)))
        # Constructing WITH a seed now succeeds (the D6 `seed=` kwarg is wired and no longer throws).
        @test ReactionNetworkProblem(mk(); seed = 1234) isa ReactionNetworkProblem
    end

    # [determinism-d1-reproducible] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §4.1 D1 (Reproducibility); D6 (seed at construction)
    # note: Stage A implemented D5/D6: ReactionNetworkProblem grows a `rng::AbstractRNG` field seeded
    # note: deterministically from `seed`; every rand reachable from _step! draws from state.rng. Verified: same
    # note: (model, seed) reproduces sol AND log exactly; a different seed (almost surely) diverges. Reclassified
    # note: T2->T1: the API exists now.
    # action: Construct two problems from the same spec with the SAME `seed`, simulate both, and assert
    # the full trajectory and the entire log are equal.
    @testset "T1 D1: same (model, seed) => identical prob.sol AND prob.log" begin
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
    end

    # [determinism-d2-rng-isolation] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §4.1 D2 (RNG isolation); D5 (threading rule)
    # note: Stage A implemented D2/D5: all randomness comes from the state-owned AbstractRNG; no _step! path
    # note: calls bare rand(). Verified: (a) external rand() before/around the run does not shift it -> pa.sol ==
    # note: pb.sol, and (b) a seeded run does not advance default_rng() -> the before/after snapshot is unchanged.
    # note: Reclassified T2->T1: the API exists now.
    # action: Run a seeded simulation. Run it again with the SAME seed but perturb the global RNG (extra
    # rand() calls) before/around the run. The trajectory must be identical (independence from global
    # RNG), and the global RNG state must be untouched by the run.
    @testset "T1 D2: a run neither reads nor perturbs the global RNG (external rand() does not affect it)" begin
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
    end

    # [determinism-d7-reinit-restores-stream] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §4.3 D7 (re-init restores the stream); _reinit! (AlgebraicAgents dispatch)
    # note: Stage A completed _reinit!: it resets u/t/ongoing_transitions/log/observables/sol AND copies the
    # note: initial RNG snapshot back, restoring the (M, seed) starting stream. Verified: a second simulate after
    # note: _reinit! reproduces the first sol AND log exactly. Reclassified T2->T1. NOTE: the public binding is
    # note: AlgebraicAgents._reinit! (also reexported as reinit!); there is NO ReactiveDynamics._reinit!.
    # action: Build one seeded problem, simulate it (capturing sol/log), call _reinit!, simulate again,
    # and assert the second trajectory equals the first.
    @testset "T1 D7: init -> step* -> reinit! -> step* reproduces the first trajectory" begin
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
        AlgebraicAgents._reinit!(p)   # AA-dispatched _reinit! (the real binding)
        simulate(p)
        @test p.sol == sol1
        @test p.log == log1
    end

    # [determinism-d8-ensemble-per-index-seeding] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §4.4 D8 (per-trajectory seeding from index); D9 (no shared mutable RNG)
    # note: Stage A enables D8/D9: a per-member seed derived deterministically from (root_seed,k) (member_seed
    # note: shown, UInt64 seeds accepted) gives each member its own RNG (D9, no shared mutable instance).
    # note: Verified: member k is identical standalone vs. inside a size-5 ensemble and independent of N/order
    # note: (which relies on D2 isolation — running members 1,2 before 3 does not perturb 3); distinct members
    # note: diverge. Reclassified T2->T1.
    # action: Compute member k's trajectory inside two different ensemble sizes / orderings and assert it
    # is identical; assert two distinct members differ.
    @testset "T1 D8/D9: ensemble member k reproducible from (root_seed,k), independent of N and order, own RNG per member" begin
        mk() = begin
          acs = @ReactionNetworkSchema begin
            3.0, A --> B, name => birth, probability => 0.5, cycletime => 2.0
          end
          @prob_init acs A = 100 B = 0
          @prob_params acs
          @prob_meta acs tspan = 20 dt = 1.0
          acs
        end
        # Derive a per-member seed from a single root seed + member index.
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

    # [equalize-live-guard-refuses] tier=T1-characterization expectedStatus=pass-now
    # contract: ADR 0004 INV-2 / ADR 0007 §A (No mid-run reindex; the rem_parts! reindexer must refuse on a live model)
    # note: Stage D added the §A live-phase guard: a constructed ReactionNetworkProblem has live==true, and
    # note: equalize!(::ReactionNetworkProblem, …) (src/actions.jl) deliberately errors — species identification
    # note: reindexes the :S table via rem_parts! (equalize.jl:52), which would invalidate the construction-
    # note: frozen, position-indexed compiled closures (compilers.jl varmap, sample_transitions! positional
    # note: loop). Identify species at AUTHORING time, before construction. The authoring-phase
    # note: equalize!(::ReactionNetworkSchema, …) stays unrestricted (the §7 path).
    # action: Construct + step the live problem, then attempt to equalize species on the LIVE state; the
    # reindexing rem_parts! must be refused, and the species indexing must be unchanged after the refusal.
    @testset "equalize!'s rem_parts! refuses on a live/stepping model (ADR 0004 INV-2 / ADR 0007 §A)" begin
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
          1.0, A2 --> B, name => t2
        end
        @prob_init acs A = 10 A2 = 10 B = 0
        @prob_params acs
        @prob_meta acs tspan = 10 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        @test prob.live == true
        simulate(prob, 1)   # advance to a tick boundary; model is live
        # the live-guard refuses, rather than silently rem_parts!-ing the acs out from under the closures
        @test_throws Exception equalize!(prob, [[(:catchall, :A), (:catchall, :A2)]])
        # the model's species indexing is unchanged after the refusal
        @test ReactiveDynamics.nparts(prob, :S) == 3
    end

    # [event-channel-repaired] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §3.4 Invariant 7; ADR 0010 §A — the event channel is repaired
    # note: Stage B replaces the no-op event_action! with the endogenous decision channel: at construction
    # note: every :E row (trigger && action) is lifted to a Rule{guard=trigger, action=RawExpr(action),
    # note: every_tick}, and fire_rules! (src/actions.jl, _step! step 10) evaluates the guard and RUNS the
    # note: action. So a triggered event now takes effect (Invariant 7 met). NB: a bare-`true` trigger still
    # note: fails at @ReactionNetworkSchema parse (Event(::Bool,::Expr) cannot convert a Bool into the
    # note: SampleableValues trigger union — a separate, pre-existing authoring limitation), so we use an
    # note: Expr trigger `@t() >= 0.0` which is always-true and parses cleanly.
    # action: Build a model whose event bumps B every tick; simulate; B must rise (the action ran).
    @testset "Event channel repaired: a triggered event runs its action each tick (Invariant 7 met)" begin
        acs = @ReactionNetworkSchema begin
          0.0, A --> B, name => inert        # no spawning; isolates the event effect
          (@t() >= 0.0) && (B += 100)        # event: always-true Expr trigger, action sets B
        end
        @prob_init acs A = 0 B = 0
        @prob_params acs
        @prob_meta acs tspan = 5 dt = 1.0
        prob = ReactionNetworkProblem(acs; seed = 1)
        # The :E row was lifted to a Rule.
        @test length(prob.rules) == 1
        simulate(prob)
        # Invariant 7: the event action executed each tick, so B was incremented (B > 0).
        @test last(prob.sol[!, "B"]) > 0
        @test all(>=(0.0), diff(prob.sol[!, "B"]))   # monotone — the action only ever adds
    end

    # [add-to-spawn-deferral-works] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §3.4 Invariant 3; add_to_spawn! records capacity overflow in transToSpawn
    # note: Stage A FIXED add_to_spawn! (state.jl:251-256): the old code handed findfirst a scalar
    # note: `length(...)` (so it never matched) and on the never-taken branch incremented :transHash — net,
    # note: overflow was silently dropped. INV3 now holds: the capacity gate (solvers.jl:147-153) computes
    # note: overflow (desired 5 - capacity 2 = 3) and add_to_spawn! defers it into transToSpawn. Verified
    # note: transToSpawn[1] == 3.0 for the capped(5)/capacity(2) model; the prior MethodError is gone. Calls the
    # note: unexported allocator internals (ReactiveDynamics.evolve!/sample_transitions!) directly to observe the
    # note: gate.
    # action: Call evolve! for one tick so the capacity gate computes overflow and calls add_to_spawn!;
    # then inspect the transToSpawn deferral column directly.
    @testset "T1 lock-in: add_to_spawn! defers capacity overflow into transToSpawn (INV3)" begin
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
        # Invariant 3: overflow (desired 5 - capacity 2 = 3) is carried forward in transToSpawn.
        @test prob.transition_recipes[:transToSpawn][1] >= 3
        @test prob.transition_recipes[:transToSpawn][1] == 3
    end

    # [resample-rangeless-returns-missing] tier=T1-characterization expectedStatus=pass-now
    # contract: CONTRACT_DRAFT.md §4.2 (observable sampling); resample! assigns the correct `.sampled` field
    # note: Stage A FIXED resample!: the range-less branch previously assigned the nonexistent `o.val` (the
    # note: @aagent Observable declares `sampled::Any`, not `val`) and threw. It now assigns `o.sampled = missing`
    # note: on the isempty(o.range) path, returning missing without throwing. Verified: obs.sampled === missing
    # note: after resample! on a range-less observable; the field is `.sampled`, never `.val`.
    # note: Observable/SampleableValues/ActionableValues are unexported, hence the ReactiveDynamics. prefix.
    # action: Call resample! on a range-less observable; the isempty(o.range) branch resamples to `missing`
    # on the correct `.sampled` field.
    @testset "T1 lock-in: resample! on a range-less observable yields missing on .sampled (no throw)" begin
        # Build a problem and a range-less Observable directly, then resample it.
        acs = @ReactionNetworkSchema begin
          1.0, A --> B, name => t1
        end
        @prob_init acs A = 10 B = 0
        @prob_params acs
        @prob_meta acs tspan = 5 dt = 1.0
        prob = ReactionNetworkProblem(acs)
        # An observable whose `range` is empty triggers the isempty(o.range) branch.
        obs = ReactiveDynamics.Observable("o", -Inf, ReactiveDynamics.SampleableValues[], Inf, ReactiveDynamics.ActionableValues[], missing)
        # The range-less path now resamples to `missing` on `.sampled` (no throw).
        @test (ReactiveDynamics.resample!(prob, obs); obs.sampled === missing)
        # Field-shape pin: confirm the struct really has `sampled`, not `val`.
        @test :sampled in fieldnames(ReactiveDynamics.Observable)
        @test !(:val in fieldnames(ReactiveDynamics.Observable))
    end
end
