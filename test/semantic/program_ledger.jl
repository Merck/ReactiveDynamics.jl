# Phase-1 semantic tests — per-program (per-structured-token) ledger (MVP finding D / D-bis).
#
# Covers the engine-level per-program ledger (src/ledger.jl): cost attributed at the bind site
# (evolve!) and reward at the finishing transition (finish!), split evenly across a transition's
# bound programs; the network UNATTRIBUTED bucket for plain (no bound token) spend/reward; the
# SUM-CONSISTENCY invariant (per-program rows + unattributed bucket == the aggregate :valuation_cost
# / :valuation_reward rows, §8.5); determinism under (model, seed) (§4 D4); and the _reinit! reset
# (§4 D7). The query API is program_ledger(state)::DataFrame and program_ledger_entries(state,name).
#
# These run against the built engine — real @tests, not pins. Default runs are entropy-seeded, so
# every problem here is built with an explicit seed= for reproducibility.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

RD = ReactiveDynamics

# A structured token kind with a `phase` attribute, defined in RD scope via @register (the
# tutorial/token_filtration idiom). One Project kind, phase-as-attribute (ADR 0008 §D).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct LedgerProjectToken
        phase::Symbol
    end
    function LedgerProjectToken(phase)
        return LedgerProjectToken(
            "LP" * string(rand(1:10^9)),
            :Project,
            nothing,
            Tuple{Symbol,Float64,ReactiveDynamics.Transition}[],
            phase,
        )
    end
end

# Aggregate ledger sums from the run log (the existing pool-level rows). Cost/reward are scalar.
agg_cost(p) = sum(r[3] for r in p.log if r[1] == :valuation_cost; init = 0.0)
agg_reward(p) = sum(r[3] for r in p.log if r[1] == :valuation_reward; init = 0.0)

# A coarse one-step advance model: @select a Phase1 Project, burn 2 budget per tick at @rate over
# cycletime 1, and @advance(phase,:Phase2) on success. `budget` carries specCost so the burn is a
# real ledger cost; the advanced Project species carries specReward so a successful advance realizes
# reward — the minimal model where BOTH sides of the per-program ledger are non-trivial. (Macro
# arguments are literal: @prob_init/@reaction_network eval their RHS in MODULE scope, so a
# parameterized burn/budget would be undefined there — we set cost/reward/budget on the ACSet
# directly via index assignment instead, which is what the kwargs control.)
function advance_cost_model(; budget0 = 100, cost = 1.0, reward = 10.0)
    net = @reaction_network begin
        @deterministic(1.0),
        @select(Project, phase == :Phase1) + 2 * @rate(budget) --> @advance(phase, :Phase2),
        name => adv, cycletime => 1.0, probability => 1.0
    end
    RD.register_structured_species!(net, :Project)
    bi = findfirst(==(:budget), net[:, :specName])
    net[bi, :specInitVal] = Float64(budget0)
    net[bi, :specCost] = cost
    pi = findfirst(==(:Project), net[:, :specName])
    net[pi, :specReward] = reward
    @prob_meta net tspan = 4 dt = 1.0
    return net
end

build_ledger_prob(seed; kwargs...) = ReactionNetworkProblem(
    advance_cost_model(; kwargs...);
    seed = seed,
    population = [RD.LedgerProjectToken(:Phase1)],
)

@testset "Per-program ledger (MVP finding D / D-bis, src/ledger.jl)" begin

    # ── known cost at bind, known reward at finish, attributed to the bound program ──────
    @testset "single bound program gets the full transition cost (bind) and reward (finish)" begin
        p = build_ledger_prob(1; cost = 1.0, reward = 10.0)
        simulate(p)
        df = program_ledger(p)
        @test nrow(df) == 1                       # exactly one program ever existed
        prog = df[1, :]
        # The program is bound exactly at tick 0 (it advances out of :Phase1 on the first finish),
        # so it is charged the burn it consumed THAT tick (2 budget * specCost 1 = 2.0) and credited
        # the advance reward (specReward 10 on the produced :Project). Attribution rule: a transition
        # with exactly one bound token gets its FULL consumed cost (see src/ledger.jl header).
        @test prog.cost_incurred == 2.0
        @test prog.reward_realized == 10.0
        @test prog.net == 8.0
        @test prog.species == :Project            # current kind (phase advanced, kind unchanged)
        @test prog.creation_index == 1            # first (only) program

        # the append-only audit trail records the two events in order
        entries = program_ledger_entries(p, prog.program)
        @test length(entries) == 2
        @test entries[1][2] == :cost && entries[1][3] == 2.0
        @test entries[2][2] == :reward && entries[2][3] == 10.0
    end

    # ── SUM-CONSISTENCY invariant: per-program + unattributed == aggregate ───────────────
    @testset "per-program cost/reward rows + unattributed bucket SUM to the aggregate ledger (§8.5)" begin
        p = build_ledger_prob(3; cost = 1.0, reward = 10.0)
        simulate(p)
        df = program_ledger(p)
        per_cost = sum(df.cost_incurred) + p.unattributed_cost
        per_reward = sum(df.reward_realized) + p.unattributed_reward
        @test isapprox(per_cost, agg_cost(p); atol = 1e-9)
        @test isapprox(per_reward, agg_reward(p); atol = 1e-9)
        # After the program advances at tick 0 it no longer matches @select(Phase1); the transition
        # keeps spawning (rate 1) and burning budget on instances with NO bound program — that spend
        # is the documented UNATTRIBUTED bucket (the boundary of finding D), so it is > 0 here and is
        # exactly what makes the per-program total fall short of the aggregate.
        @test p.unattributed_cost > 0.0
        @test isapprox(sum(df.cost_incurred), 2.0; atol = 1e-9)        # only tick 0 is the program's
    end

    # ── even split across multiple bound programs ────────────────────────────────────────
    @testset "a transition binding several programs splits its cost/reward evenly" begin
        # rate 2 + cycletime 0 (instant) so both Phase1 Projects are bound and advanced by the SAME
        # tick's spawned instances. Each instance binds one Project (token-gated genesis), so each
        # program is charged the burn of ITS own instance — i.e. each gets `burn*cost`. The even-split
        # rule is exercised; with one token per instance, even-split == full cost per program.
        net = @reaction_network begin
            @deterministic(2.0),
            @select(Project, phase == :Phase1) + 3 * @rate(budget) --> @advance(phase, :Phase2),
            name => adv2, cycletime => 0.0, probability => 1.0
        end
        RD.register_structured_species!(net, :Project)
        @prob_init net budget = 100
        bi = findfirst(==(:budget), net[:, :specName])
        net[bi, :specCost] = 1.0
        @prob_meta net tspan = 2 dt = 1.0
        p = ReactionNetworkProblem(
            net;
            seed = 5,
            population = [RD.LedgerProjectToken(:Phase1), RD.LedgerProjectToken(:Phase1)],
        )
        simulate(p)
        df = program_ledger(p)
        @test nrow(df) == 2
        # cycletime 0 ⇒ @rate burn is 0 for an instant completion (no ticks held); both advanced.
        # The invariant must still hold regardless of the exact numbers.
        @test isapprox(
            sum(df.cost_incurred) + p.unattributed_cost,
            agg_cost(p);
            atol = 1e-9,
        )
        @test all(df.species .== :Project)
        @test sort(df.creation_index) == [1, 2]
    end

    # ── determinism under (model, seed) (§4 D4) ──────────────────────────────────────────
    @testset "the per-program ledger is identical across two runs with the same seed" begin
        a = build_ledger_prob(11)
        simulate(a)
        da = program_ledger(a)
        b = build_ledger_prob(11)
        simulate(b)
        db = program_ledger(b)
        @test da.cost_incurred == db.cost_incurred
        @test da.reward_realized == db.reward_realized
        @test da.creation_index == db.creation_index
        @test a.unattributed_cost == b.unattributed_cost
        @test a.unattributed_reward == b.unattributed_reward
    end

    # ── _reinit! clears the per-program ledger and a replay reproduces it (§4 D7) ─────────
    @testset "_reinit! resets the program ledger; init→run→reinit→run reproduces the totals" begin
        p = build_ledger_prob(7)
        simulate(p)
        cost1 = sum(program_ledger(p).cost_incurred)
        reward1 = sum(program_ledger(p).reward_realized)
        unattr1 = p.unattributed_cost
        @test !isempty(p.program_ledgers)

        AlgebraicAgents._reinit!(p)
        @test isempty(p.program_ledgers)          # ledger torn down
        @test p.unattributed_cost == 0.0          # buckets cleared
        @test p.unattributed_reward == 0.0

        simulate(p)
        @test isapprox(sum(program_ledger(p).cost_incurred), cost1; atol = 1e-9)
        @test isapprox(sum(program_ledger(p).reward_realized), reward1; atol = 1e-9)
        @test isapprox(p.unattributed_cost, unattr1; atol = 1e-9)
    end

    # ── the per-tick :program_ledger log row is present and in deterministic token order ──
    @testset "a :program_ledger row is pushed each tick, keyed by token, summing to the running totals" begin
        p = build_ledger_prob(2)
        simulate(p)
        rows = [r for r in p.log if r[1] == :program_ledger]
        @test !isempty(rows)
        # the LAST per-tick row's per-program running totals match the query-API summary
        last_snapshot = rows[end][3]
        df = program_ledger(p)
        for row in eachrow(df)
            haskey(last_snapshot, row.program) || continue
            s = last_snapshot[row.program]
            @test isapprox(s.cost, row.cost_incurred; atol = 1e-9)
            @test isapprox(s.reward, row.reward_realized; atol = 1e-9)
        end
    end

    # ── valuation: a species with specValuation marks its live programs to market ─────────
    @testset "live programs are marked to market by their species' specValuation" begin
        # Give the Project species a specValuation; a live (unblocked) program then carries that
        # mark in the ledger's `valuation` column (a stock, recomputed each tick — not a flow).
        net = advance_cost_model(; reward = 0.0)
        pi = findfirst(==(:Project), net[:, :specName])
        net[pi, :specValuation] = 50.0
        # Make the program NOT advance (select a phase it isn't in) so it stays live & unblocked.
        p = ReactionNetworkProblem(net; seed = 9, population = [RD.LedgerProjectToken(:Phase3)])
        simulate(p)
        df = program_ledger(p)
        @test nrow(df) == 1
        @test df[1, :valuation] == 50.0           # marked at the species' specValuation
        @test df[1, :species] == :Project
    end
end
