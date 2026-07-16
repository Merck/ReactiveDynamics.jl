# Phase-1 Stage C semantic tests — token filtration (ADR 0008, CONTRACT §9.5).
#
# Covers: @select(Kind, clauses) predicate-based LHS binding (TokenPredicate); the
# 𝓕ₜ-measurable matches() filter in front of the unchanged priority/creation-index sort;
# phase-as-attribute as canonical (one Project kind, `phase` field); @advance(field, value)
# field write (SetField) preserving token identity; degenerate (no-predicate) backward-compat.
# These are the demo's phase-pipeline pillar (MVP_BD_DEMO.md §2/§7 / maintainer Q4).
#
# All tests run against the built engine (Stage C) — real @tests, not pins.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

# A ProjectToken structured kind with a `phase` attribute and an `npv` field, defined in RD scope
# via @register (the tutorial M2 idiom — @structured_token's bare-name doc ref fails outside RD).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct FiltProjectToken
        phase::Symbol
        npv::Float64
    end
    function FiltProjectToken(phase, npv)
        return FiltProjectToken(
            "FP" * string(rand(1:(10^9))),
            :Project,
            nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[],
            phase,
            npv,
        )
    end
end

# A coarse phase-advance model: select Phase2 projects, advance them to Phase3 on completion.
function advance_model()
    net = @reaction_network begin
        @deterministic(1.0),
            @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
            name => p2_to_p3, cycletime => 1.0, probability => 1.0
    end
    RDX.register_structured_species!(net, :Project)
    return net
end

phases(p) =
    sort(string.([t.phase for t in values(RDX.inners(RDX.getagent(p, "structured")))]))

@testset "Token filtration (@select / TokenPredicate + @advance / SetField, ADR 0008)" begin

    # ── matches(): the 𝓕ₜ-measurable predicate evaluator ────────────────────────────────
    @testset "matches() selects by attribute clause (==, >), kind-gated" begin
        net = advance_model()
        @prob_meta net tspan = 3 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1)
        t2 = RDX.FiltProjectToken(:Phase2, 100.0); add_structured_token!(p, t2)
        t1 = RDX.FiltProjectToken(:Phase1, 50.0); add_structured_token!(p, t1)
        pred_eq = RDX.TokenPredicate(:Project, [RDX.Clause(:phase, :(==), :(:Phase2))])
        @test RDX.matches(pred_eq, t2, p, nothing) == true
        @test RDX.matches(pred_eq, t1, p, nothing) == false
        # continuous-field clause: npv > 75
        pred_gt = RDX.TokenPredicate(:Project, [RDX.Clause(:npv, :(>), 75.0)])
        @test RDX.matches(pred_gt, t2, p, nothing) == true     # 100 > 75
        @test RDX.matches(pred_gt, t1, p, nothing) == false    # 50 > 75 is false
        # nothing predicate is the degenerate "any token" (backward-compat)
        @test RDX.matches(nothing, t1, p, nothing) == true
    end

    # ── @select binds only matching tokens; @advance writes the phase field ─────────────
    @testset "phase-as-attribute pipeline: @select(Phase2) --> @advance(phase,:Phase3)" begin
        net = advance_model()
        @prob_meta net tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1)
        add_structured_token!(p, RDX.FiltProjectToken(:Phase2, 100.0))
        add_structured_token!(p, RDX.FiltProjectToken(:Phase2, 200.0))
        add_structured_token!(p, RDX.FiltProjectToken(:Phase1, 50.0))
        @test phases(p) == ["Phase1", "Phase2", "Phase2"]
        simulate(p)
        # both Phase2 projects advanced to Phase3; the Phase1 project is untouched (not selected)
        @test phases(p) == ["Phase1", "Phase3", "Phase3"]
    end

    # ── @advance preserves token identity (uuid/kind/creation_index/past_bonds) ──────────
    @testset "@advance keeps token identity — only the field changes (not a new token)" begin
        net = advance_model()
        @prob_meta net tspan = 4 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1)
        tok = RDX.FiltProjectToken(:Phase2, 100.0)
        add_structured_token!(p, tok)
        name_before = AlgebraicAgents.getname(tok)
        ci_before = p.creation_index[name_before]
        simulate(p)
        @test tok.phase == :Phase3                                  # same object, advanced
        @test AlgebraicAgents.getname(tok) == name_before           # identity preserved
        @test get(p.creation_index, name_before, -1) == ci_before   # creation_index stable
        @test RDX.get_species(tok) == :Project                      # kind unchanged
    end

    # ── continuous predicate selects a subset ───────────────────────────────────────────
    @testset "@select with a continuous clause (npv > θ) advances only the qualifying subset" begin
        net = @reaction_network begin
            @deterministic(1.0),
                @select(Project, phase == :Phase2 && npv > 150.0) --> @advance(phase, :Phase3),
                name => high_npv_advance, cycletime => 1.0, probability => 1.0
        end
        RDX.register_structured_species!(net, :Project)
        @prob_meta net tspan = 5 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1)
        add_structured_token!(p, RDX.FiltProjectToken(:Phase2, 100.0))   # below θ, stays
        add_structured_token!(p, RDX.FiltProjectToken(:Phase2, 200.0))   # above θ, advances
        simulate(p)
        # exactly one advanced: the npv=200 project to Phase3; the npv=100 stays Phase2
        @test phases(p) == ["Phase2", "Phase3"]
    end

    # ── determinism: predicate-selected pipeline reproduces under (model, seed) ──────────
    @testset "predicate selection is deterministic under (model, seed)" begin
        function run_pipeline(seed)
            net = advance_model()
            @prob_meta net tspan = 5 dt = 1.0
            p = ReactionNetworkProblem(net; seed = seed)
            for ph in (:Phase2, :Phase2, :Phase1, :Phase2)
                add_structured_token!(p, RDX.FiltProjectToken(ph, 100.0))
            end
            simulate(p)
            phases(p)
        end
        @test run_pipeline(11) == run_pipeline(11)
    end

    # ── deterministic bind total order (ADR 0008 inv 3): equal-priority ties broken by ──
    # ── (species, creation_index), NOT the AA Dict / random token-name order ─────────────
    @testset "equal-priority tokens bind in deterministic creation_index order (not insertion/name)" begin
        # Two Phase2 projects with distinct npv but equal (default) priority. With rate 1 only one
        # advances per tick, so WHICH one advances first must be reproducible across runs — it is
        # the lower creation_index (the first added), regardless of the tokens' random names.
        function which_advances_first(seed)
            net = @reaction_network begin
                @deterministic(1.0),
                    @select(Project, phase == :Phase2) --> @advance(phase, :Phase3),
                    name => adv, cycletime => 1.0, probability => 1.0
            end
            RDX.register_structured_species!(net, :Project)
            @prob_meta net tspan = 2 dt = 1.0
            p = ReactionNetworkProblem(net; seed = seed)
            a = RDX.FiltProjectToken(:Phase2, 111.0); add_structured_token!(p, a)  # creation_index 1
            b = RDX.FiltProjectToken(:Phase2, 222.0); add_structured_token!(p, b)  # creation_index 2
            simulate(p, 1)   # one tick: exactly one Phase2 advances
            # return the npv of whichever advanced to Phase3 (the bound-first token)
            adv = filter(t -> t.phase == :Phase3, [a, b])
            isempty(adv) ? -1.0 : first(adv).npv
        end
        # the first-added token (creation_index 1, npv 111) advances first — every time
        @test which_advances_first(1) == 111.0
        @test which_advances_first(2) == 111.0
        @test which_advances_first(999) == 111.0
    end

    # ── degenerate (no predicate) = today's kind-only bind (backward-compat) ─────────────
    @testset "a structured LHS with no @select binds by kind only (backward-compatible)" begin
        # rate 2 (two instances/tick) + cycletime 0 (instant completion) so both distinct
        # Project tokens are bound and advanced — @advance picks bound tokens in the
        # (species, creation_index) order, one per firing instance.
        net = @reaction_network begin
            @deterministic(2.0),
                @select(Project) --> @advance(phase, :Done),
                name => any_advance, cycletime => 0.0, probability => 1.0
        end
        RDX.register_structured_species!(net, :Project)
        @prob_meta net tspan = 4 dt = 1.0
        p = ReactionNetworkProblem(net; seed = 1)
        add_structured_token!(p, RDX.FiltProjectToken(:Phase1, 10.0))
        add_structured_token!(p, RDX.FiltProjectToken(:Phase2, 20.0))
        simulate(p)
        # no predicate ⇒ any Project is bindable ⇒ both advance to :Done
        @test phases(p) == ["Done", "Done"]
    end
end
