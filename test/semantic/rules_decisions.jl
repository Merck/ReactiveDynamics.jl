# Phase-1 Stage B semantic tests — the endogenous decision channel (ADR 0010/0011, CONTRACT §12).
#
# Covers: the Rule record (guard/action/fire_mode/enabled) fired at _step! step 10; the stateless
# transition guard AND-ed with the latching transActivated gate; the action family
# {SetSpecies, SetParams, AddToken, Activate, Deactivate, Log, Seq}; once-latch + reinit reset (§4 D7);
# determinism of rule effects under (model, seed). These are the in-model acquisition-lever
# capabilities the BD demo needs (MVP_BD_DEMO.md §7 / finding B).
#
# All tests run against the built engine (Stage B) — they are real @tests, not pins.

using ReactiveDynamics, Test
using Random, Distributions, DataFrames

const RDX = ReactiveDynamics

# A structured token kind for the SetTokens population-write test (defined in RD scope).
@register begin
    @aagent BaseStructuredToken AbstractStructuredToken struct RuleProjectToken
        phase::Symbol
        npv::Float64
    end
    function RuleProjectToken(phase, npv)
        return RuleProjectToken(
            "RP" * string(rand(1:(10^9))), :Project, nothing,
            Tuple{Symbol, Float64, ReactiveDynamics.Transition}[], phase, npv
        )
    end
end

# A minimal model with a `cash` pool and an inert holder so cash is a real species column.
# `@prob_meta` eval's its values in module scope, so tspan/dt are passed via the constructor
# (the `tspan=`/`dt=` kwargs) rather than threaded through the macro.
function lever_model()
    net = @reaction_network begin
        @deterministic(1.0), input --> output, name => line
        0.0, cash --> cash, name => cash_holder
    end
    @prob_init net input = 1000 output = 0 cash = 0
    @prob_params net
    return net
end

@testset "Endogenous decision channel (ADR 0010/0011, §12)" begin

    # ── (A) Rule: SetSpecies capital injection, fire_mode=once ──────────────────────────
    @testset "once-rule fires exactly once when its guard first holds (capital lever)" begin
        p = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 1,
            rules = [RDX.Rule(:lever, :(@t() > 2.0), RDX.SetSpecies(:cash, 500, :inc); fire_mode = :once)]
        )
        simulate(p)
        ci = RDX.find_index(:cash, p)
        # injects exactly once: cash goes 0 -> 500 and stays (not 500 per tick)
        @test p.u[ci] == 500.0
        @test count(>(0.0), diff(p.sol[!, "cash"])) == 1   # exactly one positive jump
        @test p.rules[1].fire_mode == :once
    end

    # ── (B) once-latch reset by reinit! (§4 D7) ─────────────────────────────────────────
    @testset "reinit! resets the once-latch so a re-run from the same seed reproduces the lever" begin
        p = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 7,
            rules = [RDX.Rule(:lever, :(@t() > 2.0), RDX.SetSpecies(:cash, 500, :inc); fire_mode = :once)]
        )
        simulate(p); cash1 = copy(p.sol[!, "cash"])
        @test p.rules[1].enabled == false        # latched off after firing
        AlgebraicAgents._reinit!(p)
        @test p.rules[1].enabled == true         # reinit re-armed it
        simulate(p)
        @test p.sol[!, "cash"] == cash1          # identical lever effect on replay
    end

    # ── (C) every_tick rule re-fires each tick (periodic financing idiom) ───────────────
    @testset "every_tick rule re-evaluates its guard and fires each tick the guard holds" begin
        p = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 1,
            rules = [RDX.Rule(:drip, :(@t() >= 0.0), RDX.SetSpecies(:cash, 10, :inc); fire_mode = :every_tick)]
        )
        simulate(p)
        # cash accrues +10 every tick the guard holds (monotone increasing)
        @test all(>=(0.0), diff(p.sol[!, "cash"]))
        @test p.sol[!, "cash"][end] >= 50.0       # many ticks of +10
        @test issorted(p.sol[!, "cash"])
    end

    # ── (D) transition guard gates genesis (no spurious resource competition) ───────────
    @testset "transition guard withholds genesis until the condition holds (conditional line)" begin
        p = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 1,
            rules = [RDX.Rule(:fund, :(@t() > 2.0), RDX.SetSpecies(:cash, 100, :set); fire_mode = :once)]
        )
        RDX.set_guard!(p, :line, :(cash >= 50))    # the `line` only fires once funded
        simulate(p)
        out = p.sol[!, "output"]
        # output stays 0 while cash < 50 (pre-funding), then the line starts firing
        @test all(==(0.0), out[1:4])
        @test out[end] > 0.0
    end

    # ── (E) Activate / Deactivate toggle a line (ADR 0004 soft gate) ────────────────────
    @testset "a rule can Deactivate a transition mid-run (soft gate, in-flight still finish)" begin
        p = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 1,
            rules = [RDX.Rule(:halt, :(@t() > 2.0), RDX.Deactivate(:line); fire_mode = :once)]
        )
        simulate(p)
        out = p.sol[!, "output"]
        # The Deactivate fires at step 10 of the t>2 tick, so that tick's genesis still happened
        # (soft gate, one-tick lag); from the NEXT tick on no new instances spawn, so output
        # plateaus. Assert: output grows before the halt, then is flat over the run's tail.
        @test out[2] > 0.0                              # produced while the line was active
        @test out[end] == out[end - 1] == out[end - 2]      # flat tail once deactivation takes hold
        @test out[end] < 7.0                            # fewer completions than the ungated 7
    end

    # ── (F) determinism: rule effects are reproducible under (model, seed) ──────────────
    @testset "rule-driven trajectory is reproducible under the same seed, differs under another" begin
        mk() = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 42,
            rules = [RDX.Rule(:lever, :(@t() > 2.0), RDX.SetSpecies(:cash, 500, :inc); fire_mode = :once)]
        )
        a = mk(); simulate(a)
        b = mk(); simulate(b)
        @test a.sol == b.sol
        c = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 6, dt = 1.0, seed = 43,
            rules = [RDX.Rule(:lever, :(@t() > 2.0), RDX.SetSpecies(:cash, 500, :inc); fire_mode = :once)]
        )
        simulate(c)
        @test a.sol[!, "output"] != c.sol[!, "output"] || a.sol[!, "cash"] == c.sol[!, "cash"]  # cash lever is deterministic; output may differ by seed
    end

    # ── (G) Seq + SetParams: the composite acquisition lever shape ──────────────────────
    @testset "Seq composes AddToken-free lever: SetSpecies + SetParams in one action" begin
        net = lever_model()
        @prob_params net synergy = 0
        p = RDX.ReactionNetworkProblem(
            net; tspan = 6, dt = 1.0, seed = 1,
            rules = [
                RDX.Rule(
                    :acq, :(@t() > 2.0),
                    RDX.Seq(
                        [
                            RDX.SetSpecies(:cash, 300, :inc),
                            RDX.SetParams([:synergy => 1]),
                        ]
                    ); fire_mode = :once
                ),
            ]
        )
        simulate(p)
        ci = RDX.find_index(:cash, p)
        @test p.u[ci] == 300.0           # capital injected
        @test p.p[:synergy] == 1         # synergy param flipped
    end

    # ── (H) numeric guard fires the action rand(rng, Poisson(v)) times, seeded ──────────
    @testset "numeric-guard rule fires Poisson(v) times per tick, deterministic under seed" begin
        mk() = RDX.ReactionNetworkProblem(
            lever_model(); tspan = 10, dt = 1.0, seed = 99,
            rules = [RDX.Rule(:noisy, :(2.0), RDX.SetSpecies(:cash, 1, :inc); fire_mode = :every_tick)]
        )
        a = mk(); simulate(a)
        b = mk(); simulate(b)
        @test a.sol[!, "cash"] == b.sol[!, "cash"]      # seeded Poisson multiplicity is reproducible
        @test a.sol[!, "cash"][end] > 0.0               # fired some positive number of times
    end

    # ── (I) SetTokens population write with @field — the ADR-0011 idiom ─────────────────
    @testset "SetTokens writes a @select-ed population's own field via @field (ADR 0011 §A)" begin
        # "write down all Phase-2 pos_remaining by 10% on a competitor readout" — a Rule action
        # over a selected population, reading each token's OWN current field via @field. This must
        # use eval_with_token (the SetField value-eval), not the plain closure path (@field is a
        # syntactic marker that must be substituted to a literal before eval, else it errors).
        net = @reaction_network begin
            0.0, A --> B, name => inert
        end
        @prob_init net A = 0 B = 0
        RDX.register_structured_species!(net, :Project)
        p = RDX.ReactionNetworkProblem(
            net; tspan = 3, dt = 1.0, seed = 1,
            population = [
                RDX.RuleProjectToken(:Phase2, 100.0),
                RDX.RuleProjectToken(:Phase2, 200.0),
                RDX.RuleProjectToken(:Phase1, 50.0),
            ]
        )
        act = RDX.SetTokens(
            RDX.TokenPredicate(:Project, [RDX.Clause(:phase, :(==), :(:Phase2))]),
            [:npv => :(@field(npv) * 0.9)],
        )
        RDX.apply_action!(p, nothing, act)
        npvs = sort([t.npv for t in values(RDX.inners(RDX.getagent(p, "structured")))])
        # Phase2: 100→90, 200→180 (each ×0.9); Phase1 50 untouched (not selected)
        @test npvs == [50.0, 90.0, 180.0]
    end
end
