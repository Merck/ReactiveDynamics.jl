# WS-B2 — `@agentize` authoring sugar (ADR 0001 / ADR 0012).
#
# `@agentize` is a THIN macro: it expands to exactly one `ReactionNetworkProblem(acs[, u0, p]; …)`
# call and reimplements no construction logic. So the acceptance bar is behavioral equivalence to
# the constructor it lowers to, PLUS the one ergonomic win — auto-naming from a bare binding.
#
# tier=T1-characterization expectedStatus=pass-now (the macro ships with WS-B2).

using ReactiveDynamics, Test          # `simulate` is reexported from AlgebraicAgents (state.jl:1)
import AlgebraicAgents                 # `getname` is qualified, matching token_filtration.jl

# A small self-contained, fully-primed SIR used by every case below (structural conservation, so any
# seed runs). It carries init/params/meta so `ReactionNetworkProblem(build_net())` is directly
# constructable — the constructor requires a `tspan` (from `@prob_meta`), so a bare unprimed net
# would `KeyError`. `primed_net` is an alias kept for the auto-name case's bare-binding readability.
function build_net()
    net = @reaction_network begin
        α * S * I, S + I --> 2I, name => I2R
        β * I, I --> R, name => R2S
    end
    @prob_init net S = 990 I = 10 R = 0
    @prob_params net α = 0.0001 β = 0.01
    @prob_meta net tspan = 50 dt = 1.0
    return net
end

primed_net() = build_net()

@testset "@agentize is thin sugar over ReactionNetworkProblem" begin
    # [agentize-equiv] the macro-built problem matches the equivalent constructor call under a
    # fixed seed — same trajectory, same t, same sol shape. This is the "reimplements no
    # construction logic" contract: @agentize net seed=1 ≡ ReactionNetworkProblem(net; seed=1).
    @testset "macro-built problem simulates identically to the constructor call (fixed seed)" begin
        p_macro = @agentize primed_net() seed = 7
        p_ctor = ReactionNetworkProblem(primed_net(); seed = 7)
        simulate(p_macro)
        simulate(p_ctor)
        @test p_macro isa ReactionNetworkProblem
        @test size(p_macro.sol) == size(p_ctor.sol)
        @test p_macro.sol[!, "S"] == p_ctor.sol[!, "S"]
        @test p_macro.sol[!, "I"] == p_ctor.sol[!, "I"]
        @test p_macro.sol[!, "R"] == p_ctor.sol[!, "R"]
        @test p_macro.t == p_ctor.t
        @test p_macro.seed == p_ctor.seed == 7
    end

    # [agentize-autoname] a bare binding auto-names the agent after the binding symbol.
    @testset "auto-naming: `@agentize net` ⇒ name = \"net\"" begin
        net = primed_net()
        p = @agentize net
        @test AlgebraicAgents.getname(p) == "net"
    end

    # [agentize-autoname-override] an explicit name= wins over the auto-name.
    @testset "explicit `name=` overrides the auto-name" begin
        net = primed_net()
        p = @agentize net name = "custom"
        @test AlgebraicAgents.getname(p) == "custom"
    end

    # [agentize-nonsymbol] a non-symbol acs expression still works and falls back to the
    # constructor's own `name` default (no auto-name injected — hygiene-safe).
    @testset "non-symbol acs falls back to the constructor's default name" begin
        p = @agentize build_net()
        default_name = AlgebraicAgents.getname(ReactionNetworkProblem(build_net()))
        @test p isa ReactionNetworkProblem
        @test AlgebraicAgents.getname(p) == default_name
    end

    # [agentize-kwarg-hygiene] a caller-LOCAL variable as a kwarg value must resolve in the CALLER's
    # scope, not RD's module scope. Regression pin for the WS-B2 hygiene bug: `args_kwargs` esc's
    # positional args but leaves kwarg VALUES unescaped, so `@agentize net seed = my_seed` used to
    # `UndefVarError: my_seed not defined in ReactiveDynamics`. This is the concrete case the ensemble
    # runner hits (it calls the constructor with `seed = <member_seed_variable>` in a loop). Literal
    # kwarg values (`seed = 7`) masked it — hence this local-variable case.
    @testset "caller-local variable as a kwarg value resolves (macro hygiene)" begin
        my_seed = 123
        net = primed_net()
        p_macro = @agentize net seed = my_seed          # must NOT throw UndefVarError
        p_ctor = ReactionNetworkProblem(primed_net(); seed = my_seed, name = "net")
        @test p_macro.seed == my_seed
        simulate(p_macro)
        simulate(p_ctor)
        @test p_macro.sol[!, "S"] == p_ctor.sol[!, "S"]   # matches the constructor with the same local
        @test p_macro.sol[!, "I"] == p_ctor.sol[!, "I"]
        @test AlgebraicAgents.getname(p_macro) == "net"   # auto-name still applies alongside esc'd kwargs
    end

    # [agentize-positional-kwargs] positional u0/p and forwarded kwargs (seed=, tspan=) reach the
    # constructor; the run matches the equivalent explicit constructor call.
    @testset "positional u0/p + forwarded kwargs reach the constructor" begin
        u0 = Dict(:S => 990, :I => 10, :R => 0)
        prm = Dict(:α => 0.0001, :β => 0.01)
        net = build_net()
        @prob_meta net tspan = 50 dt = 1.0
        p_macro = @agentize net u0 prm seed = 3
        p_ctor = ReactionNetworkProblem(build_net(), u0, prm; seed = 3, name = "net")
        simulate(p_macro)
        simulate(p_ctor)
        @test p_macro.sol[!, "S"] == p_ctor.sol[!, "S"]
        @test p_macro.sol[!, "I"] == p_ctor.sol[!, "I"]
        @test p_macro.seed == 3
        @test AlgebraicAgents.getname(p_macro) == "net"   # auto-name still applies with positionals
    end
end
